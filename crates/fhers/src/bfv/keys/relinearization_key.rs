//! Relinearization keys for the BFV encryption scheme

use std::sync::Arc;

use super::key_switching_key::KeySwitchingKey;
use crate::bfv::{
	proto::bfv::{
		KeySwitchingKey as KeySwitchingKeyProto, RelinearizationKey as RelinearizationKeyProto,
	},
	traits::TryConvertFrom,
	BfvParameters, Ciphertext, SecretKey,
};
use crate::{Error, Result};
use fhers_traits::{DeserializeParametrized, FheParametrized, Serialize};
use math::rq::{
	switcher::Switcher, traits::TryConvertFrom as TryConvertFromPoly, Poly, Representation,
};
use protobuf::{Message, MessageField};
use zeroize::Zeroize;

/// Relinearization key for the BFV encryption scheme.
/// A relinearization key is a special type of key switching key,
/// which switch from `s^2` to `s` where `s` is the secret key.
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct RelinearizationKey {
	pub(crate) ksk: KeySwitchingKey,
}

impl RelinearizationKey {
	/// Generate a [`RelinearizationKey`] from a [`SecretKey`].
	pub fn new(sk: &SecretKey) -> Result<Self> {
		Self::new_leveled(sk, 0, 0)
	}

	#[cfg(feature = "leveled_bfv")]
	#[doc(cfg(feature = "leveled_bfv"))]
	/// Generate a [`RelinearizationKey`] from a [`SecretKey`].
	pub fn new_leveled(sk: &SecretKey, ciphertext_level: usize, key_level: usize) -> Result<Self> {
		let ctx_relin_key = sk.par.ctx_at_level(key_level)?;
		let ctx_ciphertext = sk.par.ctx_at_level(ciphertext_level)?;

		if ctx_relin_key.moduli().len() == 1 {
			return Err(Error::DefaultError(
				"These parameters do not support key switching".to_string(),
			));
		}

		let mut s = Poly::try_convert_from(
			&sk.s_coefficients as &[i64],
			ctx_ciphertext,
			false,
			Representation::PowerBasis,
		)?;
		s.change_representation(Representation::Ntt);
		let mut s2 = &s * &s;
		s2.change_representation(Representation::PowerBasis);
		let switcher_up = Switcher::new(ctx_ciphertext, ctx_relin_key)?;
		let mut s2_switched_up = s2.mod_switch_to(&switcher_up)?;
		let ksk = KeySwitchingKey::new(sk, &s2_switched_up, ciphertext_level, key_level)?;
		s2.zeroize();
		s.zeroize();
		s2_switched_up.zeroize();
		Ok(Self { ksk })
	}

	/// Relinearize an "extended" ciphertext (c0, c1, c2) into a [`Ciphertext`]
	pub fn relinearizes(&self, ct: &mut Ciphertext) -> Result<()> {
		if ct.c.len() != 3 {
			Err(Error::DefaultError(
				"Only supports relinearization of ciphertext with 3 parts".to_string(),
			))
		} else if ct.level != self.ksk.ciphertext_level {
			Err(Error::DefaultError(
				"Ciphertext has incorrect level".to_string(),
			))
		} else {
			let mut c2 = ct.c[2].clone();
			c2.change_representation(Representation::PowerBasis);
			let (mut c0, mut c1) = self.relinearizes_poly(&c2)?;

			if c0.ctx() != ct.c[0].ctx() {
				c0.change_representation(Representation::PowerBasis);
				c1.change_representation(Representation::PowerBasis);
				c0.mod_switch_down_to(ct.c[0].ctx())?;
				c1.mod_switch_down_to(ct.c[1].ctx())?;
				c0.change_representation(Representation::Ntt);
				c1.change_representation(Representation::Ntt);
			}

			ct.c[0] += &c0;
			ct.c[1] += &c1;
			ct.c.truncate(2);
			Ok(())
		}
	}

	/// Relinearize using polynomials.
	pub(crate) fn relinearizes_poly(&self, c2: &Poly) -> Result<(Poly, Poly)> {
		self.ksk.key_switch(c2)
	}
}

impl From<&RelinearizationKey> for RelinearizationKeyProto {
	fn from(value: &RelinearizationKey) -> Self {
		let mut rk = RelinearizationKeyProto::new();
		rk.ksk = MessageField::some(KeySwitchingKeyProto::from(&value.ksk));
		rk
	}
}

impl TryConvertFrom<&RelinearizationKeyProto> for RelinearizationKey {
	fn try_convert_from(value: &RelinearizationKeyProto, par: &Arc<BfvParameters>) -> Result<Self> {
		if par.moduli.len() == 1 {
			Err(Error::DefaultError(
				"Invalid parameters for a relinearization key".to_string(),
			))
		} else if value.ksk.is_some() {
			Ok(RelinearizationKey {
				ksk: KeySwitchingKey::try_convert_from(value.ksk.as_ref().unwrap(), par)?,
			})
		} else {
			Err(Error::DefaultError("Invalid serialization".to_string()))
		}
	}
}

impl Serialize for RelinearizationKey {
	fn to_bytes(&self) -> Vec<u8> {
		RelinearizationKeyProto::from(self)
			.write_to_bytes()
			.unwrap()
	}
}

impl FheParametrized for RelinearizationKey {
	type Parameters = BfvParameters;
}

impl DeserializeParametrized for RelinearizationKey {
	type Error = Error;

	fn from_bytes(bytes: &[u8], par: &Arc<Self::Parameters>) -> Result<Self> {
		let rk = RelinearizationKeyProto::parse_from_bytes(bytes);
		if let Ok(rk) = rk {
			RelinearizationKey::try_convert_from(&rk, par)
		} else {
			Err(Error::DefaultError("Invalid serialization".to_string()))
		}
	}
}

#[cfg(test)]
mod tests {
	use super::RelinearizationKey;
	use crate::bfv::{
		proto::bfv::RelinearizationKey as RelinearizationKeyProto, traits::TryConvertFrom,
		BfvParameters, Ciphertext, Encoding, SecretKey,
	};
	use fhers_traits::{FheDecoder, FheDecrypter};
	use math::rq::{traits::TryConvertFrom as TryConvertFromPoly, Poly, Representation};
	use std::{error::Error, sync::Arc};

	#[test]
	fn relinearization() -> Result<(), Box<dyn Error>> {
		for params in [Arc::new(BfvParameters::default(2, 8))] {
			for _ in 0..100 {
				let sk = SecretKey::random(&params);
				let rk = RelinearizationKey::new(&sk)?;

				let ctx = params.ctx_at_level(0)?;
				let mut s = Poly::try_convert_from(
					&sk.s_coefficients as &[i64],
					ctx,
					false,
					Representation::PowerBasis,
				)
				.map_err(crate::Error::MathError)?;
				s.change_representation(Representation::Ntt);
				let s2 = &s * &s;
				// Let's generate manually an "extended" ciphertext (c0 = e - c1 * s - c2 * s^2,
				// c1, c2) encrypting 0.
				let mut c2 = Poly::random(ctx, Representation::Ntt);
				let c1 = Poly::random(ctx, Representation::Ntt);
				let mut c0 = Poly::small(ctx, Representation::PowerBasis, 16)?;
				c0.change_representation(Representation::Ntt);
				c0 -= &(&c1 * &s);
				c0 -= &(&c2 * &s2);
				let mut ct = Ciphertext::new(vec![c0.clone(), c1.clone(), c2.clone()], &params)?;

				// Relinearize the extended ciphertext!
				rk.relinearizes(&mut ct)?;
				assert_eq!(ct.c.len(), 2);

				// Check that the relinearization by polynomials works the same way
				c2.change_representation(Representation::PowerBasis);
				let (mut c0r, mut c1r) = rk.relinearizes_poly(&c2)?;
				c0r.change_representation(Representation::PowerBasis);
				c0r.mod_switch_down_to(c0.ctx())?;
				c1r.change_representation(Representation::PowerBasis);
				c1r.mod_switch_down_to(c1.ctx())?;
				c0r.change_representation(Representation::Ntt);
				c1r.change_representation(Representation::Ntt);
				assert_eq!(ct, Ciphertext::new(vec![&c0 + &c0r, &c1 + &c1r], &params)?);

				// Print the noise and decrypt
				println!("Noise: {}", unsafe { sk.measure_noise(&ct)? });
				let pt = sk.try_decrypt(&ct)?;
				assert!(Vec::<u64>::try_decode(&pt, Encoding::poly()).is_ok_and(|v| v == &[0u64; 8]))
			}
		}
		Ok(())
	}

	#[test]
	fn relinearization_leveled() -> Result<(), Box<dyn Error>> {
		for params in [Arc::new(BfvParameters::default(5, 8))] {
			for ciphertext_level in 0..3 {
				for key_level in 0..ciphertext_level {
					for _ in 0..10 {
						let sk = SecretKey::random(&params);
						let rk = RelinearizationKey::new_leveled(&sk, ciphertext_level, key_level)?;

						let ctx = params.ctx_at_level(ciphertext_level)?;
						let mut s = Poly::try_convert_from(
							&sk.s_coefficients as &[i64],
							ctx,
							false,
							Representation::PowerBasis,
						)
						.map_err(crate::Error::MathError)?;
						s.change_representation(Representation::Ntt);
						let s2 = &s * &s;
						// Let's generate manually an "extended" ciphertext (c0 = e - c1 * s - c2 *
						// s^2, c1, c2) encrypting 0.
						let mut c2 = Poly::random(ctx, Representation::Ntt);
						let c1 = Poly::random(ctx, Representation::Ntt);
						let mut c0 = Poly::small(ctx, Representation::PowerBasis, 16)?;
						c0.change_representation(Representation::Ntt);
						c0 -= &(&c1 * &s);
						c0 -= &(&c2 * &s2);
						let mut ct =
							Ciphertext::new(vec![c0.clone(), c1.clone(), c2.clone()], &params)?;

						// Relinearize the extended ciphertext!
						rk.relinearizes(&mut ct)?;
						assert_eq!(ct.c.len(), 2);

						// Check that the relinearization by polynomials works the same way
						c2.change_representation(Representation::PowerBasis);
						let (mut c0r, mut c1r) = rk.relinearizes_poly(&c2)?;
						c0r.change_representation(Representation::PowerBasis);
						c0r.mod_switch_down_to(c0.ctx())?;
						c1r.change_representation(Representation::PowerBasis);
						c1r.mod_switch_down_to(c1.ctx())?;
						c0r.change_representation(Representation::Ntt);
						c1r.change_representation(Representation::Ntt);
						assert_eq!(ct, Ciphertext::new(vec![&c0 + &c0r, &c1 + &c1r], &params)?);

						// Print the noise and decrypt
						println!("Noise: {}", unsafe { sk.measure_noise(&ct)? });
						let pt = sk.try_decrypt(&ct)?;
						assert!(Vec::<u64>::try_decode(&pt, Encoding::poly())
							.is_ok_and(|v| v == &[0u64; 8]))
					}
				}
			}
		}
		Ok(())
	}

	#[test]
	fn proto_conversion() -> Result<(), Box<dyn Error>> {
		for params in [Arc::new(BfvParameters::default(2, 8))] {
			let sk = SecretKey::random(&params);
			let rk = RelinearizationKey::new(&sk)?;
			let proto = RelinearizationKeyProto::from(&rk);
			assert_eq!(rk, RelinearizationKey::try_convert_from(&proto, &params)?);
		}
		Ok(())
	}
}
