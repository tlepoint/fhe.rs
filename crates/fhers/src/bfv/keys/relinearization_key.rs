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
use math::rq::{traits::TryConvertFrom as TryConvertFromPoly, Poly, Representation};
use protobuf::MessageField;
use zeroize::Zeroize;

/// Relinearization key for the BFV encryption scheme.
/// A relinearization key is a special type of key switching key,
/// which switch from `s^2` to `s` where `s` is the secret key.
#[derive(Debug, PartialEq, Eq)]
pub struct RelinearizationKey {
	pub(crate) ksk: KeySwitchingKey,
}

impl RelinearizationKey {
	/// Generate a [`RelinearizationKey`] from a [`SecretKey`].
	pub fn new(sk: &SecretKey, ciphertext_level: usize, relin_key_level: usize) -> Result<Self> {
		let ctx_relin_key = sk.par.ctx_at_level(relin_key_level)?;

		if ctx_relin_key.moduli().len() == 1 {
			return Err(Error::DefaultError(
				"These parameters do not support key switching".to_string(),
			));
		}

		let mut s = Poly::try_convert_from(
			&sk.s_coefficients as &[i64],
			&ctx_relin_key,
			false,
			Representation::PowerBasis,
		)?;
		s.change_representation(Representation::Ntt);
		let mut s2 = &s * &s;
		s2.change_representation(Representation::PowerBasis);
		let ksk = KeySwitchingKey::new(sk, &s2, ciphertext_level, relin_key_level)?;
		s2.zeroize();
		s.zeroize();
		Ok(Self { ksk })
	}

	/// Relinearize an "extended" ciphertext (c0, c1, c2) into a [`Ciphertext`]
	pub fn relinearizes(&self, ct: &Ciphertext) -> Result<Ciphertext> {
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
			let (mut c0, mut c1) = self.relinearizes_with_poly(&c2)?;
			c0 += &ct.c[0];
			c1 += &ct.c[1];
			Ok(Ciphertext {
				par: ct.par.clone(),
				seed: None,
				c: vec![c0, c1],
				level: self.ksk.ciphertext_level,
			})
		}
	}

	/// Relinearize using polynomials.
	pub(crate) fn relinearizes_with_poly(&self, c2: &Poly) -> Result<(Poly, Poly)> {
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
	fn test_relinearization() -> Result<(), Box<dyn Error>> {
		for params in [Arc::new(BfvParameters::default(2))] {
			for _ in 0..100 {
				let sk = SecretKey::random(&params);
				let rk = RelinearizationKey::new(&sk, 0, 0)?;

				let ctx = params.ctx_at_level(0)?;
				let mut s = Poly::try_convert_from(
					&sk.s_coefficients as &[i64],
					&ctx,
					false,
					Representation::PowerBasis,
				)
				.map_err(crate::Error::MathError)?;
				s.change_representation(Representation::Ntt);
				let s2 = &s * &s;
				// Let's generate manually an "extended" ciphertext (c0 = e - c1 * s - c2 * s^2,
				// c1, c2) encrypting 0.
				let mut c2 = Poly::random(&ctx, Representation::Ntt);
				let c1 = Poly::random(&ctx, Representation::Ntt);
				let mut c0 = Poly::small(&ctx, Representation::PowerBasis, 16)?;
				c0.change_representation(Representation::Ntt);
				c0 -= &(&c1 * &s);
				c0 -= &(&c2 * &s2);
				let ct = Ciphertext {
					par: params.clone(),
					seed: None,
					c: vec![c0.clone(), c1.clone(), c2.clone()],
					level: 0,
				};

				// Relinearize the extended ciphertext!
				let ct_relinearized = rk.relinearizes(&ct)?;
				assert_eq!(ct_relinearized.c.len(), 2);

				// Check that the relinearization by polynomials works the same way
				c2.change_representation(Representation::PowerBasis);
				let (c0r, c1r) = rk.relinearizes_with_poly(&c2)?;
				assert_eq!(
					ct_relinearized,
					Ciphertext {
						par: params.clone(),
						seed: None,
						c: vec![&c0 + &c0r, &c1 + &c1r],
						level: 0
					}
				);

				// Print the noise and decrypt
				println!("Noise: {}", unsafe { sk.measure_noise(&ct_relinearized)? });
				let pt = sk.try_decrypt(&ct_relinearized)?;
				assert!(Vec::<u64>::try_decode(&pt, Encoding::PolyLeveled(0))
					.is_ok_and(|v| v == &[0u64; 8]))
			}
		}
		Ok(())
	}

	#[test]
	fn test_proto_conversion() -> Result<(), Box<dyn Error>> {
		for params in [Arc::new(BfvParameters::default(2))] {
			let sk = SecretKey::random(&params);
			let rk = RelinearizationKey::new(&sk, 0, 0)?;
			let proto = RelinearizationKeyProto::from(&rk);
			assert_eq!(rk, RelinearizationKey::try_convert_from(&proto, &params)?);
		}
		Ok(())
	}
}
