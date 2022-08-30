//! Galois keys for the BFV encryption scheme

use super::key_switching_key::KeySwitchingKey;
use crate::bfv::{
	proto::bfv::{GaloisKey as GaloisKeyProto, KeySwitchingKey as KeySwitchingKeyProto},
	traits::TryConvertFrom,
	BfvParameters, Ciphertext, SecretKey,
};
use crate::{Error, Result};
use math::rq::{
	switcher::Switcher, traits::TryConvertFrom as TryConvertFromPoly, Poly, Representation,
	SubstitutionExponent,
};
use protobuf::MessageField;
use std::sync::Arc;
use zeroize::Zeroize;

/// Galois key for the BFV encryption scheme.
/// A Galois key is a special type of key switching key,
/// which switch from `s(x^i)` to `s(x)` where `s(x)` is the secret key.
#[derive(Debug, PartialEq, Eq)]
pub struct GaloisKey {
	pub(crate) element: SubstitutionExponent,
	pub(crate) ksk: KeySwitchingKey,
}

impl GaloisKey {
	/// Generate a [`GaloisKey`] from a [`SecretKey`].
	pub fn new(
		sk: &SecretKey,
		exponent: usize,
		ciphertext_level: usize,
		galois_key_level: usize,
	) -> Result<Self> {
		let ctx_galois_key = sk.par.ctx_at_level(galois_key_level)?;
		let ctx_ciphertext = sk.par.ctx_at_level(ciphertext_level)?;

		let ciphertext_exponent =
			SubstitutionExponent::new(ctx_ciphertext, exponent).map_err(Error::MathError)?;

		let switcher_up = Switcher::new(ctx_ciphertext, ctx_galois_key)?;
		let mut s = Poly::try_convert_from(
			&sk.s_coefficients as &[i64],
			ctx_ciphertext,
			false,
			Representation::PowerBasis,
		)?;
		let mut s_sub = s.substitute(&ciphertext_exponent)?;
		let mut s_sub_switched_up = s_sub.mod_switch_to(&switcher_up)?;
		s_sub_switched_up.change_representation(Representation::PowerBasis);

		let ksk = KeySwitchingKey::new(sk, &s_sub_switched_up, ciphertext_level, galois_key_level)?;
		s_sub.zeroize();
		s_sub_switched_up.zeroize();
		s.zeroize();

		Ok(Self {
			element: ciphertext_exponent,
			ksk,
		})
	}

	/// Relinearize a [`Ciphertext`] using the [`GaloisKey`]
	pub fn relinearize(&self, ct: &Ciphertext) -> Result<Ciphertext> {
		// assert_eq!(ct.par, self.ksk.par);
		assert_eq!(ct.c.len(), 2);

		let mut c2 = ct.c[1].substitute(&self.element)?;
		c2.change_representation(Representation::PowerBasis);
		let (mut c0, mut c1) = self.ksk.key_switch(&c2)?;

		if c0.ctx() != ct.c[0].ctx() {
			c0.change_representation(Representation::PowerBasis);
			c1.change_representation(Representation::PowerBasis);
			c0.mod_switch_down_to(ct.c[0].ctx())?;
			c1.mod_switch_down_to(ct.c[1].ctx())?;
			c0.change_representation(Representation::Ntt);
			c1.change_representation(Representation::Ntt);
		}

		c0 += ct.c[0].substitute(&self.element)?;

		Ok(Ciphertext {
			par: ct.par.clone(),
			seed: None,
			c: vec![c0, c1],
			level: self.ksk.ciphertext_level,
		})
	}
}

impl From<&GaloisKey> for GaloisKeyProto {
	fn from(value: &GaloisKey) -> Self {
		let mut gk = GaloisKeyProto::new();
		gk.exponent = value.element.exponent as u32;
		gk.ksk = MessageField::some(KeySwitchingKeyProto::from(&value.ksk));
		gk
	}
}

impl TryConvertFrom<&GaloisKeyProto> for GaloisKey {
	fn try_convert_from(value: &GaloisKeyProto, par: &Arc<BfvParameters>) -> Result<Self> {
		if par.moduli.len() == 1 {
			Err(Error::DefaultError(
				"Invalid parameters for a relinearization key".to_string(),
			))
		} else if value.ksk.is_some() {
			let ksk = KeySwitchingKey::try_convert_from(value.ksk.as_ref().unwrap(), par)?;

			let ctx = par.ctx_at_level(ksk.ciphertext_level)?;
			let element = SubstitutionExponent::new(ctx, value.exponent as usize)
				.map_err(Error::MathError)?;

			Ok(GaloisKey { element, ksk })
		} else {
			Err(Error::DefaultError("Invalid serialization".to_string()))
		}
	}
}

#[cfg(test)]
mod tests {
	use super::GaloisKey;
	use crate::bfv::{
		proto::bfv::GaloisKey as GaloisKeyProto, traits::TryConvertFrom, BfvParameters, Encoding,
		Plaintext, SecretKey,
	};
	use fhers_traits::{FheDecoder, FheDecrypter, FheEncoder, FheEncrypter};
	use std::{error::Error, sync::Arc};

	#[test]
	fn relinearization() -> Result<(), Box<dyn Error>> {
		for params in [
			Arc::new(BfvParameters::default(2, 8)),
			Arc::new(BfvParameters::default(3, 8)),
		] {
			for _ in 0..30 {
				let sk = SecretKey::random(&params);
				let v = params.plaintext.random_vec(params.degree());
				let row_size = params.degree() >> 1;

				let pt = Plaintext::try_encode(&v as &[u64], Encoding::simd(), &params)?;
				let ct = sk.try_encrypt(&pt)?;

				for i in 1..2 * params.degree() {
					if i & 1 == 0 {
						assert!(GaloisKey::new(&sk, i, 0, 0).is_err())
					} else {
						let gk = GaloisKey::new(&sk, i, 0, 0)?;
						let ct2 = gk.relinearize(&ct)?;
						println!("Noise: {}", unsafe { sk.measure_noise(&ct2)? });

						if i == 3 {
							let pt = sk.try_decrypt(&ct2)?;

							// The expected result is rotated one on the left
							let mut expected = vec![0u64; params.degree()];
							expected[..row_size - 1].copy_from_slice(&v[1..row_size]);
							expected[row_size - 1] = v[0];
							expected[row_size..2 * row_size - 1]
								.copy_from_slice(&v[row_size + 1..]);
							expected[2 * row_size - 1] = v[row_size];
							assert_eq!(&Vec::<u64>::try_decode(&pt, Encoding::simd())?, &expected)
						} else if i == params.degree() * 2 - 1 {
							let pt = sk.try_decrypt(&ct2)?;

							// The expected result has its rows swapped
							let mut expected = vec![0u64; params.degree()];
							expected[..row_size].copy_from_slice(&v[row_size..]);
							expected[row_size..].copy_from_slice(&v[..row_size]);
							assert_eq!(&Vec::<u64>::try_decode(&pt, Encoding::simd())?, &expected)
						}
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
			let gk = GaloisKey::new(&sk, 9, 0, 0)?;
			let proto = GaloisKeyProto::from(&gk);
			assert_eq!(gk, GaloisKey::try_convert_from(&proto, &params)?);
		}
		Ok(())
	}
}
