//! Galois keys for the BFV encryption scheme

use super::key_switching_key::KeySwitchingKey;
use crate::bfv::{
	proto::bfv::{GaloisKey as GaloisKeyProto, KeySwitchingKey as KeySwitchingKeyProto},
	traits::TryConvertFrom,
	BfvParameters, Ciphertext, SecretKey,
};
use crate::{Error, Result};
use math::rq::{Poly, Representation};
use protobuf::MessageField;
use std::sync::Arc;
use zeroize::Zeroize;

/// Galois key for the BFV encryption scheme.
/// A Galois key is a special type of key switching key,
/// which switch from `s(x^i)` to `s(x)` where `s(x)` is the secret key.
#[derive(Debug, PartialEq, Eq)]
pub struct GaloisKey {
	pub(crate) exponent: usize,
	ksk: KeySwitchingKey,
}

impl GaloisKey {
	/// Generate a [`GaloisKey`] from a [`SecretKey`].
	pub fn new(sk: &SecretKey, exponent: usize) -> Result<Self> {
		let exponent = exponent % (2 * sk.par.degree());
		if exponent == 0 {
			return Err(Error::DefaultError("Invalid exponent".to_string()));
		}

		let mut s_sub = sk.s[0].substitute(exponent)?;
		s_sub.change_representation(Representation::PowerBasis);
		let ksk = KeySwitchingKey::new(sk, &s_sub)?;
		s_sub.zeroize();
		Ok(Self { exponent, ksk })
	}

	/// Relinearize a [`Ciphertext`] using the [`GaloisKey`]
	pub fn relinearize(&self, ct: &Ciphertext) -> Result<Ciphertext> {
		assert_eq!(ct.par, self.ksk.par);
		assert_eq!(ct.c.len(), 2);

		let mut c0 = ct.c[0].substitute(self.exponent)?;
		let mut c1 = Poly::zero(&self.ksk.par.ctx, Representation::Ntt);
		unsafe { c1.allow_variable_time_computations() }

		let mut c2 = ct.c[1].substitute(self.exponent)?;
		c2.change_representation(Representation::PowerBasis);
		self.ksk.key_switch(&c2, &mut c0, &mut c1)?;

		Ok(Ciphertext {
			par: ct.par.clone(),
			seed: None,
			c: vec![c0, c1],
		})
	}
}

impl From<&GaloisKey> for GaloisKeyProto {
	fn from(value: &GaloisKey) -> Self {
		let mut gk = GaloisKeyProto::new();
		gk.exponent = value.exponent as u32;
		gk.ksk = MessageField::some(KeySwitchingKeyProto::from(&value.ksk));
		gk
	}
}

impl TryConvertFrom<&GaloisKeyProto> for GaloisKey {
	fn try_convert_from(value: &GaloisKeyProto, par: &Arc<BfvParameters>) -> Result<Self> {
		let exponent = (value.exponent as usize) % (2 * par.degree());
		if exponent & 1 == 0 {
			return Err(Error::DefaultError("Invalid exponent".to_string()));
		}
		if par.ciphertext_moduli.len() == 1 {
			Err(Error::DefaultError(
				"Invalid parameters for a relinearization key".to_string(),
			))
		} else if value.ksk.is_some() {
			let ksk = KeySwitchingKey::try_convert_from(value.ksk.as_ref().unwrap(), par)?;
			Ok(GaloisKey { exponent, ksk })
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
	fn test_relinearization() -> Result<(), Box<dyn Error>> {
		for params in [Arc::new(BfvParameters::default(2))] {
			for _ in 0..50 {
				let mut sk = SecretKey::random(&params);
				let v = params.plaintext.random_vec(params.degree());
				let row_size = params.degree() >> 1;

				let pt = Plaintext::try_encode(&v as &[u64], Encoding::Simd, &params)?;
				let ct = sk.try_encrypt(&pt)?;

				for i in 1..16 {
					if i & 1 == 0 {
						assert!(GaloisKey::new(&sk, i).is_err())
					} else {
						let gk = GaloisKey::new(&sk, i)?;
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
							assert_eq!(&Vec::<u64>::try_decode(&pt, Encoding::Simd)?, &expected)
						} else if i == params.degree() * 2 - 1 {
							let pt = sk.try_decrypt(&ct2)?;

							// The expected result has its rows flipped
							let mut expected = vec![0u64; params.degree()];
							expected[..row_size].copy_from_slice(&v[row_size..]);
							expected[row_size..].copy_from_slice(&v[..row_size]);
							assert_eq!(&Vec::<u64>::try_decode(&pt, Encoding::Simd)?, &expected)
						}
					}
				}
			}
		}
		Ok(())
	}

	#[test]
	fn test_proto_conversion() -> Result<(), Box<dyn Error>> {
		for params in [Arc::new(BfvParameters::default(2))] {
			let sk = SecretKey::random(&params);
			let gk = GaloisKey::new(&sk, 9)?;
			let proto = GaloisKeyProto::from(&gk);
			assert_eq!(gk, GaloisKey::try_convert_from(&proto, &params)?);
		}
		Ok(())
	}
}
