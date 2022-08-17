//! Relinearization keys for the BFV encryption scheme

use std::sync::Arc;

use super::key_switching_key::KeySwitchingKey;
use crate::{traits::TryConvertFrom, BfvParameters, Ciphertext, Error, Result, SecretKey};
use fhers_protos::protos::bfv::{
	KeySwitchingKey as KeySwitchingKeyProto, RelinearizationKey as RelinearizationKeyProto,
};
use math::rq::{Poly, Representation};
use protobuf::MessageField;
use zeroize::Zeroize;

/// Relinearization key for the BFV encryption scheme.
/// A relinearization key is a special type of key switching key,
/// which switch from `s^2` to `s` where `s` is the secret key.
#[derive(Debug, PartialEq, Eq)]
pub struct RelinearizationKey {
	ksk: KeySwitchingKey,
}

impl RelinearizationKey {
	/// Generate a [`RelinearizationKey`] from a [`SecretKey`].
	pub fn new(sk: &SecretKey) -> Result<Self> {
		let mut s_squared = sk.s[1].clone();
		s_squared.change_representation(Representation::PowerBasis);
		let ksk = KeySwitchingKey::new(sk, &s_squared)?;
		s_squared.zeroize();
		Ok(Self { ksk })
	}

	/// Relinearize an "extended" ciphertext (c0, c1, c2) into a [`Ciphertext`]
	pub fn relinearizes(&self, ct: &Ciphertext) -> Result<Ciphertext> {
		if ct.c.len() != 3 {
			Err(Error::DefaultError(
				"Only supports relinearization of ciphertext with 3 parts".to_string(),
			))
		} else {
			let mut c2 = ct.c[2].clone();
			c2.change_representation(Representation::PowerBasis);
			let mut c0 = ct.c[0].clone();
			let mut c1 = ct.c[1].clone();
			self.relinearizes_with_poly(&c2, &mut c0, &mut c1)?;
			Ok(Ciphertext {
				par: ct.par.clone(),
				seed: None,
				c: vec![c0, c1],
			})
		}
	}

	/// Relinearize using polynomials.
	pub(crate) fn relinearizes_with_poly(
		&self,
		c2: &Poly,
		c0: &mut Poly,
		c1: &mut Poly,
	) -> Result<()> {
		if c2.representation() != &Representation::PowerBasis {
			Err(Error::DefaultError(
				"Incorrect representation for c2".to_string(),
			))
		} else {
			self.ksk.key_switch(c2, c0, c1)
		}
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
		if par.ciphertext_moduli.len() == 1 {
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
	use crate::{
		traits::{Decoder, Decryptor, TryConvertFrom},
		BfvParameters, Ciphertext, Encoding, SecretKey,
	};
	use fhers_protos::protos::bfv::RelinearizationKey as RelinearizationKeyProto;
	use math::rq::{Poly, Representation};
	use std::{error::Error, sync::Arc};

	#[test]
	fn test_relinearization() -> Result<(), Box<dyn Error>> {
		for params in [Arc::new(BfvParameters::default(2))] {
			for _ in 0..100 {
				let mut sk = SecretKey::random(&params);
				let rk = RelinearizationKey::new(&sk)?;

				// Let's generate manually an "extended" ciphertext (c0 = e - c1 * s - c2 * s^2,
				// c1, c2) encrypting 0.
				let mut c2 = Poly::random(&params.ctx, Representation::Ntt);
				let mut c1 = Poly::random(&params.ctx, Representation::Ntt);
				let mut c0 = Poly::small(&params.ctx, Representation::PowerBasis, 16)?;
				c0.change_representation(Representation::Ntt);
				c0 -= &(&c1 * &sk.s[0]);
				c0 -= &(&c2 * &sk.s[1]);
				let ct = Ciphertext {
					par: params.clone(),
					seed: None,
					c: vec![c0.clone(), c1.clone(), c2.clone()],
				};

				// Relinearize the extended ciphertext!
				let ct_relinearized = rk.relinearizes(&ct)?;
				assert_eq!(ct_relinearized.c.len(), 2);

				// Check that the relinearization by polynomials works the same way
				c2.change_representation(Representation::PowerBasis);
				rk.relinearizes_with_poly(&c2, &mut c0, &mut c1)?;
				assert_eq!(
					ct_relinearized,
					Ciphertext {
						par: params.clone(),
						seed: None,
						c: vec![c0, c1]
					}
				);

				// Print the noise and decrypt
				println!("Noise: {}", unsafe { sk.measure_noise(&ct_relinearized)? });
				let pt = sk.decrypt(&ct_relinearized)?;
				assert!(Vec::<u64>::try_decode(&pt, Encoding::Poly).is_ok_and(|v| v == &[0u64; 8]))
			}
		}
		Ok(())
	}

	#[test]
	fn test_proto_conversion() -> Result<(), Box<dyn Error>> {
		for params in [Arc::new(BfvParameters::default(2))] {
			let sk = SecretKey::random(&params);
			let rk = RelinearizationKey::new(&sk)?;
			let proto = RelinearizationKeyProto::from(&rk);
			assert_eq!(rk, RelinearizationKey::try_convert_from(&proto, &params)?);
		}
		Ok(())
	}
}
