//! Public keys for the BFV encryption scheme

use crate::bfv::traits::TryConvertFrom;
use crate::bfv::{
	proto::bfv::{Ciphertext as CiphertextProto, PublicKey as PublicKeyProto},
	BfvParameters, Ciphertext, Encoding, Plaintext,
};
use crate::{Error, Result};
use fhe_math::rq::{Poly, Representation};
use fhe_traits::{DeserializeParametrized, FheEncrypter, FheParametrized, Serialize};
use protobuf::{Message, MessageField};
use std::sync::Arc;
use zeroize::Zeroizing;

use super::SecretKey;

/// Public key for the BFV encryption scheme.
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct PublicKey {
	pub(crate) par: Arc<BfvParameters>,
	pub(crate) c: Ciphertext,
}

impl PublicKey {
	/// Generate a new [`PublicKey`] from a [`SecretKey`].
	pub fn new(sk: &SecretKey) -> Result<Self> {
		let zero = Plaintext::zero(Encoding::poly(), &sk.par)?;
		let mut c: Ciphertext = sk.try_encrypt(&zero)?;
		// The polynomials of a public key should not allow for variable time
		// computation.
		c.c.iter_mut()
			.for_each(|p| p.disallow_variable_time_computations());
		Ok(Self {
			par: sk.par.clone(),
			c,
		})
	}
}

impl FheParametrized for PublicKey {
	type Parameters = BfvParameters;
}

impl FheEncrypter<Plaintext, Ciphertext> for PublicKey {
	type Error = Error;

	fn try_encrypt(&self, pt: &Plaintext) -> Result<Ciphertext> {
		let mut ct = self.c.clone();

		#[cfg(feature = "leveled_bfv")]
		while ct.level != pt.level {
			ct.mod_switch_to_next_level();
		}

		let ctx = self.par.ctx_at_level(ct.level)?;
		let u = Zeroizing::new(
			Poly::small(ctx, Representation::Ntt, self.par.variance).map_err(Error::MathError)?,
		);
		let e1 = Zeroizing::new(
			Poly::small(ctx, Representation::Ntt, self.par.variance).map_err(Error::MathError)?,
		);
		let e2 = Zeroizing::new(
			Poly::small(ctx, Representation::Ntt, self.par.variance).map_err(Error::MathError)?,
		);

		let m = Zeroizing::new(pt.to_poly()?);
		let mut c0 = u.as_ref() * &ct.c[0];
		c0 += &e1;
		c0 += &m;
		let mut c1 = u.as_ref() * &ct.c[1];
		c1 += &e2;

		// It is now safe to enable variable time computations.
		unsafe {
			c0.allow_variable_time_computations();
			c1.allow_variable_time_computations()
		}

		Ok(Ciphertext {
			par: self.par.clone(),
			seed: None,
			c: vec![c0, c1],
			level: ct.level,
		})
	}
}

impl From<&PublicKey> for PublicKeyProto {
	fn from(pk: &PublicKey) -> Self {
		let mut proto = PublicKeyProto::new();
		proto.c = MessageField::some(CiphertextProto::from(&pk.c));
		proto
	}
}

impl Serialize for PublicKey {
	fn to_bytes(&self) -> Vec<u8> {
		PublicKeyProto::from(self).write_to_bytes().unwrap()
	}
}

impl DeserializeParametrized for PublicKey {
	type Error = Error;

	fn from_bytes(bytes: &[u8], par: &Arc<Self::Parameters>) -> Result<Self> {
		let proto =
			PublicKeyProto::parse_from_bytes(bytes).map_err(|_| Error::SerializationError)?;
		if proto.c.is_some() {
			let mut c = Ciphertext::try_convert_from(&proto.c.unwrap(), par)?;
			if c.level != 0 {
				Err(Error::SerializationError)
			} else {
				// The polynomials of a public key should not allow for variable time
				// computation.
				c.c.iter_mut()
					.for_each(|p| p.disallow_variable_time_computations());
				Ok(Self {
					par: par.clone(),
					c,
				})
			}
		} else {
			Err(Error::SerializationError)
		}
	}
}

#[cfg(test)]
mod tests {
	use super::PublicKey;
	use crate::bfv::{parameters::BfvParameters, Encoding, Plaintext, SecretKey};
	use fhe_traits::{DeserializeParametrized, FheDecrypter, FheEncoder, FheEncrypter, Serialize};
	use std::{error::Error, sync::Arc};

	#[test]
	fn keygen() -> Result<(), Box<dyn Error>> {
		let params = Arc::new(BfvParameters::default(1, 8));
		let sk = SecretKey::random(&params);
		let pk = PublicKey::new(&sk)?;
		assert_eq!(pk.par, params);
		assert_eq!(
			sk.try_decrypt(&pk.c)?,
			Plaintext::zero(Encoding::poly(), &params)?
		);
		Ok(())
	}

	#[test]
	fn encrypt_decrypt() -> Result<(), Box<dyn Error>> {
		for params in [
			Arc::new(BfvParameters::default(1, 8)),
			Arc::new(BfvParameters::default(6, 8)),
		] {
			for level in 0..params.max_level() {
				for _ in 0..20 {
					let sk = SecretKey::random(&params);
					let pk = PublicKey::new(&sk)?;

					let pt = Plaintext::try_encode(
						&params.plaintext.random_vec(params.degree()) as &[u64],
						Encoding::poly_at_level(level),
						&params,
					)?;
					let ct = pk.try_encrypt(&pt)?;
					let pt2 = sk.try_decrypt(&ct);

					println!("Noise: {}", unsafe { sk.measure_noise(&ct)? });
					assert!(pt2.is_ok_and(|pt2| pt2 == &pt));
				}
			}
		}

		Ok(())
	}

	#[test]
	fn test_serialize() -> Result<(), Box<dyn Error>> {
		for params in [
			Arc::new(BfvParameters::default(1, 8)),
			Arc::new(BfvParameters::default(6, 8)),
		] {
			let sk = SecretKey::random(&params);
			let pk = PublicKey::new(&sk)?;
			let bytes = pk.to_bytes();
			assert_eq!(pk, PublicKey::from_bytes(&bytes, &params)?);
		}
		Ok(())
	}
}
