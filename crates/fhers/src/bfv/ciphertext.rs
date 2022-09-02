//! Ciphertext type in the BFV encryption scheme.

use crate::bfv::{
	parameters::BfvParameters, proto::bfv::Ciphertext as CiphertextProto, traits::TryConvertFrom,
};
use crate::{Error, Result};
use fhers_traits::{
	DeserializeParametrized, DeserializeWithContext, FheCiphertext, FheParametrized, Serialize,
};
use math::rq::{Poly, Representation};
use protobuf::Message;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use std::sync::Arc;

/// A ciphertext encrypting a plaintext.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Ciphertext {
	/// The parameters of the underlying BFV encryption scheme.
	pub(crate) par: Arc<BfvParameters>,

	/// The seed that generated the polynomial c1 in a fresh ciphertext.
	pub(crate) seed: Option<<ChaCha8Rng as SeedableRng>::Seed>,

	/// The ciphertext elements.
	pub(crate) c: Vec<Poly>,

	/// The ciphertext level
	pub(crate) level: usize,
}

impl Ciphertext {
	#[cfg(feature = "leveled_bfv")]
	#[doc(cfg(feature = "leveled_bfv"))]
	/// Modulo switch the ciphertext to the last level.
	pub fn mod_switch_to_last_level(&mut self) {
		self.level = self.par.max_level();
		let last_ctx = self.par.ctx_at_level(self.level).unwrap();
		self.seed = None;
		self.c.iter_mut().for_each(|ci| {
			if ci.ctx() != last_ctx {
				ci.change_representation(Representation::PowerBasis);
				assert!(ci.mod_switch_down_to(last_ctx).is_ok());
				ci.change_representation(Representation::Ntt);
			}
		});
	}

	#[cfg(feature = "leveled_bfv")]
	#[doc(cfg(feature = "leveled_bfv"))]
	/// Modulo switch the ciphertext to the next level.
	pub fn mod_switch_to_next_level(&mut self) {
		if self.level < self.par.max_level() {
			self.seed = None;
			self.c.iter_mut().for_each(|ci| {
				ci.change_representation(Representation::PowerBasis);
				assert!(ci.mod_switch_down_next().is_ok());
				ci.change_representation(Representation::Ntt);
			});
			self.level += 1
		}
	}

	/// Create a ciphertext from a vector of polynomials.
	/// A ciphertext must contain at least two polynomials, and all polynomials
	/// must be in Ntt representation and with the same context.
	pub fn new(c: Vec<Poly>, par: &Arc<BfvParameters>) -> Result<Self> {
		if c.len() < 2 {
			return Err(Error::TooFewValues(c.len(), 2));
		}

		let ctx = c[0].ctx();
		let level = par.level_of_ctx(ctx)?;

		// Check that all polynomials have the expected representation and context.
		for ci in c.iter() {
			if ci.representation() != &Representation::Ntt {
				return Err(Error::MathError(math::Error::IncorrectRepresentation(
					ci.representation().clone(),
					Representation::Ntt,
				)));
			}
			if ci.ctx() != ctx {
				return Err(Error::MathError(math::Error::InvalidContext));
			}
		}

		Ok(Self {
			par: par.clone(),
			seed: None,
			c,
			level,
		})
	}

	/// Get the i-th polynomial of the ciphertext.
	pub fn get(&self, i: usize) -> Option<&Poly> {
		self.c.get(i)
	}
}

impl FheCiphertext for Ciphertext {}

impl FheParametrized for Ciphertext {
	type Parameters = BfvParameters;
}

impl Serialize for Ciphertext {
	fn to_bytes(&self) -> Vec<u8> {
		CiphertextProto::from(self).write_to_bytes().unwrap()
	}
}

impl DeserializeParametrized for Ciphertext {
	fn from_bytes(bytes: &[u8], par: &Arc<BfvParameters>) -> Result<Self> {
		if let Ok(ctp) = CiphertextProto::parse_from_bytes(bytes) {
			Ciphertext::try_convert_from(&ctp, par)
		} else {
			Err(Error::SerializationError)
		}
	}

	type Error = Error;
}

impl Ciphertext {
	/// Generate the zero ciphertext.
	pub fn zero(par: &Arc<BfvParameters>) -> Self {
		Self {
			par: par.clone(),
			seed: None,
			c: Default::default(),
			level: 0,
		}
	}
}

/// Conversions from and to protobuf.
impl From<&Ciphertext> for CiphertextProto {
	fn from(ct: &Ciphertext) -> Self {
		let mut proto = CiphertextProto::new();
		for i in 0..ct.c.len() - 1 {
			proto.c.push(ct.c[i].to_bytes())
		}
		if let Some(seed) = ct.seed {
			proto.seed = seed.to_vec()
		} else {
			proto.c.push(ct.c[ct.c.len() - 1].to_bytes())
		}
		proto.level = ct.level as u32;
		proto
	}
}

impl TryConvertFrom<&CiphertextProto> for Ciphertext {
	fn try_convert_from(value: &CiphertextProto, par: &Arc<BfvParameters>) -> Result<Self> {
		if value.c.is_empty() || (value.c.len() == 1 && value.seed.is_empty()) {
			return Err(Error::DefaultError("Not enough polynomials".to_string()));
		}

		if value.level as usize > par.max_level() {
			return Err(Error::DefaultError("Invalid level".to_string()));
		}

		#[cfg(not(feature = "leveled_bfv"))]
		if value.level != 0 {
			return Err(Error::DefaultError(
				"Invalid level: did you enable the leveled_bfv feature?".to_string(),
			));
		}

		let ctx = par.ctx_at_level(value.level as usize)?;

		let mut seed = None;

		let mut c = Vec::with_capacity(value.c.len() + 1);
		for cip in &value.c {
			c.push(Poly::from_bytes(cip, ctx)?)
		}

		if !value.seed.is_empty() {
			let try_seed = <ChaCha8Rng as SeedableRng>::Seed::try_from(value.seed.clone());
			if try_seed.is_err() {
				return Err(Error::MathError(math::Error::InvalidSeedSize(
					value.seed.len(),
					<ChaCha8Rng as SeedableRng>::Seed::default().len(),
				)));
			}
			seed = try_seed.ok();
			let mut c1 = Poly::random_from_seed(ctx, Representation::Ntt, seed.unwrap());
			unsafe { c1.allow_variable_time_computations() }
			c.push(c1)
		}

		Ok(Ciphertext {
			par: par.clone(),
			seed,
			c,
			level: value.level as usize,
		})
	}
}

#[cfg(test)]
mod tests {
	use crate::bfv::{
		proto::bfv::Ciphertext as CiphertextProto, traits::TryConvertFrom, BfvParameters,
		Ciphertext, Encoding, Plaintext, SecretKey,
	};
	use fhers_traits::{
		DeserializeParametrized, FheDecrypter, FheEncoder, FheEncrypter, Serialize,
	};
	use std::{error::Error, sync::Arc};

	#[test]
	fn proto_conversion() -> Result<(), Box<dyn Error>> {
		for params in [
			Arc::new(BfvParameters::default(1, 8)),
			Arc::new(BfvParameters::default(2, 8)),
		] {
			let sk = SecretKey::random(&params);
			let v = params.plaintext.random_vec(params.degree());
			let pt = Plaintext::try_encode(&v as &[u64], Encoding::simd(), &params)?;
			let ct = sk.try_encrypt(&pt)?;
			let ct_proto = CiphertextProto::from(&ct);
			assert_eq!(ct, Ciphertext::try_convert_from(&ct_proto, &params)?);

			let ct = &ct * &ct;
			let ct_proto = CiphertextProto::from(&ct);
			assert_eq!(ct, Ciphertext::try_convert_from(&ct_proto, &params)?)
		}
		Ok(())
	}

	#[test]
	fn serialize() -> Result<(), Box<dyn Error>> {
		for params in [
			Arc::new(BfvParameters::default(1, 8)),
			Arc::new(BfvParameters::default(2, 8)),
		] {
			let sk = SecretKey::random(&params);
			let v = params.plaintext.random_vec(params.degree());
			let pt = Plaintext::try_encode(&v as &[u64], Encoding::simd(), &params)?;
			let ct: Ciphertext = sk.try_encrypt(&pt)?;
			let ct_bytes = ct.to_bytes();
			assert_eq!(ct, Ciphertext::from_bytes(&ct_bytes, &params)?);
		}
		Ok(())
	}

	#[test]
	fn new() -> Result<(), Box<dyn Error>> {
		for params in [
			Arc::new(BfvParameters::default(1, 8)),
			Arc::new(BfvParameters::default(2, 8)),
		] {
			let sk = SecretKey::random(&params);
			let v = params.plaintext.random_vec(params.degree());
			let pt = Plaintext::try_encode(&v as &[u64], Encoding::simd(), &params)?;
			let ct: Ciphertext = sk.try_encrypt(&pt)?;
			let mut ct3 = &ct * &ct;

			let c0 = ct3.get(0).unwrap();
			let c1 = ct3.get(1).unwrap();
			let c2 = ct3.get(2).unwrap();

			assert_eq!(
				ct3,
				Ciphertext::new(vec![c0.clone(), c1.clone(), c2.clone()], &params)?
			);
			assert_eq!(ct3.level, 0);

			ct3.mod_switch_to_last_level();

			let c0 = ct3.get(0).unwrap();
			let c1 = ct3.get(1).unwrap();
			let c2 = ct3.get(2).unwrap();
			assert_eq!(
				ct3,
				Ciphertext::new(vec![c0.clone(), c1.clone(), c2.clone()], &params)?
			);
			assert_eq!(ct3.level, params.max_level());
		}

		Ok(())
	}

	#[test]
	fn mod_switch_to_last_level() -> Result<(), Box<dyn Error>> {
		for params in [
			Arc::new(BfvParameters::default(1, 8)),
			Arc::new(BfvParameters::default(2, 8)),
		] {
			let sk = SecretKey::random(&params);
			let v = params.plaintext.random_vec(params.degree());
			let pt = Plaintext::try_encode(&v as &[u64], Encoding::simd(), &params)?;
			let mut ct: Ciphertext = sk.try_encrypt(&pt)?;

			assert_eq!(ct.level, 0);
			ct.mod_switch_to_last_level();
			assert_eq!(ct.level, params.max_level());

			let decrypted = sk.try_decrypt(&ct)?;
			assert_eq!(decrypted.value, pt.value);
		}

		Ok(())
	}
}
