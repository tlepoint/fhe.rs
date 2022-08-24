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
	/// Modulo switch the ciphertext to the last level.
	#[cfg(feature = "leveled_bfv")]
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

	/// Modulo switch the ciphertext to the next level.
	#[cfg(feature = "leveled_bfv")]
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
	pub fn new(c: Vec<Poly>, par: &Arc<BfvParameters>) -> Result<Self> {
		if c.len() < 2 {
			return Err(Error::DefaultError(
				"Too few polynomials provided".to_string(),
			));
		}
		let ctx = c[0].ctx();
		let level = par.level_of_ctx(ctx)?;
		if c.iter()
			.any(|ci| ci.representation() != &Representation::Ntt)
		{
			return Err(Error::DefaultError("Invalid representation".to_string()));
		}
		if c.iter().any(|ci| ci.ctx() != ctx) {
			return Err(Error::DefaultError("Invalid context".to_string()));
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
			Err(Error::DefaultError(
				"This serialization is incorrect".to_string(),
			))
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
				return Err(Error::DefaultError("Invalid seed".to_string()));
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
		encoding::EncodingEnum, proto::bfv::Ciphertext as CiphertextProto, traits::TryConvertFrom,
		BfvParameters, Ciphertext, Encoding, Plaintext, SecretKey,
	};
	use fhers_traits::{
		DeserializeParametrized, FheDecoder, FheDecrypter, FheEncoder, FheEncrypter, Serialize,
	};
	use std::{error::Error, sync::Arc};

	#[test]
	fn test_add() {
		let ntests = 100;
		for params in [
			Arc::new(BfvParameters::default(1)),
			Arc::new(BfvParameters::default(2)),
		] {
			for _ in 0..ntests {
				let a = params.plaintext.random_vec(params.degree());
				let b = params.plaintext.random_vec(params.degree());
				let mut c = a.clone();
				params.plaintext.add_vec(&mut c, &b);

				let sk = SecretKey::random(&params);

				for encoding in [Encoding::poly(), Encoding::simd()] {
					let pt_a =
						Plaintext::try_encode(&a as &[u64], encoding.clone(), &params).unwrap();
					let pt_b =
						Plaintext::try_encode(&b as &[u64], encoding.clone(), &params).unwrap();

					let mut ct_a = sk.try_encrypt(&pt_a).unwrap();
					let ct_b = sk.try_encrypt(&pt_b).unwrap();
					let ct_c = &ct_a + &ct_b;
					ct_a += &ct_b;

					let pt_c = sk.try_decrypt(&ct_c).unwrap();
					assert_eq!(Vec::<u64>::try_decode(&pt_c, encoding.clone()).unwrap(), c);
					let pt_c = sk.try_decrypt(&ct_a).unwrap();
					assert_eq!(Vec::<u64>::try_decode(&pt_c, encoding.clone()).unwrap(), c);
				}
			}
		}
	}

	#[test]
	fn test_sub() {
		for params in [
			Arc::new(BfvParameters::default(1)),
			Arc::new(BfvParameters::default(2)),
		] {
			let ntests = 100;
			for _ in 0..ntests {
				let a = params.plaintext.random_vec(params.degree());
				let b = params.plaintext.random_vec(params.degree());
				let mut c = a.clone();
				params.plaintext.sub_vec(&mut c, &b);

				let sk = SecretKey::random(&params);

				for encoding in [Encoding::poly(), Encoding::simd()] {
					let pt_a =
						Plaintext::try_encode(&a as &[u64], encoding.clone(), &params).unwrap();
					let pt_b =
						Plaintext::try_encode(&b as &[u64], encoding.clone(), &params).unwrap();

					let mut ct_a = sk.try_encrypt(&pt_a).unwrap();
					let ct_b = sk.try_encrypt(&pt_b).unwrap();
					let ct_c = &ct_a - &ct_b;
					ct_a -= &ct_b;

					let pt_c = sk.try_decrypt(&ct_c).unwrap();
					assert_eq!(Vec::<u64>::try_decode(&pt_c, encoding.clone()).unwrap(), c);
					let pt_c = sk.try_decrypt(&ct_a).unwrap();
					assert_eq!(Vec::<u64>::try_decode(&pt_c, encoding.clone()).unwrap(), c);
				}
			}
		}
	}

	#[test]
	fn test_neg() {
		for params in [
			Arc::new(BfvParameters::default(1)),
			Arc::new(BfvParameters::default(2)),
		] {
			let ntests = 100;
			for _ in 0..ntests {
				let a = params.plaintext.random_vec(params.degree());
				let mut c = a.clone();
				params.plaintext.neg_vec(&mut c);

				let sk = SecretKey::random(&params);
				for encoding in [Encoding::poly(), Encoding::simd()] {
					let pt_a =
						Plaintext::try_encode(&a as &[u64], encoding.clone(), &params).unwrap();

					let ct_a = sk.try_encrypt(&pt_a).unwrap();
					let ct_c = -&ct_a;

					let pt_c = sk.try_decrypt(&ct_c).unwrap();
					assert_eq!(Vec::<u64>::try_decode(&pt_c, encoding.clone()).unwrap(), c);
				}
			}
		}
	}

	#[test]
	fn test_scalar_mul() {
		for params in [
			Arc::new(BfvParameters::default(1)),
			Arc::new(BfvParameters::default(2)),
		] {
			let ntests = 100;
			for _ in 0..ntests {
				let a = params.plaintext.random_vec(params.degree());
				let b = params.plaintext.random_vec(params.degree());

				let sk = SecretKey::random(&params);
				for encoding in [Encoding::poly(), Encoding::simd()] {
					let mut c = vec![0u64; params.degree()];
					match encoding.encoding {
						EncodingEnum::Poly => {
							for i in 0..params.degree() {
								for j in 0..params.degree() {
									if i + j >= params.degree() {
										c[(i + j) % params.degree()] = params.plaintext.sub(
											c[(i + j) % params.degree()],
											params.plaintext.mul(a[i], b[j]),
										);
									} else {
										c[i + j] = params
											.plaintext
											.add(c[i + j], params.plaintext.mul(a[i], b[j]));
									}
								}
							}
						}
						EncodingEnum::Simd => {
							c = a.clone();
							params.plaintext.mul_vec(&mut c, &b);
						}
					}

					let pt_a =
						Plaintext::try_encode(&a as &[u64], encoding.clone(), &params).unwrap();
					let pt_b =
						Plaintext::try_encode(&b as &[u64], encoding.clone(), &params).unwrap();

					let mut ct_a = sk.try_encrypt(&pt_a).unwrap();
					let ct_c = &ct_a * &pt_b;
					ct_a *= &pt_b;

					let pt_c = sk.try_decrypt(&ct_c).unwrap();
					assert_eq!(Vec::<u64>::try_decode(&pt_c, encoding.clone()).unwrap(), c);
					let pt_c = sk.try_decrypt(&ct_a).unwrap();
					assert_eq!(Vec::<u64>::try_decode(&pt_c, encoding.clone()).unwrap(), c);
				}
			}
		}
	}

	#[test]
	fn test_mul() -> Result<(), Box<dyn Error>> {
		let par = Arc::new(BfvParameters::default(2));
		for _ in 0..50 {
			// We will encode `values` in an Simd format, and check that the product is
			// computed correctly.
			let values = par.plaintext.random_vec(par.degree());
			let mut expected = values.clone();
			par.plaintext.mul_vec(&mut expected, &values);

			let sk = SecretKey::random(&par);
			let pt = Plaintext::try_encode(&values as &[u64], Encoding::simd(), &par)?;

			let ct1 = sk.try_encrypt(&pt)?;
			let ct2 = sk.try_encrypt(&pt)?;
			let ct3 = &ct1 * &ct2;
			let ct4 = &ct3 * &ct3;

			println!("Noise: {}", unsafe { sk.measure_noise(&ct3)? });
			let pt = sk.try_decrypt(&ct3)?;
			assert_eq!(Vec::<u64>::try_decode(&pt, Encoding::simd())?, expected);

			let e = expected.clone();
			par.plaintext.mul_vec(&mut expected, &e);
			println!("Noise: {}", unsafe { sk.measure_noise(&ct4)? });
			let pt = sk.try_decrypt(&ct4)?;
			assert_eq!(Vec::<u64>::try_decode(&pt, Encoding::simd())?, expected);
		}
		Ok(())
	}

	#[test]
	fn test_proto_conversion() -> Result<(), Box<dyn Error>> {
		for params in [
			Arc::new(BfvParameters::default(1)),
			Arc::new(BfvParameters::default(2)),
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
	fn test_serialize() -> Result<(), Box<dyn Error>> {
		for params in [
			Arc::new(BfvParameters::default(1)),
			Arc::new(BfvParameters::default(2)),
		] {
			let sk = SecretKey::random(&params);
			let v = params.plaintext.random_vec(params.degree());
			let pt = Plaintext::try_encode(&v as &[u64], Encoding::simd(), &params)?;
			let ct = sk.try_encrypt(&pt)?;
			let ct_bytes = ct.to_bytes();
			assert_eq!(ct, Ciphertext::from_bytes(&ct_bytes, &params)?);
		}
		Ok(())
	}

	#[test]
	fn test_new() -> Result<(), Box<dyn Error>> {
		for params in [
			Arc::new(BfvParameters::default(1)),
			Arc::new(BfvParameters::default(2)),
		] {
			let sk = SecretKey::random(&params);
			let v = params.plaintext.random_vec(params.degree());
			let pt = Plaintext::try_encode(&v as &[u64], Encoding::simd(), &params)?;
			let ct = sk.try_encrypt(&pt)?;
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
	fn test_mod_switch_to_last_level() -> Result<(), Box<dyn Error>> {
		for params in [
			Arc::new(BfvParameters::default(1)),
			Arc::new(BfvParameters::default(2)),
		] {
			let sk = SecretKey::random(&params);
			let v = params.plaintext.random_vec(params.degree());
			let pt = Plaintext::try_encode(&v as &[u64], Encoding::simd(), &params)?;
			let mut ct = sk.try_encrypt(&pt)?;

			assert_eq!(ct.level, 0);
			ct.mod_switch_to_last_level();
			assert_eq!(ct.level, params.max_level());

			let decrypted = sk.try_decrypt(&ct)?;
			assert_eq!(decrypted.value, pt.value);
		}

		Ok(())
	}
}
