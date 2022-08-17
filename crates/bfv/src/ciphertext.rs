//! Ciphertext type in the BFV encryption scheme.

use crate::{
	parameters::{BfvParameters, MultiplicationParameters},
	traits::{DeserializeWithParams, Serialize, TryConvertFrom},
	Error, EvaluationKey, Plaintext, Result,
};
use fhers_protos::protos::{bfv::Ciphertext as CiphertextProto, rq::Rq};
use itertools::{izip, Itertools};
use math::rq::{
	dot_product as poly_dot_product, traits::TryConvertFrom as PolyTryConvertFrom, Poly,
	Representation,
};
use num_bigint::BigUint;
use protobuf::Message;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use std::{
	cmp::min,
	ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
	sync::Arc,
};

/// A ciphertext encrypting a plaintext.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Ciphertext {
	/// The parameters of the underlying BFV encryption scheme.
	pub(crate) par: Arc<BfvParameters>,

	/// The seed that generated the polynomial c1 in a fresh ciphertext.
	pub(crate) seed: Option<<ChaCha8Rng as SeedableRng>::Seed>,

	/// The ciphertext elements.
	pub(crate) c: Vec<Poly>,
}

impl Ciphertext {
	/// Generate the zero ciphertext.
	pub fn zero(par: &Arc<BfvParameters>) -> Self {
		Self {
			par: par.clone(),
			seed: None,
			c: Default::default(),
		}
	}
}

impl Add<&Ciphertext> for &Ciphertext {
	type Output = Ciphertext;

	fn add(self, rhs: &Ciphertext) -> Ciphertext {
		let mut self_clone = self.clone();
		self_clone += rhs;
		self_clone
	}
}

impl AddAssign<&Ciphertext> for Ciphertext {
	fn add_assign(&mut self, rhs: &Ciphertext) {
		assert_eq!(self.par, rhs.par);

		if self.c.is_empty() {
			*self = rhs.clone()
		} else {
			assert_eq!(self.c.len(), rhs.c.len());
			izip!(&mut self.c, &rhs.c).for_each(|(c1i, c2i)| *c1i += c2i);
			self.seed = None
		}
	}
}

impl AddAssign<Ciphertext> for Ciphertext {
	fn add_assign(&mut self, rhs: Ciphertext) {
		assert_eq!(self.par, rhs.par);

		if self.c.is_empty() {
			*self = rhs
		} else {
			assert_eq!(self.c.len(), rhs.c.len());
			izip!(&mut self.c, rhs.c).for_each(|(c1i, c2i)| *c1i += c2i);
			self.seed = None
		}
	}
}

impl Sub<&Ciphertext> for &Ciphertext {
	type Output = Ciphertext;

	fn sub(self, rhs: &Ciphertext) -> Ciphertext {
		let mut self_clone = self.clone();
		self_clone -= rhs;
		self_clone
	}
}

impl SubAssign<&Ciphertext> for Ciphertext {
	fn sub_assign(&mut self, rhs: &Ciphertext) {
		assert_eq!(self.par, rhs.par);

		if self.c.is_empty() {
			*self = -rhs
		} else {
			assert_eq!(self.c.len(), rhs.c.len());
			izip!(&mut self.c, &rhs.c).for_each(|(c1i, c2i)| *c1i -= c2i);
			self.seed = None
		}
	}
}

impl Neg for &Ciphertext {
	type Output = Ciphertext;

	fn neg(self) -> Ciphertext {
		let c = self.c.iter().map(|c1i| -c1i).collect_vec();
		Ciphertext {
			par: self.par.clone(),
			seed: None,
			c,
		}
	}
}

impl MulAssign<&Plaintext> for Ciphertext {
	fn mul_assign(&mut self, rhs: &Plaintext) {
		assert_eq!(self.par, rhs.par);
		self.c.iter_mut().for_each(|ci| *ci *= &rhs.poly_ntt);
		self.seed = None
	}
}

impl Mul<&Plaintext> for &Ciphertext {
	type Output = Ciphertext;

	fn mul(self, rhs: &Plaintext) -> Self::Output {
		let mut self_clone = self.clone();
		self_clone *= rhs;
		self_clone
	}
}

#[allow(dead_code)]
fn print_poly(s: &str, p: &Poly) {
	println!("{} = {:?}", s, Vec::<BigUint>::from(p))
}

impl Mul<&Ciphertext> for &Ciphertext {
	type Output = Ciphertext;

	fn mul(self, rhs: &Ciphertext) -> Self::Output {
		assert_eq!(self.par, rhs.par);

		if self.c.is_empty() {
			return self.clone();
		}

		let mp = &self.par.mul_1_params;

		// Scale all ciphertexts
		// let mut now = std::time::SystemTime::now();
		let self_c = self
			.c
			.iter()
			.map(|ci| mp.extender_self.scale(ci, false).unwrap())
			.collect_vec();
		let other_c = rhs
			.c
			.iter()
			.map(|ci| mp.extender_self.scale(ci, false).unwrap())
			.collect_vec();
		// println!("Extend: {:?}", now.elapsed().unwrap());

		// Multiply
		// now = std::time::SystemTime::now();
		let mut c = vec![Poly::zero(&mp.to, Representation::Ntt); self_c.len() + other_c.len() - 1];
		for i in 0..self_c.len() {
			for j in 0..other_c.len() {
				c[i + j] += &(&self_c[i] * &other_c[j])
			}
		}
		// println!("Multiply: {:?}", now.elapsed().unwrap());

		// Scale
		// now = std::time::SystemTime::now();
		let c = c
			.iter_mut()
			.map(|ci| {
				ci.change_representation(Representation::PowerBasis);
				let mut ci = mp.down_scaler.scale(ci, false).unwrap();
				ci.change_representation(Representation::Ntt);
				ci
			})
			.collect_vec();
		// println!("Scale: {:?}", now.elapsed().unwrap());

		Ciphertext {
			par: self.par.clone(),
			seed: None,
			c,
		}
	}
}

/// Compute the dot product between an iterator of [`Ciphertext`] and an iterator of [`Plaintext`].
/// Returns an error if the iterator counts are 0, if the parameters don't match, or if the ciphertexts have different
/// number of parts.
pub fn dot_product_scalar<'a, I, J>(ct: I, pt: J) -> Result<Ciphertext>
where
	I: Iterator<Item = &'a Ciphertext> + Clone,
	J: Iterator<Item = &'a Plaintext> + Clone,
{
	let count = min(ct.clone().count(), pt.clone().count());
	if count == 0 {
		return Err(Error::DefaultError(
			"At least one iterator is empty".to_string(),
		));
	}

	let ct_first = ct.clone().next().unwrap();
	if izip!(ct.clone(), pt.clone()).any(|(cti, pti)| {
		cti.par != ct_first.par || pti.par != ct_first.par || cti.c.len() != ct_first.c.len()
	}) {
		return Err(Error::DefaultError("Mismatched parameters".to_string()));
	}
	if ct.clone().any(|cti| cti.c.len() != ct_first.c.len()) {
		return Err(Error::DefaultError(
			"Mismatched number of parts in the ciphertexts".to_string(),
		));
	}

	let c = (0..ct_first.c.len())
		.map(|i| {
			poly_dot_product(
				ct.clone().map(|cti| unsafe { cti.c.get_unchecked(i) }),
				pt.clone().map(|pti| &pti.poly_ntt),
			)
			.unwrap()
		})
		.collect_vec();

	Ok(Ciphertext {
		par: ct_first.par.clone(),
		seed: None,
		c,
	})
}

/// Multiply two ciphertext and relinearize.
fn mul_internal(
	ct0: &Ciphertext,
	ct1: &Ciphertext,
	ek: &EvaluationKey,
	mp: &MultiplicationParameters,
) -> Result<Ciphertext> {
	if !ek.supports_relinearization() {
		return Err(Error::DefaultError(
			"The evaluation key does not support relinearization".to_string(),
		));
	}
	if ct0.par != ct1.par {
		return Err(Error::DefaultError("Incompatible parameters".to_string()));
	}
	if ct0.par.ciphertext_moduli.len() == 1 {
		return Err(Error::DefaultError(
			"Parameters do not allow for multiplication".to_string(),
		));
	}
	if ct0.c.len() != 2 || ct1.c.len() != 2 {
		return Err(Error::DefaultError(
			"Multiplication can only be performed on ciphertexts of size 2".to_string(),
		));
	}

	// Extend
	// let mut now = std::time::SystemTime::now();
	let c00 = mp.extender_self.scale(&ct0.c[0], false)?;
	let c01 = mp.extender_self.scale(&ct0.c[1], false)?;
	let c10 = mp.extender_other.scale(&ct1.c[0], false)?;
	let c11 = mp.extender_other.scale(&ct1.c[1], false)?;
	// println!("Extend: {:?}", now.elapsed().unwrap());

	// Multiply
	// now = std::time::SystemTime::now();
	let mut c0 = &c00 * &c10;
	let mut c1 = &c00 * &c11;
	c1 += &(&c01 * &c10);
	let mut c2 = &c01 * &c11;
	c0.change_representation(Representation::PowerBasis);
	c1.change_representation(Representation::PowerBasis);
	c2.change_representation(Representation::PowerBasis);
	// println!("Multiply: {:?}", now.elapsed().unwrap());

	// Scale
	// TODO: This should be faster??
	// now = std::time::SystemTime::now();
	let mut c0 = mp.down_scaler.scale(&c0, false)?;
	let mut c1 = mp.down_scaler.scale(&c1, false)?;
	let c2 = mp.down_scaler.scale(&c2, false)?;
	// println!("Scale: {:?}", now.elapsed().unwrap());

	// Relinearize
	// now = std::time::SystemTime::now();
	c0.change_representation(Representation::Ntt);
	c1.change_representation(Representation::Ntt);
	ek.relinearizes_with_poly(&c2, &mut c0, &mut c1)?;
	// println!("Relinearize: {:?}", now.elapsed().unwrap());

	Ok(Ciphertext {
		par: ct0.par.clone(),
		seed: None,
		c: vec![c0, c1],
	})
}

/// Multiply two ciphertext and relinearize.
pub fn mul(ct0: &Ciphertext, ct1: &Ciphertext, ek: &EvaluationKey) -> Result<Ciphertext> {
	mul_internal(ct0, ct1, ek, &ct0.par.mul_1_params)
}

/// Multiply two ciphertext and relinearize.
pub fn mul2(ct0: &Ciphertext, ct1: &Ciphertext, ek: &EvaluationKey) -> Result<Ciphertext> {
	mul_internal(ct0, ct1, ek, &ct0.par.mul_2_params)
}

/// Conversions from and to protobuf.
impl From<&Ciphertext> for CiphertextProto {
	fn from(ct: &Ciphertext) -> Self {
		let mut proto = CiphertextProto::new();
		for i in 0..ct.c.len() - 1 {
			proto.c.push(Rq::from(&ct.c[i]))
		}
		if let Some(seed) = ct.seed {
			proto.seed = seed.to_vec()
		} else {
			proto.c.push(Rq::from(&ct.c[ct.c.len() - 1]))
		}
		proto
	}
}

impl TryConvertFrom<&CiphertextProto> for Ciphertext {
	fn try_convert_from(value: &CiphertextProto, par: &Arc<BfvParameters>) -> Result<Self> {
		if value.c.is_empty() || (value.c.len() == 1 && value.seed.is_empty()) {
			return Err(Error::DefaultError("Not enough polynomials".to_string()));
		}

		let mut seed = None;

		let mut c = Vec::with_capacity(value.c.len() + 1);
		for cip in &value.c {
			c.push(Poly::try_convert_from(cip, &par.ctx, true, None)?)
		}

		if !value.seed.is_empty() {
			let try_seed = <ChaCha8Rng as SeedableRng>::Seed::try_from(value.seed.clone());
			if try_seed.is_err() {
				return Err(Error::DefaultError("Invalid seed".to_string()));
			}
			seed = try_seed.ok();
			let mut c1 = Poly::random_from_seed(&par.ctx, Representation::Ntt, seed.unwrap());
			unsafe { c1.allow_variable_time_computations() }
			c.push(c1)
		}

		Ok(Ciphertext {
			par: par.clone(),
			seed,
			c,
		})
	}
}

impl Serialize for Ciphertext {
	fn serialize(&self) -> Vec<u8> {
		CiphertextProto::from(self).write_to_bytes().unwrap()
	}
}

impl DeserializeWithParams for Ciphertext {
	fn try_deserialize(bytes: &[u8], par: &Arc<BfvParameters>) -> Result<Self> {
		if let Ok(ctp) = CiphertextProto::parse_from_bytes(bytes) {
			Ciphertext::try_convert_from(&ctp, par)
		} else {
			Err(Error::DefaultError(
				"This serialization is incorrect".to_string(),
			))
		}
	}
}

#[cfg(test)]
mod tests {
	use super::{dot_product_scalar, mul, mul2};
	use crate::{
		traits::{Decoder, Decryptor, Encoder, Encryptor, TryConvertFrom},
		BfvParameters, Ciphertext, Encoding, EvaluationKeyBuilder, Plaintext, SecretKey,
	};
	use fhers_protos::protos::bfv::Ciphertext as CiphertextProto;
	use itertools::{izip, Itertools};
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

				let mut sk = SecretKey::random(&params);

				for encoding in [Encoding::Poly, Encoding::Simd] {
					let pt_a =
						Plaintext::try_encode(&a as &[u64], encoding.clone(), &params).unwrap();
					let pt_b =
						Plaintext::try_encode(&b as &[u64], encoding.clone(), &params).unwrap();

					let mut ct_a = sk.encrypt(&pt_a).unwrap();
					let ct_b = sk.encrypt(&pt_b).unwrap();
					let ct_c = &ct_a + &ct_b;
					ct_a += &ct_b;

					let pt_c = sk.decrypt(&ct_c).unwrap();
					assert_eq!(Vec::<u64>::try_decode(&pt_c, encoding.clone()).unwrap(), c);
					let pt_c = sk.decrypt(&ct_a).unwrap();
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

				let mut sk = SecretKey::random(&params);

				for encoding in [Encoding::Poly, Encoding::Simd] {
					let pt_a =
						Plaintext::try_encode(&a as &[u64], encoding.clone(), &params).unwrap();
					let pt_b =
						Plaintext::try_encode(&b as &[u64], encoding.clone(), &params).unwrap();

					let mut ct_a = sk.encrypt(&pt_a).unwrap();
					let ct_b = sk.encrypt(&pt_b).unwrap();
					let ct_c = &ct_a - &ct_b;
					ct_a -= &ct_b;

					let pt_c = sk.decrypt(&ct_c).unwrap();
					assert_eq!(Vec::<u64>::try_decode(&pt_c, encoding.clone()).unwrap(), c);
					let pt_c = sk.decrypt(&ct_a).unwrap();
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

				let mut sk = SecretKey::random(&params);
				for encoding in [Encoding::Poly, Encoding::Simd] {
					let pt_a =
						Plaintext::try_encode(&a as &[u64], encoding.clone(), &params).unwrap();

					let ct_a = sk.encrypt(&pt_a).unwrap();
					let ct_c = -&ct_a;

					let pt_c = sk.decrypt(&ct_c).unwrap();
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

				let mut sk = SecretKey::random(&params);
				for encoding in [Encoding::Poly, Encoding::Simd] {
					let mut c = vec![0u64; params.degree()];
					match encoding {
						Encoding::Poly => {
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
						Encoding::Simd => {
							c = a.clone();
							params.plaintext.mul_vec(&mut c, &b);
						}
					}

					let pt_a =
						Plaintext::try_encode(&a as &[u64], encoding.clone(), &params).unwrap();
					let pt_b =
						Plaintext::try_encode(&b as &[u64], encoding.clone(), &params).unwrap();

					let mut ct_a = sk.encrypt(&pt_a).unwrap();
					let ct_c = &ct_a * &pt_b;
					ct_a *= &pt_b;

					let pt_c = sk.decrypt(&ct_c).unwrap();
					assert_eq!(Vec::<u64>::try_decode(&pt_c, encoding.clone()).unwrap(), c);
					let pt_c = sk.decrypt(&ct_a).unwrap();
					assert_eq!(Vec::<u64>::try_decode(&pt_c, encoding.clone()).unwrap(), c);
				}
			}
		}
	}

	#[test]
	fn test_mul() -> Result<(), Box<dyn Error>> {
		let par = Arc::new(BfvParameters::default(2));
		for _ in 0..50 {
			// We will encode `values` in an Simd format, and check that the product is computed correctly.
			let values = par.plaintext.random_vec(par.degree());
			let mut expected = values.clone();
			par.plaintext.mul_vec(&mut expected, &values);

			let mut sk = SecretKey::random(&par);
			let pt = Plaintext::try_encode(&values as &[u64], Encoding::Simd, &par)?;

			let ct1 = sk.encrypt(&pt)?;
			let ct2 = sk.encrypt(&pt)?;
			let ct3 = &ct1 * &ct2;
			let ct4 = &ct3 * &ct3;

			println!("Noise: {}", unsafe { sk.measure_noise(&ct3)? });
			let pt = sk.decrypt(&ct3)?;
			assert_eq!(Vec::<u64>::try_decode(&pt, Encoding::Simd)?, expected);

			let e = expected.clone();
			par.plaintext.mul_vec(&mut expected, &e);
			println!("Noise: {}", unsafe { sk.measure_noise(&ct4)? });
			let pt = sk.decrypt(&ct4)?;
			assert_eq!(Vec::<u64>::try_decode(&pt, Encoding::Simd)?, expected);
		}
		Ok(())
	}

	#[test]
	fn test_mul_3() -> Result<(), Box<dyn Error>> {
		let par = Arc::new(BfvParameters::default(2));
		for _ in 0..50 {
			// We will encode `values` in an Simd format, and check that the product is computed correctly.
			let values = par.plaintext.random_vec(par.degree());
			let mut expected = values.clone();
			par.plaintext.mul_vec(&mut expected, &values);

			let mut sk = SecretKey::random(&par);
			let ek = EvaluationKeyBuilder::new(&sk)
				.enable_relinearization()?
				.build()?;
			let pt = Plaintext::try_encode(&values as &[u64], Encoding::Simd, &par)?;

			let ct1 = sk.encrypt(&pt)?;
			let ct2 = sk.encrypt(&pt)?;
			let ct3 = mul(&ct1, &ct2, &ek)?;

			println!("Noise: {}", unsafe { sk.measure_noise(&ct3)? });
			let pt = sk.decrypt(&ct3)?;
			assert_eq!(Vec::<u64>::try_decode(&pt, Encoding::Simd)?, expected);
		}
		Ok(())
	}

	#[test]
	fn test_mul2() -> Result<(), Box<dyn Error>> {
		let ntests = 100;
		let par = Arc::new(BfvParameters::default(2));
		for _ in 0..ntests {
			// We will encode `values` in an Simd format, and check that the product is computed correctly.
			let values = par.plaintext.random_vec(par.degree());
			let mut expected = values.clone();
			par.plaintext.mul_vec(&mut expected, &values);

			let mut sk = SecretKey::random(&par);
			let ek = EvaluationKeyBuilder::new(&sk)
				.enable_relinearization()?
				.build()?;
			let pt = Plaintext::try_encode(&values as &[u64], Encoding::Simd, &par)?;

			let ct1 = sk.encrypt(&pt)?;
			let ct2 = sk.encrypt(&pt)?;
			let ct3 = mul2(&ct1, &ct2, &ek)?;

			println!("Noise: {}", unsafe { sk.measure_noise(&ct3)? });
			let pt = sk.decrypt(&ct3)?;
			assert_eq!(Vec::<u64>::try_decode(&pt, Encoding::Simd)?, expected);
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
			let pt = Plaintext::try_encode(&v as &[u64], Encoding::Simd, &params)?;
			let ct = sk.encrypt(&pt)?;
			let ct_proto = CiphertextProto::from(&ct);
			assert_eq!(ct, Ciphertext::try_convert_from(&ct_proto, &params)?);

			let ct = &ct * &ct;
			let ct_proto = CiphertextProto::from(&ct);
			assert_eq!(ct, Ciphertext::try_convert_from(&ct_proto, &params)?)
		}
		Ok(())
	}

	#[test]
	fn test_dot_product_scalar() -> Result<(), Box<dyn Error>> {
		for params in [
			Arc::new(BfvParameters::default(1)),
			Arc::new(BfvParameters::default(2)),
		] {
			let sk = SecretKey::random(&params);
			for size in 1..128 {
				let ct = (0..size)
					.map(|_| {
						let v = params.plaintext.random_vec(params.degree());
						let pt =
							Plaintext::try_encode(&v as &[u64], Encoding::Simd, &params).unwrap();
						sk.encrypt(&pt).unwrap()
					})
					.collect_vec();
				let pt = (0..size)
					.map(|_| {
						let v = params.plaintext.random_vec(params.degree());
						Plaintext::try_encode(&v as &[u64], Encoding::Simd, &params).unwrap()
					})
					.collect_vec();

				let r = dot_product_scalar(ct.iter(), pt.iter())?;

				let mut expected = Ciphertext::zero(&params);
				izip!(&ct, &pt).for_each(|(cti, pti)| expected += cti * pti);
				assert_eq!(r, expected);
			}
		}
		Ok(())
	}
}
