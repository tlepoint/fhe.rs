//! Ciphertext type in the BFV encryption scheme.

use crate::{
	parameters::{BfvParameters, MultiplicationParameters},
	traits::TryConvertFrom,
	Plaintext, RelinearizationKey,
};
use fhers_protos::protos::{bfv::Ciphertext as CiphertextProto, rq::Rq};
use math::rq::{traits::TryConvertFrom as PolyTryConvertFrom, Poly, Representation};
use num_bigint::BigUint;
use protobuf::MessageField;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use std::{
	ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
	rc::Rc,
};

/// A ciphertext encrypting a plaintext.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Ciphertext {
	/// The parameters of the underlying BFV encryption scheme.
	pub(crate) par: Rc<BfvParameters>,

	/// The seed that generated the polynomial c1 in a fresh ciphertext.
	pub(crate) seed: Option<<ChaCha8Rng as SeedableRng>::Seed>,

	/// The ciphertext element c0.
	pub(crate) c0: Poly,

	/// The ciphertext element c1.
	pub(crate) c1: Poly,
}

impl Add<&Ciphertext> for &Ciphertext {
	type Output = Ciphertext;

	fn add(self, rhs: &Ciphertext) -> Ciphertext {
		debug_assert_eq!(self.par, rhs.par);

		Ciphertext {
			par: self.par.clone(),
			seed: None,
			c0: &self.c0 + &rhs.c0,
			c1: &self.c1 + &rhs.c1,
		}
	}
}

impl AddAssign<&Ciphertext> for Ciphertext {
	fn add_assign(&mut self, rhs: &Ciphertext) {
		debug_assert_eq!(self.par, rhs.par);

		self.c0 += &rhs.c0;
		self.c1 += &rhs.c1;
		self.seed = None
	}
}

impl Sub<&Ciphertext> for &Ciphertext {
	type Output = Ciphertext;

	fn sub(self, rhs: &Ciphertext) -> Ciphertext {
		assert_eq!(self.par, rhs.par);

		Ciphertext {
			par: self.par.clone(),
			seed: None,
			c0: &self.c0 - &rhs.c0,
			c1: &self.c1 - &rhs.c1,
		}
	}
}

impl SubAssign<&Ciphertext> for Ciphertext {
	fn sub_assign(&mut self, rhs: &Ciphertext) {
		debug_assert_eq!(self.par, rhs.par);

		self.c0 -= &rhs.c0;
		self.c1 -= &rhs.c1;
		self.seed = None
	}
}

impl Neg for &Ciphertext {
	type Output = Ciphertext;

	fn neg(self) -> Ciphertext {
		Ciphertext {
			par: self.par.clone(),
			seed: None,
			c0: -&self.c0,
			c1: -&self.c1,
		}
	}
}

impl MulAssign<&Plaintext> for Ciphertext {
	fn mul_assign(&mut self, rhs: &Plaintext) {
		assert_eq!(self.par, rhs.par);

		self.c0 *= &rhs.poly_ntt;
		self.c1 *= &rhs.poly_ntt;
		self.seed = None
	}
}

impl Mul<&Plaintext> for &Ciphertext {
	type Output = Ciphertext;

	fn mul(self, rhs: &Plaintext) -> Self::Output {
		Ciphertext {
			par: self.par.clone(),
			seed: None,
			c0: &self.c0 * &rhs.poly_ntt,
			c1: &self.c1 * &rhs.poly_ntt,
		}
	}
}

#[allow(dead_code)]
fn print_poly(s: &str, p: &Poly) {
	println!("{} = {:?}", s, Vec::<BigUint>::from(p))
}

/// Multiply two ciphertext and relinearize.
fn mul_internal(
	ct0: &Ciphertext,
	ct1: &Ciphertext,
	rk: &RelinearizationKey,
	mp: &MultiplicationParameters,
) -> Result<Ciphertext, String> {
	if ct0.par != ct1.par {
		return Err("Incompatible parameters".to_string());
	}
	if ct0.par.ciphertext_moduli.len() == 1 {
		return Err("Parameters do not allow for multiplication".to_string());
	}
	// Extend
	let mut now = std::time::SystemTime::now();
	let c00 = mp.extender_self.scale(&ct0.c0, false)?;
	let c01 = mp.extender_self.scale(&ct0.c1, false)?;
	let c10 = mp.extender_other.scale(&ct1.c0, false)?;
	let c11 = mp.extender_other.scale(&ct1.c1, false)?;
	println!("Extend: {:?}", now.elapsed().unwrap());

	// Multiply
	now = std::time::SystemTime::now();
	let mut c0 = &c00 * &c10;
	let mut c1 = &c00 * &c11;
	c1 += &(&c01 * &c10);
	let mut c2 = &c01 * &c11;
	c0.change_representation(Representation::PowerBasis);
	c1.change_representation(Representation::PowerBasis);
	c2.change_representation(Representation::PowerBasis);
	println!("Multiply: {:?}", now.elapsed().unwrap());

	// Scale
	// TODO: This should be faster??
	now = std::time::SystemTime::now();
	let mut c0 = mp.down_scaler.scale(&c0, false)?;
	let mut c1 = mp.down_scaler.scale(&c1, false)?;
	let c2 = mp.down_scaler.scale(&c2, false)?;
	println!("Scale: {:?}", now.elapsed().unwrap());

	// Relinearize
	now = std::time::SystemTime::now();
	c0.change_representation(Representation::Ntt);
	c1.change_representation(Representation::Ntt);
	rk.relinearize(&mut c0, &mut c1, &c2)?;
	println!("Relinearize: {:?}", now.elapsed().unwrap());

	Ok(Ciphertext {
		par: ct0.par.clone(),
		seed: None,
		c0,
		c1,
	})
}

/// Multiply two ciphertext and relinearize.
pub fn mul(
	ct0: &Ciphertext,
	ct1: &Ciphertext,
	rk: &RelinearizationKey,
) -> Result<Ciphertext, String> {
	mul_internal(ct0, ct1, rk, &ct0.par.mul_1_params)
}

/// Multiply two ciphertext and relinearize.
pub fn mul2(
	ct0: &Ciphertext,
	ct1: &Ciphertext,
	rk: &RelinearizationKey,
) -> Result<Ciphertext, String> {
	mul_internal(ct0, ct1, rk, &ct0.par.mul_2_params)
}

/// Conversions from and to protobuf.
impl From<&Ciphertext> for CiphertextProto {
	fn from(ct: &Ciphertext) -> Self {
		let mut proto = CiphertextProto::new();
		if let Some(seed) = ct.seed {
			proto.set_c1_seed(seed.to_vec())
		} else {
			proto.set_c1_poly(Rq::from(&ct.c1))
		}
		proto.c0 = MessageField::some(Rq::from(&ct.c0));
		proto
	}
}

impl TryConvertFrom<&CiphertextProto> for Ciphertext {
	type Error = String;

	fn try_convert_from(
		value: &CiphertextProto,
		par: &Rc<BfvParameters>,
	) -> Result<Self, Self::Error> {
		if value.c0.is_none() {
			return Err("Missing c0".to_string());
		}

		let mut c0 = Poly::try_convert_from(value.c0.as_ref().unwrap(), &par.ctx, None)?;
		unsafe { c0.allow_variable_time_computations() }

		let mut seed = None;
		let mut c1;
		if value.has_c1_poly() {
			c1 = Poly::try_convert_from(value.c1_poly(), &par.ctx, None)?;
		} else if value.has_c1_seed() {
			let try_seed = <ChaCha8Rng as SeedableRng>::Seed::try_from(value.c1_seed());
			if try_seed.is_err() {
				return Err("Invalid seed".to_string());
			}
			seed = Some(try_seed.unwrap());
			c1 = Poly::random_from_seed(&par.ctx, Representation::Ntt, seed.unwrap());
			unsafe { c1.allow_variable_time_computations() }
		} else {
			return Err("c1 or the seed are missing".to_string());
		}

		Ok(Ciphertext {
			par: par.clone(),
			seed,
			c0,
			c1,
		})
	}
}

#[cfg(test)]
mod tests {
	use super::{mul, mul2};
	use crate::{
		traits::{Decoder, Decryptor, Encoder, Encryptor, TryConvertFrom},
		BfvParameters, BfvParametersBuilder, Ciphertext, Encoding, Plaintext, RelinearizationKey,
		SecretKey,
	};
	use fhers_protos::protos::bfv::Ciphertext as CiphertextProto;
	use proptest::collection::vec as prop_vec;
	use proptest::prelude::any;
	use std::rc::Rc;

	proptest! {
		#[test]
		fn test_add(mut a in prop_vec(any::<u64>(), 8), mut b in prop_vec(any::<u64>(), 8)) {
			for params in [Rc::new(BfvParameters::default(1)),
						   Rc::new(BfvParameters::default(2))] {
				params.plaintext.reduce_vec(&mut a);
				params.plaintext.reduce_vec(&mut b);
				let mut c = a.clone();
				params.plaintext.add_vec(&mut c, &b);

				let sk = SecretKey::random(&params);

				for encoding in [Encoding::Poly, Encoding::Simd] {
					let pt_a = Plaintext::try_encode(&a as &[u64], encoding.clone(), &params).unwrap();
					let pt_b = Plaintext::try_encode(&b as &[u64], encoding.clone(), &params).unwrap();

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

		#[test]
		fn test_sub(mut a in prop_vec(any::<u64>(), 8), mut b in prop_vec(any::<u64>(), 8)) {
			for params in [Rc::new(BfvParameters::default(1)),
						   Rc::new(BfvParameters::default(2))] {
				params.plaintext.reduce_vec(&mut a);
				params.plaintext.reduce_vec(&mut b);
				let mut c = a.clone();
				params.plaintext.sub_vec(&mut c, &b);

				let sk = SecretKey::random(&params);

				for encoding in [Encoding::Poly, Encoding::Simd] {
					let pt_a = Plaintext::try_encode(&a as &[u64], encoding.clone(), &params).unwrap();
					let pt_b = Plaintext::try_encode(&b as &[u64], encoding.clone(), &params).unwrap();

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

		#[test]
		fn test_neg(mut a in prop_vec(any::<u64>(), 8)) {
			for params in [Rc::new(BfvParameters::default(1)),
						   Rc::new(BfvParameters::default(2))] {
				params.plaintext.reduce_vec(&mut a);
				let mut c = a.clone();
				params.plaintext.neg_vec(&mut c);

				let sk = SecretKey::random(&params);
				for encoding in [Encoding::Poly, Encoding::Simd] {
					let pt_a = Plaintext::try_encode(&a as &[u64], encoding.clone(), &params).unwrap();

					let ct_a = sk.encrypt(&pt_a).unwrap();
					let ct_c = -&ct_a;

					let pt_c = sk.decrypt(&ct_c).unwrap();
					assert_eq!(Vec::<u64>::try_decode(&pt_c, encoding.clone()).unwrap(), c);
				}
			}
		}

		#[test]
		fn test_scalar_mul(mut a in prop_vec(any::<u64>(), 8), mut b in prop_vec(any::<u64>(), 8)) {
			for params in [Rc::new(BfvParameters::default(1)),
						   Rc::new(BfvParameters::default(2))] {
				params.plaintext.reduce_vec(&mut a);
				params.plaintext.reduce_vec(&mut b);

				let sk = SecretKey::random(&params);
				for encoding in [Encoding::Poly, Encoding::Simd] {
					let mut c = vec![0u64; 8];
					match encoding {
						Encoding::Poly => {
							for i in 0..8 {
								for j in 0..8 {
									if i + j >= 8 {
										c[(i+j) % 8] = params.plaintext.sub(c[(i+j) % 8], params.plaintext.mul(a[i], b[j]));
									} else {
										c[i+j] = params.plaintext.add(c[i+j], params.plaintext.mul(a[i], b[j]));
									}
								}
							}
						}
						Encoding::Simd => {
							c = a.clone();
							params.plaintext.mul_vec(&mut c, &b);
						}
					}

					let pt_a = Plaintext::try_encode(&a as &[u64], encoding.clone(), &params).unwrap();
					let pt_b = Plaintext::try_encode(&b as &[u64], encoding.clone(), &params).unwrap();

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
	fn test_mul() -> Result<(), String> {
		let par = Rc::new(BfvParameters::default(2));
		for _ in 0..50 {
			// We will encode `values` in an Simd format, and check that the product is computed correctly.
			let values = par.plaintext.random_vec(par.polynomial_degree);
			let mut expected = values.clone();
			par.plaintext.mul_vec(&mut expected, &values);

			let sk = SecretKey::random(&par);
			let rk = RelinearizationKey::new(&sk)?;
			let pt = Plaintext::try_encode(&values as &[u64], Encoding::Simd, &par)?;

			let ct1 = sk.encrypt(&pt)?;
			let ct2 = sk.encrypt(&pt)?;
			let ct3 = mul(&ct1, &ct2, &rk)?;

			println!("Noise: {}", unsafe { sk.measure_noise(&ct3)? });
			let pt = sk.decrypt(&ct3)?;
			assert_eq!(Vec::<u64>::try_decode(&pt, Encoding::Simd)?, expected);
		}
		Ok(())
	}

	#[test]
	fn test_mul2() -> Result<(), String> {
		let par = Rc::new(BfvParameters::default(2));
		for _ in 0..100 {
			// We will encode `values` in an Simd format, and check that the product is computed correctly.
			let values = par.plaintext.random_vec(par.polynomial_degree);
			let mut expected = values.clone();
			par.plaintext.mul_vec(&mut expected, &values);

			let sk = SecretKey::random(&par);
			let rk = RelinearizationKey::new(&sk)?;
			let pt = Plaintext::try_encode(&values as &[u64], Encoding::Simd, &par)?;

			let ct1 = sk.encrypt(&pt)?;
			let ct2 = sk.encrypt(&pt)?;
			let ct3 = mul2(&ct1, &ct2, &rk)?;

			println!("Noise: {}", unsafe { sk.measure_noise(&ct3)? });
			let pt = sk.decrypt(&ct3)?;
			assert_eq!(Vec::<u64>::try_decode(&pt, Encoding::Simd)?, expected);
		}
		Ok(())
	}

	#[test]
	fn test_seq_mul_50() -> Result<(), String> {
		let par = Rc::new(
			BfvParametersBuilder::default()
				.polynomial_degree(8192)
				.plaintext_modulus(65537)
				.ciphertext_moduli_sizes(vec![50, 50, 50, 50, 50])
				.build()
				.unwrap(),
		);

		let values = par.plaintext.random_vec(par.polynomial_degree);
		let sk = SecretKey::random(&par);
		let rk = RelinearizationKey::new(&sk)?;
		let pt = Plaintext::try_encode(&values as &[u64], Encoding::Simd, &par)?;
		let mut ct1 = sk.encrypt(&pt)?;

		for _ in 0..5 {
			ct1 = mul(&ct1, &ct1, &rk)?;
			println!("Noise: {}", unsafe { sk.measure_noise(&ct1)? });
		}

		// Empirically measured.
		assert!(unsafe { sk.measure_noise(&ct1)? } <= 200);

		Ok(())
	}

	#[test]
	fn test_seq_mul_62() -> Result<(), String> {
		let par = Rc::new(
			BfvParametersBuilder::default()
				.polynomial_degree(8192)
				.plaintext_modulus(65537)
				.ciphertext_moduli_sizes(vec![62, 62, 62, 62, 62])
				.build()
				.unwrap(),
		);

		let values = par.plaintext.random_vec(par.polynomial_degree);
		let sk = SecretKey::random(&par);
		let rk = RelinearizationKey::new(&sk)?;
		let pt = Plaintext::try_encode(&values as &[u64], Encoding::Simd, &par)?;
		let mut ct1 = sk.encrypt(&pt)?;

		for _ in 0..5 {
			ct1 = mul(&ct1, &ct1, &rk)?;
			println!("Noise: {}", unsafe { sk.measure_noise(&ct1)? });
		}

		// Empirically measured.
		assert!(unsafe { sk.measure_noise(&ct1)? } <= 200);

		Ok(())
	}

	#[test]
	fn test_seq_mul2() -> Result<(), String> {
		let par = Rc::new(
			BfvParametersBuilder::default()
				.polynomial_degree(8192)
				.plaintext_modulus(65537)
				.ciphertext_moduli_sizes(vec![62, 62, 62, 62, 62])
				.build()
				.unwrap(),
		);

		let values = par.plaintext.random_vec(par.polynomial_degree);
		let sk = SecretKey::random(&par);
		let rk = RelinearizationKey::new(&sk)?;
		let pt = Plaintext::try_encode(&values as &[u64], Encoding::Simd, &par)?;
		let mut ct1 = sk.encrypt(&pt)?;

		for _ in 0..5 {
			ct1 = mul2(&ct1, &ct1, &rk)?;
			println!("Noise: {}", unsafe { sk.measure_noise(&ct1)? });
		}

		// Empirically measured.
		assert!(unsafe { sk.measure_noise(&ct1)? } <= 200);

		Ok(())
	}

	#[test]
	fn test_proto_conversion() -> Result<(), String> {
		for params in [
			Rc::new(BfvParameters::default(1)),
			Rc::new(BfvParameters::default(2)),
		] {
			let sk = SecretKey::random(&params);
			let v = params.plaintext.random_vec(8);
			let pt = Plaintext::try_encode(&v as &[u64], Encoding::Simd, &params)?;
			let ct = sk.encrypt(&pt)?;
			let ct_proto = CiphertextProto::from(&ct);
			assert_eq!(ct, Ciphertext::try_convert_from(&ct_proto, &params)?)
		}
		Ok(())
	}
}
