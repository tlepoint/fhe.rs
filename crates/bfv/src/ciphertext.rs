//! Ciphertext type in the BFV encryption scheme.

use crate::parameters::BfvParameters;
use math::rq::Poly;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use std::{
	ops::{Add, AddAssign, Neg, Sub, SubAssign},
	rc::Rc,
};

/// A ciphertext encrypting a plaintext.
#[derive(Debug, Clone, PartialEq)]
pub struct Ciphertext {
	/// The parameters of the underlying BFV encryption scheme.
	pub(crate) par: Rc<BfvParameters>,

	/// The seed that generated the polynomial c0 in a fresh ciphertext.
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

#[cfg(test)]
mod tests {
	use crate::{
		traits::{Decoder, Decryptor, Encoder, Encryptor},
		BfvParameters, Encoding, Plaintext, SecretKey,
	};
	use proptest::collection::vec as prop_vec;
	use proptest::prelude::any;
	use std::rc::Rc;

	proptest! {
		#[test]
		fn test_add(mut a in prop_vec(any::<u64>(), 8), mut b in prop_vec(any::<u64>(), 8)) {
			for params in [Rc::new(BfvParameters::default_one_modulus()),
						   Rc::new(BfvParameters::default_two_moduli())] {
				params.plaintext().reduce_vec(&mut a);
				params.plaintext().reduce_vec(&mut b);
				let mut c = a.clone();
				params.plaintext().add_vec(&mut c, &b);

				let sk = SecretKey::random(&params);

				for encoding in [Encoding::Poly, Encoding::Simd] {
					let pt_a = Plaintext::try_encode(&a, encoding.clone(), &params).unwrap();
					let pt_b = Plaintext::try_encode(&b, encoding.clone(), &params).unwrap();

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
			for params in [Rc::new(BfvParameters::default_one_modulus()),
						   Rc::new(BfvParameters::default_two_moduli())] {
				params.plaintext().reduce_vec(&mut a);
				params.plaintext().reduce_vec(&mut b);
				let mut c = a.clone();
				params.plaintext().sub_vec(&mut c, &b);

				let sk = SecretKey::random(&params);

				for encoding in [Encoding::Poly, Encoding::Simd] {
					let pt_a = Plaintext::try_encode(&a, encoding.clone(), &params).unwrap();
					let pt_b = Plaintext::try_encode(&b, encoding.clone(), &params).unwrap();

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
			for params in [Rc::new(BfvParameters::default_one_modulus()),
						   Rc::new(BfvParameters::default_two_moduli())] {
				params.plaintext().reduce_vec(&mut a);
				let mut c = a.clone();
				params.plaintext().neg_vec(&mut c);

				let sk = SecretKey::random(&params);
				for encoding in [Encoding::Poly, Encoding::Simd] {
					let pt_a = Plaintext::try_encode(&a, encoding.clone(), &params).unwrap();

					let ct_a = sk.encrypt(&pt_a).unwrap();
					let ct_c = -&ct_a;

					let pt_c = sk.decrypt(&ct_c).unwrap();
					assert_eq!(Vec::<u64>::try_decode(&pt_c, encoding.clone()).unwrap(), c);
				}
			}
		}
	}
}
