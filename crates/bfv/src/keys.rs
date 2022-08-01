//! Keys for the BFV encryption scheme

use crate::parameters::BfvParameters;
use math::rq::{traits::TryConvertFrom, Poly, Representation};
use std::rc::Rc;
use util::sample_vec_cbd;
use zeroize::{Zeroize, ZeroizeOnDrop, Zeroizing};

/// Secret key for the BFV encryption scheme.
#[derive(Debug, PartialEq)]
pub struct SecretKey {
	par: Rc<BfvParameters>,
	s: Poly,
}

impl Zeroize for SecretKey {
	fn zeroize(&mut self) {
		self.s.zeroize();
	}
}

impl ZeroizeOnDrop for SecretKey {}

impl SecretKey {
	/// Generate a random [`SecretKey`].
	pub fn random(par: &Rc<BfvParameters>) -> Self {
		let s_coefficients = Zeroizing::new(sample_vec_cbd(par.degree(), par.variance()).unwrap());
		let mut s = Poly::try_convert_from(
			s_coefficients.as_ref() as &[i64],
			par.ctx(),
			Representation::PowerBasis,
		)
		.unwrap();
		s.change_representation(Representation::NttShoup);
		Self {
			par: par.clone(),
			s,
		}
	}
}

#[cfg(test)]
mod tests {
	use super::SecretKey;
	use crate::parameters::BfvParameters;
	use math::rq::Representation;
	use std::rc::Rc;

	#[test]
	fn test_keygen() {
		let params = Rc::new(BfvParameters::default());
		let sk = SecretKey::random(&params);
		assert_eq!(sk.par, params);

		let mut s = sk.s.clone();
		s.change_representation(Representation::PowerBasis);
		let coefficients = Vec::<u64>::from(&s);
		coefficients.iter().for_each(|ci| {
			// Check that this is a small polynomial
			assert!(
				*ci <= 2 * sk.par.variance() as u64
					|| *ci >= (sk.par.moduli()[0] - 2 * sk.par.variance() as u64)
			)
		})
	}
}
