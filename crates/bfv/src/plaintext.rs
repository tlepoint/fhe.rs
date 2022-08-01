//! Plaintext type in the BFV encryption scheme.

use crate::parameters::BfvParameters;
use crate::traits::{Decoder, Encoder};
use math::rq::traits::TryConvertFrom;
use math::rq::{Poly, Representation};
use std::rc::Rc;
use zeroize::Zeroize;

/// An encoding for the plaintext.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Encoding {
	/// A Poly encoding encodes a vector as coefficients of a polynomial;
	/// homomorphic operations are therefore polynomial operations.
	Poly,
	/// A Simd encoding encodes a vector so that homomorphic operations are
	/// component-wise operations on the coefficients of the underlying vectors. The
	/// Simd encoding require that the plaintext modulus is congruent to 1 modulo
	/// the degree of the underlying polynomial.
	Simd,
}

/// A plaintext object, that encodes a vector according to a specific encoding.
/// TODO: Can the members be private?
#[derive(Debug, Clone, PartialEq)]
pub struct Plaintext {
	/// The parameters of the underlying BFV encryption scheme.
	pub par: Rc<BfvParameters>,
	/// The value after encoding.
	pub value: Vec<u64>,
}

impl Plaintext {
	/// Returns the plaintext as a polynomial
	pub fn to_poly(&self) -> Poly {
		Poly::try_convert_from(&self.value, self.par.ctx(), Representation::PowerBasis).unwrap()
	}
}

impl TryConvertFrom<&Plaintext> for Poly {
	type Error = String;

	fn try_convert_from<R>(
		pt: &Plaintext,
		ctx: &Rc<math::rq::Context>,
		_: R,
	) -> Result<Self, Self::Error>
	where
		R: Into<Option<Representation>>,
	{
		if ctx != pt.par.ctx() {
			Err("Incompatible contexts".to_string())
		} else {
			Poly::try_convert_from(&pt.value, ctx, Representation::PowerBasis)
		}
	}
}

impl Zeroize for Plaintext {
	fn zeroize(&mut self) {
		self.value.zeroize()
	}
}

impl Encoder<&[u64]> for Plaintext {
	type Error = String;

	fn try_encode(
		value: &[u64],
		encoding: Encoding,
		par: &Rc<BfvParameters>,
	) -> Result<Self, Self::Error> {
		if value.len() > par.degree() {
			return Err("There are too many values to encode".to_string());
		}

		let mut v = vec![0; par.degree()];
		v[..value.len()].copy_from_slice(value);

		match encoding {
			Encoding::Poly => Ok(Self {
				par: par.clone(),
				value: v,
			}),
			Encoding::Simd => {
				if let Some(op) = par.simd_operator() {
					op.backward(&mut v);
					Ok(Self {
						par: par.clone(),
						value: v,
					})
				} else {
					v.zeroize();
					Err("The plaintext does not allow for using the Simd encoding".to_string())
				}
			}
		}
	}
}

impl Decoder for Vec<u64> {
	type Error = String;
	fn try_decode(pt: &Plaintext, encoding: Encoding) -> Result<Vec<u64>, Self::Error> {
		let mut v = pt.value.clone();
		if encoding == Encoding::Simd {
			if let Some(op) = pt.par.simd_operator() {
				op.forward(&mut v);
				Ok(v)
			} else {
				Err("The plaintext does not allow for using the Simd encoding".to_string())
			}
		} else {
			Ok(v)
		}
	}
}

#[cfg(test)]
mod tests {
	use super::{Encoding, Plaintext};
	use crate::parameters::{BfvParameters, BfvParametersBuilder};
	use crate::traits::{Decoder, Encoder};
	use math::zq::Modulus;
	use proptest::collection::vec as prop_vec;
	use proptest::prelude::any;
	use std::rc::Rc;

	#[test]
	fn try_encode() {
		// The default test parameters support both Poly and Simd encodings
		let params = Rc::new(BfvParameters::default_one_modulus());
		let q = Modulus::new(params.plaintext()).unwrap();
		let a = q.random_vec(params.degree());

		let plaintext =
			Plaintext::try_encode(&[0, 0, 0, 0, 0, 0, 0, 0, 0], Encoding::Poly, &params);
		assert!(plaintext.is_err());

		let plaintext = Plaintext::try_encode(&a, Encoding::Poly, &params);
		assert!(plaintext.is_ok_and(|pt| pt.value == a));

		let plaintext = Plaintext::try_encode(&a, Encoding::Simd, &params);
		assert!(plaintext.is_ok());

		let plaintext = Plaintext::try_encode(&[1], Encoding::Poly, &params);
		assert!(plaintext.is_ok_and(|pt| pt.value == &[1, 0, 0, 0, 0, 0, 0, 0]));

		// The following parameters do not allow for Simd encoding
		let params = Rc::new(
			BfvParametersBuilder::default()
				.polynomial_degree(8)
				.plaintext_modulus(2)
				.ciphertext_moduli(vec![4611686018326724609])
				.build()
				.unwrap(),
		);

		let plaintext = Plaintext::try_encode(&a, Encoding::Poly, &params);
		assert!(plaintext.is_ok_and(|pt| pt.value == a));

		let plaintext = Plaintext::try_encode(&a, Encoding::Simd, &params);
		assert!(plaintext.is_err());
	}

	proptest! {
		#[test]
		fn test_encode_decode(mut a in prop_vec(any::<u64>(), 8)) {
			let params = Rc::new(BfvParameters::default_one_modulus());
			let q = Modulus::new(params.plaintext()).unwrap();
			q.reduce_vec(&mut a);

			let plaintext = Plaintext::try_encode(&a, Encoding::Poly, &params);
			assert!(plaintext.is_ok_and(|pt| pt.value == a));
			let b = Vec::<u64>::try_decode(&plaintext.unwrap(), Encoding::Poly);
			assert!(b.is_ok_and(|b| b == &a));

			let plaintext = Plaintext::try_encode(&a, Encoding::Simd, &params);
			assert!(plaintext.is_ok());
			let b = Vec::<u64>::try_decode(&plaintext.unwrap(), Encoding::Simd);
			assert!(b.is_ok_and(|b| b == &a));
		}
	}
}
