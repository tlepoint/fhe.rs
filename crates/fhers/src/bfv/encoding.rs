//! The encoding type for BFV.

use fhers_traits::FhePlaintextEncoding;

/// An encoding for the plaintext.
#[derive(Debug, Clone, Eq, PartialEq)]
pub enum Encoding {
	/// \[Advanced\] A Poly encoding embedding the level information.
	PolyLeveled(usize),
	/// \[Advanced\] A Simd encoding embedding the level information.
	SimdLeveled(usize),
}

impl Encoding {
	/// A Poly encoding encodes a vector as coefficients of a polynomial;
	/// homomorphic operations are therefore polynomial operations.
	pub fn poly() -> Self {
		Self::PolyLeveled(0)
	}

	/// A Simd encoding encodes a vector so that homomorphic operations are
	/// component-wise operations on the coefficients of the underlying vectors.
	/// The Simd encoding require that the plaintext modulus is congruent to 1
	/// modulo the degree of the underlying polynomial.
	pub fn simd() -> Self {
		Self::SimdLeveled(0)
	}
}

impl From<Encoding> for String {
	fn from(e: Encoding) -> Self {
		String::from(&e)
	}
}

impl From<&Encoding> for String {
	fn from(e: &Encoding) -> Self {
		format!("{:?}", e)
	}
}

impl FhePlaintextEncoding for Encoding {}
