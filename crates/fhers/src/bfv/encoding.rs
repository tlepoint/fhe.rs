//! The encoding type for BFV.

use fhers_traits::FhePlaintextEncoding;

/// An encoding for the plaintext.
#[derive(Debug, Clone, Eq, PartialEq)]
pub enum Encoding {
	/// A Poly encoding encodes a vector as coefficients of a polynomial;
	/// homomorphic operations are therefore polynomial operations.
	PolyLeveled(usize),
	/// A Simd encoding encodes a vector so that homomorphic operations are
	/// component-wise operations on the coefficients of the underlying vectors.
	/// The Simd encoding require that the plaintext modulus is congruent to 1
	/// modulo the degree of the underlying polynomial.
	SimdLeveled(usize),
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
