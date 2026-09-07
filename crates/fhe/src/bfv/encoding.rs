//! Interpretation of a BFV plaintext polynomial.

/// Interpretation supplied when encoding or decoding a plaintext.
/// It does not change the plaintext's modulus-switching level.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Encoding {
    /// Values are coefficients in a negacyclic polynomial ring.
    Polynomial,
    /// Values are slots with component-wise arithmetic. The current backend
    /// requires a machine-word prime plaintext modulus congruent to 1 modulo
    /// twice the polynomial degree.
    Simd,
}
