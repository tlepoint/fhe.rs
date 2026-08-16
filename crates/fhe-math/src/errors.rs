use thiserror::Error;

use crate::rq::Representation;

/// The Result type for this library.
pub type Result<T> = std::result::Result<T, Error>;

/// Enum encapsulation all the possible errors from this library.
#[derive(Debug, Error, PartialEq, Eq)]
#[expect(missing_docs, reason = "error variants are documented inline")]
#[non_exhaustive]
pub enum Error {
    /// Indicates an invalid modulus
    #[error("Invalid modulus: modulus {0} should be between 2 and (1 << 62) - 1.")]
    InvalidModulus(u64),

    /// Indicates invalid serialized polynomial data.
    #[error("Polynomial serialization error: {0}")]
    PolynomialSerialization(#[from] PolynomialSerializationError),

    /// Indicates that there is no more contexts to switch to.
    #[error("This is the last context.")]
    NoMoreContext,

    /// Indicates that a target context is not in the source context's chain.
    #[error("Target context is not reachable from the source context.")]
    ContextNotReachable,

    /// Indicates that a polynomial has the wrong context for an operation.
    #[error("Polynomial context does not match the operation context.")]
    PolynomialContextMismatch,

    /// Indicates that the seed size is incorrect.
    #[error("Invalid seed: got {0} bytes, expected {1} bytes.")]
    InvalidSeedSize(usize, usize),

    /// Indicates an invalid polynomial degree.
    #[error("Invalid polynomial degree {degree}: expected a power of two at least {minimum}.")]
    InvalidPolynomialDegree { degree: usize, minimum: usize },

    /// Indicates that an NTT operator is unavailable for a modulus and degree.
    #[error("No NTT operator for modulus {modulus} and degree {degree}.")]
    NttOperatorUnavailable { modulus: u64, degree: usize },

    /// Indicates that a context level is out of bounds.
    #[error("Context level {level} is out of bounds; maximum level is {max_level}.")]
    InvalidContextLevel { level: usize, max_level: usize },

    /// Indicates an invalid substitution exponent.
    #[error("Substitution exponent {exponent} must be odd modulo twice degree {degree}.")]
    InvalidSubstitutionExponent { exponent: usize, degree: usize },

    /// Indicates an invalid sampling variance.
    #[error("Invalid variance {variance}: expected a value in [{minimum}, {maximum}].")]
    InvalidVariance {
        variance: usize,
        minimum: usize,
        maximum: usize,
    },

    /// Indicates that a polynomial dot product received no operands.
    #[error("Polynomial dot product requires at least one operand.")]
    EmptyDotProduct,

    /// Indicates mismatched polynomial dot-product operand counts.
    #[error("Polynomial dot product length mismatch: left has {left}, right has {right}.")]
    DotProductLengthMismatch { left: usize, right: usize },

    /// Indicates an invalid number of polynomial coefficients.
    #[error(
        "Invalid coefficient count {actual} for {representation:?}: degree {degree}, moduli {moduli}."
    )]
    InvalidCoefficientCount {
        representation: Representation,
        actual: usize,
        degree: usize,
        moduli: usize,
    },

    /// Indicates an invalid two-dimensional coefficient shape.
    #[error(
        "Invalid coefficient shape ({actual_rows}, {actual_columns}); expected ({expected_rows}, {expected_columns})."
    )]
    InvalidCoefficientShape {
        actual_rows: usize,
        actual_columns: usize,
        expected_rows: usize,
        expected_columns: usize,
    },

    /// Indicates too many coefficients for a polynomial.
    #[error("Too many polynomial coefficients: found {actual}, maximum is {maximum}.")]
    TooManyCoefficients { actual: usize, maximum: usize },

    /// Indicates that polynomial coefficient storage is not contiguous.
    #[error("Polynomial coefficients are not contiguous in memory.")]
    NonContiguousCoefficients,

    /// Indicates incompatible polynomial degrees.
    #[error("Polynomial degree mismatch: found {found}, expected {expected}.")]
    DegreeMismatch { found: usize, expected: usize },

    /// Indicates that no RNS moduli were provided.
    #[error("The list of RNS moduli is empty.")]
    EmptyModuli,

    /// Indicates two RNS moduli are not coprime.
    #[error("RNS moduli {left} and {right} are not coprime.")]
    NonCoprimeModuli { left: u64, right: u64 },

    /// Indicates that a modular inverse does not exist.
    #[error("Value {value} is not invertible modulo {modulus}.")]
    NonInvertible { value: u64, modulus: u64 },
}

/// Errors arising while decoding or converting serialized polynomials.
#[derive(Debug, Error, PartialEq, Eq)]
#[expect(missing_docs, reason = "error variants are documented inline")]
#[non_exhaustive]
pub enum PolynomialSerializationError {
    /// The protobuf payload could not be decoded.
    #[error("Failed to decode polynomial serialization.")]
    Decode,

    /// The protobuf representation discriminant is invalid.
    #[error("Invalid polynomial representation value {value}.")]
    InvalidRepresentation { value: i32 },

    /// The protobuf representation was left unspecified.
    #[error("Polynomial representation is unknown.")]
    UnknownRepresentation,

    /// The serialized representation differs from the requested one.
    #[error("Incorrect representation: got {found:?}, expected {expected:?}.")]
    RepresentationMismatch {
        found: Representation,
        expected: Representation,
    },

    /// The serialized polynomial degree is invalid.
    #[error("Invalid polynomial degree {degree}.")]
    InvalidDegree { degree: usize },

    /// The serialized coefficient byte count is invalid.
    #[error("Invalid coefficient length: found {actual} bytes, expected {expected}.")]
    InvalidCoefficientCount { actual: usize, expected: usize },
}

#[cfg(test)]
mod tests {
    use crate::{Error, PolynomialSerializationError};

    #[test]
    fn error_strings() {
        assert_eq!(
            Error::InvalidModulus(0).to_string(),
            "Invalid modulus: modulus 0 should be between 2 and (1 << 62) - 1."
        );
        assert_eq!(
            Error::PolynomialSerialization(PolynomialSerializationError::Decode).to_string(),
            "Polynomial serialization error: Failed to decode polynomial serialization."
        );
        assert_eq!(
            Error::NoMoreContext.to_string(),
            "This is the last context."
        );
        assert_eq!(
            Error::ContextNotReachable.to_string(),
            "Target context is not reachable from the source context."
        );
        assert_eq!(
            Error::InvalidSeedSize(0, 1).to_string(),
            "Invalid seed: got 0 bytes, expected 1 bytes."
        );
    }
}
