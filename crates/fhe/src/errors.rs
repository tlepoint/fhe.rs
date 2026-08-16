#![expect(
    missing_docs,
    reason = "error enums rely on variant docs and error messages"
)]

use num_bigint::BigUint;
use thiserror::Error;

use crate::bfv::Encoding;

/// The Result type for this library.
pub type Result<T> = std::result::Result<T, Error>;

/// Enum encapsulating all the possible errors from this library.
#[derive(Debug, Error, PartialEq, Eq)]
#[expect(missing_docs, reason = "error variants are documented inline")]
#[non_exhaustive]
pub enum Error {
    /// An error from the underlying mathematical library.
    #[error("Math library error: {0}")]
    MathError(#[from] fhe_math::Error),

    /// Cryptographic objects were constructed with incompatible parameters.
    #[error("Parameter mismatch between {left:?} and {right:?}")]
    ParameterMismatch {
        left: ParameterSource,
        right: ParameterSource,
    },

    /// A level is outside the requested range.
    #[error("Level {level} out of bounds: valid range is [{min_level}, {max_level}]")]
    InvalidLevel {
        level: usize,
        min_level: usize,
        max_level: usize,
    },

    /// A ciphertext is structurally invalid for the requested operation.
    #[error("Ciphertext error: {0}")]
    Ciphertext(#[from] CiphertextError),

    /// A plaintext is structurally invalid for the requested operation.
    #[error("Plaintext error: {0}")]
    Plaintext(#[from] PlaintextError),

    /// A plaintext encoding is invalid or unavailable.
    #[error("Encoding error: {0}")]
    Encoding(#[from] EncodingError),

    /// A parameter error.
    #[error("Parameters error: {0}")]
    ParametersError(#[from] ParametersError),

    /// A serialization error.
    #[error("Serialization error: {0}")]
    SerializationError(#[from] SerializationError),

    /// An evaluation-key or key-switching error.
    #[error("Evaluation-key error: {0}")]
    EvaluationKey(#[from] EvaluationKeyError),

    /// A ciphertext/plaintext dot-product error.
    #[error("Dot-product error: {0}")]
    DotProduct(#[from] DotProductError),

    /// A multiparty BFV protocol error.
    #[error("Multiparty protocol error: {0}")]
    Multiparty(#[from] MultipartyError),
}

/// Identifies the role of an object in a parameter mismatch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum ParameterSource {
    Ciphertext,
    Plaintext,
    Parameters,
    SecretKey,
    InputSecretKey,
    OutputSecretKey,
    PublicKey,
    Polynomial,
    KeySwitchingKey,
    RelinearizationKey,
    Multiplicator,
}

/// Ciphertext validation failures.
#[derive(Debug, Error, PartialEq, Eq)]
#[expect(missing_docs, reason = "error variants are documented inline")]
#[non_exhaustive]
pub enum CiphertextError {
    #[error("Expected at least {minimum} polynomials, found {actual}")]
    TooFewPolynomials { actual: usize, minimum: usize },

    #[error("Polynomial context does not match ciphertext level {level}")]
    PolynomialContextMismatch { level: usize },

    #[error("{operation:?} requires {expected} polynomials, found {actual}")]
    InvalidPolynomialCount {
        operation: CiphertextOperation,
        actual: usize,
        expected: usize,
    },

    #[error("Multiplication requires {expected} polynomials per operand, found {left} and {right}")]
    MultiplicationPolynomialCount {
        left: usize,
        right: usize,
        expected: usize,
    },
}

/// Operation that imposes a particular ciphertext polynomial count.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum CiphertextOperation {
    Galois,
    EvaluationKey,
    Relinearization,
    MultipartyKeySwitch,
}

/// Plaintext validation and conversion failures.
#[derive(Debug, Error, PartialEq, Eq)]
#[expect(missing_docs, reason = "error variants are documented inline")]
#[non_exhaustive]
pub enum PlaintextError {
    #[error("Polynomial context does not match plaintext level {level}")]
    PolynomialContextMismatch { level: usize },

    #[error("No plaintext encoding was specified")]
    MissingEncoding,

    #[error("No NTT operator is available for the plaintext parameters")]
    NttOperatorUnavailable,

    #[error("Plaintext input has {actual} values but degree permits at most {maximum}")]
    TooManyValues { actual: usize, maximum: usize },

    #[error("Plaintext value does not fit in u64")]
    ValueTooLargeForU64,
}

/// Plaintext encoding failures.
#[derive(Debug, Error, PartialEq, Eq)]
#[expect(missing_docs, reason = "error variants are documented inline")]
#[non_exhaustive]
pub enum EncodingError {
    #[error("Encoding mismatch: found {found:?}, expected {expected:?}")]
    Mismatch { found: Encoding, expected: Encoding },

    #[error("SIMD encoding requires an NTT operator for the plaintext modulus")]
    SimdUnavailable,
}

/// Evaluation-key and key-switching failures.
#[derive(Debug, Error, PartialEq, Eq)]
#[expect(missing_docs, reason = "error variants are documented inline")]
#[non_exhaustive]
pub enum EvaluationKeyError {
    #[error("The parameters do not support key switching")]
    KeySwitchingNotSupported,

    #[error("Evaluation key does not support {operation:?}")]
    Unsupported { operation: EvaluationOperation },

    #[error("Evaluation key is missing {component:?}")]
    Missing { component: EvaluationKeyComponent },

    #[error("Key-switching key contains no components")]
    EmptyKeySwitchingComponents,

    #[error("Invalid rotation step {step}: must be in range [{min}, {max}]")]
    InvalidRotationStep { step: usize, min: usize, max: usize },

    #[error("Expansion size {size} must be in [1, {degree}]")]
    InvalidExpansionSize { size: usize, degree: usize },
}

/// Evaluation-key operation requested by a caller.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum EvaluationOperation {
    InnerSum,
    RowRotation,
    ColumnRotation { step: usize },
    Expansion { level: usize },
}

/// Material required by an evaluation-key operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum EvaluationKeyComponent {
    GaloisExponent { step: usize },
    GaloisKey { element: usize },
    ExpansionMonomial { level: usize },
}

/// Dot-product validation failures.
#[derive(Debug, Error, PartialEq, Eq)]
#[expect(missing_docs, reason = "error variants are documented inline")]
#[non_exhaustive]
pub enum DotProductError {
    #[error("Dot product requires at least one operand")]
    EmptyInput,

    #[error("Received {ciphertexts} ciphertexts and {plaintexts} plaintexts")]
    OperandCountMismatch {
        ciphertexts: usize,
        plaintexts: usize,
    },

    #[error("Ciphertext has {actual} polynomial parts; expected {expected}")]
    CiphertextPolynomialCountMismatch { actual: usize, expected: usize },
}

/// Multiparty BFV protocol failures.
#[derive(Debug, Error, PartialEq, Eq)]
#[expect(missing_docs, reason = "error variants are documented inline")]
#[non_exhaustive]
pub enum MultipartyError {
    #[error("At least one protocol share is required")]
    NoShares,

    #[error("Expected {expected} common random polynomials, got {actual}")]
    InvalidCommonRandomPolynomialCount { actual: usize, expected: usize },

    #[error("Round-two relinearization share is missing its round-one aggregation")]
    MissingRelinearizationRoundOneShare,
}

/// Separate enum for errors arising from serialization.
#[derive(Debug, Error, PartialEq, Eq)]
#[expect(missing_docs, reason = "error variants are documented inline")]
#[non_exhaustive]
pub enum SerializationError {
    /// A protobuf payload could not be decoded.
    #[error("Failed to decode {object:?}")]
    Decode { object: SerializedObject },

    /// A required protobuf field is absent.
    #[error("Missing required field {field:?}")]
    MissingField { field: SerializedField },

    /// A serialized polynomial collection has the wrong length.
    #[error("{component:?} has {actual} polynomials; expected {expected}")]
    WrongPolynomialCount {
        component: SerializedPolynomialComponent,
        expected: usize,
        actual: usize,
    },

    /// Key-switching levels encoded in an RGSW ciphertext disagree.
    #[error("Serialized RGSW ciphertext has inconsistent key-switching levels")]
    InconsistentKeySwitchingLevels,

    /// A public key contains a ciphertext at a nonzero level.
    #[error("Serialized public key ciphertext has level {actual}; expected {expected}")]
    InvalidPublicKeyLevel { actual: usize, expected: usize },

    /// Decomposition was encoded for non-maximal key-switching levels.
    #[error(
        "Key-switching decomposition requires level {expected}; ciphertext level is {ciphertext_level} and key level is {key_level}"
    )]
    InvalidKeySwitchingDecompositionLevels {
        ciphertext_level: usize,
        key_level: usize,
        expected: usize,
    },

    /// A serialized key-switching seed has the wrong length.
    #[error("Serialized key-switching seed has {actual} bytes; expected {expected}")]
    InvalidKeySwitchingSeedLength { actual: usize, expected: usize },

    /// A serialized secret key has the wrong coefficient count.
    #[error("Serialized secret key has {actual} coefficients; expected {expected}")]
    InvalidSecretKeyCoefficientCount { actual: usize, expected: usize },

    /// A serialized ciphertext lacks enough explicit or seeded polynomials.
    #[error(
        "Serialized ciphertext has {actual} explicit polynomials (seed present: {seed_present})"
    )]
    InvalidCiphertextPolynomialCount { actual: usize, seed_present: bool },

    /// An I/O operation failed.
    #[error("Serialization I/O error: {kind:?}")]
    Io { kind: std::io::ErrorKind },
}

impl From<std::io::Error> for SerializationError {
    fn from(error: std::io::Error) -> Self {
        SerializationError::Io { kind: error.kind() }
    }
}

/// Type of protobuf object being decoded.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum SerializedObject {
    Ciphertext,
    EvaluationKey,
    Parameters,
    PublicKey,
    RelinearizationKey,
    RgswCiphertext,
    SecretKey,
}

/// Required field in a protobuf object.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum SerializedField {
    GaloisKeySwitchingKey,
    ParametersPlaintextModulus,
    PublicKeyCiphertext,
    RelinearizationKeySwitchingKey,
    RgswKeySwitchingKey0,
    RgswKeySwitchingKey1,
}

/// Polynomial collection encoded inside a protobuf object.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum SerializedPolynomialComponent {
    KeySwitchingKeyC0,
    KeySwitchingKeyC1,
}

/// Separate enum to indicate parameters-related errors.
#[derive(Debug, Error, PartialEq, Eq)]
#[expect(missing_docs, reason = "error variants are documented inline")]
#[non_exhaustive]
pub enum ParametersError {
    /// Indicates that the degree is invalid.
    #[error("Invalid polynomial degree {degree}: must be a power of 2 between {min} and {max}")]
    InvalidDegree {
        degree: usize,
        min: usize,
        max: usize,
    },

    /// Indicates that the plaintext modulus is invalid.
    #[error("Invalid plaintext modulus {modulus}: {source}")]
    InvalidPlaintextModulus {
        modulus: u64,
        source: fhe_math::Error,
    },

    /// Indicates that a ciphertext modulus is invalid.
    #[error("Invalid ciphertext modulus at index {index}: {modulus} ({source})")]
    InvalidCiphertextModulus {
        index: usize,
        modulus: u64,
        source: fhe_math::Error,
    },

    /// Indicates that the moduli sizes are invalid.
    #[error("Invalid modulus size at index {index}: {size}, expected between {min} and {max}")]
    InvalidModulusSize {
        index: usize,
        size: usize,
        min: usize,
        max: usize,
    },

    /// Indicates that there are not enough primes of a given size
    #[error(
        "Not enough primes of size {size} for degree {degree}: need {needed}, found {available}"
    )]
    NotEnoughPrimes {
        size: usize,
        degree: usize,
        needed: usize,
        available: usize,
    },

    /// Indicates duplicate moduli
    #[error("Duplicate moduli detected: {modulus} appears at indices {indices:?}")]
    DuplicateModuli { modulus: u64, indices: Vec<usize> },

    /// Indicates moduli are not coprime
    #[error("Moduli {modulus1} and {modulus2} are not coprime (gcd = {gcd})")]
    ModuliNotCoprime {
        modulus1: u64,
        modulus2: u64,
        gcd: u64,
    },

    /// Indicates plaintext modulus is not NTT-friendly
    #[error("Plaintext modulus {modulus} is not NTT-friendly for degree {degree}")]
    PlaintextNotNttFriendly { modulus: u64, degree: usize },

    /// Indicates ciphertext modulus is not NTT-friendly
    #[error(
        "Ciphertext modulus {modulus} at index {index} is not NTT-friendly for degree {degree}"
    )]
    CiphertextModulusNotNttFriendly {
        index: usize,
        modulus: u64,
        degree: usize,
    },

    /// Indicates plaintext modulus is too large relative to ciphertext moduli
    #[error(
        "Plaintext modulus {plaintext_modulus} exceeds ciphertext modulus {ciphertext_modulus} at index {index}"
    )]
    PlaintextModulusTooLarge {
        plaintext_modulus: u64,
        ciphertext_modulus: u64,
        index: usize,
    },

    /// Indicates that the plaintext modulus is not smaller than the complete
    /// ciphertext modulus.
    #[error(
        "Plaintext modulus {plaintext_modulus} must be smaller than ciphertext modulus {ciphertext_modulus}"
    )]
    PlaintextModulusExceedsCiphertextModulus {
        plaintext_modulus: BigUint,
        ciphertext_modulus: BigUint,
    },

    /// Indicates that the plaintext modulus is not coprime with a ciphertext
    /// modulus.
    #[error(
        "Plaintext modulus {plaintext_modulus} and ciphertext modulus {ciphertext_modulus} at index {index} are not coprime (gcd = {gcd})"
    )]
    PlaintextModulusNotCoprime {
        plaintext_modulus: BigUint,
        ciphertext_modulus: u64,
        index: usize,
        gcd: u64,
    },

    /// Indicates that reducing the plaintext modulus modulo a ciphertext
    /// modulus failed.
    #[error("Failed to reduce plaintext modulus modulo ciphertext modulus {ciphertext_modulus}")]
    PlaintextReductionFailed { ciphertext_modulus: u64 },

    /// Indicates insecure parameters according to standard
    #[error(
        "Parameters provide insufficient security: estimated security level {actual} bits, minimum required {minimum} bits"
    )]
    InsufficientSecurity { actual: u32, minimum: u32 },

    /// Indicates variance parameter out of range
    #[error("Invalid variance {variance}: must be between {min} and {max}")]
    InvalidVariance {
        variance: usize,
        min: usize,
        max: usize,
    },

    /// Indicates that explicit ciphertext moduli and modulus sizes were both
    /// specified.
    #[error("Specify either ciphertext moduli or ciphertext modulus sizes, not both")]
    ConflictingCiphertextModulusSpecifications,

    /// Indicates that neither explicit ciphertext moduli nor modulus sizes
    /// were specified.
    #[error("Ciphertext moduli or ciphertext modulus sizes must be specified")]
    MissingCiphertextModulusSpecification,

    /// Indicates no default parameter set can accommodate a plaintext size.
    #[error("No default parameters support a {plaintext_bits}-bit plaintext modulus")]
    NoDefaultParameters { plaintext_bits: usize },
}

impl ParametersError {
    #[must_use]
    pub fn invalid_degree_with_bounds(degree: usize) -> Self {
        Self::InvalidDegree {
            degree,
            min: 8,
            max: 65536,
        }
    }

    #[must_use]
    pub fn insufficient_security(actual: u32) -> Self {
        Self::InsufficientSecurity {
            actual,
            minimum: 128,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{
        Error, EvaluationKeyError, EvaluationOperation, MultipartyError, ParameterSource,
        ParametersError, SerializationError, SerializedObject,
    };

    #[test]
    fn error_strings() {
        assert_eq!(
            Error::MathError(fhe_math::Error::ContextNotReachable).to_string(),
            "Math library error: Target context is not reachable from the source context."
        );
        assert_eq!(
            Error::ParameterMismatch {
                left: ParameterSource::Ciphertext,
                right: ParameterSource::Parameters,
            }
            .to_string(),
            "Parameter mismatch between Ciphertext and Parameters"
        );
        assert_eq!(
            Error::SerializationError(SerializationError::Decode {
                object: SerializedObject::Ciphertext,
            })
            .to_string(),
            "Serialization error: Failed to decode Ciphertext"
        );
        assert_eq!(
            Error::ParametersError(ParametersError::invalid_degree_with_bounds(10)).to_string(),
            "Parameters error: Invalid polynomial degree 10: must be a power of 2 between 8 and 65536"
        );
        assert_eq!(
            Error::EvaluationKey(EvaluationKeyError::Unsupported {
                operation: EvaluationOperation::ColumnRotation { step: 3 },
            })
            .to_string(),
            "Evaluation-key error: Evaluation key does not support ColumnRotation { step: 3 }"
        );
        assert_eq!(
            Error::Multiparty(MultipartyError::InvalidCommonRandomPolynomialCount {
                actual: 2,
                expected: 3,
            })
            .to_string(),
            "Multiparty protocol error: Expected 3 common random polynomials, got 2"
        );
    }
}
