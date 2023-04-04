use thiserror::Error;

/// The Result type for this library.
pub type Result<T> = std::result::Result<T, Error>;

/// Enum encapsulating all the possible errors from this library.
#[derive(Debug, Error, PartialEq, Eq)]
pub enum Error {
    /// Indicates that an error from the underlying mathematical library was
    /// encountered.
    #[error("{0}")]
    MathError(fhe_math::Error),

    /// Indicates a serialization error.
    #[error("Serialization error")]
    SerializationError,

    /// Indicates that too many values were provided.
    #[error("Too many values provided: {0} exceeds limit {1}")]
    TooManyValues(usize, usize),

    /// Indicates that too few values were provided.
    #[error("Too few values provided: {0} is below limit {1}")]
    TooFewValues(usize, usize),

    /// Indicates that an input is invalid.
    #[error("{0}")]
    UnspecifiedInput(String),

    /// Indicates a mismatch in the encodings.
    #[error("Encoding mismatch: found {0}, expected {1}")]
    EncodingMismatch(String, String),

    /// Indicates that the encoding is not supported.
    #[error("Does not support {0} encoding")]
    EncodingNotSupported(String),

    /// Indicates a parameter error.
    #[error("{0}")]
    ParametersError(ParametersError),

    /// Indicates a default error
    /// TODO: To delete eventually
    #[error("{0}")]
    DefaultError(String),
}

impl From<fhe_math::Error> for Error {
    fn from(e: fhe_math::Error) -> Self {
        Error::MathError(e)
    }
}

/// Separate enum to indicate parameters-related errors.
#[derive(Debug, Error, PartialEq, Eq)]
pub enum ParametersError {
    /// Indicates that the degree is invalid.
    #[error("Invalid degree: {0} is not a power of 2 larger than 8")]
    InvalidDegree(usize),

    /// Indicates that the moduli sizes are invalid.
    #[error("Invalid modulus size: {0}, expected an integer between {1} and {2}")]
    InvalidModulusSize(usize, usize, usize),

    /// Indicates that there exists not enough primes of this size.
    #[error("Not enough primes of size {0} for polynomials of degree {1}")]
    NotEnoughPrimes(usize, usize),

    /// Indicates that the plaintext is invalid.
    #[error("{0}")]
    InvalidPlaintext(String),

    /// Indicates that too many parameters were specified.
    #[error("{0}")]
    TooManySpecified(String),

    /// Indicates that too few parameters were specified.
    #[error("{0}")]
    TooFewSpecified(String),
}

#[cfg(test)]
mod tests {
    use crate::{Error, ParametersError};

    #[test]
    fn error_strings() {
        assert_eq!(
            Error::MathError(fhe_math::Error::InvalidContext).to_string(),
            fhe_math::Error::InvalidContext.to_string()
        );
        assert_eq!(Error::SerializationError.to_string(), "Serialization error");
        assert_eq!(
            Error::TooManyValues(20, 17).to_string(),
            "Too many values provided: 20 exceeds limit 17"
        );
        assert_eq!(
            Error::TooFewValues(10, 17).to_string(),
            "Too few values provided: 10 is below limit 17"
        );
        assert_eq!(
            Error::UnspecifiedInput("test string".to_string()).to_string(),
            "test string"
        );
        assert_eq!(
            Error::EncodingMismatch("enc1".to_string(), "enc2".to_string()).to_string(),
            "Encoding mismatch: found enc1, expected enc2"
        );
        assert_eq!(
            Error::EncodingNotSupported("test".to_string()).to_string(),
            "Does not support test encoding"
        );
        assert_eq!(
            Error::ParametersError(ParametersError::InvalidDegree(10)).to_string(),
            ParametersError::InvalidDegree(10).to_string()
        );
    }

    #[test]
    fn parameters_error_strings() {
        assert_eq!(
            ParametersError::InvalidDegree(10).to_string(),
            "Invalid degree: 10 is not a power of 2 larger than 8"
        );
        assert_eq!(
            ParametersError::InvalidModulusSize(1, 2, 3).to_string(),
            "Invalid modulus size: 1, expected an integer between 2 and 3"
        );
        assert_eq!(
            ParametersError::NotEnoughPrimes(1, 2).to_string(),
            "Not enough primes of size 1 for polynomials of degree 2"
        );
        assert_eq!(
            ParametersError::InvalidPlaintext("test".to_string()).to_string(),
            "test"
        );
        assert_eq!(
            ParametersError::TooManySpecified("test".to_string()).to_string(),
            "test"
        );
        assert_eq!(
            ParametersError::TooFewSpecified("test".to_string()).to_string(),
            "test"
        );
    }
}
