use thiserror::Error;

use crate::rq::Representation;

/// The Result type for this library.
pub type Result<T> = std::result::Result<T, Error>;

/// Enum encapsulation all the possible errors from this library.
#[derive(Debug, Error, PartialEq, Eq)]
pub enum Error {
	/// Indicates an invalid modulus
	#[error("Invalid modulus: modulus {0} should be between 2 and (1 << 62) - 1.")]
	InvalidModulus(u64),

	/// Indicates an error in the serialization / deserialization.
	#[error("{0}")]
	Serialization(String),

	/// Indicates that there is no more contexts to switch to.
	#[error("This is the last context.")]
	NoMoreContext,

	/// Indicates that the provided context is invalid.
	#[error("Invalid context provided.")]
	InvalidContext,

	/// Indicates an incorrect representation.
	#[error("Incorrect representation: got {0:?}, expected {1:?}.")]
	IncorrectRepresentation(Representation, Representation),

	/// Indicates that the seed size is incorrect.
	#[error("Invalid seed: got {0} bytes, expected {1} bytes.")]
	InvalidSeedSize(usize, usize),

	/// Indicates a default error
	/// TODO: To delete when transition is over
	#[error("{0}")]
	Default(String),
}

#[cfg(test)]
mod tests {
	use crate::{rq::Representation, Error};

	#[test]
	fn error_strings() {
		assert_eq!(
			Error::InvalidModulus(0).to_string(),
			"Invalid modulus: modulus 0 should be between 2 and (1 << 62) - 1."
		);
		assert_eq!(Error::Serialization("test".to_string()).to_string(), "test");
		assert_eq!(
			Error::NoMoreContext.to_string(),
			"This is the last context."
		);
		assert_eq!(
			Error::InvalidContext.to_string(),
			"Invalid context provided."
		);
		assert_eq!(
			Error::IncorrectRepresentation(Representation::Ntt, Representation::NttShoup)
				.to_string(),
			"Incorrect representation: got Ntt, expected NttShoup."
		);
		assert_eq!(
			Error::InvalidSeedSize(0, 1).to_string(),
			"Invalid seed: got 0 bytes, expected 1 bytes."
		);
	}
}
