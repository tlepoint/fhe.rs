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

// While switching to an Error enum, we provide a mapping to String.
// TODO: Remove when the transition is complete.
impl From<Error> for String {
	fn from(e: Error) -> Self {
		e.to_string()
	}
}

impl From<String> for Error {
	fn from(s: String) -> Self {
		Error::Default(s)
	}
}
