use thiserror::Error;

/// The Result type for this library.
pub type Result<T> = std::result::Result<T, Error>;

/// Enum encapsulation all the possible errors from this library.
#[derive(Debug, Error, PartialEq, Eq)]
pub enum Error {
	/// Indicates an invalid modulus
	#[error("Invalid modulus: modulus {0} should be between 2 and (1 << 62) - 1")]
	InvalidModulus(u64),

	/// Indicates an error in the serialization / deserialization.
	#[error("{0}")]
	Serialization(String),

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
