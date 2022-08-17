/// Serialization.
pub trait Serialize {
	/// Serialize `Self` into a vector fo bytes.
	fn serialize(&self) -> Vec<u8>;
}

/// Deserialization using the specified parameters.
pub trait DeserializeWithContext<T>
where
	Self: Sized,
{
	type Error;

	/// Attempt to deserialize from a vector of bytes
	fn try_deserialize(bytes: &[u8], ctx: T) -> Result<Self, Self::Error>;
}

/// Deserialization without context.
pub trait Deserialize
where
	Self: Sized,
{
	type Error;

	/// Attempt to deserialize from a vector of bytes
	fn try_deserialize(bytes: &[u8]) -> Result<Self, Self::Error>;
}
