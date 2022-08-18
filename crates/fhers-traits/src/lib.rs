use std::sync::Arc;

/// The homomorphic encryption parameters.
pub trait FheParameters {}

/// A plaintext which will encode one (or more) value(s).
pub trait FhePlaintext
where
	Self: Sized,
{
}

/// Encoding used when encoding a [`FhePlaintext`].
pub trait FhePlaintextEncoding {}

/// Encode a value using a specified encoding.
pub trait FheEncoder<V, E: FhePlaintextEncoding, P: FheParameters>
where
	Self: FhePlaintext,
{
	type Error;
	fn try_encode(value: V, encoding: E, par: &Arc<P>) -> Result<Self, Self::Error>;
}

// Decode the value in the plaintext with the specified (optional) encoding.
pub trait FheDecoder<E: FhePlaintextEncoding, Pt: FhePlaintext>
where
	Self: Sized,
{
	type Error;
	fn try_decode<O>(pt: &Pt, encoding: O) -> Result<Self, Self::Error>
	where
		O: Into<Option<E>>;
}

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
