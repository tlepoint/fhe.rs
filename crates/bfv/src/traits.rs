//! Traits used for the BFV homomorphic encryption scheme.

use crate::ciphertext::Ciphertext;
use crate::parameters::BfvParameters;
use crate::plaintext::{Encoding, Plaintext};
use std::sync::Arc;

/// Encode values into [`Plaintext`].
pub trait Encoder<T>
where
	Self: Sized,
{
	/// The type of errors.
	type Error;

	/// Attempt to encode the `value` with the specified [`Encoding`] and [`BfvParameters`].
	fn try_encode(
		value: T,
		encoding: Encoding,
		par: &Arc<BfvParameters>,
	) -> Result<Self, Self::Error>;
}

/// Decode [`Plaintext`] values.
pub trait Decoder
where
	Self: Sized,
{
	/// The type of errors.
	type Error;

	/// Attempt to decode the [`Plaintext`], with an optional [`Encoding`].
	fn try_decode<E>(a: &Plaintext, encoding: E) -> Result<Self, Self::Error>
	where
		E: Into<Option<Encoding>>;
}

/// Encrypt a [`Plaintext`] into a [`Ciphertext`].
pub trait Encryptor {
	/// The type of errors.
	type Error;

	/// Encrypt a [`Plaintext`].
	fn encrypt(&self, plaintext: &Plaintext) -> Result<Ciphertext, Self::Error>;
}

/// Decrypt a [`Ciphertext`] into a [`Plaintext`].
pub trait Decryptor {
	/// The type of errors.
	type Error;

	/// Decrypt a [`Ciphertext`].
	fn decrypt(&mut self, ciphertext: &Ciphertext) -> Result<Plaintext, Self::Error>;
}

/// Conversions.
///
/// We unfortunaly cannot use the `TryFrom` trait from std::convert because we
/// need to specify additional parameters, and if we try to redefine a `TryFrom` trait here,
/// we need to fully specify the trait when we use it because of the blanket implementation
/// <https://github.com/rust-lang/rust/issues/50133#issuecomment-488512355>.
pub trait TryConvertFrom<T>
where
	Self: Sized,
{
	/// The type of errors.
	type Error;

	/// Attempt to convert the `value` with a specific parameter.
	fn try_convert_from(value: T, par: &Arc<BfvParameters>) -> Result<Self, Self::Error>;
}

/// Serialization.
pub trait Serialize {
	/// Serialize `Self` into a vector fo bytes.
	fn serialize(&self) -> Vec<u8>;
}

/// Deserialization using the specified parameters.
pub trait DeserializeWithParams
where
	Self: Sized,
{
	/// The type of errors.
	type Error;

	/// Attempt to deserialize from a vector of bytes
	fn try_deserialize(bytes: &[u8], par: &Arc<BfvParameters>) -> Result<Self, Self::Error>;
}

/// Deserialization without parameters.
pub trait Deserialize
where
	Self: Sized,
{
	/// The type of errors.
	type Error;

	/// Attempt to deserialize from a vector of bytes
	fn try_deserialize(bytes: &[u8]) -> Result<Self, Self::Error>;
}
