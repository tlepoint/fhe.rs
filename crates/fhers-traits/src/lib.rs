#![crate_name = "fhers_traits"]
#![crate_type = "lib"]
#![warn(missing_docs, unused_imports)]

//! Traits for Fully Homomorphic Encryption

use std::sync::Arc;

/// The homomorphic encryption parameters.
pub trait FheParameters {}

/// Encoding used when encoding a [`FhePlaintext`].
pub trait FhePlaintextEncoding {}

/// A plaintext which will encode one (or more) value(s).
pub trait FhePlaintext
where
	Self: Sized,
{
	/// The type of the FHE parameters.
	type Parameters: FheParameters;

	/// The type of the encoding.
	type Encoding: FhePlaintextEncoding;
}

/// Encode a value using a specified encoding.
pub trait FheEncoder<V>
where
	Self: FhePlaintext,
{
	/// The type of error returned.
	type Error;

	/// Attempt to encode a value using a specified encoding.
	fn try_encode(
		value: V,
		encoding: Self::Encoding,
		par: &Arc<Self::Parameters>,
	) -> Result<Self, Self::Error>;
}

/// Decode the value in the plaintext with the specified (optional) encoding.
pub trait FheDecoder<P: FhePlaintext>
where
	Self: Sized,
{
	/// The type of error returned.
	type Error;

	/// Attempt to decode a [`FhePlaintext`] into a value, using an (optional)
	/// encoding.
	fn try_decode<O>(pt: &P, encoding: O) -> Result<Self, Self::Error>
	where
		O: Into<Option<P::Encoding>>;
}

/// A ciphertext which will encrypt a plaintext.
pub trait FheCiphertext
where
	Self: Sized + Serialize + DeserializeUsingParameters,
{
	/// The type of the FHE parameters.
	type Parameters: FheParameters;
}

/// Encrypt a plaintext into a ciphertext.
pub trait FheEncrypter<
	P: FhePlaintext<Parameters = Self::Parameters>,
	C: FheCiphertext<Parameters = Self::Parameters>,
>
{
	/// The type of error returned.
	type Error;

	/// The type of the FHE parameters used by the encrypter.
	type Parameters: FheParameters;

	/// Try to encrypt an [`FhePlaintext`] into an [`FheCiphertext`].
	fn try_encrypt(&self, pt: &P) -> Result<C, Self::Error>;
}

/// Decrypt a ciphertext into a plaintext
pub trait FheDecrypter<
	P: FhePlaintext<Parameters = Self::Parameters>,
	C: FheCiphertext<Parameters = Self::Parameters>,
>
{
	/// The type of error returned.
	type Error;

	/// The type of the FHE parameters used by the decrypter.
	type Parameters: FheParameters;

	/// Try to decrypt an [`FheCiphertext`] into an [`FhePlaintext`].
	fn try_decrypt(&mut self, ct: &C) -> Result<P, Self::Error>;
}

/// Serialization.
pub trait Serialize {
	/// Serialize `Self` into a vector fo bytes.
	fn to_bytes(&self) -> Vec<u8>;
}

/// Deserialization using the specified parameters.
pub trait DeserializeUsingParameters
where
	Self: Sized,
{
	/// The type of error returned.
	type Error;

	/// The type of the parameters used to deserialize. Note that these are not
	/// necessary FheParameters.
	type Parameters;

	/// Attempt to deserialize from a vector of bytes
	fn from_bytes(bytes: &[u8], par: &Arc<Self::Parameters>) -> Result<Self, Self::Error>;
}

/// Deserialization without context.
pub trait Deserialize
where
	Self: Sized,
{
	/// The type of error returned.
	type Error;

	/// Attempt to deserialize from a vector of bytes
	fn try_deserialize(bytes: &[u8]) -> Result<Self, Self::Error>;
}
