//! Traits used for the BFV homomorphic encryption scheme.

use crate::ciphertext::Ciphertext;
use crate::parameters::BfvParameters;
use crate::plaintext::{Encoding, Plaintext};
use std::rc::Rc;

/// Encode to create plaintexts.
pub trait Encoder<T>
where
	Self: Sized,
{
	/// The type of errors.
	type Error;

	/// Attempt to encode the value within the specified parameters.
	fn try_encode(
		value: T,
		encoding: Encoding,
		par: &Rc<BfvParameters>,
	) -> Result<Self, Self::Error>;
}

/// Decode plaintext values.
pub trait Decoder
where
	Self: Sized,
{
	/// The type of errors.
	type Error;

	/// Attempt to encode the value within the specified parameters.
	fn try_decode(a: &Plaintext, encoding: Encoding) -> Result<Self, Self::Error>;
}

/// Encrypt a plaintext into a ciphertext
pub trait Encryptor {
	/// The type of errors.
	type Error;

	/// Encrypt a plaintext.
	fn encrypt(&self, plaintext: &Plaintext) -> Result<Ciphertext, Self::Error>;
}

/// Decrypt a ciphertext into a plaintext
pub trait Decryptor {
	/// The type of errors.
	type Error;

	/// Decrypt a ciphertext
	fn decrypt(&self, ciphertext: &Ciphertext) -> Result<Plaintext, Self::Error>;
}
