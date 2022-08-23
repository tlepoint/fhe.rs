#![warn(missing_docs, unused_imports)]

//! The Brakerski-Fan-Vercauteren homomorphic encryption scheme

mod ciphertext;
mod encoding;
mod keys;
mod parameters;
mod plaintext;
mod proto;

pub mod advanced;
pub mod traits;
pub use ciphertext::Ciphertext;
pub use encoding::Encoding;
pub use keys::{EvaluationKey, EvaluationKeyBuilder, SecretKey};
pub use parameters::{BfvParameters, BfvParametersBuilder};
pub use plaintext::Plaintext;
