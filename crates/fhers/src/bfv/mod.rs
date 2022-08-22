#![warn(missing_docs, unused_imports)]

//! The Brakerski-Fan-Vercauteren homomorphic encryption scheme

mod ciphertext;
mod encoding;
mod keys;
mod parameters;
mod plaintext;
mod proto;

pub mod traits;
pub use ciphertext::{dot_product_scalar, mul, mul2, Ciphertext};
pub use encoding::Encoding;
pub use keys::{LeveledEvaluationKey, LeveledEvaluationKeyBuilder, SecretKey};
pub use parameters::{BfvParameters, BfvParametersBuilder};
pub use plaintext::{Plaintext, PlaintextVec};
