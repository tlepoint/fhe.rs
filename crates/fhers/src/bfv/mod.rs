#![warn(missing_docs, unused_imports)]

//! The Brakerski-Fan-Vercauteren homomorphic encryption scheme

mod ciphertext;
mod encoding;
mod keys;
mod leveled;
mod ops;
mod parameters;
mod plaintext;
mod proto;

pub mod traits;
pub use ciphertext::Ciphertext;
pub use encoding::Encoding;
pub use keys::{EvaluationKey, EvaluationKeyBuilder, SecretKey};
pub use parameters::{BfvParameters, BfvParametersBuilder};
pub use plaintext::Plaintext;

#[cfg(feature = "leveled_bfv")]
pub use leveled::{LeveledEvaluationKey, LeveledEvaluationKeyBuilder, PlaintextVec};

#[cfg(feature = "optimized_ops")]
pub use ops::{dot_product_scalar, mul_relin, mul_relin_2};
