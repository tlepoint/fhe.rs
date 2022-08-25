#![warn(missing_docs, unused_imports)]

//! The Brakerski-Fan-Vercauteren homomorphic encryption scheme

mod ciphertext;
mod encoding;
mod keys;
mod ops;
mod parameters;
mod plaintext;
mod plaintext_vec;
mod proto;

pub mod traits;
pub use ciphertext::Ciphertext;
pub use encoding::Encoding;
pub use keys::{EvaluationKey, EvaluationKeyBuilder, RelinearizationKey, SecretKey};
pub use ops::Multiplicator;
pub use parameters::{BfvParameters, BfvParametersBuilder};
pub use plaintext::Plaintext;
pub use plaintext_vec::PlaintextVec;

#[cfg(feature = "optimized_ops")]
pub use ops::dot_product_scalar;
