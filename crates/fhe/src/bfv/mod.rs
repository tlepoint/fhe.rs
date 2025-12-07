#![warn(missing_docs, unused_imports)]
// Allow indexing in BFV cryptographic operations for performance
#![allow(clippy::indexing_slicing)]

//! The Brakerski-Fan-Vercauteren homomorphic encryption scheme

mod ciphertext;
mod context;
mod encoding;
mod keys;
mod ops;
mod parameters;
mod plaintext;
mod plaintext_vec;
mod rgsw_ciphertext;

pub mod traits;
pub use ciphertext::Ciphertext;
pub use context::{CipherPlainContext, ContextLevel};
pub use encoding::Encoding;
pub(crate) use keys::KeySwitchingKey;
pub use keys::{EvaluationKey, EvaluationKeyBuilder, PublicKey, RelinearizationKey, SecretKey};
pub use ops::{dot_product_scalar, Multiplicator};
pub use parameters::{BfvParameters, BfvParametersBuilder};
pub use plaintext::Plaintext;
pub use plaintext_vec::PlaintextVec;
pub use rgsw_ciphertext::RGSWCiphertext;
