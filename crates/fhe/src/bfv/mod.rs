#![warn(missing_docs, unused_imports)]

//! The Brakerski-Fan-Vercauteren homomorphic encryption scheme

mod ciphertext;
mod context;
mod context_chain;
mod encoding;
mod keys;
mod level_manager;
mod ops;
mod parameters;
mod plaintext;
mod plaintext_vec;
mod rgsw_ciphertext;

pub mod traits;
pub use ciphertext::Ciphertext;
pub use context::CipherPlainContext;
pub use context_chain::ContextLevel;
pub use encoding::Encoding;
pub(crate) use keys::KeySwitchingKey;
pub use keys::{EvaluationKey, EvaluationKeyBuilder, PublicKey, RelinearizationKey, SecretKey};
pub use level_manager::LevelManager;
pub use ops::{dot_product_scalar, Multiplicator};
pub use parameters::{BfvParameters, BfvParametersBuilder};
pub use plaintext::Plaintext;
pub use plaintext_vec::PlaintextVec;
pub use rgsw_ciphertext::RGSWCiphertext;
