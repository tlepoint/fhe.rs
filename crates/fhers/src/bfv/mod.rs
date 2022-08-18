#![warn(missing_docs, unused_imports)]

//! The Brakerski-Fan-Vercauteren homomorphic encryption scheme

mod ciphertext;
mod encoding;
mod keys;
mod parameters;
mod parameters_switcher;
mod plaintext;
mod proto;

pub mod traits;
pub use ciphertext::{dot_product_scalar, mul, mul2, Ciphertext};
pub use encoding::Encoding;
pub use keys::{EvaluationKey, EvaluationKeyBuilder, SecretKey};
pub use parameters::{BfvParameters, BfvParametersBuilder};
pub use parameters_switcher::{BfvParametersSwitcher, ParametersSwitchable};
pub use plaintext::{Plaintext, VecPlaintext};
