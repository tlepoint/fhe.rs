#![crate_name = "bfv"]
#![crate_type = "lib"]
#![warn(missing_docs, unused_imports)]
#![feature(is_some_with)]
#![feature(option_get_or_insert_default)]
#![feature(int_log)]

//! The Brakerski-Fan-Vercauteren homomorphic encryption scheme

mod ciphertext;
mod keys;
mod parameters;
mod plaintext;

pub mod traits;
pub use ciphertext::{mul, mul2, Ciphertext};
pub use keys::{EvaluationKey, EvaluationKeyBuilder, SecretKey};
pub use parameters::{BfvParameters, BfvParametersBuilder};
pub use plaintext::{Encoding, Plaintext};
