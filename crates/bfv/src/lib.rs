#![crate_name = "bfv"]
#![crate_type = "lib"]
#![warn(missing_docs, unused_imports)]
#![feature(is_some_with)]
#![feature(option_get_or_insert_default)]
#![feature(int_log)]
#![feature(int_roundings)]

//! The Brakerski-Fan-Vercauteren homomorphic encryption scheme

mod ciphertext;
mod errors;
mod keys;
mod parameters;
mod parameters_switcher;
mod plaintext;

pub mod traits;
pub use ciphertext::{dot_product_scalar, mul, mul2, Ciphertext};
pub use errors::{Error, ParametersError, Result};
pub use keys::{EvaluationKey, EvaluationKeyBuilder, SecretKey};
pub use parameters::{BfvParameters, BfvParametersBuilder};
pub use parameters_switcher::{BfvParametersSwitcher, ParametersSwitchable};
pub use plaintext::{Encoding, Plaintext};
