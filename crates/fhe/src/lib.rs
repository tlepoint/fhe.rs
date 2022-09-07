#![crate_name = "fhe"]
#![crate_type = "lib"]
#![warn(missing_docs, unused_imports)]
#![feature(is_some_with, let_chains)]

//! fhers: Fully Homomorphic Encryption in Rust.

mod errors;

pub mod bfv;
pub use errors::{Error, ParametersError, Result};
