#![crate_name = "fhe"]
#![crate_type = "lib"]
#![warn(missing_docs, unused_imports)]

//! fhers: Fully Homomorphic Encryption in Rust.

mod errors;

pub mod bfv;
pub use errors::{Error, ParametersError, Result};
