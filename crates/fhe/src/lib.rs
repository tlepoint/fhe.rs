#![crate_name = "fhe"]
#![crate_type = "lib"]
#![warn(missing_docs, unused_imports)]
#![doc = include_str!("../README.md")]

mod errors;

pub mod bfv;
pub use errors::{Error, ParametersError, Result};
