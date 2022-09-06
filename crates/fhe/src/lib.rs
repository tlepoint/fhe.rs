#![crate_name = "fhe"]
#![crate_type = "lib"]
#![warn(missing_docs, unused_imports)]
#![feature(doc_cfg, int_log, int_roundings, is_some_with, let_chains)]

//! fhers: Fully Homomorphic Encryption in Rust.

mod errors;

pub mod bfv;
pub use errors::{Error, ParametersError, Result};
