#![crate_name = "fhers"]
#![crate_type = "lib"]
#![warn(missing_docs, unused_imports)]
#![feature(is_some_with)]
#![feature(option_get_or_insert_default)]
#![feature(int_log)]
#![feature(int_roundings)]
#![feature(let_chains)]

//! fhers: Fully Homomorphic Encryption in Rust.

mod errors;

pub mod bfv;
pub use errors::{Error, ParametersError, Result};
