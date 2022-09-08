#![crate_name = "fhe_math"]
#![crate_type = "lib"]
#![warn(missing_docs, unused_imports)]

//! Mathematical utilities for the fhe.rs library.

mod errors;
mod proto;

pub mod rns;
pub mod rq;
pub mod zq;

pub use errors::{Error, Result};

#[cfg(test)]
#[macro_use]
extern crate proptest;
