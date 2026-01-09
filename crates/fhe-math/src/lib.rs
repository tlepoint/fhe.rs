#![crate_name = "fhe_math"]
#![crate_type = "lib"]

//! Mathematical utilities for the fhe.rs library.

mod errors;
mod proto;

pub mod ntt;
pub mod rns;
pub mod rq;
pub mod zq;

pub use errors::{Error, Result};

#[cfg(test)]
#[macro_use]
extern crate proptest;
