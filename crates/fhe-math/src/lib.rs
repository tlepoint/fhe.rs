#![crate_name = "fhe_math"]
#![crate_type = "lib"]
#![warn(missing_docs, unused_imports)]
#![feature(bigint_helper_methods, int_log, int_roundings, is_some_with, test)]

//! Mathematical utilities for the fhe.rs library.

mod errors;
mod proto;

pub mod rns;
pub mod rq;
pub mod u256;
pub mod zq;

pub use errors::{Error, Result};

#[cfg(test)]
#[macro_use]
extern crate proptest;
