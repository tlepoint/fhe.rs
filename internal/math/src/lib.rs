#![crate_name = "math"]
#![crate_type = "lib"]
#![warn(missing_docs, unused_imports)]
#![feature(bigint_helper_methods)]
#![feature(int_roundings)]
#![feature(is_some_with)]
#![feature(test)]
#![feature(int_log)]

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