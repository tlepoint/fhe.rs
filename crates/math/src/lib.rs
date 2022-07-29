#![crate_name = "math"]
#![crate_type = "lib"]
#![warn(missing_docs, unused_imports)]
#![feature(bigint_helper_methods)]
#![feature(int_roundings)]
#![feature(is_some_with)]

//! Mathematical utilities for the fhe.rs library.

mod protos;

pub mod rns;
pub mod rq;
pub mod u256;
pub mod zq;

#[cfg(test)]
#[macro_use]
extern crate proptest;
