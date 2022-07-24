#![crate_name = "crrust_math"]
#![crate_type = "lib"]
#![warn(missing_docs, unused_imports)]
#![feature(bigint_helper_methods)]
#![feature(int_roundings)]

//! Mathematical utilities for the crrust library.

pub mod rns;
pub mod zq;
pub mod u256;

#[cfg(test)]
#[macro_use]
extern crate proptest;
