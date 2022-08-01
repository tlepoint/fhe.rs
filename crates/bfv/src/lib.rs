#![crate_name = "bfv"]
#![crate_type = "lib"]
#![warn(missing_docs, unused_imports)]
#![feature(is_some_with)]

//! The Brakerski-Fan-Vercauteren homomorphic encryption scheme

pub mod keys;
pub mod parameters;
pub mod plaintext;

#[cfg(test)]
#[macro_use]
extern crate proptest;
