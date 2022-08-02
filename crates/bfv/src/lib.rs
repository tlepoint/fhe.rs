#![crate_name = "bfv"]
#![crate_type = "lib"]
#![warn(missing_docs, unused_imports)]
#![feature(is_some_with)]

//! The Brakerski-Fan-Vercauteren homomorphic encryption scheme

mod ciphertext;
mod keys;
mod parameters;
mod plaintext;

pub mod traits;
pub use ciphertext::{mul, Ciphertext};
pub use keys::{RelinearizationKey, SecretKey};
pub use parameters::{BfvParameters, BfvParametersBuilder};
pub use plaintext::{Encoding, Plaintext};

#[cfg(test)]
#[macro_use]
extern crate proptest;
