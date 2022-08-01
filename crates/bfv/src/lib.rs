#![crate_name = "bfv"]
#![crate_type = "lib"]
#![warn(missing_docs, unused_imports)]
#![feature(is_some_with)]

//! The Brakerski-Fan-Vercauteren homomorphic encryption scheme

mod ciphertext;
mod parameters;
mod plaintext;
mod secret_key;
mod traits;

pub use ciphertext::Ciphertext;
pub use parameters::{BfvParameters, BfvParametersBuilder};
pub use plaintext::{Encoding, Plaintext};
pub use secret_key::SecretKey;

#[cfg(test)]
#[macro_use]
extern crate proptest;
