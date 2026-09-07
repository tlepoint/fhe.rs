#![crate_name = "fhe_math"]
#![crate_type = "lib"]

//! Mathematical utilities for the fhe.rs library.

mod errors;
mod proto;

pub mod ntt;
pub mod rns;
pub mod rq;
pub mod zq;

pub use errors::{Error, PolynomialSerializationError, Result};

#[cfg(test)]
#[macro_use]
extern crate proptest;

/// Explicit permissions for public variable-time operations and diagnostics.
pub use fhe_util::{PublicData, SecretDependentDiagnostics, VariableTime};
