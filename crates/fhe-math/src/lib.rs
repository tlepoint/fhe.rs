#![crate_name = "fhe_math"]
#![crate_type = "lib"]

//! Mathematical utilities for the fhe.rs library.

/// Structured errors and classifications for this crate.
pub mod error;
mod proto;

pub mod ntt;
pub mod rns;
pub mod rq;
pub mod zq;

pub use error::{Error, Result};

#[cfg(test)]
#[macro_use]
extern crate proptest;

/// Explicit permissions for public variable-time operations and diagnostics.
pub use fhe_util::{PublicData, SecretDependentDiagnostics, VariableTime};

/// Resource limits for protobuf imports.
pub use fhe_util::DecodeLimits;
