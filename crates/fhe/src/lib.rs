#![crate_name = "fhe"]
#![crate_type = "lib"]
#![doc = include_str!("../README.md")]

/// Structured errors and classifications for this crate.
pub mod error;

pub mod bfv;
#[cfg(feature = "experimental-mbfv")]
pub mod mbfv;
mod proto;
pub use error::{Error, Result};

// Test the source code included in the README.
#[macro_use]
extern crate doc_comment;
doctest!("../README.md");

/// Explicit permissions for public variable-time operations and diagnostics.
pub use fhe_util::{PublicData, SecretDependentDiagnostics, VariableTime};

/// Resource limits for protobuf imports.
pub use fhe_util::DecodeLimits;
