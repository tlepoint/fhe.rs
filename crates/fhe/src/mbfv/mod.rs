// Expect indexing in multiparty BFV cryptographic operations for performance
#![expect(
    clippy::indexing_slicing,
    reason = "performance or example code relies on validated indices"
)]

//! Experimental implementation of the Multiparty BFV scheme, as described by
//! Christian Mouchet et. al. in [Multiparty Homomorphic Encryption from
//! Ring-Learning-with-Errors](https://eprint.iacr.org/2020/304.pdf).
//!
//! # Security warning
//!
//! This module is incomplete and has not been independently audited. In
//! particular, common-reference-string generation and the noise flooding
//! required by the key-switching protocols are not fully implemented. It must
//! not be used in production or to protect sensitive data.
//!
//! The module is available only with the `experimental-mbfv` Cargo feature.

mod aggregate;
mod crp;
mod public_key_gen;
mod public_key_switch;
mod relin_key_gen;
pub mod round;
mod secret_key_switch;

pub use aggregate::{Aggregate, AggregateIter};
pub use crp::CommonRandomPoly;
pub use public_key_gen::PublicKeyShare;
pub use public_key_switch::PublicKeySwitchShare;
pub use relin_key_gen::{RelinKeyGenerator, RelinKeyShare};
pub use secret_key_switch::{DecryptionShare, SecretKeySwitchShare};
