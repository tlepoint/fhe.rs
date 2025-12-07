// Allow indexing in multiparty BFV cryptographic operations for performance
#![allow(clippy::indexing_slicing)]

//! The Multiparty BFV scheme, as described by Christian Mouchet et. al.
//! in [Multiparty Homomorphic Encryption from Ring-Learning-with-Errors](https://eprint.iacr.org/2020/304.pdf).

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
