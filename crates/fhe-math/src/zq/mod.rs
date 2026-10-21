#![warn(missing_docs, unused_imports)]
//! Ring operations for moduli up to 62 bits.
#![expect(
    clippy::indexing_slicing,
    reason = "validated arithmetic and test indices"
)]

pub(crate) mod backend;
mod modulus;
pub mod primes;
pub(crate) mod verified;

/// Structure encapsulating an integer modulus up to 62 bits.
pub type Modulus = modulus::Modulus<backend::Runtime>;

#[cfg(verus_keep_ghost)]
mod clients;
