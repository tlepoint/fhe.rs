//! Protobuf for the `fhe-math` crate.

#![allow(missing_docs)]

/// Protobuf for polynomials.
pub mod rq {
    include!(concat!(env!("OUT_DIR"), "/fhers.rq.rs"));
}
