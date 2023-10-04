//! Protobuf for the `fhe` crate.

/// Protobuf for the BFV encryption scheme.
pub mod bfv {
    #![allow(missing_docs)]
    include!(concat!(env!("OUT_DIR"), "/fhers.bfv.rs"));
}
