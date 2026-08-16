#![crate_name = "fhe"]
#![crate_type = "lib"]
#![doc = include_str!("../README.md")]

mod errors;

pub mod bfv;
#[cfg(feature = "experimental-mbfv")]
pub mod mbfv;
pub mod proto;
pub use errors::{
    CiphertextError, CiphertextOperation, DotProductError, EncodingError, Error,
    EvaluationKeyComponent, EvaluationKeyError, EvaluationOperation, MultipartyError,
    ParameterSource, ParametersError, PlaintextError, Result, SerializationError, SerializedField,
    SerializedObject, SerializedPolynomialComponent,
};

// Test the source code included in the README.
#[macro_use]
extern crate doc_comment;
doctest!("../README.md");
