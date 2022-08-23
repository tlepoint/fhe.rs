//! Advanced structs and functions for the BFV encryption scheme.

mod leveled_evaluation_key;
mod operations;
mod plaintext_vec;

pub use leveled_evaluation_key::{LeveledEvaluationKey, LeveledEvaluationKeyBuilder};
pub use operations::{dot_product_scalar, mul, mul2};
pub use plaintext_vec::PlaintextVec;
