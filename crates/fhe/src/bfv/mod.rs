#![warn(missing_docs)]
// Expect indexing in BFV cryptographic operations for performance
#![expect(
    clippy::indexing_slicing,
    reason = "performance or example code relies on validated indices"
)]

//! The Brakerski-Fan-Vercauteren homomorphic encryption scheme

mod ciphertext;
pub mod context;
mod encoding;
mod keys;
mod ops;
mod parameters;
mod plaintext;
mod plaintext_chunks;
mod rgsw_ciphertext;

mod wire;
pub use ciphertext::Ciphertext;
pub use encoding::Encoding;
#[cfg(feature = "experimental-mbfv")]
pub(crate) use keys::KeySwitchingKey;
pub use keys::{PublicKey, SecretKey};

pub use parameters::{ParameterProfile, Parameters, ParametersBuilder};
pub use plaintext::Plaintext;
mod packed_plaintext;
pub(crate) use packed_plaintext::PackedPlaintextView;

/// Evaluation keys, immutable multiplication plans, and caller-owned scratch.
/// Workspaces are independent mutable values and can be used in caller-owned
/// thread pools. This module does not start threads or retain global scratch.
pub mod evaluation {
    pub use super::keys::{
        EvaluationKey, EvaluationKeyBuilder, RelinearizationKey, RelinearizationKeyBuilder,
    };
    pub use super::ops::{
        CiphertextProductAccumulator, DotProductScalarWorkspace, MultiplicationPlan,
        MultiplicationPlanBuilder, MultiplicationScaling, PreparedMultiplicand, dot_product_scalar,
        dot_product_scalar_iter,
    };
    pub use super::rgsw_ciphertext::RgswCiphertext;
}

/// Compact, validated storage of plaintext NTT residues for repeated
/// evaluation.
pub mod packing {
    pub use super::packed_plaintext::{
        PackedPlaintext, PackedPlaintextBatch, PackedPlaintextIter, PackedPlaintextView,
    };
}
