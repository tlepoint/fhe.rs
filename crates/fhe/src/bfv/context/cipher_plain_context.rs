use fhe_math::rq::{Context, NttShoup, Poly, scaler::Scaler};
use num_bigint::BigUint;
use std::sync::Arc;

/// Stores pre-computed values relating a ciphertext and plaintext context pair.
///
/// This struct bridges the gap between ciphertext and plaintext contexts,
/// providing pre-computed values needed for efficient plaintext operations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CipherPlainContext {
    /// Scaling polynomial for the plaintext
    pub(crate) delta: Poly<NttShoup>,

    /// Q modulo the plaintext modulus
    pub(crate) q_mod_t: BigUint,

    /// Threshold for centered reduction (plaintext_modulus + 1) / 2
    pub(crate) plain_threshold: BigUint,

    /// Scaler to map a ciphertext polynomial to the plaintext context
    pub(crate) scaler: Scaler,

    /// Associated plaintext context
    pub(crate) plaintext_context: Arc<Context>,

    /// Associated ciphertext context
    pub(crate) ciphertext_context: Arc<Context>,
}

impl CipherPlainContext {
    /// Creates a new CipherPlainContext from the given plaintext and ciphertext
    /// contexts.
    pub(crate) fn new_arc(
        plaintext_context: &Arc<Context>,
        ciphertext_context: &Arc<Context>,
        delta: Poly<NttShoup>,
        q_mod_t: BigUint,
        plain_threshold: BigUint,
        scaler: Scaler,
    ) -> Arc<Self> {
        Arc::new(CipherPlainContext {
            delta,
            q_mod_t,
            plain_threshold,
            scaler,
            plaintext_context: plaintext_context.clone(),
            ciphertext_context: ciphertext_context.clone(),
        })
    }
}
