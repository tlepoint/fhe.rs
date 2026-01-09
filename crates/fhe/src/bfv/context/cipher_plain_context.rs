use fhe_math::rq::{Context, Poly, scaler::Scaler};
use std::sync::Arc;

/// Stores pre-computed values relating a ciphertext and plaintext context pair.
///
/// This struct bridges the gap between ciphertext and plaintext contexts,
/// providing pre-computed values needed for efficient plaintext operations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CipherPlainContext {
    /// Scaling polynomial for the plaintext
    pub(crate) delta: Poly,

    /// Q modulo the plaintext modulus
    pub(crate) q_mod_t: u64,

    /// Threshold for centered reduction (plaintext_modulus + 1) / 2
    pub(crate) plain_threshold: u64,

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
        delta: Poly,
        q_mod_t: u64,
        plain_threshold: u64,
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
