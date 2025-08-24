//! Context management for the BFV encryption scheme

use fhe_math::rq::{Context, Poly};
use std::sync::Arc;

/// Stores pre-computed values relating a ciphertext and plaintext context pair.
///
/// This struct bridges the gap between ciphertext and plaintext contexts,
/// providing pre-computed values needed for efficient plaintext operations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CipherPlainContext {
    /// Scaling polynomial for the plaintext
    pub delta: Poly,

    /// Q modulo the plaintext modulus
    pub q_mod_t: u64,

    /// Threshold for centered reduction (plaintext_modulus + 1) / 2
    pub plain_threshold: u64,

    /// Associated plaintext context
    pub plaintext_context: Arc<Context>,

    /// Associated ciphertext context
    pub ciphertext_context: Arc<Context>,

    /// Next context in chain (for modulus switching)
    pub next: Option<Arc<CipherPlainContext>>,
}

impl CipherPlainContext {
    /// Creates a new CipherPlainContext from the given plaintext and ciphertext
    /// contexts.
    pub fn new_arc(
        plaintext_context: &Arc<Context>,
        ciphertext_context: &Arc<Context>,
        delta: Poly,
        q_mod_t: u64,
        plain_threshold: u64,
        next: Option<Arc<CipherPlainContext>>,
    ) -> Arc<Self> {
        Arc::new(CipherPlainContext {
            delta,
            q_mod_t,
            plain_threshold,
            plaintext_context: plaintext_context.clone(),
            ciphertext_context: ciphertext_context.clone(),
            next,
        })
    }
}
