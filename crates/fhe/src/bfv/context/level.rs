use std::sync::Arc;

use fhe_math::rq::Context;

use crate::bfv::{context::CipherPlainContext, parameters::MultiplicationParameters};

/// Precomputed data for one level of the ciphertext modulus hierarchy.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContextLevel {
    /// The polynomial context at this level.
    pub poly_context: Arc<Context>,
    /// Bridge to plaintext operations.
    pub(crate) cipher_plain_context: Arc<CipherPlainContext>,
    /// Level number (0 = highest, increases as moduli are removed).
    pub(crate) level: usize,
    /// Parameters required for ciphertext-ciphertext multiplication at this
    /// level.
    pub(crate) mul_params: MultiplicationParameters,
}

impl ContextLevel {
    pub(crate) fn new(
        poly_context: Arc<Context>,
        cipher_plain_context: Arc<CipherPlainContext>,
        level: usize,
        mul_params: MultiplicationParameters,
    ) -> Self {
        Self {
            poly_context,
            cipher_plain_context,
            level,
            mul_params,
        }
    }

    /// Return this context's level number.
    #[must_use]
    pub const fn level(&self) -> usize {
        self.level
    }

    /// Check whether another modulus can be removed from this level.
    #[must_use]
    pub fn can_switch_down(&self) -> bool {
        self.poly_context.moduli().len() > 1
    }

    /// Return the deepest level reachable from this context.
    #[must_use]
    pub fn max_level(&self) -> usize {
        self.level + self.poly_context.moduli().len() - 1
    }

    /// Access multiplication parameters for this level.
    pub(crate) const fn mul_params(&self) -> &MultiplicationParameters {
        &self.mul_params
    }
}
