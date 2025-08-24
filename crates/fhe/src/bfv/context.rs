//! Unified context management for the BFV encryption scheme
//!
//! This module provides a unified context management system that consolidates
//! the scattered context management across BfvParameters, multiple Context
//! types, and encoding contexts. The design is inspired by aspire-bfv's
//! ParamContext approach.

use crate::{Error, Result};
use fhe_math::{
    ntt::NttOperator,
    rns::ScalingFactor,
    rq::{traits::TryConvertFrom, Context, Poly, Representation},
    zq::Modulus,
};
use num_traits::ToPrimitive;
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

/// Unified context manager for the BFV encryption scheme.
///
/// This struct consolidates all context management into a single, unified
/// system that provides efficient access to all the contexts and pre-computed
/// values needed for BFV operations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FheContext {
    /// The polynomial degree
    pub polynomial_degree: usize,

    /// The plaintext modulus
    pub plaintext_modulus: u64,

    /// The ciphertext moduli
    pub ciphertext_moduli: Box<[u64]>,

    /// Ciphertext contexts for each level (modulus switching chain)
    pub ciphertext_contexts: Vec<Arc<Context>>,

    /// Plaintext context
    pub plaintext_context: Arc<Context>,

    /// Bridge contexts between ciphertext and plaintext for each level
    pub cipher_plain_contexts: Vec<Arc<CipherPlainContext>>,

    /// NTT operator for SIMD plaintext operations, if possible
    pub ntt_operator: Option<Arc<NttOperator>>,

    /// Pre-computed scaling factors for each level
    pub scaling_factors: Vec<ScalingFactor>,

    /// Plaintext modulus as a Modulus type
    pub plaintext: Modulus,
}

impl FheContext {
    /// Creates a new FheContext with the given parameters.
    ///
    /// This performs all necessary pre-computations and validations to create
    /// a complete context system for BFV operations.
    pub fn new(
        polynomial_degree: usize,
        plaintext_modulus: u64,
        ciphertext_moduli: &[u64],
        _variance: usize,
    ) -> Result<Self> {
        // Validate basic parameters
        if polynomial_degree < 8 || !polynomial_degree.is_power_of_two() {
            return Err(Error::DefaultError(format!(
                "Invalid polynomial degree: {polynomial_degree}. Must be a power of 2 >= 8"
            )));
        }

        // Create plaintext modulus
        let plaintext = Modulus::new(plaintext_modulus)
            .map_err(|e| Error::DefaultError(format!("Invalid plaintext modulus: {e}")))?;

        // Create plaintext context using the first ciphertext modulus
        // This is consistent with the original implementation
        let plaintext_context = Context::new_arc(&ciphertext_moduli[..1], polynomial_degree)?;

        // Create ciphertext contexts for modulus switching chain
        let mut ciphertext_contexts = Vec::with_capacity(ciphertext_moduli.len());
        for i in 0..ciphertext_moduli.len() {
            let level_moduli = &ciphertext_moduli[..ciphertext_moduli.len() - i];
            let ctx = Context::new_arc(level_moduli, polynomial_degree)?;
            ciphertext_contexts.push(ctx);
        }

        // Create NTT operator for SIMD operations if possible
        let ntt_operator = NttOperator::new(&plaintext, polynomial_degree).map(Arc::new);

        // Create cipher-plain bridge contexts
        let mut cipher_plain_contexts = Vec::with_capacity(ciphertext_contexts.len());
        let mut next_context: Option<Arc<CipherPlainContext>> = None;

        // Build contexts in reverse order to establish the chain
        for (i, cipher_ctx) in ciphertext_contexts.iter().enumerate().rev() {
            // Compute delta (scaling polynomial)
            let level_moduli = &ciphertext_moduli[..ciphertext_moduli.len() - i];
            let mut delta_rests = vec![];
            for m in level_moduli {
                let q = Modulus::new(*m)?;
                delta_rests.push(q.inv(q.neg(plaintext_modulus)).unwrap())
            }

            // Use RnsContext to lift the delta values and create the scaling polynomial
            let rns = fhe_math::rns::RnsContext::new(level_moduli)?;
            let mut delta = Poly::try_convert_from(
                &[rns.lift((&delta_rests).into())],
                cipher_ctx,
                true,
                Representation::PowerBasis,
            )?;
            delta.change_representation(Representation::NttShoup);

            // Compute q_mod_t
            let q_mod_t = (rns.modulus() % plaintext_modulus).to_u64().unwrap();

            // Compute plain_threshold
            let plain_threshold = plaintext_modulus.div_ceil(2);

            let cipher_plain_ctx = CipherPlainContext::new_arc(
                &plaintext_context,
                cipher_ctx,
                delta,
                q_mod_t,
                plain_threshold,
                next_context.clone(),
            );

            cipher_plain_contexts.push(cipher_plain_ctx.clone());
            next_context = Some(cipher_plain_ctx);
        }

        // Reverse to get correct order (level 0 first)
        cipher_plain_contexts.reverse();

        // Create scaling factors (placeholder for now)
        let scaling_factors = vec![ScalingFactor::one(); ciphertext_contexts.len()];

        Ok(FheContext {
            polynomial_degree,
            plaintext_modulus,
            ciphertext_moduli: ciphertext_moduli.into(),
            ciphertext_contexts,
            plaintext_context,
            cipher_plain_contexts,
            ntt_operator,
            scaling_factors,
            plaintext,
        })
    }

    /// Get the context at a specific level
    pub fn context_at_level(&self, level: usize) -> Result<&Arc<Context>> {
        self.ciphertext_contexts
            .get(level)
            .ok_or_else(|| Error::DefaultError(format!("Invalid level: {level}")))
    }

    /// Get cipher-plain context for operations at a specific level
    pub fn cipher_plain_context_at_level(&self, level: usize) -> Result<&Arc<CipherPlainContext>> {
        self.cipher_plain_contexts
            .get(level)
            .ok_or_else(|| Error::DefaultError(format!("Invalid level: {level}")))
    }

    /// Get the maximum level allowed by these parameters
    pub fn max_level(&self) -> usize {
        self.ciphertext_moduli.len() - 1
    }

    /// Find the level of a given context
    pub fn level_of_context(&self, ctx: &Arc<Context>) -> Result<usize> {
        if let Some(first_ctx) = self.ciphertext_contexts.first() {
            first_ctx.niterations_to(ctx).map_err(Error::MathError)
        } else {
            Err(Error::DefaultError("No contexts available".to_string()))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fhe_context_creation() {
        let context = FheContext::new(16, 1153, &[4611686018427387617, 4611686018427387329], 10);
        assert!(context.is_ok());

        let ctx = context.unwrap();
        assert_eq!(ctx.polynomial_degree, 16);
        assert_eq!(ctx.plaintext_modulus, 1153);
        assert_eq!(ctx.ciphertext_contexts.len(), 2);
        assert_eq!(ctx.cipher_plain_contexts.len(), 2);
        assert_eq!(ctx.max_level(), 1);
    }

    #[test]
    fn test_context_access() {
        let context =
            FheContext::new(16, 1153, &[4611686018427387617, 4611686018427387329], 10).unwrap();

        // Test valid level access
        assert!(context.context_at_level(0).is_ok());
        assert!(context.context_at_level(1).is_ok());
        assert!(context.cipher_plain_context_at_level(0).is_ok());
        assert!(context.cipher_plain_context_at_level(1).is_ok());

        // Test invalid level access
        assert!(context.context_at_level(2).is_err());
        assert!(context.cipher_plain_context_at_level(2).is_err());
    }

    #[test]
    fn test_invalid_parameters() {
        // Invalid polynomial degree (not power of 2)
        assert!(FheContext::new(15, 1153, &[4611686018427387617], 10).is_err());

        // Invalid polynomial degree (too small)
        assert!(FheContext::new(4, 1153, &[4611686018427387617], 10).is_err());

        // Invalid plaintext modulus (0)
        assert!(FheContext::new(16, 0, &[4611686018427387617], 10).is_err());
    }
}
