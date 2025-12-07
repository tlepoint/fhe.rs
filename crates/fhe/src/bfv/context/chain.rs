use once_cell::sync::OnceCell;
use std::sync::{Arc, Weak};

use fhe_math::{
    rns::ScalingFactor,
    rq::{scaler::Scaler, Context},
};

use crate::bfv::{context::CipherPlainContext, parameters::MultiplicationParameters};

/// A context in the modulus switching chain
#[derive(Debug, Clone)]
pub struct ContextLevel {
    /// The polynomial context at this level
    pub poly_context: Arc<Context>,
    /// Bridge to plaintext operations
    pub(crate) cipher_plain_context: Arc<CipherPlainContext>,
    /// Level number (0 = highest, increases as we switch down)
    pub(crate) level: usize,
    /// Total number of moduli at this level
    pub(crate) num_moduli: usize,
    /// Next level in the chain (fewer moduli)
    pub next: OnceCell<Arc<ContextLevel>>,
    /// Previous level in the chain (more moduli)
    pub(crate) prev: OnceCell<Weak<ContextLevel>>,
    /// Modulus switching scaler to next level
    pub(crate) down_scaler: OnceCell<Arc<Scaler>>,
    /// Modulus switching scaler from previous level
    pub(crate) up_scaler: OnceCell<Arc<Scaler>>,
    /// Parameters required for ciphertext-ciphertext multiplication at this
    /// level
    pub(crate) mul_params: OnceCell<MultiplicationParameters>,
}

impl PartialEq for ContextLevel {
    fn eq(&self, other: &Self) -> bool {
        let Self {
            poly_context,
            cipher_plain_context,
            level,
            num_moduli,
            next: _,
            prev: _,
            down_scaler: _,
            up_scaler: _,
            mul_params: _,
        } = self;
        let Self {
            poly_context: other_poly_context,
            cipher_plain_context: other_cipher_plain_context,
            level: other_level,
            num_moduli: other_num_moduli,
            next: _,
            prev: _,
            down_scaler: _,
            up_scaler: _,
            mul_params: _,
        } = other;

        // OnceCell fields are lazily computed caching fields, not part of equality.
        level == other_level
            && num_moduli == other_num_moduli
            && poly_context == other_poly_context
            && cipher_plain_context == other_cipher_plain_context
    }
}

impl Eq for ContextLevel {}

impl ContextLevel {
    /// Create a new context level
    #[must_use]
    pub fn new(
        poly_context: Arc<Context>,
        cipher_plain_context: Arc<CipherPlainContext>,
        level: usize,
    ) -> Self {
        Self {
            num_moduli: poly_context.moduli().len(),
            poly_context,
            cipher_plain_context,
            level,
            next: OnceCell::new(),
            prev: OnceCell::new(),
            down_scaler: OnceCell::new(),
            up_scaler: OnceCell::new(),
            mul_params: OnceCell::new(),
        }
    }

    /// Chain two levels together
    pub fn chain(prev: &Arc<Self>, next: &Arc<Self>) {
        // Create scalers for modulus switching when possible
        if let Ok(ds) = Scaler::new(&prev.poly_context, &next.poly_context, ScalingFactor::one()) {
            let _ = prev.down_scaler.set(Arc::new(ds));
        }
        if let Ok(us) = Scaler::new(&next.poly_context, &prev.poly_context, ScalingFactor::one()) {
            let _ = next.up_scaler.set(Arc::new(us));
        }
        let _ = prev.next.set(next.clone());
        let _ = next.prev.set(Arc::downgrade(prev));
    }

    /// Check if this level can switch to the next
    pub fn can_switch_down(&self) -> bool {
        self.next.get().is_some()
    }

    /// Get the maximum level reachable from this context
    pub fn max_level(&self) -> usize {
        let mut max = self.level;
        let mut current = self.next.get();
        while let Some(ctx) = current {
            max = ctx.level;
            current = ctx.next.get();
        }
        max
    }

    /// Walk the entire chain without allocating a new collection
    pub fn iter_chain(&self) -> impl Iterator<Item = Arc<ContextLevel>> {
        let head = if let Some(prev) = self.prev.get().and_then(|w| w.upgrade()) {
            let mut head = prev;
            while let Some(p) = head.prev.get().and_then(|w| w.upgrade()) {
                head = p;
            }
            head
        } else {
            Arc::new(self.clone())
        };

        std::iter::successors(Some(head), |ctx| ctx.next.get().cloned())
    }

    /// Access multiplication parameters for this level
    pub(crate) fn mul_params(&self) -> &MultiplicationParameters {
        self.mul_params
            .get()
            .expect("multiplication parameters not set")
    }
}

#[cfg(test)]
mod tests {
    use crate::bfv::BfvParametersBuilder;

    #[test]
    fn chain_basics() {
        let params = BfvParametersBuilder::new()
            .set_degree(16)
            .set_plaintext_modulus(1153)
            .set_moduli_sizes(&[50, 50])
            .build()
            .unwrap();

        let head = params.context_chain();
        assert!(head.can_switch_down());
        let next = head.next.get().unwrap();
        assert!(!next.can_switch_down());
        assert_eq!(head.max_level(), 1);

        let chain: Vec<_> = head.iter_chain().collect();
        assert_eq!(chain.len(), 2);
    }
}
