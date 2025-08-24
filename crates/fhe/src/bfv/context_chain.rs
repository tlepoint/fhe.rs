use std::sync::{Arc, Weak};

use fhe_math::{
    rns::ScalingFactor,
    rq::{scaler::Scaler, Context},
};

use crate::bfv::context::CipherPlainContext;

/// A context in the modulus switching chain
#[derive(Debug, Clone)]
pub struct ContextLevel {
    /// The polynomial context at this level
    pub poly_context: Arc<Context>,
    /// Bridge to plaintext operations
    pub cipher_plain_context: Arc<CipherPlainContext>,
    /// Level number (0 = highest, increases as we switch down)
    pub level: usize,
    /// Total number of moduli at this level
    pub num_moduli: usize,
    /// Next level in the chain (fewer moduli)
    pub next: Option<Arc<ContextLevel>>,
    /// Previous level in the chain (more moduli)
    pub prev: Option<Weak<ContextLevel>>,
    /// Modulus switching scaler to next level
    pub down_scaler: Option<Arc<Scaler>>,
    /// Modulus switching scaler from previous level
    pub up_scaler: Option<Arc<Scaler>>,
}

impl PartialEq for ContextLevel {
    fn eq(&self, other: &Self) -> bool {
        self.level == other.level
            && self.num_moduli == other.num_moduli
            && self.cipher_plain_context == other.cipher_plain_context
    }
}

impl Eq for ContextLevel {}

impl ContextLevel {
    /// Create a new context level
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
            next: None,
            prev: None,
            down_scaler: None,
            up_scaler: None,
        }
    }

    /// Chain two levels together
    pub fn chain(prev: &mut Arc<Self>, next: &mut Arc<Self>) {
        // Create scalers for modulus switching when possible
        if let Ok(ds) = Scaler::new(&prev.poly_context, &next.poly_context, ScalingFactor::one()) {
            if let Some(p) = Arc::get_mut(prev) {
                p.down_scaler = Some(Arc::new(ds));
            }
        }
        if let Ok(us) = Scaler::new(&next.poly_context, &prev.poly_context, ScalingFactor::one()) {
            if let Some(n) = Arc::get_mut(next) {
                n.up_scaler = Some(Arc::new(us));
            }
        }
        if let Some(p) = Arc::get_mut(prev) {
            p.next = Some(next.clone());
        }
        if let Some(n) = Arc::get_mut(next) {
            n.prev = Some(Arc::downgrade(prev));
        }
    }

    /// Check if this level can switch to the next
    pub fn can_switch_down(&self) -> bool {
        self.next.is_some()
    }

    /// Get the maximum level reachable from this context
    pub fn max_level(&self) -> usize {
        let mut max = self.level;
        let mut current = &self.next;
        while let Some(ctx) = current {
            max = ctx.level;
            current = &ctx.next;
        }
        max
    }

    /// Walk the entire chain and collect all levels
    pub fn collect_chain(&self) -> Vec<Arc<ContextLevel>> {
        let mut chain = Vec::new();
        // Move to head
        let mut current = self.clone();
        while let Some(prev) = current.prev.as_ref().and_then(|w| w.upgrade()) {
            current = (*prev).clone();
        }
        // Collect
        let mut cur = Some(Arc::new(current));
        while let Some(node) = cur {
            chain.push(node.clone());
            cur = node.next.clone();
        }
        chain
    }
}
