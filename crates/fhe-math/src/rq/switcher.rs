#![warn(missing_docs, unused_imports)]

//! Polynomial modulus switcher.

use super::{Context, Poly, ScaleRepresentation, scaler::Scaler};
use crate::{Result, rns::ScalingFactor};
use std::sync::Arc;

/// Context switcher.
///
/// Construct this value explicitly; there is no uninitialized default.
/// ```compile_fail
/// let invalid = fhe_math::rq::switcher::Switcher::default();
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Switcher {
    pub(crate) scaler: Scaler,
}

impl Switcher {
    /// Create a switcher from a context `from` to a context `to`.
    pub fn new(from: &Arc<Context>, to: &Arc<Context>) -> Result<Self> {
        Ok(Self {
            scaler: Scaler::new(from, to, ScalingFactor::new(to.modulus(), from.modulus()))?,
        })
    }

    /// Switch a polynomial.
    pub(crate) fn switch<R: ScaleRepresentation>(&self, p: &Poly<R>) -> Result<Poly<R>> {
        self.scaler.scale(p)
    }
}
