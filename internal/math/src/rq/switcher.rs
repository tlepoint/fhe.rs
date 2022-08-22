#![warn(missing_docs, unused_imports)]

//! Polynomial modulus switcher.

use super::{scaler::Scaler, Context, Poly};
use crate::{rns::ScalingFactor, Result};
use std::sync::Arc;

/// Context switcher.
#[derive(Default, Debug, Clone, PartialEq, Eq)]
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
	pub(crate) fn switch(&self, p: &Poly) -> Result<Poly> {
		self.scaler.scale(p)
	}
}
