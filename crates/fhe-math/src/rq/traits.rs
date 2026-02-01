#![warn(missing_docs, unused_imports)]

//! Traits associated with polynomials.

use super::Context;
use crate::Result;
use std::sync::Arc;

/// Conversions to create polynomials.
///
/// We unfortunately cannot use the `TryFrom` trait from std::convert because we
/// need to specify additional parameters, and if we try to redefine a `TryFrom`
/// trait here, we need to fully specify the trait when we use it because of the
/// blanket implementation <https://github.com/rust-lang/rust/issues/50133#issuecomment-488512355>.
pub trait TryConvertFrom<T>
where
    Self: Sized,
{
    /// Attempt to convert the `value` into a polynomial with a specific
    /// context. Callers select the target representation via the `Self`
    /// type.
    fn try_convert_from(value: T, ctx: &Arc<Context>, variable_time: bool) -> Result<Self>;
}
