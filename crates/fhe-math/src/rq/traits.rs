#![warn(missing_docs, unused_imports)]

//! Traits associated with polynomials.

use super::Context;
use crate::Result;
use fhe_util::VariableTime;
use std::sync::Arc;

/// Conversions to create polynomials.
///
/// The input and a context determine the output. By default variable-time
/// processing is disabled; selecting it requires explicit public-data evidence.
/// ```compile_fail
/// use fhe_math::rq::{Context, Poly, PowerBasis, traits::TryConvertFrom};
/// use std::sync::Arc;
/// fn convert(context: &Arc<Context>) {
///     Poly::<PowerBasis>::try_convert_from(&[1_u64], context, true);
/// }
/// ```
pub trait TryConvertFrom<T>
where
    Self: Sized,
{
    /// Attempt to convert the `value` into a polynomial with a specific
    /// context. Callers select the target representation via the `Self`
    /// type.
    fn try_convert_from(value: T, ctx: &Arc<Context>) -> Result<Self> {
        Self::try_convert_from_with_timing(value, ctx, None)
    }

    /// Convert with optional permission for variable-time processing of public
    /// data.
    ///
    /// `None` keeps variable-time algorithms disabled. A supplied token asserts
    /// that the input is public; it does not establish a constant-time
    /// guarantee for arbitrary big-integer operations on the default path.
    fn try_convert_from_with_timing(
        value: T,
        ctx: &Arc<Context>,
        variable_time: Option<VariableTime>,
    ) -> Result<Self>;
}
