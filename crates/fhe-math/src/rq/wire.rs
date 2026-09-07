#![warn(missing_docs, unused_imports)]

//! Traits associated with polynomials.

use super::Context;
use crate::Result;
use fhe_util::VariableTime;
use std::sync::Arc;

/// Validated contextual decoding of private protobuf messages.
pub(crate) trait FromProto<T>
where
    Self: Sized,
{
    /// Attempt to convert the `value` into a polynomial with a specific
    /// context. Callers select the target representation via the `Self`
    /// type.
    fn from_proto(value: T, ctx: &Arc<Context>) -> Result<Self> {
        Self::from_proto_with_timing(value, ctx, None)
    }

    /// Convert with optional permission for variable-time processing of public
    /// data.
    ///
    /// `None` keeps variable-time algorithms disabled. A supplied token asserts
    /// that the input is public; it does not establish a constant-time
    /// guarantee for arbitrary big-integer operations on the default path.
    fn from_proto_with_timing(
        value: T,
        ctx: &Arc<Context>,
        variable_time: Option<VariableTime>,
    ) -> Result<Self>;
}
