#![warn(missing_docs, unused_imports)]

//! Traits associated with polynomials.

use super::{Context, Representation};
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
	/// Attempt to convert the `value` into a polynomial with a specific context
	/// and under a specific representation. The representation may optional and
	/// be specified as `None`; this is useful for example when converting from
	/// a value that encodes the representation (e.g., serialization, protobuf,
	/// etc.).
	fn try_convert_from<R>(
		value: T,
		ctx: &Arc<Context>,
		variable_time: bool,
		representation: R,
	) -> Result<Self>
	where
		R: Into<Option<Representation>>;
}

/// Unsigned trait.
pub trait Unsigned
where
	Self: Into<u64> + Copy,
{
}
