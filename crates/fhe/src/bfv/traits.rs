//! Traits used for the BFV homomorphic encryption scheme.

use crate::bfv::BfvParameters;
use crate::Result;
use std::sync::Arc;

/// Conversions.
///
/// We unfortunately cannot use the `TryFrom` trait from std::convert because we
/// need to specify additional parameters, and if we try to redefine a `TryFrom`
/// trait here, we need to fully specify the trait when we use it because of the
/// blanket implementation <https://github.com/rust-lang/rust/issues/50133#issuecomment-488512355>.
pub trait TryConvertFrom<T>
where
	Self: Sized,
{
	/// Attempt to convert the `value` with a specific parameter.
	fn try_convert_from(value: T, par: &Arc<BfvParameters>) -> Result<Self>;
}
