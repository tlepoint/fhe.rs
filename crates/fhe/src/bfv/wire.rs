//! Traits used for the BFV homomorphic encryption scheme.

use crate::Result;
use crate::bfv::Parameters;

/// Internal conversion of validated protobuf DTOs with explicit parameter
/// binding.
pub(crate) trait FromProto<T>
where
    Self: Sized,
{
    /// Attempt to convert the `value` with a specific parameter.
    fn from_proto(value: T, par: &Parameters) -> Result<Self>;
}
