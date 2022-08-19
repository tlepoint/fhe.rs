//! Implementation of serialization and deserialization.

use std::sync::Arc;

use super::{traits::TryConvertFrom, Context, Poly};
use crate::{proto::rq::Rq, Error};
use fhers_traits::{DeserializeWithContext, Serialize};
use protobuf::Message;

impl Serialize for Poly {
	fn to_bytes(&self) -> Vec<u8> {
		let rq = Rq::from(self);
		rq.write_to_bytes().unwrap()
	}
}

impl DeserializeWithContext for Poly {
	type Error = Error;
	type Context = Context;

	fn from_bytes(bytes: &[u8], ctx: &Arc<Context>) -> Result<Self, Self::Error> {
		let rq = Rq::parse_from_bytes(bytes).map_err(|e| Error::Serialization(e.to_string()))?;
		Poly::try_convert_from(&rq, ctx, false, None)
	}
}
