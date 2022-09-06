//! Implementation of serialization and deserialization.

use std::sync::Arc;

use super::{traits::TryConvertFrom, Context, Poly};
use crate::{proto::rq::Rq, Error};
use fhe_traits::{DeserializeWithContext, Serialize};
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

#[cfg(test)]
mod tests {
	use std::{error::Error, sync::Arc};

	use fhe_traits::{DeserializeWithContext, Serialize};

	use crate::rq::{Context, Poly, Representation};

	const Q: &[u64; 3] = &[
		4611686018282684417,
		4611686018326724609,
		4611686018309947393,
	];

	#[test]
	fn serialize() -> Result<(), Box<dyn Error>> {
		for qi in Q {
			let ctx = Arc::new(Context::new(&[*qi], 8)?);
			let p = Poly::random(&ctx, Representation::PowerBasis);
			assert_eq!(p, Poly::from_bytes(&p.to_bytes(), &ctx)?);
			let p = Poly::random(&ctx, Representation::Ntt);
			assert_eq!(p, Poly::from_bytes(&p.to_bytes(), &ctx)?);
			let p = Poly::random(&ctx, Representation::NttShoup);
			assert_eq!(p, Poly::from_bytes(&p.to_bytes(), &ctx)?);
		}

		let ctx = Arc::new(Context::new(Q, 8)?);
		let p = Poly::random(&ctx, Representation::PowerBasis);
		assert_eq!(p, Poly::from_bytes(&p.to_bytes(), &ctx)?);
		let p = Poly::random(&ctx, Representation::Ntt);
		assert_eq!(p, Poly::from_bytes(&p.to_bytes(), &ctx)?);
		let p = Poly::random(&ctx, Representation::NttShoup);
		assert_eq!(p, Poly::from_bytes(&p.to_bytes(), &ctx)?);

		Ok(())
	}
}
