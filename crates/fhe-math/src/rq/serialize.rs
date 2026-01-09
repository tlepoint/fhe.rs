//! Implementation of serialization and deserialization.

use std::sync::Arc;

use super::{Context, Poly, traits::TryConvertFrom};
use crate::{Error, proto::rq::Rq};
use fhe_traits::{DeserializeWithContext, Serialize};
use prost::Message;

impl Serialize for Poly {
    fn to_bytes(&self) -> Vec<u8> {
        Rq::from(self).encode_to_vec()
    }
}

impl DeserializeWithContext for Poly {
    type Error = Error;
    type Context = Context;

    fn from_bytes(bytes: &[u8], ctx: &Arc<Context>) -> Result<Self, Self::Error> {
        let rq: Rq = Message::decode(bytes).map_err(|e| Error::Serialization(e.to_string()))?;
        Poly::try_convert_from(&rq, ctx, false, None)
    }
}

#[cfg(test)]
mod tests {
    use std::{error::Error, sync::Arc};

    use fhe_traits::{DeserializeWithContext, Serialize};
    use rand::rng;

    use crate::rq::{Context, Poly, Representation};

    const Q: &[u64; 3] = &[
        4611686018282684417,
        4611686018326724609,
        4611686018309947393,
    ];

    #[test]
    fn serialize() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();

        for qi in Q {
            let ctx = Arc::new(Context::new(&[*qi], 16)?);
            let p = Poly::random(&ctx, Representation::PowerBasis, &mut rng);
            assert_eq!(p, Poly::from_bytes(&p.to_bytes(), &ctx)?);
            let p = Poly::random(&ctx, Representation::Ntt, &mut rng);
            assert_eq!(p, Poly::from_bytes(&p.to_bytes(), &ctx)?);
            let p = Poly::random(&ctx, Representation::NttShoup, &mut rng);
            assert_eq!(p, Poly::from_bytes(&p.to_bytes(), &ctx)?);
        }

        let ctx = Arc::new(Context::new(Q, 16)?);
        let p = Poly::random(&ctx, Representation::PowerBasis, &mut rng);
        assert_eq!(p, Poly::from_bytes(&p.to_bytes(), &ctx)?);
        let p = Poly::random(&ctx, Representation::Ntt, &mut rng);
        assert_eq!(p, Poly::from_bytes(&p.to_bytes(), &ctx)?);
        let p = Poly::random(&ctx, Representation::NttShoup, &mut rng);
        assert_eq!(p, Poly::from_bytes(&p.to_bytes(), &ctx)?);

        Ok(())
    }
}
