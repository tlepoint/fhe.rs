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

    use crate::proto::rq::{Representation as RepresentationProto, Rq};
    use crate::rq::{Context, Poly, Representation, traits::TryConvertFrom};
    use prost::Message;

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

    #[test]
    fn deserialize_unknown_representation_rejected() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let ctx = Arc::new(Context::new(Q, 16)?);
        let p = Poly::random(&ctx, Representation::PowerBasis, &mut rng);
        let mut proto = Rq::from(&p);
        proto.representation = RepresentationProto::Unknown as i32;
        let bytes = proto.encode_to_vec();
        let err = Poly::from_bytes(&bytes, &ctx).unwrap_err();
        assert!(err.to_string().contains("Unknown representation"));
        Ok(())
    }

    #[test]
    fn deserialize_invalid_degree_rejected() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let ctx = Arc::new(Context::new(Q, 16)?);
        let p = Poly::random(&ctx, Representation::PowerBasis, &mut rng);
        let mut proto = Rq::from(&p);
        proto.degree = 6;
        let bytes = proto.encode_to_vec();
        let err = Poly::from_bytes(&bytes, &ctx).unwrap_err();
        assert!(err.to_string().contains("Invalid degree"));
        Ok(())
    }

    #[test]
    fn deserialize_invalid_coefficients_rejected() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let ctx = Arc::new(Context::new(Q, 16)?);
        let p = Poly::random(&ctx, Representation::PowerBasis, &mut rng);
        let mut proto = Rq::from(&p);
        proto.coefficients.clear();
        let bytes = proto.encode_to_vec();
        let err = Poly::from_bytes(&bytes, &ctx).unwrap_err();
        assert!(err.to_string().contains("Invalid coefficients"));
        Ok(())
    }

    #[test]
    fn deserialize_representation_mismatch_rejected() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let ctx = Arc::new(Context::new(Q, 16)?);
        let p = Poly::random(&ctx, Representation::Ntt, &mut rng);
        let proto = Rq::from(&p);
        let err =
            Poly::try_convert_from(&proto, &ctx, false, Representation::PowerBasis).unwrap_err();
        assert!(
            err.to_string()
                .contains("representation asked for does not match")
        );
        Ok(())
    }

    #[test]
    fn deserialize_variable_time_flag_propagates() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let ctx = Arc::new(Context::new(Q, 16)?);
        let p = Poly::random(&ctx, Representation::PowerBasis, &mut rng);
        let mut proto = Rq::from(&p);
        proto.allow_variable_time = true;
        let bytes = proto.encode_to_vec();
        let decoded = Poly::from_bytes(&bytes, &ctx)?;
        assert!(decoded.allow_variable_time_computations);
        Ok(())
    }
}
