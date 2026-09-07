//! Implementation of serialization and deserialization.

use std::sync::Arc;

use super::{Context, Poly, RepresentationTag, wire::FromProto};
use crate::{Error, error::PolynomialSerializationError, proto::rq::Rq};

use prost::Message;

impl<R: RepresentationTag> Poly<R> {
    /// Serialize this polynomial using the existing protobuf wire format.
    #[must_use]
    pub fn to_bytes(&self) -> Vec<u8> {
        Rq::from(self).encode_to_vec()
    }
}

macro_rules! impl_from_bytes {
    ($($representation:ty),+ $(,)?) => {
        $(
        impl Poly<$representation> {
            /// Import validated polynomial bytes bound to `ctx`.
            /// Wire data never grants permission for variable-time computation.
            pub fn from_bytes(bytes: &[u8], ctx: &Arc<Context>) -> Result<Self, Error> {
                Self::from_bytes_with_limits(bytes, ctx, &crate::DecodeLimits::default())
            }

            /// Import with explicit bounds on wire and expanded residue storage.
            pub fn from_bytes_with_limits(bytes: &[u8], ctx: &Arc<Context>, limits: &crate::DecodeLimits) -> Result<Self, Error> {
                limits.check_context(bytes.len(), ctx.degree, ctx.moduli().len())?;
                limits.check_polynomials(1, ctx.degree, ctx.moduli().len())?;
                let rq: Rq = Message::decode(bytes).map_err(|source| PolynomialSerializationError::Decode { source })?;
                Self::from_proto(&rq, ctx)
            }
        }
        )+
    };
}
impl_from_bytes!(super::PowerBasis, super::Ntt, super::NttShoup);

#[cfg(test)]
mod tests {
    use std::{error::Error as StdError, sync::Arc};

    use rand::rng;

    use crate::rq::{Context, Ntt, NttShoup, Poly, PowerBasis, wire::FromProto};
    use crate::{
        Error,
        error::PolynomialSerializationError,
        proto::rq::{Representation as RepresentationProto, Rq},
    };
    use prost::Message;

    const Q: &[u64; 3] = &[
        4611686018282684417,
        4611686018326724609,
        4611686018309947393,
    ];

    #[test]
    fn serialization_preserves_legacy_payload_and_source() -> Result<(), Box<dyn StdError>> {
        let ctx = Context::new_arc(Q, 16)?;
        let pb = Poly::<PowerBasis>::random_from_seed(&ctx, [11; 32]);
        let ntt = pb.clone().into_ntt();
        let shoup = pb.clone().into_ntt_shoup();
        let payload: Vec<_> = pb
            .coefficients
            .outer_iter()
            .zip(ctx.q.iter())
            .flat_map(|(row, modulus)| modulus.serialize_vec(row.as_slice().unwrap()).unwrap())
            .collect();
        for (bytes, representation) in [
            (pb.to_bytes(), RepresentationProto::Powerbasis),
            (ntt.to_bytes(), RepresentationProto::Ntt),
            (shoup.to_bytes(), RepresentationProto::Nttshoup),
        ] {
            let expected = Rq {
                coefficients: payload.clone(),
                degree: 16,
                representation: representation as i32,
                allow_variable_time: false,
            };
            assert_eq!(bytes, expected.encode_to_vec());
        }
        assert_eq!(Poly::<PowerBasis>::from_bytes(&pb.to_bytes(), &ctx)?, pb);
        assert_eq!(Poly::<Ntt>::from_bytes(&ntt.to_bytes(), &ctx)?, ntt);
        assert_eq!(
            Poly::<NttShoup>::from_bytes(&shoup.to_bytes(), &ctx)?,
            shoup
        );
        assert_eq!(shoup, pb.into_ntt_shoup());
        Ok(())
    }

    #[test]
    fn serialize() -> Result<(), Box<dyn StdError>> {
        let mut rng = rng();

        for qi in Q {
            let ctx = Arc::new(Context::new(&[*qi], 16)?);
            let p = Poly::<PowerBasis>::random(&ctx, &mut rng);
            assert_eq!(p, Poly::<PowerBasis>::from_bytes(&p.to_bytes(), &ctx)?);
            let p = Poly::<Ntt>::random(&ctx, &mut rng);
            assert_eq!(p, Poly::<Ntt>::from_bytes(&p.to_bytes(), &ctx)?);
            let p = Poly::<NttShoup>::random(&ctx, &mut rng);
            assert_eq!(p, Poly::<NttShoup>::from_bytes(&p.to_bytes(), &ctx)?);
        }

        let ctx = Arc::new(Context::new(Q, 16)?);
        let p = Poly::<PowerBasis>::random(&ctx, &mut rng);
        assert_eq!(p, Poly::<PowerBasis>::from_bytes(&p.to_bytes(), &ctx)?);
        let p = Poly::<Ntt>::random(&ctx, &mut rng);
        assert_eq!(p, Poly::<Ntt>::from_bytes(&p.to_bytes(), &ctx)?);
        let p = Poly::<NttShoup>::random(&ctx, &mut rng);
        assert_eq!(p, Poly::<NttShoup>::from_bytes(&p.to_bytes(), &ctx)?);

        Ok(())
    }

    #[test]
    fn deserialize_unknown_representation_rejected() -> Result<(), Box<dyn StdError>> {
        let mut rng = rng();
        let ctx = Arc::new(Context::new(Q, 16)?);
        let p = Poly::<PowerBasis>::random(&ctx, &mut rng);
        let mut proto = Rq::from(&p);
        proto.representation = RepresentationProto::Unknown as i32;
        let bytes = proto.encode_to_vec();
        let err = Poly::<PowerBasis>::from_bytes(&bytes, &ctx).unwrap_err();
        assert_eq!(
            err,
            Error::PolynomialSerialization(PolynomialSerializationError::UnknownRepresentation)
        );
        Ok(())
    }

    #[test]
    fn deserialize_invalid_degree_rejected() -> Result<(), Box<dyn StdError>> {
        let mut rng = rng();
        let ctx = Arc::new(Context::new(Q, 16)?);
        let p = Poly::<PowerBasis>::random(&ctx, &mut rng);
        let mut proto = Rq::from(&p);
        proto.degree = 6;
        let bytes = proto.encode_to_vec();
        let err = Poly::<PowerBasis>::from_bytes(&bytes, &ctx).unwrap_err();
        assert_eq!(
            err,
            Error::PolynomialSerialization(PolynomialSerializationError::InvalidDegree {
                degree: 6
            })
        );
        Ok(())
    }

    #[test]
    fn deserialize_invalid_coefficients_rejected() -> Result<(), Box<dyn StdError>> {
        let mut rng = rng();
        let ctx = Arc::new(Context::new(Q, 16)?);
        let p = Poly::<PowerBasis>::random(&ctx, &mut rng);
        let mut proto = Rq::from(&p);
        proto.coefficients.clear();
        let bytes = proto.encode_to_vec();
        let err = Poly::<PowerBasis>::from_bytes(&bytes, &ctx).unwrap_err();
        assert!(matches!(
            err,
            Error::PolynomialSerialization(PolynomialSerializationError::InvalidCoefficientCount {
                actual: 0,
                expected: _
            })
        ));
        Ok(())
    }

    #[test]
    fn deserialize_representation_mismatch_rejected() -> Result<(), Box<dyn StdError>> {
        let mut rng = rng();
        let ctx = Arc::new(Context::new(Q, 16)?);
        let p = Poly::<Ntt>::random(&ctx, &mut rng);
        let proto = Rq::from(&p);
        let err = Poly::<PowerBasis>::from_proto(&proto, &ctx).unwrap_err();
        assert_eq!(
            err,
            Error::PolynomialSerialization(PolynomialSerializationError::RepresentationMismatch {
                found: crate::rq::Representation::Ntt,
                expected: crate::rq::Representation::PowerBasis,
            })
        );
        Ok(())
    }

    #[test]
    fn deserialize_variable_time_flag_is_ignored() -> Result<(), Box<dyn StdError>> {
        let mut rng = rng();
        let ctx = Arc::new(Context::new(Q, 16)?);
        let p = Poly::<PowerBasis>::random(&ctx, &mut rng);
        let mut proto = Rq::from(&p);
        proto.allow_variable_time = true;
        let bytes = proto.encode_to_vec();
        let decoded = Poly::<PowerBasis>::from_bytes(&bytes, &ctx)?;
        assert!(!decoded.allow_variable_time_computations);
        Ok(())
    }
}
