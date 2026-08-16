use std::{cmp::min, ops::Deref, sync::Arc};

use fhe_math::rq::{Context, Ntt, Poly, PowerBasis, traits::TryConvertFrom};
use fhe_traits::{FheEncoder, FheEncoderVariableTime, FheParametrized, FhePlaintext, VariableTime};
use num_bigint::BigUint;
use num_traits::{ToPrimitive, Zero};
use zeroize_derive::{Zeroize, ZeroizeOnDrop};

use crate::{
    Error, Result,
    bfv::{BfvParameters, Encoding, Plaintext},
};

use super::encoding::EncodingEnum;

/// A wrapper around a vector of plaintext which implements the [`FhePlaintext`]
/// trait, and therefore can be encoded to / decoded from.
#[derive(Zeroize, ZeroizeOnDrop)]
pub struct PlaintextVec(Vec<Plaintext>);

impl Deref for PlaintextVec {
    type Target = [Plaintext];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl FhePlaintext for PlaintextVec {
    type Encoding = Encoding;
}

impl FheParametrized for PlaintextVec {
    type Parameters = BfvParameters;
}

impl PlaintextVec {
    fn try_encode_with<T>(
        value: &[T],
        encoding: Encoding,
        par: &Arc<BfvParameters>,
        mut encode_chunk: impl FnMut(
            &[T],
            &Encoding,
            &Arc<BfvParameters>,
            &Arc<Context>,
        ) -> Result<Poly<Ntt>>,
    ) -> Result<Self> {
        if encoding.encoding == EncodingEnum::Simd && par.ntt_operator.is_none() {
            return Err(crate::EncodingError::SimdUnavailable.into());
        }

        let ctx = par.context_at_level(encoding.level)?;
        let num_plaintexts = value.len().div_ceil(par.degree()).max(1);
        let plaintexts = (0..num_plaintexts)
            .map(|index| {
                let start = index * par.degree();
                let end = min(value.len(), start + par.degree());
                let poly_ntt = encode_chunk(&value[start..end], &encoding, par, ctx)?;
                Ok(Plaintext {
                    par: par.clone(),
                    encoding: Some(encoding.clone()),
                    poly_ntt,
                })
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Self(plaintexts))
    }

    fn encode_u64_chunk(
        value: &[u64],
        encoding: &Encoding,
        par: &Arc<BfvParameters>,
        ctx: &Arc<Context>,
        variable_time: Option<VariableTime>,
    ) -> Result<Poly<Ntt>> {
        let mut coefficients = vec![0u64; par.degree()];
        match encoding.encoding {
            EncodingEnum::Poly => coefficients[..value.len()].copy_from_slice(value),
            EncodingEnum::Simd => {
                for (index, &coefficient) in value.iter().enumerate() {
                    coefficients[par.matrix_reps_index_map[index]] = coefficient;
                }
                let ntt_operator = par
                    .ntt_operator
                    .as_ref()
                    .ok_or(crate::PlaintextError::NttOperatorUnavailable)?;
                if variable_time.is_some() {
                    unsafe { ntt_operator.backward_vt(coefficients.as_mut_ptr()) };
                } else {
                    ntt_operator.backward(&mut coefficients);
                }
            }
        }

        let poly = if let Some(variable_time) = variable_time {
            Poly::<PowerBasis>::try_convert_from_public(&coefficients, ctx, variable_time)?
        } else {
            Poly::<PowerBasis>::try_convert_from(&coefficients, ctx, false)?
        };
        Ok(poly.into_ntt())
    }

    fn encode_biguint_chunk(
        value: &[BigUint],
        encoding: &Encoding,
        par: &Arc<BfvParameters>,
        ctx: &Arc<Context>,
    ) -> Result<Poly<Ntt>> {
        match encoding.encoding {
            EncodingEnum::Poly => {
                let mut coefficients = vec![BigUint::zero(); par.degree()];
                coefficients[..value.len()].clone_from_slice(value);
                Ok(
                    Poly::<PowerBasis>::try_convert_from(coefficients.as_slice(), ctx, false)?
                        .into_ntt(),
                )
            }
            EncodingEnum::Simd => {
                let values = value
                    .iter()
                    .map(|coefficient| {
                        coefficient
                            .to_u64()
                            .ok_or(crate::PlaintextError::ValueTooLargeForU64)
                    })
                    .collect::<std::result::Result<Vec<_>, _>>()?;
                Self::encode_u64_chunk(&values, encoding, par, ctx, None)
            }
        }
    }
}

impl FheEncoderVariableTime<&[u64]> for PlaintextVec {
    type Error = Error;

    fn try_encode_vt(
        value: &[u64],
        encoding: Encoding,
        par: &Arc<BfvParameters>,
        variable_time: VariableTime,
    ) -> Result<Self> {
        Self::try_encode_with(value, encoding, par, |value, encoding, par, ctx| {
            Self::encode_u64_chunk(value, encoding, par, ctx, Some(variable_time))
        })
    }
}

impl FheEncoder<&[BigUint]> for PlaintextVec {
    type Error = Error;
    fn try_encode(value: &[BigUint], encoding: Encoding, par: &Arc<BfvParameters>) -> Result<Self> {
        Self::try_encode_with(value, encoding, par, Self::encode_biguint_chunk)
    }
}

impl FheEncoder<&[u64]> for PlaintextVec {
    type Error = Error;
    fn try_encode(value: &[u64], encoding: Encoding, par: &Arc<BfvParameters>) -> Result<Self> {
        Self::try_encode_with(value, encoding, par, |value, encoding, par, ctx| {
            Self::encode_u64_chunk(value, encoding, par, ctx, None)
        })
    }
}

#[cfg(test)]
mod tests {
    use crate::bfv::{BfvParameters, Encoding, PlaintextVec, parameters::BfvParametersBuilder};
    use fhe_traits::{FheDecoder, FheEncoder, FheEncoderVariableTime};
    use num_bigint::BigUint;
    use num_traits::Zero;
    use rand::rng;
    use std::error::Error;

    #[test]
    fn encode_decode() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for _ in 0..20 {
            for i in 1..5 {
                let params = BfvParameters::default_arc(1, 16);
                let a = params.plaintext();
                let q = fhe_math::zq::Modulus::new(a).unwrap();
                let a_vec = q.random_vec(params.degree() * i, &mut rng);

                let plaintexts = PlaintextVec::try_encode(
                    a_vec.as_slice(),
                    Encoding::poly_at_level(0),
                    &params,
                )?;
                assert_eq!(plaintexts.0.len(), i);

                for j in 0..i {
                    let b = Vec::<u64>::try_decode(&plaintexts.0[j], Encoding::poly_at_level(0))?;
                    assert_eq!(b, &a_vec[j * params.degree()..(j + 1) * params.degree()]);
                }

                let plaintexts_vt = PlaintextVec::try_encode_vt(
                    a_vec.as_slice(),
                    Encoding::poly_at_level(0),
                    &params,
                    fhe_traits::VariableTime::new(fhe_traits::PublicData::assert_public()),
                )?;
                assert_eq!(plaintexts_vt.0.len(), i);
                for (pt, pt_vt) in plaintexts.0.iter().zip(plaintexts_vt.0.iter()) {
                    assert_eq!(pt, pt_vt);
                }

                for j in 0..i {
                    let b =
                        Vec::<u64>::try_decode(&plaintexts_vt.0[j], Encoding::poly_at_level(0))?;
                    assert_eq!(b, &a_vec[j * params.degree()..(j + 1) * params.degree()]);
                }

                let plaintexts =
                    PlaintextVec::try_encode(a_vec.as_slice(), Encoding::simd(), &params)?;
                assert_eq!(plaintexts.0.len(), i);

                for j in 0..i {
                    let b = Vec::<u64>::try_decode(&plaintexts.0[j], Encoding::simd())?;
                    assert_eq!(b, &a_vec[j * params.degree()..(j + 1) * params.degree()]);
                }

                let plaintexts_vt = PlaintextVec::try_encode_vt(
                    a_vec.as_slice(),
                    Encoding::simd(),
                    &params,
                    fhe_traits::VariableTime::new(fhe_traits::PublicData::assert_public()),
                )?;
                assert_eq!(plaintexts_vt.0.len(), i);
                for (pt, pt_vt) in plaintexts.0.iter().zip(plaintexts_vt.0.iter()) {
                    assert_eq!(pt, pt_vt);
                }

                for j in 0..i {
                    let b = Vec::<u64>::try_decode(&plaintexts_vt.0[j], Encoding::simd())?;
                    assert_eq!(b, &a_vec[j * params.degree()..(j + 1) * params.degree()]);
                }
            }
        }
        let params = BfvParametersBuilder::new()
            .set_degree(16)
            .set_plaintext_modulus(17)
            .set_moduli_sizes(&[62])
            .build_arc()?;
        let a = vec![1u64];
        assert!(matches!(
            PlaintextVec::try_encode(a.as_slice(), Encoding::simd(), &params),
            Err(crate::Error::Encoding(
                crate::EncodingError::SimdUnavailable
            ))
        ));
        assert!(matches!(
            PlaintextVec::try_encode_vt(
                a.as_slice(),
                Encoding::simd(),
                &params,
                fhe_traits::VariableTime::new(fhe_traits::PublicData::assert_public()),
            ),
            Err(crate::Error::Encoding(
                crate::EncodingError::SimdUnavailable
            ))
        ));
        Ok(())
    }

    #[test]
    fn biguint_encoding_uses_shared_chunking() -> Result<(), Box<dyn Error>> {
        let modulus = BigUint::parse_bytes(b"340282366920938463463374607431768211507", 10)
            .ok_or("invalid test modulus")?;
        let params = BfvParametersBuilder::new()
            .set_degree(16)
            .set_plaintext_modulus_biguint(modulus.clone())
            .set_moduli_sizes(&[62, 62, 62, 62, 62])
            .build_arc()?;
        let values = (0u32..20).map(BigUint::from).collect::<Vec<_>>();

        let plaintexts =
            PlaintextVec::try_encode(values.as_slice(), Encoding::poly_at_level(0), &params)?;
        assert_eq!(plaintexts.len(), 2);

        for (plaintext, chunk) in plaintexts.iter().zip(values.chunks(params.degree())) {
            let mut expected = chunk.to_vec();
            expected.resize(params.degree(), BigUint::zero());
            assert_eq!(
                Vec::<BigUint>::try_decode(plaintext, Encoding::poly_at_level(0))?,
                expected
            );
        }
        Ok(())
    }

    #[test]
    fn empty_inputs_share_zero_encoding_path() -> Result<(), Box<dyn Error>> {
        let params = BfvParameters::default_arc(1, 16);
        let encoding = Encoding::poly();
        let constant = PlaintextVec::try_encode(&[] as &[u64], encoding.clone(), &params)?;
        let big = PlaintextVec::try_encode(&[] as &[BigUint], encoding.clone(), &params)?;
        let variable = PlaintextVec::try_encode_vt(
            &[] as &[u64],
            encoding.clone(),
            &params,
            fhe_traits::VariableTime::new(fhe_traits::PublicData::assert_public()),
        )?;

        for plaintexts in [&constant, &big, &variable] {
            assert_eq!(plaintexts.len(), 1);
            assert_eq!(
                Vec::<u64>::try_decode(&plaintexts[0], encoding.clone())?,
                vec![0; params.degree()]
            );
        }
        assert!(variable[0].poly_ntt.allows_variable_time_computations());
        Ok(())
    }
}
