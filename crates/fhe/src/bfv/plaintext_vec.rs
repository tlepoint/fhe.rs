use std::{cmp::min, ops::Deref, sync::Arc};

use crate::VariableTime;
use fhe_math::rq::{Context, Ntt, Poly, PowerBasis, traits::TryConvertFrom};
use num_bigint::BigUint;
use num_traits::{ToPrimitive, Zero};
use zeroize_derive::{Zeroize, ZeroizeOnDrop};

use crate::{
    Result,
    bfv::{Encoding, Parameters, Plaintext},
};

/// Plaintexts encoded from successive degree-sized chunks of a sequence.
#[derive(Zeroize, ZeroizeOnDrop)]
pub struct PlaintextVec(Vec<Plaintext>);

impl Deref for PlaintextVec {
    type Target = [Plaintext];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl PlaintextVec {
    fn encode_with<T>(
        value: &[T],
        encoding: Encoding,
        level: usize,
        par: &Parameters,
        mut encode_chunk: impl FnMut(&[T], &Encoding, &Parameters, &Arc<Context>) -> Result<Poly<Ntt>>,
    ) -> Result<Self> {
        if encoding == Encoding::Simd && par.inner.ntt_operator.is_none() {
            return Err(crate::EncodingError::SimdUnavailable.into());
        }

        let ctx = par.context_at_level(level)?;
        let num_plaintexts = value.len().div_ceil(par.degree()).max(1);
        let plaintexts = (0..num_plaintexts)
            .map(|index| {
                let start = index * par.degree();
                let end = min(value.len(), start + par.degree());
                let poly_ntt = encode_chunk(&value[start..end], &encoding, par, ctx)?;
                Ok(Plaintext {
                    par: par.clone(),
                    poly_ntt,
                })
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Self(plaintexts))
    }

    pub(crate) fn encode_u64_chunk(
        value: &[u64],
        encoding: &Encoding,
        par: &Parameters,
        ctx: &Arc<Context>,
        variable_time: Option<VariableTime>,
    ) -> Result<Poly<Ntt>> {
        let reduced = match par.inner.plaintext.small() {
            Some(modulus) => {
                let mut reduced = value.to_vec();
                modulus.reduce_vec(&mut reduced);
                reduced
            }
            None => value
                .iter()
                .map(|value| {
                    (BigUint::from(*value) % par.plaintext_modulus())
                        .to_u64()
                        .unwrap()
                })
                .collect(),
        };
        let reduced = zeroize::Zeroizing::new(reduced);
        let value = reduced.as_slice();
        let mut coefficients = zeroize::Zeroizing::new(vec![0u64; par.degree()]);
        match *encoding {
            Encoding::Polynomial => coefficients[..value.len()].copy_from_slice(value),
            Encoding::Simd => {
                for (index, &coefficient) in value.iter().enumerate() {
                    coefficients[par.inner.matrix_reps_index_map[index]] = coefficient;
                }
                let ntt_operator = par
                    .inner
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
            Poly::<PowerBasis>::try_convert_from_public(
                coefficients.as_slice(),
                ctx,
                variable_time,
            )?
        } else {
            Poly::<PowerBasis>::try_convert_from(coefficients.as_slice(), ctx)?
        };
        Ok(poly.into_ntt())
    }

    pub(crate) fn encode_biguint_chunk(
        value: &[BigUint],
        encoding: &Encoding,
        par: &Parameters,
        ctx: &Arc<Context>,
    ) -> Result<Poly<Ntt>> {
        let reduced: Vec<_> = value
            .iter()
            .map(|value| value % par.plaintext_modulus())
            .collect();
        let value = reduced.as_slice();
        match *encoding {
            Encoding::Polynomial => {
                let mut coefficients = vec![BigUint::zero(); par.degree()];
                coefficients[..value.len()].clone_from_slice(value);
                Ok(Poly::<PowerBasis>::try_convert_from(coefficients.as_slice(), ctx)?.into_ntt())
            }
            Encoding::Simd => {
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

impl PlaintextVec {
    /// Encode a sequence in degree-sized chunks at level zero. Empty input
    /// produces one zero plaintext; the last chunk is padded with zeros.
    pub fn encode(par: &Parameters, values: &[u64], encoding: Encoding) -> Result<Self> {
        Self::encode_at_level(par, values, encoding, 0)
    }
    /// Encode degree-sized chunks at an explicit level.
    pub fn encode_at_level(
        par: &Parameters,
        values: &[u64],
        encoding: Encoding,
        level: usize,
    ) -> Result<Self> {
        Self::encode_with(
            values,
            encoding,
            level,
            par,
            |values, encoding, par, ctx| Self::encode_u64_chunk(values, encoding, par, ctx, None),
        )
    }
    /// Encode a sequence in degree-sized chunks at level zero. Empty input
    /// produces one zero plaintext; the last chunk is padded with zeros.
    pub fn encode_biguint(
        par: &Parameters,
        values: &[BigUint],
        encoding: Encoding,
    ) -> Result<Self> {
        Self::encode_biguint_at_level(par, values, encoding, 0)
    }
    /// Encode degree-sized chunks at an explicit level.
    pub fn encode_biguint_at_level(
        par: &Parameters,
        values: &[BigUint],
        encoding: Encoding,
        level: usize,
    ) -> Result<Self> {
        Self::encode_with(values, encoding, level, par, Self::encode_biguint_chunk)
    }
    /// Encode a sequence in degree-sized chunks at level zero. Empty input
    /// produces one zero plaintext; the last chunk is padded with zeros.
    pub fn encode_public(
        par: &Parameters,
        values: &[u64],
        encoding: Encoding,
        permission: VariableTime,
    ) -> Result<Self> {
        Self::encode_public_at_level(par, values, encoding, 0, permission)
    }
    /// Encode degree-sized chunks at an explicit level.
    pub fn encode_public_at_level(
        par: &Parameters,
        values: &[u64],
        encoding: Encoding,
        level: usize,
        permission: VariableTime,
    ) -> Result<Self> {
        Self::encode_with(
            values,
            encoding,
            level,
            par,
            |values, encoding, par, ctx| {
                Self::encode_u64_chunk(values, encoding, par, ctx, Some(permission))
            },
        )
    }
}

#[cfg(test)]
mod tests {
    use crate::bfv::{Encoding, Parameters, PlaintextVec, parameters::ParametersBuilder};

    use num_bigint::BigUint;
    use num_traits::Zero;
    use rand::rng;
    use std::error::Error;

    #[test]
    fn encode_decode() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for _ in 0..20 {
            for i in 1..5 {
                let params = Parameters::test_parameters(1, 16);
                let a = params.plaintext_modulus_u64().unwrap();
                let q = fhe_math::zq::Modulus::new(a).unwrap();
                let a_vec = q.random_vec(params.degree() * i, &mut rng);

                let plaintexts =
                    PlaintextVec::encode(&params, a_vec.as_slice(), Encoding::Polynomial)?;
                assert_eq!(plaintexts.0.len(), i);

                for j in 0..i {
                    let b = (plaintexts.0[j]).decode(Encoding::Polynomial)?;
                    assert_eq!(b, &a_vec[j * params.degree()..(j + 1) * params.degree()]);
                }

                let plaintexts_vt = PlaintextVec::encode_public(
                    &params,
                    a_vec.as_slice(),
                    Encoding::Polynomial,
                    crate::VariableTime::new(crate::PublicData::assert_public()),
                )?;
                assert_eq!(plaintexts_vt.0.len(), i);
                for (pt, pt_vt) in plaintexts.0.iter().zip(plaintexts_vt.0.iter()) {
                    assert_eq!(pt, pt_vt);
                }

                for j in 0..i {
                    let b = (plaintexts_vt.0[j]).decode(Encoding::Polynomial)?;
                    assert_eq!(b, &a_vec[j * params.degree()..(j + 1) * params.degree()]);
                }

                let plaintexts = PlaintextVec::encode(&params, a_vec.as_slice(), Encoding::Simd)?;
                assert_eq!(plaintexts.0.len(), i);

                for j in 0..i {
                    let b = (plaintexts.0[j]).decode(Encoding::Simd)?;
                    assert_eq!(b, &a_vec[j * params.degree()..(j + 1) * params.degree()]);
                }

                let plaintexts_vt = PlaintextVec::encode_public(
                    &params,
                    a_vec.as_slice(),
                    Encoding::Simd,
                    crate::VariableTime::new(crate::PublicData::assert_public()),
                )?;
                assert_eq!(plaintexts_vt.0.len(), i);
                for (pt, pt_vt) in plaintexts.0.iter().zip(plaintexts_vt.0.iter()) {
                    assert_eq!(pt, pt_vt);
                }

                for j in 0..i {
                    let b = (plaintexts_vt.0[j]).decode(Encoding::Simd)?;
                    assert_eq!(b, &a_vec[j * params.degree()..(j + 1) * params.degree()]);
                }
            }
        }
        let params = ParametersBuilder::new()
            .degree(16)
            .plaintext_modulus(17_u64)
            .ciphertext_modulus_bits([62])
            .build()?;
        let a = vec![1u64];
        assert!(matches!(
            PlaintextVec::encode(&params, a.as_slice(), Encoding::Simd),
            Err(crate::Error::Encoding(
                crate::EncodingError::SimdUnavailable
            ))
        ));
        assert!(matches!(
            PlaintextVec::encode_public(
                &params,
                a.as_slice(),
                Encoding::Simd,
                crate::VariableTime::new(crate::PublicData::assert_public())
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
        let params = ParametersBuilder::new()
            .degree(16)
            .plaintext_modulus(modulus.clone())
            .ciphertext_modulus_bits([62, 62, 62, 62, 62])
            .build()?;
        let values = (0u32..20).map(BigUint::from).collect::<Vec<_>>();

        let plaintexts =
            PlaintextVec::encode_biguint(&params, values.as_slice(), Encoding::Polynomial)?;
        assert_eq!(plaintexts.len(), 2);

        for (plaintext, chunk) in plaintexts.iter().zip(values.chunks(params.degree())) {
            let mut expected = chunk.to_vec();
            expected.resize(params.degree(), BigUint::zero());
            assert_eq!(plaintext.decode_biguint(Encoding::Polynomial)?, expected);
        }
        Ok(())
    }

    #[test]
    fn empty_inputs_share_zero_encoding_path() -> Result<(), Box<dyn Error>> {
        let params = Parameters::test_parameters(1, 16);
        let encoding = Encoding::Polynomial;
        let constant = PlaintextVec::encode(&params, &[] as &[u64], encoding)?;
        let big = PlaintextVec::encode_biguint(&params, &[] as &[BigUint], encoding)?;
        let variable = PlaintextVec::encode_public(
            &params,
            &[] as &[u64],
            encoding,
            crate::VariableTime::new(crate::PublicData::assert_public()),
        )?;

        for plaintexts in [&constant, &big, &variable] {
            assert_eq!(plaintexts.len(), 1);
            assert_eq!(plaintexts[0].decode(encoding)?, vec![0; params.degree()]);
        }
        assert!(variable[0].poly_ntt.allows_variable_time_computations());
        Ok(())
    }
}
