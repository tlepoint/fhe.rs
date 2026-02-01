use std::{cmp::min, ops::Deref, sync::Arc};

use fhe_math::rq::{Poly, Representation, traits::TryConvertFrom};
use fhe_traits::{FheEncoder, FheEncoderVariableTime, FheParametrized, FhePlaintext};
use num_bigint::BigUint;
use num_traits::{ToPrimitive, Zero};
use zeroize_derive::{Zeroize, ZeroizeOnDrop};

use crate::{
    Error, Result,
    bfv::{BfvParameters, Encoding, Plaintext, PlaintextValues},
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

impl FheEncoderVariableTime<&[u64]> for PlaintextVec {
    type Error = Error;

    unsafe fn try_encode_vt(
        value: &[u64],
        encoding: Encoding,
        par: &Arc<BfvParameters>,
    ) -> Result<Self> {
        if value.is_empty() {
            return Ok(PlaintextVec(vec![Plaintext::zero(encoding, par)?]));
        }
        if encoding.encoding == EncodingEnum::Simd && par.ntt_operator.is_none() {
            return Err(Error::EncodingNotSupported {
                encoding: EncodingEnum::Simd.to_string(),
                reason: "NTT operator not available".into(),
            });
        }
        let ctx = par.context_at_level(encoding.level)?;
        let num_plaintexts = value.len().div_ceil(par.degree());

        Ok(PlaintextVec(
            (0..num_plaintexts)
                .map(|i| {
                    let slice = &value[i * par.degree()..min(value.len(), (i + 1) * par.degree())];
                    let mut v = vec![0u64; par.degree()];
                    match encoding.encoding {
                        EncodingEnum::Poly => v[..slice.len()].copy_from_slice(slice),
                        EncodingEnum::Simd => {
                            for i in 0..slice.len() {
                                v[par.matrix_reps_index_map[i]] = slice[i];
                            }
                            let ntt_operator =
                                par.ntt_operator.as_ref().ok_or(Error::InvalidPlaintext {
                                    reason: "No Ntt operator".into(),
                                })?;
                            unsafe { ntt_operator.backward_vt(v.as_mut_ptr()) };
                        }
                    };

                    let mut poly =
                        Poly::try_convert_from(&v, ctx, true, Representation::PowerBasis)?;
                    poly.change_representation(Representation::Ntt);

                    let value_enum = match par.plaintext {
                        crate::bfv::PlaintextModulus::Small { .. } => {
                            PlaintextValues::Small(v.into_boxed_slice())
                        }
                        crate::bfv::PlaintextModulus::Large(_) => PlaintextValues::Large(
                            v.iter()
                                .map(|&x| BigUint::from(x))
                                .collect::<Vec<_>>()
                                .into_boxed_slice(),
                        ),
                    };

                    Ok(Plaintext {
                        par: par.clone(),
                        value: value_enum,
                        encoding: Some(encoding.clone()),
                        poly_ntt: poly,
                        level: encoding.level,
                    })
                })
                .collect::<Result<Vec<Plaintext>>>()?,
        ))
    }
}

impl FheEncoder<&[BigUint]> for PlaintextVec {
    type Error = Error;
    fn try_encode(value: &[BigUint], encoding: Encoding, par: &Arc<BfvParameters>) -> Result<Self> {
        if value.is_empty() {
            return Ok(PlaintextVec(vec![Plaintext::zero(encoding, par)?]));
        }
        if encoding.encoding == EncodingEnum::Simd && par.ntt_operator.is_none() {
            return Err(Error::EncodingNotSupported {
                encoding: EncodingEnum::Simd.to_string(),
                reason: "NTT operator not available".into(),
            });
        }
        let ctx = par.context_at_level(encoding.level)?;
        let num_plaintexts = value.len().div_ceil(par.degree());

        Ok(PlaintextVec(
            (0..num_plaintexts)
                .map(|i| {
                    let slice = &value[i * par.degree()..min(value.len(), (i + 1) * par.degree())];
                    let mut v = vec![BigUint::zero(); par.degree()];
                    match encoding.encoding {
                        EncodingEnum::Poly => v[..slice.len()].clone_from_slice(slice),
                        EncodingEnum::Simd => {
                            let mut v_u64 = vec![0u64; par.degree()];
                            for i in 0..slice.len() {
                                v_u64[par.matrix_reps_index_map[i]] =
                                    slice[i].to_u64().ok_or(Error::DefaultError(
                                        "Value too large for SIMD encoding".to_string(),
                                    ))?;
                            }
                            par.ntt_operator
                                .as_ref()
                                .ok_or(Error::InvalidPlaintext {
                                    reason: "No Ntt operator".into(),
                                })?
                                .backward(&mut v_u64);

                            v = v_u64.into_iter().map(BigUint::from).collect();
                        }
                    };

                    let mut poly = Poly::try_convert_from(
                        v.as_slice(),
                        ctx,
                        false,
                        Representation::PowerBasis,
                    )?;
                    poly.change_representation(Representation::Ntt);

                    let value_enum = match &par.plaintext {
                        crate::bfv::PlaintextModulus::Small { modulus_big, .. } => {
                            PlaintextValues::Small(
                                v.iter()
                                    .map(|x| (x % modulus_big).to_u64().unwrap())
                                    .collect::<Vec<_>>()
                                    .into_boxed_slice(),
                            )
                        }
                        crate::bfv::PlaintextModulus::Large(_) => {
                            PlaintextValues::Large(v.into_boxed_slice())
                        }
                    };

                    Ok(Plaintext {
                        par: par.clone(),
                        value: value_enum,
                        encoding: Some(encoding.clone()),
                        poly_ntt: poly,
                        level: encoding.level,
                    })
                })
                .collect::<Result<Vec<Plaintext>>>()?,
        ))
    }
}

impl FheEncoder<&[u64]> for PlaintextVec {
    type Error = Error;
    fn try_encode(value: &[u64], encoding: Encoding, par: &Arc<BfvParameters>) -> Result<Self> {
        if value.is_empty() {
            return Ok(PlaintextVec(vec![Plaintext::zero(encoding, par)?]));
        }
        if encoding.encoding == EncodingEnum::Simd && par.ntt_operator.is_none() {
            return Err(Error::EncodingNotSupported {
                encoding: EncodingEnum::Simd.to_string(),
                reason: "NTT operator not available".into(),
            });
        }
        let ctx = par.context_at_level(encoding.level)?;
        let num_plaintexts = value.len().div_ceil(par.degree());

        Ok(PlaintextVec(
            (0..num_plaintexts)
                .map(|i| {
                    let slice = &value[i * par.degree()..min(value.len(), (i + 1) * par.degree())];
                    let mut v = vec![0u64; par.degree()];
                    match encoding.encoding {
                        EncodingEnum::Poly => v[..slice.len()].copy_from_slice(slice),
                        EncodingEnum::Simd => {
                            for i in 0..slice.len() {
                                v[par.matrix_reps_index_map[i]] = slice[i];
                            }
                            par.ntt_operator
                                .as_ref()
                                .ok_or(Error::InvalidPlaintext {
                                    reason: "No Ntt operator".into(),
                                })?
                                .backward(&mut v);
                        }
                    };

                    let mut poly =
                        Poly::try_convert_from(&v, ctx, false, Representation::PowerBasis)?;
                    poly.change_representation(Representation::Ntt);

                    let value_enum = match par.plaintext {
                        crate::bfv::PlaintextModulus::Small { .. } => {
                            PlaintextValues::Small(v.into_boxed_slice())
                        }
                        crate::bfv::PlaintextModulus::Large(_) => PlaintextValues::Large(
                            v.iter()
                                .map(|&x| BigUint::from(x))
                                .collect::<Vec<_>>()
                                .into_boxed_slice(),
                        ),
                    };

                    Ok(Plaintext {
                        par: par.clone(),
                        value: value_enum,
                        encoding: Some(encoding.clone()),
                        poly_ntt: poly,
                        level: encoding.level,
                    })
                })
                .collect::<Result<Vec<Plaintext>>>()?,
        ))
    }
}

#[cfg(test)]
mod tests {
    use crate::bfv::{BfvParameters, Encoding, PlaintextVec, parameters::BfvParametersBuilder};
    use fhe_traits::{FheDecoder, FheEncoder, FheEncoderVariableTime};
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

                let plaintexts_vt = unsafe {
                    PlaintextVec::try_encode_vt(
                        a_vec.as_slice(),
                        Encoding::poly_at_level(0),
                        &params,
                    )?
                };
                assert_eq!(plaintexts_vt.0.len(), i);
                for (pt, pt_vt) in plaintexts.0.iter().zip(plaintexts_vt.0.iter()) {
                    assert_eq!(pt.value, pt_vt.value);
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

                let plaintexts_vt = unsafe {
                    PlaintextVec::try_encode_vt(a_vec.as_slice(), Encoding::simd(), &params)?
                };
                assert_eq!(plaintexts_vt.0.len(), i);
                for (pt, pt_vt) in plaintexts.0.iter().zip(plaintexts_vt.0.iter()) {
                    assert_eq!(pt.value, pt_vt.value);
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
            Err(crate::Error::EncodingNotSupported { .. })
        ));
        assert!(matches!(
            unsafe { PlaintextVec::try_encode_vt(a.as_slice(), Encoding::simd(), &params) },
            Err(crate::Error::EncodingNotSupported { .. })
        ));
        Ok(())
    }
}
