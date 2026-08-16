//! Plaintext type in the BFV encryption scheme.
use crate::{
    Error, Result,
    bfv::{BfvParameters, Encoding, PlaintextVec},
};
use fhe_math::rq::{Context, Ntt, Poly, PowerBasis, traits::TryConvertFrom};
use fhe_traits::{
    FheDecoder, FheEncoder, FheEncoderVariableTime, FheParametrized, FhePlaintext, PublicData,
    VariableTime,
};
use num_bigint::{BigInt, BigUint, Sign};
use num_traits::{ToPrimitive, Zero};
use std::sync::Arc;
use zeroize::{Zeroize, Zeroizing};

use super::encoding::EncodingEnum;

enum PlaintextCoefficients {
    Small(Vec<u64>),
    Large(Vec<BigUint>),
}

/// A plaintext object, that encodes a vector according to a specific encoding.
#[derive(Debug, Clone, Eq)]
pub struct Plaintext {
    /// The parameters of the underlying BFV encryption scheme.
    pub(crate) par: Arc<BfvParameters>,
    /// The encoding of the plaintext, if known
    pub(crate) encoding: Option<Encoding>,
    /// Canonical plaintext representation.
    pub(crate) poly_ntt: Poly<Ntt>,
}

impl Zeroize for Plaintext {
    fn zeroize(&mut self) {
        self.poly_ntt.zeroize();
    }
}

impl Drop for Plaintext {
    fn drop(&mut self) {
        self.zeroize();
    }
}

impl FheParametrized for Plaintext {
    type Parameters = BfvParameters;
}

impl FhePlaintext for Plaintext {
    type Encoding = Encoding;
}

impl Plaintext {
    #[inline]
    pub(crate) fn validate_for(&self, par: &Arc<BfvParameters>) -> Result<()> {
        if !Arc::ptr_eq(&self.par, par) {
            return Err(Error::ParameterMismatch {
                left: crate::ParameterSource::Plaintext,
                right: crate::ParameterSource::Parameters,
            });
        }
        let level = self.level();
        let expected_ctx = par.context_at_level(level)?;
        self.validate_context(level, expected_ctx)
    }

    #[inline]
    pub(crate) fn validate_for_context(
        &self,
        par: &Arc<BfvParameters>,
        expected_level: usize,
        expected_ctx: &Arc<Context>,
    ) -> Result<()> {
        if !Arc::ptr_eq(&self.par, par) {
            return Err(Error::ParameterMismatch {
                left: crate::ParameterSource::Plaintext,
                right: crate::ParameterSource::Parameters,
            });
        }
        self.validate_context(expected_level, expected_ctx)
    }

    #[inline]
    fn validate_context(&self, expected_level: usize, expected_ctx: &Arc<Context>) -> Result<()> {
        let level = self.level();
        if level != expected_level {
            return Err(Error::InvalidLevel {
                level,
                min_level: expected_level,
                max_level: expected_level,
            });
        }
        if !Arc::ptr_eq(self.poly_ntt.ctx(), expected_ctx) && self.poly_ntt.ctx() != expected_ctx {
            return Err(crate::PlaintextError::PolynomialContextMismatch {
                level: expected_level,
            }
            .into());
        }
        Ok(())
    }

    fn coefficients(&self) -> PlaintextCoefficients {
        let poly = Zeroizing::new(self.poly_ntt.clone().into_power_basis());
        match self.par.plaintext.small() {
            Some(modulus)
                if self
                    .poly_ntt
                    .ctx()
                    .moduli()
                    .first()
                    .is_some_and(|ciphertext_modulus| **modulus < *ciphertext_modulus) =>
            {
                let coefficients = Vec::<u64>::try_from(poly.as_ref()).unwrap();
                let mut values = coefficients[..self.par.degree()].to_vec();
                modulus.reduce_vec(&mut values);
                PlaintextCoefficients::Small(values)
            }
            Some(_) => {
                let mut values = Vec::<BigUint>::from(poly.as_ref());
                self.par.plaintext.reduce_vec(&mut values);
                PlaintextCoefficients::Small(
                    values
                        .into_iter()
                        .map(|value| value.to_u64().unwrap())
                        .collect(),
                )
            }
            None => {
                let mut values = Vec::<BigUint>::from(poly.as_ref());
                self.par.plaintext.reduce_vec(&mut values);
                PlaintextCoefficients::Large(values)
            }
        }
    }

    fn resolve_encoding<O>(&self, encoding: O) -> Result<Encoding>
    where
        O: Into<Option<Encoding>>,
    {
        match (self.encoding.as_ref(), encoding.into()) {
            (None, None) => Err(crate::PlaintextError::MissingEncoding.into()),
            (Some(stored), Some(provided)) if stored != &provided => {
                Err(crate::EncodingError::Mismatch {
                    found: provided,
                    expected: stored.clone(),
                }
                .into())
            }
            (Some(stored), _) => Ok(stored.clone()),
            (None, Some(provided)) => Ok(provided),
        }
    }

    fn decode_simd_u64(&self, mut values: Vec<u64>) -> Result<Vec<u64>> {
        let op = self
            .par
            .ntt_operator
            .as_ref()
            .ok_or(crate::EncodingError::SimdUnavailable)?;
        op.forward(&mut values);
        let reordered = self
            .par
            .matrix_reps_index_map
            .iter()
            .map(|&index| values[index])
            .collect();
        values.zeroize();
        Ok(reordered)
    }

    pub(crate) fn to_poly(&self) -> Poly<Ntt> {
        let ctx_lvl = self.par.context_level_at(self.level()).unwrap();
        let ctx = &ctx_lvl.poly_context;

        let m = match self.coefficients() {
            PlaintextCoefficients::Small(values) => {
                let mut values = Zeroizing::new(values);
                let Some(modulus) = self.par.plaintext.small() else {
                    unreachable!("small plaintext values require the u64 modulus fast path");
                };
                let q_mod_t = ctx_lvl.cipher_plain_context.q_mod_t.to_u64().unwrap();
                modulus.scalar_mul_vec(&mut values, q_mod_t);
                Poly::<PowerBasis>::try_convert_from(values.as_slice(), ctx, false).unwrap()
            }
            PlaintextCoefficients::Large(mut values) => {
                self.par
                    .plaintext
                    .scalar_mul_vec(&mut values, &ctx_lvl.cipher_plain_context.q_mod_t);
                Poly::<PowerBasis>::try_convert_from(values.as_slice(), ctx, false).unwrap()
            }
        };

        let mut m = m.into_ntt();
        m *= &ctx_lvl.cipher_plain_context.delta;
        m
    }

    /// Generate a zero plaintext.
    pub fn zero(encoding: Encoding, par: &Arc<BfvParameters>) -> Result<Self> {
        let ctx = par.context_at_level(encoding.level)?;
        let poly_ntt = Poly::<Ntt>::zero(ctx);
        Ok(Self {
            par: par.clone(),
            encoding: Some(encoding),
            poly_ntt,
        })
    }

    /// Returns the level of this plaintext.
    #[must_use]
    pub fn level(&self) -> usize {
        self.par.moduli().len() - self.poly_ntt.ctx().moduli().len()
    }
}

unsafe impl Send for Plaintext {}

impl PartialEq for Plaintext {
    fn eq(&self, other: &Self) -> bool {
        let Self {
            par,
            encoding,
            poly_ntt,
        } = self;
        let Self {
            par: other_par,
            encoding: other_encoding,
            poly_ntt: other_poly_ntt,
        } = other;

        let mut eq = par == other_par;
        eq &= poly_ntt == other_poly_ntt;
        if encoding.is_some() && other_encoding.is_some() {
            eq &= encoding == other_encoding;
        }
        eq
    }
}

// Conversions.
impl TryConvertFrom<&Plaintext> for Poly<PowerBasis> {
    fn try_convert_from(
        pt: &Plaintext,
        ctx: &Arc<Context>,
        variable_time: bool,
    ) -> fhe_math::Result<Self> {
        if ctx
            != pt
                .par
                .context_at_level(pt.level())
                .map_err(|_| fhe_math::Error::ContextNotReachable)?
        {
            Err(fhe_math::Error::PolynomialContextMismatch)
        } else {
            let mut poly = pt.poly_ntt.clone();
            if variable_time {
                poly.allow_variable_time_computations(VariableTime::new(
                    PublicData::assert_public(),
                ));
            } else {
                poly.disallow_variable_time_computations();
            }
            Ok(poly.into_power_basis())
        }
    }
}

// Encoding and decoding.

impl<'a, const N: usize, T> FheEncoder<&'a [T; N]> for Plaintext
where
    Plaintext: FheEncoder<&'a [T], Error = Error>,
{
    type Error = Error;
    fn try_encode(value: &'a [T; N], encoding: Encoding, par: &Arc<BfvParameters>) -> Result<Self> {
        Plaintext::try_encode(value.as_ref(), encoding, par)
    }
}

impl<'a, T> FheEncoder<&'a Vec<T>> for Plaintext
where
    Plaintext: FheEncoder<&'a [T], Error = Error>,
{
    type Error = Error;
    fn try_encode(value: &'a Vec<T>, encoding: Encoding, par: &Arc<BfvParameters>) -> Result<Self> {
        Plaintext::try_encode(value.as_ref(), encoding, par)
    }
}

impl<'a> FheEncoder<&'a [BigUint]> for Plaintext {
    type Error = Error;
    fn try_encode(
        value: &'a [BigUint],
        encoding: Encoding,
        par: &Arc<BfvParameters>,
    ) -> Result<Self> {
        if value.len() > par.degree() {
            return Err(crate::PlaintextError::TooManyValues {
                actual: value.len(),
                maximum: par.degree(),
            }
            .into());
        }

        let v = PlaintextVec::try_encode(value, encoding, par)?;
        Ok(v[0].clone())
    }
}

impl<'a> FheEncoder<&'a [u64]> for Plaintext {
    type Error = Error;
    fn try_encode(value: &'a [u64], encoding: Encoding, par: &Arc<BfvParameters>) -> Result<Self> {
        if value.len() > par.degree() {
            return Err(crate::PlaintextError::TooManyValues {
                actual: value.len(),
                maximum: par.degree(),
            }
            .into());
        }
        let v = PlaintextVec::try_encode(value, encoding, par)?;
        Ok(v[0].clone())
    }
}

impl<'a> FheEncoderVariableTime<&'a [u64]> for Plaintext {
    type Error = Error;

    fn try_encode_vt(
        value: &'a [u64],
        encoding: Encoding,
        par: &Arc<BfvParameters>,
        variable_time: VariableTime,
    ) -> Result<Self> {
        if value.len() > par.degree() {
            return Err(crate::PlaintextError::TooManyValues {
                actual: value.len(),
                maximum: par.degree(),
            }
            .into());
        }
        let v = PlaintextVec::try_encode_vt(value, encoding, par, variable_time)?;
        Ok(v[0].clone())
    }
}

impl<'a> FheEncoder<&'a [i64]> for Plaintext {
    type Error = Error;
    fn try_encode(value: &'a [i64], encoding: Encoding, par: &Arc<BfvParameters>) -> Result<Self> {
        match par.plaintext.small() {
            Some(m) => {
                let w = Zeroizing::new(m.reduce_vec_i64(value));
                Plaintext::try_encode(w.as_ref() as &[u64], encoding, par)
            }
            None => {
                let modulus_int = BigInt::from_biguint(Sign::Plus, par.plaintext_big().clone());
                let v: Vec<BigUint> = value
                    .iter()
                    .map(|&x| {
                        let mut x_int = BigInt::from(x);
                        x_int %= &modulus_int;
                        if x_int < BigInt::zero() {
                            x_int += &modulus_int;
                        }
                        x_int.to_biguint().unwrap()
                    })
                    .collect();
                Plaintext::try_encode(v.as_slice(), encoding, par)
            }
        }
    }
}

impl FheDecoder<Plaintext> for Vec<BigUint> {
    fn try_decode<O>(pt: &Plaintext, encoding: O) -> Result<Vec<BigUint>>
    where
        O: Into<Option<Encoding>>,
    {
        let encoding = pt.resolve_encoding(encoding)?;
        let values = match pt.coefficients() {
            PlaintextCoefficients::Small(values) => values.into_iter().map(BigUint::from).collect(),
            PlaintextCoefficients::Large(values) => values,
        };

        match encoding.encoding {
            EncodingEnum::Poly => Ok(values),
            EncodingEnum::Simd => {
                let values = values
                    .into_iter()
                    .map(|value| {
                        value
                            .to_u64()
                            .ok_or(crate::PlaintextError::ValueTooLargeForU64)
                    })
                    .collect::<std::result::Result<Vec<_>, _>>()?;
                Ok(pt
                    .decode_simd_u64(values)?
                    .into_iter()
                    .map(BigUint::from)
                    .collect())
            }
        }
    }
    type Error = Error;
}

impl FheDecoder<Plaintext> for Vec<u64> {
    fn try_decode<O>(pt: &Plaintext, encoding: O) -> Result<Vec<u64>>
    where
        O: Into<Option<Encoding>>,
    {
        let encoding = pt.resolve_encoding(encoding)?;
        let values = match pt.coefficients() {
            PlaintextCoefficients::Small(values) => values,
            PlaintextCoefficients::Large(values) => values
                .into_iter()
                .map(|value| {
                    value
                        .to_u64()
                        .ok_or(crate::PlaintextError::ValueTooLargeForU64.into())
                })
                .collect::<Result<Vec<_>>>()?,
        };

        match encoding.encoding {
            EncodingEnum::Poly => Ok(values),
            EncodingEnum::Simd => pt.decode_simd_u64(values),
        }
    }

    type Error = Error;
}

impl FheDecoder<Plaintext> for Vec<i64> {
    fn try_decode<E>(pt: &Plaintext, encoding: E) -> Result<Vec<i64>>
    where
        E: Into<Option<Encoding>>,
    {
        if let Some(modulus) = pt.par.plaintext.small() {
            let values = Vec::<u64>::try_decode(pt, encoding)?;
            Ok(modulus.center_vec(&values))
        } else {
            let values = Vec::<BigUint>::try_decode(pt, encoding)?;
            let modulus_big = pt.par.plaintext_big();
            let modulus_int = BigInt::from_biguint(Sign::Plus, modulus_big.clone());
            let half_modulus = modulus_big / 2u32;

            Ok(values
                .iter()
                .map(|value| {
                    if value >= &half_modulus {
                        let value_int = BigInt::from_biguint(Sign::Plus, value.clone());
                        (value_int - &modulus_int).to_i64().unwrap()
                    } else {
                        value.to_i64().unwrap()
                    }
                })
                .collect())
        }
    }

    type Error = Error;
}

#[cfg(test)]
mod tests {
    use super::{Encoding, Plaintext};
    use crate::bfv::parameters::{BfvParameters, BfvParametersBuilder};
    use fhe_math::rq::{Ntt, Poly};
    use fhe_traits::{FheDecoder, FheEncoder, FheEncoderVariableTime};
    use num_bigint::BigUint;
    use num_traits::Zero;
    use rand::rng;
    use std::error::Error;
    use zeroize::Zeroize;

    #[test]
    fn try_encode() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        // The default test parameters support both Poly and Simd encodings
        let params = BfvParameters::default_arc(1, 16);
        // random_vec returns Vec<u64>
        let a = params.plaintext();
        // use modulus directly to generate random u64s
        let q = fhe_math::zq::Modulus::new(a).unwrap();
        let a_vec = q.random_vec(params.degree(), &mut rng);

        let plaintext = Plaintext::try_encode(&[0u64; 17], Encoding::poly(), &params);
        assert!(plaintext.is_err());

        let plaintext = Plaintext::try_encode(&a_vec, Encoding::poly(), &params)?;
        assert_eq!(Vec::<u64>::try_decode(&plaintext, Encoding::poly())?, a_vec);

        let plaintext = Plaintext::try_encode(&a_vec, Encoding::simd(), &params);
        assert!(plaintext.is_ok());

        let plaintext = Plaintext::try_encode(&[1u64], Encoding::poly(), &params);
        assert!(plaintext.is_ok());

        // The following parameters do not allow for Simd encoding
        let params = BfvParametersBuilder::new()
            .set_degree(16)
            .set_plaintext_modulus(2)
            .set_moduli(&[4611686018326724609])
            .build_arc()?;

        let a = 2u64;
        let q = fhe_math::zq::Modulus::new(a).unwrap();
        let a_vec = q.random_vec(params.degree(), &mut rng);

        let plaintext = Plaintext::try_encode(&a_vec, Encoding::poly(), &params);
        assert!(plaintext.is_ok());

        let plaintext = Plaintext::try_encode(&a_vec, Encoding::simd(), &params);
        assert!(plaintext.is_err());

        Ok(())
    }

    #[test]
    fn try_encode_variable_time_marks_public_plaintexts() -> Result<(), Box<dyn Error>> {
        let params = BfvParameters::default_arc(1, 16);
        let values = [1u64, 2, 3, 4];
        let encoding = Encoding::poly();
        let variable_time = fhe_traits::VariableTime::new(fhe_traits::PublicData::assert_public());

        let constant_time = Plaintext::try_encode(&values, encoding.clone(), &params)?;
        let public =
            Plaintext::try_encode_vt(values.as_slice(), encoding.clone(), &params, variable_time)?;
        assert_eq!(constant_time, public);
        assert!(!constant_time.poly_ntt.allows_variable_time_computations());
        assert!(public.poly_ntt.allows_variable_time_computations());

        let public_zero =
            Plaintext::try_encode_vt(&[] as &[u64], encoding, &params, variable_time)?;
        assert!(public_zero.poly_ntt.allows_variable_time_computations());
        Ok(())
    }

    #[test]
    fn try_encode_big() -> Result<(), Box<dyn Error>> {
        // Test with big plaintext
        let p_val = BigUint::parse_bytes(b"340282366920938463463374607431768211507", 10).unwrap();
        let params = BfvParametersBuilder::new()
            .set_degree(16)
            .set_plaintext_modulus_biguint(p_val.clone())
            .set_moduli_sizes(&[62, 62, 62, 62, 62])
            .build_arc()?;

        let vals = vec![p_val.clone() - 1u32, BigUint::from(123u32)];
        let plaintext = Plaintext::try_encode(&vals, Encoding::poly(), &params)?;

        let decoded: Vec<BigUint> = Vec::<BigUint>::try_decode(&plaintext, Encoding::poly())?;
        assert_eq!(decoded[0], p_val - 1u32);
        assert_eq!(decoded[1], BigUint::from(123u32));
        assert_eq!(decoded[2], BigUint::zero());

        Ok(())
    }

    #[test]
    fn encode_decode() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let params = BfvParameters::default_arc(1, 16);
        let a = params.plaintext();
        let q = fhe_math::zq::Modulus::new(a).unwrap();
        let a_vec = q.random_vec(params.degree(), &mut rng);

        let plaintext = Plaintext::try_encode(&a_vec, Encoding::simd(), &params);
        assert!(plaintext.is_ok());
        let b = Vec::<u64>::try_decode(&plaintext?, Encoding::simd())?;
        assert_eq!(b, a_vec);

        // center_vec replacement logic for test
        let mut a_signed = vec![];
        for x in &a_vec {
            if *x >= a / 2 {
                a_signed.push((*x as i64) - (a as i64));
            } else {
                a_signed.push(*x as i64);
            }
        }

        let plaintext = Plaintext::try_encode(&a_signed, Encoding::poly(), &params);
        assert!(plaintext.is_ok());
        let b = Vec::<i64>::try_decode(&plaintext?, Encoding::poly())?;
        assert_eq!(b, a_signed);

        let plaintext = Plaintext::try_encode(&a_signed, Encoding::simd(), &params);
        assert!(plaintext.is_ok());
        let b = Vec::<i64>::try_decode(&plaintext?, Encoding::simd())?;
        assert_eq!(b, a_signed);

        Ok(())
    }

    #[test]
    fn partial_eq() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let params = BfvParameters::default_arc(1, 16);
        let a = params.plaintext();
        let q = fhe_math::zq::Modulus::new(a).unwrap();
        let a_vec = q.random_vec(params.degree(), &mut rng);

        let plaintext = Plaintext::try_encode(&a_vec, Encoding::poly(), &params)?;
        let mut same_plaintext = Plaintext::try_encode(&a_vec, Encoding::poly(), &params)?;
        assert_eq!(plaintext, same_plaintext);

        // Equality also holds when there is no encoding specified. In this test, we use
        // the fact that we can set it to None directly, but such a partial plaintext
        // will be created during decryption since we do not specify the encoding at the
        // time.
        same_plaintext.encoding = None;
        assert_eq!(plaintext, same_plaintext);

        Ok(())
    }

    #[test]
    fn try_decode_errors() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let params = BfvParameters::default_arc(1, 16);
        let a = params.plaintext();
        let q = fhe_math::zq::Modulus::new(a).unwrap();
        let a_vec = q.random_vec(params.degree(), &mut rng);

        let mut plaintext = Plaintext::try_encode(&a_vec, Encoding::poly(), &params)?;

        assert!(Vec::<u64>::try_decode(&plaintext, None).is_ok());
        let e = Vec::<u64>::try_decode(&plaintext, Encoding::simd());
        assert!(e.is_err());
        assert_eq!(
            e.unwrap_err(),
            crate::Error::Encoding(crate::EncodingError::Mismatch {
                found: Encoding::simd(),
                expected: Encoding::poly(),
            })
        );
        let e = Vec::<u64>::try_decode(&plaintext, Encoding::poly_at_level(1));
        assert!(e.is_err());
        assert_eq!(
            e.unwrap_err(),
            crate::Error::Encoding(crate::EncodingError::Mismatch {
                found: Encoding::poly_at_level(1),
                expected: Encoding::poly(),
            })
        );

        plaintext.encoding = None;
        let e = Vec::<u64>::try_decode(&plaintext, None);
        assert!(e.is_err());
        assert_eq!(
            e.unwrap_err(),
            crate::Error::Plaintext(crate::PlaintextError::MissingEncoding)
        );

        Ok(())
    }

    #[test]
    fn zero() -> Result<(), Box<dyn Error>> {
        let params = BfvParameters::default_arc(1, 16);
        let plaintext = Plaintext::zero(Encoding::poly(), &params)?;

        assert_eq!(
            plaintext.poly_ntt,
            Poly::<Ntt>::zero(params.context_at_level(0)?)
        );

        Ok(())
    }

    #[test]
    fn zeroize() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let params = BfvParameters::default_arc(1, 16);
        let a = params.plaintext();
        let q = fhe_math::zq::Modulus::new(a).unwrap();
        let a_vec = q.random_vec(params.degree(), &mut rng);
        let mut plaintext = Plaintext::try_encode(&a_vec, Encoding::poly(), &params)?;

        plaintext.zeroize();

        assert_eq!(plaintext, Plaintext::zero(Encoding::poly(), &params)?);

        Ok(())
    }

    #[test]
    fn try_encode_level() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        // The default test parameters support both Poly and Simd encodings
        let params = BfvParameters::default_arc(10, 16);
        let a = params.plaintext();
        let q = fhe_math::zq::Modulus::new(a).unwrap();
        let a_vec = q.random_vec(params.degree(), &mut rng);

        for level in 0..10 {
            let plaintext = Plaintext::try_encode(&a_vec, Encoding::poly_at_level(level), &params)?;
            assert_eq!(plaintext.level(), level);
            let plaintext = Plaintext::try_encode(&a_vec, Encoding::simd_at_level(level), &params)?;
            assert_eq!(plaintext.level(), level);
        }

        Ok(())
    }
}
