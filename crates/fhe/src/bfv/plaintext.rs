//! Plaintext type in the BFV encryption scheme.
use crate::{
    Error, Result,
    bfv::{Encoding, Parameters, PlaintextVec},
};
use crate::{PublicData, VariableTime};
use fhe_math::rq::{Context, Ntt, Poly, PowerBasis, traits::TryConvertFrom};
use num_bigint::{BigInt, BigUint, Sign};
use num_traits::ToPrimitive;
use std::sync::Arc;
use zeroize::{Zeroize, Zeroizing};

enum PlaintextCoefficients {
    Small(Vec<u64>),
    Large(Vec<BigUint>),
}

/// A polynomial at one modulus-switching level, with shared BFV parameters.
/// Equality compares the polynomial, level, and defining parameter settings.
/// Encoding intent and original message length are not stored; supply an
/// explicit interpretation when decoding, including after decryption.
#[derive(Clone, PartialEq, Eq)]
pub struct Plaintext {
    /// The parameters of the underlying BFV encryption scheme.
    pub(crate) par: Parameters,
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

impl Plaintext {
    #[inline]
    pub(crate) fn validate_for(&self, par: &Parameters) -> Result<()> {
        if !Parameters::compatible(&self.par, par) {
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
        par: &Parameters,
        expected_level: usize,
        expected_ctx: &Arc<Context>,
    ) -> Result<()> {
        if !Parameters::compatible(&self.par, par) {
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
        match self.par.inner.plaintext.small() {
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
                self.par.inner.plaintext.reduce_vec(&mut values);
                PlaintextCoefficients::Small(
                    values
                        .into_iter()
                        .map(|value| value.to_u64().unwrap())
                        .collect(),
                )
            }
            None => {
                let mut values = Vec::<BigUint>::from(poly.as_ref());
                self.par.inner.plaintext.reduce_vec(&mut values);
                PlaintextCoefficients::Large(values)
            }
        }
    }

    fn decode_simd_u64(&self, mut values: Vec<u64>) -> Result<Vec<u64>> {
        let op = self
            .par
            .inner
            .ntt_operator
            .as_ref()
            .ok_or(crate::EncodingError::SimdUnavailable)?;
        op.forward(&mut values);
        let reordered = self
            .par
            .inner
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
                let Some(modulus) = self.par.inner.plaintext.small() else {
                    unreachable!("small plaintext values require the u64 modulus fast path");
                };
                let q_mod_t = ctx_lvl.cipher_plain_context.q_mod_t.to_u64().unwrap();
                modulus.scalar_mul_vec(&mut values, q_mod_t);
                Poly::<PowerBasis>::try_convert_from(values.as_slice(), ctx).unwrap()
            }
            PlaintextCoefficients::Large(mut values) => {
                self.par
                    .inner
                    .plaintext
                    .scalar_mul_vec(&mut values, &ctx_lvl.cipher_plain_context.q_mod_t);
                Poly::<PowerBasis>::try_convert_from(values.as_slice(), ctx).unwrap()
            }
        };

        let mut m = m.into_ntt();
        m *= &ctx_lvl.cipher_plain_context.delta;
        m
    }

    /// Generate a zero plaintext.
    pub fn zero(par: &Parameters, level: usize) -> Result<Self> {
        let ctx = par.context_at_level(level)?;
        let poly_ntt = Poly::<Ntt>::zero(ctx);
        Ok(Self {
            par: par.clone(),
            poly_ntt,
        })
    }

    /// Returns the level of this plaintext.
    #[must_use]
    pub fn level(&self) -> usize {
        self.par.moduli().len() - self.poly_ntt.ctx().moduli().len()
    }
}

impl std::fmt::Debug for Plaintext {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Plaintext")
            .field("degree", &self.par.degree())
            .field("level", &self.level())
            .finish_non_exhaustive()
    }
}

// Conversions.
impl TryConvertFrom<&Plaintext> for Poly<PowerBasis> {
    fn try_convert_from_with_timing(
        pt: &Plaintext,
        ctx: &Arc<Context>,
        permission: Option<VariableTime>,
    ) -> fhe_math::Result<Self> {
        let variable_time = permission.is_some();
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

impl Plaintext {
    /// Shared parameters of this plaintext.
    #[must_use]
    pub fn parameters(&self) -> &Parameters {
        &self.par
    }

    /// Encode unsigned values at level zero. Inputs are reduced modulo the
    /// plaintext modulus and zero padded to the polynomial degree.
    pub fn encode(par: &Parameters, values: &[u64], encoding: Encoding) -> Result<Self> {
        Self::encode_at_level(par, values, encoding, 0)
    }

    /// Encode unsigned values at an explicit level. Increasing the level drops
    /// ciphertext moduli. More than `degree` values or an unsupported SIMD
    /// modulus returns an error.
    pub fn encode_at_level(
        par: &Parameters,
        values: &[u64],
        encoding: Encoding,
        level: usize,
    ) -> Result<Self> {
        Self::encode_with(
            par,
            values,
            encoding,
            level,
            |values, encoding, par, ctx| {
                PlaintextVec::encode_u64_chunk(values, encoding, par, ctx, None)
            },
        )
    }

    /// Encode public unsigned values at level zero with variable-time
    /// permission.
    pub fn encode_public(
        par: &Parameters,
        values: &[u64],
        encoding: Encoding,
        permission: VariableTime,
    ) -> Result<Self> {
        Self::encode_public_at_level(par, values, encoding, 0, permission)
    }

    /// Encode public unsigned values at an explicit level with variable-time
    /// permission.
    pub fn encode_public_at_level(
        par: &Parameters,
        values: &[u64],
        encoding: Encoding,
        level: usize,
        permission: VariableTime,
    ) -> Result<Self> {
        Self::encode_with(
            par,
            values,
            encoding,
            level,
            |values, encoding, par, ctx| {
                PlaintextVec::encode_u64_chunk(values, encoding, par, ctx, Some(permission))
            },
        )
    }

    /// Encode arbitrary unsigned integers, reduced modulo the plaintext
    /// modulus. Big integer input does not extend the current machine-word
    /// SIMD capability.
    pub fn encode_biguint(
        par: &Parameters,
        values: &[BigUint],
        encoding: Encoding,
    ) -> Result<Self> {
        Self::encode_biguint_at_level(par, values, encoding, 0)
    }

    /// Encode arbitrary unsigned integers at an explicit level.
    pub fn encode_biguint_at_level(
        par: &Parameters,
        values: &[BigUint],
        encoding: Encoding,
        level: usize,
    ) -> Result<Self> {
        Self::encode_with(
            par,
            values,
            encoding,
            level,
            PlaintextVec::encode_biguint_chunk,
        )
    }

    /// Encode signed integers modulo the plaintext modulus at level zero.
    pub fn encode_signed(par: &Parameters, values: &[i64], encoding: Encoding) -> Result<Self> {
        Self::encode_signed_at_level(par, values, encoding, 0)
    }

    /// Encode signed integers modulo the plaintext modulus at an explicit
    /// level.
    pub fn encode_signed_at_level(
        par: &Parameters,
        values: &[i64],
        encoding: Encoding,
        level: usize,
    ) -> Result<Self> {
        match par.inner.plaintext.small() {
            Some(m) => {
                let values = Zeroizing::new(m.reduce_vec_i64(values));
                Self::encode_at_level(par, &values, encoding, level)
            }
            None => {
                let modulus = BigInt::from_biguint(Sign::Plus, par.plaintext_modulus().clone());
                let values: Vec<BigUint> = values
                    .iter()
                    .map(|&value| {
                        let value = BigInt::from(value);
                        ((value % &modulus + &modulus) % &modulus)
                            .to_biguint()
                            .unwrap()
                    })
                    .collect();
                Self::encode_biguint_at_level(par, &values, encoding, level)
            }
        }
    }

    fn encode_with<T>(
        par: &Parameters,
        values: &[T],
        encoding: Encoding,
        level: usize,
        encode: impl FnOnce(&[T], &Encoding, &Parameters, &Arc<Context>) -> Result<Poly<Ntt>>,
    ) -> Result<Self> {
        if values.len() > par.degree() {
            return Err(crate::PlaintextError::TooManyValues {
                actual: values.len(),
                maximum: par.degree(),
            }
            .into());
        }
        if encoding == Encoding::Simd && par.inner.ntt_operator.is_none() {
            return Err(crate::EncodingError::SimdUnavailable.into());
        }
        let ctx = par.context_at_level(level)?;
        let poly_ntt = encode(values, &encoding, par, ctx)?;
        Ok(Self {
            par: par.clone(),
            poly_ntt,
        })
    }
}

impl Plaintext {
    /// Decode exactly `degree` values using the supplied interpretation.
    pub fn decode_biguint(&self, encoding: Encoding) -> Result<Vec<BigUint>> {
        let values = match self.coefficients() {
            PlaintextCoefficients::Small(values) => values.into_iter().map(BigUint::from).collect(),
            PlaintextCoefficients::Large(values) => values,
        };

        match encoding {
            Encoding::Polynomial => Ok(values),
            Encoding::Simd => {
                let values = values
                    .into_iter()
                    .map(|value| {
                        value
                            .to_u64()
                            .ok_or(crate::PlaintextError::ValueTooLargeForU64)
                    })
                    .collect::<std::result::Result<Vec<_>, _>>()?;
                Ok(self
                    .decode_simd_u64(values)?
                    .into_iter()
                    .map(BigUint::from)
                    .collect())
            }
        }
    }
}

impl Plaintext {
    /// Decode exactly `degree` values using the supplied interpretation.
    /// An output that does not fit the requested integer type returns an error.
    pub fn decode(&self, encoding: Encoding) -> Result<Vec<u64>> {
        let values = match self.coefficients() {
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

        match encoding {
            Encoding::Polynomial => Ok(values),
            Encoding::Simd => self.decode_simd_u64(values),
        }
    }
}

impl Plaintext {
    /// Decode exactly `degree` values using the supplied interpretation.
    /// Values use centered representatives; an unrepresentable output returns
    /// an error.
    pub fn decode_signed(&self, encoding: Encoding) -> Result<Vec<i64>> {
        if let Some(modulus) = self.par.inner.plaintext.small() {
            let values = self.decode(encoding)?;
            Ok(modulus.center_vec(&values))
        } else {
            let values = self.decode_biguint(encoding)?;
            let modulus_big = self.par.plaintext_modulus();
            let modulus_int = BigInt::from_biguint(Sign::Plus, modulus_big.clone());
            let half_modulus = (modulus_big + 1u32) / 2u32;

            values
                .iter()
                .map(|value| {
                    let centered = if value >= &half_modulus {
                        let value_int = BigInt::from_biguint(Sign::Plus, value.clone());
                        (value_int - &modulus_int).to_i64()
                    } else {
                        value.to_i64()
                    };
                    centered.ok_or(crate::PlaintextError::ValueTooLargeForI64.into())
                })
                .collect()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{Encoding, Plaintext};
    use crate::bfv::parameters::{Parameters, ParametersBuilder};
    use fhe_math::rq::{Ntt, Poly};

    use num_bigint::BigUint;
    use num_traits::Zero;
    use rand::rng;
    use std::error::Error;
    use zeroize::Zeroize;

    #[test]
    fn try_encode() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        // The default test parameters support both Poly and Simd encodings
        let params = Parameters::test_parameters(1, 16);
        // random_vec returns Vec<u64>
        let a = params.plaintext_modulus_u64().unwrap();
        // use modulus directly to generate random u64s
        let q = fhe_math::zq::Modulus::new(a).unwrap();
        let a_vec = q.random_vec(params.degree(), &mut rng);

        let plaintext = Plaintext::encode(&params, &[0u64; 17], Encoding::Polynomial);
        assert!(plaintext.is_err());

        let plaintext = Plaintext::encode(&params, &a_vec, Encoding::Polynomial)?;
        assert_eq!(plaintext.decode(Encoding::Polynomial)?, a_vec);

        let plaintext = Plaintext::encode(&params, &a_vec, Encoding::Simd);
        assert!(plaintext.is_ok());

        let plaintext = Plaintext::encode(&params, &[1u64], Encoding::Polynomial);
        assert!(plaintext.is_ok());

        // The following parameters do not allow for Simd encoding
        let params = ParametersBuilder::new()
            .degree(16)
            .plaintext_modulus(2_u64)
            .ciphertext_moduli([4611686018326724609])
            .build()?;

        let a = 2u64;
        let q = fhe_math::zq::Modulus::new(a).unwrap();
        let a_vec = q.random_vec(params.degree(), &mut rng);

        let plaintext = Plaintext::encode(&params, &a_vec, Encoding::Polynomial);
        assert!(plaintext.is_ok());

        let plaintext = Plaintext::encode(&params, &a_vec, Encoding::Simd);
        assert!(plaintext.is_err());

        Ok(())
    }

    #[test]
    fn try_encode_variable_time_marks_public_plaintexts() -> Result<(), Box<dyn Error>> {
        let params = Parameters::test_parameters(1, 16);
        let values = [1u64, 2, 3, 4];
        let encoding = Encoding::Polynomial;
        let variable_time = crate::VariableTime::new(crate::PublicData::assert_public());

        let constant_time = Plaintext::encode(&params, &values, encoding)?;
        let public = Plaintext::encode_public(&params, values.as_slice(), encoding, variable_time)?;
        assert_eq!(constant_time, public);
        assert!(!constant_time.poly_ntt.allows_variable_time_computations());
        assert!(public.poly_ntt.allows_variable_time_computations());

        let public_zero =
            Plaintext::encode_public(&params, &[] as &[u64], encoding, variable_time)?;
        assert!(public_zero.poly_ntt.allows_variable_time_computations());
        Ok(())
    }

    #[test]
    fn try_encode_big() -> Result<(), Box<dyn Error>> {
        // Test with big plaintext
        let p_val = BigUint::parse_bytes(b"340282366920938463463374607431768211507", 10).unwrap();
        let params = ParametersBuilder::new()
            .degree(16)
            .plaintext_modulus(p_val.clone())
            .ciphertext_modulus_bits([62, 62, 62, 62, 62])
            .build()?;

        let vals = vec![p_val.clone() - 1u32, BigUint::from(123u32)];
        let plaintext = Plaintext::encode_biguint(&params, &vals, Encoding::Polynomial)?;

        let decoded: Vec<BigUint> = plaintext.decode_biguint(Encoding::Polynomial)?;
        assert_eq!(decoded[0], p_val - 1u32);
        assert_eq!(decoded[1], BigUint::from(123u32));
        assert_eq!(decoded[2], BigUint::zero());

        Ok(())
    }

    #[test]
    fn encode_decode() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let params = Parameters::test_parameters(1, 16);
        let a = params.plaintext_modulus_u64().unwrap();
        let q = fhe_math::zq::Modulus::new(a).unwrap();
        let mut a_vec = q.random_vec(params.degree(), &mut rng);
        // Always exercise the midpoint: for odd a, floor(a / 2) stays positive.
        a_vec[params.degree() - 1] = a / 2;

        let plaintext = Plaintext::encode(&params, &a_vec, Encoding::Simd);
        assert!(plaintext.is_ok());
        let b = (plaintext?).decode(Encoding::Simd)?;
        assert_eq!(b, a_vec);

        // Center into [-a / 2, a / 2); the first negative residue is ceil(a / 2).
        let mut a_signed = vec![];
        for x in &a_vec {
            if *x >= a.div_ceil(2) {
                a_signed.push((*x as i64) - (a as i64));
            } else {
                a_signed.push(*x as i64);
            }
        }

        let plaintext = Plaintext::encode_signed(&params, &a_signed, Encoding::Polynomial);
        assert!(plaintext.is_ok());
        let b = (plaintext?).decode_signed(Encoding::Polynomial)?;
        assert_eq!(b, a_signed);

        let plaintext = Plaintext::encode_signed(&params, &a_signed, Encoding::Simd);
        assert!(plaintext.is_ok());
        let b = (plaintext?).decode_signed(Encoding::Simd)?;
        assert_eq!(b, a_signed);

        Ok(())
    }

    #[test]
    fn partial_eq() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let params = Parameters::test_parameters(1, 16);
        let a = params.plaintext_modulus_u64().unwrap();
        let q = fhe_math::zq::Modulus::new(a).unwrap();
        let a_vec = q.random_vec(params.degree(), &mut rng);

        let plaintext = Plaintext::encode(&params, &a_vec, Encoding::Polynomial)?;
        let same_plaintext = Plaintext::encode(&params, &a_vec, Encoding::Polynomial)?;
        assert_eq!(plaintext, same_plaintext);

        let sk = crate::bfv::SecretKey::generate(&params, &mut rng);
        assert_eq!(plaintext, sk.decrypt(&sk.encrypt(&plaintext, &mut rng)?)?);

        Ok(())
    }

    #[test]
    fn decoding_interpretation_is_explicit() -> Result<(), Box<dyn Error>> {
        let params = Parameters::test_parameters(2, 16);
        let plaintext = Plaintext::encode_at_level(&params, &[1_u64], Encoding::Polynomial, 1)?;
        assert_eq!(plaintext.decode(Encoding::Polynomial)?[0], 1);
        assert_eq!(plaintext.decode(Encoding::Simd)?, vec![1; params.degree()]);
        assert_eq!(plaintext.level(), 1);
        let non_simd = Parameters::builder()
            .degree(16)
            .plaintext_modulus(17_u64)
            .ciphertext_modulus_bits([62])
            .build()?;
        let plaintext = Plaintext::zero(&non_simd, 0)?;
        assert_eq!(
            plaintext.decode(Encoding::Simd),
            Err(crate::EncodingError::SimdUnavailable.into())
        );
        Ok(())
    }

    #[test]
    fn zero() -> Result<(), Box<dyn Error>> {
        let params = Parameters::test_parameters(1, 16);
        let plaintext = Plaintext::zero(&params, 0)?;

        assert_eq!(
            plaintext.poly_ntt,
            Poly::<Ntt>::zero(params.context_at_level(0)?)
        );

        Ok(())
    }

    #[test]
    fn zeroize() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let params = Parameters::test_parameters(1, 16);
        let a = params.plaintext_modulus_u64().unwrap();
        let q = fhe_math::zq::Modulus::new(a).unwrap();
        let a_vec = q.random_vec(params.degree(), &mut rng);
        let mut plaintext = Plaintext::encode(&params, &a_vec, Encoding::Polynomial)?;

        plaintext.zeroize();

        assert_eq!(plaintext, Plaintext::zero(&params, 0)?);

        Ok(())
    }

    #[test]
    fn try_encode_level() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        // The default test parameters support both Poly and Simd encodings
        let params = Parameters::test_parameters(10, 16);
        let a = params.plaintext_modulus_u64().unwrap();
        let q = fhe_math::zq::Modulus::new(a).unwrap();
        let a_vec = q.random_vec(params.degree(), &mut rng);

        for level in 0..10 {
            let plaintext =
                Plaintext::encode_at_level(&params, &a_vec, Encoding::Polynomial, level)?;
            assert_eq!(plaintext.level(), level);
            let plaintext = Plaintext::encode_at_level(&params, &a_vec, Encoding::Simd, level)?;
            assert_eq!(plaintext.level(), level);
        }

        Ok(())
    }

    #[test]
    fn signed_decoding_boundaries_and_overflow() -> Result<(), Box<dyn Error>> {
        let par = Parameters::test_parameters(2, 16);
        assert_eq!(par.plaintext_modulus_u64().unwrap(), 1153);
        for encoding in [Encoding::Polynomial, Encoding::Simd] {
            // Values just outside the centered interval wrap modulo 1153.
            let pt =
                Plaintext::encode_signed(&par, &[575i64, 576, 577, -575, -576, -577], encoding)?;
            let values = pt.decode_signed(encoding)?;
            assert_eq!(&values[..6], &[575, 576, -576, -575, -576, 576]);
        }
        for t in [
            BigUint::from(1u32) << 100usize,
            (BigUint::from(1u32) << 100usize) + 1u32,
        ] {
            let par = ParametersBuilder::new()
                .degree(16)
                .plaintext_modulus(t.clone())
                .ciphertext_modulus_bits([62, 62, 62])
                .build()?;
            let pt =
                Plaintext::encode_signed(&par, &[i64::MIN, -1, 0, i64::MAX], Encoding::Polynomial)?;
            assert_eq!(
                &pt.decode_signed(Encoding::Polynomial)?[..4],
                &[i64::MIN, -1, 0, i64::MAX]
            );
            for value in [
                BigUint::from(1u32) << 80usize,
                &t - (BigUint::from(1u32) << 80usize),
            ] {
                let pt = Plaintext::encode_biguint(&par, &[value], Encoding::Polynomial)?;
                assert!(matches!(
                    pt.decode_signed(Encoding::Polynomial),
                    Err(crate::Error::Plaintext(
                        crate::PlaintextError::ValueTooLargeForI64
                    ))
                ));
            }
        }
        Ok(())
    }
}
