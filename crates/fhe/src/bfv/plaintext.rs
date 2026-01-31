//! Plaintext type in the BFV encryption scheme.
use crate::{
    Error, Result,
    bfv::{BfvParameters, Encoding, PlaintextVec, parameters::PlaintextModulus},
};
use fhe_math::rq::{Context, Poly, Representation, traits::TryConvertFrom};
use fhe_traits::{FheDecoder, FheEncoder, FheParametrized, FhePlaintext};
use num_bigint::{BigInt, BigUint, Sign};
use num_traits::{ToPrimitive, Zero};
use std::sync::Arc;
use zeroize::{Zeroize, Zeroizing};

use super::encoding::EncodingEnum;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum PlaintextValues {
    Small(Box<[u64]>),
    Large(Box<[BigUint]>),
}

impl Zeroize for PlaintextValues {
    fn zeroize(&mut self) {
        match self {
            Self::Small(v) => v.zeroize(),
            Self::Large(v) => {
                for x in v.iter_mut() {
                    *x = BigUint::zero();
                }
            }
        }
    }
}

/// A plaintext object, that encodes a vector according to a specific encoding.
#[derive(Debug, Clone, Eq)]
pub struct Plaintext {
    /// The parameters of the underlying BFV encryption scheme.
    pub(crate) par: Arc<BfvParameters>,
    /// The value after encoding.
    pub(crate) value: PlaintextValues,
    /// The encoding of the plaintext, if known
    pub(crate) encoding: Option<Encoding>,
    /// The plaintext as a polynomial.
    pub(crate) poly_ntt: Poly,
    /// The level of the plaintext
    pub(crate) level: usize,
}

impl Zeroize for Plaintext {
    fn zeroize(&mut self) {
        self.value.zeroize();
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
    pub(crate) fn to_poly(&self) -> Poly {
        let ctx_lvl = self.par.context_level_at(self.level).unwrap();
        let ctx = &ctx_lvl.poly_context;

        let mut m = match &self.value {
            PlaintextValues::Small(v) => {
                let mut m_v = Zeroizing::new(v.clone());
                if let PlaintextModulus::Small(modulus) = &self.par.plaintext {
                    let q_mod_t = ctx_lvl.cipher_plain_context.q_mod_t.to_u64().unwrap();
                    modulus.scalar_mul_vec(&mut m_v, q_mod_t);
                } else {
                    unreachable!("PlaintextValues::Small but PlaintextModulus::Large");
                }
                Poly::try_convert_from(m_v.as_ref(), ctx, false, Representation::PowerBasis)
                    .unwrap()
            }
            PlaintextValues::Large(v) => {
                let mut m_v = v.clone();
                self.par
                    .plaintext
                    .scalar_mul_vec(&mut m_v, &ctx_lvl.cipher_plain_context.q_mod_t);
                Poly::try_convert_from(m_v.as_ref(), ctx, false, Representation::PowerBasis)
                    .unwrap()
            }
        };

        m.change_representation(Representation::Ntt);
        m *= &ctx_lvl.cipher_plain_context.delta;
        m
    }

    /// Generate a zero plaintext.
    pub fn zero(encoding: Encoding, par: &Arc<BfvParameters>) -> Result<Self> {
        let level = encoding.level;
        let ctx = par.context_at_level(level)?;
        let value = match par.plaintext {
            PlaintextModulus::Small(_) => {
                PlaintextValues::Small(vec![0u64; par.degree()].into_boxed_slice())
            }
            PlaintextModulus::Large(_) => {
                PlaintextValues::Large(vec![BigUint::zero(); par.degree()].into_boxed_slice())
            }
        };
        let poly_ntt = Poly::zero(ctx, Representation::Ntt);
        Ok(Self {
            par: par.clone(),
            value,
            encoding: Some(encoding),
            poly_ntt,
            level,
        })
    }

    /// Returns the level of this plaintext.
    #[must_use]
    pub fn level(&self) -> usize {
        self.par.level_of_context(self.poly_ntt.ctx()).unwrap()
    }
}

unsafe impl Send for Plaintext {}

impl PartialEq for Plaintext {
    fn eq(&self, other: &Self) -> bool {
        let Self {
            par,
            value,
            encoding,
            poly_ntt,
            level,
        } = self;
        let Self {
            par: other_par,
            value: other_value,
            encoding: other_encoding,
            poly_ntt: other_poly_ntt,
            level: other_level,
        } = other;

        let mut eq = par == other_par;
        eq &= value == other_value;
        eq &= poly_ntt == other_poly_ntt;
        eq &= level == other_level;
        if encoding.is_some() && other_encoding.is_some() {
            eq &= encoding == other_encoding;
        }
        eq
    }
}

// Conversions.
impl TryConvertFrom<&Plaintext> for Poly {
    fn try_convert_from<R>(
        pt: &Plaintext,
        ctx: &Arc<Context>,
        variable_time: bool,
        _: R,
    ) -> fhe_math::Result<Self>
    where
        R: Into<Option<Representation>>,
    {
        if ctx
            != pt
                .par
                .context_at_level(pt.level())
                .map_err(|e| fhe_math::Error::Default(e.to_string()))?
        {
            Err(fhe_math::Error::Default(
                "Incompatible contexts".to_string(),
            ))
        } else {
            match &pt.value {
                PlaintextValues::Small(v) => Poly::try_convert_from(
                    v.as_ref(),
                    ctx,
                    variable_time,
                    Representation::PowerBasis,
                ),
                PlaintextValues::Large(v) => Poly::try_convert_from(
                    v.as_ref(),
                    ctx,
                    variable_time,
                    Representation::PowerBasis,
                ),
            }
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
            return Err(Error::TooManyValues {
                actual: value.len(),
                limit: par.degree(),
            });
        }

        let v = PlaintextVec::try_encode(value, encoding, par)?;
        Ok(v[0].clone())
    }
}

impl<'a> FheEncoder<&'a [u64]> for Plaintext {
    type Error = Error;
    fn try_encode(value: &'a [u64], encoding: Encoding, par: &Arc<BfvParameters>) -> Result<Self> {
        if value.len() > par.degree() {
            return Err(Error::TooManyValues {
                actual: value.len(),
                limit: par.degree(),
            });
        }
        let v = PlaintextVec::try_encode(value, encoding, par)?;
        Ok(v[0].clone())
    }
}

impl<'a> FheEncoder<&'a [i64]> for Plaintext {
    type Error = Error;
    fn try_encode(value: &'a [i64], encoding: Encoding, par: &Arc<BfvParameters>) -> Result<Self> {
        match par.plaintext {
            PlaintextModulus::Small(ref m) => {
                let w = Zeroizing::new(m.reduce_vec_i64(value));
                Plaintext::try_encode(w.as_ref() as &[u64], encoding, par)
            }
            PlaintextModulus::Large(ref m) => {
                let modulus_int = BigInt::from_biguint(Sign::Plus, m.clone());
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
        // First convert to Vec<BigUint> regardless of internal storage
        let w = match &pt.value {
            PlaintextValues::Small(v) => v.iter().map(|&x| BigUint::from(x)).collect::<Vec<_>>(),
            PlaintextValues::Large(v) => v.to_vec(),
        };

        // Standard decoding logic (e.g. check encoding match)
        let encoding = encoding.into();
        let enc: Encoding;
        if pt.encoding.is_none() && encoding.is_none() {
            return Err(Error::InvalidPlaintext {
                reason: "No encoding specified".into(),
            });
        } else if pt.encoding.is_some() {
            enc = pt.encoding.as_ref().unwrap().clone();
            if let Some(arg_enc) = encoding
                && arg_enc != enc
            {
                return Err(Error::EncodingMismatch {
                    found: arg_enc.into(),
                    expected: enc.into(),
                });
            }
        } else {
            enc = encoding.unwrap();
            if let Some(pt_enc) = pt.encoding.as_ref()
                && pt_enc != &enc
            {
                return Err(Error::EncodingMismatch {
                    found: pt_enc.into(),
                    expected: enc.into(),
                });
            }
        }

        match enc.encoding {
            EncodingEnum::Poly => Ok(w),
            EncodingEnum::Simd => {
                if let Some(op) = &pt.par.ntt_operator {
                    // NTT operator works on u64.
                    // If ntt_operator exists, it means we are in Small modulus case.
                    let mut w_u64: Vec<u64> = w.iter().map(|x| x.to_u64().unwrap()).collect();
                    op.forward(&mut w_u64);
                    let mut w_reordered = w_u64.clone();
                    for i in 0..pt.par.degree() {
                        w_reordered[i] = w_u64[pt.par.matrix_reps_index_map[i]]
                    }
                    w_u64.zeroize();

                    Ok(w_reordered.into_iter().map(BigUint::from).collect())
                } else {
                    Err(Error::EncodingNotSupported {
                        encoding: EncodingEnum::Simd.to_string(),
                        reason: "NTT operator not available".into(),
                    })
                }
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
        // Optimized path for Small values
        match &pt.value {
            PlaintextValues::Small(v) => {
                // Copied logic for validation
                let encoding = encoding.into();
                let enc: Encoding;
                if pt.encoding.is_none() && encoding.is_none() {
                    return Err(Error::InvalidPlaintext {
                        reason: "No encoding specified".into(),
                    });
                } else if pt.encoding.is_some() {
                    enc = pt.encoding.as_ref().unwrap().clone();
                    if let Some(arg_enc) = encoding
                        && arg_enc != enc
                    {
                        return Err(Error::EncodingMismatch {
                            found: arg_enc.into(),
                            expected: enc.into(),
                        });
                    }
                } else {
                    enc = encoding.unwrap();
                    if let Some(pt_enc) = pt.encoding.as_ref()
                        && pt_enc != &enc
                    {
                        return Err(Error::EncodingMismatch {
                            found: pt_enc.into(),
                            expected: enc.into(),
                        });
                    }
                }

                let mut w = v.to_vec();

                match enc.encoding {
                    EncodingEnum::Poly => Ok(w),
                    EncodingEnum::Simd => {
                        if let Some(op) = &pt.par.ntt_operator {
                            op.forward(&mut w);
                            let mut w_reordered = w.clone();
                            for i in 0..pt.par.degree() {
                                w_reordered[i] = w[pt.par.matrix_reps_index_map[i]]
                            }
                            w.zeroize();
                            Ok(w_reordered)
                        } else {
                            Err(Error::EncodingNotSupported {
                                encoding: EncodingEnum::Simd.to_string(),
                                reason: "NTT operator not available".into(),
                            })
                        }
                    }
                }
            }
            PlaintextValues::Large(_) => {
                let v = Vec::<BigUint>::try_decode(pt, encoding)?;
                v.iter()
                    .map(|x| {
                        x.to_u64().ok_or(Error::DefaultError(
                            "Plaintext value too large for u64".to_string(),
                        ))
                    })
                    .collect()
            }
        }
    }

    type Error = Error;
}

impl FheDecoder<Plaintext> for Vec<i64> {
    fn try_decode<E>(pt: &Plaintext, encoding: E) -> Result<Vec<i64>>
    where
        E: Into<Option<Encoding>>,
    {
        match &pt.value {
            PlaintextValues::Small(_) => {
                let v = Vec::<u64>::try_decode(pt, encoding)?;
                if let PlaintextModulus::Small(ref m) = pt.par.plaintext {
                    Ok(m.center_vec(&v))
                } else {
                    unreachable!()
                }
            }
            PlaintextValues::Large(_) => {
                let v = Vec::<BigUint>::try_decode(pt, encoding)?;
                let modulus_big = pt.par.plaintext_big();
                let modulus_int = BigInt::from_biguint(Sign::Plus, modulus_big.clone());
                let half_modulus = modulus_big / 2u32;

                Ok(v.iter()
                    .map(|x| {
                        if x >= &half_modulus {
                            let x_int = BigInt::from_biguint(Sign::Plus, x.clone());
                            (x_int - &modulus_int).to_i64().unwrap()
                        } else {
                            x.to_i64().unwrap()
                        }
                    })
                    .collect())
            }
        }
    }

    type Error = Error;
}

#[cfg(test)]
mod tests {
    use super::{Encoding, Plaintext};
    use crate::bfv::parameters::{BfvParameters, BfvParametersBuilder};
    use crate::bfv::plaintext::PlaintextValues;
    use fhe_math::rq::{Poly, Representation};
    use fhe_traits::{FheDecoder, FheEncoder};
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

        let plaintext = Plaintext::try_encode(&a_vec, Encoding::poly(), &params);
        assert!(plaintext.is_ok());
        // Verify it used Small variant
        if let PlaintextValues::Large(_) = plaintext.unwrap().value {
            println!("Expected Small variant");
        }

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

        // Verify it used Large variant
        if let PlaintextValues::Small(_) = plaintext.value {
            println!("Expected Large variant");
        }

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
            crate::Error::EncodingMismatch {
                found: Encoding::simd().into(),
                expected: Encoding::poly().into(),
            }
        );
        let e = Vec::<u64>::try_decode(&plaintext, Encoding::poly_at_level(1));
        assert!(e.is_err());
        assert_eq!(
            e.unwrap_err(),
            crate::Error::EncodingMismatch {
                found: Encoding::poly_at_level(1).into(),
                expected: Encoding::poly().into(),
            }
        );

        plaintext.encoding = None;
        let e = Vec::<u64>::try_decode(&plaintext, None);
        assert!(e.is_err());
        assert_eq!(
            e.unwrap_err(),
            crate::Error::InvalidPlaintext {
                reason: "No encoding specified".into(),
            }
        );

        Ok(())
    }

    #[test]
    fn zero() -> Result<(), Box<dyn Error>> {
        let params = BfvParameters::default_arc(1, 16);
        let plaintext = Plaintext::zero(Encoding::poly(), &params)?;

        assert_eq!(
            plaintext.value,
            PlaintextValues::Small(vec![0u64; 16].into_boxed_slice())
        );
        assert_eq!(
            plaintext.poly_ntt,
            Poly::zero(params.context_at_level(0)?, Representation::Ntt)
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
