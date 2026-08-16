//! Leveled evaluation keys for the BFV encryption scheme.

use crate::bfv::{BfvParameters, Ciphertext, SecretKey, keys::GaloisKey, traits::TryConvertFrom};
use crate::proto::bfv::{EvaluationKey as EvaluationKeyProto, GaloisKey as GaloisKeyProto};
use crate::{Error, Result, SerializationError};
use fhe_math::rq::{NttShoup, Poly, PowerBasis};
use fhe_math::zq::Modulus;
use fhe_traits::{DeserializeParametrized, FheParametrized, Serialize};
use prost::Message;
use rand::{CryptoRng, Rng as RngCore};
use std::collections::{HashMap, HashSet};
use std::sync::Arc;
use zeroize::{Zeroize, ZeroizeOnDrop};

/// Evaluation key for the BFV encryption scheme.
///
/// An evaluation key enables one or several of the following operations:
/// - column rotation
/// - row rotation
/// - oblivious expansion
/// - inner sum
#[derive(Debug, PartialEq, Eq)]
pub struct EvaluationKey {
    par: Arc<BfvParameters>,

    ciphertext_level: usize,
    evaluation_key_level: usize,

    /// Map from Galois keys exponents to Galois keys
    gk: HashMap<usize, GaloisKey>,

    /// Map from rotation index to Galois key exponent
    rot_to_gk_exponent: HashMap<usize, usize>,

    /// Monomials used in expansion
    monomials: Vec<Poly<NttShoup>>,
}

impl EvaluationKey {
    /// Reports whether the evaluation key enables to compute an homomorphic
    /// inner sums.
    #[must_use]
    pub fn supports_inner_sum(&self) -> bool {
        let mut ret = self.gk.contains_key(&(self.par.degree() * 2 - 1));
        let mut i = 1;
        while i < self.par.degree() / 2 {
            ret &= self
                .gk
                .contains_key(self.rot_to_gk_exponent.get(&i).unwrap());
            i *= 2
        }
        ret
    }

    /// Computes the homomorphic inner sum.
    pub fn computes_inner_sum(&self, ct: &Ciphertext) -> Result<Ciphertext> {
        self.validate_ciphertext(ct)?;
        if !self.supports_inner_sum() {
            Err(crate::EvaluationKeyError::Unsupported {
                operation: crate::EvaluationOperation::InnerSum,
            }
            .into())
        } else {
            let mut out = ct.clone();
            let mut tmp = Ciphertext::zero(&ct.par);

            let mut i = 1;
            while i < ct.par.degree() / 2 {
                let exponent =
                    self.rot_to_gk_exponent
                        .get(&i)
                        .ok_or(crate::EvaluationKeyError::Missing {
                            component: crate::EvaluationKeyComponent::GaloisExponent { step: i },
                        })?;
                let gk = self
                    .gk
                    .get(exponent)
                    .ok_or(crate::EvaluationKeyError::Missing {
                        component: crate::EvaluationKeyComponent::GaloisKey { element: *exponent },
                    })?;
                gk.relinearize_into(&out, &mut tmp)?;
                out += &tmp;
                i *= 2
            }

            let row_rotation_element = self.par.degree() * 2 - 1;
            let gk =
                self.gk
                    .get(&row_rotation_element)
                    .ok_or(crate::EvaluationKeyError::Missing {
                        component: crate::EvaluationKeyComponent::GaloisKey {
                            element: row_rotation_element,
                        },
                    })?;
            gk.relinearize_into(&out, &mut tmp)?;
            out += &tmp;

            Ok(out)
        }
    }

    /// Reports whether the evaluation key enables to rotate the rows of the
    /// plaintext.
    #[must_use]
    pub fn supports_row_rotation(&self) -> bool {
        self.gk.contains_key(&(self.par.degree() * 2 - 1))
    }

    /// Homomorphically rotate the rows of the plaintext
    pub fn rotates_rows(&self, ct: &Ciphertext) -> Result<Ciphertext> {
        self.validate_ciphertext(ct)?;
        if !self.supports_row_rotation() {
            Err(crate::EvaluationKeyError::Unsupported {
                operation: crate::EvaluationOperation::RowRotation,
            }
            .into())
        } else {
            let row_rotation_element = self.par.degree() * 2 - 1;
            let gk =
                self.gk
                    .get(&row_rotation_element)
                    .ok_or(crate::EvaluationKeyError::Missing {
                        component: crate::EvaluationKeyComponent::GaloisKey {
                            element: row_rotation_element,
                        },
                    })?;
            let mut out = Ciphertext::zero(&ct.par);
            gk.relinearize_into(ct, &mut out)?;
            Ok(out)
        }
    }

    /// Reports whether the evaluation key enables to rotate the columns of the
    /// plaintext.
    #[must_use]
    pub fn supports_column_rotation_by(&self, i: usize) -> bool {
        if let Some(exp) = self.rot_to_gk_exponent.get(&i) {
            self.gk.contains_key(exp)
        } else {
            false
        }
    }

    /// Homomorphically rotate the columns of the plaintext
    pub fn rotates_columns_by(&self, ct: &Ciphertext, i: usize) -> Result<Ciphertext> {
        self.validate_ciphertext(ct)?;
        if !self.supports_column_rotation_by(i) {
            Err(crate::EvaluationKeyError::Unsupported {
                operation: crate::EvaluationOperation::ColumnRotation { step: i },
            }
            .into())
        } else {
            let exponent = self.rot_to_gk_exponent.get(&i).ok_or_else(|| {
                crate::EvaluationKeyError::InvalidRotationStep {
                    step: i,
                    min: 1,
                    max: self.par.degree() / 2 - 1,
                }
            })?;
            let gk = self
                .gk
                .get(exponent)
                .ok_or(crate::EvaluationKeyError::Missing {
                    component: crate::EvaluationKeyComponent::GaloisKey { element: *exponent },
                })?;
            let mut out = Ciphertext::zero(&ct.par);
            gk.relinearize_into(ct, &mut out)?;
            Ok(out)
        }
    }

    /// Reports whether the evaluation key supports oblivious expansion.
    #[must_use]
    pub fn supports_expansion(&self, level: usize) -> bool {
        if level == 0 {
            true
        } else if self.evaluation_key_level == self.par.moduli().len() {
            false
        } else {
            let mut ret = level <= self.par.degree().ilog2() as usize;
            for l in 0..level {
                ret &= self.gk.contains_key(&((self.par.degree() >> l) + 1));
            }
            ret
        }
    }

    /// Obliviously expands the ciphertext. Returns an error if this evaluation
    /// does not support expansion to level = ceil(log2(size)), or if the
    /// ciphertext does not have size 2. The output is a vector of `size`
    /// ciphertexts.
    pub fn expands(&self, ct: &Ciphertext, size: usize) -> Result<Vec<Ciphertext>> {
        self.validate_ciphertext(ct)?;
        if size == 0 {
            return Err(crate::EvaluationKeyError::InvalidExpansionSize {
                size,
                degree: self.par.degree(),
            }
            .into());
        }
        if size > self.par.degree() {
            return Err(crate::EvaluationKeyError::InvalidExpansionSize {
                size,
                degree: self.par.degree(),
            }
            .into());
        }

        let level = size.next_power_of_two().ilog2() as usize;
        if level == 0 {
            Ok(vec![ct.clone()])
        } else if self.supports_expansion(level) {
            let mut out = vec![Ciphertext::zero(&ct.par); 1 << level];
            out[0] = ct.clone();
            let mut sub = Ciphertext::zero(&ct.par);

            // We use the Oblivious expansion algorithm of
            // https://eprint.iacr.org/2019/1483.pdf
            for l in 0..level {
                let monomial = self
                    .monomials
                    .get(l)
                    .ok_or(crate::EvaluationKeyError::Missing {
                        component: crate::EvaluationKeyComponent::ExpansionMonomial { level: l },
                    })?;
                let element = (self.par.degree() >> l) + 1;
                let gk = self
                    .gk
                    .get(&element)
                    .ok_or(crate::EvaluationKeyError::Missing {
                        component: crate::EvaluationKeyComponent::GaloisKey { element },
                    })?;
                let step = 1 << l;
                let (low, high) = out.split_at_mut(step);
                for i in 0..step {
                    gk.relinearize_into(&low[i], &mut sub)?;
                    let j = step | i;
                    if j < size {
                        let target = &mut high[i];
                        target.clone_from(&low[i]);
                        *target -= &sub;
                        target[0] *= monomial;
                        target[1] *= monomial;
                    }
                    low[i] += &sub;
                }
            }
            out.truncate(size);
            Ok(out)
        } else {
            Err(crate::EvaluationKeyError::Unsupported {
                operation: crate::EvaluationOperation::Expansion { level },
            }
            .into())
        }
    }

    fn validate_ciphertext(&self, ct: &Ciphertext) -> Result<()> {
        ct.validate_for(&self.par)?;
        if ct.len() != 2 {
            return Err(crate::CiphertextError::InvalidPolynomialCount {
                operation: crate::CiphertextOperation::EvaluationKey,
                actual: ct.len(),
                expected: 2,
            }
            .into());
        }
        if ct.level != self.ciphertext_level {
            return Err(Error::InvalidLevel {
                level: ct.level,
                min_level: self.ciphertext_level,
                max_level: self.ciphertext_level,
            });
        }
        Ok(())
    }

    fn construct_rot_to_gk_exponent(par: &Arc<BfvParameters>) -> HashMap<usize, usize> {
        let mut m = HashMap::new();
        let q = Modulus::new(2 * par.degree() as u64).unwrap();
        for i in 1..par.degree() / 2 {
            let exp = q.pow(3, i as u64) as usize;
            m.insert(i, exp);
        }
        m
    }
}

impl FheParametrized for EvaluationKey {
    type Parameters = BfvParameters;
}

impl Serialize for EvaluationKey {
    fn to_bytes(&self) -> Vec<u8> {
        EvaluationKeyProto::from(self).encode_to_vec()
    }
}

impl DeserializeParametrized for EvaluationKey {
    type Error = Error;

    fn from_bytes(bytes: &[u8], par: &Arc<Self::Parameters>) -> Result<Self> {
        let gkp = Message::decode(bytes).map_err(|_| {
            Error::SerializationError(SerializationError::Decode {
                object: crate::SerializedObject::EvaluationKey,
            })
        })?;
        EvaluationKey::try_convert_from(&gkp, par)
    }
}

/// Builder for a leveled evaluation key from the secret key.
#[derive(Debug)]
pub struct EvaluationKeyBuilder {
    sk: SecretKey,
    ciphertext_level: usize,
    evaluation_key_level: usize,
    inner_sum: bool,
    row_rotation: bool,
    expansion_level: usize,
    column_rotation: HashSet<usize>,
    rot_to_gk_exponent: HashMap<usize, usize>,
}

impl Zeroize for EvaluationKeyBuilder {
    fn zeroize(&mut self) {
        self.sk.zeroize()
    }
}

impl ZeroizeOnDrop for EvaluationKeyBuilder {}

impl EvaluationKeyBuilder {
    /// Creates a new builder from the [`SecretKey`].
    pub fn new(sk: &SecretKey) -> Result<Self> {
        Ok(Self {
            sk: sk.clone(),
            ciphertext_level: 0,
            evaluation_key_level: 0,
            inner_sum: false,
            row_rotation: false,
            expansion_level: 0,
            column_rotation: HashSet::new(),
            rot_to_gk_exponent: EvaluationKey::construct_rot_to_gk_exponent(&sk.par),
        })
    }

    /// Creates a new builder from the [`SecretKey`], for operations on
    /// ciphertexts at level `ciphertext_level` using keys at level
    /// `evaluation_key_level`. This raises an error if the key level is larger
    /// than the ciphertext level, or if the ciphertext level is larger than the
    /// maximum level supported by these parameters.
    pub fn new_leveled(
        sk: &SecretKey,
        ciphertext_level: usize,
        evaluation_key_level: usize,
    ) -> Result<Self> {
        if ciphertext_level > sk.par.max_level() {
            return Err(Error::InvalidLevel {
                level: ciphertext_level,
                min_level: 0,
                max_level: sk.par.max_level(),
            });
        }
        if evaluation_key_level > ciphertext_level {
            return Err(Error::InvalidLevel {
                level: evaluation_key_level,
                min_level: 0,
                max_level: ciphertext_level,
            });
        }

        Ok(Self {
            sk: sk.clone(),
            ciphertext_level,
            evaluation_key_level,
            inner_sum: false,
            row_rotation: false,
            expansion_level: 0,
            column_rotation: HashSet::new(),
            rot_to_gk_exponent: EvaluationKey::construct_rot_to_gk_exponent(&sk.par),
        })
    }

    /// Allow expansion by this evaluation key.
    pub fn enable_expansion(&mut self, level: usize) -> Result<&mut Self> {
        let max_level = self.sk.par.degree().ilog2() as usize;
        if level > max_level {
            Err(Error::InvalidLevel {
                level,
                min_level: 0,
                max_level,
            })
        } else {
            self.expansion_level = level;
            Ok(self)
        }
    }

    /// Allow this evaluation key to compute homomorphic inner sums.
    pub fn enable_inner_sum(&mut self) -> Result<&mut Self> {
        self.inner_sum = true;
        Ok(self)
    }

    /// Allow this evaluation key to homomorphically rotate the plaintext rows.
    pub fn enable_row_rotation(&mut self) -> Result<&mut Self> {
        self.row_rotation = true;
        Ok(self)
    }

    /// Allow this evaluation key to homomorphically rotate the plaintext
    /// columns.
    pub fn enable_column_rotation(&mut self, i: usize) -> Result<&mut Self> {
        if let Some(exp) = self.rot_to_gk_exponent.get(&i) {
            self.column_rotation.insert(*exp);
            Ok(self)
        } else {
            Err(crate::EvaluationKeyError::InvalidRotationStep {
                step: i,
                min: 1,
                max: self.sk.par.degree() / 2 - 1,
            }
            .into())
        }
    }

    /// Build an [`EvaluationKey`] with the specified attributes.
    pub fn build<R: RngCore + CryptoRng>(&mut self, rng: &mut R) -> Result<EvaluationKey> {
        let mut ek = EvaluationKey {
            gk: HashMap::default(),
            par: self.sk.par.clone(),
            rot_to_gk_exponent: self.rot_to_gk_exponent.clone(),
            monomials: Vec::with_capacity(self.sk.par.degree().ilog2() as usize),
            ciphertext_level: self.ciphertext_level,
            evaluation_key_level: self.evaluation_key_level,
        };

        let mut indices = self.column_rotation.clone();

        if self.row_rotation {
            indices.insert(self.sk.par.degree() * 2 - 1);
        }

        if self.inner_sum {
            // Add the required indices to the set of indices
            indices.insert(self.sk.par.degree() * 2 - 1);
            let mut i = 1;
            while i < self.sk.par.degree() / 2 {
                let exponent =
                    ek.rot_to_gk_exponent
                        .get(&i)
                        .ok_or(crate::EvaluationKeyError::Missing {
                            component: crate::EvaluationKeyComponent::GaloisExponent { step: i },
                        })?;
                indices.insert(*exponent);
                i *= 2
            }
        }

        for l in 0..self.expansion_level {
            indices.insert((self.sk.par.degree() >> l) + 1);
        }

        let ciphertext_ctx = self.sk.par.context_at_level(self.ciphertext_level)?;
        for l in 0..self.sk.par.degree().ilog2() {
            let mut monomial = vec![0i64; self.sk.par.degree()];
            monomial[self.sk.par.degree() - (1 << l)] = -1;
            let monomial = Poly::<PowerBasis>::try_convert_from_public(
                &monomial,
                ciphertext_ctx,
                fhe_traits::VariableTime::new(fhe_traits::PublicData::assert_public()),
            )?;
            ek.monomials.push(monomial.into_ntt_shoup());
        }

        for index in indices {
            ek.gk.insert(
                index,
                GaloisKey::new(
                    &self.sk,
                    index,
                    self.ciphertext_level,
                    self.evaluation_key_level,
                    rng,
                )?,
            );
        }

        Ok(ek)
    }
}

impl From<&EvaluationKey> for EvaluationKeyProto {
    fn from(ek: &EvaluationKey) -> Self {
        let mut proto = EvaluationKeyProto::default();
        for gk in ek.gk.values() {
            proto.gk.push(GaloisKeyProto::from(gk))
        }
        proto.ciphertext_level = ek.ciphertext_level as u32;
        proto.evaluation_key_level = ek.evaluation_key_level as u32;
        proto
    }
}

impl TryConvertFrom<&EvaluationKeyProto> for EvaluationKey {
    fn try_convert_from(value: &EvaluationKeyProto, par: &Arc<BfvParameters>) -> Result<Self> {
        let mut gk = HashMap::new();
        for gkp in &value.gk {
            let key = GaloisKey::try_convert_from(gkp, par)?;
            if key.ksk.ciphertext_level != value.ciphertext_level as usize {
                return Err(Error::InvalidLevel {
                    level: key.ksk.ciphertext_level,
                    min_level: value.ciphertext_level as usize,
                    max_level: value.ciphertext_level as usize,
                });
            }
            if key.ksk.ksk_level != value.evaluation_key_level as usize {
                return Err(Error::InvalidLevel {
                    level: key.ksk.ksk_level,
                    min_level: value.evaluation_key_level as usize,
                    max_level: value.evaluation_key_level as usize,
                });
            }
            gk.insert(key.element.exponent, key);
        }

        let ciphertext_ctx = par.context_at_level(value.ciphertext_level as usize)?;
        let mut monomials = Vec::with_capacity(par.degree().ilog2() as usize);
        for l in 0..par.degree().ilog2() {
            let mut monomial = vec![0i64; par.degree()];
            monomial[par.degree() - (1 << l)] = -1;
            let monomial = Poly::<PowerBasis>::try_convert_from_public(
                &monomial,
                ciphertext_ctx,
                fhe_traits::VariableTime::new(fhe_traits::PublicData::assert_public()),
            )?;
            monomials.push(monomial.into_ntt_shoup());
        }

        Ok(EvaluationKey {
            gk,
            par: par.clone(),
            rot_to_gk_exponent: EvaluationKey::construct_rot_to_gk_exponent(par),
            monomials,
            ciphertext_level: value.ciphertext_level as usize,
            evaluation_key_level: value.evaluation_key_level as usize,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::{EvaluationKey, EvaluationKeyBuilder};
    use crate::bfv::{BfvParameters, Encoding, Plaintext, SecretKey, traits::TryConvertFrom};
    use crate::proto::bfv::EvaluationKey as LeveledEvaluationKeyProto;
    use fhe_traits::{
        DeserializeParametrized, FheDecoder, FheDecrypter, FheEncoder, FheEncrypter, Serialize,
    };
    use itertools::izip;
    use rand::rng;
    use std::{cmp::min, error::Error};

    #[test]
    fn builder() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let params = BfvParameters::default_arc(6, 16);
        let sk = SecretKey::random(&params, &mut rng);

        let max_level = params.max_level();
        for ciphertext_level in 0..=max_level {
            for evaluation_key_level in 0..=min(max_level, ciphertext_level) {
                let mut builder =
                    EvaluationKeyBuilder::new_leveled(&sk, ciphertext_level, evaluation_key_level)?;

                assert!(!builder.build(&mut rng)?.supports_row_rotation());
                assert!(!builder.build(&mut rng)?.supports_column_rotation_by(0));
                assert!(!builder.build(&mut rng)?.supports_column_rotation_by(1));
                assert!(!builder.build(&mut rng)?.supports_inner_sum());
                assert!(!builder.build(&mut rng)?.supports_expansion(1));
                assert!(builder.build(&mut rng)?.supports_expansion(0));
                assert!(builder.enable_column_rotation(0).is_err());
                assert!(
                    builder
                        .enable_expansion(64 - params.degree().leading_zeros() as usize)
                        .is_err()
                );

                builder.enable_column_rotation(1)?;
                assert!(builder.build(&mut rng)?.supports_column_rotation_by(1));
                assert!(!builder.build(&mut rng)?.supports_row_rotation());
                assert!(!builder.build(&mut rng)?.supports_inner_sum());
                assert!(!builder.build(&mut rng)?.supports_expansion(1));

                builder.enable_row_rotation()?;
                assert!(builder.build(&mut rng)?.supports_row_rotation());
                assert!(!builder.build(&mut rng)?.supports_inner_sum());
                assert!(!builder.build(&mut rng)?.supports_expansion(1));

                builder.enable_inner_sum()?;
                assert!(builder.build(&mut rng)?.supports_inner_sum());
                assert!(builder.build(&mut rng)?.supports_expansion(1));
                assert!(
                    !builder
                        .build(&mut rng)?
                        .supports_expansion(64 - 1 - params.degree().leading_zeros() as usize)
                );

                builder.enable_expansion(64 - 1 - params.degree().leading_zeros() as usize)?;
                assert!(
                    builder
                        .build(&mut rng)?
                        .supports_expansion(64 - 1 - params.degree().leading_zeros() as usize)
                );

                assert!(builder.build(&mut rng).is_ok());

                // Enabling inner sum enables row rotation and a few column rotations :)
                let ek = EvaluationKeyBuilder::new_leveled(&sk, 0, 0)?
                    .enable_inner_sum()?
                    .build(&mut rng)?;
                assert!(ek.supports_inner_sum());
                assert!(ek.supports_row_rotation());
                let mut i = 1;
                while i < params.degree() / 2 {
                    assert!(ek.supports_column_rotation_by(i));
                    i *= 2
                }
                assert!(!ek.supports_column_rotation_by(params.degree() / 2 - 1));
            }
        }

        let e = EvaluationKeyBuilder::new_leveled(&sk, 0, 1);
        assert!(e.is_err());
        assert_eq!(
            e.unwrap_err(),
            crate::Error::InvalidLevel {
                level: 1,
                min_level: 0,
                max_level: 0,
            }
        );

        Ok(())
    }

    #[test]
    fn inner_sum() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for params in [
            BfvParameters::default_arc(6, 16),
            BfvParameters::default_arc(5, 16),
        ] {
            for _ in 0..25 {
                for ciphertext_level in 0..=params.max_level() {
                    for evaluation_key_level in 0..=min(params.max_level() - 1, ciphertext_level) {
                        let sk = SecretKey::random(&params, &mut rng);
                        let ek = EvaluationKeyBuilder::new_leveled(
                            &sk,
                            ciphertext_level,
                            evaluation_key_level,
                        )?
                        .enable_inner_sum()?
                        .build(&mut rng)?;

                        let v = fhe_math::zq::Modulus::new(params.plaintext())
                            .unwrap()
                            .random_vec(params.degree(), &mut rng);
                        let expected = fhe_math::zq::Modulus::new(params.plaintext())
                            .unwrap()
                            .reduce_u128(v.iter().map(|vi| *vi as u128).sum());

                        let pt = Plaintext::try_encode(
                            &v,
                            Encoding::simd_at_level(ciphertext_level),
                            &params,
                        )?;
                        let ct = sk.try_encrypt(&pt, &mut rng)?;

                        let ct2 = ek.computes_inner_sum(&ct)?;
                        let pt = sk.try_decrypt(&ct2)?;
                        assert_eq!(
                            Vec::<u64>::try_decode(&pt, Encoding::simd_at_level(ciphertext_level))?,
                            vec![expected; params.degree()]
                        )
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn row_rotation() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for params in [
            BfvParameters::default_arc(6, 16),
            BfvParameters::default_arc(5, 16),
        ] {
            for _ in 0..50 {
                for ciphertext_level in 0..=params.max_level() {
                    for evaluation_key_level in 0..=min(params.max_level() - 1, ciphertext_level) {
                        let sk = SecretKey::random(&params, &mut rng);
                        let ek = EvaluationKeyBuilder::new_leveled(
                            &sk,
                            ciphertext_level,
                            evaluation_key_level,
                        )?
                        .enable_row_rotation()?
                        .build(&mut rng)?;

                        let v = fhe_math::zq::Modulus::new(params.plaintext())
                            .unwrap()
                            .random_vec(params.degree(), &mut rng);
                        let row_size = params.degree() >> 1;
                        let mut expected = vec![0u64; params.degree()];
                        expected[..row_size].copy_from_slice(&v[row_size..]);
                        expected[row_size..].copy_from_slice(&v[..row_size]);

                        let pt = Plaintext::try_encode(
                            &v,
                            Encoding::simd_at_level(ciphertext_level),
                            &params,
                        )?;
                        let ct = sk.try_encrypt(&pt, &mut rng)?;

                        let ct2 = ek.rotates_rows(&ct)?;
                        let pt = sk.try_decrypt(&ct2)?;
                        assert_eq!(
                            Vec::<u64>::try_decode(&pt, Encoding::simd_at_level(ciphertext_level))?,
                            expected
                        )
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn column_rotation() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for params in [
            BfvParameters::default_arc(6, 16),
            BfvParameters::default_arc(5, 16),
        ] {
            let row_size = params.degree() >> 1;
            for _ in 0..50 {
                for i in 1..row_size {
                    for ciphertext_level in 0..=params.max_level() {
                        for evaluation_key_level in 0..=min(params.max_level(), ciphertext_level) {
                            let sk = SecretKey::random(&params, &mut rng);
                            let ek = EvaluationKeyBuilder::new_leveled(
                                &sk,
                                ciphertext_level,
                                evaluation_key_level,
                            )?
                            .enable_column_rotation(i)?
                            .build(&mut rng)?;

                            let v = fhe_math::zq::Modulus::new(params.plaintext())
                                .unwrap()
                                .random_vec(params.degree(), &mut rng);
                            let row_size = params.degree() >> 1;
                            let mut expected = vec![0u64; params.degree()];
                            expected[..row_size - i].copy_from_slice(&v[i..row_size]);
                            expected[row_size - i..row_size].copy_from_slice(&v[..i]);
                            expected[row_size..2 * row_size - i]
                                .copy_from_slice(&v[row_size + i..]);
                            expected[2 * row_size - i..]
                                .copy_from_slice(&v[row_size..row_size + i]);

                            let pt = Plaintext::try_encode(
                                &v,
                                Encoding::simd_at_level(ciphertext_level),
                                &params,
                            )?;
                            let ct = sk.try_encrypt(&pt, &mut rng)?;

                            let ct2 = ek.rotates_columns_by(&ct, i)?;
                            let pt = sk.try_decrypt(&ct2)?;
                            assert_eq!(
                                Vec::<u64>::try_decode(
                                    &pt,
                                    Encoding::simd_at_level(ciphertext_level)
                                )?,
                                expected
                            )
                        }
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn expansion() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for params in [
            BfvParameters::default_arc(6, 16),
            BfvParameters::default_arc(5, 16),
        ] {
            let log_degree = 64 - 1 - params.degree().leading_zeros();
            for _ in 0..15 {
                for i in 1..1 + log_degree as usize {
                    for ciphertext_level in 0..=params.max_level() {
                        for evaluation_key_level in 0..=min(params.max_level(), ciphertext_level) {
                            let sk = SecretKey::random(&params, &mut rng);
                            let ek = EvaluationKeyBuilder::new_leveled(
                                &sk,
                                ciphertext_level,
                                evaluation_key_level,
                            )?
                            .enable_expansion(i)?
                            .build(&mut rng)?;

                            assert!(ek.supports_expansion(i));
                            assert!(!ek.supports_expansion(i + 1));
                            let v = fhe_math::zq::Modulus::new(params.plaintext())
                                .unwrap()
                                .random_vec(1 << i, &mut rng);
                            let pt = Plaintext::try_encode(
                                &v,
                                Encoding::poly_at_level(ciphertext_level),
                                &params,
                            )?;
                            let ct = sk.try_encrypt(&pt, &mut rng)?;

                            let ct2 = ek.expands(&ct, 1 << i)?;
                            assert_eq!(ct2.len(), 1 << i);
                            for (vi, ct2i) in izip!(&v, &ct2) {
                                let mut expected = vec![0u64; params.degree()];
                                expected[0] = fhe_math::zq::Modulus::new(params.plaintext())
                                    .unwrap()
                                    .mul(*vi, (1 << i) as u64);
                                let pt = sk.try_decrypt(ct2i)?;
                                assert_eq!(
                                    expected,
                                    Vec::<u64>::try_decode(
                                        &pt,
                                        Encoding::poly_at_level(ciphertext_level)
                                    )?
                                );
                                println!("Noise: {:?}", unsafe { sk.measure_noise(ct2i) })
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn expansion_rejects_invalid_sizes() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let params = BfvParameters::default_arc(3, 16);
        let sk = SecretKey::random(&params, &mut rng);
        let ek = EvaluationKeyBuilder::new(&sk)?
            .enable_expansion(1)?
            .build(&mut rng)?;
        let pt = Plaintext::try_encode(&[1u64][..], Encoding::poly(), &params)?;
        let ct = sk.try_encrypt(&pt, &mut rng)?;

        assert_eq!(
            ek.expands(&ct, 0),
            Err(crate::Error::EvaluationKey(
                crate::EvaluationKeyError::InvalidExpansionSize {
                    size: 0,
                    degree: params.degree(),
                }
            ))
        );
        assert_eq!(
            ek.expands(&ct, params.degree() + 1),
            Err(crate::Error::EvaluationKey(
                crate::EvaluationKeyError::InvalidExpansionSize {
                    size: params.degree() + 1,
                    degree: params.degree(),
                }
            ))
        );
        Ok(())
    }

    #[test]
    fn proto_conversion() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for params in [
            BfvParameters::default_arc(1, 16),
            BfvParameters::default_arc(6, 16),
            BfvParameters::default_arc(5, 16),
        ] {
            let sk = SecretKey::random(&params, &mut rng);

            let ek = EvaluationKeyBuilder::new_leveled(&sk, 0, 0)?.build(&mut rng)?;

            let proto = LeveledEvaluationKeyProto::from(&ek);
            assert_eq!(ek, EvaluationKey::try_convert_from(&proto, &params)?);

            let ek = EvaluationKeyBuilder::new_leveled(&sk, 0, 0)?
                .enable_row_rotation()?
                .build(&mut rng)?;

            let proto = LeveledEvaluationKeyProto::from(&ek);
            assert_eq!(ek, EvaluationKey::try_convert_from(&proto, &params)?);

            let ek = EvaluationKeyBuilder::new_leveled(&sk, 0, 0)?
                .enable_inner_sum()?
                .build(&mut rng)?;
            let proto = LeveledEvaluationKeyProto::from(&ek);
            assert_eq!(ek, EvaluationKey::try_convert_from(&proto, &params)?);

            let ek = EvaluationKeyBuilder::new_leveled(&sk, 0, 0)?
                .enable_expansion(params.degree().ilog2() as usize)?
                .build(&mut rng)?;
            let proto = LeveledEvaluationKeyProto::from(&ek);
            assert_eq!(ek, EvaluationKey::try_convert_from(&proto, &params)?);

            let ek = EvaluationKeyBuilder::new_leveled(&sk, 0, 0)?
                .enable_inner_sum()?
                .enable_expansion(params.degree().ilog2() as usize)?
                .build(&mut rng)?;
            let proto = LeveledEvaluationKeyProto::from(&ek);
            assert_eq!(ek, EvaluationKey::try_convert_from(&proto, &params)?);
        }
        Ok(())
    }

    #[test]
    fn serialize() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for params in [
            BfvParameters::default_arc(1, 16),
            BfvParameters::default_arc(6, 16),
        ] {
            let sk = SecretKey::random(&params, &mut rng);

            let ek = EvaluationKeyBuilder::new_leveled(&sk, 0, 0)?.build(&mut rng)?;
            let bytes = ek.to_bytes();
            assert_eq!(ek, EvaluationKey::from_bytes(&bytes, &params)?);

            if params.moduli.len() > 1 {
                let ek = EvaluationKeyBuilder::new_leveled(&sk, 0, 0)?
                    .enable_row_rotation()?
                    .build(&mut rng)?;
                let bytes = ek.to_bytes();
                assert_eq!(ek, EvaluationKey::from_bytes(&bytes, &params)?);

                let ek = EvaluationKeyBuilder::new_leveled(&sk, 0, 0)?
                    .enable_inner_sum()?
                    .build(&mut rng)?;
                let bytes = ek.to_bytes();
                assert_eq!(ek, EvaluationKey::from_bytes(&bytes, &params)?);

                let ek = EvaluationKeyBuilder::new_leveled(&sk, 0, 0)?
                    .enable_expansion(params.degree().ilog2() as usize)?
                    .build(&mut rng)?;
                let bytes = ek.to_bytes();
                assert_eq!(ek, EvaluationKey::from_bytes(&bytes, &params)?);

                let ek = EvaluationKeyBuilder::new_leveled(&sk, 0, 0)?
                    .enable_inner_sum()?
                    .enable_expansion(params.degree().ilog2() as usize)?
                    .build(&mut rng)?;
                let bytes = ek.to_bytes();
                assert_eq!(ek, EvaluationKey::from_bytes(&bytes, &params)?);
            }
        }
        Ok(())
    }
}
