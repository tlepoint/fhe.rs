//! Leveled evaluation keys for the BFV encryption scheme.

use crate::bfv::{Ciphertext, Parameters, SecretKey, keys::GaloisKey, wire::FromProto};
use crate::proto::bfv::{EvaluationKey as EvaluationKeyProto, GaloisKey as GaloisKeyProto};
use crate::{Error, Result, SerializationError};
use fhe_math::rq::{NttShoup, Poly, PowerBasis};
use fhe_math::zq::Modulus;

use prost::Message;
use rand::{CryptoRng, Rng as RngCore};
use std::collections::{HashMap, HashSet};

/// Evaluation key for the BFV encryption scheme.
///
/// An evaluation key enables one or several of the following operations:
/// - column rotation
/// - row rotation
/// - oblivious expansion
/// - inner sum
#[derive(Debug, PartialEq, Eq)]
pub struct EvaluationKey {
    par: Parameters,

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
    /// Start configuring the evaluation operations needed by an application.
    #[must_use]
    pub fn builder(sk: &SecretKey) -> EvaluationKeyBuilder<'_> {
        EvaluationKeyBuilder::new(sk)
    }

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
    pub fn inner_sum(&self, ct: &Ciphertext) -> Result<Ciphertext> {
        self.validate_ciphertext(ct)?;
        if !self.supports_inner_sum() {
            Err(crate::EvaluationKeyError::Unsupported {
                operation: crate::EvaluationOperation::InnerSum,
            }
            .into())
        } else {
            let mut out = ct.clone();
            let mut tmp = Ciphertext::trivial_zero(&ct.par, ct.level)?;

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
                out.add_assign(&tmp)?;
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
            out.add_assign(&tmp)?;

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
    pub fn rotate_rows(&self, ct: &Ciphertext) -> Result<Ciphertext> {
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
            let mut out = Ciphertext::trivial_zero(&ct.par, ct.level)?;
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
    pub fn rotate_columns(&self, ct: &Ciphertext, i: usize) -> Result<Ciphertext> {
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
            let mut out = Ciphertext::trivial_zero(&ct.par, ct.level)?;
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

    /// Obliviously expand the ciphertext. Returns an error if this evaluation
    /// does not support expansion to level = ceil(log2(size)), or if the
    /// ciphertext does not have size 2. The output is a vector of `size`
    /// ciphertexts.
    pub fn expand(&self, ct: &Ciphertext, size: usize) -> Result<Vec<Ciphertext>> {
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
            let mut out = Vec::with_capacity(size);
            out.push(ct.clone());
            let mut sub = Ciphertext::trivial_zero(&ct.par, ct.level)?;

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
                let step = out.len();
                for i in 0..step {
                    gk.relinearize_into(&out[i], &mut sub)?;
                    if step + i < size {
                        let mut target = out[i].subtract(&sub)?;
                        target.c[0] *= monomial;
                        target.c[1] *= monomial;
                        out.push(target);
                    }
                    out[i].add_assign(&sub)?;
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

    fn construct_rot_to_gk_exponent(par: &Parameters) -> HashMap<usize, usize> {
        let mut m = HashMap::new();
        let q = Modulus::new(2 * par.degree() as u64).unwrap();
        for i in 1..par.degree() / 2 {
            let exp = q.pow(3, i as u64) as usize;
            m.insert(i, exp);
        }
        m
    }
}

impl EvaluationKey {
    /// Serialize in the existing protobuf wire format.
    #[must_use]
    pub fn to_bytes(&self) -> Vec<u8> {
        EvaluationKeyProto::from(self).encode_to_vec()
    }
}

impl EvaluationKey {
    /// Import validated protobuf bytes, binding contextual values to the
    /// supplied parameters.
    pub fn from_bytes(bytes: &[u8], par: &Parameters) -> Result<Self> {
        let gkp = Message::decode(bytes).map_err(|_| {
            Error::SerializationError(SerializationError::Decode {
                object: crate::SerializedObject::EvaluationKey,
            })
        })?;
        EvaluationKey::from_proto(&gkp, par)
    }
}

/// Builder for a leveled evaluation key borrowing the secret key.
/// The builder retains no owned copy of the secret coefficients.
///
/// ```compile_fail
/// use fhe::bfv::{evaluation::EvaluationKeyBuilder, SecretKey};
/// fn build(sk: SecretKey) {
///     let mut builder = EvaluationKeyBuilder::new(&sk);
///     drop(sk);
///     builder.build(&mut rand::rng());
/// }
/// ```
#[derive(Clone, Debug)]
pub struct EvaluationKeyBuilder<'a> {
    sk: &'a SecretKey,
    ciphertext_level: usize,
    evaluation_key_level: usize,
    inner_sum: bool,
    row_rotation: bool,
    expansion_level: usize,
    column_rotation: HashSet<usize>,
    rot_to_gk_exponent: HashMap<usize, usize>,
}

impl<'a> EvaluationKeyBuilder<'a> {
    /// Create a builder borrowing a secret key. Configuration is validated by
    /// `build`.
    #[must_use]
    pub fn new(sk: &'a SecretKey) -> Self {
        Self {
            sk,
            ciphertext_level: 0,
            evaluation_key_level: 0,
            inner_sum: false,
            row_rotation: false,
            expansion_level: 0,
            column_rotation: HashSet::new(),
            rot_to_gk_exponent: EvaluationKey::construct_rot_to_gk_exponent(&sk.par),
        }
    }
    /// Select the level of ciphertexts that will be evaluated.
    #[must_use]
    pub fn ciphertext_level(mut self, level: usize) -> Self {
        self.ciphertext_level = level;
        self
    }
    /// Select the evaluation-key level, which must not exceed the ciphertext
    /// level.
    #[must_use]
    pub fn key_level(mut self, level: usize) -> Self {
        self.evaluation_key_level = level;
        self
    }
    /// Include keys for expansion up to `log_size` binary expansion steps.
    #[must_use]
    pub fn enable_expansion(mut self, log_size: usize) -> Self {
        self.expansion_level = log_size;
        self
    }
    /// Include keys for homomorphic inner sums.
    #[must_use]
    pub fn enable_inner_sum(mut self) -> Self {
        self.inner_sum = true;
        self
    }
    /// Include a key for swapping the two SIMD rows.
    #[must_use]
    pub fn enable_row_rotation(mut self) -> Self {
        self.row_rotation = true;
        self
    }
    /// Include a column rotation step, between one and `degree / 2 - 1`.
    #[must_use]
    pub fn enable_column_rotation(mut self, step: usize) -> Self {
        self.column_rotation.insert(step);
        self
    }

    /// Build an [`EvaluationKey`] with the specified attributes.
    pub fn build<R: RngCore + CryptoRng>(self, rng: &mut R) -> Result<EvaluationKey> {
        self.sk.par.context_at_level(self.ciphertext_level)?;
        if self.evaluation_key_level > self.ciphertext_level {
            return Err(Error::InvalidLevel {
                level: self.evaluation_key_level,
                min_level: 0,
                max_level: self.ciphertext_level,
            });
        }
        let max_expansion = self.sk.par.degree().ilog2() as usize;
        if self.expansion_level > max_expansion {
            return Err(Error::InvalidLevel {
                level: self.expansion_level,
                min_level: 0,
                max_level: max_expansion,
            });
        }
        let mut indices = self
            .column_rotation
            .iter()
            .map(|step| {
                self.rot_to_gk_exponent.get(step).copied().ok_or_else(|| {
                    crate::EvaluationKeyError::InvalidRotationStep {
                        step: *step,
                        min: 1,
                        max: self.sk.par.degree() / 2 - 1,
                    }
                    .into()
                })
            })
            .collect::<Result<HashSet<_>>>()?;
        let mut ek = EvaluationKey {
            gk: HashMap::default(),
            par: self.sk.par.clone(),
            rot_to_gk_exponent: self.rot_to_gk_exponent,
            monomials: Vec::with_capacity(self.sk.par.degree().ilog2() as usize),
            ciphertext_level: self.ciphertext_level,
            evaluation_key_level: self.evaluation_key_level,
        };

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
            let monomial = Poly::<PowerBasis>::from_signed_coefficients_with_timing(
                &monomial,
                ciphertext_ctx,
                Some(crate::VariableTime::new(crate::PublicData::assert_public())),
            )?;
            ek.monomials.push(monomial.into_ntt_shoup());
        }

        for index in indices {
            ek.gk.insert(
                index,
                GaloisKey::new(
                    self.sk,
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

impl FromProto<&EvaluationKeyProto> for EvaluationKey {
    fn from_proto(value: &EvaluationKeyProto, par: &Parameters) -> Result<Self> {
        let mut gk = HashMap::new();
        for gkp in &value.gk {
            let key = GaloisKey::from_proto(gkp, par)?;
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
            gk.insert(key.element.exponent(), key);
        }

        let ciphertext_ctx = par.context_at_level(value.ciphertext_level as usize)?;
        let mut monomials = Vec::with_capacity(par.degree().ilog2() as usize);
        for l in 0..par.degree().ilog2() {
            let mut monomial = vec![0i64; par.degree()];
            monomial[par.degree() - (1 << l)] = -1;
            let monomial = Poly::<PowerBasis>::from_signed_coefficients_with_timing(
                &monomial,
                ciphertext_ctx,
                Some(crate::VariableTime::new(crate::PublicData::assert_public())),
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
    use crate::bfv::{Encoding, Parameters, Plaintext, SecretKey, wire::FromProto};
    use crate::proto::bfv::EvaluationKey as LeveledEvaluationKeyProto;

    use itertools::izip;
    use rand::rng;
    use std::{cmp::min, error::Error};

    #[test]
    fn partial_expansion_matches_full_prefix_bytes_and_permissions() -> Result<(), Box<dyn Error>> {
        let params = Parameters::test_parameters(2, 16);
        let mut rng = rng();
        let sk = SecretKey::generate(&params, &mut rng);
        for level in 0..=params.max_level() {
            let key = EvaluationKeyBuilder::new(&sk)
                .ciphertext_level(level)
                .key_level(0)
                .enable_expansion(4)
                .build(&mut rng)?;
            let pt =
                Plaintext::encode_at_level(&params, &[1_u64, 2, 3], Encoding::Polynomial, level)?;
            let ct: crate::bfv::Ciphertext = sk.encrypt(&pt, &mut rng)?;
            for restricted in [false, true] {
                let mut input = ct.clone();
                if restricted {
                    for part in input.iter_mut() {
                        part.disallow_variable_time_computations();
                    }
                }
                for size in 1..=params.degree() {
                    let full = key.expand(&input, size.next_power_of_two())?;
                    let partial = key.expand(&input, size)?;
                    assert_eq!(partial.len(), size);
                    for (actual, expected) in partial.iter().zip(&full) {
                        assert_eq!(actual.to_bytes(), expected.to_bytes());
                        assert_eq!(actual.level(), level);
                        assert!(
                            actual
                                .iter()
                                .all(|p| p.allows_variable_time_computations() != restricted)
                        );
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn builder() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let params = Parameters::test_parameters(6, 16);
        let sk = SecretKey::generate(&params, &mut rng);

        let max_level = params.max_level();
        for ciphertext_level in 0..=max_level {
            for evaluation_key_level in 0..=min(max_level, ciphertext_level) {
                let mut builder = EvaluationKeyBuilder::new(&sk)
                    .ciphertext_level(ciphertext_level)
                    .key_level(evaluation_key_level);

                assert!(!builder.clone().build(&mut rng)?.supports_row_rotation());
                assert!(
                    !builder
                        .clone()
                        .build(&mut rng)?
                        .supports_column_rotation_by(0)
                );
                assert!(
                    !builder
                        .clone()
                        .build(&mut rng)?
                        .supports_column_rotation_by(1)
                );
                assert!(!builder.clone().build(&mut rng)?.supports_inner_sum());
                assert!(!builder.clone().build(&mut rng)?.supports_expansion(1));
                assert!(builder.clone().build(&mut rng)?.supports_expansion(0));
                assert!(
                    builder
                        .clone()
                        .enable_column_rotation(0)
                        .build(&mut rng)
                        .is_err()
                );
                assert!(
                    builder
                        .clone()
                        .enable_expansion(64 - params.degree().leading_zeros() as usize)
                        .build(&mut rng)
                        .is_err()
                );

                builder = builder.enable_column_rotation(1);
                assert!(
                    builder
                        .clone()
                        .build(&mut rng)?
                        .supports_column_rotation_by(1)
                );
                assert!(!builder.clone().build(&mut rng)?.supports_row_rotation());
                assert!(!builder.clone().build(&mut rng)?.supports_inner_sum());
                assert!(!builder.clone().build(&mut rng)?.supports_expansion(1));

                builder = builder.enable_row_rotation();
                assert!(builder.clone().build(&mut rng)?.supports_row_rotation());
                assert!(!builder.clone().build(&mut rng)?.supports_inner_sum());
                assert!(!builder.clone().build(&mut rng)?.supports_expansion(1));

                builder = builder.enable_inner_sum();
                assert!(builder.clone().build(&mut rng)?.supports_inner_sum());
                assert!(builder.clone().build(&mut rng)?.supports_expansion(1));
                assert!(
                    !builder
                        .clone()
                        .build(&mut rng)?
                        .supports_expansion(64 - 1 - params.degree().leading_zeros() as usize)
                );

                builder =
                    builder.enable_expansion(64 - 1 - params.degree().leading_zeros() as usize);
                assert!(
                    builder
                        .clone()
                        .build(&mut rng)?
                        .supports_expansion(64 - 1 - params.degree().leading_zeros() as usize)
                );

                assert!(builder.clone().build(&mut rng).is_ok());

                // Enabling inner sum enables row rotation and a few column rotations :)
                let ek = EvaluationKeyBuilder::new(&sk)
                    .ciphertext_level(0)
                    .key_level(0)
                    .enable_inner_sum()
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

        let e = EvaluationKeyBuilder::new(&sk)
            .ciphertext_level(0)
            .key_level(1)
            .build(&mut rng);
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
            Parameters::test_parameters(6, 16),
            Parameters::test_parameters(5, 16),
        ] {
            for _ in 0..25 {
                for ciphertext_level in 0..=params.max_level() {
                    for evaluation_key_level in 0..=min(params.max_level() - 1, ciphertext_level) {
                        let sk = SecretKey::generate(&params, &mut rng);
                        let ek = EvaluationKeyBuilder::new(&sk)
                            .ciphertext_level(ciphertext_level)
                            .key_level(evaluation_key_level)
                            .enable_inner_sum()
                            .build(&mut rng)?;

                        let v = fhe_math::zq::Modulus::new(params.plaintext_modulus_u64().unwrap())
                            .unwrap()
                            .random_vec(params.degree(), &mut rng);
                        let expected =
                            fhe_math::zq::Modulus::new(params.plaintext_modulus_u64().unwrap())
                                .unwrap()
                                .reduce_u128(v.iter().map(|vi| *vi as u128).sum());

                        let pt = Plaintext::encode_at_level(
                            &params,
                            &v,
                            Encoding::Simd,
                            ciphertext_level,
                        )?;
                        let ct = sk.encrypt(&pt, &mut rng)?;

                        let ct2 = ek.inner_sum(&ct)?;
                        let pt = sk.decrypt(&ct2)?;
                        assert_eq!(pt.decode(Encoding::Simd)?, vec![expected; params.degree()])
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
            Parameters::test_parameters(6, 16),
            Parameters::test_parameters(5, 16),
        ] {
            for _ in 0..50 {
                for ciphertext_level in 0..=params.max_level() {
                    for evaluation_key_level in 0..=min(params.max_level() - 1, ciphertext_level) {
                        let sk = SecretKey::generate(&params, &mut rng);
                        let ek = EvaluationKeyBuilder::new(&sk)
                            .ciphertext_level(ciphertext_level)
                            .key_level(evaluation_key_level)
                            .enable_row_rotation()
                            .build(&mut rng)?;

                        let v = fhe_math::zq::Modulus::new(params.plaintext_modulus_u64().unwrap())
                            .unwrap()
                            .random_vec(params.degree(), &mut rng);
                        let row_size = params.degree() >> 1;
                        let mut expected = vec![0u64; params.degree()];
                        expected[..row_size].copy_from_slice(&v[row_size..]);
                        expected[row_size..].copy_from_slice(&v[..row_size]);

                        let pt = Plaintext::encode_at_level(
                            &params,
                            &v,
                            Encoding::Simd,
                            ciphertext_level,
                        )?;
                        let ct = sk.encrypt(&pt, &mut rng)?;

                        let ct2 = ek.rotate_rows(&ct)?;
                        let pt = sk.decrypt(&ct2)?;
                        assert_eq!(pt.decode(Encoding::Simd)?, expected)
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
            Parameters::test_parameters(6, 16),
            Parameters::test_parameters(5, 16),
        ] {
            let row_size = params.degree() >> 1;
            for _ in 0..50 {
                for i in 1..row_size {
                    for ciphertext_level in 0..=params.max_level() {
                        for evaluation_key_level in 0..=min(params.max_level(), ciphertext_level) {
                            let sk = SecretKey::generate(&params, &mut rng);
                            let ek = EvaluationKeyBuilder::new(&sk)
                                .ciphertext_level(ciphertext_level)
                                .key_level(evaluation_key_level)
                                .enable_column_rotation(i)
                                .build(&mut rng)?;

                            let v =
                                fhe_math::zq::Modulus::new(params.plaintext_modulus_u64().unwrap())
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

                            let pt = Plaintext::encode_at_level(
                                &params,
                                &v,
                                Encoding::Simd,
                                ciphertext_level,
                            )?;
                            let ct = sk.encrypt(&pt, &mut rng)?;

                            let ct2 = ek.rotate_columns(&ct, i)?;
                            let pt = sk.decrypt(&ct2)?;
                            assert_eq!(pt.decode(Encoding::Simd)?, expected)
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
            Parameters::test_parameters(6, 16),
            Parameters::test_parameters(5, 16),
        ] {
            let log_degree = 64 - 1 - params.degree().leading_zeros();
            for _ in 0..15 {
                for i in 1..1 + log_degree as usize {
                    for ciphertext_level in 0..=params.max_level() {
                        for evaluation_key_level in 0..=min(params.max_level(), ciphertext_level) {
                            let sk = SecretKey::generate(&params, &mut rng);
                            let ek = EvaluationKeyBuilder::new(&sk)
                                .ciphertext_level(ciphertext_level)
                                .key_level(evaluation_key_level)
                                .enable_expansion(i)
                                .build(&mut rng)?;

                            assert!(ek.supports_expansion(i));
                            assert!(!ek.supports_expansion(i + 1));
                            let v =
                                fhe_math::zq::Modulus::new(params.plaintext_modulus_u64().unwrap())
                                    .unwrap()
                                    .random_vec(1 << i, &mut rng);
                            let pt = Plaintext::encode_at_level(
                                &params,
                                &v,
                                Encoding::Polynomial,
                                ciphertext_level,
                            )?;
                            let ct = sk.encrypt(&pt, &mut rng)?;

                            let ct2 = ek.expand(&ct, 1 << i)?;
                            assert_eq!(ct2.len(), 1 << i);
                            for (vi, ct2i) in izip!(&v, &ct2) {
                                let mut expected = vec![0u64; params.degree()];
                                expected[0] = fhe_math::zq::Modulus::new(
                                    params.plaintext_modulus_u64().unwrap(),
                                )
                                .unwrap()
                                .mul(*vi, (1 << i) as u64);
                                let pt = sk.decrypt(ct2i)?;
                                assert_eq!(expected, pt.decode(Encoding::Polynomial)?);
                                println!(
                                    "Noise: {:?}",
                                    sk.measure_noise_vartime(
                                        ct2i,
                                        crate::SecretDependentDiagnostics::acknowledge_leakage()
                                    )
                                )
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
        let params = Parameters::test_parameters(3, 16);
        let sk = SecretKey::generate(&params, &mut rng);
        let ek = EvaluationKeyBuilder::new(&sk)
            .enable_expansion(1)
            .build(&mut rng)?;
        let pt = Plaintext::encode(&params, &[1u64][..], Encoding::Polynomial)?;
        let ct = sk.encrypt(&pt, &mut rng)?;

        assert_eq!(
            ek.expand(&ct, 0),
            Err(crate::Error::EvaluationKey(
                crate::EvaluationKeyError::InvalidExpansionSize {
                    size: 0,
                    degree: params.degree(),
                }
            ))
        );
        assert_eq!(
            ek.expand(&ct, params.degree() + 1),
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
            Parameters::test_parameters(1, 16),
            Parameters::test_parameters(6, 16),
            Parameters::test_parameters(5, 16),
        ] {
            let sk = SecretKey::generate(&params, &mut rng);

            let ek = EvaluationKeyBuilder::new(&sk)
                .ciphertext_level(0)
                .key_level(0)
                .build(&mut rng)?;

            let proto = LeveledEvaluationKeyProto::from(&ek);
            assert_eq!(ek, EvaluationKey::from_proto(&proto, &params)?);

            let ek = EvaluationKeyBuilder::new(&sk)
                .ciphertext_level(0)
                .key_level(0)
                .enable_row_rotation()
                .build(&mut rng)?;

            let proto = LeveledEvaluationKeyProto::from(&ek);
            assert_eq!(ek, EvaluationKey::from_proto(&proto, &params)?);

            let ek = EvaluationKeyBuilder::new(&sk)
                .ciphertext_level(0)
                .key_level(0)
                .enable_inner_sum()
                .build(&mut rng)?;
            let proto = LeveledEvaluationKeyProto::from(&ek);
            assert_eq!(ek, EvaluationKey::from_proto(&proto, &params)?);

            let ek = EvaluationKeyBuilder::new(&sk)
                .ciphertext_level(0)
                .key_level(0)
                .enable_expansion(params.degree().ilog2() as usize)
                .build(&mut rng)?;
            let proto = LeveledEvaluationKeyProto::from(&ek);
            assert_eq!(ek, EvaluationKey::from_proto(&proto, &params)?);

            let ek = EvaluationKeyBuilder::new(&sk)
                .ciphertext_level(0)
                .key_level(0)
                .enable_inner_sum()
                .enable_expansion(params.degree().ilog2() as usize)
                .build(&mut rng)?;
            let proto = LeveledEvaluationKeyProto::from(&ek);
            assert_eq!(ek, EvaluationKey::from_proto(&proto, &params)?);
        }
        Ok(())
    }

    #[test]
    fn serialize() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for params in [
            Parameters::test_parameters(1, 16),
            Parameters::test_parameters(6, 16),
        ] {
            let sk = SecretKey::generate(&params, &mut rng);

            let ek = EvaluationKeyBuilder::new(&sk)
                .ciphertext_level(0)
                .key_level(0)
                .build(&mut rng)?;
            let bytes = ek.to_bytes();
            assert_eq!(ek, EvaluationKey::from_bytes(&bytes, &params)?);

            if params.inner.moduli.len() > 1 {
                let ek = EvaluationKeyBuilder::new(&sk)
                    .ciphertext_level(0)
                    .key_level(0)
                    .enable_row_rotation()
                    .build(&mut rng)?;
                let bytes = ek.to_bytes();
                assert_eq!(ek, EvaluationKey::from_bytes(&bytes, &params)?);

                let ek = EvaluationKeyBuilder::new(&sk)
                    .ciphertext_level(0)
                    .key_level(0)
                    .enable_inner_sum()
                    .build(&mut rng)?;
                let bytes = ek.to_bytes();
                assert_eq!(ek, EvaluationKey::from_bytes(&bytes, &params)?);

                let ek = EvaluationKeyBuilder::new(&sk)
                    .ciphertext_level(0)
                    .key_level(0)
                    .enable_expansion(params.degree().ilog2() as usize)
                    .build(&mut rng)?;
                let bytes = ek.to_bytes();
                assert_eq!(ek, EvaluationKey::from_bytes(&bytes, &params)?);

                let ek = EvaluationKeyBuilder::new(&sk)
                    .ciphertext_level(0)
                    .key_level(0)
                    .enable_inner_sum()
                    .enable_expansion(params.degree().ilog2() as usize)
                    .build(&mut rng)?;
                let bytes = ek.to_bytes();
                assert_eq!(ek, EvaluationKey::from_bytes(&bytes, &params)?);
            }
        }
        Ok(())
    }
}
