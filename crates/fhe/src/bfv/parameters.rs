//! Create parameters for the BFV encryption scheme

use crate::bfv::{context::CipherPlainContext, context::ContextLevel};
use crate::proto::bfv::{
    Parameters as ParametersProto, parameters::PlaintextModulus as PlaintextModulusProto,
};
use crate::{Error, Result, error::ParametersError, error::SerializationError};
use fhe_math::{
    ntt::NttOperator,
    rns::{RnsContext, ScalingFactor},
    rq::{Context, Poly, PowerBasis, scaler::Scaler},
    zq::{Modulus, primes::generate_prime},
};

use fhe_util::is_prime;
use itertools::Itertools;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use prost::Message;
use std::collections::HashMap;
use std::fmt::Debug;
use std::sync::Arc;

/// A plaintext modulus with an optional machine-word fast path.
///
/// The `BigUint` value is canonical. `small` caches the equivalent `Modulus`
/// when it fits, allowing performance-sensitive encoding and decryption to
/// keep using specialized `u64` arithmetic without exposing two independent
/// representations to the rest of the BFV implementation.
#[derive(Debug, PartialEq, Eq, Clone)]
pub(crate) struct PlaintextModulus {
    value: BigUint,
    small: Option<Modulus>,
}

impl PlaintextModulus {
    fn try_new(value: BigUint) -> Result<Self> {
        let small = value
            .to_u64()
            .map(|modulus| {
                Modulus::new(modulus).map_err(|source| {
                    Error::ParametersError(ParametersError::InvalidPlaintextModulus {
                        modulus,
                        source,
                    })
                })
            })
            .transpose()?;
        Ok(Self { value, small })
    }

    pub(crate) fn as_biguint(&self) -> &BigUint {
        &self.value
    }

    pub(crate) fn as_u64(&self) -> Option<u64> {
        self.small.as_ref().map(|modulus| **modulus)
    }

    pub(crate) fn small(&self) -> Option<&Modulus> {
        self.small.as_ref()
    }

    pub(crate) fn reduce_vec(&self, v: &mut [BigUint]) {
        v.iter_mut().for_each(|vi| *vi %= &self.value);
    }

    pub(crate) fn scalar_mul_vec(&self, a: &mut [BigUint], b: &BigUint) {
        a.iter_mut()
            .for_each(|ai| *ai = (ai as &BigUint * b) % &self.value);
    }

    fn ntt_operator(&self, degree: usize) -> Option<Arc<NttOperator>> {
        self.small
            .as_ref()
            .and_then(|modulus| NttOperator::new(modulus, degree).map(Arc::new))
    }

    fn upper_half_threshold(&self) -> BigUint {
        (&self.value + 1u32) >> 1
    }
}

/// Parameters for the BFV encryption scheme.
///
/// Cloning this immutable handle shares its precomputed contexts. Independently
/// built handles are compatible when their defining settings are equal.
#[derive(Clone)]
pub struct Parameters {
    pub(crate) inner: Arc<ParametersInner>,
}

pub(crate) struct ParametersInner {
    /// Number of coefficients in a polynomial.
    polynomial_degree: usize,

    /// Vector of coprime moduli q_i for the ciphertext.
    pub(crate) moduli: Box<[u64]>,

    /// Vector of the sized of the coprime moduli q_i for the ciphertext.
    moduli_sizes: Box<[usize]>,

    /// Error variance
    pub(crate) variance: usize,

    /// Precomputed contexts indexed by modulus-switching level.
    pub(crate) context_levels: Vec<ContextLevel>,

    /// NTT operator for SIMD plaintext operations, if possible
    pub(crate) ntt_operator: Option<Arc<NttOperator>>,

    /// Plaintext Modulus as a Modulus type or BigUint
    pub(crate) plaintext: PlaintextModulus,

    pub(crate) matrix_reps_index_map: Box<[usize]>,
}

impl Debug for Parameters {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Parameters")
            .field("polynomial_degree", &self.inner.polynomial_degree)
            .field("plaintext_modulus", &self.inner.plaintext.as_biguint())
            .field("moduli", &self.inner.moduli)
            .field("noise_variance", &self.inner.variance)
            .finish()
    }
}

impl PartialEq for Parameters {
    fn eq(&self, other: &Self) -> bool {
        self.compatible(other)
    }
}
impl Eq for Parameters {}

impl Parameters {
    /// Start configuring an immutable parameter handle.
    #[must_use]
    pub fn builder() -> ParametersBuilder {
        ParametersBuilder::default()
    }

    /// Compare the defining arithmetic and noise settings. Cloned handles share
    /// precomputation; independently built equal settings are also compatible.
    /// This does not establish that ciphertexts use the same secret key.
    #[must_use]
    pub fn compatible(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.inner, &other.inner)
            || (self.inner.polynomial_degree == other.inner.polynomial_degree
                && self.inner.moduli == other.inner.moduli
                && self.plaintext_modulus() == other.plaintext_modulus()
                && self.inner.variance == other.inner.variance)
    }

    /// Variance of the centered binomial error distribution.
    #[must_use]
    pub fn noise_variance(&self) -> usize {
        self.inner.variance
    }

    /// Returns the underlying polynomial degree
    #[must_use]
    pub fn degree(&self) -> usize {
        self.inner.polynomial_degree
    }

    /// Returns a reference to the ciphertext moduli
    #[must_use]
    pub fn moduli(&self) -> &[u64] {
        &self.inner.moduli
    }

    /// Returns a reference to the ciphertext moduli
    #[must_use]
    pub fn moduli_sizes(&self) -> &[usize] {
        &self.inner.moduli_sizes
    }

    /// Returns the plaintext modulus when it fits in a machine word.
    #[must_use]
    pub fn plaintext_modulus_u64(&self) -> Option<u64> {
        self.inner.plaintext.as_u64()
    }

    /// Returns the plaintext modulus as BigUint
    #[must_use]
    pub fn plaintext_modulus(&self) -> &BigUint {
        self.inner.plaintext.as_biguint()
    }

    /// Returns the maximum level allowed by these parameters.
    #[must_use]
    pub fn max_level(&self) -> usize {
        self.inner.context_levels.len() - 1
    }

    /// Returns the context corresponding to the level.
    pub fn context_at_level(&self, level: usize) -> Result<&Arc<Context>> {
        self.inner
            .context_levels
            .get(level)
            .map(|context_level| &context_level.poly_context)
            .ok_or_else(|| Error::InvalidLevel {
                level,
                min_level: 0,
                max_level: self.max_level(),
            })
    }

    /// Returns the level of a given context
    pub fn level_of_context(&self, ctx: &Arc<Context>) -> Result<usize> {
        let level = self
            .inner
            .moduli
            .len()
            .checked_sub(ctx.moduli().len())
            .ok_or(Error::MathError(fhe_math::Error::ContextNotReachable))?;
        let context_level = self
            .inner
            .context_levels
            .get(level)
            .ok_or(Error::MathError(fhe_math::Error::ContextNotReachable))?;
        if Arc::ptr_eq(&context_level.poly_context, ctx) || &context_level.poly_context == ctx {
            Ok(level)
        } else {
            Err(Error::MathError(fhe_math::Error::ContextNotReachable))
        }
    }

    /// Return all contexts in modulus-switching order.
    #[must_use]
    pub fn context_levels(&self) -> &[ContextLevel] {
        &self.inner.context_levels
    }

    /// Get the precomputed data for a specific modulus-switching level.
    pub fn context_level_at(&self, level: usize) -> Result<&ContextLevel> {
        self.inner
            .context_levels
            .get(level)
            .ok_or_else(|| Error::InvalidLevel {
                level,
                min_level: 0,
                max_level: self.max_level(),
            })
    }

    /// List lightweight profiles providing about 128 bits of security
    /// according to the <https://homomorphicencryption.org> standard.
    /// Filters out profiles without an appropriate plaintext prime or with an
    /// incompatible plaintext modulus. This does not build NTT or RNS contexts.
    ///
    /// Returns an error if no parameters are available after filtering.
    pub fn profiles_128(plaintext_nbits: usize) -> Result<impl Iterator<Item = ParameterProfile>> {
        if !(2..64).contains(&plaintext_nbits) {
            return Err(ParametersError::NoDefaultParameters {
                plaintext_bits: plaintext_nbits,
            }
            .into());
        }

        let mut n_and_qs = HashMap::new();
        n_and_qs.insert(1024, vec![0x7e00001]);
        n_and_qs.insert(2048, vec![0x3fffffff000001]);
        n_and_qs.insert(4096, vec![0xffffee001, 0xffffc4001, 0x1ffffe0001]);
        n_and_qs.insert(
            8192,
            vec![
                0x7fffffd8001,
                0x7fffffc8001,
                0xfffffffc001,
                0xffffff6c001,
                0xfffffebc001,
            ],
        );
        n_and_qs.insert(
            16384,
            vec![
                0xfffffffd8001,
                0xfffffffa0001,
                0xfffffff00001,
                0x1fffffff68001,
                0x1fffffff50001,
                0x1ffffffee8001,
                0x1ffffffea0001,
                0x1ffffffe88001,
                0x1ffffffe48001,
            ],
        );

        let parameters: Vec<ParameterProfile> = n_and_qs
            .into_iter()
            .sorted_by_key(|(n, _)| *n)
            .filter_map(move |(n, moduli)| {
                generate_prime(
                    plaintext_nbits,
                    2 * n as u64,
                    u64::MAX >> (64 - plaintext_nbits),
                )
                .and_then(|plaintext_modulus| {
                    // Listing must not promise a profile whose build would fail.
                    let product: BigUint = moduli.iter().map(|q| BigUint::from(*q)).product();
                    if BigUint::from(plaintext_modulus) < product
                        && moduli.iter().all(|q| plaintext_modulus % q != 0)
                    {
                        Some(ParameterProfile {
                            degree: n as usize,
                            plaintext_modulus,
                            moduli,
                        })
                    } else {
                        None
                    }
                })
            })
            .collect();

        // Check if we have any valid parameters after filtering
        if parameters.is_empty() {
            return Err(Error::ParametersError(
                ParametersError::NoDefaultParameters {
                    plaintext_bits: plaintext_nbits,
                },
            ));
        }

        Ok(parameters.into_iter())
    }

    /// Build only the preselected 128-bit profile for this degree and plaintext
    /// bit length. These retain the project's existing profile assumptions;
    /// arbitrary builder configurations are not security estimates.
    pub fn profile_128(degree: usize, plaintext_bits: usize) -> Result<Self> {
        Self::profiles_128(plaintext_bits)?
            .find(|profile| profile.degree() == degree)
            .ok_or(ParametersError::UnavailableProfile {
                degree,
                plaintext_bits,
            })?
            .build()
    }

    #[cfg(test)]
    /// Returns default parameters for tests.
    #[must_use]
    #[expect(clippy::panic, reason = "panic indicates violated internal invariant")]
    pub fn test_parameters(num_moduli: usize, degree: usize) -> Self {
        if !degree.is_power_of_two() || degree < 8 {
            panic!("Invalid degree");
        }
        ParametersBuilder::new()
            .degree(degree)
            .plaintext_modulus(1153_u64)
            .ciphertext_modulus_bits(vec![62usize; num_moduli])
            .build()
            .unwrap()
    }
}

/// A lightweight description of one of the existing preselected parameter sets.
#[derive(Clone, Debug)]
pub struct ParameterProfile {
    degree: usize,
    plaintext_modulus: u64,
    moduli: Vec<u64>,
}

impl ParameterProfile {
    /// Polynomial degree of this profile.
    #[must_use]
    pub fn degree(&self) -> usize {
        self.degree
    }
    /// Plaintext prime selected for the requested bit length.
    #[must_use]
    pub fn plaintext_modulus(&self) -> u64 {
        self.plaintext_modulus
    }
    /// Ordered ciphertext primes of this profile.
    #[must_use]
    pub fn ciphertext_moduli(&self) -> &[u64] {
        &self.moduli
    }
    /// Compute this profile's shared contexts.
    pub fn build(self) -> Result<Parameters> {
        Parameters::builder()
            .degree(self.degree)
            .plaintext_modulus(self.plaintext_modulus)
            .ciphertext_moduli(self.moduli)
            .build()
    }
}

/// Consuming builder for BFV parameters. Required fields remain unspecified
/// until set; the most recent ciphertext modulus setter replaces the previous
/// one.
#[derive(Clone, Debug)]
pub struct ParametersBuilder {
    degree: Option<usize>,
    plaintext: Option<BigUint>,
    variance: usize,
    moduli: Option<ModuliSpec>,
}

#[derive(Clone, Debug)]
enum ModuliSpec {
    Explicit(Vec<u64>),
    Bits(Vec<usize>),
}

impl Default for ParametersBuilder {
    fn default() -> Self {
        Self {
            degree: None,
            plaintext: None,
            variance: 10,
            moduli: None,
        }
    }
}

impl ParametersBuilder {
    /// Create an empty builder, with noise variance 10.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }
    /// Set a power-of-two polynomial degree between 8 and 65536.
    #[must_use]
    pub fn degree(mut self, degree: usize) -> Self {
        self.degree = Some(degree);
        self
    }
    /// Set the plaintext modulus, at least two and coprime to the ciphertext
    /// primes.
    #[must_use]
    pub fn plaintext_modulus(mut self, plaintext: impl Into<BigUint>) -> Self {
        self.plaintext = Some(plaintext.into());
        self
    }
    /// Select explicit ciphertext primes, replacing any earlier modulus
    /// selection.
    #[must_use]
    pub fn ciphertext_moduli(mut self, moduli: impl AsRef<[u64]>) -> Self {
        self.moduli = Some(ModuliSpec::Explicit(moduli.as_ref().to_vec()));
        self
    }
    /// Select ciphertext prime bit lengths, replacing any earlier modulus
    /// selection.
    #[must_use]
    pub fn ciphertext_modulus_bits(mut self, bits: impl AsRef<[usize]>) -> Self {
        self.moduli = Some(ModuliSpec::Bits(bits.as_ref().to_vec()));
        self
    }
    /// Set the centered binomial variance, between one and thirty-two.
    #[must_use]
    pub fn noise_variance(mut self, variance: usize) -> Self {
        self.variance = variance;
        self
    }
    /// Validate the settings and compute a cheap-to-clone immutable handle.
    /// Validation establishes arithmetic invariants, not a security estimate.
    pub fn build(self) -> Result<Parameters> {
        ParameterConfiguration {
            degree: self.degree.ok_or(ParametersError::MissingDegree)?,
            plaintext: self
                .plaintext
                .ok_or(ParametersError::MissingPlaintextModulus)?,
            variance: self.variance,
            moduli: self
                .moduli
                .ok_or(ParametersError::MissingCiphertextModulusSpecification)?,
        }
        .build()
    }
}

struct ParameterConfiguration {
    degree: usize,
    plaintext: BigUint,
    variance: usize,
    moduli: ModuliSpec,
}

impl ParameterConfiguration {
    const MIN_DEGREE: usize = 8;
    const MAX_DEGREE: usize = 65536;
    const MIN_VARIANCE: usize = 1;
    const MAX_VARIANCE: usize = 32;
    /// Generate ciphertext moduli with the specified sizes
    fn generate_moduli(moduli_sizes: &[usize], degree: usize) -> Result<Vec<u64>> {
        let mut moduli = vec![];
        let required_counts = moduli_sizes.iter().copied().counts();
        let mut generated_counts: HashMap<usize, usize> = HashMap::new();
        for (i, size) in moduli_sizes.iter().enumerate() {
            if *size > 62 || *size < 10 {
                return Err(Error::ParametersError(
                    ParametersError::InvalidModulusSize {
                        index: i,
                        size: *size,
                        min: 10,
                        max: 62,
                    },
                ));
            }

            let mut upper_bound = 1 << size;
            loop {
                if let Some(prime) = generate_prime(*size, 2 * degree as u64, upper_bound) {
                    if !moduli.contains(&prime) {
                        moduli.push(prime);
                        *generated_counts.entry(*size).or_insert(0) += 1;
                        break;
                    } else {
                        upper_bound = prime;
                    }
                } else {
                    let needed = *required_counts.get(size).unwrap_or(&0);
                    let available = *generated_counts.get(size).unwrap_or(&0);
                    return Err(Error::ParametersError(ParametersError::NotEnoughPrimes {
                        size: *size,
                        degree,
                        needed,
                        available,
                    }));
                }
            }
        }

        Ok(moduli)
    }

    fn gcd(mut a: u64, mut b: u64) -> u64 {
        while b != 0 {
            (a, b) = (b, a % b);
        }
        a
    }

    fn validate_configuration(&self) -> Result<()> {
        if !(Self::MIN_DEGREE..=Self::MAX_DEGREE).contains(&self.degree)
            || !self.degree.is_power_of_two()
        {
            return Err(Error::ParametersError(
                ParametersError::invalid_degree_with_bounds(self.degree),
            ));
        }

        if !(Self::MIN_VARIANCE..=Self::MAX_VARIANCE).contains(&self.variance) {
            return Err(Error::ParametersError(ParametersError::InvalidVariance {
                variance: self.variance,
                min: Self::MIN_VARIANCE,
                max: Self::MAX_VARIANCE,
            }));
        }

        Ok(())
    }

    fn validate_moduli(&self, moduli: &[u64], plaintext: &BigUint) -> Result<()> {
        for (index, modulus) in moduli.iter().copied().enumerate() {
            Modulus::new(modulus).map_err(|error| {
                Error::ParametersError(ParametersError::InvalidCiphertextModulus {
                    index,
                    modulus,
                    source: error,
                })
            })?;

            let indices = moduli
                .iter()
                .enumerate()
                .filter_map(|(i, candidate)| (*candidate == modulus).then_some(i))
                .collect_vec();
            if indices.len() > 1 {
                return Err(Error::ParametersError(ParametersError::DuplicateModuli {
                    modulus,
                    indices,
                }));
            }
        }

        for (i, modulus1) in moduli.iter().copied().enumerate() {
            for modulus2 in moduli.iter().copied().skip(i + 1) {
                let gcd = Self::gcd(modulus1, modulus2);
                if gcd != 1 {
                    return Err(Error::ParametersError(ParametersError::ModuliNotCoprime {
                        modulus1,
                        modulus2,
                        gcd,
                    }));
                }
            }
        }

        for (index, modulus) in moduli.iter().copied().enumerate() {
            if modulus % (2 * self.degree as u64) != 1 || !is_prime(modulus) {
                return Err(Error::ParametersError(
                    ParametersError::CiphertextModulusNotNttFriendly {
                        index,
                        modulus,
                        degree: self.degree,
                    },
                ));
            }
        }

        let ciphertext_modulus = moduli
            .iter()
            .map(|m| BigUint::from(*m))
            .product::<BigUint>();
        if plaintext >= &ciphertext_modulus {
            return Err(Error::ParametersError(
                ParametersError::PlaintextModulusExceedsCiphertextModulus {
                    plaintext_modulus: plaintext.clone(),
                    ciphertext_modulus,
                },
            ));
        }

        for (index, modulus) in moduli.iter().copied().enumerate() {
            let plaintext_mod_modulus = (plaintext % modulus).to_u64().ok_or({
                Error::ParametersError(ParametersError::PlaintextReductionFailed {
                    ciphertext_modulus: modulus,
                })
            })?;
            let gcd = Self::gcd(plaintext_mod_modulus, modulus);
            if gcd != 1 {
                return Err(Error::ParametersError(
                    ParametersError::PlaintextModulusNotCoprime {
                        plaintext_modulus: plaintext.clone(),
                        ciphertext_modulus: modulus,
                        index,
                        gcd,
                    },
                ));
            }
        }

        Ok(())
    }

    fn build(self) -> Result<Parameters> {
        self.validate_configuration()?;

        let plaintext_modulus_struct = PlaintextModulus::try_new(self.plaintext.clone())?;
        let plaintext_modulus = plaintext_modulus_struct.as_biguint();

        let moduli = match &self.moduli {
            ModuliSpec::Explicit(moduli) => moduli.clone(),
            ModuliSpec::Bits(bits) => Self::generate_moduli(bits, self.degree)?,
        };
        if moduli.is_empty() {
            return Err(ParametersError::MissingCiphertextModulusSpecification.into());
        }
        self.validate_moduli(&moduli, plaintext_modulus)?;

        // Recomputes the moduli sizes
        let moduli_sizes = moduli
            .iter()
            .map(|m| 64 - m.leading_zeros() as usize)
            .collect_vec();

        // Determine how many moduli needed for plaintext context
        // We need product of moduli > plaintext modulus.
        let t_bits = plaintext_modulus.bits();
        let mut accumulated_bits = 0;
        let mut plaintext_moduli_count = 0;
        for size in &moduli_sizes {
            accumulated_bits += size;
            plaintext_moduli_count += 1;
            if accumulated_bits as u64 >= t_bits + 60 {
                break;
            }
        }
        plaintext_moduli_count = std::cmp::max(plaintext_moduli_count, 1);
        plaintext_moduli_count = std::cmp::min(plaintext_moduli_count, moduli.len());

        // Create plaintext context using sufficient moduli
        let plaintext_context = Context::new_arc(&moduli[..plaintext_moduli_count], self.degree)?;

        // SIMD currently uses the cached machine-word representation.
        let ntt_operator = plaintext_modulus_struct.ntt_operator(self.degree);

        // Create cipher-plain bridge contexts
        let mut cipher_plain_contexts = Vec::with_capacity(moduli.len());

        // Build ciphertext/plaintext bridges for every modulus-switching level.
        for i in (0..moduli.len()).rev() {
            let level_moduli = &moduli[..moduli.len() - i];
            let cipher_ctx = Context::new_arc(level_moduli, self.degree)?;
            // Compute delta (scaling polynomial)
            let mut delta_rests = vec![];
            for m in level_moduli {
                let q = Modulus::new(*m)?;
                let t_mod_q = (plaintext_modulus % *m).to_u64().unwrap();
                let neg_t_mod_q = q.neg(t_mod_q);
                if let Some(inv) = q.inv(neg_t_mod_q) {
                    delta_rests.push(inv);
                } else {
                    return Err(Error::MathError(fhe_math::Error::NonInvertible {
                        value: neg_t_mod_q,
                        modulus: *m,
                    }));
                }
            }

            // Use RnsContext to lift the delta values and create the scaling polynomial
            let rns = RnsContext::new(level_moduli)?;
            let delta = Poly::<PowerBasis>::from_biguint_coefficients_with_timing(
                &[rns.lift((&delta_rests).into())],
                &cipher_ctx,
                Some(crate::VariableTime::new(crate::PublicData::assert_public())),
            )?
            .into_ntt_shoup();

            // Compute q_mod_t
            let q_mod_t = rns.modulus() % plaintext_modulus;

            // Compute plain_threshold
            let plain_threshold = plaintext_modulus_struct.upper_half_threshold();

            // Scaler from ciphertext to plaintext context
            let scaler = Scaler::new(
                &cipher_ctx,
                &plaintext_context,
                ScalingFactor::new(plaintext_modulus, rns.modulus()),
            )?;

            let cipher_plain_ctx = CipherPlainContext::new_arc(
                &plaintext_context,
                &cipher_ctx,
                delta,
                q_mod_t,
                plain_threshold,
                scaler,
            );

            cipher_plain_contexts.push(cipher_plain_ctx.clone());
        }

        // Reverse to get correct order (level 0 first)
        cipher_plain_contexts.reverse();

        // Create n+1 moduli of 62 bits for multiplication.
        let mut extended_basis = Vec::with_capacity(moduli.len() + 1);
        let mut upper_bound = 1 << 62;
        while extended_basis.len() != moduli.len() + 1 {
            upper_bound =
                generate_prime(62, 2 * self.degree as u64, upper_bound).ok_or_else(|| {
                    Error::ParametersError(ParametersError::NotEnoughPrimes {
                        size: 62,
                        degree: self.degree,
                        needed: moduli.len() + 1,
                        available: extended_basis.len(),
                    })
                })?;
            if !extended_basis.contains(&upper_bound) && !moduli.contains(&upper_bound) {
                extended_basis.push(upper_bound)
            }
        }

        // Build a fully initialized context for each level. The vector index is
        // the level, so lookups do not need to walk or initialize a chain.
        let context_levels = cipher_plain_contexts
            .into_iter()
            .enumerate()
            .map(|(level, cipher_plain_context)| {
                let poly_context = cipher_plain_context.ciphertext_context.clone();

                // For the first multiplication, extend to a context that is
                // approximately 60 bits larger.
                let modulus_size = moduli_sizes[..moduli_sizes.len() - level]
                    .iter()
                    .sum::<usize>();
                let n_moduli = (modulus_size + 60).div_ceil(62);
                let mut multiplication_moduli = moduli[..moduli_sizes.len() - level].to_vec();
                multiplication_moduli.extend_from_slice(&extended_basis[..n_moduli]);
                let multiplication_context = Context::new_arc(&multiplication_moduli, self.degree)?;
                let mul_params = MultiplicationParameters::new(
                    &poly_context,
                    &multiplication_context,
                    ScalingFactor::one(),
                    ScalingFactor::new(plaintext_modulus, poly_context.modulus()),
                )?;

                Ok(ContextLevel::new(
                    poly_context,
                    cipher_plain_context,
                    level,
                    mul_params,
                ))
            })
            .collect::<Result<Vec<_>>>()?;

        // We use the same code as SEAL
        // https://github.com/microsoft/SEAL/blob/82b07db635132e297282649e2ab5908999089ad2/native/src/seal/batchencoder.cpp
        let row_size = self.degree >> 1;
        let m = self.degree << 1;
        let generator = 3;
        let mut pos = 1;
        let mut matrix_reps_index_map = vec![0usize; self.degree];
        for i in 0..row_size {
            let index1 = (pos - 1) >> 1;
            let index2 = (m - pos - 1) >> 1;
            matrix_reps_index_map[i] = index1.reverse_bits() >> (self.degree.leading_zeros() + 1);
            matrix_reps_index_map[row_size | i] =
                index2.reverse_bits() >> (self.degree.leading_zeros() + 1);
            pos *= generator;
            pos &= m - 1;
        }

        Ok(Parameters {
            inner: Arc::new(ParametersInner {
                polynomial_degree: self.degree,
                moduli: moduli.into(),
                moduli_sizes: moduli_sizes.into(),
                variance: self.variance,
                context_levels,
                ntt_operator,
                plaintext: plaintext_modulus_struct,
                matrix_reps_index_map: matrix_reps_index_map.into(),
            }),
        })
    }
}

impl Parameters {
    /// Serialize in the existing protobuf wire format.
    #[must_use]
    pub fn to_bytes(&self) -> Vec<u8> {
        let plaintext_modulus = if let Some(plaintext_u64) = self.inner.plaintext.as_u64() {
            Some(PlaintextModulusProto::Plaintext(plaintext_u64))
        } else {
            Some(PlaintextModulusProto::PlaintextBig(
                self.inner.plaintext.as_biguint().to_bytes_le(),
            ))
        };

        ParametersProto {
            degree: self.inner.polynomial_degree as u32,
            moduli: self.inner.moduli.to_vec(),
            variance: self.inner.variance as u32,
            plaintext_modulus,
        }
        .encode_to_vec()
    }
}

impl Parameters {
    /// Import validated protobuf bytes, binding contextual values to the
    /// supplied parameters.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self> {
        Self::from_bytes_with_limits(bytes, &crate::DecodeLimits::default())
    }

    /// Import with explicit resource bounds checked before allocation.
    pub fn from_bytes_with_limits(bytes: &[u8], limits: &crate::DecodeLimits) -> Result<Self> {
        crate::bfv::wire::preflight(
            bytes,
            crate::error::SerializedObject::Parameters,
            None,
            limits,
        )?;
        let params: ParametersProto = Message::decode(bytes).map_err(|source| {
            Error::SerializationError(SerializationError::Decode {
                object: crate::error::SerializedObject::Parameters,
                source,
            })
        })?;

        let plaintext_modulus = match params.plaintext_modulus {
            Some(PlaintextModulusProto::Plaintext(value)) => BigUint::from(value),
            Some(PlaintextModulusProto::PlaintextBig(bytes)) => BigUint::from_bytes_le(&bytes),
            None => {
                return Err(Error::SerializationError(
                    SerializationError::MissingField {
                        field: crate::error::SerializedField::ParametersPlaintextModulus,
                    },
                ));
            }
        };

        ParametersBuilder::new()
            .degree(params.degree as usize)
            .plaintext_modulus(plaintext_modulus)
            .ciphertext_moduli(&params.moduli)
            .noise_variance(params.variance as usize)
            .build()
    }
}

/// Multiplication parameters
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct MultiplicationParameters {
    pub(crate) extender: Scaler,
    pub(crate) down_scaler: Scaler,
    pub(crate) from: Arc<Context>,
    pub(crate) to: Arc<Context>,
}

impl MultiplicationParameters {
    fn new(
        from: &Arc<Context>,
        to: &Arc<Context>,
        up_self_factor: ScalingFactor,
        down_factor: ScalingFactor,
    ) -> Result<Self> {
        Ok(Self {
            extender: Scaler::new(from, to, up_self_factor)?,
            down_scaler: Scaler::new(to, from, down_factor)?,
            from: from.clone(),
            to: to.clone(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::{Parameters, ParametersBuilder};
    use crate::proto::bfv::{
        Parameters as ParametersProto, parameters::PlaintextModulus as PlaintextModulusProto,
    };
    use crate::{Error as FheError, error::ParametersError};

    use num_bigint::BigUint;
    use prost::Message;
    use std::error::Error;

    #[test]
    fn default() {
        let params = Parameters::test_parameters(1, 16);
        assert_eq!(params.inner.moduli.len(), 1);
        assert_eq!(params.degree(), 16);
        assert!(params.inner.plaintext.small().is_some());
        assert_eq!(
            params.inner.plaintext.as_u64(),
            Some(params.plaintext_modulus_u64().unwrap())
        );

        let params = Parameters::test_parameters(2, 16);
        assert_eq!(params.inner.moduli.len(), 2);
        assert_eq!(params.degree(), 16);
    }

    #[test]
    fn ciphertext_moduli() -> Result<(), Box<dyn Error>> {
        let params = ParametersBuilder::new()
            .degree(16)
            .plaintext_modulus(2_u64)
            .ciphertext_modulus_bits([62, 62, 62, 61, 60, 11])
            .build()?;
        assert_eq!(
            params.inner.moduli.to_vec(),
            &[
                4611686018427387617,
                4611686018427387329,
                4611686018427387073,
                2305843009213693921,
                1152921504606845473,
                2017
            ]
        );

        let params = ParametersBuilder::new()
            .degree(16)
            .plaintext_modulus(2_u64)
            .ciphertext_moduli([
                4611686018427387617,
                4611686018427387329,
                4611686018427387073,
                2305843009213693921,
                1152921504606845473,
                2017,
            ])
            .build()?;
        assert_eq!(
            params.inner.moduli_sizes.to_vec(),
            &[62, 62, 62, 61, 60, 11]
        );

        Ok(())
    }

    #[test]
    fn big_plaintext_modulus() -> Result<(), Box<dyn Error>> {
        // Use a 128-bit prime
        let p = BigUint::parse_bytes(b"340282366920938463463374607431768211507", 10).unwrap();
        let params = ParametersBuilder::new()
            .degree(16)
            .plaintext_modulus(p.clone())
            .ciphertext_modulus_bits([62, 62, 62, 62, 62]) // Large enough for product > p
            .build()?;

        assert_eq!(params.plaintext_modulus(), &p);
        assert!(params.inner.plaintext.small().is_none());
        assert_eq!(params.inner.plaintext.as_u64(), None);
        Ok(())
    }

    #[test]
    fn serialize() -> Result<(), Box<dyn Error>> {
        let params = ParametersBuilder::new()
            .degree(16)
            .plaintext_modulus(2_u64)
            .ciphertext_modulus_bits([62, 62, 62, 61, 60, 11])
            .noise_variance(4)
            .build()?;
        let bytes = params.to_bytes();
        let proto = ParametersProto::decode(bytes.as_slice())?;
        assert!(matches!(
            proto.plaintext_modulus,
            Some(PlaintextModulusProto::Plaintext(2))
        ));
        assert_eq!(Parameters::from_bytes(&bytes)?, params);

        // Test with big plaintext
        let p = BigUint::parse_bytes(b"340282366920938463463374607431768211507", 10).unwrap();
        let params = ParametersBuilder::new()
            .degree(16)
            .plaintext_modulus(p)
            .ciphertext_modulus_bits([62, 62, 62, 62, 62])
            .noise_variance(4)
            .build()?;
        let bytes = params.to_bytes();
        let proto = ParametersProto::decode(bytes.as_slice())?;
        let proto_plaintext_bytes = match &proto.plaintext_modulus {
            Some(PlaintextModulusProto::PlaintextBig(bytes)) => bytes.as_slice(),
            _ => return Err("expected plaintext_modulus variant".into()),
        };
        assert_eq!(
            proto_plaintext_bytes,
            params.plaintext_modulus().to_bytes_le().as_slice()
        );
        let decoded = Parameters::from_bytes(&bytes)?;
        assert_eq!(decoded, params);
        assert_eq!(decoded.plaintext_modulus(), params.plaintext_modulus());

        Ok(())
    }

    #[test]
    fn deserialize_missing_plaintext_modulus() {
        let proto = ParametersProto {
            degree: 16,
            moduli: vec![4611686018427387617, 4611686018427387329],
            variance: 4,
            plaintext_modulus: None,
        };
        let bytes = proto.encode_to_vec();
        let err = Parameters::from_bytes(&bytes).unwrap_err();
        assert_eq!(
            err,
            FheError::SerializationError(crate::error::SerializationError::MissingField {
                field: crate::error::SerializedField::ParametersPlaintextModulus,
            })
        );
    }

    #[test]
    fn rejects_invalid_degree() {
        for degree in [0, 10, 131072] {
            let err = ParametersBuilder::new()
                .degree(degree)
                .plaintext_modulus(2_u64)
                .ciphertext_moduli([97])
                .build()
                .unwrap_err();
            assert!(matches!(
                err,
                FheError::ParametersError(ParametersError::InvalidDegree {
                    degree: actual,
                    min: 8,
                    max: 65536,
                }) if actual == degree
            ));
        }
    }

    #[test]
    fn validates_variance_bounds() -> Result<(), Box<dyn Error>> {
        for variance in [0, 33] {
            let err = ParametersBuilder::new()
                .degree(16)
                .plaintext_modulus(2_u64)
                .ciphertext_moduli([97])
                .noise_variance(variance)
                .build()
                .unwrap_err();
            assert!(matches!(
                err,
                FheError::ParametersError(ParametersError::InvalidVariance {
                    variance: actual,
                    min: 1,
                    max: 32,
                }) if actual == variance
            ));
        }

        for variance in [1, 32] {
            ParametersBuilder::new()
                .degree(16)
                .plaintext_modulus(2_u64)
                .ciphertext_moduli([97])
                .noise_variance(variance)
                .build()?;
        }

        Ok(())
    }

    #[test]
    fn deserialization_rejects_invalid_variance() {
        let proto = ParametersProto {
            degree: 16,
            moduli: vec![97],
            variance: 33,
            plaintext_modulus: Some(PlaintextModulusProto::Plaintext(2)),
        };
        let err = Parameters::from_bytes(&proto.encode_to_vec()).unwrap_err();
        assert!(matches!(
            err,
            FheError::ParametersError(ParametersError::InvalidVariance {
                variance: 33,
                min: 1,
                max: 32,
            })
        ));
    }

    #[test]
    fn validates_explicit_ciphertext_moduli() {
        let invalid = ParametersBuilder::new()
            .degree(16)
            .plaintext_modulus(2_u64)
            .ciphertext_moduli([1])
            .build()
            .unwrap_err();
        assert!(matches!(
            invalid,
            FheError::ParametersError(ParametersError::InvalidCiphertextModulus {
                index: 0,
                modulus: 1,
                ..
            })
        ));

        let duplicate = ParametersBuilder::new()
            .degree(16)
            .plaintext_modulus(2_u64)
            .ciphertext_moduli([97, 97])
            .build()
            .unwrap_err();
        assert!(matches!(
            duplicate,
            FheError::ParametersError(ParametersError::DuplicateModuli {
                modulus: 97,
                indices,
            }) if indices == [0, 1]
        ));

        let not_coprime = ParametersBuilder::new()
            .degree(16)
            .plaintext_modulus(2_u64)
            .ciphertext_moduli([9, 15])
            .build()
            .unwrap_err();
        assert!(matches!(
            not_coprime,
            FheError::ParametersError(ParametersError::ModuliNotCoprime {
                modulus1: 9,
                modulus2: 15,
                gcd: 3,
            })
        ));

        let not_ntt_friendly = ParametersBuilder::new()
            .degree(16)
            .plaintext_modulus(2_u64)
            .ciphertext_moduli([17])
            .build()
            .unwrap_err();
        assert!(matches!(
            not_ntt_friendly,
            FheError::ParametersError(ParametersError::CiphertextModulusNotNttFriendly {
                index: 0,
                modulus: 17,
                degree: 16,
            })
        ));
    }

    #[test]
    fn validates_plaintext_against_ciphertext_moduli() {
        let too_large = ParametersBuilder::new()
            .degree(16)
            .plaintext_modulus(98_u64)
            .ciphertext_moduli([97])
            .build()
            .unwrap_err();
        assert!(matches!(
            too_large,
            FheError::ParametersError(
                ParametersError::PlaintextModulusExceedsCiphertextModulus { .. }
            )
        ));

        let not_coprime = ParametersBuilder::new()
            .degree(16)
            .plaintext_modulus(194_u64)
            .ciphertext_moduli([97, 193])
            .build()
            .unwrap_err();
        assert!(matches!(
            not_coprime,
            FheError::ParametersError(ParametersError::PlaintextModulusNotCoprime {
                ciphertext_modulus: 97,
                index: 0,
                gcd: 97,
                ..
            })
        ));
    }

    #[test]
    fn matrix_reps_index_map_is_permutation() -> Result<(), Box<dyn Error>> {
        let params = ParametersBuilder::new()
            .degree(16)
            .plaintext_modulus(2_u64)
            .ciphertext_modulus_bits([62, 62])
            .build()?;

        let mut map = params.inner.matrix_reps_index_map.to_vec();
        assert_eq!(map.len(), params.degree());

        map.sort_unstable();
        map.dedup();
        assert_eq!(map.len(), params.degree());

        Ok(())
    }

    #[test]
    fn default_parameters_iterator() {
        let mut it = Parameters::profiles_128(20).unwrap();
        assert!(it.next().is_some());
    }

    #[test]
    fn default_parameters_filtering() {
        // Test that parameters are filtered correctly
        let params: Vec<_> = Parameters::profiles_128(20).unwrap().collect();

        // All returned parameters should have sufficient modulus bitlength
        for param in &params {
            let modulus_product_bitlength = param
                .moduli
                .iter()
                .map(|q| (64 - q.leading_zeros()) as usize)
                .sum::<usize>();
            assert!(modulus_product_bitlength >= 20);
        }

        // Test with a very small plaintext modulus for which we won't be able to
        // create any parameters
        let result = Parameters::profiles_128(10);
        assert!(result.is_err());

        assert_eq!(
            result.err(),
            Some(FheError::ParametersError(
                ParametersError::NoDefaultParameters { plaintext_bits: 10 }
            ))
        );
    }
}
