//! Create parameters for the BFV encryption scheme

use crate::bfv::{context::CipherPlainContext, context::ContextLevel};
use crate::proto::bfv::{Parameters, parameters::PlaintextModulus as PlaintextModulusProto};
use crate::{Error, ParametersError, Result, SerializationError};
use fhe_math::{
    ntt::NttOperator,
    rns::{RnsContext, ScalingFactor},
    rq::{Context, Poly, Representation, scaler::Scaler, traits::TryConvertFrom},
    zq::{Modulus, primes::generate_prime},
};
use fhe_traits::{Deserialize, FheParameters, Serialize};
use itertools::Itertools;
use num_bigint::BigUint;
use num_traits::{PrimInt as _, ToPrimitive};
use prost::Message;
use std::collections::HashMap;
use std::fmt::Debug;
use std::sync::Arc;

/// Enum to support both small (u64) and large (BigUint) plaintext moduli.
#[derive(Debug, PartialEq, Eq, Clone)]
pub(crate) enum PlaintextModulus {
    Small {
        modulus: Modulus,
        modulus_big: BigUint,
    },
    Large(BigUint),
}

impl PlaintextModulus {
    pub fn as_biguint(&self) -> &BigUint {
        match self {
            Self::Small { modulus_big, .. } => modulus_big,
            Self::Large(m) => m,
        }
    }

    pub fn as_u64(&self) -> Option<u64> {
        match self {
            Self::Small { modulus, .. } => Some(**modulus),
            Self::Large(_) => None,
        }
    }

    pub fn reduce_vec(&self, v: &mut [BigUint]) {
        match self {
            Self::Small { modulus_big, .. } => {
                v.iter_mut().for_each(|vi| *vi %= modulus_big);
            }
            Self::Large(m) => v.iter_mut().for_each(|vi| *vi %= m),
        }
    }

    // Helper to reduce BigUint vector to i64 (centered), returning as Vec<BigUint>
    // or similar? The previous implementation used center_vec_vt returning
    // Vec<i64>. If modulus is large, we can't fit in i64.

    // We need a scalar multiplication for Plaintext::to_poly
    pub fn scalar_mul_vec(&self, a: &mut [BigUint], b: &BigUint) {
        match self {
            Self::Small { modulus_big, .. } => {
                a.iter_mut()
                    .for_each(|ai| *ai = (ai as &BigUint * b) % modulus_big);
            }
            Self::Large(m) => a.iter_mut().for_each(|ai| *ai = (ai as &BigUint * b) % m),
        }
    }
}

/// Parameters for the BFV encryption scheme.
///
/// This struct consolidates all parameter-specific data and pre-computed values
/// needed for BFV operations. It contains the raw parameters as well as
/// operational contexts and pre-computed scaling factors.
#[derive(PartialEq, Eq)]
pub struct BfvParameters {
    /// Number of coefficients in a polynomial.
    polynomial_degree: usize,

    /// Vector of coprime moduli q_i for the ciphertext.
    pub(crate) moduli: Box<[u64]>,

    /// Vector of the sized of the coprime moduli q_i for the ciphertext.
    moduli_sizes: Box<[usize]>,

    /// Error variance
    pub(crate) variance: usize,

    /// Head of the context chain for modulus switching
    pub(crate) context_chain: Arc<ContextLevel>,

    /// NTT operator for SIMD plaintext operations, if possible
    pub(crate) ntt_operator: Option<Arc<NttOperator>>,

    /// Plaintext Modulus as a Modulus type or BigUint
    pub(crate) plaintext: PlaintextModulus,

    pub(crate) matrix_reps_index_map: Box<[usize]>,
}

impl Debug for BfvParameters {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BfvParameters")
            .field("polynomial_degree", &self.polynomial_degree)
            .field("plaintext_modulus", &self.plaintext.as_biguint())
            .field("moduli", &self.moduli)
            .finish()
    }
}

impl FheParameters for BfvParameters {}

unsafe impl Send for BfvParameters {}

impl BfvParameters {
    /// Returns the underlying polynomial degree
    #[must_use]
    pub const fn degree(&self) -> usize {
        self.polynomial_degree
    }

    /// Returns a reference to the ciphertext moduli
    #[must_use]
    pub fn moduli(&self) -> &[u64] {
        &self.moduli
    }

    /// Returns a reference to the ciphertext moduli
    #[must_use]
    pub fn moduli_sizes(&self) -> &[usize] {
        &self.moduli_sizes
    }

    /// Returns the plaintext modulus if it fits in u64.
    /// Panics if the modulus is too large.
    #[must_use]
    pub fn plaintext(&self) -> u64 {
        self.plaintext.as_u64().unwrap()
    }

    /// Returns the plaintext modulus as BigUint
    #[must_use]
    pub fn plaintext_big(&self) -> &BigUint {
        self.plaintext.as_biguint()
    }

    /// Returns the maximum level allowed by these parameters.
    #[must_use]
    pub fn max_level(&self) -> usize {
        self.moduli.len() - 1
    }

    /// Returns the context corresponding to the level.
    pub fn context_at_level(&self, level: usize) -> Result<&Arc<Context>> {
        let mut current: &ContextLevel = &self.context_chain;
        while current.level < level {
            current = current
                .next
                .get()
                .ok_or_else(|| Error::InvalidLevel {
                    level,
                    min_level: 0,
                    max_level: self.max_level(),
                })?
                .as_ref();
        }
        if current.level == level {
            Ok(&current.poly_context)
        } else {
            Err(Error::InvalidLevel {
                level,
                min_level: 0,
                max_level: self.max_level(),
            })
        }
    }

    /// Returns the level of a given context
    pub fn level_of_context(&self, ctx: &Arc<Context>) -> Result<usize> {
        self.context_chain
            .poly_context
            .niterations_to(ctx)
            .map_err(Error::MathError)
    }

    /// Return head of context chain
    #[must_use]
    pub fn context_chain(&self) -> Arc<ContextLevel> {
        self.context_chain.clone()
    }

    /// Get context level at a specific depth
    pub fn context_level_at(&self, level: usize) -> Result<Arc<ContextLevel>> {
        let mut current = self.context_chain.clone();
        while current.level < level {
            match current.next.get() {
                Some(n) => current = n.clone(),
                None => {
                    return Err(Error::InvalidLevel {
                        level,
                        min_level: 0,
                        max_level: self.max_level(),
                    });
                }
            }
        }
        Ok(current)
    }

    /// Iterator over default parameters providing about 128 bits of security
    /// according to the <https://homomorphicencryption.org> standard.
    /// Filters out parameters where the modulus product bitlength is smaller
    /// than the plaintext modulus bitlength.
    ///
    /// Returns an error if no parameters are available after filtering.
    pub fn default_parameters_128(
        plaintext_nbits: usize,
    ) -> Result<impl Iterator<Item = Arc<BfvParameters>>> {
        debug_assert!(plaintext_nbits < 64);

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

        let parameters: Vec<Arc<BfvParameters>> = n_and_qs
            .into_iter()
            .sorted_by_key(|(n, _)| *n)
            .filter_map(move |(n, moduli)| {
                generate_prime(
                    plaintext_nbits,
                    2 * n as u64,
                    u64::MAX >> (64 - plaintext_nbits),
                )
                .and_then(|plaintext_modulus| {
                    // Calculate the bitlength of the product of moduli
                    let modulus_product_bitlength = moduli
                        .iter()
                        .map(|&m| 64 - m.leading_zeros() as usize)
                        .sum::<usize>();

                    // Filter out parameters where modulus product bitlength < plaintext bitlength
                    if modulus_product_bitlength >= plaintext_nbits {
                        BfvParametersBuilder::new()
                            .set_degree(n as usize)
                            .set_plaintext_modulus(plaintext_modulus)
                            .set_moduli(&moduli)
                            .build_arc()
                            .ok()
                    } else {
                        None
                    }
                })
            })
            .collect();

        // Check if we have any valid parameters after filtering
        if parameters.is_empty() {
            return Err(Error::ParametersError(
                ParametersError::NoParametersAvailable {
                    reason: format!(
                        "No default parameters available for plaintext modulus of {plaintext_nbits} bits. All parameter sets have modulus product bitlength smaller than the plaintext modulus."
                    ),
                },
            ));
        }

        Ok(parameters.into_iter())
    }

    #[cfg(test)]
    /// Returns default parameters for tests.
    #[must_use]
    #[expect(clippy::panic, reason = "panic indicates violated internal invariant")]
    pub fn default_arc(num_moduli: usize, degree: usize) -> Arc<Self> {
        if !degree.is_power_of_two() || degree < 8 {
            panic!("Invalid degree");
        }
        BfvParametersBuilder::new()
            .set_degree(degree)
            .set_plaintext_modulus(1153)
            .set_moduli_sizes(&vec![62usize; num_moduli])
            .build_arc()
            .unwrap()
    }
}

/// Builder for parameters for the Bfv encryption scheme.
#[derive(Debug)]
pub struct BfvParametersBuilder {
    degree: usize,
    plaintext: BigUint,
    variance: usize,
    ciphertext_moduli: Vec<u64>,
    ciphertext_moduli_sizes: Vec<usize>,
}

impl BfvParametersBuilder {
    /// Creates a new instance of the builder
    #[expect(
        clippy::new_without_default,
        reason = "builder requires explicit configuration"
    )]
    #[must_use]
    pub fn new() -> Self {
        Self {
            degree: Default::default(),
            plaintext: Default::default(),
            variance: 10,
            ciphertext_moduli: Default::default(),
            ciphertext_moduli_sizes: Default::default(),
        }
    }

    /// Sets the polynomial degree. Returns an error if the degree is not
    /// a power of two larger or equal to 8.
    pub fn set_degree(&mut self, degree: usize) -> &mut Self {
        self.degree = degree;
        self
    }

    /// Sets the plaintext modulus.
    pub fn set_plaintext_modulus(&mut self, plaintext: u64) -> &mut Self {
        self.plaintext = BigUint::from(plaintext);
        self
    }

    /// Sets the plaintext modulus as BigUint.
    pub fn set_plaintext_modulus_biguint(&mut self, plaintext: BigUint) -> &mut Self {
        self.plaintext = plaintext;
        self
    }

    /// Sets the sizes of the ciphertext moduli.
    /// Only one of `set_moduli_sizes` and `set_moduli`
    /// can be specified.
    pub fn set_moduli_sizes(&mut self, sizes: &[usize]) -> &mut Self {
        sizes.clone_into(&mut self.ciphertext_moduli_sizes);
        self
    }

    /// Sets the ciphertext moduli to use.
    /// Only one of `set_moduli_sizes` and `set_moduli`
    /// can be specified.
    pub fn set_moduli(&mut self, moduli: &[u64]) -> &mut Self {
        moduli.clone_into(&mut self.ciphertext_moduli);
        self
    }

    /// Sets the error variance. Returns an error if the variance is not between
    /// one and sixteen.
    pub fn set_variance(&mut self, variance: usize) -> &mut Self {
        self.variance = variance;
        self
    }

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

    /// Build a new `BfvParameters` inside an `Arc`.
    pub fn build_arc(&self) -> Result<Arc<BfvParameters>> {
        self.build().map(Arc::new)
    }

    /// Build a new `BfvParameters`.
    pub fn build(&self) -> Result<BfvParameters> {
        // Check that the degree is a power of 2 (and large enough).
        if self.degree < 8 || !self.degree.is_power_of_two() {
            return Err(Error::ParametersError(
                ParametersError::invalid_degree_with_bounds(self.degree),
            ));
        }

        let plaintext_modulus_struct = if let Some(p) = self.plaintext.to_u64() {
            PlaintextModulus::Small {
                modulus: Modulus::new(p).map_err(|e| {
                    Error::ParametersError(ParametersError::InvalidPlaintextModulus {
                        modulus: p,
                        reason: e.to_string(),
                    })
                })?,
                modulus_big: BigUint::from(p),
            }
        } else {
            PlaintextModulus::Large(self.plaintext.clone())
        };
        let plaintext_big = plaintext_modulus_struct.as_biguint();

        // Check that one of `ciphertext_moduli` and `ciphertext_moduli_sizes` is
        // specified.
        if !self.ciphertext_moduli.is_empty() && !self.ciphertext_moduli_sizes.is_empty() {
            return Err(Error::ParametersError(ParametersError::ConflictingParameters {
                conflict: "Only one of `ciphertext_moduli` and `ciphertext_moduli_sizes` can be specified".into(),
            }));
        } else if self.ciphertext_moduli.is_empty() && self.ciphertext_moduli_sizes.is_empty() {
            return Err(Error::ParametersError(ParametersError::MissingParameter {
                parameter: "ciphertext_moduli or ciphertext_moduli_sizes".into(),
            }));
        }

        // Get or generate the moduli
        let mut moduli = self.ciphertext_moduli.clone();
        if !self.ciphertext_moduli_sizes.is_empty() {
            moduli = Self::generate_moduli(&self.ciphertext_moduli_sizes, self.degree)?
        }

        // Recomputes the moduli sizes
        let moduli_sizes = moduli
            .iter()
            .map(|m| 64 - m.leading_zeros() as usize)
            .collect_vec();

        // Determine how many moduli needed for plaintext context
        // We need product of moduli > plaintext modulus.
        let t_bits = plaintext_big.bits();
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

        // Create NTT operator for SIMD operations if possible
        // Only if plaintext modulus fits in u64 for now
        let ntt_operator = match &plaintext_modulus_struct {
            PlaintextModulus::Small { modulus, .. } => {
                NttOperator::new(modulus, self.degree).map(Arc::new)
            }
            PlaintextModulus::Large(_) => None,
        };

        // Create cipher-plain bridge contexts
        let mut cipher_plain_contexts = Vec::with_capacity(moduli.len());

        // Build contexts in reverse order to establish the chain
        for i in (0..moduli.len()).rev() {
            let level_moduli = &moduli[..moduli.len() - i];
            let cipher_ctx = Context::new_arc(level_moduli, self.degree)?;
            // Compute delta (scaling polynomial)
            let mut delta_rests = vec![];
            for m in level_moduli {
                let q = Modulus::new(*m)?;
                let t_mod_q = (plaintext_big % *m).to_u64().unwrap();
                let neg_t_mod_q = q.neg(t_mod_q);
                if let Some(inv) = q.inv(neg_t_mod_q) {
                    delta_rests.push(inv);
                } else {
                    Err(Error::MathError(fhe_math::Error::Default(
                        "Inverse failed".to_string(),
                    )))?;
                }
            }

            // Use RnsContext to lift the delta values and create the scaling polynomial
            let rns = RnsContext::new(level_moduli)?;
            let mut delta = Poly::try_convert_from(
                &[rns.lift((&delta_rests).into())],
                &cipher_ctx,
                true,
                Representation::PowerBasis,
            )?;
            delta.change_representation(Representation::NttShoup);

            // Compute q_mod_t
            let q_mod_t = rns.modulus() % plaintext_big;

            // Compute plain_threshold
            let plain_threshold = match &plaintext_modulus_struct {
                PlaintextModulus::Small { modulus, .. } => BigUint::from((**modulus + 1) >> 1),
                PlaintextModulus::Large(m) => (m + 1u32) >> 1,
            };

            // Scaler from ciphertext to plaintext context
            let scaler = Scaler::new(
                &cipher_ctx,
                &plaintext_context,
                ScalingFactor::new(plaintext_big, rns.modulus()),
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

        // Build linked context chain
        let nodes: Vec<Arc<ContextLevel>> = cipher_plain_contexts
            .iter()
            .enumerate()
            .map(|(lvl, cp_ctx)| {
                Arc::new(ContextLevel::new(
                    cp_ctx.ciphertext_context.clone(),
                    cp_ctx.clone(),
                    lvl,
                ))
            })
            .collect();
        for i in 0..nodes.len() - 1 {
            let (prev, rest) = nodes.split_at(i + 1);
            ContextLevel::chain(&prev[i], &rest[0]);
        }
        let context_chain = nodes.first().unwrap().clone();

        // Create n+1 moduli of 62 bits for multiplication.
        let mut extended_basis = Vec::with_capacity(moduli.len() + 1);
        let mut upper_bound = 1 << 62;
        while extended_basis.len() != moduli.len() + 1 {
            upper_bound = generate_prime(62, 2 * self.degree as u64, upper_bound).unwrap();
            if !extended_basis.contains(&upper_bound) && !moduli.contains(&upper_bound) {
                extended_basis.push(upper_bound)
            }
        }

        // Compute multiplication parameters for each level
        for (i, node) in nodes.iter().enumerate() {
            // For the first multiplication, we want to extend to a context that
            // is ~60 bits larger.
            let modulus_size = moduli_sizes[..moduli_sizes.len() - i].iter().sum::<usize>();
            let n_moduli = (modulus_size + 60).div_ceil(62);
            let mut mul_1_moduli = vec![];
            mul_1_moduli.append(&mut moduli[..moduli_sizes.len() - i].to_vec());
            mul_1_moduli.append(&mut extended_basis[..n_moduli].to_vec());
            let mul_1_ctx = Context::new_arc(&mul_1_moduli, self.degree)?;
            let mp = MultiplicationParameters::new(
                &node.poly_context,
                &mul_1_ctx,
                ScalingFactor::one(),
                ScalingFactor::new(plaintext_big, node.poly_context.modulus()),
            )?;
            node.mul_params.set(mp).unwrap();
        }

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

        Ok(BfvParameters {
            polynomial_degree: self.degree,
            moduli: moduli.into(),
            moduli_sizes: moduli_sizes.into(),
            variance: self.variance,
            context_chain,
            ntt_operator,
            plaintext: plaintext_modulus_struct,
            matrix_reps_index_map: matrix_reps_index_map.into(),
        })
    }
}

impl Serialize for BfvParameters {
    fn to_bytes(&self) -> Vec<u8> {
        let plaintext_modulus = if let Some(plaintext_u64) = self.plaintext.as_u64() {
            Some(PlaintextModulusProto::Plaintext(plaintext_u64))
        } else {
            Some(PlaintextModulusProto::PlaintextBig(
                self.plaintext.as_biguint().to_bytes_le(),
            ))
        };

        Parameters {
            degree: self.polynomial_degree as u32,
            moduli: self.moduli.to_vec(),
            variance: self.variance as u32,
            plaintext_modulus,
        }
        .encode_to_vec()
    }
}

impl Deserialize for BfvParameters {
    fn try_deserialize(bytes: &[u8]) -> Result<Self> {
        let params: Parameters = Message::decode(bytes).map_err(|_| {
            Error::SerializationError(SerializationError::ProtobufError {
                message: "Parameters decode".into(),
            })
        })?;

        let plaintext_modulus = match params.plaintext_modulus {
            Some(PlaintextModulusProto::Plaintext(value)) => BigUint::from(value),
            Some(PlaintextModulusProto::PlaintextBig(bytes)) => BigUint::from_bytes_le(&bytes),
            None => {
                return Err(Error::SerializationError(
                    SerializationError::MissingField {
                        field_name: "Parameters.plaintext_modulus".into(),
                    },
                ));
            }
        };

        BfvParametersBuilder::new()
            .set_degree(params.degree as usize)
            .set_plaintext_modulus_biguint(plaintext_modulus)
            .set_moduli(&params.moduli)
            .set_variance(params.variance as usize)
            .build()
    }
    type Error = Error;
}

/// Multiplication parameters
#[derive(Debug, Clone, PartialEq, Eq, Default)]
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
    use super::{BfvParameters, BfvParametersBuilder};
    use crate::proto::bfv::{Parameters, parameters::PlaintextModulus as PlaintextModulusProto};
    use fhe_traits::{Deserialize, Serialize};
    use num_bigint::BigUint;
    use prost::Message;
    use std::error::Error;

    #[test]
    fn default() {
        let params = BfvParameters::default_arc(1, 16);
        assert_eq!(params.moduli.len(), 1);
        assert_eq!(params.degree(), 16);

        let params = BfvParameters::default_arc(2, 16);
        assert_eq!(params.moduli.len(), 2);
        assert_eq!(params.degree(), 16);
    }

    #[test]
    fn ciphertext_moduli() -> Result<(), Box<dyn Error>> {
        let params = BfvParametersBuilder::new()
            .set_degree(16)
            .set_plaintext_modulus(2)
            .set_moduli_sizes(&[62, 62, 62, 61, 60, 11])
            .build()?;
        assert_eq!(
            params.moduli.to_vec(),
            &[
                4611686018427387617,
                4611686018427387329,
                4611686018427387073,
                2305843009213693921,
                1152921504606845473,
                2017
            ]
        );

        let params = BfvParametersBuilder::new()
            .set_degree(16)
            .set_plaintext_modulus(2)
            .set_moduli(&[
                4611686018427387617,
                4611686018427387329,
                4611686018427387073,
                2305843009213693921,
                1152921504606845473,
                2017,
            ])
            .build()?;
        assert_eq!(params.moduli_sizes.to_vec(), &[62, 62, 62, 61, 60, 11]);

        Ok(())
    }

    #[test]
    fn big_plaintext_modulus() -> Result<(), Box<dyn Error>> {
        // Use a 128-bit prime
        let p = BigUint::parse_bytes(b"340282366920938463463374607431768211507", 10).unwrap();
        let params = BfvParametersBuilder::new()
            .set_degree(16)
            .set_plaintext_modulus_biguint(p.clone())
            .set_moduli_sizes(&[62, 62, 62, 62, 62]) // Large enough for product > p
            .build()?;

        assert_eq!(params.plaintext_big(), &p);
        Ok(())
    }

    #[test]
    fn serialize() -> Result<(), Box<dyn Error>> {
        let params = BfvParametersBuilder::new()
            .set_degree(16)
            .set_plaintext_modulus(2)
            .set_moduli_sizes(&[62, 62, 62, 61, 60, 11])
            .set_variance(4)
            .build()?;
        let bytes = params.to_bytes();
        let proto = Parameters::decode(bytes.as_slice())?;
        assert!(matches!(
            proto.plaintext_modulus,
            Some(PlaintextModulusProto::Plaintext(2))
        ));
        assert_eq!(BfvParameters::try_deserialize(&bytes)?, params);

        // Test with big plaintext
        let p = BigUint::parse_bytes(b"340282366920938463463374607431768211507", 10).unwrap();
        let params = BfvParametersBuilder::new()
            .set_degree(16)
            .set_plaintext_modulus_biguint(p)
            .set_moduli_sizes(&[62, 62, 62, 62, 62])
            .set_variance(4)
            .build()?;
        let bytes = params.to_bytes();
        let proto = Parameters::decode(bytes.as_slice())?;
        let proto_plaintext_bytes = match &proto.plaintext_modulus {
            Some(PlaintextModulusProto::PlaintextBig(bytes)) => bytes.as_slice(),
            _ => return Err("expected plaintext_big variant".into()),
        };
        assert_eq!(
            proto_plaintext_bytes,
            params.plaintext_big().to_bytes_le().as_slice()
        );
        let decoded = BfvParameters::try_deserialize(&bytes)?;
        assert_eq!(decoded, params);
        assert_eq!(decoded.plaintext_big(), params.plaintext_big());

        Ok(())
    }

    #[test]
    fn deserialize_missing_plaintext_modulus() {
        let proto = Parameters {
            degree: 16,
            moduli: vec![4611686018427387617, 4611686018427387329],
            variance: 4,
            plaintext_modulus: None,
        };
        let bytes = proto.encode_to_vec();
        let err = BfvParameters::try_deserialize(&bytes).unwrap_err();
        assert!(format!("{err}").contains("Missing required field"));
    }

    #[test]
    fn matrix_reps_index_map_is_permutation() -> Result<(), Box<dyn Error>> {
        let params = BfvParametersBuilder::new()
            .set_degree(16)
            .set_plaintext_modulus(2)
            .set_moduli_sizes(&[62, 62])
            .build()?;

        let mut map = params.matrix_reps_index_map.to_vec();
        assert_eq!(map.len(), params.degree());

        map.sort_unstable();
        map.dedup();
        assert_eq!(map.len(), params.degree());

        Ok(())
    }

    #[test]
    fn default_parameters_iterator() {
        let mut it = BfvParameters::default_parameters_128(20).unwrap();
        assert!(it.next().is_some());
    }

    #[test]
    fn default_parameters_filtering() {
        // Test that parameters are filtered correctly
        let params: Vec<_> = BfvParameters::default_parameters_128(20).unwrap().collect();

        // All returned parameters should have sufficient modulus bitlength
        for param in &params {
            let modulus_product_bitlength = param.moduli_sizes.iter().sum::<usize>();
            assert!(modulus_product_bitlength >= 20);
        }

        // Test with a very small plaintext modulus for which we won't be able to
        // create any parameters
        let result = BfvParameters::default_parameters_128(10);
        assert!(result.is_err());

        #[expect(clippy::panic, reason = "panic indicates violated internal invariant")]
        match result {
            Err(e) => {
                let error_string = format!("{e}");
                assert!(error_string.contains("No parameters available"));
                assert!(error_string.contains("10 bits"));
            }
            Ok(_) => panic!("Expected error"),
        }
    }
}
