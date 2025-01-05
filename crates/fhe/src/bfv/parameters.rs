//! Create parameters for the BFV encryption scheme

use crate::proto::bfv::Parameters;
use crate::{Error, ParametersError, Result};
use fhe_math::{
    ntt::NttOperator,
    rns::{RnsContext, ScalingFactor},
    rq::{scaler::Scaler, traits::TryConvertFrom, Context, Poly, Representation},
    zq::{primes::generate_prime, Modulus},
};
use fhe_traits::{Deserialize, FheParameters, Serialize};
use itertools::Itertools;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use prost::Message;
use std::collections::HashMap;
use std::fmt::Debug;
use std::sync::Arc;

/// Parameters for the BFV encryption scheme.
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(PartialEq, Eq)]
pub struct BfvParameters {
    /// Number of coefficients in a polynomial.
    polynomial_degree: usize,

    /// Modulus of the plaintext.
    plaintext_modulus: u64,

    /// Vector of coprime moduli q_i for the ciphertext.
    /// One and only one of `ciphertext_moduli` or `ciphertext_moduli_sizes`
    /// must be specified.
    pub(crate) moduli: Box<[u64]>,

    /// Vector of the sized of the coprime moduli q_i for the ciphertext.
    /// One and only one of `ciphertext_moduli` or `ciphertext_moduli_sizes`
    /// must be specified.
    moduli_sizes: Box<[usize]>,

    /// Error variance
    pub(crate) variance: usize,

    /// Context for the underlying polynomials
    pub(crate) ctx: Vec<Arc<Context>>,

    /// Ntt operator for the SIMD plaintext, if possible.
    pub(crate) op: Option<Arc<NttOperator>>,

    /// Scaling polynomial for the plaintext
    pub(crate) delta: Box<[Poly]>,

    /// Q modulo the plaintext modulus
    pub(crate) q_mod_t: Box<[u64]>,

    /// Down scaler for the plaintext
    pub(crate) scalers: Box<[Scaler]>,

    /// Plaintext Modulus
    pub(crate) plaintext: Modulus,

    // Parameters for the multiplications
    pub(crate) mul_params: Box<[MultiplicationParameters]>,

    pub(crate) matrix_reps_index_map: Box<[usize]>,
}

impl Debug for BfvParameters {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BfvParameters")
            .field("polynomial_degree", &self.polynomial_degree)
            .field("plaintext_modulus", &self.plaintext_modulus)
            .field("moduli", &self.moduli)
            // .field("moduli_sizes", &self.moduli_sizes)
            // .field("variance", &self.variance)
            // .field("ctx", &self.ctx)
            // .field("op", &self.op)
            // .field("delta", &self.delta)
            // .field("q_mod_t", &self.q_mod_t)
            // .field("scaler", &self.scaler)
            // .field("plaintext", &self.plaintext)
            // .field("mul_params", &self.mul_params)
            // .field("matrix_reps_index_map", &self.matrix_reps_index_map)
            .finish()
    }
}

impl FheParameters for BfvParameters {}

unsafe impl Send for BfvParameters {}

impl BfvParameters {
    /// Returns the underlying polynomial degree
    pub const fn degree(&self) -> usize {
        self.polynomial_degree
    }

    /// Returns a reference to the ciphertext moduli
    pub fn moduli(&self) -> &[u64] {
        &self.moduli
    }

    /// Returns a reference to the ciphertext moduli
    pub fn moduli_sizes(&self) -> &[usize] {
        &self.moduli_sizes
    }

    /// Returns the plaintext modulus
    pub const fn plaintext(&self) -> u64 {
        self.plaintext_modulus
    }

    /// Returns the maximum level allowed by these parameters.
    pub fn max_level(&self) -> usize {
        self.moduli.len() - 1
    }

    /// Returns the context corresponding to the level.
    pub(crate) fn ctx_at_level(&self, level: usize) -> Result<&Arc<Context>> {
        self.ctx
            .get(level)
            .ok_or_else(|| Error::DefaultError("No context".to_string()))
    }

    /// Returns the level of a given context
    pub(crate) fn level_of_ctx(&self, ctx: &Arc<Context>) -> Result<usize> {
        self.ctx[0].niterations_to(ctx).map_err(Error::MathError)
    }

    /// Vector of default parameters providing about 128 bits of security
    /// according to the <https://homomorphicencryption.org> standard.
    pub fn default_parameters_128(plaintext_nbits: usize) -> Vec<Arc<BfvParameters>> {
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

        let mut params = vec![];

        for n in n_and_qs.keys().sorted() {
            let moduli = n_and_qs.get(n).unwrap();
            if let Some(plaintext_modulus) = generate_prime(
                plaintext_nbits,
                2 * *n as u64,
                u64::MAX >> (64 - plaintext_nbits),
            ) {
                params.push(
                    BfvParametersBuilder::new()
                        .set_degree(*n as usize)
                        .set_plaintext_modulus(plaintext_modulus)
                        .set_moduli(moduli)
                        .build_arc()
                        .unwrap(),
                )
            }
        }

        params
    }

    #[cfg(test)]
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
    plaintext: u64,
    variance: usize,
    ciphertext_moduli: Vec<u64>,
    ciphertext_moduli_sizes: Vec<usize>,
}

impl BfvParametersBuilder {
    /// Creates a new instance of the builder
    #[allow(clippy::new_without_default)]
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

    /// Sets the plaintext modulus. Returns an error if the plaintext is not
    /// between 2 and 2^62 - 1.
    pub fn set_plaintext_modulus(&mut self, plaintext: u64) -> &mut Self {
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
        for size in moduli_sizes {
            if *size > 62 || *size < 10 {
                return Err(Error::ParametersError(ParametersError::InvalidModulusSize(
                    *size, 10, 62,
                )));
            }

            let mut upper_bound = 1 << size;
            loop {
                if let Some(prime) = generate_prime(*size, 2 * degree as u64, upper_bound) {
                    if !moduli.contains(&prime) {
                        moduli.push(prime);
                        break;
                    } else {
                        upper_bound = prime;
                    }
                } else {
                    return Err(Error::ParametersError(ParametersError::NotEnoughPrimes(
                        *size, degree,
                    )));
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
            return Err(Error::ParametersError(ParametersError::InvalidDegree(
                self.degree,
            )));
        }

        // This checks that the plaintext modulus is valid.
        // TODO: Check bound on the plaintext modulus.
        let plaintext_modulus = Modulus::new(self.plaintext).map_err(|e| {
            Error::ParametersError(ParametersError::InvalidPlaintext(e.to_string()))
        })?;

        // Check that one of `ciphertext_moduli` and `ciphertext_moduli_sizes` is
        // specified.
        if !self.ciphertext_moduli.is_empty() && !self.ciphertext_moduli_sizes.is_empty() {
            return Err(Error::ParametersError(ParametersError::TooManySpecified(
                "Only one of `ciphertext_moduli` and `ciphertext_moduli_sizes` can be specified"
                    .to_string(),
            )));
        } else if self.ciphertext_moduli.is_empty() && self.ciphertext_moduli_sizes.is_empty() {
            return Err(Error::ParametersError(ParametersError::TooFewSpecified(
                "One of `ciphertext_moduli` and `ciphertext_moduli_sizes` must be specified"
                    .to_string(),
            )));
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

        // Create n+1 moduli of 62 bits for multiplication.
        let mut extended_basis = Vec::with_capacity(moduli.len() + 1);
        let mut upper_bound = 1 << 62;
        while extended_basis.len() != moduli.len() + 1 {
            upper_bound = generate_prime(62, 2 * self.degree as u64, upper_bound).unwrap();
            if !extended_basis.contains(&upper_bound) && !moduli.contains(&upper_bound) {
                extended_basis.push(upper_bound)
            }
        }

        let op = NttOperator::new(&plaintext_modulus, self.degree);

        let plaintext_ctx = Context::new_arc(&moduli[..1], self.degree)?;

        let mut delta_rests = vec![];
        for m in &moduli {
            let q = Modulus::new(*m)?;
            delta_rests.push(q.inv(q.neg(*plaintext_modulus)).unwrap())
        }

        let mut ctx = Vec::with_capacity(moduli.len());
        let mut delta = Vec::with_capacity(moduli.len());
        let mut q_mod_t = Vec::with_capacity(moduli.len());
        let mut scalers = Vec::with_capacity(moduli.len());
        let mut mul_params = Vec::with_capacity(moduli.len());
        for i in 0..moduli.len() {
            let rns = RnsContext::new(&moduli[..moduli.len() - i])?;
            let ctx_i = Context::new_arc(&moduli[..moduli.len() - i], self.degree)?;
            let mut p = Poly::try_convert_from(
                &[rns.lift((&delta_rests).into())],
                &ctx_i,
                true,
                Representation::PowerBasis,
            )?;
            p.change_representation(Representation::NttShoup);
            delta.push(p);

            q_mod_t.push((rns.modulus() % *plaintext_modulus).to_u64().unwrap());

            scalers.push(Scaler::new(
                &ctx_i,
                &plaintext_ctx,
                ScalingFactor::new(&BigUint::from(*plaintext_modulus), rns.modulus()),
            )?);

            // For the first multiplication, we want to extend to a context that
            // is ~60 bits larger.
            let modulus_size = moduli_sizes[..moduli_sizes.len() - i].iter().sum::<usize>();
            let n_moduli = (modulus_size + 60).div_ceil(62);
            let mut mul_1_moduli = vec![];
            mul_1_moduli.append(&mut moduli[..moduli_sizes.len() - i].to_vec());
            mul_1_moduli.append(&mut extended_basis[..n_moduli].to_vec());
            let mul_1_ctx = Context::new_arc(&mul_1_moduli, self.degree)?;
            mul_params.push(MultiplicationParameters::new(
                &ctx_i,
                &mul_1_ctx,
                ScalingFactor::one(),
                ScalingFactor::new(&BigUint::from(*plaintext_modulus), ctx_i.modulus()),
            )?);

            ctx.push(ctx_i);
        }

        // We use the same code as SEAL
        // https://github.com/microsoft/SEAL/blob/82b07db635132e297282649e2ab5908999089ad2/native/src/seal/batchencoder.cpp
        let row_size = self.degree >> 1;
        let m = self.degree << 1;
        let gen = 3;
        let mut pos = 1;
        let mut matrix_reps_index_map = vec![0usize; self.degree];
        for i in 0..row_size {
            let index1 = (pos - 1) >> 1;
            let index2 = (m - pos - 1) >> 1;
            matrix_reps_index_map[i] = index1.reverse_bits() >> (self.degree.leading_zeros() + 1);
            matrix_reps_index_map[row_size | i] =
                index2.reverse_bits() >> (self.degree.leading_zeros() + 1);
            pos *= gen;
            pos &= m - 1;
        }

        Ok(BfvParameters {
            polynomial_degree: self.degree,
            plaintext_modulus: self.plaintext,
            moduli: moduli.into(),
            moduli_sizes: moduli_sizes.into(),
            variance: self.variance,
            ctx,
            op: op.map(Arc::new),
            delta: delta.into(),
            q_mod_t: q_mod_t.into(),
            scalers: scalers.into(),
            plaintext: plaintext_modulus,
            mul_params: mul_params.into(),
            matrix_reps_index_map: matrix_reps_index_map.into(),
        })
    }
}

impl Serialize for BfvParameters {
    fn to_bytes(&self) -> Vec<u8> {
        Parameters {
            degree: self.polynomial_degree as u32,
            plaintext: self.plaintext_modulus,
            moduli: self.moduli.to_vec(),
            variance: self.variance as u32,
        }
        .encode_to_vec()
    }
}

impl Deserialize for BfvParameters {
    fn try_deserialize(bytes: &[u8]) -> Result<Self> {
        let params: Parameters = Message::decode(bytes).map_err(|_| Error::SerializationError)?;
        BfvParametersBuilder::new()
            .set_degree(params.degree as usize)
            .set_plaintext_modulus(params.plaintext)
            .set_moduli(&params.moduli)
            .set_variance(params.variance as usize)
            .build()
    }
    type Error = Error;
}

/// Multiplication parameters
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, PartialEq, Eq, Default)]
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
    use fhe_traits::{Deserialize, Serialize};
    use std::error::Error;

    // TODO: To fix when errors handling is fixed.
    // #[test]
    // fn builder()  -> Result<(), Box<dyn Error>> {
    // 	let params = BfvParametersBuilder::new().build();
    // 	assert!(params.is_err_and(|e| e.to_string() == "Unspecified degree"));

    // 	assert!(BfvParametersBuilder::new()
    // 		.set_degree(7)
    // 		.build()
    // 		.is_err_and(
    // 			|e| e.to_string() == "The degree should be a power of two larger or equal to
    // 8" 		));

    // 	assert!(BfvParametersBuilder::new()
    // 		.set_degree(1023)
    // 		.build()
    // 		.is_err_and(
    // 			|e| e.to_string() == "The degree should be a power of two larger or equal to
    // 8" 		));

    // 	let params = BfvParametersBuilder::new().set_degree(1024).build();
    // 	assert!(params.is_err_and(|e| e.to_string() == "Unspecified plaintext
    // modulus"));

    // 	assert!(BfvParametersBuilder::new()
    // 		.set_degree(1024)
    // 		.set_plaintext_modulus(0)
    // 		.build()
    // 		.is_err_and(|e| e.to_string() == "modulus should be between 2 and 2^62-1"));

    // 	let params = BfvParametersBuilder::new()
    // 		.set_degree(1024)
    // 		.set_plaintext_modulus(2)
    // 		.build();
    // 	assert!(params.is_err_and(|e| e.to_string() == "Unspecified ciphertext
    // moduli"));

    // 	assert!(BfvParametersBuilder::new()
    // 		.set_degree(1024)
    // 		.set_plaintext_modulus(2)
    // 		.set_moduli(&[])
    // 		.build()
    // 		.is_err_and(|e| e.to_string() == "Unspecified ciphertext moduli"));

    // 	assert!(BfvParametersBuilder::new()
    // 		.set_degree(1024)
    // 		.set_plaintext_modulus(2)
    // 		.set_moduli(&[1153])
    // 		.set_moduli_sizes(&[62])
    // 		.build()
    // 		.is_err_and(|e| e.to_string() == "The set of ciphertext moduli is already
    // specified"));

    // 	assert!(BfvParametersBuilder::new()
    // 		.set_degree(8)
    // 		.set_plaintext_modulus(2)
    // 		.set_moduli(&[1])
    // 		.build()
    // 		.is_err_and(|e| e.to_string() == "modulus should be between 2 and 2^62-1"));

    // 	let params = BfvParametersBuilder::new()
    // 		.set_degree(8)
    // 		.set_plaintext_modulus(2)
    // 		.set_moduli(&[2])
    // 		.build();
    // 	assert!(params.is_err_and(|e| e.to_string() == "Impossible to construct a Ntt
    // operator"));

    // 	let params = BfvParametersBuilder::new()
    // 		.set_degree(8)
    // 		.set_plaintext_modulus(2)
    // 		.set_moduli(&[1153])
    // 		.build();
    // 	assert!(params.is_ok());

    // 	let params = params.unwrap();
    // 	assert_eq!(params.ciphertext_moduli, vec![1153]);
    // 	assert_eq!(params.moduli(), vec![1153]);
    // 	assert_eq!(params.plaintext_modulus, 2);
    // 	assert_eq!(params.polynomial_degree, 16);
    // 	assert_eq!(params.degree(), 16);
    // 	assert_eq!(params.variance, 1);
    // 	assert!(params.op.is_none());

    // 	Ok(())
    // }

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
    fn serialize() -> Result<(), Box<dyn Error>> {
        let params = BfvParametersBuilder::new()
            .set_degree(16)
            .set_plaintext_modulus(2)
            .set_moduli_sizes(&[62, 62, 62, 61, 60, 11])
            .set_variance(4)
            .build()?;
        let bytes = params.to_bytes();
        assert_eq!(BfvParameters::try_deserialize(&bytes)?, params);
        Ok(())
    }
}
