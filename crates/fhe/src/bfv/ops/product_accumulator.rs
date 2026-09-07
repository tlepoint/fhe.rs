use fhe_math::rq::{Ntt, Poly};
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use zeroize::Zeroizing;

use crate::{
    CiphertextError, DotProductError, Result,
    bfv::{Ciphertext, Parameters},
};

/// Accumulate products of two-part BFV ciphertexts, then scale the sum once.
///
/// This computes a ciphertext dot product with fewer basis conversions than
/// multiplying each pair separately. Products are accumulated in the extended
/// multiplication basis. [`Self::finish`] applies the BFV scaling factor and
/// returns a three-part ciphertext, ready for optional relinearization.
///
/// Rounding once changes ciphertext coefficients and rounding noise relative to
/// a sum of separately rounded products; the results need not be
/// byte-identical. The usual BFV noise-budget requirement still applies. This
/// uses the same fixed-point scaler as ordinary multiplication, including its
/// precision limits.
///
/// The accumulator checks a conservative coefficient bound to prevent wrapping
/// in the extended basis. Its scratch is zeroized when finished or dropped.
pub struct CiphertextProductAccumulator {
    par: Parameters,
    level: usize,
    products: usize,
    max_products: usize,
    c: Zeroizing<Vec<Poly<Ntt>>>,
}

// A lifted input coefficient has magnitude at most Q (a conservative bound
// allowing for either centered representative at a rounding boundary). Each
// output coefficient sums at most 2*N products per ciphertext pair. Thus k
// pairs are unambiguous in the centered extended basis M if 4*k*N*Q^2 < M.
fn product_limit(q: &BigUint, m: &BigUint, degree: usize) -> usize {
    ((m - 1u32) / (q * q * degree * 4u32))
        .to_usize()
        .unwrap_or(usize::MAX)
}

impl CiphertextProductAccumulator {
    /// Create an empty accumulator at a fixed ciphertext level.
    pub fn new(par: &Parameters, level: usize) -> Result<Self> {
        let mp = par.context_level_at(level)?.mul_params();
        let max_products = product_limit(mp.from.modulus(), mp.to.modulus(), par.degree());
        let mut zero = Poly::zero(&mp.to);
        zero.allow_variable_time_computations(crate::VariableTime::new(
            crate::PublicData::assert_public(),
        ));
        Ok(Self {
            par: par.clone(),
            level,
            products: 0,
            max_products,
            c: Zeroizing::new(vec![zero; 3]),
        })
    }

    /// Maximum number of pairs permitted by the extended-basis coefficient
    /// bound.
    #[must_use]
    pub fn max_products(&self) -> usize {
        self.max_products
    }

    /// Add one product without BFV downscaling or relinearization.
    ///
    /// Both inputs must have two polynomial parts and match this accumulator's
    /// parameters and level. Invalid inputs or a full accumulator return an
    /// error before changing the accumulated sum. No input-dependent shortcut
    /// is used, and variable-time permission propagates from all input parts.
    pub fn add_product(&mut self, lhs: &Ciphertext, rhs: &Ciphertext) -> Result<()> {
        let mp = self.par.context_level_at(self.level)?.mul_params();
        lhs.validate_for_context(&self.par, self.level, &mp.from)?;
        rhs.validate_for_context(&self.par, self.level, &mp.from)?;
        if lhs.len() != 2 || rhs.len() != 2 {
            return Err(CiphertextError::MultiplicationPolynomialCount {
                left: lhs.len(),
                right: rhs.len(),
                expected: 2,
            }
            .into());
        }
        if self.products == self.max_products {
            return Err(DotProductError::TooManyProducts {
                maximum: self.max_products,
            }
            .into());
        }

        let mut left = Zeroizing::new(
            lhs.iter()
                .map(|p| p.scale(&mp.extender))
                .collect::<fhe_math::Result<Vec<_>>>()?,
        );
        let mut right = Zeroizing::new(
            rhs.iter()
                .map(|p| p.scale(&mp.extender))
                .collect::<fhe_math::Result<Vec<_>>>()?,
        );
        if !lhs
            .iter()
            .chain(rhs.iter())
            .all(Poly::allows_variable_time_computations)
        {
            for p in left.iter_mut().chain(right.iter_mut()) {
                p.disallow_variable_time_computations();
            }
        }
        let product = Zeroizing::new(super::tensor::product(&left, &right));
        for (sum, term) in self.c.iter_mut().zip(product.iter()) {
            *sum += term;
        }
        self.products += 1;
        Ok(())
    }

    /// Scale the sum and return a three-part ciphertext at the original level.
    /// Returns an error if no products were added. Consumes the accumulator.
    pub fn finish(self) -> Result<Ciphertext> {
        if self.products == 0 {
            return Err(DotProductError::EmptyInput.into());
        }
        let mp = self.par.context_level_at(self.level)?.mul_params();
        let c = self
            .c
            .iter()
            .map(|p| p.scale(&mp.down_scaler))
            .collect::<fhe_math::Result<Vec<_>>>()?;
        Ciphertext::from_components(c, &self.par)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bfv::{Encoding, ParametersBuilder, Plaintext, RelinearizationKey, SecretKey};
    use fhe_math::rq::traits::TryConvertFrom;

    use num_bigint::BigInt;
    use num_traits::{Signed, Zero};
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    #[test]
    fn exact_negacyclic_product_sum() -> Result<()> {
        let par = ParametersBuilder::new()
            .degree(16)
            .plaintext_modulus(17_u64)
            .ciphertext_modulus_bits([12, 12])
            .build()?;
        for level in 0..=par.max_level() {
            let ctx = par.context_at_level(level)?;
            let q = BigInt::from(ctx.modulus().clone());
            let half = (&q >> 1usize) - 2u32;
            for count in [1, 2, 17, 81] {
                let mut accumulator = CiphertextProductAccumulator::new(&par, level)?;
                let mut exact = vec![vec![BigInt::zero(); par.degree()]; 3];
                for pair in 0..count {
                    // Include near-maximal centered coefficients, both signs,
                    // cancellation and wraparound in X^N + 1.
                    let operands: Vec<Vec<Vec<BigInt>>> = (0..2)
                        .map(|side| {
                            (0..2)
                                .map(|part| {
                                    (0..par.degree())
                                        .map(|i| match (i + pair + 2 * side + part) % 5 {
                                            0 => half.clone(),
                                            1 => -&half,
                                            2 => BigInt::from(1),
                                            3 => BigInt::from(-1),
                                            _ => BigInt::zero(),
                                        })
                                        .collect()
                                })
                                .collect()
                        })
                        .collect();
                    let ciphertexts = operands
                        .iter()
                        .map(|parts| {
                            let c = parts
                                .iter()
                                .map(|coefficients| {
                                    let residues: Vec<_> = coefficients
                                        .iter()
                                        .map(|x| ((x + &q) % &q).to_biguint().unwrap())
                                        .collect();
                                    Poly::<Ntt>::try_convert_from(residues.as_slice(), ctx)
                                })
                                .collect::<fhe_math::Result<Vec<_>>>()?;
                            Ciphertext::from_components(c, &par)
                        })
                        .collect::<Result<Vec<_>>>()?;
                    accumulator.add_product(&ciphertexts[0], &ciphertexts[1])?;
                    for a in 0..2 {
                        for b in 0..2 {
                            for i in 0..par.degree() {
                                for j in 0..par.degree() {
                                    let term = &operands[0][a][i] * &operands[1][b][j];
                                    let coefficient = &mut exact[a + b][(i + j) % par.degree()];
                                    if i + j < par.degree() {
                                        *coefficient += term;
                                    } else {
                                        *coefficient -= term;
                                    }
                                }
                            }
                        }
                    }
                }
                let result = accumulator.finish()?;
                for (part, expected) in result.iter().zip(exact) {
                    let expected: Vec<_> = expected
                        .iter()
                        .map(|x| {
                            let numerator = x * BigInt::from(par.plaintext_modulus_u64().unwrap());
                            let rounded =
                                numerator.signum() * ((numerator.abs() + (&q >> 1usize)) / &q);
                            (((rounded % &q) + &q) % &q).to_biguint().unwrap()
                        })
                        .collect();
                    assert_eq!(Vec::<BigUint>::from(&part.to_power_basis()), expected);
                }
            }
        }
        Ok(())
    }

    #[test]
    fn decrypts_product_sums_at_multiple_levels() -> Result<()> {
        let par = Parameters::test_parameters(3, 16);
        let mut rng = ChaCha8Rng::seed_from_u64(0xf053d);
        let sk = SecretKey::generate(&par, &mut rng);
        for level in 0..=par.max_level() {
            let encoding = Encoding::Simd;
            let rk = RelinearizationKey::new_leveled(
                &sk,
                level,
                level.min(par.max_level() - 1),
                &mut rng,
            )?;
            for count in [1, 2, 17, 81] {
                let mut accumulator = CiphertextProductAccumulator::new(&par, level)?;
                let mut expected = vec![0u64; par.degree()];
                let mut separate: Option<Ciphertext> = None;
                for pair in 0..count {
                    let values: Vec<_> = (0..2)
                        .map(|side| {
                            (0..par.degree())
                                .map(|i| {
                                    (17 * i + 5 * pair + 31 * side) as u64
                                        % par.plaintext_modulus_u64().unwrap()
                                })
                                .collect::<Vec<_>>()
                        })
                        .collect();
                    let ct = values
                        .iter()
                        .map(|v| {
                            let pt =
                                Plaintext::encode_at_level(&par, v.as_slice(), encoding, level)?;
                            sk.encrypt(&pt, &mut rng)
                        })
                        .collect::<Result<Vec<Ciphertext>>>()?;
                    accumulator.add_product(&ct[0], &ct[1])?;
                    let product = ct[0].multiply(&ct[1]).unwrap();
                    if let Some(sum) = separate.as_mut() {
                        sum.add_assign(&product).unwrap();
                    } else {
                        separate = Some(product);
                    }
                    for (i, expected) in expected.iter_mut().enumerate() {
                        *expected = (*expected + values[0][i] * values[1][i])
                            % par.plaintext_modulus_u64().unwrap();
                    }
                }
                let separate = separate.unwrap();
                let mut result = accumulator.finish()?;
                if count == 1 {
                    assert_eq!(result, separate);
                }
                for ciphertext in [&result, &separate] {
                    assert_eq!(sk.decrypt(ciphertext)?.decode(encoding)?, expected);
                }
                rk.relinearize(&mut result)?;
                assert_eq!(sk.decrypt(&result)?.decode(encoding)?, expected);
            }
        }
        Ok(())
    }

    #[test]
    fn validates_inputs_and_preserves_sum_on_errors() -> Result<()> {
        let par = Parameters::test_parameters(2, 16);
        let other_par = Parameters::test_parameters(3, 16);
        let mut rng = ChaCha8Rng::seed_from_u64(31);
        let ctx = par.context_at_level(0)?;
        let ct = Ciphertext::from_components(vec![Poly::random(ctx, &mut rng); 2], &par)?;
        assert!(CiphertextProductAccumulator::new(&par, 2).is_err());
        assert_eq!(
            CiphertextProductAccumulator::new(&par, 0)?.finish(),
            Err(DotProductError::EmptyInput.into())
        );
        let mut accumulator = CiphertextProductAccumulator::new(&par, 0)?;
        accumulator.add_product(&ct, &ct)?;
        let saved = accumulator.c.clone();
        let wrong_level =
            Ciphertext::from_components(vec![Poly::zero(par.context_at_level(1)?); 2], &par)?;
        let wrong_par = Ciphertext::from_components(
            vec![Poly::zero(other_par.context_at_level(0)?); 2],
            &other_par,
        )?;
        let three_parts = Ciphertext::from_components(vec![Poly::zero(ctx); 3], &par)?;
        let mut wrong_context = ct.clone();
        wrong_context.c[1] = Poly::zero(par.context_at_level(1)?);
        for invalid in [
            Ciphertext::invalid_empty(&par),
            wrong_level,
            wrong_par,
            three_parts,
            wrong_context,
        ] {
            assert!(accumulator.add_product(&invalid, &ct).is_err());
            assert!(accumulator.add_product(&ct, &invalid).is_err());
            assert_eq!(*accumulator.c, *saved);
            assert_eq!(accumulator.products, 1);
        }
        // Exercise the full-accumulator guard without allocating an infeasibly
        // large input. The coefficient-bound arithmetic is checked separately.
        accumulator.max_products = 1;
        assert_eq!(
            accumulator.add_product(&ct, &ct),
            Err(DotProductError::TooManyProducts { maximum: 1 }.into())
        );
        assert_eq!(*accumulator.c, *saved);
        assert_eq!(accumulator.finish()?, ct.multiply(&ct).unwrap());
        Ok(())
    }

    #[test]
    fn coefficient_bound_is_strict_and_capped() {
        let q = BigUint::from(17u32);
        let denominator = &q * &q * 16u32 * 4u32;
        for count in [0usize, 1, 2, 81, 65536] {
            let boundary = &denominator * count;
            assert_eq!(product_limit(&q, &(&boundary + 1u32), 16), count);
            if count != 0 {
                assert_eq!(product_limit(&q, &boundary, 16), count - 1);
            }
        }
        assert_eq!(
            product_limit(&q, &(BigUint::from(1u32) << 256), 16),
            usize::MAX
        );
    }

    #[test]
    fn timing_permission_includes_every_part_and_pair() -> Result<()> {
        let par = Parameters::test_parameters(2, 16);
        let ctx = par.context_at_level(0)?;
        let mut p = Poly::zero(ctx);
        p.allow_variable_time_computations(crate::VariableTime::new(
            crate::PublicData::assert_public(),
        ));
        let public = Ciphertext::from_components(vec![p; 2], &par)?;
        for secret_pair in 0..=2 {
            for secret_part in 0..4 {
                let mut accumulator = CiphertextProductAccumulator::new(&par, 0)?;
                for pair in 0..2 {
                    let mut operands = [public.clone(), public.clone()];
                    if pair == secret_pair {
                        operands[secret_part / 2].c[secret_part % 2]
                            .disallow_variable_time_computations();
                    }
                    accumulator.add_product(&operands[0], &operands[1])?;
                }
                assert!(
                    accumulator
                        .finish()?
                        .iter()
                        .all(|p| { p.allows_variable_time_computations() == (secret_pair == 2) })
                );
            }
        }
        Ok(())
    }
}
