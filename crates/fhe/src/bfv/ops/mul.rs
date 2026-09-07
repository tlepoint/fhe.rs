use std::sync::Arc;

use fhe_math::{
    rns::ScalingFactor,
    rq::{Context, Ntt, NttShoup, Poly, scaler::Scaler},
    zq::primes::generate_prime,
};

use zeroize::Zeroizing;

use super::tensor::{self, Scratch};
use crate::{
    Error, ParametersError, Result,
    bfv::{Ciphertext, Parameters, keys::RelinearizationKey},
};

/// MultiplicationPlan that implements a strategy for multiplying. In
/// particular, the following information can be specified:
/// - Whether `lhs` must be scaled;
/// - Whether `rhs` must be scaled;
/// - The basis at which the multiplication will occur;
/// - The scaling factor after multiplication;
/// - Whether relinearization should be used.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MultiplicationPlan {
    par: Parameters,
    pub(crate) extender_lhs: Scaler,
    pub(crate) extender_rhs: Scaler,
    pub(crate) down_scaler: Scaler,
    pub(crate) base_ctx: Arc<Context>,
    pub(crate) mul_ctx: Arc<Context>,
    rk: Option<RelinearizationKey>,
    mod_switch: bool,
    level: usize,
    symmetric: bool,
}

impl MultiplicationPlan {
    /// Construct a multiplicator using custom scaling factors and extended
    /// basis.
    pub fn new(
        lhs_scaling_factor: ScalingFactor,
        rhs_scaling_factor: ScalingFactor,
        extended_basis: &[u64],
        post_mul_scaling_factor: ScalingFactor,
        par: &Parameters,
    ) -> Result<Self> {
        Self::new_leveled_internal(
            lhs_scaling_factor,
            rhs_scaling_factor,
            extended_basis,
            post_mul_scaling_factor,
            0,
            par,
        )
    }

    /// Construct a multiplicator using custom scaling factors and extended
    /// basis at a given level.
    pub fn new_leveled(
        lhs_scaling_factor: ScalingFactor,
        rhs_scaling_factor: ScalingFactor,
        extended_basis: &[u64],
        post_mul_scaling_factor: ScalingFactor,
        level: usize,
        par: &Parameters,
    ) -> Result<Self> {
        Self::new_leveled_internal(
            lhs_scaling_factor,
            rhs_scaling_factor,
            extended_basis,
            post_mul_scaling_factor,
            level,
            par,
        )
    }

    fn new_leveled_internal(
        lhs_scaling_factor: ScalingFactor,
        rhs_scaling_factor: ScalingFactor,
        extended_basis: &[u64],
        post_mul_scaling_factor: ScalingFactor,
        level: usize,
        par: &Parameters,
    ) -> Result<Self> {
        let base_ctx = par.context_at_level(level)?;
        let mul_ctx = Arc::new(Context::new(extended_basis, par.degree())?);
        let symmetric = lhs_scaling_factor == rhs_scaling_factor;
        let extender_lhs = Scaler::new(base_ctx, &mul_ctx, lhs_scaling_factor)?;
        let extender_rhs = Scaler::new(base_ctx, &mul_ctx, rhs_scaling_factor)?;
        let down_scaler = Scaler::new(&mul_ctx, base_ctx, post_mul_scaling_factor)?;
        Ok(Self {
            par: par.clone(),
            extender_lhs,
            extender_rhs,
            down_scaler,
            base_ctx: base_ctx.clone(),
            mul_ctx,
            rk: None,
            mod_switch: false,
            level,
            symmetric,
        })
    }

    /// Default multiplication strategy using relinearization.
    pub fn with_relinearization(rk: &RelinearizationKey) -> Result<Self> {
        let ctx = rk.ksk.par.context_at_level(rk.ksk.ciphertext_level)?;

        let modulus_size = rk.ksk.par.moduli_sizes()[..ctx.moduli().len()]
            .iter()
            .sum::<usize>();
        let n_moduli = (modulus_size + 60).div_ceil(62);

        let mut extended_basis = Vec::with_capacity(ctx.moduli().len() + n_moduli);
        extended_basis.append(&mut ctx.moduli().to_vec());
        let mut upper_bound = 1 << 62;
        while extended_basis.len() != ctx.moduli().len() + n_moduli {
            upper_bound = generate_prime(62, 2 * rk.ksk.par.degree() as u64, upper_bound)
                .ok_or_else(|| {
                    Error::ParametersError(ParametersError::NotEnoughPrimes {
                        size: 62,
                        degree: rk.ksk.par.degree(),
                        needed: n_moduli,
                        available: extended_basis.len() - ctx.moduli().len(),
                    })
                })?;
            if !extended_basis.contains(&upper_bound) && !ctx.moduli().contains(&upper_bound) {
                extended_basis.push(upper_bound)
            }
        }

        let mut multiplicator = Self::new_leveled_internal(
            ScalingFactor::one(),
            ScalingFactor::one(),
            &extended_basis,
            ScalingFactor::new(rk.ksk.par.plaintext_modulus(), ctx.modulus()),
            rk.ksk.ciphertext_level,
            &rk.ksk.par,
        )?;

        multiplicator.enable_relinearization(rk)?;
        Ok(multiplicator)
    }

    /// Enable relinearization after multiplication.
    pub fn enable_relinearization(&mut self, rk: &RelinearizationKey) -> Result<()> {
        if !Parameters::compatible(&self.par, &rk.ksk.par)
            || rk.ksk.ciphertext_level != self.level
            || rk.ksk.ctx_ciphertext != self.base_ctx
        {
            return Err(Error::ParameterMismatch {
                left: crate::ParameterSource::RelinearizationKey,
                right: crate::ParameterSource::MultiplicationPlan,
            });
        }
        self.rk = Some(rk.clone());
        Ok(())
    }

    /// Enable modulus switching after multiplication (and relinearization, if
    /// applicable).
    pub fn enable_mod_switching(&mut self) -> Result<()> {
        if self.par.context_at_level(self.par.max_level())? == &self.base_ctx {
            Err(fhe_math::Error::NoMoreContext.into())
        } else {
            self.mod_switch = true;
            Ok(())
        }
    }

    /// Use the parameters' precomputed multiplication basis at `level`, without
    /// relinearization or modulus switching. These can be enabled afterwards.
    pub fn without_relinearization(par: &Parameters, level: usize) -> Result<Self> {
        let mp = par.context_level_at(level)?.mul_params();
        Ok(Self {
            par: par.clone(),
            extender_lhs: mp.extender.clone(),
            extender_rhs: mp.extender.clone(),
            down_scaler: mp.down_scaler.clone(),
            base_ctx: mp.from.clone(),
            mul_ctx: mp.to.clone(),
            rk: None,
            mod_switch: false,
            level,
            symmetric: true,
        })
    }

    fn validate_operand(&self, ct: &Ciphertext) -> Result<()> {
        ct.validate_for_context(&self.par, self.level, &self.base_ctx)?;
        if ct.len() != 2 {
            return Err(crate::CiphertextError::MultiplicationPolynomialCount {
                left: ct.len(),
                right: ct.len(),
                expected: 2,
            }
            .into());
        }
        Ok(())
    }

    /// Multiply two two-part ciphertexts using the defined strategy.
    pub fn multiply(&self, lhs: &Ciphertext, rhs: &Ciphertext) -> Result<Ciphertext> {
        if std::ptr::eq(lhs, rhs) && self.symmetric {
            return self.square(lhs);
        }
        self.validate_operand(lhs)?;
        self.validate_operand(rhs)?;
        let left = Scratch([
            lhs.c[0].scale(&self.extender_lhs)?,
            lhs.c[1].scale(&self.extender_lhs)?,
        ]);
        let right = Scratch([
            rhs.c[0].scale(&self.extender_rhs)?,
            rhs.c[1].scale(&self.extender_rhs)?,
        ]);
        self.finish_product(tensor::product(&left.0, &right.0))
    }

    /// Square a two-part ciphertext and apply the configured relinearization
    /// and modulus switching. Symmetric strategies reuse one basis extension
    /// and compute the cross product once. Asymmetric custom scaling factors
    /// retain their separate left and right extensions.
    pub fn square(&self, ct: &Ciphertext) -> Result<Ciphertext> {
        self.validate_operand(ct)?;
        if !self.symmetric {
            return self.multiply(ct, ct);
        }
        let input = Scratch([
            ct.c[0].scale(&self.extender_lhs)?,
            ct.c[1].scale(&self.extender_lhs)?,
        ]);
        self.finish_product(tensor::square(&input.0))
    }

    /// Prepare a reusable left operand for this multiplication strategy.
    ///
    /// Caches its two extended NTT parts and their sum with Shoup quotients.
    /// Subsequent products skip the left basis extension and use cheaper
    /// multiplication by fixed polynomials. Preparation is useful when the
    /// same ciphertext participates in several products.
    ///
    /// The result owns a snapshot of the operand and borrows this strategy,
    /// preventing changes to its level, scalers, or relinearization settings
    /// while the prepared value is in use. Cached coefficients are zeroized
    /// on drop. Their storage is approximately `6*N*L` words, where `L` is the
    /// number of primes in the extended multiplication basis.
    pub fn prepare_lhs(&self, ct: &Ciphertext) -> Result<PreparedMultiplicand<'_>> {
        self.validate_operand(ct)?;
        let mut c0 = ct.c[0].scale(&self.extender_lhs)?;
        let mut c1 = ct.c[1].scale(&self.extender_lhs)?;
        if !ct.iter().all(Poly::allows_variable_time_computations) {
            c0.disallow_variable_time_computations();
            c1.disallow_variable_time_computations();
        }
        let sum = &c0 + &c1;
        Ok(PreparedMultiplicand {
            multiplicator: self,
            c: Zeroizing::new([
                c0.into_ntt_shoup(),
                c1.into_ntt_shoup(),
                sum.into_ntt_shoup(),
            ]),
        })
    }

    fn finish_product(&self, product: Vec<Poly<Ntt>>) -> Result<Ciphertext> {
        let product = Scratch(product);
        let mut c = product
            .0
            .iter()
            .map(|p| p.scale(&self.down_scaler))
            .collect::<fhe_math::Result<Vec<_>>>()?;

        // Relinearize
        if let Some(rk) = self.rk.as_ref() {
            let c2_pb = c[2].clone().into_power_basis();
            let (mut c0r, mut c1r) = rk.relinearizes_poly(&c2_pb)?;

            if c0r.ctx() != c[0].ctx() {
                let mut c0r_pb = c0r.into_power_basis();
                let mut c1r_pb = c1r.into_power_basis();
                c0r_pb.switch_down_to(c[0].ctx())?;
                c1r_pb.switch_down_to(c[1].ctx())?;
                c0r = c0r_pb.into_ntt();
                c1r = c1r_pb.into_ntt();
            }

            c[0] += &c0r;
            c[1] += &c1r;
            c.truncate(2);
        }

        let mut c = Ciphertext {
            par: self.par.clone(),
            seed: None,
            c,
            level: self.level,
        };

        if self.mod_switch {
            c.switch_down()?;
        }

        Ok(c)
    }
}

/// A reusable, basis-extended left operand tied to a [`MultiplicationPlan`].
/// Construct with [`MultiplicationPlan::prepare_lhs`].
///
/// This is an owned snapshot, so later changes to the source ciphertext do not
/// affect its products. Each multiplication still rounds independently; use
/// [`crate::bfv::CiphertextProductAccumulator`] to round a sum of products
/// once.
pub struct PreparedMultiplicand<'a> {
    multiplicator: &'a MultiplicationPlan,
    c: Zeroizing<[Poly<NttShoup>; 3]>,
}

impl PreparedMultiplicand<'_> {
    /// Multiply the prepared left operand by `rhs`, including the strategy's
    /// optional relinearization and modulus switching. `rhs` must have two
    /// parts and use compatible parameter settings and the strategy's input
    /// level.
    pub fn multiply(&self, rhs: &Ciphertext) -> Result<Ciphertext> {
        let strategy = self.multiplicator;
        strategy.validate_operand(rhs)?;
        let mut right = Scratch([
            rhs.c[0].scale(&strategy.extender_rhs)?,
            rhs.c[1].scale(&strategy.extender_rhs)?,
        ]);
        if !rhs.iter().all(Poly::allows_variable_time_computations) {
            for p in right.0.iter_mut() {
                p.disallow_variable_time_computations();
            }
        }
        let c0 = &right.0[0] * &self.c[0];
        let c2 = &right.0[1] * &self.c[1];
        let mut c1 = &right.0[0] + &right.0[1];
        c1 *= &self.c[2];
        c1 -= &c0;
        c1 -= &c2;
        strategy.finish_product(vec![c0, c1, c2])
    }
}

#[cfg(test)]
mod tests {
    use crate::bfv::{Ciphertext, Encoding, Parameters, Plaintext, RelinearizationKey, SecretKey};
    use fhe_math::{
        rns::{RnsContext, ScalingFactor},
        zq::primes::generate_prime,
    };

    use num_bigint::BigUint;
    use rand::rng;
    use std::error::Error;

    use super::MultiplicationPlan;

    #[test]
    fn prepared_products_and_squares_at_every_level() -> Result<(), Box<dyn Error>> {
        use rand::SeedableRng;
        use rand_chacha::ChaCha8Rng;
        use zeroize::Zeroize;

        let par = Parameters::test_parameters(3, 16);
        let mut rng = ChaCha8Rng::seed_from_u64(0xcac4e);
        let sk = SecretKey::generate(&par, &mut rng);
        for level in 0..=par.max_level() {
            let encoding = Encoding::Simd;
            let left: Vec<_> = (1..=par.degree() as u64).collect();
            let right: Vec<_> = left.iter().map(|x| 2 * x + 1).collect();
            let pt = Plaintext::encode_at_level(&par, &left, encoding, level)?;
            let other_pt = Plaintext::encode_at_level(&par, &right, encoding, level)?;
            let ct: Ciphertext = sk.encrypt(&pt, &mut rng)?;
            let other: Ciphertext = sk.encrypt(&other_pt, &mut rng)?;
            let base = MultiplicationPlan::without_relinearization(&par, level)?;
            let raw_product = ct.multiply(&other)?;
            let raw_square = ct.multiply(&ct.clone())?;
            assert_eq!(raw_square, ct.square()?);
            assert_eq!(base.multiply(&ct, &other)?, raw_product);
            for key_level in [0, level.min(par.max_level() - 1)] {
                for mode in 0..3 {
                    if mode == 2 && level == par.max_level() {
                        continue;
                    }
                    let mut strategy = base.clone();
                    let rk = RelinearizationKey::new_leveled(&sk, level, key_level, &mut rng)?;
                    let mut expected = raw_product.clone();
                    let mut expected_square = raw_square.clone();
                    if mode > 0 {
                        strategy.enable_relinearization(&rk)?;
                        rk.relinearize(&mut expected)?;
                        rk.relinearize(&mut expected_square)?;
                    }
                    if mode == 2 {
                        strategy.enable_mod_switching()?;
                        expected.switch_down()?;
                        expected_square.switch_down()?;
                    }
                    let mut snapshot_source = ct.clone();
                    let prepared = strategy.prepare_lhs(&snapshot_source)?;
                    snapshot_source.c[0].zeroize();
                    for _ in 0..2 {
                        let result = prepared.multiply(&other)?;
                        assert_eq!(result, expected);
                        assert_eq!(strategy.square(&ct)?, expected_square);
                        assert_eq!(strategy.multiply(&ct, &ct)?, expected_square);
                        assert_eq!(prepared.multiply(&ct)?, expected_square);
                        assert!(result.seed.is_none());
                        assert_eq!(
                            sk.decrypt(&result)?.decode(Encoding::Simd)?,
                            left.iter()
                                .zip(&right)
                                .map(|(x, y)| x * y % par.plaintext_modulus_u64().unwrap())
                                .collect::<Vec<_>>()
                        );
                        assert_eq!(
                            sk.decrypt(&expected_square)?.decode(Encoding::Simd)?,
                            left.iter()
                                .map(|x| x * x % par.plaintext_modulus_u64().unwrap())
                                .collect::<Vec<_>>()
                        );
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn prepared_products_preserve_all_timing_restrictions() -> Result<(), Box<dyn Error>> {
        use fhe_math::rq::{Ntt, Poly};
        use rand::SeedableRng;
        use rand_chacha::ChaCha8Rng;

        let par = Parameters::test_parameters(3, 16);
        let mut rng = ChaCha8Rng::seed_from_u64(0xc07);
        for level in 0..=par.max_level() {
            let ctx = par.context_at_level(level)?;
            let strategy = MultiplicationPlan::without_relinearization(&par, level)?;
            for mask in 0..16 {
                let mut parts: Vec<_> =
                    (0..4).map(|_| Poly::<Ntt>::random(ctx, &mut rng)).collect();
                for (i, p) in parts.iter_mut().enumerate() {
                    if mask & (1 << i) != 0 {
                        p.allow_variable_time_computations(crate::VariableTime::new(
                            crate::PublicData::assert_public(),
                        ));
                    }
                }
                let right = Ciphertext::from_components(parts.split_off(2), &par)?;
                let left = Ciphertext::from_components(parts, &par)?;
                let prepared = strategy.prepare_lhs(&left)?;
                let result = prepared.multiply(&right)?;
                assert_eq!(result, left.multiply(&right)?);
                assert!(
                    result
                        .iter()
                        .all(|p| p.allows_variable_time_computations() == (mask == 15))
                );
                assert_eq!(prepared.multiply(&left)?, strategy.square(&left)?);
            }
        }
        Ok(())
    }

    #[test]
    fn checked_multiplication_rejects_invalid_inputs() -> Result<(), Box<dyn Error>> {
        use fhe_math::rq::Poly;
        let par = Parameters::test_parameters(3, 16);
        let mut rng = rng();
        let sk = SecretKey::generate(&par, &mut rng);
        let pt = Plaintext::encode(&par, &[2u64], Encoding::Polynomial)?;
        let ct: Ciphertext = sk.encrypt(&pt, &mut rng)?;
        let strategy = MultiplicationPlan::without_relinearization(&par, 0)?;
        let prepared = strategy.prepare_lhs(&ct)?;
        let expected = prepared.multiply(&ct)?;

        let mut wrong_params = ct.clone();
        wrong_params.par = Parameters::builder()
            .degree(16)
            .plaintext_modulus(1153_u64)
            .ciphertext_moduli(par.moduli())
            .noise_variance(11)
            .build()?;
        let mut wrong_level = ct.clone();
        wrong_level.switch_down()?;
        let mut wrong_context = ct.clone();
        wrong_context.c[1] = Poly::zero(par.context_at_level(1)?);
        let mut one_part = ct.clone();
        one_part.c.truncate(1);
        let empty = Ciphertext::invalid_empty(&par);
        for invalid in [
            &wrong_params,
            &wrong_level,
            &wrong_context,
            &one_part,
            &empty,
        ] {
            let saved = invalid.clone();
            assert!(ct.multiply(invalid).is_err());
            assert!(strategy.multiply(&ct, invalid).is_err());
            assert!(strategy.multiply(invalid, &ct).is_err());
            assert!(strategy.prepare_lhs(invalid).is_err());
            assert!(strategy.square(invalid).is_err());
            assert!(prepared.multiply(invalid).is_err());
            assert_eq!(invalid, &saved);
        }
        for invalid in [&wrong_context, &one_part, &empty] {
            assert!(invalid.square().is_err());
        }
        assert!(strategy.prepare_lhs(&expected).is_err());
        assert!(prepared.multiply(&expected).is_err());
        assert_eq!(prepared.multiply(&ct)?, expected);
        assert!(MultiplicationPlan::without_relinearization(&par, par.max_level() + 1).is_err());
        // Identical moduli are insufficient: an evaluation key must belong
        // to compatible BFV parameter settings and the same ciphertext level.
        let foreign_sk = SecretKey::generate(&wrong_params.par, &mut rng);
        let foreign_rk = RelinearizationKey::new(&foreign_sk, &mut rng)?;
        let leveled_rk = RelinearizationKey::new_leveled(&sk, 1, 0, &mut rng)?;
        let mut strategy = strategy.clone();
        assert!(strategy.enable_relinearization(&foreign_rk).is_err());
        assert!(strategy.enable_relinearization(&leveled_rk).is_err());
        assert_eq!(strategy.square(&ct)?, expected);
        Ok(())
    }

    #[test]
    fn square_supports_unrelinearized_ciphertexts() -> Result<(), Box<dyn Error>> {
        let par = Parameters::test_parameters(3, 16);
        let mut rng = rng();
        let sk = SecretKey::generate(&par, &mut rng);
        let pt = Plaintext::encode(&par, &[2u64, 3], Encoding::Polynomial)?;
        let ct: Ciphertext = sk.encrypt(&pt, &mut rng)?;
        let product = ct.square()?;
        assert_eq!(product.len(), 3);
        let fourth_power = product.square()?;
        assert_eq!(fourth_power.len(), 5);
        assert_eq!(fourth_power, product.multiply(&product.clone())?);
        let mut expected = vec![0u64; par.degree()];
        expected[..5].copy_from_slice(&[16, 96, 216, 216, 81]);
        assert_eq!(
            sk.decrypt(&fourth_power)?.decode(Encoding::Polynomial)?,
            expected
        );
        Ok(())
    }

    #[test]
    fn mul() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let par = Parameters::test_parameters(3, 16);
        let q = fhe_math::zq::Modulus::new(par.plaintext_modulus_u64().unwrap()).unwrap();
        for _ in 0..30 {
            // We will encode `values` in an Simd format, and check that the product is
            // computed correctly.
            let values = q.random_vec(par.degree(), &mut rng);
            let mut expected = values.clone();
            q.mul_vec(&mut expected, &values);

            let sk = SecretKey::generate(&par, &mut rng);
            let rk = RelinearizationKey::new(&sk, &mut rng)?;
            let pt = Plaintext::encode(&par, &values, Encoding::Simd)?;
            let ct1 = sk.encrypt(&pt, &mut rng)?;
            let ct2 = sk.encrypt(&pt, &mut rng)?;

            let mut multiplicator = MultiplicationPlan::with_relinearization(&rk)?;
            let ct3 = multiplicator.multiply(&ct1, &ct2)?;
            println!(
                "Noise: {}",
                sk.measure_noise_vartime(
                    &ct3,
                    crate::SecretDependentDiagnostics::acknowledge_leakage()
                )?
            );
            let pt = sk.decrypt(&ct3)?;
            assert_eq!(pt.decode(Encoding::Simd)?, expected);

            multiplicator.enable_mod_switching()?;
            let ct3 = multiplicator.multiply(&ct1, &ct2)?;
            assert_eq!(ct3.level, 1);
            println!(
                "Noise: {}",
                sk.measure_noise_vartime(
                    &ct3,
                    crate::SecretDependentDiagnostics::acknowledge_leakage()
                )?
            );
            let pt = sk.decrypt(&ct3)?;
            assert_eq!(pt.decode(Encoding::Simd)?, expected);
        }
        Ok(())
    }

    #[test]
    fn mul_at_level() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let par = Parameters::test_parameters(3, 16);
        let q = fhe_math::zq::Modulus::new(par.plaintext_modulus_u64().unwrap()).unwrap();
        for _ in 0..15 {
            for level in 0..2 {
                let values = q.random_vec(par.degree(), &mut rng);
                let mut expected = values.clone();
                q.mul_vec(&mut expected, &values);

                let sk = SecretKey::generate(&par, &mut rng);
                let rk = RelinearizationKey::new_leveled(&sk, level, level, &mut rng)?;
                let pt = Plaintext::encode_at_level(&par, &values, Encoding::Simd, level)?;
                let ct1: Ciphertext = sk.encrypt(&pt, &mut rng)?;
                let ct2: Ciphertext = sk.encrypt(&pt, &mut rng)?;
                assert_eq!(ct1.level, level);
                assert_eq!(ct2.level, level);

                let mut multiplicator = MultiplicationPlan::with_relinearization(&rk).unwrap();
                let ct3 = multiplicator.multiply(&ct1, &ct2).unwrap();
                println!(
                    "Noise: {}",
                    sk.measure_noise_vartime(
                        &ct3,
                        crate::SecretDependentDiagnostics::acknowledge_leakage()
                    )?
                );
                let pt = sk.decrypt(&ct3)?;
                assert_eq!(pt.decode(Encoding::Simd)?, expected);

                multiplicator.enable_mod_switching()?;
                let ct3 = multiplicator.multiply(&ct1, &ct2)?;
                assert_eq!(ct3.level, level + 1);
                println!(
                    "Noise: {}",
                    sk.measure_noise_vartime(
                        &ct3,
                        crate::SecretDependentDiagnostics::acknowledge_leakage()
                    )?
                );
                let pt = sk.decrypt(&ct3)?;
                assert_eq!(pt.decode(Encoding::Simd)?, expected);
            }
        }
        Ok(())
    }

    #[test]
    fn mul_no_relin() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let par = Parameters::test_parameters(6, 16);
        let q = fhe_math::zq::Modulus::new(par.plaintext_modulus_u64().unwrap()).unwrap();
        for _ in 0..30 {
            // We will encode `values` in an Simd format, and check that the product is
            // computed correctly.
            let values = q.random_vec(par.degree(), &mut rng);
            let mut expected = values.clone();
            q.mul_vec(&mut expected, &values);

            let sk = SecretKey::generate(&par, &mut rng);
            let rk = RelinearizationKey::new(&sk, &mut rng)?;
            let pt = Plaintext::encode(&par, &values, Encoding::Simd)?;
            let ct1 = sk.encrypt(&pt, &mut rng)?;
            let ct2 = sk.encrypt(&pt, &mut rng)?;

            let mut multiplicator = MultiplicationPlan::with_relinearization(&rk)?;
            // Remove the relinearization key.
            multiplicator.rk = None;
            let ct3 = multiplicator.multiply(&ct1, &ct2)?;
            println!(
                "Noise: {}",
                sk.measure_noise_vartime(
                    &ct3,
                    crate::SecretDependentDiagnostics::acknowledge_leakage()
                )?
            );
            let pt = sk.decrypt(&ct3)?;
            assert_eq!(pt.decode(Encoding::Simd)?, expected);

            multiplicator.enable_mod_switching()?;
            let ct3 = multiplicator.multiply(&ct1, &ct2)?;
            assert_eq!(ct3.level, 1);
            println!(
                "Noise: {}",
                sk.measure_noise_vartime(
                    &ct3,
                    crate::SecretDependentDiagnostics::acknowledge_leakage()
                )?
            );
            let pt = sk.decrypt(&ct3)?;
            assert_eq!(pt.decode(Encoding::Simd)?, expected);
        }
        Ok(())
    }

    #[test]
    fn different_mul_strategy() -> Result<(), Box<dyn Error>> {
        // Implement the second multiplication strategy from <https://eprint.iacr.org/2021/204>

        let mut rng = rng();
        let par = Parameters::test_parameters(3, 16);
        let q = fhe_math::zq::Modulus::new(par.plaintext_modulus_u64().unwrap()).unwrap();
        let mut extended_basis = par.moduli().to_vec();
        extended_basis
            .push(generate_prime(62, 2 * par.degree() as u64, extended_basis[2]).unwrap());
        extended_basis
            .push(generate_prime(62, 2 * par.degree() as u64, extended_basis[3]).unwrap());
        extended_basis
            .push(generate_prime(62, 2 * par.degree() as u64, extended_basis[4]).unwrap());
        let rns = RnsContext::new(&extended_basis[3..])?;

        for _ in 0..30 {
            // We will encode `values` in an Simd format, and check that the product is
            // computed correctly.
            let values = q.random_vec(par.degree(), &mut rng);
            let mut expected = values.clone();
            q.mul_vec(&mut expected, &values);

            let sk = SecretKey::generate(&par, &mut rng);
            let pt = Plaintext::encode(&par, &values, Encoding::Simd)?;
            let ct1 = sk.encrypt(&pt, &mut rng)?;
            let ct2 = sk.encrypt(&pt, &mut rng)?;

            let mut multiplicator = MultiplicationPlan::new(
                ScalingFactor::one(),
                ScalingFactor::new(rns.modulus(), par.context_at_level(0)?.modulus()),
                &extended_basis,
                ScalingFactor::new(
                    &BigUint::from(par.plaintext_modulus_u64().unwrap()),
                    rns.modulus(),
                ),
                &par,
            )?;

            let ct3 = multiplicator.multiply(&ct1, &ct2)?;
            assert_eq!(multiplicator.prepare_lhs(&ct1)?.multiply(&ct2)?, ct3);
            let squared = multiplicator.square(&ct1)?;
            assert_eq!(squared, multiplicator.multiply(&ct1, &ct1.clone())?);
            assert_eq!(squared, multiplicator.prepare_lhs(&ct1)?.multiply(&ct1)?);
            assert_eq!(sk.decrypt(&squared)?.decode(Encoding::Simd)?, expected);
            println!(
                "Noise: {}",
                sk.measure_noise_vartime(
                    &ct3,
                    crate::SecretDependentDiagnostics::acknowledge_leakage()
                )?
            );
            let pt = sk.decrypt(&ct3)?;
            assert_eq!(pt.decode(Encoding::Simd)?, expected);

            multiplicator.enable_mod_switching()?;
            let ct3 = multiplicator.multiply(&ct1, &ct2)?;
            assert_eq!(ct3.level, 1);
            assert_eq!(multiplicator.prepare_lhs(&ct1)?.multiply(&ct2)?, ct3);
            let squared = multiplicator.square(&ct1)?;
            assert_eq!(squared, multiplicator.multiply(&ct1, &ct1.clone())?);
            assert_eq!(squared, multiplicator.prepare_lhs(&ct1)?.multiply(&ct1)?);
            assert_eq!(sk.decrypt(&squared)?.decode(Encoding::Simd)?, expected);
            println!(
                "Noise: {}",
                sk.measure_noise_vartime(
                    &ct3,
                    crate::SecretDependentDiagnostics::acknowledge_leakage()
                )?
            );
            let pt = sk.decrypt(&ct3)?;
            assert_eq!(pt.decode(Encoding::Simd)?, expected);
        }

        Ok(())
    }
}
