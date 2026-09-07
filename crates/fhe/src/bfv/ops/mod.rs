//! Operations over ciphertexts

mod dot_product;
pub use dot_product::{DotProductScalarWorkspace, dot_product_scalar};

mod mul;
pub use mul::{MultiplicationPlan, PreparedMultiplicand};

mod tensor;

mod product_accumulator;
pub use product_accumulator::CiphertextProductAccumulator;

use super::{Ciphertext, Plaintext};
use crate::Result;
use std::ops::Neg;
use tensor::Scratch;

impl Neg for &Ciphertext {
    type Output = Ciphertext;

    fn neg(self) -> Ciphertext {
        let c = self.iter().map(|c1i| -c1i).collect::<Vec<_>>();
        Ciphertext {
            par: self.par.clone(),
            seed: None,
            c,
            level: self.level,
        }
    }
}

impl Neg for Ciphertext {
    type Output = Ciphertext;

    fn neg(mut self) -> Ciphertext {
        self.iter_mut().for_each(|c1i| *c1i = -&*c1i);
        self.seed = None;
        self
    }
}

impl Ciphertext {
    fn validate_binary(&self, rhs: &Self) -> Result<()> {
        self.validate_for(&self.par)?;
        rhs.validate_for_context(
            &self.par,
            self.level,
            self.par.context_at_level(self.level)?,
        )?;
        if self.c.len() != rhs.c.len() {
            return Err(crate::CiphertextError::ComponentCountMismatch {
                left: self.c.len(),
                right: rhs.c.len(),
            }
            .into());
        }
        Ok(())
    }

    fn validate_plaintext(&self, rhs: &Plaintext) -> Result<()> {
        self.validate_for(&self.par)?;
        rhs.validate_for_context(
            &self.par,
            self.level,
            self.par.context_at_level(self.level)?,
        )
    }
    /// Add a ciphertext with compatible parameters, level, and component count.
    pub fn add(&self, rhs: &Self) -> Result<Self> {
        let mut out = self.clone();
        out.add_assign(rhs)?;
        Ok(out)
    }
    /// Add in place. A validation error leaves this ciphertext unchanged.
    pub fn add_assign(&mut self, rhs: &Self) -> Result<()> {
        self.validate_binary(rhs)?;
        for (left, right) in self.c.iter_mut().zip(&rhs.c) {
            *left += right;
        }
        self.seed = None;
        Ok(())
    }
    /// Subtract a ciphertext with compatible parameters, level, and component
    /// count.
    pub fn subtract(&self, rhs: &Self) -> Result<Self> {
        let mut out = self.clone();
        out.subtract_assign(rhs)?;
        Ok(out)
    }
    /// Subtract in place. A validation error leaves this ciphertext unchanged.
    pub fn subtract_assign(&mut self, rhs: &Self) -> Result<()> {
        self.validate_binary(rhs)?;
        for (left, right) in self.c.iter_mut().zip(&rhs.c) {
            *left -= right;
        }
        self.seed = None;
        Ok(())
    }
    /// Add a plaintext with compatible parameters at the same level.
    pub fn add_plaintext(&self, rhs: &Plaintext) -> Result<Self> {
        let mut out = self.clone();
        out.add_plaintext_assign(rhs)?;
        Ok(out)
    }
    /// Add a plaintext in place. Errors leave this ciphertext unchanged.
    pub fn add_plaintext_assign(&mut self, rhs: &Plaintext) -> Result<()> {
        self.validate_plaintext(rhs)?;
        let poly = rhs.to_poly();
        self.c[0] += &poly;
        self.seed = None;
        Ok(())
    }
    /// Subtract a plaintext with compatible parameters at the same level.
    pub fn subtract_plaintext(&self, rhs: &Plaintext) -> Result<Self> {
        let mut out = self.clone();
        out.subtract_plaintext_assign(rhs)?;
        Ok(out)
    }
    /// Subtract a plaintext in place. Errors leave this ciphertext unchanged.
    pub fn subtract_plaintext_assign(&mut self, rhs: &Plaintext) -> Result<()> {
        self.validate_plaintext(rhs)?;
        let poly = rhs.to_poly();
        self.c[0] -= &poly;
        self.seed = None;
        Ok(())
    }
    /// Multiply a plaintext with compatible parameters at the same level.
    pub fn multiply_plaintext(&self, rhs: &Plaintext) -> Result<Self> {
        let mut out = self.clone();
        out.multiply_plaintext_assign(rhs)?;
        Ok(out)
    }
    /// Multiply a plaintext in place. Errors leave this ciphertext unchanged.
    pub fn multiply_plaintext_assign(&mut self, rhs: &Plaintext) -> Result<()> {
        self.validate_plaintext(rhs)?;
        for left in &mut self.c {
            *left *= &rhs.poly_ntt;
        }
        self.seed = None;
        Ok(())
    }
    /// Multiply without relinearization and replace this ciphertext only on
    /// success.
    pub fn multiply_assign(&mut self, rhs: &Self) -> Result<()> {
        *self = self.multiply(rhs)?;
        Ok(())
    }
    /// Square without relinearization and replace this ciphertext only on
    /// success.
    pub fn square_assign(&mut self) -> Result<()> {
        *self = self.square()?;
        Ok(())
    }
}

impl Ciphertext {
    /// Multiply without relinearization, returning validation errors instead
    /// of panicking. Both ciphertexts must use the compatible parameter
    /// settings and level and contain at least two polynomial parts.
    ///
    /// Multiplying `m` and `n` components produces `m + n - 1` components
    /// at the same level. Relinearization and switching are separate
    /// operations. The two-part case uses three pointwise products
    /// Karatsuba. Multiplying an object by itself uses [`Self::square`]
    /// without comparing coefficients.
    pub fn multiply(&self, rhs: &Ciphertext) -> Result<Ciphertext> {
        if std::ptr::eq(self, rhs) {
            return self.square();
        }
        self.validate_for(&self.par)?;
        let mp = self.par.context_level_at(self.level)?.mul_params();
        rhs.validate_for_context(&self.par, self.level, &mp.from)?;
        let left = Scratch(
            self.iter()
                .map(|p| p.scale(&mp.extender))
                .collect::<fhe_math::Result<Vec<_>>>()?,
        );
        let right = Scratch(
            rhs.iter()
                .map(|p| p.scale(&mp.extender))
                .collect::<fhe_math::Result<Vec<_>>>()?,
        );
        let product = Scratch(tensor::product(&left.0, &right.0));
        let c = product
            .0
            .iter()
            .map(|p| p.scale(&mp.down_scaler))
            .collect::<fhe_math::Result<_>>()?;
        Ciphertext::from_components(c, &self.par)
    }

    /// Square without relinearization. Each input part is extended once and
    /// each off-diagonal product is computed once, then doubled before BFV
    /// rounding. An input with `k` parts produces `2*k - 1` parts.
    ///
    /// Returns an error for an invalid ciphertext. Use
    /// [`MultiplicationPlan::square`] to also relinearize or switch
    /// down according to a configured strategy.
    pub fn square(&self) -> Result<Ciphertext> {
        self.validate_for(&self.par)?;
        let mp = self.par.context_level_at(self.level)?.mul_params();
        let input = Scratch(
            self.iter()
                .map(|p| p.scale(&mp.extender))
                .collect::<fhe_math::Result<Vec<_>>>()?,
        );
        let product = Scratch(tensor::square(&input.0));
        let c = product
            .0
            .iter()
            .map(|p| p.scale(&mp.down_scaler))
            .collect::<fhe_math::Result<_>>()?;
        Ciphertext::from_components(c, &self.par)
    }
}

#[cfg(test)]
mod tests {
    use crate::bfv::{Ciphertext, Encoding, Parameters, Plaintext, SecretKey};

    use rand::rng;
    use std::error::Error;

    #[test]
    fn add() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();

        for params in [
            Parameters::test_parameters(1, 16),
            Parameters::test_parameters(6, 16),
        ] {
            let zero = Ciphertext::trivial_zero(&params, 0)?;
            let q = fhe_math::zq::Modulus::new(params.plaintext_modulus_u64().unwrap()).unwrap();
            for _ in 0..50 {
                let a = q.random_vec(params.degree(), &mut rng);
                let b = q.random_vec(params.degree(), &mut rng);
                let mut c = a.clone();
                q.add_vec(&mut c, &b);

                let sk = SecretKey::generate(&params, &mut rng);

                for encoding in [Encoding::Polynomial, Encoding::Simd] {
                    let pt_a = Plaintext::encode(&params, &a, encoding)?;
                    let pt_b = Plaintext::encode(&params, &b, encoding)?;

                    let mut ct_a: Ciphertext = sk.encrypt(&pt_a, &mut rng)?;
                    assert_eq!(ct_a, ct_a.add(&zero).unwrap());
                    assert_eq!(ct_a, zero.add(&ct_a).unwrap());
                    let ct_b: Ciphertext = sk.encrypt(&pt_b, &mut rng)?;
                    let ct_c = ct_a.add(&ct_b).unwrap();
                    let ct_c_owned = (ct_a.clone()).add(&ct_b).unwrap();
                    ct_a.add_assign(&ct_b).unwrap();

                    let pt_c = sk.decrypt(&ct_c)?;
                    assert_eq!(pt_c.decode(encoding)?, c);
                    assert_eq!(ct_c_owned, ct_c);
                    let pt_c = sk.decrypt(&ct_a)?;
                    assert_eq!(pt_c.decode(encoding)?, c);
                }
            }
        }

        Ok(())
    }

    #[test]
    fn add_scalar() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();

        for params in [
            Parameters::test_parameters(1, 16),
            Parameters::test_parameters(6, 16),
        ] {
            let q = fhe_math::zq::Modulus::new(params.plaintext_modulus_u64().unwrap()).unwrap();
            for _ in 0..50 {
                let a = q.random_vec(params.degree(), &mut rng);
                let b = q.random_vec(params.degree(), &mut rng);
                let mut c = a.clone();
                q.add_vec(&mut c, &b);

                let sk = SecretKey::generate(&params, &mut rng);

                for encoding in [Encoding::Polynomial, Encoding::Simd] {
                    let zero = Plaintext::zero(&params, 0)?;
                    let pt_a = Plaintext::encode(&params, &a, encoding)?;
                    let pt_b = Plaintext::encode(&params, &b, encoding)?;

                    let mut ct_a: Ciphertext = sk.encrypt(&pt_a, &mut rng)?;
                    assert_eq!(
                        sk.decrypt(&(ct_a.add_plaintext(&zero).unwrap()))?
                            .decode(encoding)?,
                        a
                    );
                    assert_eq!(
                        sk.decrypt(&(ct_a.add_plaintext(&zero).unwrap()))?
                            .decode(encoding)?,
                        a
                    );
                    let ct_c = ct_a.add_plaintext(&pt_b).unwrap();
                    let ct_c_owned = (ct_a.clone()).add_plaintext(&pt_b).unwrap();
                    ct_a.add_plaintext_assign(&pt_b).unwrap();

                    let pt_c = sk.decrypt(&ct_c)?;
                    assert_eq!(pt_c.decode(encoding)?, c);
                    assert_eq!(ct_c_owned, ct_c);
                    let pt_c = sk.decrypt(&ct_a)?;
                    assert_eq!(pt_c.decode(encoding)?, c);
                }
            }
        }

        Ok(())
    }

    #[test]
    fn sub() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for params in [
            Parameters::test_parameters(1, 16),
            Parameters::test_parameters(6, 16),
        ] {
            let zero = Ciphertext::trivial_zero(&params, 0)?;
            let q = fhe_math::zq::Modulus::new(params.plaintext_modulus_u64().unwrap()).unwrap();
            for _ in 0..50 {
                let a = q.random_vec(params.degree(), &mut rng);
                let mut a_neg = a.clone();
                q.neg_vec(&mut a_neg);
                let b = q.random_vec(params.degree(), &mut rng);
                let mut c = a.clone();
                q.sub_vec(&mut c, &b);

                let sk = SecretKey::generate(&params, &mut rng);

                for encoding in [Encoding::Polynomial, Encoding::Simd] {
                    let pt_a = Plaintext::encode(&params, &a, encoding)?;
                    let pt_b = Plaintext::encode(&params, &b, encoding)?;

                    let mut ct_a: Ciphertext = sk.encrypt(&pt_a, &mut rng)?;
                    assert_eq!(ct_a, ct_a.subtract(&zero).unwrap());
                    assert_eq!(
                        sk.decrypt(&(zero.subtract(&ct_a).unwrap()))?
                            .decode(encoding)?,
                        a_neg
                    );
                    let ct_b: Ciphertext = sk.encrypt(&pt_b, &mut rng)?;
                    let ct_c = ct_a.subtract(&ct_b).unwrap();
                    let ct_c_owned = (ct_a.clone()).subtract(&ct_b).unwrap();
                    ct_a.subtract_assign(&ct_b).unwrap();

                    let pt_c = sk.decrypt(&ct_c)?;
                    assert_eq!(pt_c.decode(encoding)?, c);
                    assert_eq!(ct_c_owned, ct_c);
                    let pt_c = sk.decrypt(&ct_a)?;
                    assert_eq!(pt_c.decode(encoding)?, c);
                }
            }
        }

        Ok(())
    }

    #[test]
    fn sub_scalar() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for params in [
            Parameters::test_parameters(1, 16),
            Parameters::test_parameters(6, 16),
        ] {
            let q = fhe_math::zq::Modulus::new(params.plaintext_modulus_u64().unwrap()).unwrap();
            for _ in 0..50 {
                let a = q.random_vec(params.degree(), &mut rng);
                let mut a_neg = a.clone();
                q.neg_vec(&mut a_neg);
                let b = q.random_vec(params.degree(), &mut rng);
                let mut c = a.clone();
                q.sub_vec(&mut c, &b);

                let sk = SecretKey::generate(&params, &mut rng);

                for encoding in [Encoding::Polynomial, Encoding::Simd] {
                    let zero = Plaintext::zero(&params, 0)?;
                    let pt_a = Plaintext::encode(&params, &a, encoding)?;
                    let pt_b = Plaintext::encode(&params, &b, encoding)?;

                    let mut ct_a: Ciphertext = sk.encrypt(&pt_a, &mut rng)?;
                    assert_eq!(
                        (sk.decrypt(&(ct_a.subtract_plaintext(&zero).unwrap()))?)
                            .decode(encoding)?,
                        a
                    );
                    assert_eq!(
                        (sk.decrypt(&(-(ct_a.subtract_plaintext(&zero).unwrap())))?)
                            .decode(encoding)?,
                        a_neg
                    );
                    let ct_c = ct_a.subtract_plaintext(&pt_b).unwrap();
                    let ct_c_owned = (ct_a.clone()).subtract_plaintext(&pt_b).unwrap();
                    ct_a.subtract_plaintext_assign(&pt_b).unwrap();

                    let pt_c = sk.decrypt(&ct_c)?;
                    assert_eq!(pt_c.decode(encoding)?, c);
                    assert_eq!(ct_c_owned, ct_c);
                    let pt_c = sk.decrypt(&ct_a)?;
                    assert_eq!(pt_c.decode(encoding)?, c);
                }
            }
        }

        Ok(())
    }

    #[test]
    fn neg() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for params in [
            Parameters::test_parameters(1, 16),
            Parameters::test_parameters(6, 16),
        ] {
            let q = fhe_math::zq::Modulus::new(params.plaintext_modulus_u64().unwrap()).unwrap();
            for _ in 0..50 {
                let a = q.random_vec(params.degree(), &mut rng);
                let mut c = a.clone();
                q.neg_vec(&mut c);

                let sk = SecretKey::generate(&params, &mut rng);
                for encoding in [Encoding::Polynomial, Encoding::Simd] {
                    let pt_a = Plaintext::encode(&params, &a, encoding)?;

                    let ct_a: Ciphertext = sk.encrypt(&pt_a, &mut rng)?;

                    let ct_c = -&ct_a;
                    let pt_c = sk.decrypt(&ct_c)?;
                    assert_eq!(pt_c.decode(encoding)?, c);

                    let ct_c = -ct_a;
                    let pt_c = sk.decrypt(&ct_c)?;
                    assert_eq!(pt_c.decode(encoding)?, c);
                }
            }
        }

        Ok(())
    }

    #[test]
    fn mul_scalar() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();

        for params in [
            Parameters::test_parameters(1, 16),
            Parameters::test_parameters(6, 16),
        ] {
            let q = fhe_math::zq::Modulus::new(params.plaintext_modulus_u64().unwrap()).unwrap();
            for _ in 0..50 {
                let a = q.random_vec(params.degree(), &mut rng);
                let b = q.random_vec(params.degree(), &mut rng);

                let sk = SecretKey::generate(&params, &mut rng);
                for encoding in [Encoding::Polynomial, Encoding::Simd] {
                    let mut c = vec![0u64; params.degree()];
                    match encoding {
                        Encoding::Polynomial => {
                            for i in 0..params.degree() {
                                for j in 0..params.degree() {
                                    if i + j >= params.degree() {
                                        c[(i + j) % params.degree()] =
                                            q.sub(c[(i + j) % params.degree()], q.mul(a[i], b[j]));
                                    } else {
                                        c[i + j] = q.add(c[i + j], q.mul(a[i], b[j]));
                                    }
                                }
                            }
                        }
                        Encoding::Simd => {
                            c.clone_from(&a);
                            q.mul_vec(&mut c, &b);
                        }
                    }

                    let pt_a = Plaintext::encode(&params, &a, encoding)?;
                    let pt_b = Plaintext::encode(&params, &b, encoding)?;

                    let mut ct_a: Ciphertext = sk.encrypt(&pt_a, &mut rng)?;
                    let ct_c = ct_a.multiply_plaintext(&pt_b).unwrap();
                    let ct_c_owned = (ct_a.clone()).multiply_plaintext(&pt_b).unwrap();
                    ct_a.multiply_plaintext_assign(&pt_b).unwrap();

                    let pt_c = sk.decrypt(&ct_c)?;
                    assert_eq!(pt_c.decode(encoding)?, c);
                    assert_eq!(ct_c_owned, ct_c);
                    let pt_c = sk.decrypt(&ct_a)?;
                    assert_eq!(pt_c.decode(encoding)?, c);
                }
            }
        }

        Ok(())
    }

    #[test]
    fn mul() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for par in [
            Parameters::test_parameters(2, 16),
            Parameters::test_parameters(8, 16),
        ] {
            let q = fhe_math::zq::Modulus::new(par.plaintext_modulus_u64().unwrap()).unwrap();
            for _ in 0..1 {
                // We will encode `values` in an Simd format, and check that the product is
                // computed correctly.
                let v1 = q.random_vec(par.degree(), &mut rng);
                let v2 = q.random_vec(par.degree(), &mut rng);
                let mut expected = v1.clone();
                q.mul_vec(&mut expected, &v2);

                let sk = SecretKey::generate(&par, &mut rng);
                let pt1 = Plaintext::encode(&par, &v1, Encoding::Simd)?;
                let pt2 = Plaintext::encode(&par, &v2, Encoding::Simd)?;

                let ct1: Ciphertext = sk.encrypt(&pt1, &mut rng)?;
                let ct2: Ciphertext = sk.encrypt(&pt2, &mut rng)?;
                let ct3 = ct1.multiply(&ct2).unwrap();
                let ct4 = ct3.multiply(&ct3).unwrap();
                assert!(
                    ct3.iter()
                        .chain(ct4.iter())
                        .all(|poly| poly.allows_variable_time_computations())
                );

                let mut mixed = ct2.clone();
                mixed.c[0].disallow_variable_time_computations();
                let mixed_product = ct1.multiply(&mixed).unwrap();
                assert!(
                    mixed_product
                        .iter()
                        .all(|poly| !poly.allows_variable_time_computations())
                );

                println!(
                    "Noise: {}",
                    sk.measure_noise_vartime(
                        &ct3,
                        crate::SecretDependentDiagnostics::acknowledge_leakage()
                    )?
                );
                let pt = sk.decrypt(&ct3)?;
                assert_eq!(pt.decode(Encoding::Simd)?, expected);

                let e = expected.clone();
                q.mul_vec(&mut expected, &e);
                println!(
                    "Noise: {}",
                    sk.measure_noise_vartime(
                        &ct4,
                        crate::SecretDependentDiagnostics::acknowledge_leakage()
                    )?
                );
                let pt = sk.decrypt(&ct4)?;
                assert_eq!(pt.decode(Encoding::Simd)?, expected);
            }
        }
        Ok(())
    }

    #[test]
    fn square() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let par = Parameters::test_parameters(6, 16);
        let q = fhe_math::zq::Modulus::new(par.plaintext_modulus_u64().unwrap()).unwrap();
        for _ in 0..20 {
            // We will encode `values` in an Simd format, and check that the product is
            // computed correctly.
            let v = q.random_vec(par.degree(), &mut rng);
            let mut expected = v.clone();
            q.mul_vec(&mut expected, &v);

            let sk = SecretKey::generate(&par, &mut rng);
            let pt = Plaintext::encode(&par, &v, Encoding::Simd)?;

            let ct1: Ciphertext = sk.encrypt(&pt, &mut rng)?;
            let ct2 = ct1.multiply(&ct1).unwrap();

            println!(
                "Noise: {}",
                sk.measure_noise_vartime(
                    &ct2,
                    crate::SecretDependentDiagnostics::acknowledge_leakage()
                )?
            );
            let pt = sk.decrypt(&ct2)?;
            assert_eq!(pt.decode(Encoding::Simd)?, expected);
        }
        Ok(())
    }

    #[test]
    fn zero_multiplication_is_symmetric_at_every_level() -> Result<(), Box<dyn Error>> {
        let params = Parameters::test_parameters(2, 16);
        let mut rng = rng();
        let sk = SecretKey::generate(&params, &mut rng);
        for level in 0..=params.max_level() {
            let zero = Ciphertext::trivial_zero(&params, level)?;
            let encoding = Encoding::Polynomial;
            let pt = Plaintext::encode_at_level(&params, &[3u64], encoding, level)?;
            let ct: Ciphertext = sk.encrypt(&pt, &mut rng)?;
            for operand in [zero.clone(), ct.clone(), ct.multiply(&ct).unwrap()] {
                let left = zero.multiply(&operand).unwrap();
                let right = operand.multiply(&zero).unwrap();
                assert_eq!(left, right);
                assert_eq!(left.component_count(), operand.component_count() + 1);
                assert_eq!(left.level(), level);
                assert_eq!(
                    sk.decrypt(&left)?.decode(encoding)?,
                    vec![0; params.degree()],
                );
            }
        }
        Ok(())
    }
}
