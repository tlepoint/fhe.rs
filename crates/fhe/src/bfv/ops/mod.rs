//! Operations over ciphertexts

mod dot_product;
pub use dot_product::{DotProductScalarWorkspace, dot_product_scalar};

mod mul;
pub use mul::{Multiplicator, PreparedMultiplicand};

mod tensor;

mod product_accumulator;
pub use product_accumulator::CiphertextProductAccumulator;

use super::{Ciphertext, Plaintext};
use crate::Result;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use std::sync::Arc;
use tensor::Scratch;

impl Add<&Ciphertext> for &Ciphertext {
    type Output = Ciphertext;

    fn add(self, rhs: &Ciphertext) -> Ciphertext {
        assert!(Arc::ptr_eq(&self.par, &rhs.par));

        assert_eq!(self.level, rhs.level);
        assert_eq!(self.len(), rhs.len());

        let c = self
            .iter()
            .zip(rhs.iter())
            .map(|(c1i, c2i)| c1i + c2i)
            .collect::<Vec<_>>();
        Ciphertext {
            par: self.par.clone(),
            seed: None,
            c,
            level: self.level,
        }
    }
}

impl Add<&Ciphertext> for Ciphertext {
    type Output = Ciphertext;

    fn add(mut self, rhs: &Ciphertext) -> Ciphertext {
        self += rhs;
        self
    }
}

impl AddAssign<&Ciphertext> for Ciphertext {
    fn add_assign(&mut self, rhs: &Ciphertext) {
        assert!(Arc::ptr_eq(&self.par, &rhs.par));

        assert_eq!(self.level, rhs.level);
        assert_eq!(self.len(), rhs.len());
        self.iter_mut()
            .zip(rhs.iter())
            .for_each(|(c1i, c2i)| *c1i += c2i);
        self.seed = None
    }
}

impl Add<&Plaintext> for &Ciphertext {
    type Output = Ciphertext;

    fn add(self, rhs: &Plaintext) -> Ciphertext {
        let mut self_clone = self.clone();
        self_clone += rhs;
        self_clone
    }
}

impl Add<&Ciphertext> for &Plaintext {
    type Output = Ciphertext;

    fn add(self, rhs: &Ciphertext) -> Ciphertext {
        rhs + self
    }
}

impl AddAssign<&Plaintext> for Ciphertext {
    fn add_assign(&mut self, rhs: &Plaintext) {
        assert!(Arc::ptr_eq(&self.par, &rhs.par));
        assert_eq!(self.level, rhs.level());

        let poly = rhs.to_poly();
        self.c[0] += &poly;
        self.seed = None
    }
}

impl Add<&Plaintext> for Ciphertext {
    type Output = Ciphertext;

    fn add(mut self, rhs: &Plaintext) -> Ciphertext {
        self += rhs;
        self
    }
}

impl Sub<&Ciphertext> for &Ciphertext {
    type Output = Ciphertext;

    fn sub(self, rhs: &Ciphertext) -> Ciphertext {
        assert!(Arc::ptr_eq(&self.par, &rhs.par));

        assert_eq!(self.level, rhs.level);
        assert_eq!(self.len(), rhs.len());

        let c = self
            .iter()
            .zip(rhs.iter())
            .map(|(c1i, c2i)| c1i - c2i)
            .collect::<Vec<_>>();
        Ciphertext {
            par: self.par.clone(),
            seed: None,
            c,
            level: self.level,
        }
    }
}

impl Sub<&Ciphertext> for Ciphertext {
    type Output = Ciphertext;

    fn sub(mut self, rhs: &Ciphertext) -> Ciphertext {
        self -= rhs;
        self
    }
}

impl SubAssign<&Ciphertext> for Ciphertext {
    fn sub_assign(&mut self, rhs: &Ciphertext) {
        assert!(Arc::ptr_eq(&self.par, &rhs.par));

        assert_eq!(self.level, rhs.level);
        assert_eq!(self.len(), rhs.len());
        self.iter_mut()
            .zip(rhs.iter())
            .for_each(|(c1i, c2i)| *c1i -= c2i);
        self.seed = None
    }
}

impl Sub<&Plaintext> for &Ciphertext {
    type Output = Ciphertext;

    fn sub(self, rhs: &Plaintext) -> Ciphertext {
        let mut self_clone = self.clone();
        self_clone -= rhs;
        self_clone
    }
}

impl Sub<&Ciphertext> for &Plaintext {
    type Output = Ciphertext;

    fn sub(self, rhs: &Ciphertext) -> Ciphertext {
        -(rhs - self)
    }
}

impl SubAssign<&Plaintext> for Ciphertext {
    fn sub_assign(&mut self, rhs: &Plaintext) {
        assert!(Arc::ptr_eq(&self.par, &rhs.par));
        assert_eq!(self.level, rhs.level());

        let poly = rhs.to_poly();
        self.c[0] -= &poly;
        self.seed = None
    }
}

impl Sub<&Plaintext> for Ciphertext {
    type Output = Ciphertext;

    fn sub(mut self, rhs: &Plaintext) -> Ciphertext {
        self -= rhs;
        self
    }
}

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

impl MulAssign<&Plaintext> for Ciphertext {
    fn mul_assign(&mut self, rhs: &Plaintext) {
        assert!(Arc::ptr_eq(&self.par, &rhs.par));
        assert_eq!(self.level, rhs.level());
        self.iter_mut().for_each(|ci| *ci *= &rhs.poly_ntt);
        self.seed = None
    }
}

impl Mul<&Plaintext> for &Ciphertext {
    type Output = Ciphertext;

    fn mul(self, rhs: &Plaintext) -> Ciphertext {
        let mut self_clone = self.clone();
        self_clone *= rhs;
        self_clone
    }
}

impl Mul<&Plaintext> for Ciphertext {
    type Output = Ciphertext;

    fn mul(mut self, rhs: &Plaintext) -> Ciphertext {
        self *= rhs;
        self
    }
}

impl Mul<&Ciphertext> for &Ciphertext {
    type Output = Ciphertext;

    fn mul(self, rhs: &Ciphertext) -> Ciphertext {
        assert!(Arc::ptr_eq(&self.par, &rhs.par));
        self.try_mul(rhs).unwrap()
    }
}

impl Ciphertext {
    /// Multiply without relinearization, returning validation errors instead
    /// of panicking. Both ciphertexts must use the same parameter instance and
    /// level and contain at least two polynomial parts.
    ///
    /// The two-part case uses three pointwise products (Karatsuba). Multiplying
    /// an object by itself uses [`Self::square`] without comparing
    /// coefficients.
    pub fn try_mul(&self, rhs: &Ciphertext) -> Result<Ciphertext> {
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
    /// [`Multiplicator::square`] to also relinearize or switch
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
    use crate::bfv::{
        BfvParameters, Ciphertext, Encoding, Plaintext, SecretKey, encoding::EncodingEnum,
    };
    use fhe_traits::{FheDecoder, FheDecrypter, FheEncoder, FheEncrypter};
    use rand::rng;
    use std::error::Error;

    #[test]
    fn add() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();

        for params in [
            BfvParameters::default_arc(1, 16),
            BfvParameters::default_arc(6, 16),
        ] {
            let zero = Ciphertext::trivial_zero(&params, 0)?;
            let q = fhe_math::zq::Modulus::new(params.plaintext()).unwrap();
            for _ in 0..50 {
                let a = q.random_vec(params.degree(), &mut rng);
                let b = q.random_vec(params.degree(), &mut rng);
                let mut c = a.clone();
                q.add_vec(&mut c, &b);

                let sk = SecretKey::random(&params, &mut rng);

                for encoding in [Encoding::poly(), Encoding::simd()] {
                    let pt_a = Plaintext::try_encode(&a, encoding.clone(), &params)?;
                    let pt_b = Plaintext::try_encode(&b, encoding.clone(), &params)?;

                    let mut ct_a: Ciphertext = sk.try_encrypt(&pt_a, &mut rng)?;
                    assert_eq!(ct_a, &ct_a + &zero);
                    assert_eq!(ct_a, &zero + &ct_a);
                    let ct_b: Ciphertext = sk.try_encrypt(&pt_b, &mut rng)?;
                    let ct_c = &ct_a + &ct_b;
                    let ct_c_owned = ct_a.clone() + &ct_b;
                    ct_a += &ct_b;

                    let pt_c = sk.try_decrypt(&ct_c)?;
                    assert_eq!(Vec::<u64>::try_decode(&pt_c, encoding.clone())?, c);
                    assert_eq!(ct_c_owned, ct_c);
                    let pt_c = sk.try_decrypt(&ct_a)?;
                    assert_eq!(Vec::<u64>::try_decode(&pt_c, encoding.clone())?, c);
                }
            }
        }

        Ok(())
    }

    #[test]
    fn add_scalar() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();

        for params in [
            BfvParameters::default_arc(1, 16),
            BfvParameters::default_arc(6, 16),
        ] {
            let q = fhe_math::zq::Modulus::new(params.plaintext()).unwrap();
            for _ in 0..50 {
                let a = q.random_vec(params.degree(), &mut rng);
                let b = q.random_vec(params.degree(), &mut rng);
                let mut c = a.clone();
                q.add_vec(&mut c, &b);

                let sk = SecretKey::random(&params, &mut rng);

                for encoding in [Encoding::poly(), Encoding::simd()] {
                    let zero = Plaintext::zero(encoding.clone(), &params)?;
                    let pt_a = Plaintext::try_encode(&a, encoding.clone(), &params)?;
                    let pt_b = Plaintext::try_encode(&b, encoding.clone(), &params)?;

                    let mut ct_a: Ciphertext = sk.try_encrypt(&pt_a, &mut rng)?;
                    assert_eq!(
                        Vec::<u64>::try_decode(
                            &sk.try_decrypt(&(&ct_a + &zero))?,
                            encoding.clone()
                        )?,
                        a
                    );
                    assert_eq!(
                        Vec::<u64>::try_decode(
                            &sk.try_decrypt(&(&zero + &ct_a))?,
                            encoding.clone()
                        )?,
                        a
                    );
                    let ct_c = &ct_a + &pt_b;
                    let ct_c_owned = ct_a.clone() + &pt_b;
                    ct_a += &pt_b;

                    let pt_c = sk.try_decrypt(&ct_c)?;
                    assert_eq!(Vec::<u64>::try_decode(&pt_c, encoding.clone())?, c);
                    assert_eq!(ct_c_owned, ct_c);
                    let pt_c = sk.try_decrypt(&ct_a)?;
                    assert_eq!(Vec::<u64>::try_decode(&pt_c, encoding.clone())?, c);
                }
            }
        }

        Ok(())
    }

    #[test]
    fn sub() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for params in [
            BfvParameters::default_arc(1, 16),
            BfvParameters::default_arc(6, 16),
        ] {
            let zero = Ciphertext::trivial_zero(&params, 0)?;
            let q = fhe_math::zq::Modulus::new(params.plaintext()).unwrap();
            for _ in 0..50 {
                let a = q.random_vec(params.degree(), &mut rng);
                let mut a_neg = a.clone();
                q.neg_vec(&mut a_neg);
                let b = q.random_vec(params.degree(), &mut rng);
                let mut c = a.clone();
                q.sub_vec(&mut c, &b);

                let sk = SecretKey::random(&params, &mut rng);

                for encoding in [Encoding::poly(), Encoding::simd()] {
                    let pt_a = Plaintext::try_encode(&a, encoding.clone(), &params)?;
                    let pt_b = Plaintext::try_encode(&b, encoding.clone(), &params)?;

                    let mut ct_a: Ciphertext = sk.try_encrypt(&pt_a, &mut rng)?;
                    assert_eq!(ct_a, &ct_a - &zero);
                    assert_eq!(
                        Vec::<u64>::try_decode(
                            &sk.try_decrypt(&(&zero - &ct_a))?,
                            encoding.clone()
                        )?,
                        a_neg
                    );
                    let ct_b: Ciphertext = sk.try_encrypt(&pt_b, &mut rng)?;
                    let ct_c = &ct_a - &ct_b;
                    let ct_c_owned = ct_a.clone() - &ct_b;
                    ct_a -= &ct_b;

                    let pt_c = sk.try_decrypt(&ct_c)?;
                    assert_eq!(Vec::<u64>::try_decode(&pt_c, encoding.clone())?, c);
                    assert_eq!(ct_c_owned, ct_c);
                    let pt_c = sk.try_decrypt(&ct_a)?;
                    assert_eq!(Vec::<u64>::try_decode(&pt_c, encoding.clone())?, c);
                }
            }
        }

        Ok(())
    }

    #[test]
    fn sub_scalar() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for params in [
            BfvParameters::default_arc(1, 16),
            BfvParameters::default_arc(6, 16),
        ] {
            let q = fhe_math::zq::Modulus::new(params.plaintext()).unwrap();
            for _ in 0..50 {
                let a = q.random_vec(params.degree(), &mut rng);
                let mut a_neg = a.clone();
                q.neg_vec(&mut a_neg);
                let b = q.random_vec(params.degree(), &mut rng);
                let mut c = a.clone();
                q.sub_vec(&mut c, &b);

                let sk = SecretKey::random(&params, &mut rng);

                for encoding in [Encoding::poly(), Encoding::simd()] {
                    let zero = Plaintext::zero(encoding.clone(), &params)?;
                    let pt_a = Plaintext::try_encode(&a, encoding.clone(), &params)?;
                    let pt_b = Plaintext::try_encode(&b, encoding.clone(), &params)?;

                    let mut ct_a: Ciphertext = sk.try_encrypt(&pt_a, &mut rng)?;
                    assert_eq!(
                        Vec::<u64>::try_decode(
                            &sk.try_decrypt(&(&ct_a - &zero))?,
                            encoding.clone()
                        )?,
                        a
                    );
                    assert_eq!(
                        Vec::<u64>::try_decode(
                            &sk.try_decrypt(&(&zero - &ct_a))?,
                            encoding.clone()
                        )?,
                        a_neg
                    );
                    let ct_c = &ct_a - &pt_b;
                    let ct_c_owned = ct_a.clone() - &pt_b;
                    ct_a -= &pt_b;

                    let pt_c = sk.try_decrypt(&ct_c)?;
                    assert_eq!(Vec::<u64>::try_decode(&pt_c, encoding.clone())?, c);
                    assert_eq!(ct_c_owned, ct_c);
                    let pt_c = sk.try_decrypt(&ct_a)?;
                    assert_eq!(Vec::<u64>::try_decode(&pt_c, encoding.clone())?, c);
                }
            }
        }

        Ok(())
    }

    #[test]
    fn neg() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for params in [
            BfvParameters::default_arc(1, 16),
            BfvParameters::default_arc(6, 16),
        ] {
            let q = fhe_math::zq::Modulus::new(params.plaintext()).unwrap();
            for _ in 0..50 {
                let a = q.random_vec(params.degree(), &mut rng);
                let mut c = a.clone();
                q.neg_vec(&mut c);

                let sk = SecretKey::random(&params, &mut rng);
                for encoding in [Encoding::poly(), Encoding::simd()] {
                    let pt_a = Plaintext::try_encode(&a, encoding.clone(), &params)?;

                    let ct_a: Ciphertext = sk.try_encrypt(&pt_a, &mut rng)?;

                    let ct_c = -&ct_a;
                    let pt_c = sk.try_decrypt(&ct_c)?;
                    assert_eq!(Vec::<u64>::try_decode(&pt_c, encoding.clone())?, c);

                    let ct_c = -ct_a;
                    let pt_c = sk.try_decrypt(&ct_c)?;
                    assert_eq!(Vec::<u64>::try_decode(&pt_c, encoding.clone())?, c);
                }
            }
        }

        Ok(())
    }

    #[test]
    fn mul_scalar() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();

        for params in [
            BfvParameters::default_arc(1, 16),
            BfvParameters::default_arc(6, 16),
        ] {
            let q = fhe_math::zq::Modulus::new(params.plaintext()).unwrap();
            for _ in 0..50 {
                let a = q.random_vec(params.degree(), &mut rng);
                let b = q.random_vec(params.degree(), &mut rng);

                let sk = SecretKey::random(&params, &mut rng);
                for encoding in [Encoding::poly(), Encoding::simd()] {
                    let mut c = vec![0u64; params.degree()];
                    match encoding.encoding {
                        EncodingEnum::Poly => {
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
                        EncodingEnum::Simd => {
                            c.clone_from(&a);
                            q.mul_vec(&mut c, &b);
                        }
                    }

                    let pt_a = Plaintext::try_encode(&a, encoding.clone(), &params)?;
                    let pt_b = Plaintext::try_encode(&b, encoding.clone(), &params)?;

                    let mut ct_a: Ciphertext = sk.try_encrypt(&pt_a, &mut rng)?;
                    let ct_c = &ct_a * &pt_b;
                    let ct_c_owned = ct_a.clone() * &pt_b;
                    ct_a *= &pt_b;

                    let pt_c = sk.try_decrypt(&ct_c)?;
                    assert_eq!(Vec::<u64>::try_decode(&pt_c, encoding.clone())?, c);
                    assert_eq!(ct_c_owned, ct_c);
                    let pt_c = sk.try_decrypt(&ct_a)?;
                    assert_eq!(Vec::<u64>::try_decode(&pt_c, encoding.clone())?, c);
                }
            }
        }

        Ok(())
    }

    #[test]
    fn mul() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for par in [
            BfvParameters::default_arc(2, 16),
            BfvParameters::default_arc(8, 16),
        ] {
            let q = fhe_math::zq::Modulus::new(par.plaintext()).unwrap();
            for _ in 0..1 {
                // We will encode `values` in an Simd format, and check that the product is
                // computed correctly.
                let v1 = q.random_vec(par.degree(), &mut rng);
                let v2 = q.random_vec(par.degree(), &mut rng);
                let mut expected = v1.clone();
                q.mul_vec(&mut expected, &v2);

                let sk = SecretKey::random(&par, &mut rng);
                let pt1 = Plaintext::try_encode(&v1, Encoding::simd(), &par)?;
                let pt2 = Plaintext::try_encode(&v2, Encoding::simd(), &par)?;

                let ct1: Ciphertext = sk.try_encrypt(&pt1, &mut rng)?;
                let ct2: Ciphertext = sk.try_encrypt(&pt2, &mut rng)?;
                let ct3 = &ct1 * &ct2;
                let ct4 = &ct3 * &ct3;
                assert!(
                    ct3.iter()
                        .chain(ct4.iter())
                        .all(|poly| poly.allows_variable_time_computations())
                );

                let mut mixed = ct2.clone();
                mixed.c[0].disallow_variable_time_computations();
                let mixed_product = &ct1 * &mixed;
                assert!(
                    mixed_product
                        .iter()
                        .all(|poly| !poly.allows_variable_time_computations())
                );

                println!(
                    "Noise: {}",
                    sk.measure_noise_vartime(
                        &ct3,
                        fhe_traits::SecretDependentDiagnostics::acknowledge_leakage()
                    )?
                );
                let pt = sk.try_decrypt(&ct3)?;
                assert_eq!(Vec::<u64>::try_decode(&pt, Encoding::simd())?, expected);

                let e = expected.clone();
                q.mul_vec(&mut expected, &e);
                println!(
                    "Noise: {}",
                    sk.measure_noise_vartime(
                        &ct4,
                        fhe_traits::SecretDependentDiagnostics::acknowledge_leakage()
                    )?
                );
                let pt = sk.try_decrypt(&ct4)?;
                assert_eq!(Vec::<u64>::try_decode(&pt, Encoding::simd())?, expected);
            }
        }
        Ok(())
    }

    #[test]
    fn square() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let par = BfvParameters::default_arc(6, 16);
        let q = fhe_math::zq::Modulus::new(par.plaintext()).unwrap();
        for _ in 0..20 {
            // We will encode `values` in an Simd format, and check that the product is
            // computed correctly.
            let v = q.random_vec(par.degree(), &mut rng);
            let mut expected = v.clone();
            q.mul_vec(&mut expected, &v);

            let sk = SecretKey::random(&par, &mut rng);
            let pt = Plaintext::try_encode(&v, Encoding::simd(), &par)?;

            let ct1: Ciphertext = sk.try_encrypt(&pt, &mut rng)?;
            let ct2 = &ct1 * &ct1;

            println!(
                "Noise: {}",
                sk.measure_noise_vartime(
                    &ct2,
                    fhe_traits::SecretDependentDiagnostics::acknowledge_leakage()
                )?
            );
            let pt = sk.try_decrypt(&ct2)?;
            assert_eq!(Vec::<u64>::try_decode(&pt, Encoding::simd())?, expected);
        }
        Ok(())
    }

    #[test]
    fn zero_multiplication_is_symmetric_at_every_level() -> Result<(), Box<dyn Error>> {
        let params = BfvParameters::default_arc(2, 16);
        let mut rng = rng();
        let sk = SecretKey::random(&params, &mut rng);
        for level in 0..=params.max_level() {
            let zero = Ciphertext::trivial_zero(&params, level)?;
            let encoding = Encoding::poly_at_level(level);
            let pt = Plaintext::try_encode(&[3u64], encoding.clone(), &params)?;
            let ct: Ciphertext = sk.try_encrypt(&pt, &mut rng)?;
            for operand in [zero.clone(), ct.clone(), &ct * &ct] {
                let left = &zero * &operand;
                let right = &operand * &zero;
                assert_eq!(left, right);
                assert_eq!(left.component_count(), operand.component_count() + 1);
                assert_eq!(left.level(), level);
                assert_eq!(
                    Vec::<u64>::try_decode(&sk.try_decrypt(&left)?, encoding.clone())?,
                    vec![0; params.degree()],
                );
            }
        }
        Ok(())
    }
}
