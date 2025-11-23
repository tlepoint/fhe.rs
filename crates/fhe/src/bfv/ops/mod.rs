//! Operations over ciphertexts

mod dot_product;
pub use dot_product::dot_product_scalar;

mod mul;
pub use mul::Multiplicator;

use super::{Ciphertext, Plaintext};
use crate::{Error, Result};
use fhe_math::rq::{Poly, Representation};
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use std::sync::Arc;

impl Add<&Ciphertext> for &Ciphertext {
    type Output = Ciphertext;

    fn add(self, rhs: &Ciphertext) -> Ciphertext {
        assert!(Arc::ptr_eq(&self.par, &rhs.par));

        if self.is_empty() {
            return rhs.clone();
        }
        if rhs.is_empty() {
            return self.clone();
        }

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

        if self.is_empty() {
            *self = rhs.clone()
        } else if !rhs.is_empty() {
            assert_eq!(self.level, rhs.level);
            assert_eq!(self.len(), rhs.len());
            self.iter_mut()
                .zip(rhs.iter())
                .for_each(|(c1i, c2i)| *c1i += c2i);
            self.seed = None
        }
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
        assert!(!self.is_empty());
        assert_eq!(self.level, rhs.level);

        let poly = rhs.to_poly();
        self[0] += &poly;
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

        if self.is_empty() {
            return -rhs.clone();
        }
        if rhs.is_empty() {
            return self.clone();
        }

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

        if self.is_empty() {
            *self = -rhs
        } else if !rhs.is_empty() {
            assert_eq!(self.level, rhs.level);
            assert_eq!(self.len(), rhs.len());
            self.iter_mut()
                .zip(rhs.iter())
                .for_each(|(c1i, c2i)| *c1i -= c2i);
            self.seed = None
        }
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
        assert!(!self.is_empty());
        assert_eq!(self.level, rhs.level);

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
        if !self.is_empty() {
            assert_eq!(self.level, rhs.level);
            self.iter_mut().for_each(|ci| *ci *= &rhs.poly_ntt);
        }
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
        if self.is_empty() {
            return self.clone();
        }

        if rhs == self {
            // Squaring operation
            let ctx_lvl = self.par.context_level_at(self.level).unwrap();
            let mp = ctx_lvl.mul_params();

            // Scale all ciphertexts
            let self_c = self
                .iter()
                .map(|ci| ci.scale(&mp.extender).map_err(Error::MathError))
                .collect::<Result<Vec<Poly>>>()
                .unwrap();

            // Multiply
            let mut c = vec![Poly::zero(&mp.to, Representation::Ntt); 2 * self_c.len() - 1];
            for i in 0..self_c.len() {
                for j in 0..self_c.len() {
                    c[i + j] += &(&self_c[i] * &self_c[j])
                }
            }

            // Scale
            let c = c
                .iter_mut()
                .map(|ci| {
                    ci.change_representation(Representation::PowerBasis);
                    let mut ci = ci.scale(&mp.down_scaler).map_err(Error::MathError)?;
                    ci.change_representation(Representation::Ntt);
                    Ok(ci)
                })
                .collect::<Result<Vec<Poly>>>()
                .unwrap();

            Ciphertext {
                par: self.par.clone(),
                seed: None,
                c,
                level: rhs.level,
            }
        } else {
            assert!(Arc::ptr_eq(&self.par, &rhs.par));
            assert_eq!(self.level, rhs.level);

            let ctx_lvl = self.par.context_level_at(self.level).unwrap();
            let mp = ctx_lvl.mul_params();

            // Scale all ciphertexts
            let self_c = self
                .iter()
                .map(|ci| ci.scale(&mp.extender).map_err(Error::MathError))
                .collect::<Result<Vec<Poly>>>()
                .unwrap();
            let other_c = rhs
                .iter()
                .map(|ci| ci.scale(&mp.extender).map_err(Error::MathError))
                .collect::<Result<Vec<Poly>>>()
                .unwrap();

            // Multiply
            let mut c =
                vec![Poly::zero(&mp.to, Representation::Ntt); self_c.len() + other_c.len() - 1];
            for i in 0..self_c.len() {
                for j in 0..other_c.len() {
                    c[i + j] += &(&self_c[i] * &other_c[j])
                }
            }

            // Scale
            let c = c
                .iter_mut()
                .map(|ci| {
                    ci.change_representation(Representation::PowerBasis);
                    let mut ci = ci.scale(&mp.down_scaler).map_err(Error::MathError)?;
                    ci.change_representation(Representation::Ntt);
                    Ok(ci)
                })
                .collect::<Result<Vec<Poly>>>()
                .unwrap();

            Ciphertext {
                par: self.par.clone(),
                seed: None,
                c,
                level: rhs.level,
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::bfv::{
        encoding::EncodingEnum, BfvParameters, Ciphertext, Encoding, Plaintext, SecretKey,
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
            let zero = Ciphertext::zero(&params);
            for _ in 0..50 {
                let a = params.plaintext.random_vec(params.degree(), &mut rng);
                let b = params.plaintext.random_vec(params.degree(), &mut rng);
                let mut c = a.clone();
                params.plaintext.add_vec(&mut c, &b);

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
            for _ in 0..50 {
                let a = params.plaintext.random_vec(params.degree(), &mut rng);
                let b = params.plaintext.random_vec(params.degree(), &mut rng);
                let mut c = a.clone();
                params.plaintext.add_vec(&mut c, &b);

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
            let zero = Ciphertext::zero(&params);
            for _ in 0..50 {
                let a = params.plaintext.random_vec(params.degree(), &mut rng);
                let mut a_neg = a.clone();
                params.plaintext.neg_vec(&mut a_neg);
                let b = params.plaintext.random_vec(params.degree(), &mut rng);
                let mut c = a.clone();
                params.plaintext.sub_vec(&mut c, &b);

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
            for _ in 0..50 {
                let a = params.plaintext.random_vec(params.degree(), &mut rng);
                let mut a_neg = a.clone();
                params.plaintext.neg_vec(&mut a_neg);
                let b = params.plaintext.random_vec(params.degree(), &mut rng);
                let mut c = a.clone();
                params.plaintext.sub_vec(&mut c, &b);

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
            for _ in 0..50 {
                let a = params.plaintext.random_vec(params.degree(), &mut rng);
                let mut c = a.clone();
                params.plaintext.neg_vec(&mut c);

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
            for _ in 0..50 {
                let a = params.plaintext.random_vec(params.degree(), &mut rng);
                let b = params.plaintext.random_vec(params.degree(), &mut rng);

                let sk = SecretKey::random(&params, &mut rng);
                for encoding in [Encoding::poly(), Encoding::simd()] {
                    let mut c = vec![0u64; params.degree()];
                    match encoding.encoding {
                        EncodingEnum::Poly => {
                            for i in 0..params.degree() {
                                for j in 0..params.degree() {
                                    if i + j >= params.degree() {
                                        c[(i + j) % params.degree()] = params.plaintext.sub(
                                            c[(i + j) % params.degree()],
                                            params.plaintext.mul(a[i], b[j]),
                                        );
                                    } else {
                                        c[i + j] = params
                                            .plaintext
                                            .add(c[i + j], params.plaintext.mul(a[i], b[j]));
                                    }
                                }
                            }
                        }
                        EncodingEnum::Simd => {
                            c.clone_from(&a);
                            params.plaintext.mul_vec(&mut c, &b);
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
            for _ in 0..1 {
                // We will encode `values` in an Simd format, and check that the product is
                // computed correctly.
                let v1 = par.plaintext.random_vec(par.degree(), &mut rng);
                let v2 = par.plaintext.random_vec(par.degree(), &mut rng);
                let mut expected = v1.clone();
                par.plaintext.mul_vec(&mut expected, &v2);

                let sk = SecretKey::random(&par, &mut rng);
                let pt1 = Plaintext::try_encode(&v1, Encoding::simd(), &par)?;
                let pt2 = Plaintext::try_encode(&v2, Encoding::simd(), &par)?;

                let ct1: Ciphertext = sk.try_encrypt(&pt1, &mut rng)?;
                let ct2: Ciphertext = sk.try_encrypt(&pt2, &mut rng)?;
                let ct3 = &ct1 * &ct2;
                let ct4 = &ct3 * &ct3;

                println!("Noise: {}", unsafe { sk.measure_noise(&ct3)? });
                let pt = sk.try_decrypt(&ct3)?;
                assert_eq!(Vec::<u64>::try_decode(&pt, Encoding::simd())?, expected);

                let e = expected.clone();
                par.plaintext.mul_vec(&mut expected, &e);
                println!("Noise: {}", unsafe { sk.measure_noise(&ct4)? });
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
        for _ in 0..20 {
            // We will encode `values` in an Simd format, and check that the product is
            // computed correctly.
            let v = par.plaintext.random_vec(par.degree(), &mut rng);
            let mut expected = v.clone();
            par.plaintext.mul_vec(&mut expected, &v);

            let sk = SecretKey::random(&par, &mut rng);
            let pt = Plaintext::try_encode(&v, Encoding::simd(), &par)?;

            let ct1: Ciphertext = sk.try_encrypt(&pt, &mut rng)?;
            let ct2 = &ct1 * &ct1;

            println!("Noise: {}", unsafe { sk.measure_noise(&ct2)? });
            let pt = sk.try_decrypt(&ct2)?;
            assert_eq!(Vec::<u64>::try_decode(&pt, Encoding::simd())?, expected);
        }
        Ok(())
    }
}
