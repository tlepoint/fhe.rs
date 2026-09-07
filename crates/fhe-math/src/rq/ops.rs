//! Implementation of operations over polynomials.

use super::{Ntt, NttShoup, Poly, PowerBasis};
use itertools::izip;
use num_bigint::BigUint;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

// These operations are shared by PowerBasis and Ntt only. Do not add
// NttShoup: mutating its coefficients would invalidate its cached quotients.
macro_rules! impl_coefficient_ops {
    ($repr:ty) => {
        impl AddAssign<&Poly<$repr>> for Poly<$repr> {
            fn add_assign(&mut self, p: &Poly<$repr>) {
                assert!(!self.has_lazy_coefficients && !p.has_lazy_coefficients);
                debug_assert_eq!(self.ctx, p.ctx, "Incompatible contexts");

                self.allow_variable_time_computations &= p.allow_variable_time_computations;
                if self.allow_variable_time_computations {
                    izip!(
                        self.coefficients.outer_iter_mut(),
                        p.coefficients.outer_iter(),
                        self.ctx.q.iter()
                    )
                    .for_each(|(mut v1, v2, qi)| unsafe {
                        qi.add_vec_vt(v1.as_slice_mut().unwrap(), v2.as_slice().unwrap())
                    });
                } else {
                    izip!(
                        self.coefficients.outer_iter_mut(),
                        p.coefficients.outer_iter(),
                        self.ctx.q.iter()
                    )
                    .for_each(|(mut v1, v2, qi)| {
                        qi.add_vec(v1.as_slice_mut().unwrap(), v2.as_slice().unwrap())
                    });
                }
            }
        }
        impl Add<&Poly<$repr>> for &Poly<$repr> {
            type Output = Poly<$repr>;
            fn add(self, p: &Poly<$repr>) -> Poly<$repr> {
                let mut q = self.clone();
                q += p;
                q
            }
        }
        impl Add for Poly<$repr> {
            type Output = Poly<$repr>;
            fn add(self, mut p: Poly<$repr>) -> Poly<$repr> {
                p += &self;
                p
            }
        }
        impl SubAssign<&Poly<$repr>> for Poly<$repr> {
            fn sub_assign(&mut self, p: &Poly<$repr>) {
                assert!(!self.has_lazy_coefficients && !p.has_lazy_coefficients);
                debug_assert_eq!(self.ctx, p.ctx, "Incompatible contexts");

                self.allow_variable_time_computations &= p.allow_variable_time_computations;
                if self.allow_variable_time_computations {
                    izip!(
                        self.coefficients.outer_iter_mut(),
                        p.coefficients.outer_iter(),
                        self.ctx.q.iter()
                    )
                    .for_each(|(mut v1, v2, qi)| unsafe {
                        qi.sub_vec_vt(v1.as_slice_mut().unwrap(), v2.as_slice().unwrap())
                    });
                } else {
                    izip!(
                        self.coefficients.outer_iter_mut(),
                        p.coefficients.outer_iter(),
                        self.ctx.q.iter()
                    )
                    .for_each(|(mut v1, v2, qi)| {
                        qi.sub_vec(v1.as_slice_mut().unwrap(), v2.as_slice().unwrap())
                    });
                }
            }
        }
        impl Sub<&Poly<$repr>> for &Poly<$repr> {
            type Output = Poly<$repr>;
            fn sub(self, p: &Poly<$repr>) -> Poly<$repr> {
                let mut q = self.clone();
                q -= p;
                q
            }
        }
        impl MulAssign<&BigUint> for Poly<$repr> {
            fn mul_assign(&mut self, p: &BigUint) {
                let scalar_crt = self.ctx.rns.project(p);

                if self.allow_variable_time_computations {
                    unsafe {
                        izip!(
                            self.coefficients.outer_iter_mut(),
                            scalar_crt.iter(),
                            self.ctx.q.iter()
                        )
                        .for_each(|(mut v1, scalar_qi, qi)| {
                            qi.scalar_mul_vec_vt(v1.as_slice_mut().unwrap(), *scalar_qi)
                        });
                    }
                } else {
                    izip!(
                        self.coefficients.outer_iter_mut(),
                        scalar_crt.iter(),
                        self.ctx.q.iter()
                    )
                    .for_each(|(mut v1, scalar_qi, qi)| {
                        qi.scalar_mul_vec(v1.as_slice_mut().unwrap(), *scalar_qi)
                    });
                }
            }
        }
        impl Neg for &Poly<$repr> {
            type Output = Poly<$repr>;
            fn neg(self) -> Self::Output {
                -self.clone()
            }
        }
        impl Neg for Poly<$repr> {
            type Output = Poly<$repr>;

            fn neg(mut self) -> Poly<$repr> {
                assert!(!self.has_lazy_coefficients);
                if self.allow_variable_time_computations {
                    izip!(self.coefficients.outer_iter_mut(), self.ctx.q.iter()).for_each(
                        |(mut v1, qi)| unsafe { qi.neg_vec_vt(v1.as_slice_mut().unwrap()) },
                    );
                } else {
                    izip!(self.coefficients.outer_iter_mut(), self.ctx.q.iter())
                        .for_each(|(mut v1, qi)| qi.neg_vec(v1.as_slice_mut().unwrap()));
                }
                self
            }
        }
    };
}
impl_coefficient_ops!(PowerBasis);
impl_coefficient_ops!(Ntt);

impl MulAssign<&Poly<Ntt>> for Poly<Ntt> {
    fn mul_assign(&mut self, p: &Poly<Ntt>) {
        assert!(!p.has_lazy_coefficients);
        assert!(
            !self.has_lazy_coefficients,
            "Cannot multiply lazy coefficients by an Ntt polynomial"
        );
        debug_assert_eq!(self.ctx, p.ctx, "Incompatible contexts");
        self.allow_variable_time_computations &= p.allow_variable_time_computations;

        if self.allow_variable_time_computations {
            unsafe {
                izip!(
                    self.coefficients.outer_iter_mut(),
                    p.coefficients.outer_iter(),
                    self.ctx.q.iter()
                )
                .for_each(|(mut v1, v2, qi)| {
                    qi.mul_vec_vt(v1.as_slice_mut().unwrap(), v2.as_slice().unwrap());
                });
            }
        } else {
            izip!(
                self.coefficients.outer_iter_mut(),
                p.coefficients.outer_iter(),
                self.ctx.q.iter()
            )
            .for_each(|(mut v1, v2, qi)| {
                qi.mul_vec(v1.as_slice_mut().unwrap(), v2.as_slice().unwrap())
            });
        }
    }
}

impl MulAssign<&Poly<NttShoup>> for Poly<Ntt> {
    fn mul_assign(&mut self, p: &Poly<NttShoup>) {
        assert!(!p.has_lazy_coefficients);
        debug_assert_eq!(self.ctx, p.ctx, "Incompatible contexts");
        self.allow_variable_time_computations &= p.allow_variable_time_computations;

        if self.allow_variable_time_computations {
            izip!(
                self.coefficients.outer_iter_mut(),
                p.coefficients.outer_iter(),
                p.coefficients_shoup.as_ref().unwrap().outer_iter(),
                self.ctx.q.iter()
            )
            .for_each(|(mut v1, v2, v2_shoup, qi)| unsafe {
                qi.mul_shoup_vec_vt(
                    v1.as_slice_mut().unwrap(),
                    v2.as_slice().unwrap(),
                    v2_shoup.as_slice().unwrap(),
                )
            });
        } else {
            izip!(
                self.coefficients.outer_iter_mut(),
                p.coefficients.outer_iter(),
                p.coefficients_shoup.as_ref().unwrap().outer_iter(),
                self.ctx.q.iter()
            )
            .for_each(|(mut v1, v2, v2_shoup, qi)| {
                qi.mul_shoup_vec(
                    v1.as_slice_mut().unwrap(),
                    v2.as_slice().unwrap(),
                    v2_shoup.as_slice().unwrap(),
                )
            });
        }
        self.has_lazy_coefficients = false;
    }
}

impl Mul<&Poly<Ntt>> for &Poly<Ntt> {
    type Output = Poly<Ntt>;
    fn mul(self, p: &Poly<Ntt>) -> Poly<Ntt> {
        let mut q = self.clone();
        q *= p;
        q
    }
}

impl Mul<&Poly<NttShoup>> for &Poly<Ntt> {
    type Output = Poly<Ntt>;
    fn mul(self, p: &Poly<NttShoup>) -> Poly<Ntt> {
        let mut q = self.clone();
        q *= p;
        q
    }
}

impl Mul<&BigUint> for &Poly<Ntt> {
    type Output = Poly<Ntt>;
    fn mul(self, p: &BigUint) -> Poly<Ntt> {
        let mut q = self.clone();
        q *= p;
        q
    }
}

impl Mul<&BigUint> for &Poly<PowerBasis> {
    type Output = Poly<PowerBasis>;
    fn mul(self, p: &BigUint) -> Poly<PowerBasis> {
        let mut q = self.clone();
        q *= p;
        q
    }
}

impl Mul<&Poly<Ntt>> for &BigUint {
    type Output = Poly<Ntt>;
    fn mul(self, p: &Poly<Ntt>) -> Poly<Ntt> {
        p * self
    }
}

impl Mul<&Poly<PowerBasis>> for &BigUint {
    type Output = Poly<PowerBasis>;
    fn mul(self, p: &Poly<PowerBasis>) -> Poly<PowerBasis> {
        p * self
    }
}

#[cfg(test)]
mod tests {
    use itertools::{Itertools, izip};
    use num_bigint::BigUint;
    use rand::rng;

    use crate::rq::dot_product;
    use crate::{
        rq::{Context, Ntt, NttShoup, Poly, PowerBasis},
        zq::Modulus,
    };
    use std::{error::Error, sync::Arc};

    static MODULI: &[u64; 3] = &[1153, 4611686018326724609, 4611686018309947393];

    #[test]
    fn shared_arithmetic_preserves_values_and_timing_policy() {
        let ctx = crate::rq::Context::new_arc(&[1153, 2017], 16).unwrap();
        macro_rules! check {
            ($repr:ty) => {
                for left_public in [false, true] {
                    for right_public in [false, true] {
                        let mut left = crate::rq::Poly::<$repr>::random_from_seed(&ctx, [1; 32]);
                        let mut right = crate::rq::Poly::<$repr>::random_from_seed(&ctx, [2; 32]);
                        left.allow_variable_time_computations = left_public;
                        right.allow_variable_time_computations = right_public;
                        let sum = &left + &right;
                        let difference = &left - &right;
                        let negated = -&left;
                        let mut scaled = left.clone();
                        scaled *= &BigUint::from(7u64);
                        assert_eq!(
                            sum.allows_variable_time_computations(),
                            left_public && right_public
                        );
                        assert_eq!(
                            difference.allows_variable_time_computations(),
                            left_public && right_public
                        );
                        assert_eq!(negated.allows_variable_time_computations(), left_public);
                        assert_eq!(scaled.allows_variable_time_computations(), left_public);
                        assert_eq!(left.clone() + right.clone(), sum);
                        assert_eq!(-left.clone(), negated);
                        for (row, modulus) in ctx.q.iter().enumerate() {
                            for column in 0..ctx.degree {
                                let index = [row, column];
                                let a = left.coefficients[index];
                                let b = right.coefficients[index];
                                assert_eq!(sum.coefficients[index], modulus.add(a, b));
                                assert_eq!(difference.coefficients[index], modulus.sub(a, b));
                                assert_eq!(negated.coefficients[index], modulus.neg(a));
                                assert_eq!(scaled.coefficients[index], modulus.mul(a, 7));
                            }
                        }
                        left.has_lazy_coefficients = true;
                        assert!(std::panic::catch_unwind(|| &left + &right).is_err());
                        assert!(std::panic::catch_unwind(|| &left - &right).is_err());
                        assert!(std::panic::catch_unwind(|| -&left).is_err());
                    }
                }
            };
        }
        check!(crate::rq::PowerBasis);
        check!(crate::rq::Ntt);
    }

    #[test]
    fn add() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let n = 16;
        for _ in 0..100 {
            for modulus in MODULI {
                let ctx = Arc::new(Context::new(&[*modulus], n)?);
                let m = Modulus::new(*modulus).unwrap();

                let p = Poly::<PowerBasis>::random(&ctx, &mut rng);
                let q = Poly::<PowerBasis>::random(&ctx, &mut rng);
                let r = &p + &q;
                let mut a = Vec::<u64>::try_from(&p).unwrap();
                m.add_vec(&mut a, &Vec::<u64>::try_from(&q).unwrap());
                assert_eq!(Vec::<u64>::try_from(&r).unwrap(), a);

                let p = Poly::<Ntt>::random(&ctx, &mut rng);
                let q = Poly::<Ntt>::random(&ctx, &mut rng);
                let r = &p + &q;
                let mut a = Vec::<u64>::try_from(&p).unwrap();
                m.add_vec(&mut a, &Vec::<u64>::try_from(&q).unwrap());
                assert_eq!(Vec::<u64>::try_from(&r).unwrap(), a);
            }

            let ctx = Arc::new(Context::new(MODULI, 16)?);
            let p = Poly::<PowerBasis>::random(&ctx, &mut rng);
            let q = Poly::<PowerBasis>::random(&ctx, &mut rng);
            let mut a = Vec::<u64>::try_from(&p).unwrap();
            let b = Vec::<u64>::try_from(&q).unwrap();
            for i in 0..MODULI.len() {
                let m = Modulus::new(MODULI[i]).unwrap();
                m.add_vec(&mut a[i * 16..(i + 1) * 16], &b[i * 16..(i + 1) * 16])
            }
            let r = &p + &q;
            assert_eq!(Vec::<u64>::try_from(&r).unwrap(), a);
        }
        Ok(())
    }

    #[test]
    fn sub() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for _ in 0..100 {
            for modulus in MODULI {
                let ctx = Arc::new(Context::new(&[*modulus], 16)?);
                let m = Modulus::new(*modulus).unwrap();

                let p = Poly::<PowerBasis>::random(&ctx, &mut rng);
                let q = Poly::<PowerBasis>::random(&ctx, &mut rng);
                let r = &p - &q;
                let mut a = Vec::<u64>::try_from(&p).unwrap();
                m.sub_vec(&mut a, &Vec::<u64>::try_from(&q).unwrap());
                assert_eq!(Vec::<u64>::try_from(&r).unwrap(), a);

                let p = Poly::<Ntt>::random(&ctx, &mut rng);
                let q = Poly::<Ntt>::random(&ctx, &mut rng);
                let r = &p - &q;
                let mut a = Vec::<u64>::try_from(&p).unwrap();
                m.sub_vec(&mut a, &Vec::<u64>::try_from(&q).unwrap());
                assert_eq!(Vec::<u64>::try_from(&r).unwrap(), a);
            }

            let ctx = Arc::new(Context::new(MODULI, 16)?);
            let p = Poly::<PowerBasis>::random(&ctx, &mut rng);
            let q = Poly::<PowerBasis>::random(&ctx, &mut rng);
            let mut a = Vec::<u64>::try_from(&p).unwrap();
            let b = Vec::<u64>::try_from(&q).unwrap();
            for i in 0..MODULI.len() {
                let m = Modulus::new(MODULI[i]).unwrap();
                m.sub_vec(&mut a[i * 16..(i + 1) * 16], &b[i * 16..(i + 1) * 16])
            }
            let r = &p - &q;
            assert_eq!(Vec::<u64>::try_from(&r).unwrap(), a);
        }
        Ok(())
    }

    #[test]
    fn mul() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for _ in 0..100 {
            for modulus in MODULI {
                let ctx = Arc::new(Context::new(&[*modulus], 16)?);
                let m = Modulus::new(*modulus).unwrap();

                let p = Poly::<Ntt>::random(&ctx, &mut rng);
                let q = Poly::<Ntt>::random(&ctx, &mut rng);
                let r = &p * &q;
                let mut a = Vec::<u64>::try_from(&p).unwrap();
                m.mul_vec(&mut a, &Vec::<u64>::try_from(&q).unwrap());
                assert_eq!(Vec::<u64>::try_from(&r).unwrap(), a);
            }

            let ctx = Arc::new(Context::new(MODULI, 16)?);
            let p = Poly::<Ntt>::random(&ctx, &mut rng);
            let q = Poly::<Ntt>::random(&ctx, &mut rng);
            let mut a = Vec::<u64>::try_from(&p).unwrap();
            let b = Vec::<u64>::try_from(&q).unwrap();
            for i in 0..MODULI.len() {
                let m = Modulus::new(MODULI[i]).unwrap();
                m.mul_vec(&mut a[i * 16..(i + 1) * 16], &b[i * 16..(i + 1) * 16])
            }
            let r = &p * &q;
            assert_eq!(Vec::<u64>::try_from(&r).unwrap(), a);
        }
        Ok(())
    }

    #[test]
    fn mul_shoup() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for _ in 0..100 {
            for modulus in MODULI {
                let ctx = Arc::new(Context::new(&[*modulus], 16)?);
                let m = Modulus::new(*modulus).unwrap();

                let p = Poly::<Ntt>::random(&ctx, &mut rng);
                let q = Poly::<NttShoup>::random(&ctx, &mut rng);
                let r = &p * &q;
                let mut a = Vec::<u64>::try_from(&p).unwrap();
                m.mul_vec(&mut a, &Vec::<u64>::try_from(&q).unwrap());
                assert_eq!(Vec::<u64>::try_from(&r).unwrap(), a);
            }

            let ctx = Arc::new(Context::new(MODULI, 16)?);
            let p = Poly::<Ntt>::random(&ctx, &mut rng);
            let q = Poly::<NttShoup>::random(&ctx, &mut rng);
            let mut a = Vec::<u64>::try_from(&p).unwrap();
            let b = Vec::<u64>::try_from(&q).unwrap();
            for i in 0..MODULI.len() {
                let m = Modulus::new(MODULI[i]).unwrap();
                m.mul_vec(&mut a[i * 16..(i + 1) * 16], &b[i * 16..(i + 1) * 16])
            }
            let r = &p * &q;
            assert_eq!(Vec::<u64>::try_from(&r).unwrap(), a);
        }
        Ok(())
    }

    #[test]
    fn neg() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for _ in 0..100 {
            for modulus in MODULI {
                let ctx = Arc::new(Context::new(&[*modulus], 16)?);
                let m = Modulus::new(*modulus).unwrap();

                let p = Poly::<PowerBasis>::random(&ctx, &mut rng);
                let r = -&p;
                let mut a = Vec::<u64>::try_from(&p).unwrap();
                m.neg_vec(&mut a);
                assert_eq!(Vec::<u64>::try_from(&r).unwrap(), a);

                let p = Poly::<Ntt>::random(&ctx, &mut rng);
                let r = -&p;
                let mut a = Vec::<u64>::try_from(&p).unwrap();
                m.neg_vec(&mut a);
                assert_eq!(Vec::<u64>::try_from(&r).unwrap(), a);
            }

            let ctx = Arc::new(Context::new(MODULI, 16)?);
            let p = Poly::<PowerBasis>::random(&ctx, &mut rng);
            let mut a = Vec::<u64>::try_from(&p).unwrap();
            for i in 0..MODULI.len() {
                let m = Modulus::new(MODULI[i]).unwrap();
                m.neg_vec(&mut a[i * 16..(i + 1) * 16])
            }
            let r = -&p;
            assert_eq!(Vec::<u64>::try_from(&r).unwrap(), a);

            let r = -p;
            assert_eq!(Vec::<u64>::try_from(&r).unwrap(), a);
        }
        Ok(())
    }

    #[test]
    fn test_dot_product() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for _ in 0..20 {
            for modulus in MODULI {
                let ctx = Arc::new(Context::new(&[*modulus], 16)?);

                for len in 1..50 {
                    let p = (0..len)
                        .map(|_| Poly::<Ntt>::random(&ctx, &mut rng))
                        .collect_vec();
                    let q = (0..len)
                        .map(|_| Poly::<Ntt>::random(&ctx, &mut rng))
                        .collect_vec();
                    let r = dot_product(p.iter(), q.iter())?;

                    let mut expected = Poly::<Ntt>::zero(&ctx);
                    izip!(&p, &q).for_each(|(pi, qi)| expected += &(pi * qi));
                    assert_eq!(r, expected);
                }
            }

            let ctx = Arc::new(Context::new(MODULI, 16)?);
            for len in 1..50 {
                let p = (0..len)
                    .map(|_| Poly::<Ntt>::random(&ctx, &mut rng))
                    .collect_vec();
                let q = (0..len)
                    .map(|_| Poly::<Ntt>::random(&ctx, &mut rng))
                    .collect_vec();
                let r = dot_product(p.iter(), q.iter())?;

                let mut expected = Poly::<Ntt>::zero(&ctx);
                izip!(&p, &q).for_each(|(pi, qi)| expected += &(pi * qi));
                assert_eq!(r, expected);
            }
        }
        Ok(())
    }

    #[test]
    fn dot_product_requires_all_operands_to_allow_variable_time() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let ctx = Arc::new(Context::new(&MODULI[..1], 16)?);
        let variable_time = fhe_util::VariableTime::new(fhe_util::PublicData::assert_public());
        let mut p = (0..2)
            .map(|_| Poly::<Ntt>::random(&ctx, &mut rng))
            .collect_vec();
        let mut q = (0..2)
            .map(|_| Poly::<Ntt>::random(&ctx, &mut rng))
            .collect_vec();

        p.iter_mut()
            .chain(q.iter_mut())
            .for_each(|poly| poly.allow_variable_time_computations(variable_time));
        let all_public = dot_product(p.iter(), q.iter())?;
        assert!(all_public.allows_variable_time_computations());

        q[1].disallow_variable_time_computations();
        let mixed = dot_product(p.iter(), q.iter())?;
        assert!(!mixed.allows_variable_time_computations());
        Ok(())
    }

    #[test]
    fn dot_product_rejects_mismatched_inputs() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let ctx = Arc::new(Context::new(&MODULI[..1], 16)?);
        let other_ctx = Arc::new(Context::new(&MODULI[1..2], 16)?);
        let p = [Poly::<Ntt>::random(&ctx, &mut rng)];
        let q = [
            Poly::<Ntt>::random(&ctx, &mut rng),
            Poly::<Ntt>::random(&ctx, &mut rng),
        ];

        assert!(matches!(
            dot_product(p.iter(), q.iter()),
            Err(crate::Error::DotProductLengthMismatch { left: 1, right: 2 })
        ));

        let q = [Poly::<Ntt>::random(&other_ctx, &mut rng)];
        assert_eq!(
            dot_product(p.iter(), q.iter()),
            Err(crate::Error::PolynomialContextMismatch)
        );
        Ok(())
    }

    #[test]
    fn mul_scalar() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for _ in 0..100 {
            for modulus in MODULI {
                let ctx = Arc::new(Context::new(&[*modulus], 16)?);
                let m = Modulus::new(*modulus).unwrap();

                // Test with PowerBasis representation
                let p = Poly::<PowerBasis>::random(&ctx, &mut rng);
                let scalar = BigUint::from(42u64);
                let r = &p * &scalar;
                let mut expected = Vec::<u64>::try_from(&p).unwrap();
                m.scalar_mul_vec(&mut expected, 42u64);
                assert_eq!(Vec::<u64>::try_from(&r).unwrap(), expected);

                // Test with NTT representation
                let p = Poly::<Ntt>::random(&ctx, &mut rng);
                let scalar = BigUint::from(123u64);
                let r = &p * &scalar;
                let mut expected = Vec::<u64>::try_from(&p).unwrap();
                m.scalar_mul_vec(&mut expected, 123u64);
                assert_eq!(Vec::<u64>::try_from(&r).unwrap(), expected);
            }

            let ctx = Arc::new(Context::new(MODULI, 16)?);

            // Test with PowerBasis representation
            let p = Poly::<PowerBasis>::random(&ctx, &mut rng);
            let scalar = BigUint::from(99u64);
            let r = &p * &scalar;
            let mut expected = Vec::<u64>::try_from(&p).unwrap();
            for i in 0..MODULI.len() {
                let m = Modulus::new(MODULI[i]).unwrap();
                m.scalar_mul_vec(&mut expected[i * 16..(i + 1) * 16], 99u64)
            }
            assert_eq!(Vec::<u64>::try_from(&r).unwrap(), expected);

            // Test with NTT representation
            let p = Poly::<Ntt>::random(&ctx, &mut rng);
            let scalar = BigUint::from(77u64);
            let r = &p * &scalar;
            let mut expected = Vec::<u64>::try_from(&p).unwrap();
            for i in 0..MODULI.len() {
                let m = Modulus::new(MODULI[i]).unwrap();
                m.scalar_mul_vec(&mut expected[i * 16..(i + 1) * 16], 77u64)
            }
            assert_eq!(Vec::<u64>::try_from(&r).unwrap(), expected);
        }
        Ok(())
    }

    #[test]
    fn mul_scalar_large_crt() -> Result<(), Box<dyn Error>> {
        let ctx = Arc::new(Context::new(MODULI, 16)?);

        // Create a large scalar that exceeds the max modulus
        let q_prod = MODULI.iter().fold(BigUint::from(1u64), |acc, &m| acc * m);
        let large_scalar = &q_prod + BigUint::from(12345u64);

        let p = Poly::<Ntt>::random(&ctx, &mut rng());
        let r = &p * &large_scalar;

        // Verify by computing the expected result manually for each modulus
        let mut expected = Vec::<u64>::try_from(&p).unwrap();
        for i in 0..MODULI.len() {
            let m = Modulus::new(MODULI[i]).unwrap();
            // Reduce the large scalar modulo this prime
            let scalar_mod_qi = (&large_scalar % MODULI[i]).to_u64_digits()[0];
            m.scalar_mul_vec(&mut expected[i * 16..(i + 1) * 16], scalar_mod_qi)
        }
        assert_eq!(Vec::<u64>::try_from(&r).unwrap(), expected);

        Ok(())
    }

    #[test]
    fn mul_scalar_ntt_shoup() {
        let ctx = Arc::new(Context::new(MODULI, 16).unwrap());
        let p = Poly::<NttShoup>::random(&ctx, &mut rng());
        let mut p_ntt = p.clone().into_ntt();
        let scalar = BigUint::from(42u64);

        let mut p_ntt_scaled = p_ntt.clone();
        p_ntt_scaled *= &scalar;

        p_ntt *= &scalar;
        assert_eq!(p_ntt_scaled, p_ntt);
    }
}
