//! Implementation of operations over polynomials.

use super::{traits::TryConvertFrom, Poly, Representation};
use crate::{Error, Result};
use itertools::{izip, Itertools};
use ndarray::Array2;
use num_bigint::BigUint;
use std::{
    cmp::min,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};
use zeroize::Zeroize;

impl AddAssign<&Poly> for Poly {
    fn add_assign(&mut self, p: &Poly) {
        assert!(!self.has_lazy_coefficients && !p.has_lazy_coefficients);
        assert_ne!(
            self.representation,
            Representation::NttShoup,
            "Cannot add to a polynomial in NttShoup representation"
        );
        assert_eq!(
            self.representation, p.representation,
            "Incompatible representations"
        );
        debug_assert_eq!(self.ctx, p.ctx, "Incompatible contexts");
        self.allow_variable_time_computations |= p.allow_variable_time_computations;
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

impl Add<&Poly> for &Poly {
    type Output = Poly;
    fn add(self, p: &Poly) -> Poly {
        let mut q = self.clone();
        q += p;
        q
    }
}

impl Add for Poly {
    type Output = Poly;
    fn add(self, mut p: Poly) -> Poly {
        p += &self;
        p
    }
}

impl SubAssign<&Poly> for Poly {
    fn sub_assign(&mut self, p: &Poly) {
        assert!(!self.has_lazy_coefficients && !p.has_lazy_coefficients);
        assert_ne!(
            self.representation,
            Representation::NttShoup,
            "Cannot subtract from a polynomial in NttShoup representation"
        );
        assert_eq!(
            self.representation, p.representation,
            "Incompatible representations"
        );
        debug_assert_eq!(self.ctx, p.ctx, "Incompatible contexts");
        self.allow_variable_time_computations |= p.allow_variable_time_computations;
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

impl Sub<&Poly> for &Poly {
    type Output = Poly;
    fn sub(self, p: &Poly) -> Poly {
        let mut q = self.clone();
        q -= p;
        q
    }
}

impl MulAssign<&Poly> for Poly {
    fn mul_assign(&mut self, p: &Poly) {
        assert!(!p.has_lazy_coefficients);
        assert_ne!(
            self.representation,
            Representation::NttShoup,
            "Cannot multiply to a polynomial in NttShoup representation"
        );
        if self.has_lazy_coefficients && self.representation == Representation::Ntt {
            assert!(
				p.representation == Representation::NttShoup,
				"Can only multiply a polynomial with lazy coefficients by an NttShoup representation."
			);
        } else {
            assert_eq!(
                self.representation,
                Representation::Ntt,
                "Multiplication requires an Ntt representation."
            );
        }
        debug_assert_eq!(self.ctx, p.ctx, "Incompatible contexts");
        self.allow_variable_time_computations |= p.allow_variable_time_computations;

        match p.representation {
            Representation::Ntt => {
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
            Representation::NttShoup => {
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
                self.has_lazy_coefficients = false
            }
            _ => {
                panic!("Multiplication requires a multipliand in Ntt or NttShoup representation.")
            }
        }
    }
}

impl MulAssign<&BigUint> for Poly {
    fn mul_assign(&mut self, p: &BigUint) {
        let v: Vec<BigUint> = vec![p.clone()];
        let mut q = Poly::try_convert_from(
            v.as_ref() as &[BigUint],
            &self.ctx,
            self.allow_variable_time_computations,
            self.representation.clone(),
        )
        .unwrap();
        q.change_representation(Representation::Ntt);
        if self.allow_variable_time_computations {
            unsafe {
                izip!(
                    self.coefficients.outer_iter_mut(),
                    q.coefficients.outer_iter(),
                    self.ctx.q.iter()
                )
                .for_each(|(mut v1, v2, qi)| {
                    qi.mul_vec_vt(v1.as_slice_mut().unwrap(), v2.as_slice().unwrap())
                });
            }
        } else {
            izip!(
                self.coefficients.outer_iter_mut(),
                q.coefficients.outer_iter(),
                self.ctx.q.iter()
            )
            .for_each(|(mut v1, v2, qi)| {
                qi.mul_vec(v1.as_slice_mut().unwrap(), v2.as_slice().unwrap())
            });
        }
    }
}

impl Mul<&Poly> for &Poly {
    type Output = Poly;
    fn mul(self, p: &Poly) -> Poly {
        match self.representation {
            Representation::NttShoup => {
                // TODO: To test, and do the same thing for add, sub, and neg
                let mut q = p.clone();
                if q.representation == Representation::NttShoup {
                    q.coefficients_shoup
                        .as_mut()
                        .unwrap()
                        .as_slice_mut()
                        .unwrap()
                        .zeroize();
                    unsafe { q.override_representation(Representation::Ntt) }
                }
                q *= self;
                q
            }
            _ => {
                let mut q = self.clone();
                q *= p;
                q
            }
        }
    }
}

impl Mul<&BigUint> for &Poly {
    type Output = Poly;
    fn mul(self, p: &BigUint) -> Poly {
        let mut q = self.clone();
        q *= p;
        q
    }
}

impl Mul<&Poly> for &BigUint {
    type Output = Poly;
    fn mul(self, p: &Poly) -> Poly {
        p * self
    }
}

impl Neg for &Poly {
    type Output = Poly;

    fn neg(self) -> Poly {
        assert!(!self.has_lazy_coefficients);
        let mut out = self.clone();
        if self.allow_variable_time_computations {
            izip!(out.coefficients.outer_iter_mut(), out.ctx.q.iter())
                .for_each(|(mut v1, qi)| unsafe { qi.neg_vec_vt(v1.as_slice_mut().unwrap()) });
        } else {
            izip!(out.coefficients.outer_iter_mut(), out.ctx.q.iter())
                .for_each(|(mut v1, qi)| qi.neg_vec(v1.as_slice_mut().unwrap()));
        }
        out
    }
}

impl Neg for Poly {
    type Output = Poly;

    fn neg(mut self) -> Poly {
        assert!(!self.has_lazy_coefficients);
        if self.allow_variable_time_computations {
            izip!(self.coefficients.outer_iter_mut(), self.ctx.q.iter())
                .for_each(|(mut v1, qi)| unsafe { qi.neg_vec_vt(v1.as_slice_mut().unwrap()) });
        } else {
            izip!(self.coefficients.outer_iter_mut(), self.ctx.q.iter())
                .for_each(|(mut v1, qi)| qi.neg_vec(v1.as_slice_mut().unwrap()));
        }
        self
    }
}

/// Computes the Fused-Mul-Add operation `out[i] += x[i] * y[i]`
unsafe fn fma(out: &mut [u128], x: &[u64], y: &[u64]) {
    let n = out.len();
    assert_eq!(x.len(), n);
    assert_eq!(y.len(), n);

    macro_rules! fma_at {
        ($idx:expr) => {
            *out.get_unchecked_mut($idx) +=
                (*x.get_unchecked($idx) as u128) * (*y.get_unchecked($idx) as u128);
        };
    }

    let r = n / 16;
    for i in 0..r {
        fma_at!(16 * i);
        fma_at!(16 * i + 1);
        fma_at!(16 * i + 2);
        fma_at!(16 * i + 3);
        fma_at!(16 * i + 4);
        fma_at!(16 * i + 5);
        fma_at!(16 * i + 6);
        fma_at!(16 * i + 7);
        fma_at!(16 * i + 8);
        fma_at!(16 * i + 9);
        fma_at!(16 * i + 10);
        fma_at!(16 * i + 11);
        fma_at!(16 * i + 12);
        fma_at!(16 * i + 13);
        fma_at!(16 * i + 14);
        fma_at!(16 * i + 15);
    }

    for i in 0..n % 16 {
        fma_at!(16 * r + i);
    }
}

/// Compute the dot product between two iterators of polynomials.
/// Returna an error if the iterator counts are 0, or if any of the polynomial
/// is not in Ntt or NttShoup representation.
pub fn dot_product<'a, 'b, I, J>(p: I, q: J) -> Result<Poly>
where
    I: Iterator<Item = &'a Poly> + Clone,
    J: Iterator<Item = &'b Poly> + Clone,
{
    debug_assert!(!p
        .clone()
        .any(|pi| pi.representation == Representation::PowerBasis));
    debug_assert!(!q
        .clone()
        .any(|qi| qi.representation == Representation::PowerBasis));

    let count = min(p.clone().count(), q.clone().count());
    if count == 0 {
        return Err(Error::Default("At least one iterator is empty".to_string()));
    }

    let p_first = p.clone().next().unwrap();

    // Initialize the accumulator
    let mut acc: Array2<u128> = Array2::zeros((p_first.ctx.q.len(), p_first.ctx.degree));
    let acc_ptr = acc.as_mut_ptr();

    // Current number of products accumulated
    let mut num_acc = vec![1u128; p_first.ctx.q.len()];
    let num_acc_ptr = num_acc.as_mut_ptr();

    // Maximum number of products that can be accumulated
    let max_acc = p_first
        .ctx
        .q
        .iter()
        .map(|qi| 1u128 << (2 * (*qi).leading_zeros()))
        .collect_vec();
    let max_acc_ptr = max_acc.as_ptr();

    let q_ptr = p_first.ctx.q.as_ptr();
    let degree = p_first.ctx.degree as isize;

    let min_of_max = max_acc.iter().min().unwrap();

    let out_slice = acc.as_slice_mut().unwrap();
    if count as u128 > *min_of_max {
        for (pi, qi) in izip!(p, q) {
            let pij = pi.coefficients();
            let qij = qi.coefficients();
            let pi_slice = pij.as_slice().unwrap();
            let qi_slice = qij.as_slice().unwrap();
            unsafe {
                fma(out_slice, pi_slice, qi_slice);

                for j in 0..p_first.ctx.q.len() as isize {
                    let qj = &*q_ptr.offset(j);
                    *num_acc_ptr.offset(j) += 1;
                    if *num_acc_ptr.offset(j) == *max_acc_ptr.offset(j) {
                        if p_first.allow_variable_time_computations {
                            for i in j * degree..(j + 1) * degree {
                                *acc_ptr.offset(i) = qj.reduce_u128_vt(*acc_ptr.offset(i)) as u128;
                            }
                        } else {
                            for i in j * degree..(j + 1) * degree {
                                *acc_ptr.offset(i) = qj.reduce_u128(*acc_ptr.offset(i)) as u128;
                            }
                        }
                        *num_acc_ptr.offset(j) = 1;
                    }
                }
            }
        }
    } else {
        // We don't need to check the condition on the max, it should shave off a few
        // cycles.
        for (pi, qi) in izip!(p, q) {
            let pij = pi.coefficients();
            let qij = qi.coefficients();
            let pi_slice = pij.as_slice().unwrap();
            let qi_slice = qij.as_slice().unwrap();
            unsafe { fma(out_slice, pi_slice, qi_slice) }
        }
    }
    // Last reduction to create the coefficients
    let mut coeffs: Array2<u64> = Array2::zeros((p_first.ctx.q.len(), p_first.ctx.degree));
    izip!(
        coeffs.outer_iter_mut(),
        acc.outer_iter(),
        p_first.ctx.q.iter()
    )
    .for_each(|(mut coeffsj, accj, m)| {
        if p_first.allow_variable_time_computations {
            izip!(coeffsj.iter_mut(), accj.iter())
                .for_each(|(cj, accjk)| *cj = unsafe { m.reduce_u128_vt(*accjk) });
        } else {
            izip!(coeffsj.iter_mut(), accj.iter())
                .for_each(|(cj, accjk)| *cj = m.reduce_u128(*accjk));
        }
    });

    Ok(Poly {
        ctx: p_first.ctx.clone(),
        representation: Representation::Ntt,
        allow_variable_time_computations: p_first.allow_variable_time_computations,
        coefficients: coeffs,
        coefficients_shoup: None,
        has_lazy_coefficients: false,
    })
}

#[cfg(test)]
mod tests {
    use itertools::{izip, Itertools};
    use rand::rng;

    use super::dot_product;
    use crate::{
        rq::{Context, Poly, Representation},
        zq::Modulus,
    };
    use std::{error::Error, sync::Arc};

    static MODULI: &[u64; 3] = &[1153, 4611686018326724609, 4611686018309947393];

    #[test]
    fn add() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let n = 16;
        for _ in 0..100 {
            for modulus in MODULI {
                let ctx = Arc::new(Context::new(&[*modulus], n)?);
                let m = Modulus::new(*modulus).unwrap();

                let p = Poly::random(&ctx, Representation::PowerBasis, &mut rng);
                let q = Poly::random(&ctx, Representation::PowerBasis, &mut rng);
                let r = &p + &q;
                assert_eq!(r.representation, Representation::PowerBasis);
                let mut a = Vec::<u64>::from(&p);
                m.add_vec(&mut a, &Vec::<u64>::from(&q));
                assert_eq!(Vec::<u64>::from(&r), a);

                let p = Poly::random(&ctx, Representation::Ntt, &mut rng);
                let q = Poly::random(&ctx, Representation::Ntt, &mut rng);
                let r = &p + &q;
                assert_eq!(r.representation, Representation::Ntt);
                let mut a = Vec::<u64>::from(&p);
                m.add_vec(&mut a, &Vec::<u64>::from(&q));
                assert_eq!(Vec::<u64>::from(&r), a);
            }

            let ctx = Arc::new(Context::new(MODULI, 16)?);
            let p = Poly::random(&ctx, Representation::PowerBasis, &mut rng);
            let q = Poly::random(&ctx, Representation::PowerBasis, &mut rng);
            let mut a = Vec::<u64>::from(&p);
            let b = Vec::<u64>::from(&q);
            for i in 0..MODULI.len() {
                let m = Modulus::new(MODULI[i]).unwrap();
                m.add_vec(&mut a[i * 16..(i + 1) * 16], &b[i * 16..(i + 1) * 16])
            }
            let r = &p + &q;
            assert_eq!(r.representation, Representation::PowerBasis);
            assert_eq!(Vec::<u64>::from(&r), a);
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

                let p = Poly::random(&ctx, Representation::PowerBasis, &mut rng);
                let q = Poly::random(&ctx, Representation::PowerBasis, &mut rng);
                let r = &p - &q;
                assert_eq!(r.representation, Representation::PowerBasis);
                let mut a = Vec::<u64>::from(&p);
                m.sub_vec(&mut a, &Vec::<u64>::from(&q));
                assert_eq!(Vec::<u64>::from(&r), a);

                let p = Poly::random(&ctx, Representation::Ntt, &mut rng);
                let q = Poly::random(&ctx, Representation::Ntt, &mut rng);
                let r = &p - &q;
                assert_eq!(r.representation, Representation::Ntt);
                let mut a = Vec::<u64>::from(&p);
                m.sub_vec(&mut a, &Vec::<u64>::from(&q));
                assert_eq!(Vec::<u64>::from(&r), a);
            }

            let ctx = Arc::new(Context::new(MODULI, 16)?);
            let p = Poly::random(&ctx, Representation::PowerBasis, &mut rng);
            let q = Poly::random(&ctx, Representation::PowerBasis, &mut rng);
            let mut a = Vec::<u64>::from(&p);
            let b = Vec::<u64>::from(&q);
            for i in 0..MODULI.len() {
                let m = Modulus::new(MODULI[i]).unwrap();
                m.sub_vec(&mut a[i * 16..(i + 1) * 16], &b[i * 16..(i + 1) * 16])
            }
            let r = &p - &q;
            assert_eq!(r.representation, Representation::PowerBasis);
            assert_eq!(Vec::<u64>::from(&r), a);
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

                let p = Poly::random(&ctx, Representation::Ntt, &mut rng);
                let q = Poly::random(&ctx, Representation::Ntt, &mut rng);
                let r = &p * &q;
                assert_eq!(r.representation, Representation::Ntt);
                let mut a = Vec::<u64>::from(&p);
                m.mul_vec(&mut a, &Vec::<u64>::from(&q));
                assert_eq!(Vec::<u64>::from(&r), a);
            }

            let ctx = Arc::new(Context::new(MODULI, 16)?);
            let p = Poly::random(&ctx, Representation::Ntt, &mut rng);
            let q = Poly::random(&ctx, Representation::Ntt, &mut rng);
            let mut a = Vec::<u64>::from(&p);
            let b = Vec::<u64>::from(&q);
            for i in 0..MODULI.len() {
                let m = Modulus::new(MODULI[i]).unwrap();
                m.mul_vec(&mut a[i * 16..(i + 1) * 16], &b[i * 16..(i + 1) * 16])
            }
            let r = &p * &q;
            assert_eq!(r.representation, Representation::Ntt);
            assert_eq!(Vec::<u64>::from(&r), a);
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

                let p = Poly::random(&ctx, Representation::Ntt, &mut rng);
                let q = Poly::random(&ctx, Representation::NttShoup, &mut rng);
                let r = &p * &q;
                assert_eq!(r.representation, Representation::Ntt);
                let mut a = Vec::<u64>::from(&p);
                m.mul_vec(&mut a, &Vec::<u64>::from(&q));
                assert_eq!(Vec::<u64>::from(&r), a);
            }

            let ctx = Arc::new(Context::new(MODULI, 16)?);
            let p = Poly::random(&ctx, Representation::Ntt, &mut rng);
            let q = Poly::random(&ctx, Representation::NttShoup, &mut rng);
            let mut a = Vec::<u64>::from(&p);
            let b = Vec::<u64>::from(&q);
            for i in 0..MODULI.len() {
                let m = Modulus::new(MODULI[i]).unwrap();
                m.mul_vec(&mut a[i * 16..(i + 1) * 16], &b[i * 16..(i + 1) * 16])
            }
            let r = &p * &q;
            assert_eq!(r.representation, Representation::Ntt);
            assert_eq!(Vec::<u64>::from(&r), a);
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

                let p = Poly::random(&ctx, Representation::PowerBasis, &mut rng);
                let r = -&p;
                assert_eq!(r.representation, Representation::PowerBasis);
                let mut a = Vec::<u64>::from(&p);
                m.neg_vec(&mut a);
                assert_eq!(Vec::<u64>::from(&r), a);

                let p = Poly::random(&ctx, Representation::Ntt, &mut rng);
                let r = -&p;
                assert_eq!(r.representation, Representation::Ntt);
                let mut a = Vec::<u64>::from(&p);
                m.neg_vec(&mut a);
                assert_eq!(Vec::<u64>::from(&r), a);
            }

            let ctx = Arc::new(Context::new(MODULI, 16)?);
            let p = Poly::random(&ctx, Representation::PowerBasis, &mut rng);
            let mut a = Vec::<u64>::from(&p);
            for i in 0..MODULI.len() {
                let m = Modulus::new(MODULI[i]).unwrap();
                m.neg_vec(&mut a[i * 16..(i + 1) * 16])
            }
            let r = -&p;
            assert_eq!(r.representation, Representation::PowerBasis);
            assert_eq!(Vec::<u64>::from(&r), a);

            let r = -p;
            assert_eq!(r.representation, Representation::PowerBasis);
            assert_eq!(Vec::<u64>::from(&r), a);
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
                        .map(|_| Poly::random(&ctx, Representation::Ntt, &mut rng))
                        .collect_vec();
                    let q = (0..len)
                        .map(|_| Poly::random(&ctx, Representation::Ntt, &mut rng))
                        .collect_vec();
                    let r = dot_product(p.iter(), q.iter())?;

                    let mut expected = Poly::zero(&ctx, Representation::Ntt);
                    izip!(&p, &q).for_each(|(pi, qi)| expected += &(pi * qi));
                    assert_eq!(r, expected);
                }
            }

            let ctx = Arc::new(Context::new(MODULI, 16)?);
            for len in 1..50 {
                let p = (0..len)
                    .map(|_| Poly::random(&ctx, Representation::Ntt, &mut rng))
                    .collect_vec();
                let q = (0..len)
                    .map(|_| Poly::random(&ctx, Representation::Ntt, &mut rng))
                    .collect_vec();
                let r = dot_product(p.iter(), q.iter())?;

                let mut expected = Poly::zero(&ctx, Representation::Ntt);
                izip!(&p, &q).for_each(|(pi, qi)| expected += &(pi * qi));
                assert_eq!(r, expected);
            }
        }
        Ok(())
    }
}
