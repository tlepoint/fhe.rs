//! Modulus switching while retaining the surviving rows in NTT form.

use super::{Context, Ntt, Poly};
use crate::{Error, Result};
use itertools::izip;
use ndarray::Axis;
use std::sync::Arc;
use zeroize::{Zeroize, Zeroizing};

impl Poly<Ntt> {
    /// Divide and round each power-basis coefficient by the last modulus,
    /// dropping that modulus while retaining NTT representation.
    ///
    /// This produces the same polynomial as converting to power basis, calling
    /// [`Poly::<super::PowerBasis>::switch_down`], and converting back. Only
    /// the dropped row needs an inverse NTT; the surviving rows subtract a
    /// forward transform of its rounding correction. Timing permission is
    /// preserved. Returns an error without changing the polynomial if no
    /// next context exists.
    pub fn switch_down(&mut self) -> Result<()> {
        let next_context = self.ctx.next_context.as_ref().ok_or(Error::NoMoreContext)?;
        debug_assert!(!self.has_lazy_coefficients);

        let last_index = self.ctx.q.len() - 1;
        let q_last = &self.ctx.q[last_index];
        let half = **q_last / 2;
        let (mut rows, mut last_row) = self.coefficients.view_mut().split_at(Axis(0), last_index);
        let last = last_row.as_slice_mut().unwrap();
        let last_op = &self.ctx.ops[last_index];
        if self.allow_variable_time_computations {
            unsafe { last_op.backward_vt(last.as_mut_ptr()) };
        } else {
            last_op.backward(last);
        }
        for x in last.iter_mut() {
            *x = q_last.add(*x, half);
        }

        let mut correction = Zeroizing::new(vec![0u64; self.ctx.degree]);
        for (mut row, qi, op, inv, inv_shoup) in izip!(
            rows.outer_iter_mut(),
            self.ctx.q.iter(),
            self.ctx.ops.iter(),
            self.ctx.inv_last_qi_mod_qj.iter(),
            self.ctx.inv_last_qi_mod_qj_shoup.iter(),
        ) {
            // r = ((x mod q_last + floor(q_last/2)) mod q_last)
            //     - floor(q_last/2). Then round(x/q_last) = (x-r)/q_last.
            let half_mod_qi = qi.reduce(half);
            for (r, x) in correction.iter_mut().zip(last.iter()) {
                *r = qi.sub(qi.reduce(*x), half_mod_qi);
            }
            if self.allow_variable_time_computations {
                unsafe { op.forward_vt(correction.as_mut_ptr()) };
            } else {
                op.forward(&mut correction);
            }
            for (x, r) in row.iter_mut().zip(correction.iter()) {
                // Both NTT inputs are canonical, so x + q_i - r < 2*q_i.
                *x = qi.mul_shoup(*x + **qi - *r, *inv, *inv_shoup);
            }
        }

        if !self.allow_variable_time_computations {
            last.zeroize();
        }
        self.coefficients.remove_index(Axis(0), last_index);
        self.ctx = next_context.clone();
        Ok(())
    }

    /// Modulus switch to a descendant context, staying in NTT representation.
    /// An unreachable context returns an error before changing the polynomial.
    pub fn switch_down_to(&mut self, context: &Arc<Context>) -> Result<()> {
        let iterations = self.ctx.niterations_to(context)?;
        if iterations == 1 {
            self.switch_down()?;
        } else if iterations > 1 {
            // A single bulk round trip is cheaper for several drops: applying
            // the NTT correction repeatedly would transform rows that a later
            // step discards.
            let mut power_basis = self.to_power_basis();
            power_basis.switch_down_to(context)?;
            *self = power_basis.into_ntt();
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rq::{PowerBasis, traits::TryConvertFrom};
    use crate::zq::primes::generate_prime;
    use num_bigint::BigUint;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    fn check_chain(p: Poly<PowerBasis>) -> Result<()> {
        for public in [false, true] {
            let mut expected = p.clone();
            expected.allow_variable_time_computations = public;
            let mut actual = expected.clone().into_ntt();
            while let Some(next) = expected.ctx.next_context.as_ref() {
                let last = *expected.ctx.moduli.last().unwrap();
                let integer_reference: Vec<_> = Vec::<BigUint>::from(&expected)
                    .into_iter()
                    .map(|x| ((x + last / 2) / last) % next.modulus())
                    .collect();
                expected.switch_down()?;
                actual.switch_down()?;
                assert_eq!(actual, expected.clone().into_ntt());
                assert_eq!(
                    Vec::<BigUint>::from(&actual.to_power_basis()),
                    integer_reference
                );
                assert_eq!(actual.allows_variable_time_computations(), public);
            }
            let saved = actual.clone();
            assert_eq!(actual.switch_down(), Err(Error::NoMoreContext));
            assert_eq!(actual, saved);
        }
        Ok(())
    }

    #[test]
    fn ntt_switch_down_exhausts_small_coefficients() -> Result<()> {
        let ctx = Context::new_arc(&[97, 193], 16)?;
        for start in (0..97 * 193).step_by(ctx.degree) {
            let coefficients: Vec<_> = (start..start + ctx.degree as u64)
                .map(|x| BigUint::from(x % (97 * 193)))
                .collect();
            check_chain(Poly::try_convert_from(
                coefficients.as_slice(),
                &ctx,
                false,
            )?)?;
        }
        Ok(())
    }

    #[test]
    fn ntt_switch_down_matches_rounding_boundaries_and_random_inputs() -> Result<()> {
        let mut rng = ChaCha8Rng::seed_from_u64(0x517c4);
        for degree in [16, 1024] {
            // Include the PIR widths and a descending chain with large ratios
            // between moduli. Reversing each chain tests both ratio directions.
            for widths in [[50, 55, 55], [36, 36, 37], [62, 40, 20]] {
                let mut moduli = Vec::new();
                for width in widths {
                    let mut bound = 1u64 << width;
                    loop {
                        let prime = generate_prime(width, 2 * degree as u64, bound).unwrap();
                        if !moduli.contains(&prime) {
                            moduli.push(prime);
                            break;
                        }
                        bound = prime;
                    }
                }
                for reverse in [false, true] {
                    if reverse {
                        moduli.reverse();
                    }
                    let ctx = Context::new_arc(&moduli, degree)?;
                    let q = ctx.modulus();
                    let last = *moduli.last().unwrap();
                    let edges = [
                        BigUint::from(0u32),
                        BigUint::from(1u32),
                        q - 1u32,
                        (q >> 1usize) - 1u32,
                        q >> 1usize,
                        (q >> 1usize) + 1u32,
                        BigUint::from(last / 2),
                        BigUint::from(last / 2 + 1),
                        BigUint::from(last - 1),
                        BigUint::from(last),
                        BigUint::from(last + 1),
                        q - last + last / 2,
                        q - last + last / 2 + 1u32,
                    ];
                    let coefficients: Vec<_> = edges.iter().cycle().take(degree).cloned().collect();
                    check_chain(Poly::try_convert_from(
                        coefficients.as_slice(),
                        &ctx,
                        false,
                    )?)?;
                    for _ in 0..3 {
                        check_chain(Poly::random(&ctx, &mut rng))?;
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn ntt_switch_down_to_validates_context_before_mutating() -> Result<()> {
        let ctx = Context::new_arc(&[97, 193, 257, 353], 16)?;
        let mut rng = ChaCha8Rng::seed_from_u64(41);
        let original = Poly::<PowerBasis>::random(&ctx, &mut rng);
        for level in 0..ctx.moduli.len() {
            let target = ctx.context_at_level(level)?;
            let mut expected = original.clone();
            expected.switch_down_to(&target)?;
            let mut actual = original.clone().into_ntt();
            actual.switch_down_to(&target)?;
            assert_eq!(actual, expected.into_ntt());
        }
        for invalid in [
            Context::new_arc(&[97, 257], 16)?,
            Context::new_arc(&[193], 16)?,
            Context::new_arc(&[97], 8)?,
        ] {
            let mut actual = original.clone().into_ntt();
            let saved = actual.clone();
            assert_eq!(
                actual.switch_down_to(&invalid),
                Err(Error::ContextNotReachable)
            );
            assert_eq!(actual, saved);
        }
        Ok(())
    }
}
