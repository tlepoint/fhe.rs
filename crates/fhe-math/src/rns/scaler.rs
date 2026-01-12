#![warn(missing_docs, unused_imports)]
// Expect indexing in this performance-critical RNS scaler implementation
#![expect(
    clippy::indexing_slicing,
    reason = "performance or example code relies on validated indices"
)]

//! RNS scaler inspired from Remark 3.2 of <https://eprint.iacr.org/2021/204.pdf>.

use super::RnsContext;
use ethnum::{U256, u256};
use itertools::{Itertools, izip};
use ndarray::{ArrayView1, ArrayViewMut1};
use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};
use std::{cmp::min, sync::Arc};

/// Scaling factor when performing a RNS scaling.
#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub struct ScalingFactor {
    numerator: BigUint,
    denominator: BigUint,
    pub(crate) is_one: bool,
}

impl ScalingFactor {
    /// Create a new scaling factor. Aborts if the denominator is 0.
    #[must_use]
    pub fn new(numerator: &BigUint, denominator: &BigUint) -> Self {
        assert_ne!(denominator, &BigUint::zero());
        Self {
            numerator: numerator.clone(),
            denominator: denominator.clone(),
            is_one: numerator == denominator,
        }
    }

    /// Returns the identity element of `Self`.
    #[must_use]
    pub fn one() -> Self {
        Self {
            numerator: BigUint::one(),
            denominator: BigUint::one(),
            is_one: true,
        }
    }
}

/// Scaler for a RNS context.
/// This is a helper struct to perform RNS scaling.
#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub struct RnsScaler {
    from: Arc<RnsContext>,
    to: Arc<RnsContext>,
    scaling_factor: ScalingFactor,

    gamma: Box<[u64]>,
    gamma_shoup: Box<[u64]>,
    theta_gamma_lo: u64,
    theta_gamma_hi: u64,
    theta_gamma_sign: bool,

    omega: Box<[Box<[u64]>]>,
    omega_shoup: Box<[Box<[u64]>]>,
    theta_omega_lo: Box<[u64]>,
    theta_omega_hi: Box<[u64]>,
    theta_omega_sign: Box<[bool]>,

    theta_garner_lo: Box<[u64]>,
    theta_garner_hi: Box<[u64]>,
    theta_garner_shift: usize,
}

impl RnsScaler {
    /// Create a RNS scaler by numerator / denominator.
    ///
    /// Aborts if denominator is equal to 0.
    #[must_use]
    pub fn new(
        from: &Arc<RnsContext>,
        to: &Arc<RnsContext>,
        scaling_factor: ScalingFactor,
    ) -> Self {
        // Let's define gamma = round(numerator * from.product / denominator)
        let (gamma, theta_gamma_lo, theta_gamma_hi, theta_gamma_sign) =
            Self::extract_projection_and_theta(
                to,
                &from.product,
                &scaling_factor.numerator,
                &scaling_factor.denominator,
                false,
            );
        let gamma_shoup = izip!(&gamma, &to.moduli)
            .map(|(wi, q)| q.shoup(*wi))
            .collect_vec();

        // Let's define omega_i = round(from.garner_i * numerator / denominator)
        let mut omega = vec![vec![0u64; from.moduli.len()].into_boxed_slice(); to.moduli.len()];
        let mut omega_shoup =
            vec![vec![0u64; from.moduli.len()].into_boxed_slice(); to.moduli.len()];
        let (omegas_i, theta_omega_lo, theta_omega_hi, theta_omega_sign): (
            Vec<Vec<u64>>,
            Vec<u64>,
            Vec<u64>,
            Vec<bool>,
        ) = from
            .garner
            .iter()
            .map(|garner_i| {
                Self::extract_projection_and_theta(
                    to,
                    garner_i,
                    &scaling_factor.numerator,
                    &scaling_factor.denominator,
                    true,
                )
            })
            .multiunzip();

        for (i, omega_i) in omegas_i.iter().enumerate() {
            for j in 0..to.moduli.len() {
                let qj = &to.moduli[j];
                omega[j][i] = qj.reduce(omega_i[j]);
                omega_shoup[j][i] = qj.shoup(omega[j][i]);
            }
        }

        // Determine the shift so that the sum of the scaled theta_garner fit on an U192
        // (shift + 1) + log(q * n) <= 192
        let theta_garner_shift = min(
            from.moduli_u64
                .iter()
                .map(|qi| {
                    192 - 1
                        - ((*qi as u128) * (from.moduli_u64.len() as u128))
                            .next_power_of_two()
                            .ilog2()
                })
                .min()
                .unwrap(),
            127,
        );
        // Finally, define theta_garner_i = from.garner_i / product, also scaled by
        // 2^127.
        let (theta_garner_lo, theta_garner_hi): (Vec<u64>, Vec<u64>) = from
            .garner
            .iter()
            .map(|garner_i| {
                let mut theta: BigUint =
                    ((garner_i << theta_garner_shift) + (&from.product >> 1)) / &from.product;
                let theta_hi: BigUint = &theta >> 64;
                theta -= &theta_hi << 64;
                (theta.to_u64().unwrap(), theta_hi.to_u64().unwrap())
            })
            .unzip();

        Self {
            from: from.clone(),
            to: to.clone(),
            scaling_factor,
            gamma: gamma.into_boxed_slice(),
            gamma_shoup: gamma_shoup.into_boxed_slice(),
            theta_gamma_lo,
            theta_gamma_hi,
            theta_gamma_sign,
            omega: omega.into_boxed_slice(),
            omega_shoup: omega_shoup.into_boxed_slice(),
            theta_omega_lo: theta_omega_lo.into_boxed_slice(),
            theta_omega_hi: theta_omega_hi.into_boxed_slice(),
            theta_omega_sign: theta_omega_sign.into_boxed_slice(),
            theta_garner_lo: theta_garner_lo.into_boxed_slice(),
            theta_garner_hi: theta_garner_hi.into_boxed_slice(),
            theta_garner_shift: theta_garner_shift as usize,
        }
    }

    // Let's define gamma = round(numerator * input / denominator)
    // and theta_gamma such that theta_gamma = numerator * input / denominator -
    // gamma. This function projects gamma in the RNS context, and scales
    // theta_gamma by 2**127 and rounds. It outputs the projection of gamma in the
    // RNS context, and theta_lo, theta_hi, theta_sign such that theta_gamma =
    // (-1)**theta_sign * (theta_lo + 2^64 * theta_hi).
    fn extract_projection_and_theta(
        ctx: &RnsContext,
        input: &BigUint,
        numerator: &BigUint,
        denominator: &BigUint,
        round_up: bool,
    ) -> (Vec<u64>, u64, u64, bool) {
        let gamma = (numerator * input + (denominator >> 1)) / denominator;
        let projected = ctx.project(&gamma);

        let mut theta = (numerator * input) % denominator;
        let mut theta_sign = false;
        if denominator > &BigUint::one() {
            // If denominator is odd, flip theta if theta > (denominator >> 1)
            if denominator & BigUint::one() == BigUint::one() {
                if theta > (denominator >> 1) {
                    theta_sign = true;
                    theta = denominator - theta;
                }
            } else {
                // denominator is even, flip if theta >= (denominator >> 1)
                if theta >= (denominator >> 1) {
                    theta_sign = true;
                    theta = denominator - theta;
                }
            }
        }
        // theta = ((theta << 127) + (denominator >> 1)) / denominator;
        // We can now split theta into two u64 words.
        if round_up {
            if theta_sign {
                theta = (theta << 127) / denominator;
            } else {
                theta = ((theta << 127) + denominator - BigUint::one()) / denominator;
            }
        } else if theta_sign {
            theta = ((theta << 127) + denominator - BigUint::one()) / denominator;
        } else {
            theta = (theta << 127) / denominator;
        }
        let theta_hi_biguint: BigUint = &theta >> 64;
        theta -= &theta_hi_biguint << 64;
        let theta_lo = theta.to_u64().unwrap();
        let theta_hi = theta_hi_biguint.to_u64().unwrap();

        (projected, theta_lo, theta_hi, theta_sign)
    }

    /// Output the RNS representation of the rests scaled by numerator *
    /// denominator, and either rounded or floored.
    ///
    /// Aborts if the number of rests is different than the number of moduli in
    /// debug mode, or if the size is not in [1, ..., rests.len()].
    #[must_use]
    pub fn scale_new(&self, rests: ArrayView1<u64>, size: usize) -> Vec<u64> {
        let mut out = vec![0; size];
        self.scale(rests, (&mut out).into(), 0);
        out
    }

    /// Compute the RNS representation of the rests scaled by numerator *
    /// denominator, and either rounded or floored, and store the result in
    /// `out`.
    ///
    /// Aborts if the number of rests is different than the number of moduli in
    /// debug mode, or if the size of out is not in [1, ..., rests.len()].
    pub fn scale(
        &self,
        rests: ArrayView1<u64>,
        mut out: ArrayViewMut1<u64>,
        starting_index: usize,
    ) {
        debug_assert_eq!(rests.len(), self.from.moduli_u64.len());
        debug_assert!(!out.is_empty());
        debug_assert!(starting_index + out.len() <= self.to.moduli_u64.len());

        // First, let's compute the inner product of the rests with theta_omega.
        let mut sum_theta_garner = u256::ZERO;
        for (thetag_lo, thetag_hi, ri) in izip!(
            self.theta_garner_lo.iter(),
            self.theta_garner_hi.iter(),
            rests
        ) {
            sum_theta_garner = sum_theta_garner.wrapping_add(
                U256::from(*ri) * U256::from((*thetag_lo as u128) | ((*thetag_hi as u128) << 64)),
            );
        }
        // Let's compute v = round(sum_theta_garner / 2^theta_garner_shift)
        sum_theta_garner >>= self.theta_garner_shift - 1;
        let v = sum_theta_garner.as_u128().div_ceil(2);

        // If the scaling factor is not 1, compute the inner product with the
        // theta_omega
        let mut w_sign = false;
        let mut w = 0u128;
        if !self.scaling_factor.is_one {
            let mut sum_theta_omega = u256::ZERO;
            for (thetao_lo, thetao_hi, thetao_sign, ri) in izip!(
                self.theta_omega_lo.iter(),
                self.theta_omega_hi.iter(),
                self.theta_omega_sign.iter(),
                rests
            ) {
                let product = U256::from(*ri)
                    * U256::from((*thetao_lo as u128) | ((*thetao_hi as u128) << 64));
                if *thetao_sign {
                    sum_theta_omega = sum_theta_omega.wrapping_sub(product);
                } else {
                    sum_theta_omega = sum_theta_omega.wrapping_add(product);
                }
            }

            // Let's subtract v * theta_gamma to sum_theta_omega.
            let v_theta_gamma = U256::from(v)
                * U256::from((self.theta_gamma_lo as u128) | ((self.theta_gamma_hi as u128) << 64));
            if self.theta_gamma_sign {
                sum_theta_omega = sum_theta_omega.wrapping_add(v_theta_gamma);
            } else {
                sum_theta_omega = sum_theta_omega.wrapping_sub(v_theta_gamma);
            }

            // Let's compute w = round(sum_theta_omega / 2^(192)).
            w_sign = (sum_theta_omega >> (63 + 128)) > u256::ZERO;

            if w_sign {
                w = ((!sum_theta_omega) >> 126isize).as_u128() + 1;
                w /= 2;
            } else {
                w = (sum_theta_omega >> 126isize).as_u128();
                w = w.div_ceil(2)
            }
        }

        unsafe {
            for i in 0..out.len() {
                debug_assert!(starting_index + i < self.to.moduli.len());
                debug_assert!(starting_index + i < self.omega.len());
                debug_assert!(starting_index + i < self.omega_shoup.len());
                debug_assert!(starting_index + i < self.gamma.len());
                debug_assert!(starting_index + i < self.gamma_shoup.len());
                let out_i = out.get_mut(i).unwrap();
                let qi = self.to.moduli.get_unchecked(starting_index + i);
                let omega_i = self.omega.get_unchecked(starting_index + i);
                let omega_shoup_i = self.omega_shoup.get_unchecked(starting_index + i);
                let gamma_i = self.gamma.get_unchecked(starting_index + i);
                let gamma_shoup_i = self.gamma_shoup.get_unchecked(starting_index + i);

                let mut yi = (**qi * 2
                    - qi.lazy_mul_shoup(qi.reduce_u128(v), *gamma_i, *gamma_shoup_i))
                    as u128;

                if !self.scaling_factor.is_one {
                    let wi = qi.lazy_reduce_u128(w);
                    yi += if w_sign { **qi * 2 - wi } else { wi } as u128;
                }

                debug_assert!(rests.len() <= omega_i.len());
                debug_assert!(rests.len() <= omega_shoup_i.len());
                for j in 0..rests.len() {
                    yi += qi.lazy_mul_shoup(
                        *rests.get(j).unwrap(),
                        *omega_i.get_unchecked(j),
                        *omega_shoup_i.get_unchecked(j),
                    ) as u128;
                }

                *out_i = qi.reduce_u128(yi)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use std::{error::Error, panic::catch_unwind, sync::Arc};

    use super::RnsScaler;
    use crate::rns::{RnsContext, scaler::ScalingFactor};
    use ndarray::ArrayView1;
    use num_bigint::BigUint;
    use num_traits::{ToPrimitive, Zero};
    use rand::{RngCore, rng};

    #[test]
    fn constructor() -> Result<(), Box<dyn Error>> {
        let q = Arc::new(RnsContext::new(&[4, 4611686018326724609, 1153])?);

        let scaler = RnsScaler::new(&q, &q, ScalingFactor::one());
        assert_eq!(scaler.from, q);

        assert!(
            catch_unwind(|| ScalingFactor::new(&BigUint::from(1u64), &BigUint::zero())).is_err()
        );
        Ok(())
    }

    #[test]
    fn scale_same_context() -> Result<(), Box<dyn Error>> {
        let ntests = 1000;
        let q = Arc::new(RnsContext::new(&[4u64, 4611686018326724609, 1153])?);
        let mut rng = rng();

        for numerator in &[1u64, 2, 3, 100, 1000, 4611686018326724610] {
            for denominator in &[1u64, 2, 3, 4, 100, 101, 1000, 1001, 4611686018326724610] {
                let n = BigUint::from(*numerator);
                let d = BigUint::from(*denominator);
                let scaler = RnsScaler::new(&q, &q, ScalingFactor::new(&n, &d));

                for _ in 0..ntests {
                    let x = vec![
                        rng.next_u64() % q.moduli_u64[0],
                        rng.next_u64() % q.moduli_u64[1],
                        rng.next_u64() % q.moduli_u64[2],
                    ];
                    let mut x_lift = q.lift(ArrayView1::from(&x));
                    let x_sign = x_lift >= (q.modulus() >> 1);
                    if x_sign {
                        x_lift = q.modulus() - x_lift;
                    }

                    let z = scaler.scale_new((&x).into(), x.len());
                    let x_scaled_round = if x_sign {
                        if d.to_u64().unwrap() % 2 == 0 {
                            q.modulus()
                                - (&(&x_lift * &n + ((&d >> 1usize) - 1u64)) / &d) % q.modulus()
                        } else {
                            q.modulus() - (&(&x_lift * &n + (&d >> 1)) / &d) % q.modulus()
                        }
                    } else {
                        &(&x_lift * &n + (&d >> 1)) / &d
                    };
                    assert_eq!(z, q.project(&x_scaled_round));
                }
            }
        }
        Ok(())
    }

    #[test]
    fn scale_different_contexts() -> Result<(), Box<dyn Error>> {
        let ntests = 100;
        let q = Arc::new(RnsContext::new(&[4u64, 4611686018326724609, 1153])?);
        let r = Arc::new(RnsContext::new(&[
            4u64,
            4611686018326724609,
            1153,
            4611686018309947393,
            4611686018282684417,
            4611686018257518593,
            4611686018232352769,
            4611686018171535361,
            4611686018106523649,
            4611686018058289153,
        ])?);
        let mut rng = rng();

        for numerator in &[1u64, 2, 3, 100, 1000, 4611686018326724610] {
            for denominator in &[1u64, 2, 3, 4, 100, 101, 1000, 1001, 4611686018326724610] {
                let n = BigUint::from(*numerator);
                let d = BigUint::from(*denominator);
                let scaler = RnsScaler::new(&q, &r, ScalingFactor::new(&n, &d));
                for _ in 0..ntests {
                    let x = vec![
                        rng.next_u64() % q.moduli_u64[0],
                        rng.next_u64() % q.moduli_u64[1],
                        rng.next_u64() % q.moduli_u64[2],
                    ];

                    let mut x_lift = q.lift(ArrayView1::from(&x));
                    let x_sign = x_lift >= (q.modulus() >> 1);
                    if x_sign {
                        x_lift = q.modulus() - x_lift;
                    }

                    let y = scaler.scale_new((&x).into(), r.moduli.len());
                    let x_scaled_round = if x_sign {
                        if d.to_u64().unwrap() % 2 == 0 {
                            r.modulus()
                                - (&(&x_lift * &n + ((&d >> 1usize) - 1u64)) / &d) % r.modulus()
                        } else {
                            r.modulus() - (&(&x_lift * &n + (&d >> 1)) / &d) % r.modulus()
                        }
                    } else {
                        &(&x_lift * &n + (&d >> 1)) / &d
                    };
                    assert_eq!(y, r.project(&x_scaled_round));
                }
            }
        }
        Ok(())
    }
}
