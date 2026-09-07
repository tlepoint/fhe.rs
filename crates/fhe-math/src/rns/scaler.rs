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
///
/// Construct this value explicitly; there is no uninitialized default.
/// ```compile_fail
/// let invalid = fhe_math::rns::ScalingFactor::default();
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
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
///
/// Fractional corrections use fixed-point approximations. For large contexts,
/// results extremely close to centering or rounding boundaries can differ from
/// exact rational arithmetic.
///
/// Construct this value explicitly; there is no uninitialized default.
/// ```compile_fail
/// let invalid = fhe_math::rns::RnsScaler::default();
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
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
    theta_omega: Box<[RoundingTerm]>,

    theta_garner_lo: Box<[u64]>,
    theta_garner_hi: Box<[u64]>,
    theta_garner_shift: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct RoundingTerm {
    index: usize,
    lo: u64,
    hi: u64,
    negative: bool,
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

        // Integral projections need no rounding correction. In particular,
        // when BFV scales from Q*P by t/Q, the extra P residues have integral
        // Garner projections. Their zero terms can be omitted for every
        // coefficient. This schedule depends only on the public contexts and
        // scaling factor, never on the residues being scaled.
        let theta_omega = izip!(theta_omega_lo, theta_omega_hi, theta_omega_sign)
            .enumerate()
            .filter_map(|(index, (lo, hi, negative))| {
                ((lo | hi) != 0).then_some(RoundingTerm {
                    index,
                    lo,
                    hi,
                    negative,
                })
            })
            .collect();

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
            theta_omega,
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
    /// Panics if the input length differs from the source modulus count, or
    /// `size` is zero or exceeds the destination modulus count.
    #[must_use]
    pub fn scale_new(&self, rests: ArrayView1<u64>, size: usize) -> Vec<u64> {
        assert_eq!(rests.len(), self.from.moduli_u64.len());
        assert!(size > 0 && size <= self.to.moduli_u64.len());
        let mut out = vec![0; size];
        self.scale(rests, (&mut out).into(), 0);
        out
    }

    /// Compute the RNS representation of the rests scaled by numerator *
    /// denominator, and either rounded or floored, and store the result in
    /// `out`.
    ///
    /// Panics if the input length differs from the source modulus count, or
    /// the nonempty output range lies outside the destination moduli.
    pub fn scale(&self, rests: ArrayView1<u64>, out: ArrayViewMut1<u64>, starting_index: usize) {
        assert_eq!(rests.len(), self.from.moduli_u64.len());
        assert!(!out.is_empty());
        assert!(starting_index <= self.to.moduli_u64.len());
        assert!(out.len() <= self.to.moduli_u64.len() - starting_index);

        // Specialize once on the public scaling factor, so basis conversion
        // does not carry the fractional-correction branches and scratch.
        if self.scaling_factor.is_one {
            self.scale_inner::<true>(rests, out, starting_index);
        } else {
            self.scale_inner::<false>(rests, out, starting_index);
        }
    }

    fn scale_inner<const IDENTITY: bool>(
        &self,
        rests: ArrayView1<u64>,
        mut out: ArrayViewMut1<u64>,
        starting_index: usize,
    ) {
        // First, let's compute the inner product of the rests with theta_garner.
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
        let v = sum_theta_garner.as_u128();
        let v = (v >> 1) + (v & 1);

        // If the scaling factor is not 1, compute the inner product with the
        // theta_omega
        let mut w_sign = 0u64;
        let mut w = 0u128;
        if !IDENTITY {
            let mut sum_theta_omega = u256::ZERO;
            for term in &self.theta_omega {
                let product = U256::from(rests[term.index])
                    * U256::from((term.lo as u128) | ((term.hi as u128) << 64));
                if term.negative {
                    sum_theta_omega = sum_theta_omega.wrapping_sub(product);
                } else {
                    sum_theta_omega = sum_theta_omega.wrapping_add(product);
                }
            }

            // Let's subtract v * theta_gamma to sum_theta_omega.
            // This correction also vanishes exactly for BFV's t/Q scaling
            // from Q*P, since t*P is integral.
            if (self.theta_gamma_lo | self.theta_gamma_hi) != 0 {
                let v_theta_gamma = U256::from(v)
                    * U256::from(
                        (self.theta_gamma_lo as u128) | ((self.theta_gamma_hi as u128) << 64),
                    );
                if self.theta_gamma_sign {
                    sum_theta_omega = sum_theta_omega.wrapping_add(v_theta_gamma);
                } else {
                    sum_theta_omega = sum_theta_omega.wrapping_sub(v_theta_gamma);
                }
            }

            // Let's compute w = round(sum_theta_omega / 2^(192)).
            let sign_bits = (sum_theta_omega >> 191isize).as_u128();
            // Branch-free nonzero test of the upper bits, preserving the
            // existing signed rounding convention.
            w_sign = ((sign_bits | sign_bits.wrapping_neg()) >> 127) as u64;

            // Both rounding candidates are computed without branching on the
            // secret-dependent sign. Wrapping handles the unused candidate.
            let positive = (sum_theta_omega >> 126isize).as_u128();
            let positive = (positive >> 1) + (positive & 1);
            let negative = ((!sum_theta_omega) >> 126isize).as_u128().wrapping_add(1) >> 1;
            let mask = 0u128.wrapping_sub(w_sign as u128);
            w = positive ^ ((positive ^ negative) & mask);
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

                if !IDENTITY {
                    let wi = qi.lazy_reduce_u128(w);
                    let mask = 0u64.wrapping_sub(w_sign);
                    let signed_wi = wi ^ ((wi ^ (**qi * 2 - wi)) & mask);
                    yi += signed_wi as u128;
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
    use rand::{Rng as RngCore, rng};

    #[test]
    fn rejects_invalid_dimensions_in_all_profiles() -> Result<(), Box<dyn Error>> {
        let ctx = Arc::new(RnsContext::new(&[17])?);
        let scaler = RnsScaler::new(&ctx, &ctx, ScalingFactor::one());
        for rests in [vec![], vec![1, 2]] {
            assert!(catch_unwind(|| scaler.scale_new(rests.as_slice().into(), 1)).is_err());
        }
        for size in [0, 2, usize::MAX] {
            assert!(catch_unwind(|| scaler.scale_new((&[1u64][..]).into(), size)).is_err());
        }
        for start in [1, 2, usize::MAX] {
            assert!(
                catch_unwind(|| {
                    let mut out = [0u64];
                    scaler.scale((&[1u64][..]).into(), (&mut out[..]).into(), start);
                })
                .is_err()
            );
        }
        Ok(())
    }

    #[test]
    fn signed_rounding_matches_integer_reference() -> Result<(), Box<dyn Error>> {
        // Exhaust both signs and rounding ties, with odd and even source products.
        for moduli in [&[17u64, 19][..], &[4u64, 17][..]] {
            let from = Arc::new(RnsContext::new(moduli)?);
            let to = Arc::new(RnsContext::new(&[37, 41])?);
            let modulus = from.modulus().to_i64().unwrap();
            for (num, den) in [(1i64, 2i64), (2, 3), (5, 7), (17, modulus)] {
                let scaler = RnsScaler::new(
                    &from,
                    &to,
                    ScalingFactor::new(&BigUint::from(num as u64), &BigUint::from(den as u64)),
                );
                for x in 0..modulus {
                    let centered = if x >= (modulus + 1) / 2 {
                        x - modulus
                    } else {
                        x
                    };
                    let product = centered * num;
                    let rounded = if product < 0 {
                        -((-product + (den - 1) / 2) / den)
                    } else {
                        (product + den / 2) / den
                    };
                    let rests = from.project(&BigUint::from(x as u64));
                    let actual = scaler.scale_new(rests.as_slice().into(), 2);
                    assert_eq!(
                        actual,
                        vec![rounded.rem_euclid(37) as u64, rounded.rem_euclid(41) as u64],
                        "x={x}, factor={num}/{den}"
                    );
                    let mut suffix = [0u64];
                    scaler.scale(rests.as_slice().into(), (&mut suffix[..]).into(), 1);
                    assert_eq!(suffix[0], actual[1]);
                }
            }
        }
        Ok(())
    }

    fn exact_scaled_residues(
        x: &BigUint,
        from: &RnsContext,
        to: &RnsContext,
        numerator: &BigUint,
        denominator: &BigUint,
    ) -> Vec<u64> {
        let negative = x >= &((from.modulus() + 1u64) >> 1usize);
        let magnitude = if negative {
            from.modulus() - x
        } else {
            x.clone()
        };
        // Signed ties round toward positive infinity.
        let bias = if negative {
            (denominator - 1u64) >> 1usize
        } else {
            denominator >> 1usize
        };
        let rounded = ((magnitude * numerator + bias) / denominator) % to.modulus();
        let reduced = if negative && !rounded.is_zero() {
            to.modulus() - rounded
        } else {
            rounded
        };
        to.project(&reduced)
    }

    #[test]
    fn sparse_rounding_matches_exact_arithmetic_with_strided_views() -> Result<(), Box<dyn Error>> {
        use ndarray::s;

        let to = Arc::new(RnsContext::new(&[37, 41, 43])?);
        for (moduli, denominator, indices) in [
            ([5u64, 7, 11], 35u64, [0, 1]),
            ([11, 5, 7], 35, [1, 2]),
            ([4, 9, 5], 36, [0, 1]),
        ] {
            let from = Arc::new(RnsContext::new(&moduli)?);
            for (num, den) in [(0u64, 1u64), (2, 1), (2, 2), (3, denominator), (7, 13)] {
                let numerator = BigUint::from(num);
                let denominator = BigUint::from(den);
                let scaler =
                    RnsScaler::new(&from, &to, ScalingFactor::new(&numerator, &denominator));
                if num == 3 {
                    assert_eq!(
                        scaler
                            .theta_omega
                            .iter()
                            .map(|term| term.index)
                            .collect::<Vec<_>>(),
                        indices
                    );
                    assert_eq!(scaler.theta_gamma_lo | scaler.theta_gamma_hi, 0);
                } else if num == 7 {
                    assert_eq!(scaler.theta_omega.len(), moduli.len());
                    assert_ne!(scaler.theta_gamma_lo | scaler.theta_gamma_hi, 0);
                } else {
                    assert!(scaler.theta_omega.is_empty());
                }
                for x in 0..from.modulus().to_u64().unwrap() {
                    let x = BigUint::from(x);
                    let expected = exact_scaled_residues(&x, &from, &to, &numerator, &denominator);
                    let rests = from.project(&x);
                    assert_eq!(scaler.scale_new(rests.as_slice().into(), 3), expected);

                    let interleaved: Vec<_> = rests.iter().flat_map(|r| [*r, u64::MAX]).collect();
                    let input = ArrayView1::from(&interleaved);
                    let mut output = [u64::MAX; 4];
                    scaler.scale(
                        input.slice(s![..;2]),
                        ndarray::ArrayViewMut1::from(&mut output[..]).slice_move(s![1..;2]),
                        1,
                    );
                    assert_eq!(output, [u64::MAX, expected[1], u64::MAX, expected[2]]);
                }
            }
        }
        Ok(())
    }

    #[test]
    fn bfv_sparse_rounding_preserves_dense_boundaries_and_exact_random_results()
    -> Result<(), Box<dyn Error>> {
        use rand::SeedableRng;
        use rand_chacha::ChaCha8Rng;

        let moduli = [
            562949954093057u64,
            4611686018326724609,
            4611686018309947393,
            4611686018282684417,
            4611686018257518593,
        ];
        let from = Arc::new(RnsContext::new(&moduli)?);
        let to = Arc::new(RnsContext::new(&moduli[..2])?);
        let denominator = to.modulus();
        let mut rng = ChaCha8Rng::seed_from_u64(42);
        let mut inputs = vec![
            BigUint::from(0u64),
            BigUint::from(1u64),
            from.modulus() - 1u64,
        ];
        let mut boundaries = Vec::new();
        for center in [from.modulus() >> 1usize, denominator >> 1usize] {
            for x in [&center - 1u64, center.clone(), &center + 1u64] {
                boundaries.push(from.modulus() - &x);
                boundaries.push(x);
            }
        }
        for _ in 0..128 {
            let residues: Vec<_> = moduli.iter().map(|q| rng.next_u64() % q).collect();
            inputs.push(from.lift(residues.as_slice().into()));
        }
        for num in [1u64, 2, 2056193] {
            let numerator = BigUint::from(num);
            let scaler = RnsScaler::new(&from, &to, ScalingFactor::new(&numerator, denominator));
            assert_eq!(scaler.theta_omega.len(), 2);
            assert_eq!(scaler.theta_gamma_lo | scaler.theta_gamma_hi, 0);
            // Reconstruct the original dense schedule, including zero terms.
            // For large contexts the existing fixed-point approximation can
            // disagree with exact arithmetic extremely close to centering or
            // rounding ties. Preserve its behavior at those boundaries; this
            // optimization changes only the work performed, not the precision.
            let mut dense = scaler.clone();
            dense.theta_omega = from
                .garner
                .iter()
                .enumerate()
                .map(|(index, garner)| {
                    let (_, lo, hi, negative) = RnsScaler::extract_projection_and_theta(
                        &to,
                        garner,
                        &numerator,
                        denominator,
                        true,
                    );
                    super::RoundingTerm {
                        index,
                        lo,
                        hi,
                        negative,
                    }
                })
                .collect();
            for x in inputs.iter().chain(&boundaries) {
                let residues = from.project(x);
                assert_eq!(
                    scaler.scale_new(residues.as_slice().into(), 2),
                    dense.scale_new(residues.as_slice().into(), 2)
                );
            }
            for x in &inputs {
                let residues = from.project(x);
                assert_eq!(
                    scaler.scale_new(residues.as_slice().into(), 2),
                    exact_scaled_residues(x, &from, &to, &numerator, denominator),
                    "num={num}, x={x}"
                );
            }
        }
        Ok(())
    }

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
