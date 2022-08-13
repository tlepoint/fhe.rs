#![warn(missing_docs, unused_imports)]

//! RNS scaler inspired from Remark 3.2 of <https://eprint.iacr.org/2021/204.pdf>.

use super::RnsContext;
use crate::u256::U256;
use itertools::{izip, Itertools};
use ndarray::{ArrayView1, ArrayViewMut1};
use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};
use std::sync::Arc;

/// Scaling factor when performing a RNS scaling.
#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub struct ScalingFactor {
	numerator: BigUint,
	denominator: BigUint,
	pub(crate) is_one: bool,
}

impl ScalingFactor {
	/// Create a new scaling factor. Aborts if the denominator is 0.
	pub fn new(numerator: &BigUint, denominator: &BigUint) -> Self {
		assert_ne!(denominator, &BigUint::zero());
		Self {
			numerator: numerator.clone(),
			denominator: denominator.clone(),
			is_one: numerator == denominator,
		}
	}

	/// Returns the identity element of `Self`.
	pub fn one() -> Self {
		Self {
			numerator: BigUint::one(),
			denominator: BigUint::one(),
			is_one: true,
		}
	}
}

/// Scaler in RNS basis.
#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub struct RnsScaler {
	from: Arc<RnsContext>,
	to: Arc<RnsContext>,
	scaling_factor: ScalingFactor,

	gamma: Vec<u64>,
	gamma_shoup: Vec<u64>,
	theta_gamma_lo: u64,
	theta_gamma_hi: u64,
	theta_gamma_sign: bool,

	omega: Vec<Vec<u64>>,
	omega_shoup: Vec<Vec<u64>>,
	theta_omega_lo: Vec<u64>,
	theta_omega_hi: Vec<u64>,
	theta_omega_sign: Vec<bool>,

	theta_garner_lo: Vec<u64>,
	theta_garner_hi: Vec<u64>,
}

impl RnsScaler {
	/// Create a RNS scaler by numerator / denominator.
	///
	/// Aborts if denominator is equal to 0.
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
		let mut omega = Vec::with_capacity(to.moduli.len());
		let mut omega_shoup = Vec::with_capacity(to.moduli.len());
		for _ in &to.moduli {
			omega.push(vec![0u64; from.moduli.len()]);
			omega_shoup.push(vec![0u64; from.moduli.len()]);
		}
		let mut theta_omega_lo = Vec::with_capacity(from.garner.len());
		let mut theta_omega_hi = Vec::with_capacity(from.garner.len());
		let mut theta_omega_sign = Vec::with_capacity(from.garner.len());
		for i in 0..from.garner.len() {
			let (omega_i, theta_omega_i_lo, theta_omega_i_hi, theta_omega_i_sign) =
				Self::extract_projection_and_theta(
					to,
					&from.garner[i],
					&scaling_factor.numerator,
					&scaling_factor.denominator,
					true,
				);
			for j in 0..to.moduli.len() {
				let qj = &to.moduli[j];
				omega[j][i] = qj.reduce(omega_i[j]);
				omega_shoup[j][i] = qj.shoup(omega[j][i]);
			}
			theta_omega_lo.push(theta_omega_i_lo);
			theta_omega_hi.push(theta_omega_i_hi);
			theta_omega_sign.push(theta_omega_i_sign);
		}

		// Finally, define theta_garner_i = from.garner_i / product, also scaled by 2^127.
		let mut theta_garner_lo = Vec::with_capacity(from.garner.len());
		let mut theta_garner_hi = Vec::with_capacity(from.garner.len());
		for garner_i in &from.garner {
			let mut theta: BigUint = ((garner_i << 127) + (&from.product >> 1)) / &from.product;
			let theta_hi: BigUint = &theta >> 64;
			theta -= &theta_hi << 64;
			theta_garner_lo.push(theta.to_u64().unwrap());
			theta_garner_hi.push(theta_hi.to_u64().unwrap());
		}

		Self {
			from: from.clone(),
			to: to.clone(),
			scaling_factor,
			gamma,
			gamma_shoup,
			theta_gamma_lo,
			theta_gamma_hi,
			theta_gamma_sign,
			omega,
			omega_shoup,
			theta_omega_lo,
			theta_omega_hi,
			theta_omega_sign,
			theta_garner_lo,
			theta_garner_hi,
		}
	}

	// Let's define gamma = round(numerator * input / denominator)
	// and theta_gamma such that theta_gamma = numerator * input / denominator - gamma.
	// This function projects gamma in the RNS context, and scales theta_gamma by 2**127 and rounds.
	// It outputs the projection of gamma in the RNS context,
	// and theta_lo, theta_hi, theta_sign such that theta_gamma = (-1)**theta_sign * (theta_lo + 2^64 * theta_hi).
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

	/// Output the RNS representation of the rests scaled by numerator * denominator,
	/// and either rounded or floored.
	///
	/// Aborts if the number of rests is different than the number of moduli in debug mode, or
	/// if the size is not in [1, ..., rests.len()].
	pub fn scale_new(&self, rests: &ArrayView1<u64>, size: usize, floor: bool) -> Vec<u64> {
		let mut out = vec![0; size];
		self.scale(rests, &mut (&mut out).into(), 0, floor);
		out
	}

	/// Compute the RNS representation of the rests scaled by numerator * denominator,
	/// and either rounded or floored, and store the result in `out`.
	///
	/// Aborts if the number of rests is different than the number of moduli in debug mode, or
	/// if the size of out is not in [1, ..., rests.len()].
	pub fn scale(
		&self,
		rests: &ArrayView1<u64>,
		out: &mut ArrayViewMut1<u64>,
		starting_index: usize,
		floor: bool,
	) {
		debug_assert_eq!(rests.len(), self.from.moduli_u64.len());
		debug_assert!(!out.is_empty() && starting_index + out.len() <= self.to.moduli_u64.len());

		// First, let's compute the inner product of the rests with theta_omega.
		let mut sum_theta_garner = U256::zero();
		for (thetag_lo, thetag_hi, ri) in izip!(&self.theta_garner_lo, &self.theta_garner_hi, rests)
		{
			let lo = (*ri as u128) * (*thetag_lo as u128);
			let hi = (*ri as u128) * (*thetag_hi as u128) + (lo >> 64);
			sum_theta_garner.overflowing_add(U256::from([
				lo as u64,
				hi as u64,
				(hi >> 64) as u64,
				0,
			]));
		}
		// Let's compute v = round(sum_theta_garner / 2^127)
		sum_theta_garner >>= 126;
		let v = sum_theta_garner.as_u128();
		let v = (v & 1) + (v >> 1);

		// If the scaling factor is not 1, compute the inner product with the theta_omega
		let mut w_sign = false;
		let mut w = 0u128;
		if !self.scaling_factor.is_one {
			let mut sum_theta_omega = U256::zero();
			for (thetao_lo, thetao_hi, thetao_sign, ri) in izip!(
				&self.theta_omega_lo,
				&self.theta_omega_hi,
				&self.theta_omega_sign,
				rests
			) {
				let lo = (*ri as u128) * (*thetao_lo as u128);
				let hi = (*ri as u128) * (*thetao_hi as u128) + (lo >> 64);
				if *thetao_sign {
					sum_theta_omega.overflowing_sub(U256::from([
						lo as u64,
						hi as u64,
						(hi >> 64) as u64,
						0,
					]));
				} else {
					sum_theta_omega.overflowing_add(U256::from([
						lo as u64,
						hi as u64,
						(hi >> 64) as u64,
						0,
					]));
				}
			}

			// Let's substract v * theta_gamma to sum_theta_omega.
			let vt_lo_lo = ((v as u64) as u128) * (self.theta_gamma_lo as u128);
			let vt_lo_hi = ((v as u64) as u128) * (self.theta_gamma_hi as u128);
			let vt_hi_lo = ((v >> 64) as u128) * (self.theta_gamma_lo as u128);
			let vt_hi_hi = ((v >> 64) as u128) * (self.theta_gamma_hi as u128);
			let vt_mi =
				(vt_lo_lo >> 64) + ((vt_lo_hi as u64) as u128) + ((vt_hi_lo as u64) as u128);
			let vt_hi = (vt_lo_hi >> 64) + (vt_mi >> 64) + ((vt_hi_hi as u64) as u128);
			if self.theta_gamma_sign {
				sum_theta_omega.overflowing_add(U256::from([
					vt_lo_lo as u64,
					vt_mi as u64,
					vt_hi as u64,
					0,
				]))
			} else {
				sum_theta_omega.overflowing_sub(U256::from([
					vt_lo_lo as u64,
					vt_mi as u64,
					vt_hi as u64,
					0,
				]))
			}

			// Let's compute w = round(sum_theta_omega / 2^127).
			w_sign = sum_theta_omega.msb() > 0;

			if w_sign {
				w = ((!sum_theta_omega) >> 126).as_u128() + 1;
				if !floor {
					w = w.div_ceil(2) - (w & 1);
				} else {
					w = w.div_ceil(2);
				}
			} else {
				w = (sum_theta_omega >> 126).as_u128();
				if !floor {
					w = w.div_floor(2) + (w & 1)
				} else {
					w = w.div_floor(2)
				}
			}
		}

		unsafe {
			for i in 0..out.len() {
				debug_assert!(starting_index + i <= self.to.moduli.len());
				debug_assert!(starting_index + i <= self.omega.len());
				debug_assert!(starting_index + i <= self.omega_shoup.len());
				debug_assert!(starting_index + i <= self.gamma.len());
				debug_assert!(starting_index + i <= self.gamma_shoup.len());
				let out_i = out.get_mut(i).unwrap();
				let qi = self.to.moduli.get_unchecked(starting_index + i);
				let omega_i = self.omega.get_unchecked(starting_index + i);
				let omega_shoup_i = self.omega_shoup.get_unchecked(starting_index + i);
				let gamma_i = self.gamma.get_unchecked(starting_index + i);
				let gamma_shoup_i = self.gamma_shoup.get_unchecked(starting_index + i);

				let mut yi = (qi.modulus() * 2
					- qi.lazy_mul_shoup(qi.lazy_reduce_u128(v), *gamma_i, *gamma_shoup_i))
					as u128;

				if !self.scaling_factor.is_one {
					let wi = qi.lazy_reduce_u128(w);
					yi += if w_sign { qi.modulus() * 2 - wi } else { wi } as u128;
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
	use std::sync::Arc;

	use super::RnsScaler;
	use crate::rns::{scaler::ScalingFactor, RnsContext};
	use ndarray::ArrayView1;
	use num_bigint::BigUint;
	use num_traits::{ToPrimitive, Zero};
	use rand::{thread_rng, RngCore};
	use util::catch_unwind;

	#[test]
	fn test_constructor() -> Result<(), String> {
		let q = Arc::new(RnsContext::new(&[4, 4611686018326724609, 1153])?);

		let scaler = RnsScaler::new(&q, &q, ScalingFactor::one());
		assert_eq!(scaler.from, q);

		assert!(
			catch_unwind(|| ScalingFactor::new(&BigUint::from(1u64), &BigUint::zero())).is_err()
		);
		Ok(())
	}

	#[test]
	fn test_scale_same_context() -> Result<(), String> {
		let ntests = 1000;
		let q = Arc::new(RnsContext::new(&[4u64, 4611686018326724609, 1153])?);
		let mut rng = thread_rng();

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
					let mut x_lift = q.lift(&ArrayView1::from(&x));
					let x_sign = x_lift >= (q.modulus() >> 1);
					if x_sign {
						x_lift = q.modulus() - x_lift;
					}

					let y = scaler.scale_new(&(&x).into(), x.len(), true);
					let x_scaled_floor = if x_sign {
						q.modulus() - (&(&x_lift * &n + &d - 1u64) / &d) % q.modulus()
					} else {
						((&x_lift * &n) / &d) % q.modulus()
					};
					assert_eq!(y, q.project(&x_scaled_floor));

					let z = scaler.scale_new(&(&x).into(), x.len(), false);
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
	fn test_scale_different_contexts() -> Result<(), String> {
		let ntests = 10;
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
		let mut rng = thread_rng();

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

					let mut x_lift = q.lift(&ArrayView1::from(&x));
					let x_sign = x_lift >= (q.modulus() >> 1);
					if x_sign {
						x_lift = q.modulus() - x_lift;
					}

					let y = scaler.scale_new(&(&x).into(), r.moduli.len(), true);
					let x_scaled_floor = if x_sign {
						&r.product - (&(&x_lift * &n + &d - 1u64) / &d) % &r.product
					} else {
						((&x_lift * &n) / &d) % &r.product
					};
					assert_eq!(y, r.project(&x_scaled_floor));
				}
			}
		}
		Ok(())
	}
}
