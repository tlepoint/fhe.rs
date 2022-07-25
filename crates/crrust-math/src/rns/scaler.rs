#![warn(missing_docs, unused_imports)]

//! RNS scaler inspired from Remark 3.2 of <https://eprint.iacr.org/2021/204.pdf>.

use super::RnsContext;
use crate::u256::U256;
use itertools::{izip, Itertools};
use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};

/// Scaler in RNS basis.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RnsScaler {
	ctx: RnsContext,
	numerator: BigUint,
	denominator: BigUint,

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
	pub fn new(ctx: &RnsContext, numerator: &BigUint, denominator: &BigUint) -> Self {
		assert_ne!(denominator, &BigUint::zero());

		// Let's define gamma = round(numerator * product / denominator)
		let (gamma, theta_gamma_lo, theta_gamma_hi, theta_gamma_sign) =
			Self::extract_projection_and_theta(ctx, &ctx.product, numerator, denominator, false);
		let gamma_shoup = izip!(&gamma, &ctx.moduli)
			.map(|(wi, q)| q.shoup(*wi))
			.collect_vec();

		// Let's define omega_i = round(garner_i * numerator / denominator)
		let mut omega = vec![];
		let mut omega_shoup = vec![];
		for _ in &ctx.moduli {
			omega.push(vec![0u64; ctx.moduli.len()]);
			omega_shoup.push(vec![0u64; ctx.moduli.len()]);
		}
		let mut theta_omega_lo = vec![];
		let mut theta_omega_hi = vec![];
		let mut theta_omega_sign = vec![];
		for i in 0..ctx.garner.len() {
			let (omega_i, theta_omega_i_lo, theta_omega_i_hi, theta_omega_i_sign) =
				Self::extract_projection_and_theta(
					ctx,
					&ctx.garner[i],
					numerator,
					denominator,
					true,
				);
			for j in 0..ctx.moduli.len() {
				let qj = &ctx.moduli[j];
				omega[j][i] = qj.reduce(omega_i[j]);
				omega_shoup[j][i] = qj.shoup(omega[j][i]);
			}
			theta_omega_lo.push(theta_omega_i_lo);
			theta_omega_hi.push(theta_omega_i_hi);
			theta_omega_sign.push(theta_omega_i_sign);
		}

		// Finally, define theta_garner_i = garner_i / product, also scaled by 2^127.
		let mut theta_garner_lo = vec![];
		let mut theta_garner_hi = vec![];
		for garner_i in &ctx.garner {
			let mut theta = ((garner_i << 127) + (&ctx.product >> 1)) / &ctx.product;
			let theta_hi = &theta >> 64;
			theta -= &theta_hi << 64;
			theta_garner_lo.push(theta.to_u64().unwrap());
			theta_garner_hi.push(theta_hi.to_u64().unwrap());
		}

		Self {
			ctx: ctx.clone(),
			numerator: numerator.clone(),
			denominator: denominator.clone(),
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
	// This function projects gamma in the RNS context, and scales theta_gamma by 2**128 and round.
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
		let theta_hi_biguint = &theta >> 64;
		theta -= &theta_hi_biguint << 64;
		let theta_lo = theta.to_u64().unwrap();
		let theta_hi = theta_hi_biguint.to_u64().unwrap();

		(projected, theta_lo, theta_hi, theta_sign)
	}

	/// Output the RNS representation of the underlying value scaled by numerator * denominator,
	/// and either rounded or floored.
	///
	///
	/// Aborts if the number of rests is different than the number of moduli in debug mode, or
	/// if the size is not in [1, ..., rests.len()].
	pub fn scale(&self, rests: &[u64], size: usize, floor: bool) -> Vec<u64> {
		debug_assert_eq!(rests.len(), self.ctx.moduli_u64.len());
		debug_assert!(size >= 1 && size <= rests.len());

		let mut out = Vec::with_capacity(size);

		// First, let's compute the inner product of the rests with theta_omega and theta_garner.
		let mut sum_theta_garner = U256::zero();
		let mut sum_theta_omega = U256::zero();
		for (thetag_lo, thetag_hi, thetao_lo, thetao_hi, thetao_sign, ri) in izip!(
			&self.theta_garner_lo,
			&self.theta_garner_hi,
			&self.theta_omega_lo,
			&self.theta_omega_hi,
			&self.theta_omega_sign,
			rests
		) {
			let mut lo = (*ri as u128) * (*thetag_lo as u128);
			let mut hi = (*ri as u128) * (*thetag_hi as u128) + (lo >> 64);
			// sum_theta_garner.add(lo as u64, hi as u64, (hi >> 64) as u64, false);
			sum_theta_garner.overflowing_add(U256::from([
				lo as u64,
				hi as u64,
				(hi >> 64) as u64,
				0,
			]));
			lo = (*ri as u128) * (*thetao_lo as u128);
			hi = (*ri as u128) * (*thetao_hi as u128) + (lo >> 64);
			// sum_theta_omega.add(lo as u64, hi as u64, (hi >> 64) as u64, *thetao_sign);
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

		// Let's compute v = floor(sum_theta_garner / 2^127)
		sum_theta_garner >>= 127;
		let v = sum_theta_garner.as_u128();

		// Let's substract v * theta_gamma to sum_theta_omega.
		let vt_lo_lo = ((v as u64) as u128) * (self.theta_gamma_lo as u128);
		let vt_lo_hi = ((v as u64) as u128) * (self.theta_gamma_hi as u128);
		let vt_hi_lo = ((v >> 64) as u128) * (self.theta_gamma_lo as u128);
		let vt_hi_hi = ((v >> 64) as u128) * (self.theta_gamma_hi as u128);
		let vt_mi = (vt_lo_lo >> 64) + ((vt_lo_hi as u64) as u128) + ((vt_hi_lo as u64) as u128);
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
		let w_sign = sum_theta_omega.msb() > 0;
		let mut w: u128;
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

		// Let's compute [ sum(r_j * omega_j) - v * gamma + w] mod q_i
		for i in 0..size {
			let qi = &self.ctx.moduli[i];
			let vi = qi.lazy_reduce_u128(v);
			let wi = qi.lazy_reduce_u128(w);
			let mut yi = if w_sign { qi.modulus() * 2 - wi } else { wi } as u128;

			yi += (qi.modulus() * 2 - qi.lazy_mul_shoup(vi, self.gamma[i], self.gamma_shoup[i]))
				as u128;

			for (ri, omega_i, omega_shoup_i) in izip!(rests, &self.omega[i], &self.omega_shoup[i]) {
				yi += qi.lazy_mul_shoup(*ri, *omega_i, *omega_shoup_i) as u128;
			}
			out.push(qi.reduce_u128(yi));
		}

		out
	}
}

#[cfg(test)]
mod tests {
	use super::RnsScaler;
	use crate::rns::RnsContext;
	use num_bigint::BigUint;
	use num_traits::Zero;
	use rand::{thread_rng, RngCore};
	use std::panic::UnwindSafe;

	// Redefine catch_unwind to silence the panic.
	pub fn catch_unwind<F, R>(f: F) -> std::thread::Result<R>
	where
		F: FnOnce() -> R + UnwindSafe,
	{
		let prev_hook = std::panic::take_hook();
		std::panic::set_hook(Box::new(|_| {}));
		let r = std::panic::catch_unwind(f);
		std::panic::set_hook(prev_hook);
		r
	}

	#[test]
	fn test_constructor() {
		let q = RnsContext::new(&[4, 4611686018326724609, 1153]).unwrap();

		let scaler = RnsScaler::new(&q, &BigUint::from(1u64), &BigUint::from(1u64));
		assert_eq!(scaler.ctx, q);

		assert!(
			catch_unwind(|| RnsScaler::new(&q, &BigUint::from(1u64), &BigUint::zero())).is_err()
		);
	}

	#[test]
	fn test_identity() {
		let ntests = 100;
		let q = RnsContext::new(&[4, 4611686018326724609, 1153]).unwrap();
		let identity = RnsScaler::new(&q, &BigUint::from(1u64), &BigUint::from(1u64));

		let mut rng = thread_rng();
		for _ in 0..ntests {
			let x = vec![
				rng.next_u64() % q.moduli_u64[0],
				rng.next_u64() % q.moduli_u64[1],
				rng.next_u64() % q.moduli_u64[2],
			];
			let y = identity.scale(&x, x.len(), true);
			let z = identity.scale(&x, x.len(), true);
			assert_eq!(&x, &y);
			assert_eq!(&x, &z);
		}
	}

	#[test]
	fn test_scale_up() {
		let ntests = 100;
		let q = RnsContext::new(&[4u64, 4611686018326724609, 1153]).unwrap();
		let mut rng = thread_rng();

		for numerator in &[1u64, 2, 3, 100, 1000, 4611686018326724610] {
			let scaler = RnsScaler::new(&q, &BigUint::from(*numerator), &BigUint::from(1u64));

			for _ in 0..ntests {
				let x = vec![
					rng.next_u64() % q.moduli_u64[0],
					rng.next_u64() % q.moduli_u64[1],
					rng.next_u64() % q.moduli_u64[2],
				];
				let y = scaler.scale(&x, x.len(), true);
				let z = scaler.scale(&x, x.len(), false);

				let x_lift = q.lift(&x);
				let x_scaled = &x_lift * *numerator;
				assert_eq!(y, q.project(&x_scaled));
				assert_eq!(z, q.project(&x_scaled));
			}
		}
	}

	#[test]
	fn test_scale() {
		let ntests = 100;
		let q = RnsContext::new(&[4u64, 4611686018326724609, 1153]).unwrap();
		let mut rng = thread_rng();

		for numerator in &[1u64, 2, 3, 100, 1000, 4611686018326724610] {
			for denominator in &[1u64, 2, 3, 101, 1001, 4611686018326724610] {
				let n = BigUint::from(*numerator);
				let d = BigUint::from(*denominator);
				let scaler = RnsScaler::new(&q, &n, &d);

				for _ in 0..ntests {
					let x = vec![
						rng.next_u64() % q.moduli_u64[0],
						rng.next_u64() % q.moduli_u64[1],
						rng.next_u64() % q.moduli_u64[2],
					];
					let x_lift = q.lift(&x);

					let y = scaler.scale(&x, x.len(), true);
					let x_scaled_floor = (&x_lift * &n) / &d;
					assert_eq!(y, q.project(&x_scaled_floor));

					let z = scaler.scale(&x, x.len(), false);
					let x_scaled_round = (&x_lift * &n + (&d >> 1)) / &d;
					assert_eq!(z, q.project(&x_scaled_round));
				}
			}
		}
	}
}
