#![warn(missing_docs, unused_imports)]

//! Converter from one RNS basis to another RNS basis. The converter uses the
//! RNS Basis Extension algorithm described in Section 2.2 of <https://eprint.iacr.org/2018/117.pdf>.

use super::RnsContext;
use crate::{u256::U256, zq::Modulus};
use itertools::izip;
use num_bigint::BigUint;
use num_traits::ToPrimitive;

/// Converter from one RNS basis to another.
#[derive(Debug, Clone, PartialEq)]
pub struct RnsConverter {
	from: RnsContext, // Moduli q
	to: RnsContext,   // Moduli p
	q_star_mod_p: Vec<Vec<u64>>,
	q_star_mod_p_shoup: Vec<Vec<u64>>,
	theta_lo: Vec<u64>,
	theta_hi: Vec<u64>,
	q_mod_p: Vec<u64>,
	q_mod_p_shoup: Vec<u64>,
}

impl RnsConverter {
	/// Create a RNS basis converter between two RNS contexts.
	pub fn new(from: &RnsContext, to: &RnsContext) -> Self {
		// Store the product of the moduli in `from` modulo the moduli of `to`.
		let mut q_mod_p = Vec::with_capacity(to.moduli_u64.len());
		let mut q_mod_p_shoup = Vec::with_capacity(to.moduli_u64.len());
		to.moduli_u64
			.iter()
			.for_each(|modulus| q_mod_p.push((&from.product % *modulus).to_u64().unwrap()));
		izip!(&q_mod_p, &to.moduli_u64)
			.for_each(|(qi, p)| q_mod_p_shoup.push(Modulus::new(*p).unwrap().shoup(*qi)));

		// Store the q_star from `from` modulo the moduli of `to`.
		let mut q_star_mod_p = Vec::with_capacity(to.moduli_u64.len());
		let mut q_star_mod_p_shoup = Vec::with_capacity(to.moduli_u64.len());
		for p in &to.moduli_u64 {
			let mut q_star_mod_p_i = Vec::with_capacity(from.moduli_u64.len());
			from.q_star
				.iter()
				.for_each(|q_star| q_star_mod_p_i.push((q_star % *p).to_u64().unwrap()));
			q_star_mod_p_shoup.push(Modulus::new(*p).unwrap().shoup_vec(&q_star_mod_p_i));
			q_star_mod_p.push(q_star_mod_p_i);
		}

		// Define theta = floor(2^128 / q_i) = theta_lo + 2^64 * theta_hi.
		let mut theta_lo = Vec::with_capacity(from.moduli_u64.len());
		let mut theta_hi = Vec::with_capacity(from.moduli_u64.len());
		for q in &from.moduli_u64 {
			let mut theta_lo_i: BigUint = (BigUint::from(1u64) << 128) / *q;
			let theta_hi_i: BigUint = &theta_lo_i >> 64;
			theta_lo_i -= &theta_hi_i << 64;
			theta_lo.push(theta_lo_i.to_u64().unwrap());
			theta_hi.push(theta_hi_i.to_u64().unwrap());
		}

		Self {
			from: from.clone(),
			to: to.clone(),
			q_star_mod_p,
			q_star_mod_p_shoup,
			theta_lo,
			theta_hi,
			q_mod_p,
			q_mod_p_shoup,
		}
	}

	/// Convert a RNS representation in context `from` into a representation in context `to`.
	///
	/// Aborts if the number of rests is different than the number of moduli in debug mode.
	pub fn convert(&self, rests_from: &[u64]) -> Vec<u64> {
		debug_assert_eq!(rests_from.len(), self.from.moduli_u64.len());

		let mut rests_to = Vec::with_capacity(self.to.moduli_u64.len());

		let mut y = Vec::with_capacity(rests_from.len());
		let mut sum = U256::zero();
		for (rests_from_i, theta_lo, theta_hi, q_tilde, q_tilde_shoup, qi) in izip!(
			rests_from,
			&self.theta_lo,
			&self.theta_hi,
			&self.from.q_tilde,
			&self.from.q_tilde_shoup,
			&self.from.moduli,
		) {
			let yi = qi.mul_shoup(*rests_from_i, *q_tilde, *q_tilde_shoup);
			y.push(yi);
			// Compute yi * theta
			let lo = (*theta_lo as u128) * (yi as u128);
			let hi = (*theta_hi as u128) * (yi as u128) + (lo >> 64);
			sum.overflowing_add(U256::from([lo as u64, hi as u64, (hi >> 64) as u64, 0]));
		}
		sum >>= 128;
		let value2 = sum.as_u64();

		for (q_star_mod_p_j, q_star_mod_p_shoup_j, p_j, q_mod_p_j, q_mod_p_shoup_j) in izip!(
			&self.q_star_mod_p,
			&self.q_star_mod_p_shoup,
			&self.to.moduli,
			&self.q_mod_p,
			&self.q_mod_p_shoup,
		) {
			let mut x = (2 * p_j.modulus()
				- p_j.lazy_mul_shoup(value2, *q_mod_p_j, *q_mod_p_shoup_j)) as u128;
			izip!(&y, q_star_mod_p_j, q_star_mod_p_shoup_j).for_each(
				|(yi, q_star_mod_p_j_i, q_star_mod_p_shoup_j_i)| {
					x += p_j.lazy_mul_shoup(*yi, *q_star_mod_p_j_i, *q_star_mod_p_shoup_j_i) as u128
				},
			);
			rests_to.push(p_j.reduce_u128(x));
		}

		rests_to
	}
}

#[cfg(test)]
mod tests {
	use super::RnsConverter;
	use crate::rns::RnsContext;
	use num_bigint::BigUint;
	use rand::{thread_rng, RngCore};

	#[test]
	fn test_constructor() {
		let q = RnsContext::new(&[4, 15, 1153]).unwrap();
		let p = RnsContext::new(&[7, 13, 907]).unwrap();
		let converter = RnsConverter::new(&q, &p);
		assert_eq!(converter.from, q);
		assert_eq!(converter.to, p);
	}

	#[test]
	fn test_converter() {
		let ntests = 100;
		let q = RnsContext::new(&[4, 15, 1153]).unwrap();
		let q_product = 4u64 * 15 * 1153;
		let p = RnsContext::new(&[7, 13, 907]).unwrap();
		let converter = RnsConverter::new(&q, &p);

		let a = converter.convert(&[0, 0, 0]);
		assert_eq!(a, &[0, 0, 0]);

		let mut rng = thread_rng();
		for _ in 0..ntests {
			let u = BigUint::from(rng.next_u64() % q_product);
			let a = converter.convert(&q.project(&u));
			let b = p.project(&u);
			assert_eq!(&a, &b);
		}
	}
}
