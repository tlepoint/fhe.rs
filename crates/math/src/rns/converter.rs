#![warn(missing_docs, unused_imports)]

//! Converter from one RNS basis to another RNS basis. The converter uses the
//! RNS Basis Extension algorithm described in Section 2.2 of <https://eprint.iacr.org/2018/117.pdf>.

use super::RnsContext;
use crate::{u256::U256, zq::Modulus};
use itertools::izip;
use ndarray::{ArrayView1, ArrayViewMut1};
use num_bigint::BigUint;
use num_traits::ToPrimitive;

/// Converter from one RNS basis to another.
#[derive(Default, Debug, Clone, PartialEq)]
pub struct RnsConverter {
	from: RnsContext, // Moduli q
	to: RnsContext,   // Moduli p
	garner_mod_p: Vec<Vec<u64>>,
	garner_mod_p_shoup: Vec<Vec<u64>>,
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

		// Store the garner from `from` modulo the moduli of `to`.
		let mut garner_mod_p = Vec::with_capacity(to.moduli_u64.len());
		let mut garner_mod_p_shoup = Vec::with_capacity(to.moduli_u64.len());
		for p in &to.moduli_u64 {
			let mut garner_mod_p_i = Vec::with_capacity(from.moduli_u64.len());
			from.garner
				.iter()
				.for_each(|garner| garner_mod_p_i.push((garner % *p).to_u64().unwrap()));
			garner_mod_p_shoup.push(Modulus::new(*p).unwrap().shoup_vec(&garner_mod_p_i));
			garner_mod_p.push(garner_mod_p_i);
		}

		// Define theta = ceil(2^127 * qtilde_i / q_i) = theta_lo + 2^64 * theta_hi.
		let mut theta_lo = Vec::with_capacity(from.moduli_u64.len());
		let mut theta_hi = Vec::with_capacity(from.moduli_u64.len());
		for (q, qt) in izip!(&from.moduli_u64, &from.q_tilde) {
			let mut theta_lo_i: BigUint = ((BigUint::from(*qt) << 127) + (*q >> 1)) / *q;
			let theta_hi_i: BigUint = &theta_lo_i >> 64;
			theta_lo_i -= &theta_hi_i << 64;
			theta_lo.push(theta_lo_i.to_u64().unwrap());
			theta_hi.push(theta_hi_i.to_u64().unwrap());
		}

		Self {
			from: from.clone(),
			to: to.clone(),
			garner_mod_p,
			garner_mod_p_shoup,
			theta_lo,
			theta_hi,
			q_mod_p,
			q_mod_p_shoup,
		}
	}

	/// Converts  a RNS representation in context `from` and returns a representation in context `to`.
	///
	/// Aborts if the number of rests is different than the number of moduli in debug mode.
	pub fn convert_new(&self, rests_from: &ArrayView1<u64>) -> Vec<u64> {
		let mut rests_to = vec![0; self.to.moduli.len()];
		self.convert(rests_from, &mut (&mut rests_to).into());
		rests_to
	}

	/// Convert a RNS representation in context `from` into a representation in context `to`.
	///
	/// Aborts if the number of rests is different than the number of moduli in debug mode.
	pub fn convert(&self, rests_from: &ArrayView1<u64>, rests_to: &mut ArrayViewMut1<u64>) {
		debug_assert_eq!(rests_from.len(), self.from.moduli_u64.len());

		let mut sum = U256::zero();
		for (rests_from_i, theta_lo, theta_hi) in izip!(rests_from, &self.theta_lo, &self.theta_hi,)
		{
			// Compute rests_from_i * theta
			let lo = (*rests_from_i as u128) * (*theta_lo as u128);
			let mi = (*rests_from_i as u128) * (*theta_hi as u128) + (lo >> 64);
			let hi = mi >> 64;
			sum.overflowing_add(U256::from([lo as u64, mi as u64, hi as u64, 0]));
		}
		sum >>= 126;
		let value = sum.as_u128();
		let value = (value & 1) + (value >> 1);

		for (to, garner_mod_p_j, garner_mod_p_shoup_j, p_j, q_mod_p_j, q_mod_p_shoup_j) in izip!(
			rests_to.iter_mut(),
			&self.garner_mod_p,
			&self.garner_mod_p_shoup,
			&self.to.moduli,
			&self.q_mod_p,
			&self.q_mod_p_shoup,
		) {
			let mut x = (2 * p_j.modulus()
				- p_j.lazy_mul_shoup(p_j.reduce_u128(value), *q_mod_p_j, *q_mod_p_shoup_j))
				as u128;
			for (rests_from_i, garner_mod_p_j_i, garner_mod_p_shoup_j_i) in
				izip!(rests_from, garner_mod_p_j, garner_mod_p_shoup_j)
			{
				x += p_j.lazy_mul_shoup(*rests_from_i, *garner_mod_p_j_i, *garner_mod_p_shoup_j_i)
					as u128
			}
			*to = p_j.reduce_u128(x);
		}
	}
}

#[cfg(test)]
mod tests {
	use std::vec;

	use super::RnsConverter;
	use crate::rns::RnsContext;
	use rand::{thread_rng, RngCore};

	static Q: &[u64; 3] = &[
		4611686018282684417,
		4611686018326724609,
		4611686018309947393,
	];
	static P: &[u64; 5] = &[
		1153,
		4611686018257518593,
		4611686018232352769,
		4611686018171535361,
		4611686018106523649,
	];

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
		let q = RnsContext::new(Q).unwrap();
		let p = RnsContext::new(P).unwrap();
		let converter = RnsConverter::new(&q, &p);

		let a = converter.convert_new(&(&[0, 0, 0]).into());
		assert_eq!(a, &[0, 0, 0, 0, 0]);
		let mut c = vec![1; 5];
		converter.convert(&(&[0, 0, 0]).into(), &mut (&mut c).into());
		assert_eq!(c, &[0, 0, 0, 0, 0]);

		let mut rng = thread_rng();
		for _ in 0..ntests {
			let mut u = vec![];
			for qi in Q {
				u.push(rng.next_u64() % *qi)
			}
			let u_biguint = q.lift(&(&u).into());
			let a = converter.convert_new(&(&u).into());
			let b = if &u_biguint > &(q.modulus() >> 1u64) {
				p.project(&(p.modulus() - q.modulus() + &u_biguint))
			} else {
				p.project(&u_biguint)
			};
			assert_eq!(&a, &b);
			converter.convert(&(&u).into(), &mut (&mut c).into());
			assert_eq!(&c, &b);
		}
	}
}
