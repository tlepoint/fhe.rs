#![warn(missing_docs, unused_imports)]

//! Residue-Number System operations.

use crate::zq::Modulus;
use itertools::izip;
use num_bigint::{BigInt, BigUint, ExtendedGcd, ModInverse};
use num_traits::{cast::ToPrimitive, One, Zero};
use std::cmp::Ordering;

mod converter;
mod i193;
mod scaler;

pub use converter::RnsConverter;
pub use scaler::RnsScaler;

/// Context for a Residue Number System.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RnsContext {
	moduli_u64: Vec<u64>,
	moduli: Vec<Modulus>,
	q_tilde: Vec<u64>,
	q_tilde_shoup: Vec<u64>,
	q_star: Vec<BigUint>,
	garner: Vec<BigUint>,
	product: BigUint,
}

impl RnsContext {
	/// Create a RNS context from a list of moduli.
	///
	/// Returns None if the list is empty, or if the moduli are no coprime.
	pub fn new(moduli_u64: &[u64]) -> std::option::Option<Self> {
		if moduli_u64.is_empty() {
			None
		} else {
			let mut moduli = Vec::with_capacity(moduli_u64.len());
			let mut q_tilde = Vec::with_capacity(moduli_u64.len());
			let mut q_tilde_shoup = Vec::with_capacity(moduli_u64.len());
			let mut q_star = Vec::with_capacity(moduli_u64.len());
			let mut garner = Vec::with_capacity(moduli_u64.len());
			let mut product = BigUint::one();

			for i in 0..moduli_u64.len() {
				// Return None if the moduli are not coprime.
				for j in 0..moduli_u64.len() {
					if i != j {
						let (d, _, _) = BigUint::from(moduli_u64[i])
							.extended_gcd(&BigUint::from(moduli_u64[j]));
						if d.cmp(&BigInt::from(1)) != Ordering::Equal {
							return None;
						}
					}
				}

				product *= &BigUint::from(moduli_u64[i]);
			}

			for modulus in moduli_u64 {
				moduli.push(Modulus::new(*modulus).unwrap());
				// q* = product / modulus
				let q_star_i = &product / modulus;
				// q~ = (product / modulus) ^ (-1) % modulus
				let q_tilde_i = (&product / modulus)
					.mod_inverse(&BigUint::from(*modulus))
					.unwrap();
				// garner = (q*) * (q~)
				let garner_i = (&q_star_i * &q_tilde_i).clone().to_biguint().unwrap();
				q_tilde.push(q_tilde_i.to_u64().unwrap());
				garner.push(garner_i);
				q_star.push(q_star_i);
				q_tilde_shoup.push(
					Modulus::new(*modulus)
						.unwrap()
						.shoup(q_tilde_i.to_u64().unwrap()),
				);
			}

			Some(Self {
				moduli_u64: moduli_u64.to_owned(),
				moduli,
				q_tilde,
				q_tilde_shoup,
				q_star,
				garner,
				product,
			})
		}
	}

	/// Returns the product of the moduli used when creating the RNS context.
	pub fn modulus(&self) -> &BigUint {
		&self.product
	}

	/// Project a BigUint into its rests.
	pub fn project(&self, a: &BigUint) -> Vec<u64> {
		let mut rests = Vec::with_capacity(self.moduli_u64.len());
		for modulus in &self.moduli_u64 {
			rests.push((a % modulus).to_u64().unwrap())
		}
		rests
	}

	/// Lift rests into a BigUint.
	///
	/// Aborts if the number of rests is different than the number of moduli in debug mode.
	pub fn lift(&self, rests: &[u64]) -> BigUint {
		let mut result = BigUint::zero();
		izip!(rests.iter(), self.garner.iter())
			.for_each(|(r_i, garner_i)| result += garner_i * *r_i);
		result % &self.product
	}
}

#[cfg(test)]
mod tests {

	use super::RnsContext;
	use num_bigint::BigUint;
	use rand::RngCore;

	#[test]
	pub fn test_constructor() {
		assert!(RnsContext::new(&[2]).is_some());
		assert!(RnsContext::new(&[2, 3]).is_some());
		assert!(RnsContext::new(&[4, 15, 1153]).is_some());

		assert!(RnsContext::new(&[]).is_none());
		assert!(RnsContext::new(&[2, 4]).is_none());
		assert!(RnsContext::new(&[2, 3, 5, 30]).is_none());
	}

	#[test]
	pub fn test_modulus() {
		let mut rns = RnsContext::new(&[2]).unwrap();
		debug_assert_eq!(rns.modulus(), &BigUint::from(2u64));

		rns = RnsContext::new(&[2, 5]).unwrap();
		debug_assert_eq!(rns.modulus(), &BigUint::from(2u64 * 5));

		rns = RnsContext::new(&[4, 15, 1153]).unwrap();
		debug_assert_eq!(rns.modulus(), &BigUint::from(4u64 * 15 * 1153));
	}

	#[test]
	pub fn test_project_lift() {
		let ntests = 100;
		let rns = RnsContext::new(&[4, 15, 1153]).unwrap();
		let product = 4u64 * 15 * 1153;

		let mut rests = rns.project(&BigUint::from(0u64));
		assert_eq!(&rests, &[0u64, 0, 0]);
		assert_eq!(rns.lift(&rests), BigUint::from(0u64));

		rests = rns.project(&BigUint::from(4u64));
		assert_eq!(&rests, &[0u64, 4, 4]);
		assert_eq!(rns.lift(&rests), BigUint::from(4u64));

		rests = rns.project(&BigUint::from(15u64));
		assert_eq!(&rests, &[3u64, 0, 15]);
		assert_eq!(rns.lift(&rests), BigUint::from(15u64));

		rests = rns.project(&BigUint::from(1153u64));
		assert_eq!(&rests, &[1u64, 13, 0]);
		assert_eq!(rns.lift(&rests), BigUint::from(1153u64));

		rests = rns.project(&BigUint::from(product - 1));
		assert_eq!(&rests, &[3u64, 14, 1152]);
		assert_eq!(rns.lift(&rests), BigUint::from(product - 1));

		let mut rng = rand::thread_rng();

		for _ in 0..ntests {
			let b = BigUint::from(rng.next_u64() % product);
			rests = rns.project(&b);
			assert_eq!(rns.lift(&rests), b);
		}
	}
}
