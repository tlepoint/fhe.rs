#![warn(missing_docs, unused_imports)]

//! Polynomials in R_q\[x\] = (ZZ_q1 x ... x ZZ_qn)\[x\] where the qi's are prime moduli in zq.

mod convert;
mod ops;

pub mod extender;
pub mod scaler;
pub mod traits;

use crate::{
	rns::RnsContext,
	zq::{ntt::NttOperator, Modulus},
};
use itertools::izip;
use ndarray::{s, Array2, ArrayView2};
use num_bigint::BigUint;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use std::sync::Arc;
use traits::TryConvertFrom;
use util::sample_vec_cbd;
use zeroize::Zeroize;

/// Struct that holds the context associated with elements in rq.
#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub struct Context {
	q: Vec<Modulus>,
	rns: Arc<RnsContext>,
	ops: Vec<NttOperator>,
	degree: usize,
}

impl Context {
	/// Creates a context from a list of moduli and a polynomial degree.
	///
	/// Returns an error if the moduli are not primes less than 62 bits which supports the NTT of size `degree`.
	pub fn new(moduli: &[u64], degree: usize) -> Result<Self, String> {
		if !degree.is_power_of_two() || degree < 8 {
			Err("The degree is not a power of two larger or equal to 8".to_string())
		} else {
			let mut q = Vec::with_capacity(moduli.len());
			let rns = Arc::new(RnsContext::new(moduli)?);
			let mut ops = Vec::with_capacity(moduli.len());
			for modulus in moduli {
				let qi = Modulus::new(*modulus)?;
				if let Some(op) = NttOperator::new(&qi, degree) {
					q.push(qi);
					ops.push(op);
				} else {
					return Err("Impossible to construct a Ntt operator".to_string());
				}
			}

			Ok(Self {
				q,
				rns,
				ops,
				degree,
			})
		}
	}

	/// Returns the modulus as a BigUint
	pub fn modulus(&self) -> &BigUint {
		self.rns.modulus()
	}
}

/// Possible representations of the underlying polynomial.
#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub enum Representation {
	/// This is the list of coefficients ci, such that the polynomial is c0 + c1 * x + ... + c_(degree - 1) * x^(degree - 1)
	#[default]
	PowerBasis,
	/// This is the NTT representation of the PowerBasis representation.
	Ntt,
	/// This is a "Shoup" representation of the Ntt representation used for faster multiplication.
	NttShoup,
}

/// Struct that holds a polynomial for a specific context.
#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub struct Poly {
	ctx: Arc<Context>,
	representation: Representation,
	allow_variable_time_computations: bool,
	coefficients: Array2<u64>,
	coefficients_shoup: Option<Array2<u64>>,
}

impl Poly {
	/// Creates a polynomial holding the constant 0.
	pub fn zero(ctx: &Arc<Context>, representation: Representation) -> Self {
		Self {
			ctx: ctx.clone(),
			representation: representation.clone(),
			allow_variable_time_computations: false,
			coefficients: Array2::zeros((ctx.q.len(), ctx.degree)),
			coefficients_shoup: if representation == Representation::NttShoup {
				Some(Array2::zeros((ctx.q.len(), ctx.degree)))
			} else {
				None
			},
		}
	}

	/// Enable variable time computations when this polynomial is involved.
	///
	/// # Safety
	///
	/// By default, this is marked as unsafe, but is usually safe when only public data is processed.
	pub unsafe fn allow_variable_time_computations(&mut self) {
		self.allow_variable_time_computations = true
	}

	/// Disable variable time computations when this polynomial is involved.
	pub fn disallow_variable_time_computations(&mut self) {
		self.allow_variable_time_computations = false
	}

	/// Change the representation of the underlying polynomial.
	pub fn change_representation(&mut self, to: Representation) {
		match self.representation {
			Representation::PowerBasis => {
				match to {
					Representation::Ntt => self.ntt_forward(),
					Representation::NttShoup => {
						self.ntt_forward();
						self.compute_coefficients_shoup();
					}
					Representation::PowerBasis => {} // no-op
				}
			}
			Representation::Ntt => {
				match to {
					Representation::PowerBasis => self.ntt_backward(),
					Representation::NttShoup => self.compute_coefficients_shoup(),
					Representation::Ntt => {} // no-op
				}
			}
			Representation::NttShoup => {
				if to != Representation::NttShoup {
					// We are not sure whether this polynomial was sensitive or not,
					// so for security, we zeroize the Shoup coefficients.
					self.coefficients_shoup
						.as_mut()
						.unwrap()
						.as_slice_mut()
						.unwrap()
						.zeroize();
					self.coefficients_shoup = None
				}
				match to {
					Representation::PowerBasis => self.ntt_backward(),
					Representation::Ntt => {}      // no-op
					Representation::NttShoup => {} // no-op
				}
			}
		}

		self.representation = to;
	}

	/// Compute the Shoup representation of the coefficients.
	fn compute_coefficients_shoup(&mut self) {
		let mut coefficients_shoup = Array2::zeros((self.ctx.q.len(), self.ctx.degree));
		izip!(
			coefficients_shoup.outer_iter_mut(),
			self.coefficients.outer_iter(),
			&self.ctx.q
		)
		.for_each(|(mut v_shoup, v, qi)| {
			v_shoup
				.as_slice_mut()
				.unwrap()
				.copy_from_slice(&qi.shoup_vec(v.as_slice().unwrap()))
		});
		self.coefficients_shoup = Some(coefficients_shoup)
	}

	/// Override the internal representation to a given representation.
	///
	/// # Safety
	///
	/// Prefer the `change_representation` function to safely modify the polynomial representation.
	/// If the `to` representation is NttShoup, the coefficients are still computed correctly to
	/// avoid being in an unstable state. Similarly, if we override a representation which was
	/// NttShoup, we zeroize the existing Shoup coefficients.
	pub unsafe fn override_representation(&mut self, to: Representation) {
		if to == Representation::NttShoup {
			self.compute_coefficients_shoup()
		} else if self.coefficients_shoup.is_some() {
			self.coefficients_shoup
				.as_mut()
				.unwrap()
				.as_slice_mut()
				.unwrap()
				.zeroize();
			self.coefficients_shoup = None
		}
		self.representation = to;
	}

	/// Generate a random polynomial.
	pub fn random(ctx: &Arc<Context>, representation: Representation) -> Self {
		let mut p = Poly::zero(ctx, representation);
		izip!(p.coefficients.outer_iter_mut(), &ctx.q).for_each(|(mut v, qi)| {
			v.as_slice_mut()
				.unwrap()
				.copy_from_slice(&qi.random_vec(ctx.degree))
		});
		if p.representation == Representation::NttShoup {
			p.compute_coefficients_shoup()
		}
		p
	}

	/// Generate a random polynomial deterministically from a seed.
	pub fn random_from_seed(
		ctx: &Arc<Context>,
		representation: Representation,
		seed: <ChaCha8Rng as SeedableRng>::Seed,
	) -> Self {
		let mut rng = ChaCha8Rng::from_seed(seed);
		let mut p = Poly::zero(ctx, representation);
		izip!(p.coefficients.outer_iter_mut(), &ctx.q).for_each(|(mut v, qi)| {
			let mut seed_for_vec = <ChaCha8Rng as SeedableRng>::Seed::default();
			rng.fill(&mut seed_for_vec);
			v.as_slice_mut()
				.unwrap()
				.copy_from_slice(&qi.random_vec_from_seed(ctx.degree, seed_for_vec))
		});
		if p.representation == Representation::NttShoup {
			p.compute_coefficients_shoup()
		}
		p
	}

	/// Generate a small polynomial and convert into the specified representation.
	///
	/// Returns an error if the variance does not belong to [1, ..., 16].
	pub fn small(
		ctx: &Arc<Context>,
		representation: Representation,
		variance: usize,
	) -> Result<Self, String> {
		if !(1..=16).contains(&variance) {
			Err("The variance should be an integer between 1 and 16".to_string())
		} else {
			let mut coeffs = sample_vec_cbd(ctx.degree, variance)?;
			let mut p =
				Poly::try_convert_from(coeffs.as_ref() as &[i64], ctx, Representation::PowerBasis)?;
			if representation != Representation::PowerBasis {
				p.change_representation(representation);
			}
			coeffs.zeroize();
			Ok(p)
		}
	}

	/// Access the polynomial coefficients in RNS representation.
	pub fn coefficients(&self) -> ArrayView2<u64> {
		self.coefficients.view()
	}

	/// Computes the forward Ntt on the coefficients
	fn ntt_forward(&mut self) {
		if self.allow_variable_time_computations {
			izip!(self.coefficients.outer_iter_mut(), &self.ctx.ops)
				.for_each(|(mut v, op)| unsafe { op.forward_vt(v.as_slice_mut().unwrap()) });
		} else {
			izip!(self.coefficients.outer_iter_mut(), &self.ctx.ops)
				.for_each(|(mut v, op)| op.forward(v.as_slice_mut().unwrap()));
		}
	}

	/// Computes the backward Ntt on the coefficients
	fn ntt_backward(&mut self) {
		if self.allow_variable_time_computations {
			izip!(self.coefficients.outer_iter_mut(), &self.ctx.ops)
				.for_each(|(mut v, op)| unsafe { op.backward_vt(v.as_slice_mut().unwrap()) });
		} else {
			izip!(self.coefficients.outer_iter_mut(), &self.ctx.ops)
				.for_each(|(mut v, op)| op.backward(v.as_slice_mut().unwrap()));
		}
	}

	/// Substitute x by x^i in a polynomial.
	/// In PowerBasis representation, i can be any integer that is not a multiple of 2 * degree.
	/// In Ntt and NttShoup representation, i can be any odd integer that is not a multiple of 2 * degree.
	pub fn substitute(&self, i: usize) -> Result<Poly, String> {
		let degree = self.ctx.degree as u32;
		let exponent = (i as u32) % (2 * degree);
		if exponent == 0 {
			return Err("The exponent is a multiple of 2 * degree".to_string());
		}

		let mut q = Poly::zero(&self.ctx, self.representation.clone());
		if self.allow_variable_time_computations {
			unsafe { q.allow_variable_time_computations() }
		}
		let mask = (self.ctx.degree - 1) as u32;
		match self.representation {
			Representation::Ntt => {
				if exponent & 1 == 0 {
					return Err("The exponent should be odd modulo 2 * degree".to_string());
				}

				let mut power = (exponent - 1) / 2;
				for j in 0..degree {
					let j_bitrev = (j.reverse_bits() >> (degree.leading_zeros() + 1)) as usize;
					let power_bitrev =
						((power & mask).reverse_bits() >> (degree.leading_zeros() + 1)) as usize;
					q.coefficients
						.slice_mut(s![.., j_bitrev])
						.assign(&self.coefficients.slice(s![.., power_bitrev]));
					power += exponent;
				}
			}
			Representation::NttShoup => {
				if exponent & 1 == 0 {
					return Err("The exponent should be odd modulo 2 * degree".to_string());
				}

				let mut power = (exponent - 1) / 2;
				for j in 0..degree {
					let j_bitrev = (j.reverse_bits() >> (degree.leading_zeros() + 1)) as usize;
					let power_bitrev =
						(power.reverse_bits() >> (degree.leading_zeros() + 1)) as usize;
					q.coefficients
						.slice_mut(s![.., j_bitrev])
						.assign(&self.coefficients.slice(s![.., power_bitrev]));
					q.coefficients_shoup
						.as_mut()
						.unwrap()
						.slice_mut(s![.., j_bitrev])
						.assign(
							&self
								.coefficients_shoup
								.as_ref()
								.unwrap()
								.slice(s![.., power_bitrev]),
						);
					power += exponent;
				}
			}
			Representation::PowerBasis => {
				let mut power = 0u32;
				for j in 0..degree {
					izip!(
						&self.ctx.q,
						q.coefficients.slice_mut(s![.., (power & mask) as usize]),
						self.coefficients.slice(s![.., j as usize])
					)
					.for_each(|(qi, qij, pij)| {
						if (power / degree) & 1 == 1 {
							*qij = qi.sub(*qij, *pij)
						} else {
							*qij = qi.add(*qij, *pij)
						}
					});
					power += exponent
				}
			}
		}

		Ok(q)
	}
}

impl Zeroize for Poly {
	fn zeroize(&mut self) {
		self.coefficients.as_slice_mut().unwrap().zeroize();
		if let Some(s) = self.coefficients_shoup.as_mut() {
			s.as_slice_mut().unwrap().zeroize();
		}
	}
}

#[cfg(test)]
mod tests {
	extern crate test;
	use super::{Context, Poly, Representation};
	use crate::zq::{ntt::supports_ntt, Modulus};
	use itertools::Itertools;
	use num_bigint::BigUint;
	use num_traits::{One, Zero};
	use rand::{thread_rng, Rng, SeedableRng};
	use rand_chacha::ChaCha8Rng;
	use std::sync::Arc;

	// Moduli to be used in tests.
	static MODULI: &[u64; 3] = &[1153, 4611686018326724609, 4611686018309947393];

	#[test]
	fn test_context_constructor() {
		for modulus in MODULI {
			assert!(Context::new(&[*modulus], 8).is_ok());
			if supports_ntt(*modulus, 128) {
				assert!(Context::new(&[*modulus], 128).is_ok());
			} else {
				assert!(Context::new(&[*modulus], 128).is_err());
			}
		}

		assert!(Context::new(MODULI, 8).is_ok());
		assert!(Context::new(MODULI, 128).is_err());
	}

	#[test]
	fn test_poly_zero() -> Result<(), String> {
		let reference = &[
			BigUint::zero(),
			BigUint::zero(),
			BigUint::zero(),
			BigUint::zero(),
			BigUint::zero(),
			BigUint::zero(),
			BigUint::zero(),
			BigUint::zero(),
		];

		for modulus in MODULI {
			let ctx = Arc::new(Context::new(&[*modulus], 8)?);
			let p = Poly::zero(&ctx, Representation::PowerBasis);
			let q = Poly::zero(&ctx, Representation::Ntt);
			assert_ne!(p, q);
			assert_eq!(Vec::<u64>::from(&p), &[0; 8]);
			assert_eq!(Vec::<u64>::from(&q), &[0; 8]);
		}

		let ctx = Arc::new(Context::new(MODULI, 8)?);
		let p = Poly::zero(&ctx, Representation::PowerBasis);
		let q = Poly::zero(&ctx, Representation::Ntt);
		assert_ne!(p, q);
		assert_eq!(Vec::<u64>::from(&p), [0; 24]);
		assert_eq!(Vec::<u64>::from(&q), [0; 24]);
		assert_eq!(Vec::<BigUint>::from(&p), reference);
		assert_eq!(Vec::<BigUint>::from(&q), reference);

		Ok(())
	}

	#[test]
	fn test_random() -> Result<(), String> {
		for _ in 0..100 {
			let mut seed = <ChaCha8Rng as SeedableRng>::Seed::default();
			thread_rng().fill(&mut seed);

			for modulus in MODULI {
				let ctx = Arc::new(Context::new(&[*modulus], 8)?);
				let p = Poly::random_from_seed(&ctx, Representation::Ntt, seed.clone());
				let q = Poly::random_from_seed(&ctx, Representation::Ntt, seed.clone());
				assert_eq!(p, q);
			}

			let ctx = Arc::new(Context::new(MODULI, 8)?);
			let p = Poly::random_from_seed(&ctx, Representation::Ntt, seed.clone());
			let q = Poly::random_from_seed(&ctx, Representation::Ntt, seed.clone());
			assert_eq!(p, q);

			thread_rng().fill(&mut seed);
			let p = Poly::random_from_seed(&ctx, Representation::Ntt, seed);
			assert_ne!(p, q);

			let r = Poly::random(&ctx, Representation::Ntt);
			assert_ne!(p, r);
			assert_ne!(q, r);
		}
		Ok(())
	}

	#[test]
	fn test_coefficients() -> Result<(), String> {
		for _ in 0..50 {
			for modulus in MODULI {
				let ctx = Arc::new(Context::new(&[*modulus], 8)?);
				let p = Poly::random(&ctx, Representation::Ntt);
				let p_coefficients = Vec::<u64>::from(&p);
				assert_eq!(p_coefficients, p.coefficients().as_slice().unwrap())
			}

			let ctx = Arc::new(Context::new(MODULI, 8)?);
			let p = Poly::random(&ctx, Representation::Ntt);
			let p_coefficients = Vec::<u64>::from(&p);
			assert_eq!(p_coefficients, p.coefficients().as_slice().unwrap())
		}
		Ok(())
	}

	#[test]
	fn test_modulus() -> Result<(), String> {
		for modulus in MODULI {
			let modulus_biguint = BigUint::from(*modulus);
			let ctx = Arc::new(Context::new(&[*modulus], 8)?);
			assert_eq!(ctx.modulus(), &modulus_biguint)
		}

		let mut modulus_biguint = BigUint::one();
		MODULI.iter().for_each(|m| modulus_biguint *= *m);
		let ctx = Arc::new(Context::new(MODULI, 8)?);
		assert_eq!(ctx.modulus(), &modulus_biguint);

		Ok(())
	}

	#[test]
	fn test_allow_variable_time_computations() -> Result<(), String> {
		for modulus in MODULI {
			let ctx = Arc::new(Context::new(&[*modulus], 8)?);
			let mut p = Poly::random(&ctx, Representation::default());
			assert!(!p.allow_variable_time_computations);

			unsafe { p.allow_variable_time_computations() }
			assert!(p.allow_variable_time_computations);

			let q = p.clone();
			assert!(q.allow_variable_time_computations);

			p.disallow_variable_time_computations();
			assert!(!p.allow_variable_time_computations);
		}

		let ctx = Arc::new(Context::new(MODULI, 8)?);
		let mut p = Poly::random(&ctx, Representation::default());
		assert!(!p.allow_variable_time_computations);

		unsafe { p.allow_variable_time_computations() }
		assert!(p.allow_variable_time_computations);

		let q = p.clone();
		assert!(q.allow_variable_time_computations);

		// Allowing variable time propagates.
		let mut p = Poly::random(&ctx, Representation::Ntt);
		unsafe { p.allow_variable_time_computations() }
		let mut q = Poly::random(&ctx, Representation::Ntt);

		assert!(!q.allow_variable_time_computations);
		q *= &p;
		assert!(q.allow_variable_time_computations);

		q.disallow_variable_time_computations();
		q += &p;
		assert!(q.allow_variable_time_computations);

		q.disallow_variable_time_computations();
		q -= &p;
		assert!(q.allow_variable_time_computations);

		q = -&p;
		assert!(q.allow_variable_time_computations);

		Ok(())
	}

	#[test]
	fn test_change_representation() -> Result<(), String> {
		let ctx = Arc::new(Context::new(MODULI, 8)?);

		let mut p = Poly::random(&ctx, Representation::default());
		assert_eq!(p.representation, Representation::default());

		p.change_representation(Representation::PowerBasis);
		assert_eq!(p.representation, Representation::PowerBasis);
		assert!(p.coefficients_shoup.is_none());
		let q = p.clone();

		p.change_representation(Representation::Ntt);
		assert_eq!(p.representation, Representation::Ntt);
		assert_ne!(p.coefficients, q.coefficients);
		assert!(p.coefficients_shoup.is_none());
		let q_ntt = p.clone();

		p.change_representation(Representation::NttShoup);
		assert_eq!(p.representation, Representation::NttShoup);
		assert_ne!(p.coefficients, q.coefficients);
		assert!(p.coefficients_shoup.is_some());
		let q_ntt_shoup = p.clone();

		p.change_representation(Representation::PowerBasis);
		assert_eq!(p, q);

		p.change_representation(Representation::NttShoup);
		assert_eq!(p, q_ntt_shoup);

		p.change_representation(Representation::Ntt);
		assert_eq!(p, q_ntt);

		p.change_representation(Representation::PowerBasis);
		assert_eq!(p, q);

		Ok(())
	}

	#[test]
	fn test_override_representation() -> Result<(), String> {
		let ctx = Arc::new(Context::new(MODULI, 8)?);

		let mut p = Poly::random(&ctx, Representation::PowerBasis);
		let q = p.clone();

		unsafe { p.override_representation(Representation::Ntt) }
		assert_eq!(p.representation, Representation::Ntt);
		assert_eq!(p.coefficients, q.coefficients);
		assert!(p.coefficients_shoup.is_none());

		unsafe { p.override_representation(Representation::NttShoup) }
		assert_eq!(p.representation, Representation::NttShoup);
		assert_eq!(p.coefficients, q.coefficients);
		assert!(p.coefficients_shoup.is_some());

		unsafe { p.override_representation(Representation::PowerBasis) }
		assert_eq!(p, q);

		unsafe { p.override_representation(Representation::NttShoup) }
		assert!(p.coefficients_shoup.is_some());

		unsafe { p.override_representation(Representation::Ntt) }
		assert!(p.coefficients_shoup.is_none());

		Ok(())
	}

	#[test]
	fn test_small() -> Result<(), String> {
		for modulus in MODULI {
			let ctx = Arc::new(Context::new(&[*modulus], 8)?);
			let q = Modulus::new(*modulus).unwrap();
			assert!(Poly::small(&ctx, Representation::PowerBasis, 0).is_err_and(
				|e| e.to_string() == "The variance should be an integer between 1 and 16"
			));
			assert!(
				Poly::small(&ctx, Representation::PowerBasis, 17).is_err_and(
					|e| e.to_string() == "The variance should be an integer between 1 and 16"
				)
			);
			for i in 1..=16 {
				let p = Poly::small(&ctx, Representation::PowerBasis, i)?;
				let coefficients = p.coefficients().to_slice().unwrap();
				let v = unsafe { q.center_vec_vt(coefficients) }
					.iter()
					.map(|ai| *ai as f64)
					.collect_vec();
				let s = test::stats::Summary::new(&v);
				assert!(s.min >= -2.0 * (i as f64));
				assert!(s.max <= 2.0 * (i as f64));
			}
		}

		// TODO: This is *way* too slow??
		// // Generate a very large polynomial to check the variance (here equal to 8).
		// let ctx = Arc::new(Context::new(&[4611686018326724609], 1 << 13)?);
		// let q = Modulus::new(4611686018326724609).unwrap();
		// let p = Poly::small(&ctx, Representation::PowerBasis, 8)?;
		// let coefficients = p.coefficients().to_slice().unwrap();
		// let v = unsafe { q.center_vec_vt(coefficients) }
		// 	.iter()
		// 	.map(|ai| *ai as f64)
		// 	.collect_vec();
		// let s = test::stats::Summary::new(&v);
		// assert!(s.min >= -16.0);
		// assert!(s.max <= 16.0);
		// assert_eq!(s.var.round(), 8.0);

		Ok(())
	}

	#[test]
	fn test_substitute() -> Result<(), String> {
		for modulus in MODULI {
			let ctx = Arc::new(Context::new(&[*modulus], 8)?);
			let p = Poly::random(&ctx, Representation::PowerBasis);
			let mut p_ntt = p.clone();
			p_ntt.change_representation(Representation::Ntt);
			let mut p_ntt_shoup = p.clone();
			p_ntt_shoup.change_representation(Representation::NttShoup);
			let p_coeffs = Vec::<u64>::from(&p);

			// Substitution by a multiple of 2 * degree should fail
			assert!(
				p.substitute(0).is_err()
					&& p_ntt.substitute(0).is_err()
					&& p_ntt_shoup.substitute(0).is_err()
			);
			assert!(
				p.substitute(16).is_err()
					&& p_ntt.substitute(16).is_err()
					&& p_ntt_shoup.substitute(16).is_err()
			);

			// Substitution by 1 leaves the polynomials unchanged
			assert_eq!(p, p.substitute(1)?);
			assert_eq!(p_ntt, p_ntt.substitute(1)?);
			assert_eq!(p_ntt_shoup, p_ntt_shoup.substitute(1)?);

			// Substitution by an even number is only possible in PowerBasis representation.
			assert!(p_ntt.substitute(2).is_err() && p_ntt_shoup.substitute(2).is_err());
			let mut v = vec![0u64; 8];
			for i in 0..8 {
				let vi = if 2 * i >= 8 && p_coeffs[i] > 0 {
					*modulus - p_coeffs[i]
				} else {
					p_coeffs[i]
				};
				v[(2 * i) % 8] = (v[(2 * i) % 8] + vi) % *modulus
			}
			assert_eq!(&Vec::<u64>::from(&p.substitute(2)?), &v);

			// In Ntt and NttShoup representations, the exponent must be odd
			let mut q = p.substitute(3)?;
			let mut v = vec![0u64; 8];
			for i in 0..8 {
				v[(3 * i) % 8] = if ((3 * i) / 8) & 1 == 1 && p_coeffs[i] > 0 {
					*modulus - p_coeffs[i]
				} else {
					p_coeffs[i]
				};
			}
			assert_eq!(&Vec::<u64>::from(&q), &v);

			let q_ntt = p_ntt.substitute(3)?;
			q.change_representation(Representation::Ntt);
			assert_eq!(q, q_ntt);

			let q_ntt_shoup = p_ntt_shoup.substitute(3)?;
			q.change_representation(Representation::NttShoup);
			assert_eq!(q, q_ntt_shoup);

			// 11 = 3^(-1) % 16
			assert_eq!(p, p.substitute(3)?.substitute(11)?);
			assert_eq!(p_ntt, p_ntt.substitute(3)?.substitute(11)?);
			assert_eq!(p_ntt_shoup, p_ntt_shoup.substitute(3)?.substitute(11)?);
		}

		let ctx = Arc::new(Context::new(MODULI, 8)?);
		let p = Poly::random(&ctx, Representation::PowerBasis);
		let mut p_ntt = p.clone();
		p_ntt.change_representation(Representation::Ntt);
		let mut p_ntt_shoup = p.clone();
		p_ntt_shoup.change_representation(Representation::NttShoup);

		assert_eq!(p, p.substitute(3)?.substitute(11)?);
		assert_eq!(p_ntt, p_ntt.substitute(3)?.substitute(11)?);
		assert_eq!(p_ntt_shoup, p_ntt_shoup.substitute(3)?.substitute(11)?);

		Ok(())
	}
}
