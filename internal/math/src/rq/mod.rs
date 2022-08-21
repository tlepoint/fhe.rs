#![warn(missing_docs, unused_imports)]

//! Polynomials in R_q\[x\] = (ZZ_q1 x ... x ZZ_qn)\[x\] where the qi's are
//! prime moduli in zq.

mod context;
mod convert;
mod ops;
mod serialize;

pub mod scaler;
pub mod traits;
pub use context::Context;
pub use ops::dot_product;

use crate::{Error, Result};
use itertools::{izip, Itertools};
use ndarray::{s, Array2, ArrayView2, Axis};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use std::sync::Arc;
use traits::TryConvertFrom;
use util::sample_vec_cbd;
use zeroize::Zeroize;

/// Possible representations of the underlying polynomial.
#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub enum Representation {
	/// This is the list of coefficients ci, such that the polynomial is c0 + c1
	/// * x + ... + c_(degree - 1) * x^(degree - 1)
	#[default]
	PowerBasis,
	/// This is the NTT representation of the PowerBasis representation.
	Ntt,
	/// This is a "Shoup" representation of the Ntt representation used for
	/// faster multiplication.
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
	has_lazy_coefficients: bool,
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
			has_lazy_coefficients: false,
		}
	}

	/// Enable variable time computations when this polynomial is involved.
	///
	/// # Safety
	///
	/// By default, this is marked as unsafe, but is usually safe when only
	/// public data is processed.
	pub unsafe fn allow_variable_time_computations(&mut self) {
		self.allow_variable_time_computations = true
	}

	/// Disable variable time computations when this polynomial is involved.
	pub fn disallow_variable_time_computations(&mut self) {
		self.allow_variable_time_computations = false
	}

	/// Current representation of the polynomial.
	/// TODO: To test
	pub const fn representation(&self) -> &Representation {
		&self.representation
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
	/// Prefer the `change_representation` function to safely modify the
	/// polynomial representation. If the `to` representation is NttShoup, the
	/// coefficients are still computed correctly to avoid being in an unstable
	/// state. Similarly, if we override a representation which was NttShoup, we
	/// zeroize the existing Shoup coefficients.
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

	/// Generate a small polynomial and convert into the specified
	/// representation.
	///
	/// Returns an error if the variance does not belong to [1, ..., 16].
	pub fn small(
		ctx: &Arc<Context>,
		representation: Representation,
		variance: usize,
	) -> Result<Self> {
		if !(1..=16).contains(&variance) {
			Err(Error::Default(
				"The variance should be an integer between 1 and 16".to_string(),
			))
		} else {
			let mut coeffs =
				sample_vec_cbd(ctx.degree, variance).map_err(|e| Error::Default(e.to_string()))?;
			let mut p = Poly::try_convert_from(
				coeffs.as_ref() as &[i64],
				ctx,
				false,
				Representation::PowerBasis,
			)?;
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
				.for_each(|(mut v, op)| unsafe { op.forward_vt(v.as_mut_ptr()) });
		} else {
			izip!(self.coefficients.outer_iter_mut(), &self.ctx.ops)
				.for_each(|(mut v, op)| op.forward(v.as_slice_mut().unwrap()));
		}
	}

	/// Computes the backward Ntt on the coefficients
	fn ntt_backward(&mut self) {
		if self.allow_variable_time_computations {
			izip!(self.coefficients.outer_iter_mut(), &self.ctx.ops)
				.for_each(|(mut v, op)| unsafe { op.backward_vt(v.as_mut_ptr()) });
		} else {
			izip!(self.coefficients.outer_iter_mut(), &self.ctx.ops)
				.for_each(|(mut v, op)| op.backward(v.as_slice_mut().unwrap()));
		}
	}

	/// Substitute x by x^i in a polynomial.
	/// In PowerBasis representation, i can be any integer that is not a
	/// multiple of 2 * degree. In Ntt and NttShoup representation, i can be any
	/// odd integer that is not a multiple of 2 * degree.
	pub fn substitute(&self, i: usize) -> Result<Poly> {
		let degree = self.ctx.degree as u32;
		let exponent = (i as u32) % (2 * degree);
		if exponent == 0 {
			return Err(Error::Default(
				"The exponent is a multiple of 2 * degree".to_string(),
			));
		}

		let mut q = Poly::zero(&self.ctx, self.representation.clone());
		if self.allow_variable_time_computations {
			unsafe { q.allow_variable_time_computations() }
		}
		let mask = (self.ctx.degree - 1) as u32;
		match self.representation {
			Representation::Ntt => {
				if exponent & 1 == 0 {
					return Err(Error::Default(
						"The exponent should be odd modulo 2 * degree".to_string(),
					));
				}
				let mut power = (exponent - 1) / 2;
				let power_bitrev = (0..degree)
					.map(|_| {
						let r = ((power & mask).reverse_bits() >> (degree.leading_zeros() + 1))
							as usize;
						power += exponent;
						r
					})
					.collect_vec();
				izip!(
					q.coefficients.outer_iter_mut(),
					self.coefficients.outer_iter()
				)
				.for_each(|(mut q_row, p_row)| {
					for (j, k) in izip!(&self.ctx.bitrev, &power_bitrev) {
						q_row[*j] = p_row[*k]
					}
				});
			}
			Representation::NttShoup => {
				if exponent & 1 == 0 {
					return Err(Error::Default(
						"The exponent should be odd modulo 2 * degree".to_string(),
					));
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
						if power & degree != 0 {
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

	/// Create a polynomial which can only be multiplied by a polynomial in
	/// NttShoup representation. All other operations may panic.
	///
	/// # Safety
	/// This operation also creates a polynomial that allows variable time
	/// operations.
	pub unsafe fn create_constant_ntt_polynomial_with_lazy_coefficients_and_variable_time(
		power_basis_coefficients: &[u64],
		ctx: &Arc<Context>,
	) -> Self {
		let mut coefficients = Array2::zeros((ctx.q.len(), ctx.degree));
		izip!(coefficients.outer_iter_mut(), &ctx.q, &ctx.ops).for_each(|(mut p, qi, op)| {
			p.as_slice_mut()
				.unwrap()
				.clone_from_slice(power_basis_coefficients);
			qi.reduce_vec(p.as_slice_mut().unwrap());
			op.forward_vt_lazy(p.as_mut_ptr());
		});
		Self {
			ctx: ctx.clone(),
			representation: Representation::Ntt,
			allow_variable_time_computations: true,
			coefficients,
			coefficients_shoup: None,
			has_lazy_coefficients: true,
		}
	}

	pub fn ctx(&self) -> &Arc<Context> {
		&self.ctx
	}

	/// Modulus switch down the polynomial by dividing and rounding each
	/// coefficient by the last modulus in the chain, then drops the last
	/// modulus, as described in Algorithm 2 of <https://eprint.iacr.org/2018/931.pdf>.
	///
	/// Returns an error if there is no next context or if the representation
	/// is not PowerBasis.
	pub fn mod_switch_down_next(&mut self) -> Result<()> {
		if self.ctx.next_context.is_none() {
			return Err(Error::NoMoreContext);
		}

		if self.representation != Representation::PowerBasis {
			return Err(Error::IncorrectRepresentation(
				self.representation.clone(),
				Representation::PowerBasis,
			));
		}

		// Unwrap the next_context.
		let next_context = self.ctx.next_context.as_ref().unwrap();

		let q_len = self.ctx.q.len();
		let q_last = self.ctx.q.last().unwrap();
		let q_last_div_2 = q_last.modulus() / 2;

		// Add (q_last - 1) / 2 to change from flooring to rounding
		let mut q_last_poly = self.coefficients.slice_mut(s![q_len - 1.., ..]);
		q_last_poly
			.iter_mut()
			.for_each(|coeff| *coeff = q_last.add(*coeff, q_last_div_2));

		let (mut q_new_polys, q_last_poly) =
			self.coefficients.view_mut().split_at(Axis(0), q_len - 1);

		izip!(
			q_new_polys.outer_iter_mut(),
			&self.ctx.q,
			&self.ctx.inv_last_qi_mod_qj,
			&self.ctx.inv_last_qi_mod_qj_shoup
		)
		.for_each(|(coeffs, qi, inv, inv_shoup)| {
			let q_last_div_2_mod_qi = qi.reduce(q_last_div_2);
			for (coeff, q_last_coeff) in izip!(coeffs, &q_last_poly) {
				// (x mod q_last - q_L/2) mod q_i
				let mut tmp = qi.reduce(*q_last_coeff);
				tmp = qi.sub(tmp, q_last_div_2_mod_qi);

				// ((x mod q_i) - (x mod q_last) + (q_L/2 mod q_i)) mod q_i
				// = (x - x mod q_last + q_L/2) mod q_i
				*coeff = qi.sub(*coeff, tmp);

				// q_last^{-1} * (x - x mod q_last) mod q_i
				*coeff = qi.mul_shoup(*coeff, *inv, *inv_shoup);
			}
		});

		// Remove the last row, and update the context.
		self.coefficients.remove_index(Axis(0), q_len - 1);
		self.ctx = next_context.clone();

		Ok(())
	}

	/// Modulo switch to a smaller context.
	///
	/// Returns an error if there is the provided context is not a child of the
	/// current context, or if the polynomial is not in PowerBasis
	/// representation.
	pub fn mod_switch_down_to(&mut self, context: &Arc<Context>) -> Result<()> {
		let niterations = self.ctx.niterations_to(context)?;
		for _ in 0..niterations {
			self.mod_switch_down_next()?;
		}
		assert_eq!(&self.ctx, context);
		Ok(())
	}

	pub fn multiply_inverse_power_of_x(&mut self, power: usize) -> Result<()> {
		let shift = ((self.ctx.degree << 1) - power) % (self.ctx.degree << 1);
		let mask = self.ctx.degree - 1;
		let original_coefficients = self.coefficients.clone();
		izip!(
			self.coefficients.outer_iter_mut(),
			original_coefficients.outer_iter(),
			&self.ctx.q
		)
		.for_each(|(mut coeffs, orig_coeffs, qi)| {
			for k in 0..self.ctx.degree {
				let index = shift + k;
				if index & self.ctx.degree == 0 {
					coeffs[index & mask] = orig_coeffs[k];
				} else {
					coeffs[index & mask] = qi.neg(orig_coeffs[k]);
				}
			}
		});
		Ok(())
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
	use crate::zq::Modulus;
	use itertools::Itertools;
	use num_bigint::BigUint;
	use num_traits::{One, Zero};
	use rand::{thread_rng, Rng, SeedableRng};
	use rand_chacha::ChaCha8Rng;
	use std::{error::Error, sync::Arc};

	// Moduli to be used in tests.
	const MODULI: &[u64; 5] = &[
		1153,
		4611686018326724609,
		4611686018309947393,
		4611686018232352769,
		4611686018171535361,
	];

	#[test]
	fn test_poly_zero() -> Result<(), Box<dyn Error>> {
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
		assert_eq!(Vec::<u64>::from(&p), [0; 8 * MODULI.len()]);
		assert_eq!(Vec::<u64>::from(&q), [0; 8 * MODULI.len()]);
		assert_eq!(Vec::<BigUint>::from(&p), reference);
		assert_eq!(Vec::<BigUint>::from(&q), reference);

		Ok(())
	}

	#[test]
	fn test_random() -> Result<(), Box<dyn Error>> {
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
	fn test_coefficients() -> Result<(), Box<dyn Error>> {
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
	fn test_modulus() -> Result<(), Box<dyn Error>> {
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
	fn test_allow_variable_time_computations() -> Result<(), Box<dyn Error>> {
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
	fn test_change_representation() -> Result<(), Box<dyn Error>> {
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
	fn test_override_representation() -> Result<(), Box<dyn Error>> {
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
	fn test_small() -> Result<(), Box<dyn Error>> {
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

		// Generate a very large polynomial to check the variance (here equal to 8).
		let ctx = Arc::new(Context::new(&[4611686018326724609], 1 << 18)?);
		let q = Modulus::new(4611686018326724609).unwrap();
		let p = Poly::small(&ctx, Representation::PowerBasis, 8)?;
		let coefficients = p.coefficients().to_slice().unwrap();
		let v = unsafe { q.center_vec_vt(coefficients) }
			.iter()
			.map(|ai| *ai as f64)
			.collect_vec();
		let s = test::stats::Summary::new(&v);
		assert!(s.min >= -16.0);
		assert!(s.max <= 16.0);
		assert_eq!(s.var.round(), 8.0);

		Ok(())
	}

	#[test]
	fn test_substitute() -> Result<(), Box<dyn Error>> {
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

	#[test]
	fn test_mod_switch_down_next() -> Result<(), Box<dyn Error>> {
		let ntests = 100;
		let ctx = Arc::new(Context::new(MODULI, 8)?);

		for _ in 0..ntests {
			// If the polynomial has incorrect representation, an error is returned
			assert!(Poly::random(&ctx, Representation::Ntt)
				.mod_switch_down_next()
				.is_err_and(|e| e
					== &crate::Error::IncorrectRepresentation(
						Representation::Ntt,
						Representation::PowerBasis
					)));

			// Otherwise, no error happens and the coefficients evolve as expected.
			let mut p = Poly::random(&ctx, Representation::PowerBasis);
			let mut reference = Vec::<BigUint>::from(&p);
			let mut current_ctx = ctx.clone();
			assert_eq!(p.ctx, current_ctx);
			while current_ctx.next_context.is_some() {
				let denominator = current_ctx.modulus().clone();
				current_ctx = current_ctx.next_context.as_ref().unwrap().clone();
				let numerator = current_ctx.modulus().clone();
				assert!(p.mod_switch_down_next().is_ok());
				assert_eq!(p.ctx, current_ctx);
				let p_biguint = Vec::<BigUint>::from(&p);
				assert_eq!(
					p_biguint,
					reference
						.iter()
						.map(|b| ((b * &numerator) + (&denominator >> 1)) / &denominator)
						.collect_vec()
				);
				reference = p_biguint.clone();
			}
		}
		Ok(())
	}

	#[test]
	fn test_mod_switch_down_to() -> Result<(), Box<dyn Error>> {
		let ntests = 100;
		let ctx1 = Arc::new(Context::new(MODULI, 8)?);
		let ctx2 = Arc::new(Context::new(&MODULI[..2], 8)?);

		for _ in 0..ntests {
			let mut p = Poly::random(&ctx1, Representation::PowerBasis);
			let reference = Vec::<BigUint>::from(&p);

			p.mod_switch_down_to(&ctx2)?;

			assert_eq!(p.ctx, ctx2);
			assert_eq!(
				Vec::<BigUint>::from(&p),
				reference
					.iter()
					.map(|b| ((b * ctx2.modulus()) + (ctx1.modulus() >> 1)) / ctx1.modulus())
					.collect_vec()
			);
		}

		Ok(())
	}
}
