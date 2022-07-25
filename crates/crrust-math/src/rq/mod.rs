#![warn(missing_docs, unused_imports)]

//! Polynomials in R_q[x] = (ZZ_q1 x ... x ZZ_qn)[x] where the qi's are prime moduli in zq.

use crate::zq::{ntt::NttOperator, Modulus};
use itertools::izip;
use ndarray::Array2;
use std::{
	ops::{Add, AddAssign},
	rc::Rc,
};

/// Struct that holds the context associated with elements in rq.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Context {
	q: Vec<Rc<Modulus>>,
	ops: Vec<NttOperator>,
	degree: usize,
}

impl Context {
	/// Creates a context from a list of moduli and a polynomial degree.
	///
	/// Returns None if the moduli are not primes less than 62 bits which supports the NTT of size `degree`.
	pub fn new(moduli: &[u64], degree: usize) -> Option<Self> {
		if moduli.is_empty() || !degree.is_power_of_two() || degree < 8 {
			None
		} else {
			let mut q = Vec::with_capacity(moduli.len());
			let mut ops = Vec::with_capacity(moduli.len());
			for modulus in moduli {
				let qi = Modulus::new(*modulus);
				qi.as_ref()?;
				let rc = Rc::new(qi.unwrap());
				let op = NttOperator::new(&rc, degree);
				op.as_ref()?;
				q.push(rc);
				ops.push(op.unwrap());
			}

			Some(Self { q, ops, degree })
		}
	}
}

/// Possible representations of the underlying polynomial.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Representation {
	/// This is the list of coefficients ci, such that the polynomial is c0 + c1 * x + ... + c_(degree - 1) * x^(degree - 1)
	PowerBasis,
	/// This is the NTT representation of the PowerBasis representation.
	Ntt,
}

/// Struct that holds a polynomial for a specific context.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Poly {
	ctx: Rc<Context>,
	coefficients: Array2<u64>,
	representation: Representation,
}

impl Poly {
	/// Creates a polynomial holding the constant 0.
	pub fn zero(ctx: &Rc<Context>, representation: Representation) -> Self {
		Self {
			ctx: ctx.clone(),
			coefficients: Array2::zeros((ctx.q.len(), ctx.degree)),
			representation,
		}
	}

	/// Change the representation of the underlying polynomial.
	///
	/// Panics if the change of representation is illegal.
	pub fn change_representation(&mut self, to: Representation) {
		if self.representation == Representation::PowerBasis && to == Representation::Ntt {
			izip!(self.coefficients.outer_iter_mut(), &self.ctx.ops)
				.for_each(|(mut v, op)| op.forward(v.as_slice_mut().unwrap()));
		} else if self.representation == Representation::Ntt && to == Representation::PowerBasis {
			izip!(self.coefficients.outer_iter_mut(), &self.ctx.ops)
				.for_each(|(mut v, op)| op.backward(v.as_slice_mut().unwrap()));
		} else {
			panic!(
				"Invalid change of representation from {:?} to {:?}",
				self.representation, to
			)
		}
		self.representation = to;
	}
}

// Implement our own version of TryFrom.

pub trait TryFrom<T>
where
	Self: Sized,
{
	type Error;

	fn try_from(
		value: T,
		ctx: &Rc<Context>,
		representation: Representation,
	) -> Result<Self, Self::Error>;
}

impl TryFrom<u64> for Poly {
	type Error = &'static str;

	fn try_from(
		value: u64,
		ctx: &Rc<Context>,
		representation: Representation,
	) -> Result<Self, Self::Error> {
		match representation {
			Representation::PowerBasis => {
				<Poly as crate::rq::TryFrom<&[u64]>>::try_from(&[value], ctx, representation)
			}
			_ => {
				Err("Converting from constant values is only possible in PowerBasis representation")
			}
		}
	}
}

impl TryFrom<Vec<u64>> for Poly {
	type Error = &'static str;

	fn try_from(
		v: Vec<u64>,
		ctx: &Rc<Context>,
		representation: Representation,
	) -> Result<Self, Self::Error> {
		match representation {
			Representation::Ntt => {
				if let Ok(coefficients) = Array2::from_shape_vec((ctx.q.len(), ctx.degree), v) {
					Ok(Self {
						ctx: ctx.clone(),
						representation,
						coefficients,
					})
				} else {
					Err("In Ntt representation, all coefficients must be specified")
				}
			}
			Representation::PowerBasis => {
				if v.len() == ctx.q.len() * ctx.degree {
					let coefficients =
						Array2::from_shape_vec((ctx.q.len(), ctx.degree), v).unwrap();
					Ok(Self {
						ctx: ctx.clone(),
						representation,
						coefficients,
					})
				} else if v.len() <= ctx.degree {
					let mut out = Self::zero(ctx, representation);
					izip!(out.coefficients.outer_iter_mut(), &ctx.q).for_each(|(mut w, qi)| {
						let wi = w.as_slice_mut().unwrap();
						wi[..v.len()].copy_from_slice(&v);
						qi.reduce_vec(wi);
					});
					Ok(out)
				} else {
					Err("In PowerBasis representation, either all coefficients must be specified, or only coefficients up to the degree")
				}
			}
		}
	}
}

impl TryFrom<&[u64]> for Poly {
	type Error = &'static str;

	fn try_from(
		v: &[u64],
		ctx: &Rc<Context>,
		representation: Representation,
	) -> Result<Self, Self::Error> {
		let v_clone: Vec<u64> = v.to_vec();
		<Poly as crate::rq::TryFrom<Vec<u64>>>::try_from(v_clone, ctx, representation)
	}
}

impl AddAssign<&Poly> for Poly {
	fn add_assign(&mut self, p: &Poly) {
		assert_eq!(
			self.representation, p.representation,
			"Incompatible representations"
		);
		assert_eq!(self.ctx, p.ctx, "Incompatible contexts");
		izip!(
			self.coefficients.outer_iter_mut(),
			p.coefficients.outer_iter(),
			&self.ctx.q
		)
		.for_each(|(mut v1, v2, qi)| {
			qi.add_vec(v1.as_slice_mut().unwrap(), v2.as_slice().unwrap())
		});
	}
}

impl Add<&Poly> for Poly {
	type Output = Self;
	fn add(self, p: &Poly) -> Self {
		let mut q = self;
		q += p;
		q
	}
}

#[cfg(test)]
mod tests {
	use super::{Context, Poly, Representation, TryFrom};
	use crate::zq::ntt::supports_ntt;
	use std::rc::Rc;

	// Moduli to be used in tests.
	static MODULI: &[u64; 3] = &[1153, 4611686018326724609, 4611686018309947393];

	// Strategy for short and long polynomials.

	#[test]
	fn test_context_constructor() {
		for modulus in MODULI {
			assert!(Context::new(&[*modulus], 8).is_some());
			if supports_ntt(*modulus, 128) {
				assert!(Context::new(&[*modulus], 128).is_some());
			} else {
				assert!(Context::new(&[*modulus], 128).is_none());
			}
		}

		assert!(Context::new(MODULI, 8).is_some());
		assert!(Context::new(MODULI, 128).is_none());
	}

	#[test]
	fn test_poly_zero() {
		for modulus in MODULI {
			let ctx = Rc::new(Context::new(&[*modulus], 8).unwrap());
			let p = Poly::zero(&ctx, Representation::PowerBasis);
			let q = Poly::zero(&ctx, Representation::Ntt);
			assert_ne!(p, q);
			assert_eq!(
				p.coefficients.as_slice().unwrap(),
				&[0, 0, 0, 0, 0, 0, 0, 0]
			);
			assert_eq!(
				q.coefficients.as_slice().unwrap(),
				&[0, 0, 0, 0, 0, 0, 0, 0]
			);
		}

		let ctx = Rc::new(Context::new(MODULI, 8).unwrap());
		let p = Poly::zero(&ctx, Representation::PowerBasis);
		let q = Poly::zero(&ctx, Representation::Ntt);
		assert_ne!(p, q);
		assert_eq!(
			p.coefficients.as_slice().unwrap(),
			&[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		);
		assert_eq!(
			q.coefficients.as_slice().unwrap(),
			&[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		);
	}

	#[test]
	fn test_try_from_slice_zero() {
		for modulus in MODULI {
			let ctx = Rc::new(Context::new(&[*modulus], 8).unwrap());
			let mut p = <Poly as TryFrom<&[u64]>>::try_from(&[], &ctx, Representation::PowerBasis);
			assert_eq!(
				p.ok().unwrap(),
				Poly::zero(&ctx, Representation::PowerBasis)
			);
			p = <Poly as TryFrom<&[u64]>>::try_from(&[], &ctx, Representation::Ntt);
			assert!(p.is_err());

			p = <Poly as TryFrom<&[u64]>>::try_from(&[0], &ctx, Representation::PowerBasis);
			assert_eq!(
				p.ok().unwrap(),
				Poly::zero(&ctx, Representation::PowerBasis)
			);
			p = <Poly as TryFrom<&[u64]>>::try_from(&[0], &ctx, Representation::Ntt);
			assert!(p.is_err());

			p = <Poly as TryFrom<&[u64]>>::try_from(
				&[0, 0, 0, 0, 0, 0, 0, 0],
				&ctx,
				Representation::PowerBasis,
			);
			assert_eq!(
				p.ok().unwrap(),
				Poly::zero(&ctx, Representation::PowerBasis)
			);
			p = <Poly as TryFrom<&[u64]>>::try_from(
				&[0, 0, 0, 0, 0, 0, 0, 0],
				&ctx,
				Representation::Ntt,
			);
			assert_eq!(p.ok().unwrap(), Poly::zero(&ctx, Representation::Ntt));

			p = <Poly as TryFrom<&[u64]>>::try_from(
				&[0, 0, 0, 0, 0, 0, 0, 0, 0],
				&ctx,
				Representation::PowerBasis,
			);
			assert!(p.is_err());
			p = <Poly as TryFrom<&[u64]>>::try_from(
				&[0, 0, 0, 0, 0, 0, 0, 0, 0],
				&ctx,
				Representation::Ntt,
			);
			assert!(p.is_err());
		}

		let ctx = Rc::new(Context::new(MODULI, 8).unwrap());
		let mut p = <Poly as TryFrom<&[u64]>>::try_from(&[], &ctx, Representation::PowerBasis);
		assert_eq!(
			p.ok().unwrap(),
			Poly::zero(&ctx, Representation::PowerBasis)
		);
		p = <Poly as TryFrom<&[u64]>>::try_from(&[], &ctx, Representation::Ntt);
		assert!(p.is_err());

		p = <Poly as TryFrom<&[u64]>>::try_from(&[0], &ctx, Representation::PowerBasis);
		assert_eq!(
			p.ok().unwrap(),
			Poly::zero(&ctx, Representation::PowerBasis)
		);
		p = <Poly as TryFrom<&[u64]>>::try_from(&[0], &ctx, Representation::Ntt);
		assert!(p.is_err());

		p = <Poly as TryFrom<&[u64]>>::try_from(
			&[0, 0, 0, 0, 0, 0, 0, 0],
			&ctx,
			Representation::PowerBasis,
		);
		assert_eq!(
			p.ok().unwrap(),
			Poly::zero(&ctx, Representation::PowerBasis)
		);
		p = <Poly as TryFrom<&[u64]>>::try_from(
			&[0, 0, 0, 0, 0, 0, 0, 0],
			&ctx,
			Representation::Ntt,
		);
		assert!(p.is_err());

		p = <Poly as TryFrom<&[u64]>>::try_from(
			&[0, 0, 0, 0, 0, 0, 0, 0, 0],
			&ctx,
			Representation::PowerBasis,
		);
		assert!(p.is_err());
		p = <Poly as TryFrom<&[u64]>>::try_from(
			&[0, 0, 0, 0, 0, 0, 0, 0, 0],
			&ctx,
			Representation::Ntt,
		);
		assert!(p.is_err());

		p = <Poly as TryFrom<&[u64]>>::try_from(
			&[
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			],
			&ctx,
			Representation::PowerBasis,
		);
		assert_eq!(
			p.ok().unwrap(),
			Poly::zero(&ctx, Representation::PowerBasis)
		);
		p = <Poly as TryFrom<&[u64]>>::try_from(
			&[
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			],
			&ctx,
			Representation::Ntt,
		);
		assert_eq!(p.ok().unwrap(), Poly::zero(&ctx, Representation::Ntt));
	}

	#[test]
	fn test_try_from_vec_zero() {
		for modulus in MODULI {
			let ctx = Rc::new(Context::new(&[*modulus], 8).unwrap());
			let mut p =
				<Poly as TryFrom<Vec<u64>>>::try_from(vec![], &ctx, Representation::PowerBasis);
			assert_eq!(
				p.ok().unwrap(),
				Poly::zero(&ctx, Representation::PowerBasis)
			);
			p = <Poly as TryFrom<Vec<u64>>>::try_from(vec![], &ctx, Representation::Ntt);
			assert!(p.is_err());

			p = <Poly as TryFrom<Vec<u64>>>::try_from(vec![0], &ctx, Representation::PowerBasis);
			assert_eq!(
				p.ok().unwrap(),
				Poly::zero(&ctx, Representation::PowerBasis)
			);
			p = <Poly as TryFrom<Vec<u64>>>::try_from(vec![0], &ctx, Representation::Ntt);
			assert!(p.is_err());

			p = <Poly as TryFrom<Vec<u64>>>::try_from(
				vec![0, 0, 0, 0, 0, 0, 0, 0],
				&ctx,
				Representation::PowerBasis,
			);
			assert_eq!(
				p.ok().unwrap(),
				Poly::zero(&ctx, Representation::PowerBasis)
			);
			p = <Poly as TryFrom<Vec<u64>>>::try_from(
				vec![0, 0, 0, 0, 0, 0, 0, 0],
				&ctx,
				Representation::Ntt,
			);
			assert_eq!(p.ok().unwrap(), Poly::zero(&ctx, Representation::Ntt));

			p = <Poly as TryFrom<Vec<u64>>>::try_from(
				vec![0, 0, 0, 0, 0, 0, 0, 0, 0],
				&ctx,
				Representation::PowerBasis,
			);
			assert!(p.is_err());
			p = <Poly as TryFrom<Vec<u64>>>::try_from(
				vec![0, 0, 0, 0, 0, 0, 0, 0, 0],
				&ctx,
				Representation::Ntt,
			);
			assert!(p.is_err());
		}

		let ctx = Rc::new(Context::new(MODULI, 8).unwrap());
		let mut p = <Poly as TryFrom<Vec<u64>>>::try_from(vec![], &ctx, Representation::PowerBasis);
		assert_eq!(
			p.ok().unwrap(),
			Poly::zero(&ctx, Representation::PowerBasis)
		);
		p = <Poly as TryFrom<Vec<u64>>>::try_from(vec![], &ctx, Representation::Ntt);
		assert!(p.is_err());

		p = <Poly as TryFrom<Vec<u64>>>::try_from(vec![0], &ctx, Representation::PowerBasis);
		assert_eq!(
			p.ok().unwrap(),
			Poly::zero(&ctx, Representation::PowerBasis)
		);
		p = <Poly as TryFrom<Vec<u64>>>::try_from(vec![0], &ctx, Representation::Ntt);
		assert!(p.is_err());

		p = <Poly as TryFrom<Vec<u64>>>::try_from(
			vec![0, 0, 0, 0, 0, 0, 0, 0],
			&ctx,
			Representation::PowerBasis,
		);
		assert_eq!(
			p.ok().unwrap(),
			Poly::zero(&ctx, Representation::PowerBasis)
		);
		p = <Poly as TryFrom<Vec<u64>>>::try_from(
			vec![0, 0, 0, 0, 0, 0, 0, 0],
			&ctx,
			Representation::Ntt,
		);
		assert!(p.is_err());

		p = <Poly as TryFrom<Vec<u64>>>::try_from(
			vec![0, 0, 0, 0, 0, 0, 0, 0, 0],
			&ctx,
			Representation::PowerBasis,
		);
		assert!(p.is_err());
		p = <Poly as TryFrom<Vec<u64>>>::try_from(
			vec![0, 0, 0, 0, 0, 0, 0, 0, 0],
			&ctx,
			Representation::Ntt,
		);
		assert!(p.is_err());

		p = <Poly as TryFrom<Vec<u64>>>::try_from(
			vec![
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			],
			&ctx,
			Representation::PowerBasis,
		);
		assert_eq!(
			p.ok().unwrap(),
			Poly::zero(&ctx, Representation::PowerBasis)
		);
		p = <Poly as TryFrom<Vec<u64>>>::try_from(
			vec![
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			],
			&ctx,
			Representation::Ntt,
		);
		assert_eq!(p.ok().unwrap(), Poly::zero(&ctx, Representation::Ntt));
	}

	#[test]
	fn test_try_from_u64_zero() {
		for modulus in MODULI {
			let ctx = Rc::new(Context::new(&[*modulus], 8).unwrap());
			let mut p = <Poly as TryFrom<u64>>::try_from(0, &ctx, Representation::PowerBasis);
			assert_eq!(
				p.ok().unwrap(),
				Poly::zero(&ctx, Representation::PowerBasis)
			);
			p = <Poly as TryFrom<u64>>::try_from(0, &ctx, Representation::Ntt);
			assert!(p.is_err());
		}

		let ctx = Rc::new(Context::new(MODULI, 8).unwrap());
		let mut p = <Poly as TryFrom<u64>>::try_from(0, &ctx, Representation::PowerBasis);
		assert_eq!(
			p.ok().unwrap(),
			Poly::zero(&ctx, Representation::PowerBasis)
		);
		p = <Poly as TryFrom<u64>>::try_from(0, &ctx, Representation::Ntt);
		assert!(p.is_err());
	}
}
