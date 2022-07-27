#![warn(missing_docs, unused_imports)]

//! Polynomials in R_q[x] = (ZZ_q1 x ... x ZZ_qn)[x] where the qi's are prime moduli in zq.

use crate::zq::{ntt::NttOperator, Modulus};
use itertools::izip;
use ndarray::{Array2, Axis};
use std::{
	ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign},
	rc::Rc,
};

/// Struct that holds the context associated with elements in rq.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Context {
	q: Vec<Modulus>,
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
				let qi = qi.unwrap();
				let op = NttOperator::new(&qi, degree);
				op.as_ref()?;
				q.push(qi);
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
	representation: Representation,
	coefficients: Array2<u64>,
}

impl Poly {
	/// Creates a polynomial holding the constant 0.
	pub fn zero(ctx: &Rc<Context>, representation: Representation) -> Self {
		Self {
			ctx: ctx.clone(),
			representation,
			coefficients: Array2::zeros((ctx.q.len(), ctx.degree)),
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

impl From<&Poly> for Vec<u64> {
	fn from(p: &Poly) -> Vec<u64> {
		p.coefficients.as_slice().unwrap().to_vec()
	}
}

impl AddAssign<&Poly> for Poly {
	fn add_assign(&mut self, p: &Poly) {
		debug_assert_eq!(
			self.representation, p.representation,
			"Incompatible representations"
		);
		debug_assert_eq!(self.ctx, p.ctx, "Incompatible contexts");
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

impl Add<&Poly> for &Poly {
	type Output = Poly;
	fn add(self, p: &Poly) -> Poly {
		let mut q = self.clone();
		q += p;
		q
	}
}

impl SubAssign<&Poly> for Poly {
	fn sub_assign(&mut self, p: &Poly) {
		assert_eq!(
			self.representation, p.representation,
			"Incompatible representations"
		);
		debug_assert_eq!(self.ctx, p.ctx, "Incompatible contexts");
		izip!(
			self.coefficients.outer_iter_mut(),
			p.coefficients.outer_iter(),
			&self.ctx.q
		)
		.for_each(|(mut v1, v2, qi)| {
			qi.sub_vec(v1.as_slice_mut().unwrap(), v2.as_slice().unwrap())
		});
	}
}

impl Sub<&Poly> for &Poly {
	type Output = Poly;
	fn sub(self, p: &Poly) -> Poly {
		let mut q = self.clone();
		q -= p;
		q
	}
}

impl MulAssign<&Poly> for Poly {
	fn mul_assign(&mut self, p: &Poly) {
		assert_eq!(
			self.representation,
			Representation::Ntt,
			"Multiplication requires an Ntt representation."
		);
		assert_eq!(
			p.representation,
			Representation::Ntt,
			"Multiplication requires an Ntt representation."
		);
		debug_assert_eq!(self.ctx, p.ctx, "Incompatible contexts");
		izip!(
			self.coefficients.outer_iter_mut(),
			p.coefficients.outer_iter(),
			&self.ctx.q
		)
		.for_each(|(mut v1, v2, qi)| {
			qi.mul_vec(v1.as_slice_mut().unwrap(), v2.as_slice().unwrap())
		});
	}
}

impl Mul<&Poly> for &Poly {
	type Output = Poly;
	fn mul(self, p: &Poly) -> Poly {
		let mut q = self.clone();
		q *= p;
		q
	}
}

#[cfg(test)]
mod tests {
	use super::{Context, Poly, Representation, TryFrom};
	use crate::zq::{ntt::supports_ntt, Modulus};
	use proptest::collection::vec as prop_vec;
	use proptest::prelude::{any, ProptestConfig};
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

	proptest! {
		#![proptest_config(ProptestConfig::with_cases(100))]
		#[test]
		fn test_add(a in prop_vec(any::<u64>(), 8), b in prop_vec(any::<u64>(), 8)) {
			for modulus in MODULI {
				let ctx = Rc::new(Context::new(&[*modulus], 8).unwrap());
				let m = Modulus::new(*modulus).unwrap();
				let mut c = m.reduce_vec_new(&a);
				let d = m.reduce_vec_new(&b);

				let p = <Poly as TryFrom<&[u64]>>::try_from(&c, &ctx, Representation::PowerBasis).ok().unwrap();
				let q = <Poly as TryFrom<&[u64]>>::try_from(&d, &ctx, Representation::PowerBasis).ok().unwrap();
				let r = &p + &q;
				prop_assert_eq!(r.representation, Representation::PowerBasis);
				m.add_vec(&mut c, &d);
				prop_assert_eq!(r.coefficients.as_slice().unwrap(), &c);

				let p = <Poly as TryFrom<&[u64]>>::try_from(&c, &ctx, Representation::Ntt).ok().unwrap();
				let q = <Poly as TryFrom<&[u64]>>::try_from(&d, &ctx, Representation::Ntt).ok().unwrap();
				let r = &p + &q;
				prop_assert_eq!(r.representation, Representation::Ntt);
				m.add_vec(&mut c, &d);
				prop_assert_eq!(r.coefficients.as_slice().unwrap(), &c);
			}

			let mut reference = vec![];
			for modulus in MODULI {
				let m = Modulus::new(*modulus).unwrap();
				let mut c = m.reduce_vec_new(&a);
				let d = m.reduce_vec_new(&b);
				m.add_vec(&mut c, &d);
				reference.extend(c)
			}
			let ctx = Rc::new(Context::new(MODULI, 8).unwrap());
			let p = <Poly as TryFrom<&[u64]>>::try_from(&a, &ctx, Representation::PowerBasis).ok().unwrap();
			let q = <Poly as TryFrom<&[u64]>>::try_from(&b, &ctx, Representation::PowerBasis).ok().unwrap();
			let r = &p + &q;
			prop_assert_eq!(r.representation, Representation::PowerBasis);
			prop_assert_eq!(r.coefficients.as_slice().unwrap(), &reference);
		}

		#[test]
		fn test_sub(a in prop_vec(any::<u64>(), 8), b in prop_vec(any::<u64>(), 8)) {
			for modulus in MODULI {
				let ctx = Rc::new(Context::new(&[*modulus], 8).unwrap());
				let m = Modulus::new(*modulus).unwrap();
				let mut c = m.reduce_vec_new(&a);
				let d = m.reduce_vec_new(&b);

				let p = <Poly as TryFrom<&[u64]>>::try_from(&c, &ctx, Representation::PowerBasis).ok().unwrap();
				let q = <Poly as TryFrom<&[u64]>>::try_from(&d, &ctx, Representation::PowerBasis).ok().unwrap();
				let r = &p - &q;
				prop_assert_eq!(r.representation, Representation::PowerBasis);
				m.sub_vec(&mut c, &d);
				prop_assert_eq!(r.coefficients.as_slice().unwrap(), &c);

				let p = <Poly as TryFrom<&[u64]>>::try_from(&c, &ctx, Representation::Ntt).ok().unwrap();
				let q = <Poly as TryFrom<&[u64]>>::try_from(&d, &ctx, Representation::Ntt).ok().unwrap();
				let r = &p - &q;
				prop_assert_eq!(r.representation, Representation::Ntt);
				m.sub_vec(&mut c, &d);
				prop_assert_eq!(r.coefficients.as_slice().unwrap(), &c);
			}

			let mut reference = vec![];
			for modulus in MODULI {
				let m = Modulus::new(*modulus).unwrap();
				let mut c = m.reduce_vec_new(&a);
				let d = m.reduce_vec_new(&b);
				m.sub_vec(&mut c, &d);
				reference.extend(c)
			}
			let ctx = Rc::new(Context::new(MODULI, 8).unwrap());
			let p = <Poly as TryFrom<&[u64]>>::try_from(&a, &ctx, Representation::PowerBasis).ok().unwrap();
			let q = <Poly as TryFrom<&[u64]>>::try_from(&b, &ctx, Representation::PowerBasis).ok().unwrap();
			let r = &p - &q;
			prop_assert_eq!(r.representation, Representation::PowerBasis);
			prop_assert_eq!(r.coefficients.as_slice().unwrap(), &reference);
		}

		#[test]
		fn test_mul(a in prop_vec(any::<u64>(), 8), b in prop_vec(any::<u64>(), 8), mut a2 in prop_vec(any::<u64>(), 24), mut b2 in prop_vec(any::<u64>(), 24)) {
			for modulus in MODULI {
				let ctx = Rc::new(Context::new(&[*modulus], 8).unwrap());
				let m = Modulus::new(*modulus).unwrap();
				let mut c = m.reduce_vec_new(&a);
				let d = m.reduce_vec_new(&b);

				let p = <Poly as TryFrom<&[u64]>>::try_from(&c, &ctx, Representation::Ntt).ok().unwrap();
				let q = <Poly as TryFrom<&[u64]>>::try_from(&d, &ctx, Representation::Ntt).ok().unwrap();
				let r = &p * &q;
				prop_assert_eq!(r.representation, Representation::Ntt);
				m.mul_vec(&mut c, &d);
				prop_assert_eq!(r.coefficients.as_slice().unwrap(), &c);
			}

			let mut reference = vec![];
			for i in 0..MODULI.len() {
				let m = Modulus::new(MODULI[i]).unwrap();
				m.reduce_vec(&mut a2[i * 8..(i+1)*8]);
				m.reduce_vec(&mut b2[i * 8..(i+1)*8]);
				let mut a3 = a2[i * 8..(i+1)*8].to_vec();
				m.mul_vec(&mut a3, &b2[i * 8..(i+1)*8]);
				reference.extend(a3)
			}
			let ctx = Rc::new(Context::new(MODULI, 8).unwrap());
			let p = <Poly as TryFrom<&[u64]>>::try_from(&a2, &ctx, Representation::Ntt).ok().unwrap();
			let q = <Poly as TryFrom<&[u64]>>::try_from(&b2, &ctx, Representation::Ntt).ok().unwrap();
			let r = &p * &q;
			prop_assert_eq!(r.representation, Representation::Ntt);
			prop_assert_eq!(r.coefficients.as_slice().unwrap(), &reference);
		}

	}
}
