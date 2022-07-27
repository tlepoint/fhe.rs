#![warn(missing_docs, unused_imports)]

//! Polynomials in R_q[x] = (ZZ_q1 x ... x ZZ_qn)[x] where the qi's are prime moduli in zq.

use crate::zq::{ntt::NttOperator, Modulus};
use itertools::izip;
use ndarray::Array2;
use std::{
	ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
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

/// Conversions to create polynomials. We unfortunaly cannot use the `TryFrom` trait from std::convert
/// because we need to specify additional parameters, and if we try to redefine a `TryFrom` trait here,
/// we need to fully specify the trait when we use it because of the blanket implementation
/// <https://github.com/rust-lang/rust/issues/50133#issuecomment-488512355>.
pub trait TryConvertFrom<T>
where
	Self: Sized,
{
	/// The type of errors.
	type Error;

	/// Attempt to convert the `value` into a polynomial with a specific context and under a specific representation.
	fn try_convert_from(
		value: T,
		ctx: &Rc<Context>,
		representation: Representation,
	) -> Result<Self, Self::Error>;
}

impl TryConvertFrom<u64> for Poly {
	type Error = &'static str;

	fn try_convert_from(
		value: u64,
		ctx: &Rc<Context>,
		representation: Representation,
	) -> Result<Self, Self::Error> {
		match representation {
			Representation::PowerBasis => Poly::try_convert_from(&[value], ctx, representation),
			_ => {
				Err("Converting from constant values is only possible in PowerBasis representation")
			}
		}
	}
}

impl TryConvertFrom<Vec<u64>> for Poly {
	type Error = &'static str;

	fn try_convert_from(
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

impl TryConvertFrom<Array2<u64>> for Poly {
	type Error = &'static str;

	fn try_convert_from(
		a: Array2<u64>,
		ctx: &Rc<Context>,
		representation: Representation,
	) -> Result<Self, Self::Error> {
		if a.shape() != [ctx.q.len(), ctx.degree] {
			Err("The array of coefficient does not have the correct shape")
		} else {
			Ok(Self {
				ctx: ctx.clone(),
				representation,
				coefficients: a,
			})
		}
	}
}

impl TryConvertFrom<&[u64]> for Poly {
	type Error = &'static str;

	fn try_convert_from(
		v: &[u64],
		ctx: &Rc<Context>,
		representation: Representation,
	) -> Result<Self, Self::Error> {
		let v_clone: Vec<u64> = v.to_vec();
		Poly::try_convert_from(v_clone, ctx, representation)
	}
}

impl TryConvertFrom<&Vec<u64>> for Poly {
	type Error = &'static str;

	fn try_convert_from(
		v: &Vec<u64>,
		ctx: &Rc<Context>,
		representation: Representation,
	) -> Result<Self, Self::Error> {
		let v_clone: Vec<u64> = v.clone();
		Poly::try_convert_from(v_clone, ctx, representation)
	}
}

impl<const N: usize> TryConvertFrom<&[u64; N]> for Poly {
	type Error = &'static str;

	fn try_convert_from(
		v: &[u64; N],
		ctx: &Rc<Context>,
		representation: Representation,
	) -> Result<Self, Self::Error> {
		let v_clone: Vec<u64> = v.to_vec();
		Poly::try_convert_from(v_clone, ctx, representation)
	}
}

impl From<&Poly> for Vec<u64> {
	fn from(p: &Poly) -> Vec<u64> {
		p.coefficients.as_slice().unwrap().to_vec()
	}
}

impl AddAssign<&Poly> for Poly {
	fn add_assign(&mut self, p: &Poly) {
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

impl Neg for Poly {
	type Output = Poly;

	fn neg(self) -> Poly {
		let mut out = self;
		izip!(out.coefficients.outer_iter_mut(), &out.ctx.q)
			.for_each(|(mut v1, qi)| qi.neg_vec(v1.as_slice_mut().unwrap()));
		out
	}
}

impl Neg for &Poly {
	type Output = Poly;

	fn neg(self) -> Poly {
		let mut out = self.clone();
		izip!(out.coefficients.outer_iter_mut(), &out.ctx.q)
			.for_each(|(mut v1, qi)| qi.neg_vec(v1.as_slice_mut().unwrap()));
		out
	}
}

#[cfg(test)]
mod tests {
	use super::{Context, Poly, Representation, TryConvertFrom};
	use crate::zq::{ntt::supports_ntt, Modulus};
	use proptest::collection::vec as prop_vec;
	use proptest::prelude::{any, ProptestConfig};
	use std::rc::Rc;

	// Moduli to be used in tests.
	static MODULI: &[u64; 3] = &[1153, 4611686018326724609, 4611686018309947393];

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
			assert_eq!(Vec::from(&p), &[0, 0, 0, 0, 0, 0, 0, 0]);
			assert_eq!(Vec::from(&q), &[0, 0, 0, 0, 0, 0, 0, 0]);
		}

		let ctx = Rc::new(Context::new(MODULI, 8).unwrap());
		let p = Poly::zero(&ctx, Representation::PowerBasis);
		let q = Poly::zero(&ctx, Representation::Ntt);
		assert_ne!(p, q);
		assert_eq!(
			Vec::from(&p),
			&[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		);
		assert_eq!(
			Vec::from(&q),
			&[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		);
	}

	#[test]
	fn test_try_convert_from_slice_zero() {
		for modulus in MODULI {
			let ctx = Rc::new(Context::new(&[*modulus], 8).unwrap());
			let mut p = Poly::try_convert_from(&[], &ctx, Representation::PowerBasis);
			assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
			p = Poly::try_convert_from(&[], &ctx, Representation::Ntt);
			assert!(p.is_err());

			p = Poly::try_convert_from(&[0], &ctx, Representation::PowerBasis);
			assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
			p = Poly::try_convert_from(&[0], &ctx, Representation::Ntt);
			assert!(p.is_err());

			p = Poly::try_convert_from(&[0, 0, 0, 0, 0, 0, 0, 0], &ctx, Representation::PowerBasis);
			assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
			p = Poly::try_convert_from(&[0, 0, 0, 0, 0, 0, 0, 0], &ctx, Representation::Ntt);
			assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::Ntt)));

			p = Poly::try_convert_from(
				&[0, 0, 0, 0, 0, 0, 0, 0, 0],
				&ctx,
				Representation::PowerBasis,
			);
			assert!(p.is_err());
			p = Poly::try_convert_from(&[0, 0, 0, 0, 0, 0, 0, 0, 0], &ctx, Representation::Ntt);
			assert!(p.is_err());
		}

		let ctx = Rc::new(Context::new(MODULI, 8).unwrap());
		let mut p = Poly::try_convert_from(&[], &ctx, Representation::PowerBasis);
		assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
		p = Poly::try_convert_from(&[], &ctx, Representation::Ntt);
		assert!(p.is_err());

		p = Poly::try_convert_from(&[0], &ctx, Representation::PowerBasis);
		assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
		p = Poly::try_convert_from(&[0], &ctx, Representation::Ntt);
		assert!(p.is_err());

		p = Poly::try_convert_from(&[0, 0, 0, 0, 0, 0, 0, 0], &ctx, Representation::PowerBasis);
		assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
		p = Poly::try_convert_from(&[0, 0, 0, 0, 0, 0, 0, 0], &ctx, Representation::Ntt);
		assert!(p.is_err());

		p = Poly::try_convert_from(
			&[0, 0, 0, 0, 0, 0, 0, 0, 0],
			&ctx,
			Representation::PowerBasis,
		);
		assert!(p.is_err());
		p = Poly::try_convert_from(&[0, 0, 0, 0, 0, 0, 0, 0, 0], &ctx, Representation::Ntt);
		assert!(p.is_err());

		p = Poly::try_convert_from(
			&[
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			],
			&ctx,
			Representation::PowerBasis,
		);
		assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
		p = Poly::try_convert_from(
			&[
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			],
			&ctx,
			Representation::Ntt,
		);
		assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::Ntt)));
	}

	#[test]
	fn test_try_convert_from_vec_zero() {
		for modulus in MODULI {
			let ctx = Rc::new(Context::new(&[*modulus], 8).unwrap());
			let mut p = Poly::try_convert_from(vec![], &ctx, Representation::PowerBasis);
			assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
			p = Poly::try_convert_from(vec![], &ctx, Representation::Ntt);
			assert!(p.is_err());

			p = Poly::try_convert_from(vec![0], &ctx, Representation::PowerBasis);
			assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
			p = Poly::try_convert_from(vec![0], &ctx, Representation::Ntt);
			assert!(p.is_err());

			p = Poly::try_convert_from(
				vec![0, 0, 0, 0, 0, 0, 0, 0],
				&ctx,
				Representation::PowerBasis,
			);
			assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
			p = Poly::try_convert_from(vec![0, 0, 0, 0, 0, 0, 0, 0], &ctx, Representation::Ntt);
			assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::Ntt)));

			p = Poly::try_convert_from(
				vec![0, 0, 0, 0, 0, 0, 0, 0, 0],
				&ctx,
				Representation::PowerBasis,
			);
			assert!(p.is_err());
			p = Poly::try_convert_from(vec![0, 0, 0, 0, 0, 0, 0, 0, 0], &ctx, Representation::Ntt);
			assert!(p.is_err());
		}

		let ctx = Rc::new(Context::new(MODULI, 8).unwrap());
		let mut p = Poly::try_convert_from(vec![], &ctx, Representation::PowerBasis);
		assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
		p = Poly::try_convert_from(vec![], &ctx, Representation::Ntt);
		assert!(p.is_err());

		p = Poly::try_convert_from(vec![0], &ctx, Representation::PowerBasis);
		assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
		p = Poly::try_convert_from(vec![0], &ctx, Representation::Ntt);
		assert!(p.is_err());

		p = Poly::try_convert_from(
			vec![0, 0, 0, 0, 0, 0, 0, 0],
			&ctx,
			Representation::PowerBasis,
		);
		assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
		p = Poly::try_convert_from(vec![0, 0, 0, 0, 0, 0, 0, 0], &ctx, Representation::Ntt);
		assert!(p.is_err());

		p = Poly::try_convert_from(
			vec![0, 0, 0, 0, 0, 0, 0, 0, 0],
			&ctx,
			Representation::PowerBasis,
		);
		assert!(p.is_err());
		p = Poly::try_convert_from(vec![0, 0, 0, 0, 0, 0, 0, 0, 0], &ctx, Representation::Ntt);
		assert!(p.is_err());

		p = Poly::try_convert_from(
			vec![
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			],
			&ctx,
			Representation::PowerBasis,
		);
		assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
		p = Poly::try_convert_from(
			vec![
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			],
			&ctx,
			Representation::Ntt,
		);
		assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::Ntt)));
	}

	#[test]
	fn test_try_convert_from_u64_zero() {
		for modulus in MODULI {
			let ctx = Rc::new(Context::new(&[*modulus], 8).unwrap());
			let mut p = <Poly as TryConvertFrom<u64>>::try_convert_from(
				0,
				&ctx,
				Representation::PowerBasis,
			);
			assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
			p = <Poly as TryConvertFrom<u64>>::try_convert_from(0, &ctx, Representation::Ntt);
			assert!(p.is_err());
		}

		let ctx = Rc::new(Context::new(MODULI, 8).unwrap());
		let mut p =
			<Poly as TryConvertFrom<u64>>::try_convert_from(0, &ctx, Representation::PowerBasis);
		assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
		p = <Poly as TryConvertFrom<u64>>::try_convert_from(0, &ctx, Representation::Ntt);
		assert!(p.is_err());
	}

	proptest! {
		#![proptest_config(ProptestConfig::with_cases(64))]
		#[test]
		fn test_add(a in prop_vec(any::<u64>(), 8), b in prop_vec(any::<u64>(), 8)) {
			for modulus in MODULI {
				let ctx = Rc::new(Context::new(&[*modulus], 8).unwrap());
				let m = Modulus::new(*modulus).unwrap();
				let mut c = m.reduce_vec_new(&a);
				let d = m.reduce_vec_new(&b);

				let p = Poly::try_convert_from(&c, &ctx, Representation::PowerBasis).ok().unwrap();
				let q = Poly::try_convert_from(&d, &ctx, Representation::PowerBasis).ok().unwrap();
				let r = &p + &q;
				prop_assert_eq!(&r.representation, &Representation::PowerBasis);
				m.add_vec(&mut c, &d);
				prop_assert_eq!(&Vec::from(&r), &c);

				let p = Poly::try_convert_from(&c, &ctx, Representation::Ntt).ok().unwrap();
				let q = Poly::try_convert_from(&d, &ctx, Representation::Ntt).ok().unwrap();
				let r = &p + &q;
				prop_assert_eq!(&r.representation, &Representation::Ntt);
				m.add_vec(&mut c, &d);
				prop_assert_eq!(&Vec::from(&r), &c);
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
			let p = Poly::try_convert_from(&a, &ctx, Representation::PowerBasis).ok().unwrap();
			let q = Poly::try_convert_from(&b, &ctx, Representation::PowerBasis).ok().unwrap();
			let r = &p + &q;
			prop_assert_eq!(&r.representation, &Representation::PowerBasis);
			prop_assert_eq!(&Vec::from(&r), &reference);
		}

		#[test]
		fn test_sub(a in prop_vec(any::<u64>(), 8), b in prop_vec(any::<u64>(), 8)) {
			for modulus in MODULI {
				let ctx = Rc::new(Context::new(&[*modulus], 8).unwrap());
				let m = Modulus::new(*modulus).unwrap();
				let mut c = m.reduce_vec_new(&a);
				let d = m.reduce_vec_new(&b);

				let p = Poly::try_convert_from(&c, &ctx, Representation::PowerBasis).ok().unwrap();
				let q = Poly::try_convert_from(&d, &ctx, Representation::PowerBasis).ok().unwrap();
				let r = &p - &q;
				prop_assert_eq!(&r.representation, &Representation::PowerBasis);
				m.sub_vec(&mut c, &d);
				prop_assert_eq!(&Vec::from(&r), &c);

				let p = Poly::try_convert_from(&c, &ctx, Representation::Ntt).ok().unwrap();
				let q = Poly::try_convert_from(&d, &ctx, Representation::Ntt).ok().unwrap();
				let r = &p - &q;
				prop_assert_eq!(&r.representation, &Representation::Ntt);
				m.sub_vec(&mut c, &d);
				prop_assert_eq!(&Vec::from(&r), &c);
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
			let p = Poly::try_convert_from(&a, &ctx, Representation::PowerBasis).ok().unwrap();
			let q = Poly::try_convert_from(&b, &ctx, Representation::PowerBasis).ok().unwrap();
			let r = &p - &q;
			prop_assert_eq!(&r.representation, &Representation::PowerBasis);
			prop_assert_eq!(&Vec::from(&r), &reference);
		}

		#[test]
		fn test_mul(a in prop_vec(any::<u64>(), 8), b in prop_vec(any::<u64>(), 8), mut a2 in prop_vec(any::<u64>(), 24), mut b2 in prop_vec(any::<u64>(), 24)) {
			for modulus in MODULI {
				let ctx = Rc::new(Context::new(&[*modulus], 8).unwrap());
				let m = Modulus::new(*modulus).unwrap();
				let mut c = m.reduce_vec_new(&a);
				let d = m.reduce_vec_new(&b);

				let p = Poly::try_convert_from(&c, &ctx, Representation::Ntt).ok().unwrap();
				let q = Poly::try_convert_from(&d, &ctx, Representation::Ntt).ok().unwrap();
				let r = &p * &q;
				prop_assert_eq!(&r.representation, &Representation::Ntt);
				m.mul_vec(&mut c, &d);
				prop_assert_eq!(&Vec::from(&r), &c);
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
			let p = Poly::try_convert_from(&a2, &ctx, Representation::Ntt).ok().unwrap();
			let q = Poly::try_convert_from(&b2, &ctx, Representation::Ntt).ok().unwrap();
			let r = &p * &q;
			prop_assert_eq!(&r.representation, &Representation::Ntt);
			prop_assert_eq!(&Vec::from(&r), &reference);
		}

		#[test]
		fn test_neg(a in prop_vec(any::<u64>(), 8)) {
			for modulus in MODULI {
				let ctx = Rc::new(Context::new(&[*modulus], 8).unwrap());
				let m = Modulus::new(*modulus).unwrap();
				let mut c = m.reduce_vec_new(&a);

				let p = Poly::try_convert_from(&c, &ctx, Representation::PowerBasis).ok().unwrap();
				let r = -p;
				prop_assert_eq!(&r.representation, &Representation::PowerBasis);
				m.neg_vec(&mut c);
				prop_assert_eq!(&Vec::from(&r), &c);

				let p = Poly::try_convert_from(&c, &ctx, Representation::Ntt).ok().unwrap();
				let r = -p;
				prop_assert_eq!(&r.representation, &Representation::Ntt);
				m.neg_vec(&mut c);
				prop_assert_eq!(&Vec::from(&r), &c);
			}

			let mut reference = vec![];
			for modulus in MODULI {
				let m = Modulus::new(*modulus).unwrap();
				let mut c = m.reduce_vec_new(&a);
				m.neg_vec(&mut c);
				reference.extend(c)
			}
			let ctx = Rc::new(Context::new(MODULI, 8).unwrap());
			let p = Poly::try_convert_from(&a, &ctx, Representation::PowerBasis).ok().unwrap();
			let r = -&p;
			prop_assert_eq!(&r.representation, &Representation::PowerBasis);
			prop_assert_eq!(&Vec::from(&r), &reference);
			let r = -p;
			prop_assert_eq!(&r.representation, &Representation::PowerBasis);
			prop_assert_eq!(&Vec::from(&r), &reference);
		}
	}
}
