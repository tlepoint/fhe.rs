//! Implementation of operations over polynomials.

use super::{traits::TryConvertFrom, Poly, Representation};
use itertools::izip;
use num_bigint::BigUint;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use zeroize::Zeroize;

impl AddAssign<&Poly> for Poly {
	fn add_assign(&mut self, p: &Poly) {
		assert_ne!(
			self.representation,
			Representation::NttShoup,
			"Cannot add to a polynomial in NttShoup representation"
		);
		assert_eq!(
			self.representation, p.representation,
			"Incompatible representations"
		);
		debug_assert_eq!(self.ctx, p.ctx, "Incompatible contexts");
		if self.allow_variable_time_computations || p.allow_variable_time_computations {
			izip!(
				self.coefficients.outer_iter_mut(),
				p.coefficients.outer_iter(),
				&self.ctx.q
			)
			.for_each(|(mut v1, v2, qi)| unsafe {
				qi.add_vec_vt(v1.as_slice_mut().unwrap(), v2.as_slice().unwrap())
			});
			self.allow_variable_time_computations = true
		} else {
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
		assert_ne!(
			self.representation,
			Representation::NttShoup,
			"Cannot subtract from a polynomial in NttShoup representation"
		);
		assert_eq!(
			self.representation, p.representation,
			"Incompatible representations"
		);
		debug_assert_eq!(self.ctx, p.ctx, "Incompatible contexts");
		if self.allow_variable_time_computations || p.allow_variable_time_computations {
			izip!(
				self.coefficients.outer_iter_mut(),
				p.coefficients.outer_iter(),
				&self.ctx.q
			)
			.for_each(|(mut v1, v2, qi)| unsafe {
				qi.sub_vec_vt(v1.as_slice_mut().unwrap(), v2.as_slice().unwrap())
			});
			self.allow_variable_time_computations = true
		} else {
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
		assert_ne!(
			self.representation,
			Representation::NttShoup,
			"Cannot multiply to a polynomial in NttShoup representation"
		);
		assert_eq!(
			self.representation,
			Representation::Ntt,
			"Multiplication requires an Ntt representation."
		);
		debug_assert_eq!(self.ctx, p.ctx, "Incompatible contexts");

		match p.representation {
			Representation::Ntt => {
				if self.allow_variable_time_computations || p.allow_variable_time_computations {
					unsafe {
						izip!(
							self.coefficients.outer_iter_mut(),
							p.coefficients.outer_iter(),
							&self.ctx.q
						)
						.for_each(|(mut v1, v2, qi)| {
							qi.mul_vec_vt(v1.as_slice_mut().unwrap(), v2.as_slice().unwrap());
						});
					}
					self.allow_variable_time_computations = true
				} else {
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
			Representation::NttShoup => {
				if self.allow_variable_time_computations || p.allow_variable_time_computations {
					izip!(
						self.coefficients.outer_iter_mut(),
						p.coefficients.outer_iter(),
						p.coefficients_shoup.as_ref().unwrap().outer_iter(),
						&self.ctx.q
					)
					.for_each(|(mut v1, v2, v2_shoup, qi)| unsafe {
						qi.mul_shoup_vec_vt(
							v1.as_slice_mut().unwrap(),
							v2.as_slice().unwrap(),
							v2_shoup.as_slice().unwrap(),
						)
					});
					self.allow_variable_time_computations = true
				} else {
					izip!(
						self.coefficients.outer_iter_mut(),
						p.coefficients.outer_iter(),
						p.coefficients_shoup.as_ref().unwrap().outer_iter(),
						&self.ctx.q
					)
					.for_each(|(mut v1, v2, v2_shoup, qi)| {
						qi.mul_shoup_vec(
							v1.as_slice_mut().unwrap(),
							v2.as_slice().unwrap(),
							v2_shoup.as_slice().unwrap(),
						)
					});
				}
			}
			_ => {
				panic!("Multiplication requires a multipliand in Ntt or NttShoup representation.")
			}
		}
	}
}

impl MulAssign<&BigUint> for Poly {
	fn mul_assign(&mut self, p: &BigUint) {
		let v: Vec<BigUint> = vec![p.clone()];
		let mut q = Poly::try_convert_from(
			v.as_ref() as &[BigUint],
			&self.ctx,
			self.representation.clone(),
		)
		.unwrap();
		q.change_representation(Representation::Ntt);
		if self.allow_variable_time_computations {
			unsafe {
				izip!(
					self.coefficients.outer_iter_mut(),
					q.coefficients.outer_iter(),
					&self.ctx.q
				)
				.for_each(|(mut v1, v2, qi)| {
					qi.mul_vec_vt(v1.as_slice_mut().unwrap(), v2.as_slice().unwrap())
				});
			}
		} else {
			izip!(
				self.coefficients.outer_iter_mut(),
				q.coefficients.outer_iter(),
				&self.ctx.q
			)
			.for_each(|(mut v1, v2, qi)| {
				qi.mul_vec(v1.as_slice_mut().unwrap(), v2.as_slice().unwrap())
			});
		}
	}
}

impl Mul<&Poly> for &Poly {
	type Output = Poly;
	fn mul(self, p: &Poly) -> Poly {
		match self.representation {
			Representation::NttShoup => {
				// TODO: To test, and do the same thing for add, sub, and neg
				let mut q = p.clone();
				if q.representation == Representation::NttShoup {
					q.coefficients_shoup
						.as_mut()
						.unwrap()
						.as_slice_mut()
						.unwrap()
						.zeroize();
					unsafe { q.override_representation(Representation::Ntt) }
				}
				q *= self;
				q
			}
			_ => {
				let mut q = self.clone();
				q *= p;
				q
			}
		}
	}
}

impl Mul<&BigUint> for &Poly {
	type Output = Poly;
	fn mul(self, p: &BigUint) -> Poly {
		let mut q = self.clone();
		q *= p;
		q
	}
}

impl Mul<&Poly> for &BigUint {
	type Output = Poly;
	fn mul(self, p: &Poly) -> Poly {
		p * self
	}
}

impl Neg for &Poly {
	type Output = Poly;

	fn neg(self) -> Poly {
		let mut out = self.clone();
		if self.allow_variable_time_computations {
			izip!(out.coefficients.outer_iter_mut(), &out.ctx.q)
				.for_each(|(mut v1, qi)| unsafe { qi.neg_vec_vt(v1.as_slice_mut().unwrap()) });
		} else {
			izip!(out.coefficients.outer_iter_mut(), &out.ctx.q)
				.for_each(|(mut v1, qi)| qi.neg_vec(v1.as_slice_mut().unwrap()));
		}
		out
	}
}

#[cfg(test)]
mod tests {
	use crate::{
		rq::{Context, Poly, Representation},
		zq::Modulus,
	};
	use std::rc::Rc;

	static MODULI: &[u64; 3] = &[1153, 4611686018326724609, 4611686018309947393];

	#[test]
	fn test_add() -> Result<(), String> {
		for _ in 0..100 {
			for modulus in MODULI {
				let ctx = Rc::new(Context::new(&[*modulus], 8)?);
				let m = Modulus::new(*modulus).unwrap();

				let p = Poly::random(&ctx, Representation::PowerBasis);
				let q = Poly::random(&ctx, Representation::PowerBasis);
				let r = &p + &q;
				assert_eq!(r.representation, Representation::PowerBasis);
				let mut a = Vec::<u64>::from(&p);
				m.add_vec(&mut a, &Vec::<u64>::from(&q));
				assert_eq!(Vec::<u64>::from(&r), a);

				let p = Poly::random(&ctx, Representation::Ntt);
				let q = Poly::random(&ctx, Representation::Ntt);
				let r = &p + &q;
				assert_eq!(r.representation, Representation::Ntt);
				let mut a = Vec::<u64>::from(&p);
				m.add_vec(&mut a, &Vec::<u64>::from(&q));
				assert_eq!(Vec::<u64>::from(&r), a);
			}

			let ctx = Rc::new(Context::new(MODULI, 8)?);
			let p = Poly::random(&ctx, Representation::PowerBasis);
			let q = Poly::random(&ctx, Representation::PowerBasis);
			let mut a = Vec::<u64>::from(&p);
			let b = Vec::<u64>::from(&q);
			for i in 0..MODULI.len() {
				let m = Modulus::new(MODULI[i]).unwrap();
				m.add_vec(&mut a[i * 8..(i + 1) * 8], &b[i * 8..(i + 1) * 8])
			}
			let r = &p + &q;
			assert_eq!(r.representation, Representation::PowerBasis);
			assert_eq!(Vec::<u64>::from(&r), a);
		}
		Ok(())
	}

	#[test]
	fn test_sub() -> Result<(), String> {
		for _ in 0..100 {
			for modulus in MODULI {
				let ctx = Rc::new(Context::new(&[*modulus], 8)?);
				let m = Modulus::new(*modulus).unwrap();

				let p = Poly::random(&ctx, Representation::PowerBasis);
				let q = Poly::random(&ctx, Representation::PowerBasis);
				let r = &p - &q;
				assert_eq!(r.representation, Representation::PowerBasis);
				let mut a = Vec::<u64>::from(&p);
				m.sub_vec(&mut a, &Vec::<u64>::from(&q));
				assert_eq!(Vec::<u64>::from(&r), a);

				let p = Poly::random(&ctx, Representation::Ntt);
				let q = Poly::random(&ctx, Representation::Ntt);
				let r = &p - &q;
				assert_eq!(r.representation, Representation::Ntt);
				let mut a = Vec::<u64>::from(&p);
				m.sub_vec(&mut a, &Vec::<u64>::from(&q));
				assert_eq!(Vec::<u64>::from(&r), a);
			}

			let ctx = Rc::new(Context::new(MODULI, 8)?);
			let p = Poly::random(&ctx, Representation::PowerBasis);
			let q = Poly::random(&ctx, Representation::PowerBasis);
			let mut a = Vec::<u64>::from(&p);
			let b = Vec::<u64>::from(&q);
			for i in 0..MODULI.len() {
				let m = Modulus::new(MODULI[i]).unwrap();
				m.sub_vec(&mut a[i * 8..(i + 1) * 8], &b[i * 8..(i + 1) * 8])
			}
			let r = &p - &q;
			assert_eq!(r.representation, Representation::PowerBasis);
			assert_eq!(Vec::<u64>::from(&r), a);
		}
		Ok(())
	}

	#[test]
	fn test_mul() -> Result<(), String> {
		for _ in 0..100 {
			for modulus in MODULI {
				let ctx = Rc::new(Context::new(&[*modulus], 8)?);
				let m = Modulus::new(*modulus).unwrap();

				let p = Poly::random(&ctx, Representation::Ntt);
				let q = Poly::random(&ctx, Representation::Ntt);
				let r = &p * &q;
				assert_eq!(r.representation, Representation::Ntt);
				let mut a = Vec::<u64>::from(&p);
				m.mul_vec(&mut a, &Vec::<u64>::from(&q));
				assert_eq!(Vec::<u64>::from(&r), a);
			}

			let ctx = Rc::new(Context::new(MODULI, 8)?);
			let p = Poly::random(&ctx, Representation::Ntt);
			let q = Poly::random(&ctx, Representation::Ntt);
			let mut a = Vec::<u64>::from(&p);
			let b = Vec::<u64>::from(&q);
			for i in 0..MODULI.len() {
				let m = Modulus::new(MODULI[i]).unwrap();
				m.mul_vec(&mut a[i * 8..(i + 1) * 8], &b[i * 8..(i + 1) * 8])
			}
			let r = &p * &q;
			assert_eq!(r.representation, Representation::Ntt);
			assert_eq!(Vec::<u64>::from(&r), a);
		}
		Ok(())
	}

	#[test]
	fn test_mul_shoup() -> Result<(), String> {
		for _ in 0..100 {
			for modulus in MODULI {
				let ctx = Rc::new(Context::new(&[*modulus], 8)?);
				let m = Modulus::new(*modulus).unwrap();

				let p = Poly::random(&ctx, Representation::Ntt);
				let q = Poly::random(&ctx, Representation::NttShoup);
				let r = &p * &q;
				assert_eq!(r.representation, Representation::Ntt);
				let mut a = Vec::<u64>::from(&p);
				m.mul_vec(&mut a, &Vec::<u64>::from(&q));
				assert_eq!(Vec::<u64>::from(&r), a);
			}

			let ctx = Rc::new(Context::new(MODULI, 8)?);
			let p = Poly::random(&ctx, Representation::Ntt);
			let q = Poly::random(&ctx, Representation::NttShoup);
			let mut a = Vec::<u64>::from(&p);
			let b = Vec::<u64>::from(&q);
			for i in 0..MODULI.len() {
				let m = Modulus::new(MODULI[i]).unwrap();
				m.mul_vec(&mut a[i * 8..(i + 1) * 8], &b[i * 8..(i + 1) * 8])
			}
			let r = &p * &q;
			assert_eq!(r.representation, Representation::Ntt);
			assert_eq!(Vec::<u64>::from(&r), a);
		}
		Ok(())
	}

	#[test]
	fn test_neg() -> Result<(), String> {
		for _ in 0..100 {
			for modulus in MODULI {
				let ctx = Rc::new(Context::new(&[*modulus], 8)?);
				let m = Modulus::new(*modulus).unwrap();

				let p = Poly::random(&ctx, Representation::PowerBasis);
				let r = -&p;
				assert_eq!(r.representation, Representation::PowerBasis);
				let mut a = Vec::<u64>::from(&p);
				m.neg_vec(&mut a);
				assert_eq!(Vec::<u64>::from(&r), a);

				let p = Poly::random(&ctx, Representation::Ntt);
				let r = -&p;
				assert_eq!(r.representation, Representation::Ntt);
				let mut a = Vec::<u64>::from(&p);
				m.neg_vec(&mut a);
				assert_eq!(Vec::<u64>::from(&r), a);
			}

			let ctx = Rc::new(Context::new(MODULI, 8)?);
			let p = Poly::random(&ctx, Representation::PowerBasis);
			let mut a = Vec::<u64>::from(&p);
			for i in 0..MODULI.len() {
				let m = Modulus::new(MODULI[i]).unwrap();
				m.neg_vec(&mut a[i * 8..(i + 1) * 8])
			}
			let r = -&p;
			assert_eq!(r.representation, Representation::PowerBasis);
			assert_eq!(Vec::<u64>::from(&r), a);
		}
		Ok(())
	}
}
