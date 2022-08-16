//! Implementation of operations over polynomials.

use super::{traits::TryConvertFrom, Poly, Representation};
use itertools::{izip, Itertools};
use ndarray::Array2;
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

impl AddAssign<Poly> for Poly {
	fn add_assign(&mut self, p: Poly) {
		*self += &p
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

impl Add for Poly {
	type Output = Poly;
	fn add(self, mut p: Poly) -> Poly {
		p += self;
		p
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

/// Compute the dot product between two iterators of polynomials.
/// Returna an error if the iterator counts are 0, or if any of the polynomial is not
/// in Ntt or NttShoup representation.
pub fn dot_product<'a, 'b, I, J>(p: I, q: J) -> Result<Poly, String>
where
	I: Iterator<Item = &'a Poly> + Clone,
	J: Iterator<Item = &'b Poly> + Clone,
{
	let p_count = p.clone().count();
	let q_count = q.clone().count();
	if p_count == 0 || q_count == 0 {
		return Err("At least one iterator is empty".to_string());
	}
	if p.clone()
		.any(|pi| pi.representation == Representation::PowerBasis)
	{
		return Err("All polynomials in p should be in Ntt representation".to_string());
	}
	if q.clone()
		.any(|pi| pi.representation == Representation::PowerBasis)
	{
		return Err("All polynomials in q should be in Ntt representation".to_string());
	}

	let p_first = p.clone().next().unwrap();

	// Initialize the accumulator
	let mut acc: Array2<u128> = Array2::zeros((p_first.ctx.q.len(), p_first.ctx.degree));
	// Current number of products accumulated
	let mut num_acc = vec![1u128; p_first.ctx.q.len()];
	// Maximum number of products that can be accumulated
	let mut max_acc = p_first
		.ctx
		.q
		.iter()
		.map(|qi| 1u128 << (2 * qi.modulus().leading_zeros()))
		.collect_vec();

	let min_of_max = max_acc.iter().min().unwrap();
	if p_count as u128 > *min_of_max || q_count as u128 > *min_of_max {
		for (pi, qi) in izip!(p, q) {
			izip!(
				acc.outer_iter_mut(),
				pi.coefficients().outer_iter(),
				qi.coefficients().outer_iter(),
				num_acc.iter_mut(),
				max_acc.iter_mut(),
				&p_first.ctx.q,
			)
			.for_each(|(mut accj, pij, qij, na, ma, m)| {
				izip!(accj.iter_mut(), pij.iter(), qij.iter())
					.for_each(|(accjk, pijk, qijk)| *accjk += (*pijk as u128) * (*qijk as u128));
				*na += 1;
				if *na == *ma {
					if p_first.allow_variable_time_computations {
						accj.iter_mut()
							.for_each(|accjk| *accjk = m.reduce_u128(*accjk) as u128);
					} else {
						accj.iter_mut()
							.for_each(|accjk| *accjk = unsafe { m.reduce_u128_vt(*accjk) as u128 });
					}
					*na = 1
				}
			});
		}
	} else {
		// We don't need to check the condition on the max, it should shave off a few cycles.
		for (pi, qi) in izip!(p, q) {
			izip!(
				acc.outer_iter_mut(),
				pi.coefficients().outer_iter(),
				qi.coefficients().outer_iter(),
				num_acc.iter_mut(),
			)
			.for_each(|(mut accj, pij, qij, na)| {
				izip!(accj.iter_mut(), pij.iter(), qij.iter())
					.for_each(|(accjk, pijk, qijk)| *accjk += (*pijk as u128) * (*qijk as u128));
				*na += 1;
			});
		}
	}

	// Last (conditional) reduction to create the coefficients
	let mut coeffs = Array2::zeros((p_first.ctx.q.len(), p_first.ctx.degree));
	izip!(
		coeffs.outer_iter_mut(),
		acc.outer_iter(),
		num_acc.iter(),
		&p_first.ctx.q,
	)
	.for_each(|(mut coeffsj, accj, na, m)| {
		if *na > 1 {
			if p_first.allow_variable_time_computations {
				izip!(coeffsj.iter_mut(), accj.iter())
					.for_each(|(cj, accjk)| *cj = m.reduce_u128(*accjk));
			} else {
				izip!(coeffsj.iter_mut(), accj.iter())
					.for_each(|(cj, accjk)| *cj = unsafe { m.reduce_u128_vt(*accjk) });
			}
		} else {
			izip!(coeffsj.iter_mut(), accj.iter()).for_each(|(cj, accjk)| *cj = *accjk as u64);
		}
	});

	Ok(Poly {
		ctx: p_first.ctx.clone(),
		representation: Representation::Ntt,
		allow_variable_time_computations: p_first.allow_variable_time_computations,
		coefficients: coeffs,
		coefficients_shoup: None,
	})
}

#[cfg(test)]
mod tests {
	use itertools::{izip, Itertools};

	use super::dot_product;
	use crate::{
		rq::{Context, Poly, Representation},
		zq::Modulus,
	};
	use std::sync::Arc;

	static MODULI: &[u64; 3] = &[1153, 4611686018326724609, 4611686018309947393];

	#[test]
	fn test_add() -> Result<(), String> {
		for _ in 0..100 {
			for modulus in MODULI {
				let ctx = Arc::new(Context::new(&[*modulus], 8)?);
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

			let ctx = Arc::new(Context::new(MODULI, 8)?);
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
				let ctx = Arc::new(Context::new(&[*modulus], 8)?);
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

			let ctx = Arc::new(Context::new(MODULI, 8)?);
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
				let ctx = Arc::new(Context::new(&[*modulus], 8)?);
				let m = Modulus::new(*modulus).unwrap();

				let p = Poly::random(&ctx, Representation::Ntt);
				let q = Poly::random(&ctx, Representation::Ntt);
				let r = &p * &q;
				assert_eq!(r.representation, Representation::Ntt);
				let mut a = Vec::<u64>::from(&p);
				m.mul_vec(&mut a, &Vec::<u64>::from(&q));
				assert_eq!(Vec::<u64>::from(&r), a);
			}

			let ctx = Arc::new(Context::new(MODULI, 8)?);
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
				let ctx = Arc::new(Context::new(&[*modulus], 8)?);
				let m = Modulus::new(*modulus).unwrap();

				let p = Poly::random(&ctx, Representation::Ntt);
				let q = Poly::random(&ctx, Representation::NttShoup);
				let r = &p * &q;
				assert_eq!(r.representation, Representation::Ntt);
				let mut a = Vec::<u64>::from(&p);
				m.mul_vec(&mut a, &Vec::<u64>::from(&q));
				assert_eq!(Vec::<u64>::from(&r), a);
			}

			let ctx = Arc::new(Context::new(MODULI, 8)?);
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
				let ctx = Arc::new(Context::new(&[*modulus], 8)?);
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

			let ctx = Arc::new(Context::new(MODULI, 8)?);
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

	#[test]
	fn test_dot_product() -> Result<(), String> {
		for _ in 0..20 {
			for modulus in MODULI {
				let ctx = Arc::new(Context::new(&[*modulus], 8)?);

				for len in 1..50 {
					let p = (0..len)
						.map(|_| Poly::random(&ctx, Representation::Ntt))
						.collect_vec();
					let q = (0..len)
						.map(|_| Poly::random(&ctx, Representation::Ntt))
						.collect_vec();
					let r = dot_product(p.iter(), q.iter())?;

					let mut expected = Poly::zero(&ctx, Representation::Ntt);
					izip!(&p, &q).for_each(|(pi, qi)| expected += pi * qi);
					assert_eq!(r, expected);
				}
			}

			let ctx = Arc::new(Context::new(MODULI, 8)?);
			for len in 1..50 {
				let p = (0..len)
					.map(|_| Poly::random(&ctx, Representation::Ntt))
					.collect_vec();
				let q = (0..len)
					.map(|_| Poly::random(&ctx, Representation::Ntt))
					.collect_vec();
				let r = dot_product(p.iter(), q.iter())?;

				let mut expected = Poly::zero(&ctx, Representation::Ntt);
				izip!(&p, &q).for_each(|(pi, qi)| expected += pi * qi);
				assert_eq!(r, expected);
			}
		}
		Ok(())
	}
}
