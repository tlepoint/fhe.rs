#![warn(missing_docs, unused_imports)]

//! Polynomial scaler.

use super::{Context, Poly, Representation};
use crate::rns::RnsScaler;
use itertools::izip;
use ndarray::{Array2, Axis};
use num_bigint::BigUint;
use std::rc::Rc;

/// Context extender.
#[derive(Default, Debug, Clone, PartialEq)]
pub struct Scaler {
	from: Rc<Context>,
	to: Rc<Context>,
	scaler: RnsScaler,
}

impl Scaler {
	/// Create a scaler from a context `from` to a context `to`.
	pub fn new(
		from: &Rc<Context>,
		to: &Rc<Context>,
		numerator: &BigUint,
		denominator: &BigUint,
	) -> Result<Self, String> {
		let scaler = RnsScaler::new(&from.rns, &to.rns, numerator, denominator);

		Ok(Self {
			from: from.clone(),
			to: to.clone(),
			scaler,
		})
	}
}

impl Scaler {
	/// Scale a polynomial
	pub fn scale(&self, p: &Poly, floor: bool) -> Result<Poly, String> {
		if p.ctx.as_ref() != self.from.as_ref() {
			Err("The input polynomial does not have the correct context".to_string())
		} else if p.representation != Representation::PowerBasis {
			Err("The input polynomial should be in power basis representation".to_string())
		} else {
			let mut new_coefficients = Array2::<u64>::zeros((self.to.q.len(), self.to.degree));
			izip!(
				new_coefficients.axis_iter_mut(Axis(1)),
				p.coefficients.axis_iter(Axis(1))
			)
			.for_each(|(mut new_column, column)| {
				self.scaler.scale(&column, &mut new_column, floor)
			});
			Ok(Poly {
				ctx: self.to.clone(),
				representation: Representation::PowerBasis,
				allow_variable_time_computations: p.allow_variable_time_computations,
				coefficients: new_coefficients,
				coefficients_shoup: None,
			})
		}
	}
}

#[cfg(test)]
mod tests {
	use super::Scaler;
	use crate::rq::{Context, Poly, Representation};
	use itertools::Itertools;
	use num_bigint::BigUint;
	use std::rc::Rc;

	// Moduli to be used in tests.
	static Q: &[u64; 3] = &[
		4611686018282684417,
		4611686018326724609,
		4611686018309947393,
	];

	static P: &[u64; 3] = &[
		4611686018309947393,
		4611686018282684417,
		4611686018257518593,
	];

	#[test]
	fn test_scaler() -> Result<(), String> {
		let ntests = 100;
		let from = Rc::new(Context::new(Q, 8)?);
		let to = Rc::new(Context::new(P, 8)?);

		for numerator in &[1u64, 2, 3, 100, 1000, 4611686018326724610] {
			for denominator in &[1u64, 2, 3, 4, 100, 101, 1000, 1001, 4611686018326724610] {
				let n = BigUint::from(*numerator);
				let d = BigUint::from(*denominator);

				let scaler = Scaler::new(&from, &to, &n, &d)?;

				for _ in 0..ntests {
					let poly = Poly::random(&from, Representation::PowerBasis);
					let poly_biguint = Vec::<BigUint>::from(&poly);

					let scaled_poly = scaler.scale(&poly, true)?;
					let scaled_biguint = Vec::<BigUint>::from(&scaled_poly);

					let expected = poly_biguint
						.iter()
						.map(|i| {
							if i >= &(from.modulus() >> 1usize) {
								to.modulus()
									- (&(&(from.modulus() - i) * &n + &d - 1u64) / &d)
										% to.modulus()
							} else {
								((i * &n) / &d) % to.modulus()
							}
						})
						.collect_vec();
					assert_eq!(expected, scaled_biguint);
				}
			}
		}

		Ok(())
	}
}
