//! Create parameters for the BFV encryption scheme

use derive_builder::Builder;
use math::{
	rns::RnsContext,
	rq::{
		extender::Extender, scaler::Scaler, traits::TryConvertFrom, Context, Poly, Representation,
	},
	zq::{nfl::generate_opt_prime, ntt::NttOperator, Modulus},
};
use ndarray::ArrayView1;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use std::rc::Rc;

/// Parameters for the BFV encryption scheme.
#[derive(Debug, Builder, PartialEq)]
#[builder(build_fn(private, name = "fallible_build"))]
pub struct BfvParameters {
	/// Number of coefficients in a polynomial.
	pub(crate) polynomial_degree: usize,

	/// Modulus of the plaintext.
	plaintext_modulus: u64,

	/// Vector of coprime moduli q_i for the ciphertext.
	/// One and only one of `ciphertext_moduli` or `ciphertext_moduli_sizes` must be specified.
	pub(crate) ciphertext_moduli: Vec<u64>,

	/// Vector of the sized of the coprime moduli q_i for the ciphertext.
	/// One and only one of `ciphertext_moduli` or `ciphertext_moduli_sizes` must be specified.
	ciphertext_moduli_sizes: Vec<usize>,

	/// Error variance
	pub(crate) variance: usize,

	/// Context for the underlying polynomials
	#[builder(setter(skip))]
	pub(crate) ctx: Rc<Context>,

	/// Ntt operator for the SIMD plaintext, if possible.
	#[builder(setter(skip))]
	pub(crate) op: Option<Rc<NttOperator>>,

	/// Scaling polynomial for the plaintext
	#[builder(setter(skip))]
	pub(crate) delta: Poly,

	/// Q modulo the plaintext modulus
	#[builder(setter(skip))]
	pub(crate) q_mod_t: u64,

	/// Down scaler for the plaintext
	#[builder(setter(skip))]
	pub(crate) scaler: Scaler,

	/// Plaintext Modulus
	// #[builder(setter(skip))] // TODO: How can we handle this?
	pub(crate) plaintext: Modulus,

	/// Polynomial extender used in the homomorphic multiplication
	#[builder(setter(skip))]
	pub(crate) extender: Extender,

	/// Scaler used in the homomorphic multiplication
	#[builder(setter(skip))]
	pub(crate) rounder: Scaler,
}

impl BfvParameters {
	/// Returns the underlying polynomial degree
	pub fn degree(&self) -> usize {
		self.polynomial_degree
	}

	/// Returns a reference to the ciphertext moduli
	pub fn moduli(&self) -> &[u64] {
		&self.ciphertext_moduli
	}

	#[cfg(test)]
	pub fn default_one_modulus() -> Self {
		BfvParametersBuilder::default()
			.polynomial_degree(8)
			.plaintext_modulus(1153)
			.ciphertext_moduli(vec![4611686018326724609])
			.build()
			.unwrap()
	}

	#[cfg(test)]
	pub fn default_two_moduli() -> Self {
		BfvParametersBuilder::default()
			.polynomial_degree(8)
			.plaintext_modulus(1153)
			.ciphertext_moduli(vec![4611686018326724609, 4611686018309947393])
			.build()
			.unwrap()
	}
}

impl BfvParametersBuilder {
	/// Build a new `BfvParameters`.
	pub fn build(&self) -> Result<BfvParameters, BfvParametersBuilderError> {
		// Check the polynomial degree
		if self.polynomial_degree.is_none() {
			return Err(BfvParametersBuilderError::UninitializedField(
				"polynomial_degree",
			));
		}
		let polynomial_degree = self.polynomial_degree.unwrap();
		if polynomial_degree < 8 || !polynomial_degree.is_power_of_two() {
			return Err(BfvParametersBuilderError::ValidationError(
				"`polynomial_degree` must be a power of two larger or equal to 8".to_string(),
			));
		}

		// Check the plaintext modulus
		if self.plaintext_modulus.is_none() {
			return Err(BfvParametersBuilderError::UninitializedField(
				"plaintext_modulus",
			));
		}
		let plaintext_modulus = Modulus::new(self.plaintext_modulus.unwrap());
		if plaintext_modulus.is_none() {
			// TODO: More checks needed
			return Err(BfvParametersBuilderError::ValidationError(
				"`plaintext_modulus` must be larger or equal to 2".to_string(),
			));
		}
		let plaintext_modulus = plaintext_modulus.unwrap();

		// Check the ciphertext moduli
		if self.ciphertext_moduli.is_none() && self.ciphertext_moduli.is_none() {
			return Err(
				BfvParametersBuilderError::ValidationError("One and only one of `ciphertext_moduli` or `ciphertext_moduli_sizes` must be specified"
					.to_string())
			);
		}
		// TODO: This needs to be fixed when only sizes are specified
		assert!(self.ciphertext_moduli.is_some());
		let ciphertext_moduli = self.ciphertext_moduli.clone().unwrap();

		let variance = if let Some(var) = self.variance {
			var
		} else {
			1
		};
		if !(1..=16).contains(&variance) {
			return Err(BfvParametersBuilderError::ValidationError(
				"The variance should be an integer between 1 and 16".to_string(),
			));
		}

		let op = NttOperator::new(&plaintext_modulus, polynomial_degree);

		// Compute the scaling factors for the plaintext
		let rns = RnsContext::new(&ciphertext_moduli)?;

		let ctx = Rc::new(Context::new(&ciphertext_moduli, polynomial_degree)?);
		let plaintext_ctx = Rc::new(Context::new(&ciphertext_moduli[..1], polynomial_degree)?);
		let scaler = Scaler::new(
			&ctx,
			&plaintext_ctx,
			&BigUint::from(plaintext_modulus.modulus()),
			rns.modulus(),
		)?;

		// Compute the NttShoup representation of delta = -1/t mod Q
		let mut delta_rests = vec![];
		for m in &ciphertext_moduli {
			let q = Modulus::new(*m).unwrap();
			delta_rests.push(q.inv(q.neg(plaintext_modulus.modulus())).unwrap())
		}
		let delta = rns.lift(&ArrayView1::from(&delta_rests)); // -1/t mod Q
		let mut delta_poly = Poly::try_convert_from(&[delta], &ctx, Representation::PowerBasis)?;
		delta_poly.change_representation(Representation::NttShoup);

		// Compute Q mod t
		let q_mod_t = (rns.modulus() % plaintext_modulus.modulus())
			.to_u64()
			.unwrap();

		let mut extended_moduli = Vec::with_capacity(2 * ciphertext_moduli.len() + 1);
		for m in &ciphertext_moduli {
			extended_moduli.push(*m);
		}
		let mut upper_bound = u64::MAX >> 2;
		while extended_moduli.len() != 2 * ciphertext_moduli.len() + 1 {
			upper_bound =
				generate_opt_prime(62, 2 * polynomial_degree as u64, upper_bound).unwrap();
			if !ciphertext_moduli.contains(&upper_bound) {
				extended_moduli.push(upper_bound)
			}
		}
		let to_ctx = Rc::new(Context::new(&extended_moduli, polynomial_degree)?);
		let extender = Extender::new(&ctx, &to_ctx)?;
		let rounder = Scaler::new(
			&to_ctx,
			&ctx,
			&BigUint::from(plaintext_modulus.modulus()),
			rns.modulus(),
		)?;

		Ok(BfvParameters {
			polynomial_degree,
			plaintext_modulus: plaintext_modulus.modulus(),
			ciphertext_moduli,
			ciphertext_moduli_sizes: vec![],
			variance,
			ctx,
			op: op.map(Rc::new),
			delta: delta_poly,
			q_mod_t,
			scaler,
			plaintext: plaintext_modulus,
			extender,
			rounder,
		})
	}
}

#[cfg(test)]
mod tests {
	use super::{BfvParameters, BfvParametersBuilder};

	#[test]
	fn test_builder() {
		let params = BfvParametersBuilder::default().build();
		assert!(params.is_err_and(|e| e.to_string() == "`polynomial_degree` must be initialized"));

		let params = BfvParametersBuilder::default().polynomial_degree(7).build();
		assert!(params
			.is_err_and(|e| e.to_string()
				== "`polynomial_degree` must be a power of two larger or equal to 8"));

		let params = BfvParametersBuilder::default()
			.polynomial_degree(1023)
			.build();
		assert!(params
			.is_err_and(|e| e.to_string()
				== "`polynomial_degree` must be a power of two larger or equal to 8"));

		let params = BfvParametersBuilder::default()
			.polynomial_degree(1024)
			.build();
		assert!(params.is_err_and(|e| e.to_string() == "`plaintext_modulus` must be initialized"));

		let params = BfvParametersBuilder::default()
			.polynomial_degree(1024)
			.plaintext_modulus(0)
			.build();
		assert!(params
			.is_err_and(|e| e.to_string() == "`plaintext_modulus` must be larger or equal to 2"));

		let params = BfvParametersBuilder::default()
			.polynomial_degree(1024)
			.plaintext_modulus(2)
			.build();
		assert!(params
			.is_err_and(|e| e.to_string() == "One and only one of `ciphertext_moduli` or `ciphertext_moduli_sizes` must be specified"));

		let params = BfvParametersBuilder::default()
			.polynomial_degree(1024)
			.plaintext_modulus(2)
			.ciphertext_moduli(vec![])
			.build();
		assert!(params.is_err_and(|e| e.to_string() == "The list of moduli is empty"));

		let params = BfvParametersBuilder::default()
			.polynomial_degree(8)
			.plaintext_modulus(2)
			.ciphertext_moduli(vec![1])
			.build();
		assert!(params.is_err_and(|e| e.to_string() == "The modulus is invalid"));

		let params = BfvParametersBuilder::default()
			.polynomial_degree(8)
			.plaintext_modulus(2)
			.ciphertext_moduli(vec![2])
			.build();
		assert!(params.is_err_and(|e| e.to_string() == "Impossible to construct a Ntt operator"));

		let params = BfvParametersBuilder::default()
			.polynomial_degree(8)
			.plaintext_modulus(2)
			.ciphertext_moduli(vec![1153])
			.build();
		assert!(params.is_ok());

		let params = params.unwrap();
		assert_eq!(params.ciphertext_moduli, vec![1153]);
		assert_eq!(params.moduli(), vec![1153]);
		assert_eq!(params.plaintext_modulus, 2);
		assert_eq!(params.polynomial_degree, 8);
		assert_eq!(params.degree(), 8);
		assert_eq!(params.variance, 1);
		assert!(params.op.is_none());
	}

	#[test]
	fn test_default() {
		let params = BfvParameters::default_one_modulus();
		assert_eq!(params.ciphertext_moduli.len(), 1);
		assert!(params.op.is_some());

		let params = BfvParameters::default_two_moduli();
		assert_eq!(params.ciphertext_moduli.len(), 2);
		assert!(params.op.is_some());
	}
}
