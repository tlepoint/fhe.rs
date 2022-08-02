//! Create parameters for the BFV encryption scheme

use derive_builder::Builder;
use math::{
	rns::RnsContext,
	rq::{scaler::Scaler, traits::TryConvertFrom, Context, Poly, Representation},
	zq::{ntt::NttOperator, Modulus},
};
use num_bigint::BigUint;
use std::rc::Rc;

/// Parameters for the BFV encryption scheme.
#[derive(Debug, Builder, PartialEq)]
#[builder(build_fn(private, name = "fallible_build"))]
pub struct BfvParameters {
	/// Number of coefficients in a polynomial.
	polynomial_degree: usize,

	/// Modulus of the plaintext.
	plaintext_modulus: u64,

	/// Vector of coprime moduli q_i for the ciphertext.
	/// One and only one of `ciphertext_moduli` or `ciphertext_moduli_sizes` must be specified.
	ciphertext_moduli: Vec<u64>,

	/// Vector of the sized of the coprime moduli q_i for the ciphertext.
	/// One and only one of `ciphertext_moduli` or `ciphertext_moduli_sizes` must be specified.
	ciphertext_moduli_sizes: Vec<usize>,

	/// Error variance
	variance: usize,

	/// Context for the underlying polynomials
	#[builder(setter(skip))]
	ctx: Rc<Context>,

	/// Ntt operator for the SIMD plaintext, if possible.
	#[builder(setter(skip))]
	op: Option<Rc<NttOperator>>,

	/// Scaling polynomial for the plaintext
	#[builder(setter(skip))]
	delta: Poly,

	/// Down scaler for the plaintext
	#[builder(setter(skip))]
	scaler: Scaler,

	/// Plaintext Modulus
	// #[builder(setter(skip))] // TODO: How can we handle this?
	plaintext: Modulus,
}

impl BfvParameters {
	/// Returns the underlying polynomial degree
	pub fn degree(&self) -> usize {
		self.polynomial_degree
	}

	/// Returns the underlying plaintext modulus
	pub fn plaintext(&self) -> &Modulus {
		&self.plaintext
	}

	/// Returns the error variance
	pub fn variance(&self) -> usize {
		self.variance
	}

	/// Returns the underlying polynomial context.
	pub fn ctx(&self) -> &Rc<Context> {
		&self.ctx
	}

	/// Returns a reference to the ciphertext moduli
	pub fn moduli(&self) -> &[u64] {
		&self.ciphertext_moduli
	}

	/// Returns a reference to the Ntt operator if it exists
	pub fn simd_operator(&self) -> &Option<Rc<NttOperator>> {
		&self.op
	}

	/// Returns the scaling constant polynomial in NttShoup representation
	pub fn delta(&self) -> &Poly {
		&self.delta
	}

	/// Returns a reference to the scaler that enables to multiply by plaintext / ciphertext_modulus
	pub fn scaler(&self) -> &Scaler {
		&self.scaler
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

		let ctx = Context::new(&ciphertext_moduli, polynomial_degree);
		if ctx.is_none() {
			return Err(BfvParametersBuilderError::ValidationError(
				"The polynomial context could not be generated".to_string(),
			));
		}

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

		let plaintext_modulus = plaintext_modulus.unwrap();
		let op = NttOperator::new(&plaintext_modulus, polynomial_degree);

		// Compute the scaling factors for the plaintext
		let rns = RnsContext::new(&ciphertext_moduli).unwrap();
		let delta = rns.modulus() / plaintext_modulus.modulus();

		let ctx = Rc::new(ctx.unwrap());
		let scaler = Scaler::new(
			&ctx,
			&BigUint::from(plaintext_modulus.modulus()),
			rns.modulus(),
		)?;

		let mut delta_vec = Vec::with_capacity(polynomial_degree);
		for _ in 0..polynomial_degree {
			delta_vec.push(delta.clone())
		}
		let delta_poly =
			Poly::try_convert_from(delta_vec.as_slice(), &ctx, Representation::NttShoup)?;

		Ok(BfvParameters {
			polynomial_degree,
			plaintext_modulus: plaintext_modulus.modulus(),
			ciphertext_moduli,
			ciphertext_moduli_sizes: vec![],
			variance,
			ctx,
			op: op.map(Rc::new),
			delta: delta_poly,
			scaler,
			plaintext: plaintext_modulus,
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
		assert!(
			params.is_err_and(|e| e.to_string() == "The polynomial context could not be generated")
		);

		let params = BfvParametersBuilder::default()
			.polynomial_degree(8)
			.plaintext_modulus(2)
			.ciphertext_moduli(vec![1])
			.build();
		assert!(
			params.is_err_and(|e| e.to_string() == "The polynomial context could not be generated")
		);

		let params = BfvParametersBuilder::default()
			.polynomial_degree(8)
			.plaintext_modulus(2)
			.ciphertext_moduli(vec![2])
			.build();
		assert!(
			params.is_err_and(|e| e.to_string() == "The polynomial context could not be generated")
		);

		let params = BfvParametersBuilder::default()
			.polynomial_degree(8)
			.plaintext_modulus(2)
			.ciphertext_moduli(vec![1153])
			.build();
		assert!(params.is_ok());

		let params = params.unwrap();
		assert_eq!(params.ciphertext_moduli, vec![1153]);
		assert_eq!(params.ciphertext_moduli_sizes, vec![]); // TODO: Should be fixed
		assert_eq!(params.plaintext_modulus, 2);
		assert_eq!(params.polynomial_degree, 8);
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
