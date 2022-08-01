//! Create parameters for the BFV encryption scheme

use derive_builder::Builder;
use math::{
	rq::Context,
	zq::{ntt::NttOperator, Modulus},
};
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
	op: Option<Rc<NttOperator>>,
}

impl BfvParameters {
	/// Returns the underlying polynomial degree
	pub fn degree(&self) -> usize {
		self.polynomial_degree
	}

	/// Returns the underlying plaintext modulus
	pub fn plaintext(&self) -> u64 {
		self.plaintext_modulus
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

	#[cfg(test)]
	pub fn default() -> Self {
		BfvParametersBuilder::default()
			.polynomial_degree(8)
			.plaintext_modulus(1153)
			.ciphertext_moduli(vec![4611686018326724609])
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

		Ok(BfvParameters {
			polynomial_degree,
			plaintext_modulus: plaintext_modulus.modulus(),
			ciphertext_moduli,
			ciphertext_moduli_sizes: vec![],
			variance,
			ctx: Rc::new(ctx.unwrap()),
			op: op.map(Rc::new),
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

		let params = BfvParameters::default();
		assert!(params.op.is_some());
	}
}
