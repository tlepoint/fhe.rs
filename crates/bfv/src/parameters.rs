//! Create parameters for the BFV encryption scheme

use derive_builder::Builder;
use math::rns::RnsContext;

/// Parameters for the BFV encryption scheme.
#[derive(Debug, Builder)]
#[builder(build_fn(validate = "Self::validate"))]
pub struct BfvParameters {
	/// Number of coefficients in a polynomial.
	pub polynomial_degree: usize,
	/// Modulus of the plaintext.
	pub plaintext_modulus: u64,
	/// Vector of coprime moduli q_i for the ciphertext.
	/// One and only one of `ciphertext_moduli` or `ciphertext_moduli_sizes` must be specified.
	#[builder(setter(strip_option), default)]
	pub ciphertext_moduli: Option<Vec<u64>>,
	/// Vector of the sized of the coprime moduli q_i for the ciphertext.
	/// One and only one of `ciphertext_moduli` or `ciphertext_moduli_sizes` must be specified.
	#[builder(setter(strip_option), default)]
	pub ciphertext_moduli_sizes: Option<Vec<usize>>,
}

impl BfvParametersBuilder {
	fn validate(&self) -> Result<(), String> {
		println!("{:?}", self.plaintext_modulus);
		if let Some(degree) = self.polynomial_degree {
			if degree < 8 || !degree.is_power_of_two() {
				return Err(
					"`polynomial_degree` must be a power of two larger or equal to 8".to_string(),
				);
			}
		}

		if let Some(plaintext) = self.plaintext_modulus {
			if plaintext < 2 {
				return Err("`plaintext_modulus` must be larger or equal to 2".to_string());
			}
		}

		if self.ciphertext_moduli.is_none() && self.ciphertext_moduli.is_none() {
			return Err(
				"One and only one of `ciphertext_moduli` or `ciphertext_moduli_sizes` must be specified"
					.to_string(),
			);
		}

		if self.ciphertext_moduli.is_some() {
			let rns = RnsContext::new(&self.ciphertext_moduli.clone().unwrap().unwrap());
			if rns.is_none() {
				return Err(
					"The list of moduli is either empty, or the moduli are not coprime integers between 2 and 2^62-1".to_string(),
				);
			}
		}

		Ok(())
	}
}

#[cfg(test)]
mod tests {
	use super::BfvParametersBuilder;

	#[test]
	fn test_builder() {
		let params = BfvParametersBuilder::default().build();
		println!("{:?}", params);
		assert!(params.is_err_and(|e| {
			e.to_string()
			== "One and only one of `ciphertext_moduli` or `ciphertext_moduli_sizes` must be specified"
		}));

		let params = BfvParametersBuilder::default()
			.ciphertext_moduli(vec![0u64])
			.build();
		println!("{:?}", params);
		assert!(params.is_err_and(|e| e.to_string() == "The list of moduli is either empty, or the moduli are not coprime integers between 2 and 2^62-1"));

		let params = BfvParametersBuilder::default()
			.ciphertext_moduli(vec![0u64])
			.polynomial_degree(7)
			.build();
		println!("{:?}", params);
		assert!(params
			.is_err_and(|e| e.to_string()
				== "`polynomial_degree` must be a power of two larger or equal to 8"));

		let params = BfvParametersBuilder::default()
			.ciphertext_moduli(vec![0u64])
			.polynomial_degree(1023)
			.build();
		println!("{:?}", params);
		assert!(params
			.is_err_and(|e| e.to_string()
				== "`polynomial_degree` must be a power of two larger or equal to 8"));

		let params = BfvParametersBuilder::default()
			.ciphertext_moduli(vec![0u64])
			.polynomial_degree(1024usize)
			.build();
		println!("{:?}", params);
		assert!(params.is_err_and(|e| e.to_string() == "The list of moduli is either empty, or the moduli are not coprime integers between 2 and 2^62-1"));

		let params = BfvParametersBuilder::default()
			.ciphertext_moduli(vec![2u64])
			.polynomial_degree(1024usize)
			.build();
		println!("{:?}", params);
		assert!(params.is_err_and(|e| e.to_string() == "`plaintext_modulus` must be initialized"));

		let params = BfvParametersBuilder::default()
			.ciphertext_moduli(vec![2u64])
			.polynomial_degree(1024usize)
			.plaintext_modulus(0)
			.build();
		assert!(params
			.is_err_and(|e| e.to_string() == "`plaintext_modulus` must be larger or equal to 2"));

		let params = BfvParametersBuilder::default()
			.ciphertext_moduli(vec![2u64])
			.polynomial_degree(1024usize)
			.plaintext_modulus(11)
			.build();
		assert!(params.is_ok())
	}
}
