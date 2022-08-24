//! Evaluation keys for the BFV encryption scheme.

use fhers_traits::{DeserializeParametrized, FheParametrized, Serialize};

use crate::{
	bfv::{
		leveled::{LeveledEvaluationKey, LeveledEvaluationKeyBuilder},
		BfvParameters, Ciphertext,
	},
	Error, Result,
};

use super::SecretKey;

/// Evaluation key for the BFV encryption scheme.
/// An evaluation key enables one or several of the following operations:
/// - column rotation
/// - row rotation
/// - relinearization
/// - oblivious expansion
/// - inner sum
#[derive(Debug, PartialEq, Eq)]
pub struct EvaluationKey {
	leveled_ek: LeveledEvaluationKey,
}

impl EvaluationKey {
	/// Reports whether the evaluation key enables to compute an homomorphic
	/// inner sums.
	pub fn supports_inner_sum(&self) -> bool {
		self.leveled_ek.supports_inner_sum()
	}

	/// Computes the homomorphic inner sum.
	pub fn computes_inner_sum(&self, ct: &Ciphertext) -> Result<Ciphertext> {
		self.leveled_ek.computes_inner_sum(ct)
	}

	/// Reports whether the evaluation key enables to rotate the rows of the
	/// plaintext.
	pub fn supports_row_rotation(&self) -> bool {
		self.leveled_ek.supports_row_rotation()
	}

	/// Homomorphically rotate the rows of the plaintext
	pub fn rotates_row(&self, ct: &Ciphertext) -> Result<Ciphertext> {
		self.leveled_ek.rotates_row(ct)
	}

	/// Reports whether the evaluation key enables to rotate the columns of the
	/// plaintext.
	pub fn supports_column_rotation_by(&self, i: usize) -> bool {
		self.leveled_ek.supports_column_rotation_by(i)
	}

	/// Homomorphically rotate the columns of the plaintext
	pub fn rotates_column_by(&self, ct: &Ciphertext, i: usize) -> Result<Ciphertext> {
		self.leveled_ek.rotates_column_by(ct, i)
	}

	/// Reports whether the evaluation key supports oblivious expansion.
	pub fn supports_expansion(&self, level: usize) -> bool {
		self.leveled_ek.supports_expansion(level)
	}

	/// Obliviously expands the ciphertext. Returns an error if this evaluation
	/// does not support expansion to level = ceil(log2(size)), or if the
	/// ciphertext does not have size 2. The output is a vector of `size`
	/// ciphertexts.
	pub fn expands(&self, ct: &Ciphertext, size: usize) -> Result<Vec<Ciphertext>> {
		self.leveled_ek.expands(ct, size)
	}
}

/// Builder for an evaluation key from the secret key.
pub struct EvaluationKeyBuilder {
	builder: LeveledEvaluationKeyBuilder,
}

impl EvaluationKeyBuilder {
	/// Creates a new builder from the [`SecretKey`], for operations on
	/// ciphertexts.
	pub fn new(sk: &SecretKey) -> Self {
		Self {
			builder: LeveledEvaluationKeyBuilder::new(sk, 0, 0).unwrap(),
		}
	}

	/// Allow expansion by this evaluation key.
	#[allow(unused_must_use)]
	pub fn enable_expansion(&mut self, level: usize) -> Result<&mut Self> {
		self.builder.enable_expansion(level)?;
		Ok(self)
	}

	/// Allow this evaluation key to compute homomorphic inner sums.
	#[allow(unused_must_use)]
	pub fn enable_inner_sum(&mut self) -> Result<&mut Self> {
		self.builder.enable_inner_sum()?;
		Ok(self)
	}

	/// Allow this evaluation key to homomorphically rotate the plaintext rows.
	#[allow(unused_must_use)]
	pub fn enable_row_rotation(&mut self) -> Result<&mut Self> {
		self.builder.enable_row_rotation()?;
		Ok(self)
	}

	/// Allow this evaluation key to homomorphically rotate the plaintext
	/// columns.
	#[allow(unused_must_use)]
	pub fn enable_column_rotation(&mut self, i: usize) -> Result<&mut Self> {
		self.builder.enable_column_rotation(i)?;
		Ok(self)
	}

	/// Build an [`EvaluationKey`] with the specified attributes.
	pub fn build(&mut self) -> Result<EvaluationKey> {
		Ok(EvaluationKey {
			leveled_ek: self.builder.build()?,
		})
	}
}

impl Serialize for EvaluationKey {
	fn to_bytes(&self) -> Vec<u8> {
		self.leveled_ek.to_bytes()
	}
}

impl FheParametrized for EvaluationKey {
	type Parameters = BfvParameters;
}

impl DeserializeParametrized for EvaluationKey {
	type Error = Error;

	fn from_bytes(bytes: &[u8], par: &std::sync::Arc<Self::Parameters>) -> Result<Self> {
		Ok(EvaluationKey {
			leveled_ek: LeveledEvaluationKey::from_bytes(bytes, par)?,
		})
	}
}

#[cfg(test)]
mod tests {
	use crate::bfv::{
		BfvParameters, Encoding, EvaluationKey, EvaluationKeyBuilder, Plaintext, SecretKey,
	};
	use fhers_traits::{
		DeserializeParametrized, FheDecoder, FheDecrypter, FheEncoder, FheEncrypter, Serialize,
	};
	use itertools::izip;
	use std::{error::Error, sync::Arc};

	#[test]
	fn test_builder() -> Result<(), Box<dyn Error>> {
		let params = Arc::new(BfvParameters::default(2));
		let sk = SecretKey::random(&params);
		let mut builder = EvaluationKeyBuilder::new(&sk);

		assert!(!builder.build()?.supports_row_rotation());
		assert!(!builder.build()?.supports_column_rotation_by(0));
		assert!(!builder.build()?.supports_column_rotation_by(1));
		assert!(!builder.build()?.supports_inner_sum());
		assert!(!builder.build()?.supports_expansion(1));
		assert!(builder.build()?.supports_expansion(0));
		assert!(builder.enable_column_rotation(0).is_err());
		assert!(builder
			.enable_expansion(64 - params.degree().leading_zeros() as usize)
			.is_err());

		builder.enable_column_rotation(1)?;
		assert!(builder.build()?.supports_column_rotation_by(1));
		assert!(!builder.build()?.supports_row_rotation());
		assert!(!builder.build()?.supports_inner_sum());
		assert!(!builder.build()?.supports_expansion(1));

		builder.enable_row_rotation()?;
		assert!(builder.build()?.supports_row_rotation());
		assert!(!builder.build()?.supports_inner_sum());
		assert!(!builder.build()?.supports_expansion(1));

		builder.enable_inner_sum()?;
		assert!(builder.build()?.supports_inner_sum());
		assert!(builder.build()?.supports_expansion(1));
		assert!(!builder
			.build()?
			.supports_expansion(64 - 1 - params.degree().leading_zeros() as usize));

		builder.enable_expansion(64 - 1 - params.degree().leading_zeros() as usize)?;
		assert!(builder
			.build()?
			.supports_expansion(64 - 1 - params.degree().leading_zeros() as usize));

		assert!(builder.build().is_ok());

		// Enabling inner sum enables row rotation and a few column rotations :)
		let ek = EvaluationKeyBuilder::new(&sk).enable_inner_sum()?.build()?;
		assert!(ek.supports_inner_sum());
		assert!(ek.supports_row_rotation());
		let mut i = 1;
		while i < params.degree() / 2 {
			assert!(ek.supports_column_rotation_by(i));
			i *= 2
		}
		assert!(!ek.supports_column_rotation_by(params.degree() / 2 - 1));

		Ok(())
	}

	#[test]
	fn test_inner_sum() -> Result<(), Box<dyn Error>> {
		for params in [
			Arc::new(BfvParameters::default(2)),
			Arc::new(BfvParameters::default(5)),
		] {
			for _ in 0..25 {
				let sk = SecretKey::random(&params);
				let ek = EvaluationKeyBuilder::new(&sk).enable_inner_sum()?.build()?;

				let v = params.plaintext.random_vec(params.degree());
				let expected = params
					.plaintext
					.reduce_u128(v.iter().map(|vi| *vi as u128).sum());

				let pt = Plaintext::try_encode(&v as &[u64], Encoding::simd(), &params)?;
				let ct = sk.try_encrypt(&pt)?;

				let ct2 = ek.computes_inner_sum(&ct)?;
				let pt = sk.try_decrypt(&ct2)?;
				assert_eq!(
					Vec::<u64>::try_decode(&pt, Encoding::simd())?,
					vec![expected; params.degree()]
				)
			}
		}
		Ok(())
	}

	#[test]
	fn test_row_rotation() -> Result<(), Box<dyn Error>> {
		for params in [
			Arc::new(BfvParameters::default(2)),
			Arc::new(BfvParameters::default(5)),
		] {
			for _ in 0..50 {
				let sk = SecretKey::random(&params);
				let ek = EvaluationKeyBuilder::new(&sk)
					.enable_row_rotation()?
					.build()?;

				let v = params.plaintext.random_vec(params.degree());
				let row_size = params.degree() >> 1;
				let mut expected = vec![0u64; params.degree()];
				expected[..row_size].copy_from_slice(&v[row_size..]);
				expected[row_size..].copy_from_slice(&v[..row_size]);

				let pt = Plaintext::try_encode(&v as &[u64], Encoding::simd(), &params)?;
				let ct = sk.try_encrypt(&pt)?;

				let ct2 = ek.rotates_row(&ct)?;
				let pt = sk.try_decrypt(&ct2)?;
				assert_eq!(Vec::<u64>::try_decode(&pt, Encoding::simd())?, expected)
			}
		}
		Ok(())
	}

	#[test]
	fn test_column_rotation() -> Result<(), Box<dyn Error>> {
		for params in [
			Arc::new(BfvParameters::default(2)),
			Arc::new(BfvParameters::default(5)),
		] {
			let row_size = params.degree() >> 1;
			for _ in 0..50 {
				for i in 1..row_size {
					let sk = SecretKey::random(&params);
					let ek = EvaluationKeyBuilder::new(&sk)
						.enable_column_rotation(i)?
						.build()?;

					let v = params.plaintext.random_vec(params.degree());
					let row_size = params.degree() >> 1;
					let mut expected = vec![0u64; params.degree()];
					expected[..row_size - i].copy_from_slice(&v[i..row_size]);
					expected[row_size - i..row_size].copy_from_slice(&v[..i]);
					expected[row_size..2 * row_size - i].copy_from_slice(&v[row_size + i..]);
					expected[2 * row_size - i..].copy_from_slice(&v[row_size..row_size + i]);

					let pt = Plaintext::try_encode(&v as &[u64], Encoding::simd(), &params)?;
					let ct = sk.try_encrypt(&pt)?;

					let ct2 = ek.rotates_column_by(&ct, i)?;
					let pt = sk.try_decrypt(&ct2)?;
					assert_eq!(Vec::<u64>::try_decode(&pt, Encoding::simd())?, expected)
				}
			}
		}
		Ok(())
	}

	#[test]
	fn test_expansion() -> Result<(), Box<dyn Error>> {
		for params in [
			Arc::new(BfvParameters::default(2)),
			Arc::new(BfvParameters::default(5)),
		] {
			let log_degree = 64 - 1 - params.degree().leading_zeros();
			for _ in 0..15 {
				for i in 1..1 + log_degree as usize {
					let sk = SecretKey::random(&params);
					let ek = EvaluationKeyBuilder::new(&sk)
						.enable_expansion(i)?
						.build()?;

					assert!(ek.supports_expansion(i));
					assert!(!ek.supports_expansion(i + 1));
					let v = params.plaintext.random_vec(1 << i);
					let pt = Plaintext::try_encode(&v as &[u64], Encoding::poly(), &params)?;
					let ct = sk.try_encrypt(&pt)?;

					let ct2 = ek.expands(&ct, 1 << i)?;
					assert_eq!(ct2.len(), 1 << i);
					for (vi, ct2i) in izip!(&v, &ct2) {
						let mut expected = vec![0u64; params.degree()];
						expected[0] = params.plaintext.mul(*vi, (1 << i) as u64);
						let pt = sk.try_decrypt(ct2i)?;
						assert_eq!(expected, Vec::<u64>::try_decode(&pt, Encoding::poly())?);
						println!("Noise: {:?}", unsafe { sk.measure_noise(ct2i) })
					}
				}
			}
		}
		Ok(())
	}

	#[test]
	fn test_serialize() -> Result<(), Box<dyn Error>> {
		for params in [
			Arc::new(BfvParameters::default(1)),
			Arc::new(BfvParameters::default(2)),
		] {
			let sk = SecretKey::random(&params);

			let ek = EvaluationKeyBuilder::new(&sk).build()?;
			let bytes = ek.to_bytes();
			assert_eq!(ek, EvaluationKey::from_bytes(&bytes, &params)?);

			if params.moduli.len() > 1 {
				let ek = EvaluationKeyBuilder::new(&sk)
					.enable_row_rotation()?
					.build()?;
				let bytes = ek.to_bytes();
				assert_eq!(ek, EvaluationKey::from_bytes(&bytes, &params)?);

				let ek = EvaluationKeyBuilder::new(&sk).enable_inner_sum()?.build()?;
				let bytes = ek.to_bytes();
				assert_eq!(ek, EvaluationKey::from_bytes(&bytes, &params)?);

				let ek = EvaluationKeyBuilder::new(&sk)
					.enable_expansion(params.degree().ilog2() as usize)?
					.build()?;
				let bytes = ek.to_bytes();
				assert_eq!(ek, EvaluationKey::from_bytes(&bytes, &params)?);

				let ek = EvaluationKeyBuilder::new(&sk)
					.enable_inner_sum()?
					.enable_expansion(params.degree().ilog2() as usize)?
					.build()?;
				let bytes = ek.to_bytes();
				assert_eq!(ek, EvaluationKey::from_bytes(&bytes, &params)?);
			}
		}
		Ok(())
	}
}
