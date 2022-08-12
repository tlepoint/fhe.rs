//! Evaluation keys for the BFV encryption scheme.

use crate::{traits::TryConvertFrom, Ciphertext};

use super::{GaloisKey, RelinearizationKey, SecretKey};
use fhers_protos::protos::bfv::{
	EvaluationKey as EvaluationKeyProto, GaloisKey as GaloisKeyProto,
	RelinearizationKey as RelinearizationKeyProto,
};
use math::rq::Poly;
use math::zq::Modulus;
use protobuf::MessageField;
use std::collections::{HashMap, HashSet};
use std::rc::Rc;
use zeroize::{Zeroize, ZeroizeOnDrop};

/// Evaluation key for the BFV encryption scheme.
/// An evaluation key enables one or several of the following operations:
/// - column rotation
/// - row rotation
/// - relinearization
/// - oblivious expansion
/// - inner sum
#[derive(Debug, PartialEq, Eq)]
pub struct EvaluationKey {
	inner_sum: bool,

	/// Relinearization key
	rk: Option<RelinearizationKey>,

	/// Galois keys
	gk: HashMap<usize, GaloisKey>,
}

impl EvaluationKey {
	/// Reports whether the evaluation key enables to compute an homomorphic inner sums.
	pub fn supports_inner_sum(&self) -> bool {
		self.inner_sum
	}

	/// Computes the homomorphic inner sum.
	pub fn computes_inner_sum(&self, ct: &Ciphertext) -> Result<Ciphertext, String> {
		if !self.supports_inner_sum() {
			Err("This key does not support the inner sum functionality".to_string())
		} else {
			let mut out = ct.clone();

			let mut i = 1;
			while i < ct.par.degree() / 2 {
				let gk = self.gk.get(&i).unwrap();
				out += &gk.relinearize(&out)?;
				i *= 2
			}

			let gk = self.gk.get(&0).unwrap();
			out += &gk.relinearize(&out)?;

			Ok(out)
		}
	}

	/// Reports whether the evaluation key enables to rotate the rows of the plaintext.
	pub fn supports_row_rotation(&self) -> bool {
		self.gk.contains_key(&0)
	}

	/// Homomorphically rotate the rows of the plaintext
	pub fn rotates_row(&self, ct: &Ciphertext) -> Result<Ciphertext, String> {
		if !self.supports_row_rotation() {
			Err("This key does not support the row rotation functionality".to_string())
		} else {
			let gk = self.gk.get(&0).unwrap();
			gk.relinearize(ct)
		}
	}

	/// Reports whether the evaluation key enables to rotate the columns of the plaintext.
	pub fn supports_column_rotation_by(&self, i: usize) -> bool {
		if i == 0 {
			false
		} else {
			self.gk.contains_key(&i)
		}
	}

	/// Homomorphically rotate the columns of the plaintext
	pub fn rotates_column_by(&self, ct: &Ciphertext, i: usize) -> Result<Ciphertext, String> {
		if !self.supports_column_rotation_by(i) {
			Err("This key does not support rotating the columns by this index".to_string())
		} else {
			let gk = self.gk.get(&i).unwrap();
			gk.relinearize(ct)
		}
	}

	/// Reports whether the evaluation key enable to perform relinearizations
	pub fn supports_relinearization(&self) -> bool {
		self.rk.is_some()
	}

	/// Relinearizes the expanded ciphertext
	pub fn relinearizes(&self, c0: &mut Poly, c1: &mut Poly, c2: &Poly) -> Result<(), String> {
		if !self.supports_relinearization() {
			Err("This key does not support relinearization".to_string())
		} else {
			self.rk.as_ref().unwrap().relinearize(c0, c1, c2)
		}
	}
}

/// Builder for an evaluation key from the secret key.
pub struct EvaluationKeyBuilder {
	relin: bool,
	inner_sum: bool,
	row_rotation: bool,
	column_rotation: HashSet<usize>,
	sk: SecretKey,
}

impl Zeroize for EvaluationKeyBuilder {
	fn zeroize(&mut self) {
		self.sk.zeroize()
	}
}

impl ZeroizeOnDrop for EvaluationKeyBuilder {}

impl EvaluationKeyBuilder {
	/// Creates a new builder form the [`SecretKey`]
	pub fn new(sk: &SecretKey) -> Self {
		Self {
			sk: sk.clone(),
			relin: false,
			inner_sum: false,
			row_rotation: false,
			column_rotation: HashSet::default(),
		}
	}

	/// Allow relinearizations by this evaluation key.
	#[allow(unused_must_use)]
	pub fn enable_relinearization(&mut self) -> &mut Self {
		self.relin = true;
		self
	}

	/// Allow this evaluation to compute homomorphic inner sums.
	#[allow(unused_must_use)]
	pub fn enable_inner_sum(&mut self) -> &mut Self {
		self.inner_sum = true;
		self
	}

	/// Allow this evaluation to homomorphically rotate the plaintext rows.
	#[allow(unused_must_use)]
	pub fn enable_row_rotation(&mut self) -> &mut Self {
		self.row_rotation = true;
		self
	}

	/// Allow this evaluation to homomorphically rotate the plaintext columns.
	#[allow(unused_must_use)]
	pub fn enable_column_rotation(&mut self, i: usize) -> Result<&mut Self, String> {
		if i == 0 || i >= self.sk.par.degree() {
			Err("Invalid column index".to_string())
		} else {
			self.column_rotation.insert(i);
			Ok(self)
		}
	}

	/// Build an [`EvaluationKey`] with the specified attributes.
	pub fn build(&self) -> Result<EvaluationKey, String> {
		let mut ek = EvaluationKey {
			rk: None,
			gk: HashMap::default(),
			inner_sum: false,
		};

		let mut indices = self.column_rotation.clone();

		if self.relin {
			ek.rk = Some(RelinearizationKey::new(&self.sk)?)
		}

		if self.row_rotation {
			indices.insert(0);
		}

		if self.inner_sum {
			ek.inner_sum = true;

			// Add the required indices to the set of indices
			indices.insert(0);
			let mut i = 1;
			while i < self.sk.par.degree() / 2 {
				indices.insert(i);
				i *= 2
			}
		}

		let q = Modulus::new(2 * self.sk.par.degree() as u64)?;
		for index in indices {
			if index == 0 {
				ek.gk
					.insert(0, GaloisKey::new(&self.sk, 2 * self.sk.par.degree() - 1)?);
			} else {
				let exp = q.pow(3, index as u64);
				ek.gk.insert(index, GaloisKey::new(&self.sk, exp as usize)?);
			}
		}

		Ok(ek)
	}
}

impl From<&EvaluationKey> for EvaluationKeyProto {
	fn from(ek: &EvaluationKey) -> Self {
		let mut proto = EvaluationKeyProto::new();
		for (_, gk) in ek.gk.iter() {
			proto.gk.push(GaloisKeyProto::from(gk))
		}
		if let Some(rk) = &ek.rk {
			proto.rk = MessageField::some(RelinearizationKeyProto::from(rk))
		}
		proto
	}
}

impl TryConvertFrom<&EvaluationKeyProto> for EvaluationKey {
	type Error = String;

	fn try_convert_from(
		value: &EvaluationKeyProto,
		par: &Rc<crate::BfvParameters>,
	) -> Result<Self, Self::Error> {
		let mut rk = None;
		if value.rk.is_some() {
			rk = Some(RelinearizationKey::try_convert_from(
				value.rk.as_ref().unwrap(),
				par,
			)?);
		}

		let mut gk = HashMap::new();
		for gkp in &value.gk {
			let key = GaloisKey::try_convert_from(gkp, par)?;
			if key.exponent == 2 * par.degree() - 1 {
				// row rotation key
				gk.insert(0, key);
			} else {
				let q = Modulus::new(2 * par.degree() as u64)?;
				let mut i = 1usize;
				while i < par.degree() / 2 {
					let e = q.pow(3, i as u64);
					if key.exponent == e as usize {
						gk.insert(i, key);
						break;
					}
					i *= 2
				}
			}
		}

		let mut inner_sum = gk.contains_key(&0);
		let mut i = 1usize;
		while i < par.degree() / 2 {
			inner_sum &= gk.contains_key(&i);
			i *= 2
		}

		Ok(EvaluationKey { inner_sum, rk, gk })
	}
}

#[cfg(test)]
mod tests {
	use super::{EvaluationKey, EvaluationKeyBuilder};
	use crate::{
		traits::{Decoder, Decryptor, Encoder, Encryptor, TryConvertFrom},
		BfvParameters, Encoding, Plaintext, SecretKey,
	};
	use fhers_protos::protos::bfv::EvaluationKey as EvaluationKeyProto;
	use std::rc::Rc;

	#[test]
	fn test_builder() -> Result<(), String> {
		let params = Rc::new(BfvParameters::default(2));
		let sk = SecretKey::random(&params);
		let mut builder = EvaluationKeyBuilder::new(&sk);

		assert!(!builder.build()?.supports_row_rotation());
		assert!(!builder.build()?.supports_column_rotation_by(0));
		assert!(!builder.build()?.supports_column_rotation_by(1));
		assert!(!builder.build()?.supports_inner_sum());
		assert!(builder.enable_column_rotation(0).is_err());

		builder.enable_column_rotation(1)?;
		assert!(builder.build()?.supports_column_rotation_by(1));
		assert!(!builder.build()?.supports_row_rotation());
		assert!(!builder.build()?.supports_inner_sum());

		builder.enable_row_rotation();
		assert!(builder.build()?.supports_row_rotation());
		assert!(!builder.build()?.supports_inner_sum());

		builder.enable_inner_sum();
		assert!(builder.build()?.supports_inner_sum());

		assert!(builder.build().is_ok());

		// Enabling inner sum enables row rotation and a few column rotations :)
		let ek = EvaluationKeyBuilder::new(&sk).enable_inner_sum().build()?;
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
	fn test_inner_sum() -> Result<(), String> {
		for params in [Rc::new(BfvParameters::default(2))] {
			for _ in 0..50 {
				let mut sk = SecretKey::random(&params);
				let ek = EvaluationKeyBuilder::new(&sk).enable_inner_sum().build()?;

				let v = params.plaintext.random_vec(params.degree());
				let expected = params
					.plaintext
					.reduce_u128(v.iter().map(|vi| *vi as u128).sum());

				let pt = Plaintext::try_encode(&v as &[u64], Encoding::Simd, &params)?;
				let mut ct = sk.encrypt(&pt)?;

				let ct2 = ek.computes_inner_sum(&mut ct)?;
				let pt = sk.decrypt(&ct2)?;
				assert_eq!(
					Vec::<u64>::try_decode(&pt, Encoding::Simd)?,
					vec![expected; params.degree()]
				)
			}
		}
		Ok(())
	}

	#[test]
	fn test_row_rotation() -> Result<(), String> {
		for params in [Rc::new(BfvParameters::default(2))] {
			for _ in 0..50 {
				let mut sk = SecretKey::random(&params);
				let ek = EvaluationKeyBuilder::new(&sk)
					.enable_row_rotation()
					.build()?;

				let v = params.plaintext.random_vec(params.degree());
				let row_size = params.degree() >> 1;
				let mut expected = vec![0u64; params.degree()];
				expected[..row_size].copy_from_slice(&v[row_size..]);
				expected[row_size..].copy_from_slice(&v[..row_size]);

				let pt = Plaintext::try_encode(&v as &[u64], Encoding::Simd, &params)?;
				let mut ct = sk.encrypt(&pt)?;

				let ct2 = ek.rotates_row(&mut ct)?;
				let pt = sk.decrypt(&ct2)?;
				assert_eq!(Vec::<u64>::try_decode(&pt, Encoding::Simd)?, expected)
			}
		}
		Ok(())
	}

	#[test]
	fn test_column_rotation() -> Result<(), String> {
		for params in [Rc::new(BfvParameters::default(2))] {
			let row_size = params.degree() >> 1;
			for _ in 0..50 {
				for i in 1..row_size {
					let mut sk = SecretKey::random(&params);
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

					let pt = Plaintext::try_encode(&v as &[u64], Encoding::Simd, &params)?;
					let mut ct = sk.encrypt(&pt)?;

					let ct2 = ek.rotates_column_by(&mut ct, i)?;
					let pt = sk.decrypt(&ct2)?;
					assert_eq!(Vec::<u64>::try_decode(&pt, Encoding::Simd)?, expected)
				}
			}
		}
		Ok(())
	}

	#[test]
	fn test_proto_conversion() -> Result<(), String> {
		for params in [
			Rc::new(BfvParameters::default(1)),
			Rc::new(BfvParameters::default(2)),
		] {
			let sk = SecretKey::random(&params);

			let ek = EvaluationKeyBuilder::new(&sk)
				.enable_row_rotation()
				.build()?;
			let proto = EvaluationKeyProto::from(&ek);
			assert_eq!(ek, EvaluationKey::try_convert_from(&proto, &params)?);

			let ek = EvaluationKeyBuilder::new(&sk).build()?;
			let proto = EvaluationKeyProto::from(&ek);
			assert_eq!(ek, EvaluationKey::try_convert_from(&proto, &params)?);

			let ek = EvaluationKeyBuilder::new(&sk)
				.enable_inner_sum()
				.enable_relinearization()
				.build()?;
			let proto = EvaluationKeyProto::from(&ek);
			assert_eq!(ek, EvaluationKey::try_convert_from(&proto, &params)?);
		}
		Ok(())
	}
}
