//! Evaluation keys for the BFV encryption scheme.

use super::{GaloisKey, RelinearizationKey, SecretKey};
use crate::bfv::{
	proto::bfv::{
		EvaluationKey as EvaluationKeyProto, GaloisKey as GaloisKeyProto,
		RelinearizationKey as RelinearizationKeyProto,
	},
	traits::TryConvertFrom,
	BfvParameters, Ciphertext,
};
use crate::{Error, Result};
use fhers_traits::{DeserializeWithContext, Serialize};
use math::rq::{traits::TryConvertFrom as TryConvertFromPoly, Poly, Representation};
use math::zq::Modulus;
use protobuf::{Message, MessageField};
use std::collections::{HashMap, HashSet};
use std::sync::Arc;
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
	par: Arc<BfvParameters>,

	/// Relinearization key
	rk: Option<RelinearizationKey>,

	/// Map from Galois keys exponents to Galois keys
	gk: HashMap<usize, GaloisKey>,

	/// Map from rotation index to Galois key exponent
	rot_to_gk_exponent: HashMap<usize, usize>,

	/// Monomials used in expansion
	monomials: Vec<Poly>,
}

impl EvaluationKey {
	/// Reports whether the evaluation key enables to compute an homomorphic
	/// inner sums.
	pub fn supports_inner_sum(&self) -> bool {
		if self.par.ciphertext_moduli.len() == 1 {
			false
		} else {
			let mut ret = self.gk.contains_key(&(self.par.degree() * 2 - 1));
			let mut i = 1;
			while i < self.par.degree() / 2 {
				ret &= self
					.gk
					.contains_key(self.rot_to_gk_exponent.get(&i).unwrap());
				i *= 2
			}
			ret
		}
	}

	/// Computes the homomorphic inner sum.
	pub fn computes_inner_sum(&self, ct: &Ciphertext) -> Result<Ciphertext> {
		if !self.supports_inner_sum() {
			Err(Error::DefaultError(
				"This key does not support the inner sum functionality".to_string(),
			))
		} else {
			let mut out = ct.clone();

			let mut i = 1;
			while i < ct.par.degree() / 2 {
				let gk = self
					.gk
					.get(self.rot_to_gk_exponent.get(&i).unwrap())
					.unwrap();
				out += &gk.relinearize(&out)?;
				i *= 2
			}

			let gk = self.gk.get(&(self.par.degree() * 2 - 1)).unwrap();
			out += &gk.relinearize(&out)?;

			Ok(out)
		}
	}

	/// Reports whether the evaluation key enables to rotate the rows of the
	/// plaintext.
	pub fn supports_row_rotation(&self) -> bool {
		if self.par.ciphertext_moduli.len() == 1 {
			false
		} else {
			self.gk.contains_key(&(self.par.degree() * 2 - 1))
		}
	}

	/// Homomorphically rotate the rows of the plaintext
	pub fn rotates_row(&self, ct: &Ciphertext) -> Result<Ciphertext> {
		if !self.supports_row_rotation() {
			Err(Error::DefaultError(
				"This key does not support the row rotation functionality".to_string(),
			))
		} else {
			let gk = self.gk.get(&(self.par.degree() * 2 - 1)).unwrap();
			gk.relinearize(ct)
		}
	}

	/// Reports whether the evaluation key enables to rotate the columns of the
	/// plaintext.
	pub fn supports_column_rotation_by(&self, i: usize) -> bool {
		if self.par.ciphertext_moduli.len() == 1 {
			false
		} else if let Some(exp) = self.rot_to_gk_exponent.get(&i) {
			self.gk.contains_key(exp)
		} else {
			false
		}
	}

	/// Homomorphically rotate the columns of the plaintext
	pub fn rotates_column_by(&self, ct: &Ciphertext, i: usize) -> Result<Ciphertext> {
		if !self.supports_column_rotation_by(i) {
			Err(Error::DefaultError(
				"This key does not support rotating the columns by this index".to_string(),
			))
		} else {
			let gk = self
				.gk
				.get(self.rot_to_gk_exponent.get(&i).unwrap())
				.unwrap();
			gk.relinearize(ct)
		}
	}

	/// Reports whether the evaluation key enable to perform relinearizations
	pub fn supports_relinearization(&self) -> bool {
		self.rk.is_some()
	}

	/// Relinearizes the expanded ciphertext
	pub fn relinearizes_new(&self, ct: &Ciphertext) -> Result<Ciphertext> {
		if !self.supports_relinearization() {
			Err(Error::DefaultError(
				"This key does not support relinearization".to_string(),
			))
		} else {
			self.rk.as_ref().unwrap().relinearizes(ct)
		}
	}

	/// Relinearize a 3-part ciphertext in place.
	pub fn relinearizes(&self, ct: &mut Ciphertext) -> Result<()> {
		if !self.supports_relinearization() {
			Err(Error::DefaultError(
				"This key does not support relinearization".to_string(),
			))
		} else if ct.c.len() != 3 {
			Err(Error::DefaultError(
				"The ciphertext does not have 3 parts".to_string(),
			))
		} else {
			let mut c2 = ct.c[2].clone();
			c2.change_representation(Representation::PowerBasis);
			let mut c1 = ct.c[1].clone();
			self.rk
				.as_ref()
				.unwrap()
				.relinearizes_with_poly(&c2, &mut ct.c[0], &mut c1)?;
			ct.c[1] = c1;
			ct.c.truncate(2);
			Ok(())
		}
	}

	/// Relinearize using polynomials.
	pub(crate) fn relinearizes_with_poly(
		&self,
		c2: &Poly,
		c0: &mut Poly,
		c1: &mut Poly,
	) -> Result<()> {
		if !self.supports_relinearization() {
			Err(Error::DefaultError(
				"This key does not support relinearization".to_string(),
			))
		} else {
			self.rk.as_ref().unwrap().relinearizes_with_poly(c2, c0, c1)
		}
	}

	/// Reports whether the evaluation key supports oblivious expansion.
	pub fn supports_expansion(&self, level: usize) -> bool {
		if self.par.ciphertext_moduli.len() == 1 {
			false
		} else {
			let mut ret = level < self.par.degree().leading_zeros() as usize;
			for l in 0..level {
				ret &= self.gk.contains_key(&((self.par.degree() >> l) + 1));
			}
			ret
		}
	}

	/// Obliviously expands the ciphertext. Returns an error if this evaluation
	/// does not support expansion to this level, or if the ciphertext does not
	/// have size 2. The output is a vector of 2^level ciphertexts.
	pub fn expands(&self, ct: &Ciphertext, level: usize) -> Result<Vec<Ciphertext>> {
		if ct.c.len() != 2 {
			Err(Error::DefaultError(
				"The ciphertext is not of size 2".to_string(),
			))
		} else if level == 0 {
			Ok(vec![ct.clone()])
		} else if self.supports_expansion(level) {
			let mut out = vec![Ciphertext::zero(&ct.par); 1 << level];
			out[0] = ct.clone();

			// We use the Oblivious expansion algorithm of
			// https://eprint.iacr.org/2019/1483.pdf
			for l in 0..level {
				let monomial = &self.monomials[l];
				let gk = self.gk.get(&((self.par.degree() >> l) + 1)).unwrap();
				for i in 0..(1 << l) {
					let sub = gk.relinearize(&out[i])?;
					out[(1 << l) | i] = &out[i] - &sub;
					out[(1 << l) | i].c[0] *= monomial;
					out[(1 << l) | i].c[1] *= monomial;
					out[i] += &sub;
				}
			}

			Ok(out)
		} else {
			Err(Error::DefaultError(
				"This key does not support expansion at this level".to_string(),
			))
		}
	}

	fn construct_rot_to_gk_exponent(par: &Arc<BfvParameters>) -> HashMap<usize, usize> {
		let mut m = HashMap::new();
		let q = Modulus::new(2 * par.degree() as u64).unwrap();
		for i in 1..par.degree() / 2 {
			let exp = q.pow(3, i as u64) as usize;
			m.insert(i, exp);
		}
		m
	}
}

impl Serialize for EvaluationKey {
	fn serialize(&self) -> Vec<u8> {
		let ekp = EvaluationKeyProto::from(self);
		ekp.write_to_bytes().unwrap()
	}
}

impl DeserializeWithContext<&Arc<BfvParameters>> for EvaluationKey {
	fn try_deserialize(bytes: &[u8], par: &Arc<BfvParameters>) -> Result<Self> {
		let gkp = EvaluationKeyProto::parse_from_bytes(bytes);
		if let Ok(gkp) = gkp {
			EvaluationKey::try_convert_from(&gkp, par)
		} else {
			Err(Error::DefaultError("Invalid serialization".to_string()))
		}
	}
	type Error = Error;
}

/// Builder for an evaluation key from the secret key.
pub struct EvaluationKeyBuilder {
	relin: bool,
	inner_sum: bool,
	row_rotation: bool,
	expansion_level: usize,
	column_rotation: HashSet<usize>,
	sk: SecretKey,
	rot_to_gk_exponent: HashMap<usize, usize>,
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
			expansion_level: 0,
			column_rotation: HashSet::default(),
			rot_to_gk_exponent: EvaluationKey::construct_rot_to_gk_exponent(&sk.par),
		}
	}

	/// Allow relinearizations by this evaluation key.
	#[allow(unused_must_use)]
	pub fn enable_relinearization(&mut self) -> Result<&mut Self> {
		if self.sk.par.ciphertext_moduli.len() == 1 {
			Err(Error::DefaultError(
				"Not enough moduli to enable relinearization".to_string(),
			))
		} else {
			self.relin = true;
			Ok(self)
		}
	}

	/// Allow relinearizations by this evaluation key.
	#[allow(unused_must_use)]
	pub fn enable_expansion(&mut self, level: usize) -> Result<&mut Self> {
		if self.sk.par.ciphertext_moduli.len() == 1 {
			Err(Error::DefaultError(
				"Not enough moduli to enable expansion".to_string(),
			))
		} else if level >= 64 - self.sk.par.degree().leading_zeros() as usize {
			Err(Error::DefaultError("Invalid level 2".to_string()))
		} else {
			self.expansion_level = level;
			Ok(self)
		}
	}

	/// Allow this evaluation key to compute homomorphic inner sums.
	#[allow(unused_must_use)]
	pub fn enable_inner_sum(&mut self) -> Result<&mut Self> {
		if self.sk.par.ciphertext_moduli.len() == 1 {
			Err(Error::DefaultError(
				"Not enough moduli to enable relinearization".to_string(),
			))
		} else {
			self.inner_sum = true;
			Ok(self)
		}
	}

	/// Allow this evaluation key to homomorphically rotate the plaintext rows.
	#[allow(unused_must_use)]
	pub fn enable_row_rotation(&mut self) -> Result<&mut Self> {
		if self.sk.par.ciphertext_moduli.len() == 1 {
			Err(Error::DefaultError(
				"Not enough moduli to enable relinearization".to_string(),
			))
		} else {
			self.row_rotation = true;
			Ok(self)
		}
	}

	/// Allow this evaluation key to homomorphically rotate the plaintext
	/// columns.
	#[allow(unused_must_use)]
	pub fn enable_column_rotation(&mut self, i: usize) -> Result<&mut Self> {
		if let Some(exp) = self.rot_to_gk_exponent.get(&i) {
			self.column_rotation.insert(*exp);
			Ok(self)
		} else {
			Err(Error::DefaultError("Invalid column index".to_string()))
		}
	}

	/// Build an [`EvaluationKey`] with the specified attributes.
	pub fn build(&self) -> Result<EvaluationKey> {
		let mut ek = EvaluationKey {
			rk: None,
			gk: HashMap::default(),
			par: self.sk.par.clone(),
			rot_to_gk_exponent: EvaluationKey::construct_rot_to_gk_exponent(&self.sk.par),
			monomials: vec![],
		};

		let mut indices = self.column_rotation.clone();

		if self.relin {
			ek.rk = Some(RelinearizationKey::new(&self.sk)?)
		}

		if self.row_rotation {
			indices.insert(self.sk.par.degree() * 2 - 1);
		}

		if self.inner_sum {
			// Add the required indices to the set of indices
			indices.insert(self.sk.par.degree() * 2 - 1);
			let mut i = 1;
			while i < self.sk.par.degree() / 2 {
				indices.insert(*ek.rot_to_gk_exponent.get(&i).unwrap());
				i *= 2
			}
		}

		for l in 0..self.expansion_level {
			indices.insert((self.sk.par.degree() >> l) + 1);
		}

		for l in 0..self.sk.par.degree().ilog2() {
			let mut monomial = vec![0i64; self.sk.par.degree()];
			monomial[self.sk.par.degree() - (1 << l)] = -1;
			let mut monomial = Poly::try_convert_from(
				&monomial as &[i64],
				&self.sk.par.ctx,
				true,
				Representation::PowerBasis,
			)?;
			unsafe { monomial.allow_variable_time_computations() }
			monomial.change_representation(Representation::Ntt);
			ek.monomials.push(monomial);
		}

		for index in indices {
			ek.gk.insert(index, GaloisKey::new(&self.sk, index)?);
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
	fn try_convert_from(value: &EvaluationKeyProto, par: &Arc<BfvParameters>) -> Result<Self> {
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
			gk.insert(key.exponent, key);
		}

		let mut monomials = Vec::with_capacity(par.degree().ilog2() as usize);
		for l in 0..par.degree().ilog2() {
			let mut monomial = vec![0i64; par.degree()];
			monomial[par.degree() - (1 << l)] = -1;
			let mut monomial = Poly::try_convert_from(
				&monomial as &[i64],
				&par.ctx,
				true,
				Representation::PowerBasis,
			)?;
			unsafe { monomial.allow_variable_time_computations() }
			monomial.change_representation(Representation::Ntt);
			monomials.push(monomial);
		}

		Ok(EvaluationKey {
			rk,
			gk,
			par: par.clone(),
			rot_to_gk_exponent: EvaluationKey::construct_rot_to_gk_exponent(par),
			monomials,
		})
	}
}

#[cfg(test)]
mod tests {
	use super::{EvaluationKey, EvaluationKeyBuilder};
	use crate::bfv::{
		proto::bfv::EvaluationKey as EvaluationKeyProto,
		traits::{Decoder, Decryptor, Encoder, Encryptor, TryConvertFrom},
		BfvParameters, Encoding, Plaintext, SecretKey,
	};
	use fhers_traits::{DeserializeWithContext, Serialize};
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
		for params in [Arc::new(BfvParameters::default(2))] {
			for _ in 0..50 {
				let mut sk = SecretKey::random(&params);
				let ek = EvaluationKeyBuilder::new(&sk).enable_inner_sum()?.build()?;

				let v = params.plaintext.random_vec(params.degree());
				let expected = params
					.plaintext
					.reduce_u128(v.iter().map(|vi| *vi as u128).sum());

				let pt = Plaintext::try_encode(&v as &[u64], Encoding::Simd, &params)?;
				let ct = sk.encrypt(&pt)?;

				let ct2 = ek.computes_inner_sum(&ct)?;
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
	fn test_row_rotation() -> Result<(), Box<dyn Error>> {
		for params in [Arc::new(BfvParameters::default(2))] {
			for _ in 0..50 {
				let mut sk = SecretKey::random(&params);
				let ek = EvaluationKeyBuilder::new(&sk)
					.enable_row_rotation()?
					.build()?;

				let v = params.plaintext.random_vec(params.degree());
				let row_size = params.degree() >> 1;
				let mut expected = vec![0u64; params.degree()];
				expected[..row_size].copy_from_slice(&v[row_size..]);
				expected[row_size..].copy_from_slice(&v[..row_size]);

				let pt = Plaintext::try_encode(&v as &[u64], Encoding::Simd, &params)?;
				let ct = sk.encrypt(&pt)?;

				let ct2 = ek.rotates_row(&ct)?;
				let pt = sk.decrypt(&ct2)?;
				assert_eq!(Vec::<u64>::try_decode(&pt, Encoding::Simd)?, expected)
			}
		}
		Ok(())
	}

	#[test]
	fn test_column_rotation() -> Result<(), Box<dyn Error>> {
		for params in [Arc::new(BfvParameters::default(2))] {
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
					let ct = sk.encrypt(&pt)?;

					let ct2 = ek.rotates_column_by(&ct, i)?;
					let pt = sk.decrypt(&ct2)?;
					assert_eq!(Vec::<u64>::try_decode(&pt, Encoding::Simd)?, expected)
				}
			}
		}
		Ok(())
	}

	#[test]
	fn test_expansion() -> Result<(), Box<dyn Error>> {
		for params in [Arc::new(BfvParameters::default(2))] {
			let log_degree = 64 - 1 - params.degree().leading_zeros();
			for _ in 0..1 {
				for i in 1..1 + log_degree as usize {
					let mut sk = SecretKey::random(&params);
					let ek = EvaluationKeyBuilder::new(&sk)
						.enable_expansion(i)?
						.build()?;

					assert!(ek.supports_expansion(i));
					assert!(!ek.supports_expansion(i + 1));
					let v = params.plaintext.random_vec(1 << i);
					let pt = Plaintext::try_encode(&v as &[u64], Encoding::Poly, &params)?;
					let ct = sk.encrypt(&pt)?;

					let ct2 = ek.expands(&ct, i)?;
					assert_eq!(ct2.len(), 1 << i);
					for (vi, ct2i) in izip!(&v, &ct2) {
						let mut expected = vec![0u64; params.degree()];
						expected[0] = params.plaintext.mul(*vi, (1 << i) as u64);
						let pt = sk.decrypt(ct2i)?;
						assert_eq!(expected, Vec::<u64>::try_decode(&pt, Encoding::Poly)?);
					}
				}
			}
		}
		Ok(())
	}

	#[test]
	fn test_proto_conversion() -> Result<(), Box<dyn Error>> {
		for params in [
			Arc::new(BfvParameters::default(1)),
			Arc::new(BfvParameters::default(2)),
		] {
			let sk = SecretKey::random(&params);

			let ek = EvaluationKeyBuilder::new(&sk).build()?;
			let proto = EvaluationKeyProto::from(&ek);
			assert_eq!(ek, EvaluationKey::try_convert_from(&proto, &params)?);

			if params.ciphertext_moduli.len() > 1 {
				let ek = EvaluationKeyBuilder::new(&sk)
					.enable_row_rotation()?
					.build()?;
				let proto = EvaluationKeyProto::from(&ek);
				assert_eq!(ek, EvaluationKey::try_convert_from(&proto, &params)?);

				let ek = EvaluationKeyBuilder::new(&sk)
					.enable_inner_sum()?
					.enable_relinearization()?
					.build()?;
				let proto = EvaluationKeyProto::from(&ek);
				assert_eq!(ek, EvaluationKey::try_convert_from(&proto, &params)?);

				let ek = EvaluationKeyBuilder::new(&sk)
					.enable_expansion(params.degree().ilog2() as usize)?
					.build()?;
				let proto = EvaluationKeyProto::from(&ek);
				assert_eq!(ek, EvaluationKey::try_convert_from(&proto, &params)?);

				let ek = EvaluationKeyBuilder::new(&sk)
					.enable_inner_sum()?
					.enable_relinearization()?
					.enable_expansion(params.degree().ilog2() as usize)?
					.build()?;
				let proto = EvaluationKeyProto::from(&ek);
				assert_eq!(ek, EvaluationKey::try_convert_from(&proto, &params)?);
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
			let bytes = ek.serialize();
			assert_eq!(ek, EvaluationKey::try_deserialize(&bytes, &params)?);

			if params.ciphertext_moduli.len() > 1 {
				let ek = EvaluationKeyBuilder::new(&sk)
					.enable_row_rotation()?
					.build()?;
				let bytes = ek.serialize();
				assert_eq!(ek, EvaluationKey::try_deserialize(&bytes, &params)?);

				let ek = EvaluationKeyBuilder::new(&sk)
					.enable_inner_sum()?
					.enable_relinearization()?
					.build()?;
				let bytes = ek.serialize();
				assert_eq!(ek, EvaluationKey::try_deserialize(&bytes, &params)?);

				let ek = EvaluationKeyBuilder::new(&sk)
					.enable_expansion(params.degree().ilog2() as usize)?
					.build()?;
				let bytes = ek.serialize();
				assert_eq!(ek, EvaluationKey::try_deserialize(&bytes, &params)?);

				let ek = EvaluationKeyBuilder::new(&sk)
					.enable_inner_sum()?
					.enable_relinearization()?
					.enable_expansion(params.degree().ilog2() as usize)?
					.build()?;
				let bytes = ek.serialize();
				assert_eq!(ek, EvaluationKey::try_deserialize(&bytes, &params)?);
			}
		}
		Ok(())
	}
}
