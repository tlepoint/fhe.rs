//! Plaintext type in the BFV encryption scheme.
use crate::parameters::BfvParameters;
use crate::traits::{Decoder, Encoder};
use math::rq::{traits::TryConvertFrom, Context, Poly, Representation};
use std::cmp::min;
use std::sync::Arc;
use zeroize::{Zeroize, ZeroizeOnDrop};

/// An encoding for the plaintext.
#[derive(Debug, Clone, Eq, PartialEq)]
pub enum Encoding {
	/// A Poly encoding encodes a vector as coefficients of a polynomial;
	/// homomorphic operations are therefore polynomial operations.
	Poly,
	/// A Simd encoding encodes a vector so that homomorphic operations are
	/// component-wise operations on the coefficients of the underlying vectors. The
	/// Simd encoding require that the plaintext modulus is congruent to 1 modulo
	/// the degree of the underlying polynomial.
	Simd,
}

/// A plaintext object, that encodes a vector according to a specific encoding.
#[derive(Debug, Clone)]
pub struct Plaintext {
	/// The parameters of the underlying BFV encryption scheme.
	pub(crate) par: Arc<BfvParameters>,
	/// The value after encoding.
	pub(crate) value: Vec<i64>,
	/// The encoding of the plaintext, if known
	pub(crate) encoding: Option<Encoding>,
	/// The plaintext as a polynomial.
	pub(crate) poly_ntt: Poly,
}

impl Plaintext {
	pub(crate) fn encode(&self) -> Result<Poly, String> {
		let mut m_v = self.par.plaintext.reduce_vec_i64(&self.value);

		self.par
			.plaintext
			.scalar_mul_vec(&mut m_v, self.par.q_mod_t);
		let mut m = Poly::try_convert_from(&m_v, &self.par.ctx, Representation::PowerBasis)?;
		m.change_representation(Representation::Ntt);
		m *= &self.par.delta;
		Ok(m)
	}

	/// Generate a zero plaintext.
	pub fn zero(encoding: Encoding, par: &Arc<BfvParameters>) -> Self {
		let value = vec![0i64; par.degree()];
		let poly_ntt = Poly::zero(&par.ctx, Representation::Ntt);
		Self {
			par: par.clone(),
			value,
			encoding: Some(encoding),
			poly_ntt,
		}
	}
}

unsafe impl Send for Plaintext {}

// Equality.
impl PartialEq for Plaintext {
	fn eq(&self, other: &Self) -> bool {
		let mut eq = self.par == other.par && self.value == other.value;
		if self.encoding.is_some() && other.encoding.is_some() {
			eq &= self.encoding.as_ref().unwrap() == other.encoding.as_ref().unwrap()
		}
		eq
	}
}

impl Eq for Plaintext {}

// Conversions.
impl TryConvertFrom<&Plaintext> for Poly {
	type Error = String;

	fn try_convert_from<R>(pt: &Plaintext, ctx: &Arc<Context>, _: R) -> Result<Self, Self::Error>
	where
		R: Into<Option<Representation>>,
	{
		if ctx != &pt.par.ctx {
			Err("Incompatible contexts".to_string())
		} else {
			Poly::try_convert_from(&pt.value as &[i64], &pt.par.ctx, Representation::PowerBasis)
		}
	}
}

// Zeroizing of plaintexts.

impl Zeroize for Plaintext {
	fn zeroize(&mut self) {
		self.value.zeroize()
	}
}

impl ZeroizeOnDrop for Plaintext {}

// Encoding and decoding.

impl Encoder<&[u64]> for Plaintext {
	type Error = String;

	fn try_encode(
		value: &[u64],
		encoding: Encoding,
		par: &Arc<BfvParameters>,
	) -> Result<Self, Self::Error> {
		if value.len() > par.degree() {
			return Err("There are too many values to encode".to_string());
		}
		let v = Vec::<Plaintext>::try_encode(value, encoding, par)?;
		Ok(v[0].clone())
	}
}

impl Encoder<&[u64]> for Vec<Plaintext> {
	type Error = String;

	fn try_encode(
		value: &[u64],
		encoding: Encoding,
		par: &Arc<BfvParameters>,
	) -> Result<Self, Self::Error> {
		if value.is_empty() {
			return Ok(vec![Plaintext::zero(encoding, par)]);
		}
		let encoding = &encoding;

		let num_plaintexts = value.len().div_ceil(par.degree());
		let mut out = Vec::with_capacity(num_plaintexts);
		for i in 0..num_plaintexts {
			let slice = &value[i * par.degree()..min((i + 1) * par.degree(), value.len())];
			let mut v = vec![0u64; par.degree()];

			match *encoding {
				Encoding::Poly => v[..slice.len()].copy_from_slice(slice),
				Encoding::Simd => {
					if let Some(op) = &par.op {
						for i in 0..slice.len() {
							v[par.matrix_reps_index_map[i]] = slice[i];
						}
						op.backward(&mut v);
					} else {
						return Err(
							"The plaintext does not allow for using the Simd encoding".to_string()
						);
					}
				}
			}

			let w = unsafe { par.plaintext.center_vec_vt(&v) };
			v.zeroize();

			let mut poly =
				Poly::try_convert_from(&w as &[i64], &par.ctx, Representation::PowerBasis)?;
			poly.change_representation(Representation::Ntt);

			out.push(Plaintext {
				par: par.clone(),
				value: w,
				encoding: Some(encoding.clone()),
				poly_ntt: poly,
			})
		}
		Ok(out)
	}
}

impl Encoder<&[i64]> for Plaintext {
	type Error = String;

	fn try_encode(
		value: &[i64],
		encoding: Encoding,
		par: &Arc<BfvParameters>,
	) -> Result<Self, Self::Error> {
		if value.len() > par.degree() {
			return Err("There are too many values to encode".to_string());
		}
		let mut v = vec![0u64; par.degree()];

		if !value.is_empty() {
			let mut value_u64 = par.plaintext.reduce_vec_i64(value);
			match encoding {
				Encoding::Poly => v[..value_u64.len()].copy_from_slice(&value_u64),
				Encoding::Simd => {
					if let Some(op) = &par.op {
						for i in 0..value_u64.len() {
							v[par.matrix_reps_index_map[i]] = value_u64[i];
						}
						op.backward(&mut v);
					} else {
						return Err(
							"The plaintext does not allow for using the Simd encoding".to_string()
						);
					}
				}
			}
			value_u64.zeroize()
		}

		let w = unsafe { par.plaintext.center_vec_vt(&v) };
		v.zeroize();

		let mut poly = Poly::try_convert_from(&w as &[i64], &par.ctx, Representation::PowerBasis)?;
		poly.change_representation(Representation::Ntt);

		Ok(Self {
			par: par.clone(),
			value: w,
			encoding: Some(encoding),
			poly_ntt: poly,
		})
	}
}

impl Decoder for Vec<u64> {
	type Error = String;

	fn try_decode<E>(pt: &Plaintext, encoding: E) -> Result<Vec<u64>, Self::Error>
	where
		E: Into<Option<Encoding>>,
	{
		let encoding = encoding.into();
		let enc: Encoding;
		if pt.encoding.is_none() && encoding.is_none() {
			return Err("No encoding specified".to_string());
		} else if pt.encoding.is_some() {
			enc = pt.encoding.as_ref().unwrap().clone();
			if let Some(arg_enc) = encoding && arg_enc != enc {
				return Err("Mismatched encodings".to_string())
			}
		} else {
			enc = encoding.unwrap();
			if let Some(pt_enc) = pt.encoding.as_ref() && pt_enc != &enc {
				return Err("Mismatched encodings".to_string())
			}
		}

		let mut w = pt.par.plaintext.reduce_vec_i64(&pt.value);

		match enc {
			Encoding::Poly => Ok(w),
			Encoding::Simd => {
				if let Some(op) = &pt.par.op {
					op.forward(&mut w);
					let mut w_reordered = w.clone();
					for i in 0..pt.par.degree() {
						w_reordered[i] = w[pt.par.matrix_reps_index_map[i]]
					}
					w.zeroize();
					Ok(w_reordered)
				} else {
					Err("The plaintext does not allow for using the Simd encoding".to_string())
				}
			}
		}
	}
}

impl Decoder for Vec<i64> {
	type Error = String;

	fn try_decode<E>(pt: &Plaintext, encoding: E) -> Result<Vec<i64>, Self::Error>
	where
		E: Into<Option<Encoding>>,
	{
		let v = Vec::<u64>::try_decode(pt, encoding)?;
		Ok(unsafe { pt.par.plaintext.center_vec_vt(&v) })
	}
}

#[cfg(test)]
mod tests {
	use super::{Encoding, Plaintext};
	use crate::parameters::{BfvParameters, BfvParametersBuilder};
	use crate::traits::{Decoder, Encoder};
	use std::sync::Arc;

	#[test]
	fn try_encode() -> Result<(), String> {
		// The default test parameters support both Poly and Simd encodings
		let params = Arc::new(BfvParameters::default(1));
		let a = params.plaintext.random_vec(params.degree());

		let plaintext = Plaintext::try_encode(&[0u64; 9] as &[u64], Encoding::Poly, &params);
		assert!(plaintext.is_err());

		let plaintext = Plaintext::try_encode(&a as &[u64], Encoding::Poly, &params);
		assert!(plaintext.is_ok());

		let plaintext = Plaintext::try_encode(&a as &[u64], Encoding::Simd, &params);
		assert!(plaintext.is_ok());

		let plaintext = Plaintext::try_encode(&[1u64] as &[u64], Encoding::Poly, &params);
		assert!(plaintext.is_ok());

		// The following parameters do not allow for Simd encoding
		let params = Arc::new(
			BfvParametersBuilder::new()
				.set_degree(8)?
				.set_plaintext_modulus(2)?
				.set_ciphertext_moduli(&[4611686018326724609])?
				.build()?,
		);

		let a = params.plaintext.random_vec(params.degree());

		let plaintext = Plaintext::try_encode(&a as &[u64], Encoding::Poly, &params);
		assert!(plaintext.is_ok());

		let plaintext = Plaintext::try_encode(&a as &[u64], Encoding::Simd, &params);
		assert!(plaintext.is_err());

		Ok(())
	}

	#[test]
	fn test_encode_decode() {
		(0..40).for_each(|_| {
			for i in 1..5 {
				let params = Arc::new(BfvParameters::default(1));
				let a = params.plaintext.random_vec(params.degree() * i);

				let plaintexts =
					Vec::<Plaintext>::try_encode(&a as &[u64], Encoding::Poly, &params);
				assert!(plaintexts.is_ok());
				let plaintexts = plaintexts.unwrap();
				assert_eq!(plaintexts.len(), i);
				for j in 0..i {
					let b = Vec::<u64>::try_decode(&plaintexts[j], Encoding::Poly);
					assert!(
						b.is_ok_and(|b| b == &a[j * params.degree()..(j + 1) * params.degree()])
					);
				}

				let a = a[..params.degree()].to_vec();

				let plaintext = Plaintext::try_encode(&a as &[u64], Encoding::Simd, &params);
				assert!(plaintext.is_ok());
				let b = Vec::<u64>::try_decode(&plaintext.unwrap(), Encoding::Simd);
				assert!(b.is_ok_and(|b| b == &a));

				let a = unsafe { params.plaintext.center_vec_vt(&a) };
				let plaintext = Plaintext::try_encode(&a as &[i64], Encoding::Poly, &params);
				assert!(plaintext.is_ok());
				let b = Vec::<i64>::try_decode(&plaintext.unwrap(), Encoding::Poly);
				assert!(b.is_ok_and(|b| b == &a));

				let plaintext = Plaintext::try_encode(&a as &[i64], Encoding::Simd, &params);
				assert!(plaintext.is_ok());
				let b = Vec::<i64>::try_decode(&plaintext.unwrap(), Encoding::Simd);
				assert!(b.is_ok_and(|b| b == &a));
			}
		})
	}

	#[test]
	fn test_partial_eq() -> Result<(), String> {
		let params = Arc::new(BfvParameters::default(1));
		let a = params.plaintext.random_vec(params.degree());

		let plaintext = Plaintext::try_encode(&a as &[u64], Encoding::Poly, &params)?;
		let mut same_plaintext = Plaintext::try_encode(&a as &[u64], Encoding::Poly, &params)?;
		assert_eq!(plaintext, same_plaintext);

		// Equality also holds when there is no encoding specified. In this test, we use the fact that
		// we can set it to None directly, but such a partial plaintext will be created during
		// decryption since we do not specify the encoding at the time.
		same_plaintext.encoding = None;
		assert_eq!(plaintext, same_plaintext);

		Ok(())
	}
}
