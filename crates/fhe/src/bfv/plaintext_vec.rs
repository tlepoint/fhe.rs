use std::{cmp::min, sync::Arc};

use fhe_math::rq::{traits::TryConvertFrom, Poly, Representation};
use fhe_traits::{FheEncoder, FheEncoderVariableTime, FheParametrized, FhePlaintext};
use fhe_util::div_ceil;
use zeroize::{Zeroize, ZeroizeOnDrop};

use crate::{
	bfv::{BfvParameters, Encoding, Plaintext},
	Error, Result,
};

use super::encoding::EncodingEnum;

/// A wrapper around a vector of plaintext which implements the [`FhePlaintext`]
/// trait, and therefore can be encoded to / decoded from.
pub struct PlaintextVec(pub Vec<Plaintext>);

impl FhePlaintext for PlaintextVec {
	type Encoding = Encoding;
}

impl FheParametrized for PlaintextVec {
	type Parameters = BfvParameters;
}

impl Zeroize for PlaintextVec {
	fn zeroize(&mut self) {
		self.0.zeroize()
	}
}

impl ZeroizeOnDrop for PlaintextVec {}

impl FheEncoderVariableTime<&[u64]> for PlaintextVec {
	type Error = Error;

	unsafe fn try_encode_vt(
		value: &[u64],
		encoding: Encoding,
		par: &Arc<BfvParameters>,
	) -> Result<Self> {
		if value.is_empty() {
			return Ok(PlaintextVec(vec![Plaintext::zero(encoding, par)?]));
		}
		if encoding.encoding == EncodingEnum::Simd && par.op.is_none() {
			return Err(Error::EncodingNotSupported(EncodingEnum::Simd.to_string()));
		}
		let ctx = par.ctx_at_level(encoding.level)?;
		let num_plaintexts = div_ceil(value.len(), par.degree());

		Ok(PlaintextVec(
			(0..num_plaintexts)
				.map(|i| {
					let slice = &value[i * par.degree()..min(value.len(), (i + 1) * par.degree())];
					let mut v = vec![0u64; par.degree()];
					match encoding.encoding {
						EncodingEnum::Poly => v[..slice.len()].copy_from_slice(slice),
						EncodingEnum::Simd => {
							for i in 0..slice.len() {
								v[par.matrix_reps_index_map[i]] = slice[i];
							}
							par.op.as_ref().unwrap().backward_vt(v.as_mut_ptr());
						}
					};

					let mut poly = Poly::try_convert_from(
						&v as &[u64],
						ctx,
						true,
						Representation::PowerBasis,
					)?;
					poly.change_representation(Representation::Ntt);

					Ok(Plaintext {
						par: par.clone(),
						value: v.into_boxed_slice(),
						encoding: Some(encoding.clone()),
						poly_ntt: poly,
						level: encoding.level,
					})
				})
				.collect::<Result<Vec<Plaintext>>>()?,
		))
	}
}

impl FheEncoder<&[u64]> for PlaintextVec {
	type Error = Error;
	fn try_encode(value: &[u64], encoding: Encoding, par: &Arc<BfvParameters>) -> Result<Self> {
		if value.is_empty() {
			return Ok(PlaintextVec(vec![Plaintext::zero(encoding, par)?]));
		}
		if encoding.encoding == EncodingEnum::Simd && par.op.is_none() {
			return Err(Error::EncodingNotSupported(EncodingEnum::Simd.to_string()));
		}
		let ctx = par.ctx_at_level(encoding.level)?;
		let num_plaintexts = div_ceil(value.len(), par.degree());

		Ok(PlaintextVec(
			(0..num_plaintexts)
				.map(|i| {
					let slice = &value[i * par.degree()..min(value.len(), (i + 1) * par.degree())];
					let mut v = vec![0u64; par.degree()];
					match encoding.encoding {
						EncodingEnum::Poly => v[..slice.len()].copy_from_slice(slice),
						EncodingEnum::Simd => {
							for i in 0..slice.len() {
								v[par.matrix_reps_index_map[i]] = slice[i];
							}
							par.op.as_ref().unwrap().backward(&mut v);
						}
					};

					let mut poly = Poly::try_convert_from(
						&v as &[u64],
						ctx,
						false,
						Representation::PowerBasis,
					)?;
					poly.change_representation(Representation::Ntt);

					Ok(Plaintext {
						par: par.clone(),
						value: v.into_boxed_slice(),
						encoding: Some(encoding.clone()),
						poly_ntt: poly,
						level: encoding.level,
					})
				})
				.collect::<Result<Vec<Plaintext>>>()?,
		))
	}
}

#[cfg(test)]
mod tests {
	use std::{error::Error, sync::Arc};

	use fhe_traits::{FheDecoder, FheEncoder, FheEncoderVariableTime};

	use crate::bfv::{BfvParameters, Encoding, PlaintextVec};

	#[test]
	fn encode_decode() -> Result<(), Box<dyn Error>> {
		for _ in 0..20 {
			for i in 1..5 {
				let params = Arc::new(BfvParameters::default(1, 8));
				let a = params.plaintext.random_vec(params.degree() * i);

				let plaintexts =
					PlaintextVec::try_encode(&a as &[u64], Encoding::poly_at_level(0), &params)?;
				assert_eq!(plaintexts.0.len(), i);

				for j in 0..i {
					let b = Vec::<u64>::try_decode(&plaintexts.0[j], Encoding::poly_at_level(0));
					assert!(
						b.is_ok_and(|b| b == &a[j * params.degree()..(j + 1) * params.degree()])
					);
				}

				let plaintexts = unsafe {
					PlaintextVec::try_encode_vt(&a as &[u64], Encoding::poly_at_level(0), &params)?
				};
				assert_eq!(plaintexts.0.len(), i);

				for j in 0..i {
					let b = Vec::<u64>::try_decode(&plaintexts.0[j], Encoding::poly_at_level(0));
					assert!(
						b.is_ok_and(|b| b == &a[j * params.degree()..(j + 1) * params.degree()])
					);
				}
			}
		}
		Ok(())
	}
}
