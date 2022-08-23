use std::cmp::min;

use itertools::{izip, Itertools};
use math::rq::{dot_product as poly_dot_product, Representation};

use crate::{
	bfv::{
		advanced::LeveledEvaluationKey, parameters::MultiplicationParameters, Ciphertext, Plaintext,
	},
	Error, Result,
};

/// Multiply two ciphertext and relinearize.
fn mul_internal(
	ct0: &Ciphertext,
	ct1: &Ciphertext,
	ek: &LeveledEvaluationKey,
	mp: &MultiplicationParameters,
) -> Result<Ciphertext> {
	if !ek.supports_relinearization() {
		return Err(Error::DefaultError(
			"The evaluation key does not support relinearization".to_string(),
		));
	}
	if ct0.par != ct1.par {
		return Err(Error::DefaultError("Incompatible parameters".to_string()));
	}
	if ct0.level != ct1.level {
		return Err(Error::DefaultError("Incompatible levels".to_string()));
	}
	if ct0.par.moduli.len() == 1 {
		return Err(Error::DefaultError(
			"Parameters do not allow for multiplication".to_string(),
		));
	}
	if ct0.c.len() != 2 || ct1.c.len() != 2 {
		return Err(Error::DefaultError(
			"Multiplication can only be performed on ciphertexts of size 2".to_string(),
		));
	}

	// Extend
	// let mut now = std::time::SystemTime::now();
	let c00 = ct0.c[0].scale(&mp.extender_self)?;
	let c01 = ct0.c[1].scale(&mp.extender_self)?;
	let c10 = ct1.c[0].scale(&mp.extender_other)?;
	let c11 = ct1.c[1].scale(&mp.extender_other)?;
	// println!("Extend: {:?}", now.elapsed().unwrap());

	// Multiply
	// now = std::time::SystemTime::now();
	let mut c0 = &c00 * &c10;
	let mut c1 = &c00 * &c11;
	c1 += &(&c01 * &c10);
	let mut c2 = &c01 * &c11;
	c0.change_representation(Representation::PowerBasis);
	c1.change_representation(Representation::PowerBasis);
	c2.change_representation(Representation::PowerBasis);
	// println!("Multiply: {:?}", now.elapsed().unwrap());

	// Scale
	// TODO: This should be faster??
	// now = std::time::SystemTime::now();
	let mut c0 = c0.scale(&mp.down_scaler)?;
	let mut c1 = c1.scale(&mp.down_scaler)?;
	let c2 = c2.scale(&mp.down_scaler)?;
	// println!("Scale: {:?}", now.elapsed().unwrap());

	// Relinearize
	// now = std::time::SystemTime::now();
	c0.change_representation(Representation::Ntt);
	c1.change_representation(Representation::Ntt);
	let (c0r, c1r) = ek.relinearizes_with_poly(&c2)?;
	c0 += &c0r;
	c1 += &c1r;
	// println!("Relinearize: {:?}", now.elapsed().unwrap());

	Ok(Ciphertext {
		par: ct0.par.clone(),
		seed: None,
		c: vec![c0, c1],
		level: ct0.level,
	})
}

/// Multiply two ciphertext and relinearize.
pub fn mul(ct0: &Ciphertext, ct1: &Ciphertext, ek: &LeveledEvaluationKey) -> Result<Ciphertext> {
	mul_internal(ct0, ct1, ek, &ct0.par.mul_1_params[ct0.level])
}

/// Multiply two ciphertext and relinearize.
pub fn mul2(ct0: &Ciphertext, ct1: &Ciphertext, ek: &LeveledEvaluationKey) -> Result<Ciphertext> {
	mul_internal(ct0, ct1, ek, &ct0.par.mul_1_params[ct0.level])
}

/// Compute the dot product between an iterator of [`Ciphertext`] and an
/// iterator of [`Plaintext`]. Returns an error if the iterator counts are 0, if
/// the parameters don't match, or if the ciphertexts have different
/// number of parts.
pub fn dot_product_scalar<'a, I, J>(ct: I, pt: J) -> Result<Ciphertext>
where
	I: Iterator<Item = &'a Ciphertext> + Clone,
	J: Iterator<Item = &'a Plaintext> + Clone,
{
	let count = min(ct.clone().count(), pt.clone().count());
	if count == 0 {
		return Err(Error::DefaultError(
			"At least one iterator is empty".to_string(),
		));
	}

	let ct_first = ct.clone().next().unwrap();
	if izip!(ct.clone(), pt.clone()).any(|(cti, pti)| {
		cti.par != ct_first.par || pti.par != ct_first.par || cti.c.len() != ct_first.c.len()
	}) {
		return Err(Error::DefaultError("Mismatched parameters".to_string()));
	}
	if ct.clone().any(|cti| cti.c.len() != ct_first.c.len()) {
		return Err(Error::DefaultError(
			"Mismatched number of parts in the ciphertexts".to_string(),
		));
	}

	let c = (0..ct_first.c.len())
		.map(|i| {
			poly_dot_product(
				ct.clone().map(|cti| unsafe { cti.c.get_unchecked(i) }),
				pt.clone().map(|pti| &pti.poly_ntt),
			)
			.unwrap()
		})
		.collect_vec();

	Ok(Ciphertext {
		par: ct_first.par.clone(),
		seed: None,
		c,
		level: ct_first.level,
	})
}

#[cfg(test)]
mod tests {
	use super::{dot_product_scalar, mul, mul2};
	use crate::bfv::{
		advanced::LeveledEvaluationKeyBuilder, BfvParameters, Ciphertext, Encoding, Plaintext,
		SecretKey,
	};
	use fhers_traits::{FheDecoder, FheDecrypter, FheEncoder, FheEncrypter};
	use itertools::{izip, Itertools};
	use std::{error::Error, sync::Arc};

	#[test]
	fn test_mul() -> Result<(), Box<dyn Error>> {
		let par = Arc::new(BfvParameters::default(2));
		for _ in 0..50 {
			// We will encode `values` in an Simd format, and check that the product is
			// computed correctly.
			let values = par.plaintext.random_vec(par.degree());
			let mut expected = values.clone();
			par.plaintext.mul_vec(&mut expected, &values);

			let sk = SecretKey::random(&par);
			let ek = LeveledEvaluationKeyBuilder::new(&sk, 0, 0)?
				.enable_relinearization()?
				.build()?;
			let pt = Plaintext::try_encode(&values as &[u64], Encoding::SimdLeveled(0), &par)?;

			let ct1 = sk.try_encrypt(&pt)?;
			let ct2 = sk.try_encrypt(&pt)?;
			let ct3 = mul(&ct1, &ct2, &ek)?;

			println!("Noise: {}", unsafe { sk.measure_noise(&ct3)? });
			let pt = sk.try_decrypt(&ct3)?;
			assert_eq!(
				Vec::<u64>::try_decode(&pt, Encoding::SimdLeveled(0))?,
				expected
			);
		}
		Ok(())
	}

	#[test]
	fn test_mul2() -> Result<(), Box<dyn Error>> {
		let ntests = 100;
		let par = Arc::new(BfvParameters::default(2));
		for _ in 0..ntests {
			// We will encode `values` in an Simd format, and check that the product is
			// computed correctly.
			let values = par.plaintext.random_vec(par.degree());
			let mut expected = values.clone();
			par.plaintext.mul_vec(&mut expected, &values);

			let sk = SecretKey::random(&par);
			let ek = LeveledEvaluationKeyBuilder::new(&sk, 0, 0)?
				.enable_relinearization()?
				.build()?;
			let pt = Plaintext::try_encode(&values as &[u64], Encoding::SimdLeveled(0), &par)?;

			let ct1 = sk.try_encrypt(&pt)?;
			let ct2 = sk.try_encrypt(&pt)?;
			let ct3 = mul2(&ct1, &ct2, &ek)?;

			println!("Noise: {}", unsafe { sk.measure_noise(&ct3)? });
			let pt = sk.try_decrypt(&ct3)?;
			assert_eq!(
				Vec::<u64>::try_decode(&pt, Encoding::SimdLeveled(0))?,
				expected
			);
		}
		Ok(())
	}

	#[test]
	fn test_dot_product_scalar() -> Result<(), Box<dyn Error>> {
		for params in [
			Arc::new(BfvParameters::default(1)),
			Arc::new(BfvParameters::default(2)),
		] {
			let sk = SecretKey::random(&params);
			for size in 1..128 {
				let ct = (0..size)
					.map(|_| {
						let v = params.plaintext.random_vec(params.degree());
						let pt =
							Plaintext::try_encode(&v as &[u64], Encoding::SimdLeveled(0), &params)
								.unwrap();
						sk.try_encrypt(&pt).unwrap()
					})
					.collect_vec();
				let pt = (0..size)
					.map(|_| {
						let v = params.plaintext.random_vec(params.degree());
						Plaintext::try_encode(&v as &[u64], Encoding::SimdLeveled(0), &params)
							.unwrap()
					})
					.collect_vec();

				let r = dot_product_scalar(ct.iter(), pt.iter())?;

				let mut expected = Ciphertext::zero(&params);
				izip!(&ct, &pt).for_each(|(cti, pti)| expected += cti * pti);
				assert_eq!(r, expected);
			}
		}
		Ok(())
	}
}
