use std::cmp::min;

use itertools::{izip, Itertools};
use math::rq::dot_product as poly_dot_product;

use crate::{
	bfv::{Ciphertext, Plaintext},
	Error, Result,
};

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
	use super::dot_product_scalar;
	use crate::bfv::{BfvParameters, Ciphertext, Encoding, Plaintext, SecretKey};
	use fhers_traits::{FheEncoder, FheEncrypter};
	use itertools::{izip, Itertools};
	use std::{error::Error, sync::Arc};

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
							Plaintext::try_encode(&v as &[u64], Encoding::simd(), &params).unwrap();
						sk.try_encrypt(&pt).unwrap()
					})
					.collect_vec();
				let pt = (0..size)
					.map(|_| {
						let v = params.plaintext.random_vec(params.degree());
						Plaintext::try_encode(&v as &[u64], Encoding::simd(), &params).unwrap()
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
