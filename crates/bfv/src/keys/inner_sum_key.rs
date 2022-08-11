//! Inner Sum keys for the BFV encryption scheme

use std::{collections::HashMap, rc::Rc};

use crate::{traits::TryConvertFrom, BfvParameters, Ciphertext, GaloisKey, SecretKey};
use fhers_protos::protos::bfv::{GaloisKey as GaloisKeyProto, InnerSumKey as InnerSumKeyProto};
use math::zq::Modulus;

/// Inner-sum key for the BFV encryption scheme, which enables to
/// compute the sum of all the slots.
#[derive(Debug, PartialEq, Eq)]
pub struct InnerSumKey {
	gk: HashMap<usize, GaloisKey>,
}

impl InnerSumKey {
	/// Generate a [`GaloisKey`] from a [`SecretKey`].
	pub fn new(sk: &SecretKey) -> Result<Self, String> {
		let q = Modulus::new(2 * sk.par.degree() as u64)?;
		let mut gk = HashMap::new();
		let mut i = 1;
		while i < sk.par.degree() {
			let exp = q.pow(3, i as u64);
			gk.insert(i, GaloisKey::new(sk, exp as usize)?);
			i *= 2
		}
		gk.insert(0, GaloisKey::new(sk, sk.par.degree() * 2 - 1)?);
		Ok(Self { gk })
	}

	/// Compute the inner sum of a [`Ciphertext`]
	pub fn inner_sum(&self, ct: &Ciphertext) -> Result<Ciphertext, String> {
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

impl From<&InnerSumKey> for InnerSumKeyProto {
	fn from(value: &InnerSumKey) -> Self {
		let mut isk = InnerSumKeyProto::new();
		for (_, gk) in value.gk.iter() {
			isk.gk.push(GaloisKeyProto::from(gk));
		}
		isk
	}
}

impl TryConvertFrom<&InnerSumKeyProto> for InnerSumKey {
	type Error = String;

	fn try_convert_from(
		value: &InnerSumKeyProto,
		par: &Rc<BfvParameters>,
	) -> Result<Self, Self::Error> {
		let mut gkmap = HashMap::new();
		let q = Modulus::new(2 * par.degree() as u64)?;

		for gkp in value.gk.iter() {
			let gk = GaloisKey::try_convert_from(gkp, par)?;

			if gk.exponent == 2 * par.degree() - 1 {
				gkmap.insert(0, gk);
			} else {
				let mut i = 1;
				while i < par.degree() {
					let exp = q.pow(3, i as u64);
					if exp as usize == gk.exponent {
						gkmap.insert(i, gk);
						break;
					}
					i *= 2
				}
			}
		}

		// Verify that we have all the necessary keys
		let mut all_keys_present = gkmap.contains_key(&0);
		let mut i = 1;
		while i < par.degree() {
			all_keys_present &= gkmap.contains_key(&i);
			i *= 2;
		}

		if !all_keys_present {
			Err("Some Galois keys are missing".to_string())
		} else {
			Ok(Self { gk: gkmap })
		}
	}
}

#[cfg(test)]
mod tests {
	use crate::{
		traits::{Decoder, Decryptor, Encoder, Encryptor, TryConvertFrom},
		BfvParameters, Encoding, InnerSumKey, Plaintext, SecretKey,
	};
	use fhers_protos::protos::bfv::InnerSumKey as InnerSumKeyProto;
	use std::rc::Rc;

	#[test]
	fn test_inner_sum() -> Result<(), String> {
		for params in [Rc::new(BfvParameters::default(2))] {
			for _ in 0..50 {
				let sk = SecretKey::random(&params);
				let isk = InnerSumKey::new(&sk)?;

				let v = params.plaintext.random_vec(params.degree());
				let expected = params
					.plaintext
					.reduce_u128(v.iter().map(|vi| *vi as u128).sum());

				let pt = Plaintext::try_encode(&v as &[u64], Encoding::Simd, &params)?;
				let mut ct = sk.encrypt(&pt)?;

				let ct2 = isk.inner_sum(&mut ct)?;
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
	fn test_proto_conversion() -> Result<(), String> {
		for params in [
			Rc::new(BfvParameters::default(1)),
			Rc::new(BfvParameters::default(2)),
		] {
			let sk = SecretKey::random(&params);
			let isk = InnerSumKey::new(&sk)?;
			let proto = InnerSumKeyProto::from(&isk);
			assert_eq!(isk, InnerSumKey::try_convert_from(&proto, &params)?);
		}
		Ok(())
	}
}
