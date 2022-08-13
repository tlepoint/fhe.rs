//! Relinearization keys for the BFV encryption scheme

use std::sync::Arc;

use super::key_switching_key::KeySwitchingKey;
use crate::{traits::TryConvertFrom, BfvParameters, SecretKey};
use fhers_protos::protos::bfv::{
	KeySwitchingKey as KeySwitchingKeyProto, RelinearizationKey as RelinearizationKeyProto,
};
use math::rq::{Poly, Representation};
use protobuf::MessageField;
use zeroize::Zeroize;

/// Relinearization key for the BFV encryption scheme.
/// A relinearization key is a special type of key switching key,
/// which switch from `s^2` to `s` where `s` is the secret key.
#[derive(Debug, PartialEq, Eq)]
pub struct RelinearizationKey {
	ksk: KeySwitchingKey,
}

impl RelinearizationKey {
	/// Generate a [`RelinearizationKey`] from a [`SecretKey`].
	pub fn new(sk: &SecretKey) -> Result<Self, String> {
		let mut s_squared = sk.s[1].clone();
		s_squared.change_representation(Representation::PowerBasis);
		let ksk = KeySwitchingKey::new(sk, &s_squared)?;
		s_squared.zeroize();
		Ok(Self { ksk })
	}

	/// Relinearize an "extended" ciphertext (c0, c1, c2) into a [`Ciphertext`]
	pub fn relinearize(&self, c0: &mut Poly, c1: &mut Poly, c2: &Poly) -> Result<(), String> {
		self.ksk.key_switch(c2, c0, c1)
	}
}

impl From<&RelinearizationKey> for RelinearizationKeyProto {
	fn from(value: &RelinearizationKey) -> Self {
		let mut rk = RelinearizationKeyProto::new();
		rk.ksk = MessageField::some(KeySwitchingKeyProto::from(&value.ksk));
		rk
	}
}

impl TryConvertFrom<&RelinearizationKeyProto> for RelinearizationKey {
	type Error = String;

	fn try_convert_from(
		value: &RelinearizationKeyProto,
		par: &Arc<BfvParameters>,
	) -> Result<Self, Self::Error> {
		if value.ksk.is_some() {
			Ok(RelinearizationKey {
				ksk: KeySwitchingKey::try_convert_from(value.ksk.as_ref().unwrap(), par)?,
			})
		} else {
			Err("Invalid serialization".to_string())
		}
	}
}

#[cfg(test)]
mod tests {
	use super::RelinearizationKey;
	use crate::{
		traits::{Decoder, Decryptor, TryConvertFrom},
		BfvParameters, Ciphertext, Encoding, SecretKey,
	};
	use fhers_protos::protos::bfv::RelinearizationKey as RelinearizationKeyProto;
	use math::rq::{Poly, Representation};
	use std::sync::Arc;

	#[test]
	fn test_relinearization() -> Result<(), String> {
		for params in [Arc::new(BfvParameters::default(2))] {
			for _ in 0..100 {
				let mut sk = SecretKey::random(&params);
				let rk = RelinearizationKey::new(&sk)?;

				// Let's generate manually an "extended" ciphertext (c0 = e - c1 * s - c2 * s^2, c1, c2) encrypting 0.
				let mut c2 = Poly::random(&params.ctx, Representation::Ntt);
				let mut c1 = Poly::random(&params.ctx, Representation::Ntt);
				let mut c0 = Poly::small(&params.ctx, Representation::PowerBasis, 16)?;
				c0.change_representation(Representation::Ntt);
				c0 -= &(&c1 * &sk.s[0]);
				c0 -= &(&c2 * &sk.s[1]);

				// Relinearize the extended ciphertext!
				c2.change_representation(Representation::PowerBasis);
				rk.relinearize(&mut c0, &mut c1, &c2)?;

				let ct = Ciphertext {
					par: params.clone(),
					seed: None,
					c: vec![c0, c1],
				};

				// Print the noise and decrypt
				println!("Noise: {}", unsafe { sk.measure_noise(&ct)? });
				let pt = sk.decrypt(&ct)?;
				assert!(Vec::<u64>::try_decode(&pt, Encoding::Poly).is_ok_and(|v| v == &[0u64; 8]))
			}
		}
		Ok(())
	}

	#[test]
	fn test_proto_conversion() -> Result<(), String> {
		for params in [
			Arc::new(BfvParameters::default(1)),
			Arc::new(BfvParameters::default(2)),
		] {
			let sk = SecretKey::random(&params);
			let rk = RelinearizationKey::new(&sk)?;
			let proto = RelinearizationKeyProto::from(&rk);
			assert_eq!(rk, RelinearizationKey::try_convert_from(&proto, &params)?);
		}
		Ok(())
	}
}
