use std::ops::Mul;

use fhers_traits::{
	DeserializeParametrized, FheCiphertext, FheEncrypter, FheParametrized, Serialize,
};
use math::rq::{traits::TryConvertFrom as TryConvertFromPoly, Poly, Representation};
use protobuf::{Message, MessageField};

use crate::{
	bfv::proto::bfv::{
		KeySwitchingKey as KeySwitchingKeyProto, RGSWCiphertext as RGSWCiphertextProto,
	},
	Error, Result,
};

use super::{
	keys::KeySwitchingKey, traits::TryConvertFrom, BfvParameters, Ciphertext, Plaintext, SecretKey,
};

#[doc(cfg(feature = "rgsw"))]
/// A RGSW ciphertext encrypting a plaintext.
#[derive(Debug, PartialEq, Eq)]
pub struct RGSWCiphertext {
	ksk0: KeySwitchingKey,
	ksk1: KeySwitchingKey,
}

impl FheParametrized for RGSWCiphertext {
	type Parameters = BfvParameters;
}

impl From<&RGSWCiphertext> for RGSWCiphertextProto {
	fn from(ct: &RGSWCiphertext) -> Self {
		let mut proto = RGSWCiphertextProto::new();
		proto.ksk0 = MessageField::some(KeySwitchingKeyProto::from(&ct.ksk0));
		proto.ksk1 = MessageField::some(KeySwitchingKeyProto::from(&ct.ksk1));
		proto
	}
}

impl TryConvertFrom<&RGSWCiphertextProto> for RGSWCiphertext {
	fn try_convert_from(
		value: &RGSWCiphertextProto,
		par: &std::sync::Arc<BfvParameters>,
	) -> Result<Self> {
		if value.ksk0.is_none() || value.ksk1.is_none() {
			return Err(Error::SerializationError);
		}

		let ksk0 = KeySwitchingKey::try_convert_from(value.ksk0.as_ref().unwrap(), par)?;
		let ksk1 = KeySwitchingKey::try_convert_from(value.ksk1.as_ref().unwrap(), par)?;
		if ksk0.ksk_level != ksk0.ciphertext_level
			|| ksk0.ciphertext_level != ksk1.ciphertext_level
			|| ksk1.ciphertext_level != ksk1.ksk_level
		{
			return Err(Error::SerializationError);
		}

		Ok(Self { ksk0, ksk1 })
	}
}

impl DeserializeParametrized for RGSWCiphertext {
	type Error = Error;

	fn from_bytes(bytes: &[u8], par: &std::sync::Arc<Self::Parameters>) -> Result<Self> {
		let proto =
			RGSWCiphertextProto::parse_from_bytes(bytes).map_err(|_| Error::SerializationError)?;
		RGSWCiphertext::try_convert_from(&proto, par)
	}
}

impl Serialize for RGSWCiphertext {
	fn to_bytes(&self) -> Vec<u8> {
		RGSWCiphertextProto::from(self).write_to_bytes().unwrap()
	}
}

impl FheCiphertext for RGSWCiphertext {}

impl FheEncrypter<Plaintext, RGSWCiphertext> for SecretKey {
	type Error = Error;

	fn try_encrypt(&self, pt: &Plaintext) -> Result<RGSWCiphertext> {
		let level = pt.level;
		let ctx = self.par.ctx_at_level(level)?;

		let mut m = pt.poly_ntt.clone();
		let mut m_s = Poly::try_convert_from(
			&self.s_coefficients as &[i64],
			ctx,
			false,
			Representation::PowerBasis,
		)?;
		m_s.change_representation(Representation::Ntt);
		m_s *= &m;
		m_s.change_representation(Representation::PowerBasis);
		m.change_representation(Representation::PowerBasis);

		let ksk0 = KeySwitchingKey::new(self, &m, pt.level, pt.level)?;
		let ksk1 = KeySwitchingKey::new(self, &m_s, pt.level, pt.level)?;

		Ok(RGSWCiphertext { ksk0, ksk1 })
	}
}

impl Mul<&RGSWCiphertext> for &Ciphertext {
	type Output = Ciphertext;

	fn mul(self, rhs: &RGSWCiphertext) -> Self::Output {
		assert_eq!(
			self.par, rhs.ksk0.par,
			"Ciphertext and RGSWCiphertext must have the same parameters"
		);
		assert_eq!(
			self.level, rhs.ksk0.ciphertext_level,
			"Ciphertext and RGSWCiphertext must have the same level"
		);
		assert_eq!(self.c.len(), 2, "Ciphertext must have two parts");

		let mut ct0 = self.c[0].clone();
		let mut ct1 = self.c[1].clone();
		ct0.change_representation(Representation::PowerBasis);
		ct1.change_representation(Representation::PowerBasis);

		let (c0, c1) = rhs.ksk0.key_switch(&ct0).unwrap();
		let (c0p, c1p) = rhs.ksk1.key_switch(&ct1).unwrap();

		Ciphertext {
			par: self.par.clone(),
			seed: None,
			c: vec![&c0 + &c0p, &c1 + &c1p],
			level: self.level,
		}
	}
}

impl Mul<&Ciphertext> for &RGSWCiphertext {
	type Output = Ciphertext;

	fn mul(self, rhs: &Ciphertext) -> Self::Output {
		rhs * self
	}
}

#[cfg(test)]
mod tests {
	use std::{error::Error, sync::Arc};

	use fhers_traits::{
		DeserializeParametrized, FheDecoder, FheDecrypter, FheEncoder, FheEncrypter, Serialize,
	};

	use crate::bfv::{BfvParameters, Ciphertext, Encoding, Plaintext, SecretKey};

	use super::RGSWCiphertext;

	#[test]
	fn external_product() -> Result<(), Box<dyn Error>> {
		for params in [
			Arc::new(BfvParameters::default(6, 8)),
			Arc::new(BfvParameters::default(8, 8)),
		] {
			let sk = SecretKey::random(&params);
			let v1 = params.plaintext.random_vec(params.degree());
			let v2 = params.plaintext.random_vec(params.degree());
			let mut expected = v1.clone();
			params.plaintext.mul_vec(&mut expected, &v2);

			let pt1 = Plaintext::try_encode(&v1 as &[u64], Encoding::simd(), &params)?;
			let pt2 = Plaintext::try_encode(&v2 as &[u64], Encoding::simd(), &params)?;

			let ct: Ciphertext = sk.try_encrypt(&pt1)?;
			let ct2: RGSWCiphertext = sk.try_encrypt(&pt2)?;

			let ct3 = &ct * &ct2;
			let ct4 = &ct2 * &ct;
			let pt3 = sk.try_decrypt(&ct3)?;
			let pt4 = sk.try_decrypt(&ct4)?;

			println!("Noise: {:?}", unsafe { sk.measure_noise(&ct3) });
			println!("Noise: {:?}", unsafe { sk.measure_noise(&ct4) });
			assert_eq!(expected, Vec::<u64>::try_decode(&pt3, Encoding::simd())?);
			assert_eq!(expected, Vec::<u64>::try_decode(&pt4, Encoding::simd())?);
		}
		Ok(())
	}

	#[test]
	fn serialize() -> Result<(), Box<dyn Error>> {
		for params in [
			Arc::new(BfvParameters::default(6, 8)),
			Arc::new(BfvParameters::default(5, 8)),
		] {
			let sk = SecretKey::random(&params);
			let v = params.plaintext.random_vec(params.degree());
			let pt = Plaintext::try_encode(&v as &[u64], Encoding::simd(), &params)?;
			let ct: RGSWCiphertext = sk.try_encrypt(&pt)?;

			let bytes = ct.to_bytes();
			assert_eq!(RGSWCiphertext::from_bytes(&bytes, &params)?, ct);
		}

		Ok(())
	}
}
