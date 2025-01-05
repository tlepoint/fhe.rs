use std::ops::Mul;

use crate::proto::bfv::{
    KeySwitchingKey as KeySwitchingKeyProto, RgswCiphertext as RGSWCiphertextProto,
};
use crate::{Error, Result};
use fhe_math::rq::{traits::TryConvertFrom as TryConvertFromPoly, Poly, Representation};
use fhe_traits::{
    DeserializeParametrized, FheCiphertext, FheEncrypter, FheParametrized, Serialize,
};
use prost::Message;
use rand::{CryptoRng, RngCore};
use zeroize::Zeroizing;

use super::{
    keys::KeySwitchingKey, traits::TryConvertFrom, BfvParameters, Ciphertext, Plaintext, SecretKey,
};

/// A RGSW ciphertext encrypting a plaintext.
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
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
        RGSWCiphertextProto {
            ksk0: Some(KeySwitchingKeyProto::from(&ct.ksk0)),
            ksk1: Some(KeySwitchingKeyProto::from(&ct.ksk1)),
        }
    }
}

impl TryConvertFrom<&RGSWCiphertextProto> for RGSWCiphertext {
    fn try_convert_from(
        value: &RGSWCiphertextProto,
        par: &std::sync::Arc<BfvParameters>,
    ) -> Result<Self> {
        let ksk0 = KeySwitchingKey::try_convert_from(
            value.ksk0.as_ref().ok_or(Error::SerializationError)?,
            par,
        )?;
        let ksk1 = KeySwitchingKey::try_convert_from(
            value.ksk1.as_ref().ok_or(Error::SerializationError)?,
            par,
        )?;
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
        let proto = Message::decode(bytes).map_err(|_| Error::SerializationError)?;
        RGSWCiphertext::try_convert_from(&proto, par)
    }
}

impl Serialize for RGSWCiphertext {
    fn to_bytes(&self) -> Vec<u8> {
        RGSWCiphertextProto::from(self).encode_to_vec()
    }
}

impl FheCiphertext for RGSWCiphertext {}

impl FheEncrypter<Plaintext, RGSWCiphertext> for SecretKey {
    type Error = Error;

    fn try_encrypt<R: RngCore + CryptoRng>(
        &self,
        pt: &Plaintext,
        rng: &mut R,
    ) -> Result<RGSWCiphertext> {
        let level = pt.level;
        let ctx = self.par.ctx_at_level(level)?;

        let mut m = Zeroizing::new(pt.poly_ntt.clone());
        let mut m_s = Zeroizing::new(Poly::try_convert_from(
            self.coeffs.as_ref(),
            ctx,
            false,
            Representation::PowerBasis,
        )?);
        m_s.change_representation(Representation::Ntt);
        *m_s.as_mut() *= m.as_ref();
        m_s.change_representation(Representation::PowerBasis);
        m.change_representation(Representation::PowerBasis);

        let ksk0 = KeySwitchingKey::new(self, &m, pt.level, pt.level, rng)?;
        let ksk1 = KeySwitchingKey::new(self, &m_s, pt.level, pt.level, rng)?;

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
        assert_eq!(self.len(), 2, "Ciphertext must have two parts");

        let mut ct0 = self[0].clone();
        let mut ct1 = self[1].clone();
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
    use std::error::Error;

    use crate::bfv::{BfvParameters, Ciphertext, Encoding, Plaintext, SecretKey};
    use fhe_traits::{DeserializeParametrized, FheDecrypter, FheEncoder, FheEncrypter, Serialize};
    use rand::thread_rng;

    use super::RGSWCiphertext;

    #[test]
    fn external_product() -> Result<(), Box<dyn Error>> {
        let mut rng = thread_rng();
        for params in [
            BfvParameters::default_arc(2, 16),
            BfvParameters::default_arc(8, 16),
        ] {
            let sk = SecretKey::random(&params, &mut rng);
            let v1 = params.plaintext.random_vec(params.degree(), &mut rng);
            let v2 = params.plaintext.random_vec(params.degree(), &mut rng);

            let pt1 = Plaintext::try_encode(&v1, Encoding::simd(), &params)?;
            let pt2 = Plaintext::try_encode(&v2, Encoding::simd(), &params)?;

            let ct1: Ciphertext = sk.try_encrypt(&pt1, &mut rng)?;
            let ct2: Ciphertext = sk.try_encrypt(&pt2, &mut rng)?;
            let ct2_rgsw: RGSWCiphertext = sk.try_encrypt(&pt2, &mut rng)?;

            let product = &ct1 * &ct2;
            let expected = sk.try_decrypt(&product)?;

            let ct3 = &ct1 * &ct2_rgsw;
            let ct4 = &ct2_rgsw * &ct1;

            println!("Noise 1: {:?}", unsafe { sk.measure_noise(&ct3) });
            println!("Noise 2: {:?}", unsafe { sk.measure_noise(&ct4) });
            assert_eq!(expected, sk.try_decrypt(&ct3)?);
            assert_eq!(expected, sk.try_decrypt(&ct4)?);
        }
        Ok(())
    }

    #[test]
    fn serialize() -> Result<(), Box<dyn Error>> {
        let mut rng = thread_rng();
        for params in [
            BfvParameters::default_arc(6, 16),
            BfvParameters::default_arc(5, 16),
        ] {
            let sk = SecretKey::random(&params, &mut rng);
            let v = params.plaintext.random_vec(params.degree(), &mut rng);
            let pt = Plaintext::try_encode(&v, Encoding::simd(), &params)?;
            let ct: RGSWCiphertext = sk.try_encrypt(&pt, &mut rng)?;

            let bytes = ct.to_bytes();
            assert_eq!(RGSWCiphertext::from_bytes(&bytes, &params)?, ct);
        }

        Ok(())
    }
}
