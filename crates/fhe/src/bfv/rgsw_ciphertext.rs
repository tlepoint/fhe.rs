use crate::proto::bfv::{
    KeySwitchingKey as KeySwitchingKeyProto, RgswCiphertext as RgswCiphertextProto,
};
use crate::{Error, Result, error::SerializationError};
use fhe_math::rq::{Ntt, Poly, PowerBasis};

use prost::Message;
use rand::{CryptoRng, Rng as RngCore};
use zeroize::Zeroizing;

use super::{Ciphertext, Parameters, Plaintext, SecretKey, keys::KeySwitchingKey, wire::FromProto};

/// A RGSW ciphertext encrypting a plaintext.
#[derive(Debug, PartialEq, Eq)]
pub struct RgswCiphertext {
    ksk0: KeySwitchingKey,
    ksk1: KeySwitchingKey,
}

impl From<&RgswCiphertext> for RgswCiphertextProto {
    fn from(ct: &RgswCiphertext) -> Self {
        RgswCiphertextProto {
            ksk0: Some(KeySwitchingKeyProto::from(&ct.ksk0)),
            ksk1: Some(KeySwitchingKeyProto::from(&ct.ksk1)),
        }
    }
}

impl FromProto<&RgswCiphertextProto> for RgswCiphertext {
    fn from_proto(
        value: &RgswCiphertextProto,
        par: &Parameters,
        limits: &crate::DecodeLimits,
    ) -> Result<Self> {
        let ksk0 = KeySwitchingKey::from_proto(
            value.ksk0.as_ref().ok_or(Error::SerializationError(
                SerializationError::MissingField {
                    field: crate::error::SerializedField::RgswKeySwitchingKey0,
                },
            ))?,
            par,
            limits,
        )?;
        let ksk1 = KeySwitchingKey::from_proto(
            value.ksk1.as_ref().ok_or(Error::SerializationError(
                SerializationError::MissingField {
                    field: crate::error::SerializedField::RgswKeySwitchingKey1,
                },
            ))?,
            par,
            limits,
        )?;
        if ksk0.ksk_level != ksk0.ciphertext_level
            || ksk0.ciphertext_level != ksk1.ciphertext_level
            || ksk1.ciphertext_level != ksk1.ksk_level
        {
            return Err(Error::SerializationError(
                SerializationError::InconsistentKeySwitchingLevels,
            ));
        }

        Ok(Self { ksk0, ksk1 })
    }
}

impl RgswCiphertext {
    /// Import validated protobuf bytes, binding contextual values to the
    /// supplied parameters.
    pub fn from_bytes(bytes: &[u8], par: &Parameters) -> Result<Self> {
        Self::from_bytes_with_limits(bytes, par, &crate::DecodeLimits::default())
    }

    /// Import with explicit resource bounds checked before allocation.
    pub fn from_bytes_with_limits(
        bytes: &[u8],
        par: &Parameters,
        limits: &crate::DecodeLimits,
    ) -> Result<Self> {
        crate::bfv::wire::preflight(
            bytes,
            crate::error::SerializedObject::RgswCiphertext,
            Some(par),
            limits,
        )?;
        let proto = Message::decode(bytes).map_err(|source| {
            Error::SerializationError(SerializationError::Decode {
                object: crate::error::SerializedObject::RgswCiphertext,
                source,
            })
        })?;
        RgswCiphertext::from_proto(&proto, par, limits)
    }
}

impl RgswCiphertext {
    /// Serialize in the existing protobuf wire format.
    #[must_use]
    pub fn to_bytes(&self) -> Vec<u8> {
        RgswCiphertextProto::from(self).encode_to_vec()
    }
}

impl SecretKey {
    /// Encrypt a plaintext as an RGSW ciphertext using caller-owned
    /// cryptographic randomness.
    pub fn encrypt_rgsw<R: RngCore + CryptoRng>(
        &self,
        pt: &Plaintext,
        rng: &mut R,
    ) -> Result<RgswCiphertext> {
        pt.validate_for(&self.par)?;
        let level = pt.level();
        let ctx = self.par.context_at_level(level)?;

        let m = Zeroizing::new(pt.poly_ntt.clone().into_power_basis());
        let mut m_s = Zeroizing::new(
            Poly::<PowerBasis>::from_signed_coefficients(self.coeffs.as_ref(), ctx)?.into_ntt(),
        );
        *m_s.as_mut() *= pt.poly_ntt.as_ref();
        let ctx = m_s.ctx().clone();
        let m_s_inner = std::mem::replace(m_s.as_mut(), Poly::<Ntt>::zero(&ctx));
        let m_s = Zeroizing::new(m_s_inner.into_power_basis());

        let ksk0 = KeySwitchingKey::new(self, &m, level, level, rng)?;
        let ksk1 = KeySwitchingKey::new(self, &m_s, level, level, rng)?;

        Ok(RgswCiphertext { ksk0, ksk1 })
    }
}

impl Ciphertext {
    /// Replace this ciphertext with its RGSW external product only on success.
    pub fn multiply_rgsw_assign(&mut self, rhs: &RgswCiphertext) -> Result<()> {
        *self = self.multiply_rgsw(rhs)?;
        Ok(())
    }

    /// Compute the RGSW external product at the same level. The BFV input must
    /// contain exactly two components and use compatible parameters.
    pub fn multiply_rgsw(&self, rhs: &RgswCiphertext) -> Result<Self> {
        self.validate_for_context(
            &rhs.ksk0.par,
            rhs.ksk0.ciphertext_level,
            &rhs.ksk0.ctx_ciphertext,
        )?;
        if self.len() != 2 {
            return Err(crate::error::CiphertextError::InvalidPolynomialCount {
                operation: crate::error::CiphertextOperation::RgswProduct,
                actual: self.len(),
                expected: 2,
            }
            .into());
        }
        let ct0 = self.c[0].clone().into_power_basis();
        let ct1 = self.c[1].clone().into_power_basis();

        let mut c0 = Poly::<Ntt>::zero(&rhs.ksk0.ctx_ksk);
        let mut c1 = Poly::<Ntt>::zero(&rhs.ksk0.ctx_ksk);
        rhs.ksk0.key_switch_assign(&ct0, &mut c0, &mut c1)?;

        let mut c0p = Poly::<Ntt>::zero(&rhs.ksk1.ctx_ksk);
        let mut c1p = Poly::<Ntt>::zero(&rhs.ksk1.ctx_ksk);
        rhs.ksk1.key_switch_assign(&ct1, &mut c0p, &mut c1p)?;

        c0 += &c0p;
        c1 += &c1p;
        Ciphertext::from_components(vec![c0, c1], &self.par)
    }
}

#[cfg(test)]
mod tests {
    use std::error::Error;

    use crate::bfv::{Ciphertext, Encoding, Parameters, Plaintext, SecretKey};

    use rand::rng;

    use super::RgswCiphertext;

    #[test]
    fn import_rejects_different_switching_and_ciphertext_levels() -> crate::Result<()> {
        use crate::bfv::keys::KeySwitchingKey;
        use fhe_math::rq::{Poly, PowerBasis};
        let par = Parameters::test_parameters(3, 16);
        let mut rng = rng();
        let sk = SecretKey::generate(&par, &mut rng);
        let pt = Plaintext::encode(&par, &[2, 3], Encoding::Simd)?;
        let secret = Poly::<PowerBasis>::from_signed_coefficients(
            sk.coeffs.as_ref(),
            par.context_at_level(0)?,
        )?
        .into_ntt();
        let m_s = (&pt.poly_ntt * &secret).into_power_basis();
        let m = pt.poly_ntt.clone().into_power_basis();
        let rgsw = RgswCiphertext {
            ksk0: KeySwitchingKey::new(&sk, &m, 1, 0, &mut rng)?,
            ksk1: KeySwitchingKey::new(&sk, &m_s, 1, 0, &mut rng)?,
        };
        assert!(matches!(
            RgswCiphertext::from_bytes(&rgsw.to_bytes(), &par),
            Err(crate::Error::SerializationError(
                crate::error::SerializationError::InconsistentKeySwitchingLevels
            ))
        ));
        Ok(())
    }

    #[test]
    fn external_product() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for params in [
            Parameters::test_parameters(2, 16),
            Parameters::test_parameters(8, 16),
        ] {
            let sk = SecretKey::generate(&params, &mut rng);
            let v1 = fhe_math::zq::Modulus::new(params.plaintext_modulus_u64().unwrap())
                .unwrap()
                .random_vec(params.degree(), &mut rng);
            let v2 = fhe_math::zq::Modulus::new(params.plaintext_modulus_u64().unwrap())
                .unwrap()
                .random_vec(params.degree(), &mut rng);

            let pt1 = Plaintext::encode(&params, &v1, Encoding::Simd)?;
            let pt2 = Plaintext::encode(&params, &v2, Encoding::Simd)?;

            let ct1: Ciphertext = sk.encrypt(&pt1, &mut rng)?;
            let ct2: Ciphertext = sk.encrypt(&pt2, &mut rng)?;
            let ct2_rgsw: RgswCiphertext = sk.encrypt_rgsw(&pt2, &mut rng)?;

            let product = ct1.multiply(&ct2).unwrap();
            let expected = sk.decrypt(&product)?;

            let ct3 = ct1.multiply_rgsw(&ct2_rgsw).unwrap();
            let ct4 = ct1.multiply_rgsw(&ct2_rgsw).unwrap();

            println!(
                "Noise 1: {:?}",
                sk.measure_noise_vartime(
                    &ct3,
                    crate::SecretDependentDiagnostics::acknowledge_leakage()
                )
            );
            println!(
                "Noise 2: {:?}",
                sk.measure_noise_vartime(
                    &ct4,
                    crate::SecretDependentDiagnostics::acknowledge_leakage()
                )
            );
            assert_eq!(expected, sk.decrypt(&ct3)?);
            assert_eq!(expected, sk.decrypt(&ct4)?);
        }
        Ok(())
    }

    #[test]
    fn encryption_rejects_mismatched_parameters() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let params = Parameters::test_parameters(1, 16);
        let other_params = Parameters::test_parameters(1, 32);
        let sk = SecretKey::generate(&params, &mut rng);
        let pt = Plaintext::encode(&other_params, &[1u64][..], Encoding::Polynomial)?;
        let encrypted: crate::Result<RgswCiphertext> = sk.encrypt_rgsw(&pt, &mut rng);

        assert!(encrypted.is_err());
        Ok(())
    }

    #[test]
    fn serialize() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for params in [
            Parameters::test_parameters(6, 16),
            Parameters::test_parameters(5, 16),
        ] {
            let sk = SecretKey::generate(&params, &mut rng);
            let v = fhe_math::zq::Modulus::new(params.plaintext_modulus_u64().unwrap())
                .unwrap()
                .random_vec(params.degree(), &mut rng);
            let pt = Plaintext::encode(&params, &v, Encoding::Simd)?;
            let ct: RgswCiphertext = sk.encrypt_rgsw(&pt, &mut rng)?;

            let bytes = ct.to_bytes();
            assert_eq!(RgswCiphertext::from_bytes(&bytes, &params)?, ct);
        }

        Ok(())
    }
}
