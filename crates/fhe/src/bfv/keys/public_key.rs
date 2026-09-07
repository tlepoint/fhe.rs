//! Public keys for the BFV encryption scheme

use crate::bfv::wire::FromProto;
use crate::bfv::{Ciphertext, Parameters, Plaintext};
use crate::proto::bfv::{Ciphertext as CiphertextProto, PublicKey as PublicKeyProto};
use crate::{Error, Result, error::SerializationError};
use fhe_math::rq::{Ntt, Poly};

use prost::Message;
use rand::{CryptoRng, Rng as RngCore};
use std::borrow::Cow;
use zeroize::Zeroizing;

use super::SecretKey;

/// Public key for the BFV encryption scheme.
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct PublicKey {
    pub(crate) par: Parameters,
    pub(crate) c: Ciphertext,
}

impl PublicKey {
    /// Generate a new [`PublicKey`] from a [`SecretKey`].
    pub fn from_secret_key<R: RngCore + CryptoRng>(sk: &SecretKey, rng: &mut R) -> Self {
        let zero = Plaintext::zero(&sk.par, 0).unwrap();
        let mut c: Ciphertext = sk.encrypt(&zero, rng).unwrap();
        // The polynomials of a public key should not allow for variable time
        // computation. Only timing metadata changes, so the seed remains valid.
        c.c.iter_mut()
            .for_each(|p| p.disallow_variable_time_computations());
        Self {
            par: sk.par.clone(),
            c,
        }
    }
}

impl PublicKey {
    /// Encrypt a plaintext with compatible parameters using caller-owned
    /// cryptographic randomness.
    pub fn encrypt<R: RngCore + CryptoRng>(
        &self,
        pt: &Plaintext,
        rng: &mut R,
    ) -> Result<Ciphertext> {
        pt.validate_for(&self.par)?;
        self.c.validate_for(&self.par)?;
        let plaintext_level = pt.level();
        if plaintext_level < self.c.level {
            return Err(Error::InvalidLevel {
                level: plaintext_level,
                min_level: self.c.level,
                max_level: self.par.max_level(),
            });
        }

        let needs_switch = self.c.level != plaintext_level;
        let ct: Cow<'_, Ciphertext> = if needs_switch {
            let mut owned = self.c.clone();
            while owned.level != plaintext_level {
                owned.switch_down()?;
            }
            Cow::Owned(owned)
        } else {
            Cow::Borrowed(&self.c)
        };

        let ctx = self.par.context_at_level(ct.level)?;
        let u = Zeroizing::new(Poly::<Ntt>::small(ctx, self.par.inner.variance, rng)?);
        let e1 = Zeroizing::new(Poly::<Ntt>::small(ctx, self.par.inner.variance, rng)?);
        let e2 = Zeroizing::new(Poly::<Ntt>::small(ctx, self.par.inner.variance, rng)?);

        let m = Zeroizing::new(pt.to_poly());
        let mut c0 = u.as_ref() * &ct.c[0];
        c0 += &e1;
        c0 += &m;
        let mut c1 = u.as_ref() * &ct.c[1];
        c1 += &e2;

        // It is now safe to enable variable time computations.
        let variable_time = crate::VariableTime::new(crate::PublicData::assert_public());
        c0.allow_variable_time_computations(variable_time);
        c1.allow_variable_time_computations(variable_time);

        Ok(Ciphertext {
            par: self.par.clone(),
            seed: None,
            c: vec![c0, c1],
            level: ct.level,
        })
    }
}

impl From<&PublicKey> for PublicKeyProto {
    fn from(pk: &PublicKey) -> Self {
        PublicKeyProto {
            c: Some(CiphertextProto::from(&pk.c)),
        }
    }
}

impl PublicKey {
    /// Serialize in the existing protobuf wire format.
    #[must_use]
    pub fn to_bytes(&self) -> Vec<u8> {
        PublicKeyProto::from(self).encode_to_vec()
    }
}

impl PublicKey {
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
            crate::error::SerializedObject::PublicKey,
            Some(par),
            limits,
        )?;
        let proto: PublicKeyProto = Message::decode(bytes).map_err(|source| {
            Error::SerializationError(SerializationError::Decode {
                object: crate::error::SerializedObject::PublicKey,
                source,
            })
        })?;
        if let Some(proto_c) = &proto.c {
            let mut c = Ciphertext::from_proto(proto_c, par, limits)?;
            if c.level != 0 {
                Err(Error::SerializationError(
                    SerializationError::InvalidPublicKeyLevel {
                        actual: c.level,
                        expected: 0,
                    },
                ))
            } else {
                // The polynomials of a public key should not allow for variable time
                // computation. Only timing metadata changes, so the seed remains valid.
                c.c.iter_mut()
                    .for_each(|p| p.disallow_variable_time_computations());
                Ok(Self {
                    par: par.clone(),
                    c,
                })
            }
        } else {
            Err(Error::SerializationError(
                SerializationError::MissingField {
                    field: crate::error::SerializedField::PublicKeyCiphertext,
                },
            ))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::PublicKey;
    use crate::bfv::{Encoding, Plaintext, SecretKey, parameters::Parameters};

    use rand::rng;
    use std::error::Error;

    #[test]
    fn keygen() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let params = Parameters::test_parameters(1, 16);
        let sk = SecretKey::generate(&params, &mut rng);
        let pk = PublicKey::from_secret_key(&sk, &mut rng);
        assert_eq!(pk.par, params);
        assert_eq!(
            sk.decrypt(&pk.c)?.poly_ntt,
            Plaintext::zero(&params, 0)?.poly_ntt
        );
        Ok(())
    }

    #[test]
    fn encrypt_decrypt() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for params in [
            Parameters::test_parameters(1, 16),
            Parameters::test_parameters(6, 16),
        ] {
            for level in 0..params.max_level() {
                for _ in 0..20 {
                    let sk = SecretKey::generate(&params, &mut rng);
                    let pk = PublicKey::from_secret_key(&sk, &mut rng);

                    let pt = Plaintext::encode_at_level(
                        &params,
                        &fhe_math::zq::Modulus::new(params.plaintext_modulus_u64().unwrap())
                            .unwrap()
                            .random_vec(params.degree(), &mut rng),
                        Encoding::Polynomial,
                        level,
                    )?;
                    let ct = pk.encrypt(&pt, &mut rng)?;
                    let pt2 = sk.decrypt(&ct)?;

                    println!(
                        "Noise: {}",
                        sk.measure_noise_vartime(
                            &ct,
                            crate::SecretDependentDiagnostics::acknowledge_leakage()
                        )?
                    );
                    assert_eq!(pt2.poly_ntt, pt.poly_ntt);
                }
            }
        }

        Ok(())
    }

    #[test]
    fn encrypt_rejects_mismatched_parameters() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let params = Parameters::test_parameters(1, 16);
        let other_params = Parameters::test_parameters(1, 32);
        let sk = SecretKey::generate(&params, &mut rng);
        let pk = PublicKey::from_secret_key(&sk, &mut rng);
        let pt = Plaintext::encode(&other_params, &[1u64][..], Encoding::Polynomial)?;

        assert!(pk.encrypt(&pt, &mut rng).is_err());
        Ok(())
    }

    #[test]
    fn test_serialize() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for params in [
            Parameters::test_parameters(1, 16),
            Parameters::test_parameters(6, 16),
        ] {
            let sk = SecretKey::generate(&params, &mut rng);
            let pk = PublicKey::from_secret_key(&sk, &mut rng);
            let bytes = pk.to_bytes();
            assert_eq!(pk, PublicKey::from_bytes(&bytes, &params)?);
        }
        Ok(())
    }
    #[test]
    fn timing_policy_preserves_seeded_serialization() -> Result<(), Box<dyn Error>> {
        let params = Parameters::test_parameters(2, 16);
        let mut rng = rng();
        let sk = SecretKey::generate(&params, &mut rng);
        let pk = PublicKey::from_secret_key(&sk, &mut rng);
        assert!(pk.c.seed.is_some());
        let bytes = pk.to_bytes();
        let restored = PublicKey::from_bytes(&bytes, &params)?;
        assert_eq!(restored, pk);
        assert!(restored.c.seed.is_some());
        assert!(
            restored
                .c
                .iter()
                .all(|poly| !poly.allows_variable_time_computations())
        );
        let mut unseeded = pk;
        unseeded.c.seed = None;
        assert!(bytes.len() < unseeded.to_bytes().len());
        Ok(())
    }
}
