//! Public keys for the BFV encryption scheme

use crate::bfv::traits::TryConvertFrom;
use crate::bfv::{BfvParameters, Ciphertext, Encoding, Plaintext};
use crate::proto::bfv::{Ciphertext as CiphertextProto, PublicKey as PublicKeyProto};
use crate::{Error, Result, SerializationError};
use fhe_math::rq::{Ntt, Poly};
use fhe_traits::{DeserializeParametrized, FheEncrypter, FheParametrized, Serialize};
use prost::Message;
use rand::{CryptoRng, Rng as RngCore};
use std::borrow::Cow;
use std::sync::Arc;
use zeroize::Zeroizing;

use super::SecretKey;

/// Public key for the BFV encryption scheme.
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct PublicKey {
    pub(crate) par: Arc<BfvParameters>,
    pub(crate) c: Ciphertext,
}

impl PublicKey {
    /// Generate a new [`PublicKey`] from a [`SecretKey`].
    pub fn new<R: RngCore + CryptoRng>(sk: &SecretKey, rng: &mut R) -> Self {
        let zero = Plaintext::zero(Encoding::poly(), &sk.par).unwrap();
        let mut c: Ciphertext = sk.try_encrypt(&zero, rng).unwrap();
        // The polynomials of a public key should not allow for variable time
        // computation.
        c.iter_mut()
            .for_each(|p| p.disallow_variable_time_computations());
        Self {
            par: sk.par.clone(),
            c,
        }
    }
}

impl FheParametrized for PublicKey {
    type Parameters = BfvParameters;
}

impl FheEncrypter<Plaintext, Ciphertext> for PublicKey {
    type Error = Error;

    fn try_encrypt<R: RngCore + CryptoRng>(
        &self,
        pt: &Plaintext,
        rng: &mut R,
    ) -> Result<Ciphertext> {
        pt.validate_for(&self.par)?;
        self.c.validate_for(&self.par)?;
        if pt.level < self.c.level {
            return Err(Error::InvalidLevel {
                level: pt.level,
                min_level: self.c.level,
                max_level: self.par.max_level(),
            });
        }

        let needs_switch = self.c.level != pt.level;
        let ct: Cow<'_, Ciphertext> = if needs_switch {
            let mut owned = self.c.clone();
            while owned.level != pt.level {
                owned.switch_down()?;
            }
            Cow::Owned(owned)
        } else {
            Cow::Borrowed(&self.c)
        };

        let ctx = self.par.context_at_level(ct.level)?;
        let u = Zeroizing::new(Poly::<Ntt>::small(ctx, self.par.variance, rng)?);
        let e1 = Zeroizing::new(Poly::<Ntt>::small(ctx, self.par.variance, rng)?);
        let e2 = Zeroizing::new(Poly::<Ntt>::small(ctx, self.par.variance, rng)?);

        let m = Zeroizing::new(pt.to_poly());
        let mut c0 = u.as_ref() * &ct[0];
        c0 += &e1;
        c0 += &m;
        let mut c1 = u.as_ref() * &ct[1];
        c1 += &e2;

        // It is now safe to enable variable time computations.
        let variable_time = fhe_traits::VariableTime::new(fhe_traits::PublicData::assert_public());
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

impl Serialize for PublicKey {
    fn to_bytes(&self) -> Vec<u8> {
        PublicKeyProto::from(self).encode_to_vec()
    }
}

impl DeserializeParametrized for PublicKey {
    type Error = Error;

    fn from_bytes(bytes: &[u8], par: &Arc<Self::Parameters>) -> Result<Self> {
        let proto: PublicKeyProto = Message::decode(bytes).map_err(|_| {
            Error::SerializationError(SerializationError::Decode {
                object: crate::SerializedObject::PublicKey,
            })
        })?;
        if let Some(proto_c) = &proto.c {
            let mut c = Ciphertext::try_convert_from(proto_c, par)?;
            if c.level != 0 {
                Err(Error::SerializationError(
                    SerializationError::InvalidPublicKeyLevel {
                        actual: c.level,
                        expected: 0,
                    },
                ))
            } else {
                // The polynomials of a public key should not allow for variable time
                // computation.
                c.iter_mut()
                    .for_each(|p| p.disallow_variable_time_computations());
                Ok(Self {
                    par: par.clone(),
                    c,
                })
            }
        } else {
            Err(Error::SerializationError(
                SerializationError::MissingField {
                    field: crate::SerializedField::PublicKeyCiphertext,
                },
            ))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::PublicKey;
    use crate::bfv::{Encoding, Plaintext, SecretKey, parameters::BfvParameters};
    use fhe_traits::{DeserializeParametrized, FheDecrypter, FheEncoder, FheEncrypter, Serialize};
    use rand::rng;
    use std::error::Error;

    #[test]
    fn keygen() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let params = BfvParameters::default_arc(1, 16);
        let sk = SecretKey::random(&params, &mut rng);
        let pk = PublicKey::new(&sk, &mut rng);
        assert_eq!(pk.par, params);
        assert_eq!(
            sk.try_decrypt(&pk.c)?,
            Plaintext::zero(Encoding::poly(), &params)?
        );
        Ok(())
    }

    #[test]
    fn encrypt_decrypt() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for params in [
            BfvParameters::default_arc(1, 16),
            BfvParameters::default_arc(6, 16),
        ] {
            for level in 0..params.max_level() {
                for _ in 0..20 {
                    let sk = SecretKey::random(&params, &mut rng);
                    let pk = PublicKey::new(&sk, &mut rng);

                    let pt = Plaintext::try_encode(
                        &fhe_math::zq::Modulus::new(params.plaintext())
                            .unwrap()
                            .random_vec(params.degree(), &mut rng),
                        Encoding::poly_at_level(level),
                        &params,
                    )?;
                    let ct = pk.try_encrypt(&pt, &mut rng)?;
                    let pt2 = sk.try_decrypt(&ct)?;

                    println!("Noise: {}", unsafe { sk.measure_noise(&ct)? });
                    assert_eq!(pt2, pt);
                }
            }
        }

        Ok(())
    }

    #[test]
    fn encrypt_rejects_mismatched_parameters() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let params = BfvParameters::default_arc(1, 16);
        let other_params = BfvParameters::default_arc(1, 16);
        let sk = SecretKey::random(&params, &mut rng);
        let pk = PublicKey::new(&sk, &mut rng);
        let pt = Plaintext::try_encode(&[1u64][..], Encoding::poly(), &other_params)?;

        assert!(pk.try_encrypt(&pt, &mut rng).is_err());
        Ok(())
    }

    #[test]
    fn test_serialize() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for params in [
            BfvParameters::default_arc(1, 16),
            BfvParameters::default_arc(6, 16),
        ] {
            let sk = SecretKey::random(&params, &mut rng);
            let pk = PublicKey::new(&sk, &mut rng);
            let bytes = pk.to_bytes();
            assert_eq!(pk, PublicKey::from_bytes(&bytes, &params)?);
        }
        Ok(())
    }
}
