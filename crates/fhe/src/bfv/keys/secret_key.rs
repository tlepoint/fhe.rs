//! Secret keys for the BFV encryption scheme

use crate::bfv::{
    BfvParameters, Ciphertext, Plaintext, parameters::PlaintextModulus, plaintext::PlaintextValues,
};
use crate::proto::bfv::SecretKey as SecretKeyProto;
use crate::{Error, Result, SerializationError};
use fhe_math::{
    rq::{Poly, Representation, traits::TryConvertFrom},
    zq::Modulus,
};
use fhe_traits::{DeserializeParametrized, FheDecrypter, FheEncrypter, FheParametrized, Serialize};
use fhe_util::sample_vec_cbd;
use itertools::Itertools;
use num_bigint::BigUint;
use prost::Message;
use rand::{CryptoRng, Rng, RngCore, SeedableRng};
use rand_chacha::ChaCha8Rng;
use std::sync::Arc;
use zeroize::{Zeroize, Zeroizing};

/// Secret key for the BFV encryption scheme.
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct SecretKey {
    pub(crate) par: Arc<BfvParameters>,
    pub(crate) coeffs: Box<[i64]>,
}

impl Zeroize for SecretKey {
    fn zeroize(&mut self) {
        // Only zeroize the sensitive coefficients field
        self.coeffs.zeroize();
    }
}

impl Drop for SecretKey {
    fn drop(&mut self) {
        self.zeroize();
    }
}

impl SecretKey {
    /// Generate a random [`SecretKey`].
    pub fn random<R: RngCore + CryptoRng>(par: &Arc<BfvParameters>, rng: &mut R) -> Self {
        let s_coefficients = sample_vec_cbd(par.degree(), par.variance, rng).unwrap();
        Self::new(s_coefficients, par)
    }

    /// Generate a [`SecretKey`] from its coefficients.
    pub(crate) fn new(coeffs: Vec<i64>, par: &Arc<BfvParameters>) -> Self {
        Self {
            par: par.to_owned(),
            coeffs: coeffs.into_boxed_slice(),
        }
    }

    /// Measure the noise in a [`Ciphertext`].
    ///
    /// # Safety
    ///
    /// This operations may run in a variable time depending on the value of the
    /// noise.
    pub unsafe fn measure_noise(&self, ct: &Ciphertext) -> Result<usize> {
        let plaintext = Zeroizing::new(self.try_decrypt(ct)?);
        let m = Zeroizing::new(plaintext.to_poly());

        // Let's create a secret key with the ciphertext context
        let mut s = Zeroizing::new(Poly::try_convert_from(
            self.coeffs.as_ref(),
            ct[0].ctx(),
            false,
            Representation::PowerBasis,
        )?);
        s.change_representation(Representation::Ntt);
        let mut si = s.clone();

        // Let's disable variable time computations
        let mut c = Zeroizing::new(ct[0].clone());
        c.disallow_variable_time_computations();

        for i in 1..ct.len() {
            let mut cis = Zeroizing::new(ct[i].clone());
            cis.disallow_variable_time_computations();
            *cis.as_mut() *= si.as_ref();
            *c.as_mut() += &cis;
            *si.as_mut() *= s.as_ref();
        }
        *c.as_mut() -= &m;
        c.change_representation(Representation::PowerBasis);

        let ciphertext_modulus = ct[0].ctx().modulus();
        let mut noise = 0usize;
        for coeff in Vec::<BigUint>::from(c.as_ref()) {
            noise = std::cmp::max(
                noise,
                std::cmp::min(coeff.bits(), (ciphertext_modulus - &coeff).bits()) as usize,
            )
        }

        Ok(noise)
    }

    pub(crate) fn encrypt_poly<R: RngCore + CryptoRng>(
        &self,
        p: &Poly,
        rng: &mut R,
    ) -> Result<Ciphertext> {
        assert_eq!(p.representation(), &Representation::Ntt);

        let level = self.par.level_of_context(p.ctx())?;

        let mut seed = <ChaCha8Rng as SeedableRng>::Seed::default();
        rand::rng().fill(&mut seed);

        // Let's create a secret key with the ciphertext context
        let mut s = Zeroizing::new(Poly::try_convert_from(
            self.coeffs.as_ref(),
            p.ctx(),
            false,
            Representation::PowerBasis,
        )?);
        s.change_representation(Representation::Ntt);

        let mut a = Poly::random_from_seed(p.ctx(), Representation::Ntt, seed);
        let a_s = Zeroizing::new(&a * s.as_ref());

        let mut b = Poly::small(p.ctx(), Representation::Ntt, self.par.variance, rng)
            .map_err(Error::MathError)?;
        b -= &a_s;
        b += p;

        // It is now safe to enable variable time computations.
        unsafe {
            a.allow_variable_time_computations();
            b.allow_variable_time_computations()
        }

        Ok(Ciphertext {
            par: self.par.clone(),
            seed: Some(seed),
            c: vec![b, a],
            level,
        })
    }
}

impl From<&SecretKey> for SecretKeyProto {
    fn from(sk: &SecretKey) -> Self {
        Self {
            coeffs: sk.coeffs.to_vec(),
        }
    }
}

impl Serialize for SecretKey {
    fn to_bytes(&self) -> Vec<u8> {
        SecretKeyProto::from(self).encode_to_vec()
    }
}

impl DeserializeParametrized for SecretKey {
    type Error = Error;

    fn from_bytes(bytes: &[u8], par: &Arc<Self::Parameters>) -> Result<Self> {
        let proto: SecretKeyProto = Message::decode(bytes).map_err(|_| {
            Error::SerializationError(SerializationError::ProtobufError {
                message: "SecretKey decode".into(),
            })
        })?;

        if proto.coeffs.len() != par.degree() {
            return Err(Error::SerializationError(
                SerializationError::InvalidFormat {
                    reason: "SecretKey coeffs length and parameters degree mismatch".into(),
                },
            ));
        }

        Ok(Self {
            par: par.clone(),
            coeffs: proto.coeffs.into_boxed_slice(),
        })
    }
}

impl FheParametrized for SecretKey {
    type Parameters = BfvParameters;
}

impl FheEncrypter<Plaintext, Ciphertext> for SecretKey {
    type Error = Error;

    fn try_encrypt<R: RngCore + CryptoRng>(
        &self,
        pt: &Plaintext,
        rng: &mut R,
    ) -> Result<Ciphertext> {
        assert!(Arc::ptr_eq(&self.par, &pt.par));
        let m = Zeroizing::new(pt.to_poly());
        self.encrypt_poly(m.as_ref(), rng)
    }
}

impl FheDecrypter<Plaintext, Ciphertext> for SecretKey {
    type Error = Error;

    fn try_decrypt(&self, ct: &Ciphertext) -> Result<Plaintext> {
        if !Arc::ptr_eq(&self.par, &ct.par) {
            Err(Error::DefaultError(
                "Incompatible BFV parameters".to_string(),
            ))
        } else {
            // Let's create a secret key with the ciphertext context
            let mut s = Zeroizing::new(Poly::try_convert_from(
                self.coeffs.as_ref(),
                ct[0].ctx(),
                false,
                Representation::PowerBasis,
            )?);
            s.change_representation(Representation::Ntt);
            let mut si = s.clone();

            let mut c = Zeroizing::new(ct[0].clone());
            c.disallow_variable_time_computations();

            // Compute the phase c0 + c1*s + c2*s^2 + ... where the secret power
            // s^k is computed on-the-fly
            for i in 1..ct.len() {
                let mut cis = Zeroizing::new(ct[i].clone());
                cis.disallow_variable_time_computations();
                *cis.as_mut() *= si.as_ref();
                *c.as_mut() += &cis;
                if i + 1 < ct.len() {
                    *si.as_mut() *= s.as_ref();
                }
            }
            c.change_representation(Representation::PowerBasis);

            let ctx_lvl = self.par.context_level_at(ct.level).unwrap();
            let d = Zeroizing::new(c.scale(&ctx_lvl.cipher_plain_context.scaler)?);

            let value = match self.par.plaintext {
                PlaintextModulus::Small(_) => {
                    let mut v = Vec::<u64>::try_from(d.as_ref())?;
                    let plaintext_modulus = self.par.plaintext();
                    v.iter_mut().for_each(|vi| *vi += plaintext_modulus);
                    let mut w = v[..self.par.degree()].to_vec();

                    let q = Modulus::new(self.par.moduli[0]).map_err(Error::MathError)?;
                    q.reduce_vec(&mut w);
                    if let PlaintextModulus::Small(ref m) = self.par.plaintext {
                        m.reduce_vec(&mut w);
                    }
                    PlaintextValues::Small(w.into_boxed_slice())
                }
                PlaintextModulus::Large(_) => {
                    let v: Vec<BigUint> = Vec::<BigUint>::from(d.as_ref())
                        .into_iter()
                        .map(|vi| vi + self.par.plaintext_big())
                        .collect_vec();

                    let mut w = v[..self.par.degree()].to_vec();
                    let q_poly = d.as_ref().ctx().modulus();
                    w.iter_mut().for_each(|wi| *wi %= q_poly);

                    self.par.plaintext.reduce_vec(&mut w);
                    PlaintextValues::Large(w.into_boxed_slice())
                }
            };

            let mut poly = match &value {
                PlaintextValues::Small(v) => Poly::try_convert_from(
                    v.as_ref(),
                    ct[0].ctx(),
                    false,
                    Representation::PowerBasis,
                )?,
                PlaintextValues::Large(v) => Poly::try_convert_from(
                    v.as_ref(),
                    ct[0].ctx(),
                    false,
                    Representation::PowerBasis,
                )?,
            };

            poly.change_representation(Representation::Ntt);

            let pt = Plaintext {
                par: self.par.clone(),
                value,
                encoding: None,
                poly_ntt: poly,
                level: ct.level,
            };

            Ok(pt)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::SecretKey;
    use crate::bfv::{Encoding, Plaintext, parameters::BfvParameters};
    use crate::proto::bfv::SecretKey as SecretKeyProto;
    use fhe_traits::{DeserializeParametrized, FheDecrypter, FheEncoder, FheEncrypter, Serialize};
    use prost::Message;
    use rand::rng;
    use std::error::Error;

    #[test]
    fn keygen() {
        let mut rng = rng();
        let params = BfvParameters::default_arc(1, 16);
        let sk = SecretKey::random(&params, &mut rng);
        assert_eq!(sk.par, params);

        sk.coeffs.iter().for_each(|ci| {
            // Check that this is a small polynomial
            assert!((*ci).abs() <= 2 * sk.par.variance as i64)
        })
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
                    let q = fhe_math::zq::Modulus::new(params.plaintext()).unwrap();

                    let pt = Plaintext::try_encode(
                        &q.random_vec(params.degree(), &mut rng),
                        Encoding::poly_at_level(level),
                        &params,
                    )?;
                    let ct = sk.try_encrypt(&pt, &mut rng)?;
                    let pt2 = sk.try_decrypt(&ct)?;

                    println!("Noise: {}", unsafe { sk.measure_noise(&ct)? });
                    assert_eq!(pt2, pt);
                }
            }
        }

        Ok(())
    }

    #[test]
    fn serialize_roundtrip() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let params = BfvParameters::default_arc(2, 16);
        let sk = SecretKey::random(&params, &mut rng);

        let bytes = sk.to_bytes();
        let decoded = SecretKey::from_bytes(&bytes, &params)?;

        assert_eq!(decoded, sk);
        Ok(())
    }

    #[test]
    fn deserialize_invalid_length() {
        let params = BfvParameters::default_arc(1, 16);
        let mut proto = SecretKeyProto {
            coeffs: vec![0; params.degree()],
        };
        proto.coeffs.pop();

        let bytes = proto.encode_to_vec();
        let err = SecretKey::from_bytes(&bytes, &params).unwrap_err();

        assert!(matches!(
            err,
            crate::Error::SerializationError(crate::SerializationError::InvalidFormat { .. })
        ));
    }
}
