//! Secret keys for the BFV encryption scheme

use crate::bfv::{Ciphertext, Parameters, Plaintext};
use crate::proto::bfv::SecretKey as SecretKeyProto;
use crate::{Error, Result, SerializationError};
use fhe_math::{
    rq::{Ntt, Poly, PowerBasis, traits::TryConvertFrom},
    zq::Modulus,
};

use fhe_util::sample_vec_cbd;
use itertools::Itertools;
use num_bigint::BigUint;
use prost::Message;
use rand::{CryptoRng, Rng as RngCore, RngExt, SeedableRng};
use rand_chacha::ChaCha8Rng;
use zeroize::{Zeroize, Zeroizing};

/// Secret key for the BFV encryption scheme.
#[derive(PartialEq, Eq, Clone)]
pub struct SecretKey {
    pub(crate) par: Parameters,
    pub(crate) coeffs: Box<[i64]>,
}

impl std::fmt::Debug for SecretKey {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SecretKey")
            .field("degree", &self.par.degree())
            .finish_non_exhaustive()
    }
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
    pub fn generate<R: RngCore + CryptoRng>(par: &Parameters, rng: &mut R) -> Self {
        let s_coefficients = sample_vec_cbd(par.degree(), par.inner.variance, rng).unwrap();
        Self::new(s_coefficients, par)
    }

    /// Generate a [`SecretKey`] from its coefficients.
    pub(crate) fn new(coeffs: Vec<i64>, par: &Parameters) -> Self {
        Self {
            par: par.to_owned(),
            coeffs: coeffs.into_boxed_slice(),
        }
    }

    /// Measure the noise in a [`Ciphertext`].
    ///
    /// ```compile_fail
    /// use fhe::bfv::{SecretKey, Ciphertext};
    /// fn diagnostic(sk: &SecretKey, ct: &Ciphertext) {
    ///     sk.measure_noise_vartime(ct);
    /// }
    /// ```
    ///
    /// Both the result and the running time reveal secret-dependent noise.
    /// Use only in a trusted diagnostic setting, never as an oracle for
    /// untrusted callers. The acknowledgment concerns information leakage,
    /// not Rust memory safety.
    pub fn measure_noise_vartime(
        &self,
        ct: &Ciphertext,
        _diagnostics: crate::SecretDependentDiagnostics,
    ) -> Result<usize> {
        let plaintext = Zeroizing::new(self.decrypt(ct)?);
        let m = Zeroizing::new(plaintext.to_poly());

        // Let's create a secret key with the ciphertext context
        let s = Zeroizing::new(
            Poly::<PowerBasis>::try_convert_from(self.coeffs.as_ref(), ct.c[0].ctx())?.into_ntt(),
        );
        let mut si = s.clone();

        // Let's disable variable time computations
        let mut c = Zeroizing::new(ct.c[0].clone());
        c.disallow_variable_time_computations();

        for i in 1..ct.len() {
            let mut cis = Zeroizing::new(ct.c[i].clone());
            cis.disallow_variable_time_computations();
            *cis.as_mut() *= si.as_ref();
            *c.as_mut() += &cis;
            *si.as_mut() *= s.as_ref();
        }
        *c.as_mut() -= &m;
        let ctx = c.ctx().clone();
        let c_inner = std::mem::replace(c.as_mut(), Poly::<Ntt>::zero(&ctx));
        let c = Zeroizing::new(c_inner.into_power_basis());

        let ciphertext_modulus = ct.c[0].ctx().modulus();
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
        p: &Poly<Ntt>,
        rng: &mut R,
    ) -> Result<Ciphertext> {
        let level = self.par.level_of_context(p.ctx())?;

        let mut seed = <ChaCha8Rng as SeedableRng>::Seed::default();
        rng.fill(&mut seed);

        // Let's create a secret key with the ciphertext context
        let s = Zeroizing::new(
            Poly::<PowerBasis>::try_convert_from(self.coeffs.as_ref(), p.ctx())?.into_ntt(),
        );

        let mut a = Poly::<Ntt>::random_from_seed(p.ctx(), seed);
        let a_s = Zeroizing::new(&a * s.as_ref());

        let mut b =
            Poly::<Ntt>::small(p.ctx(), self.par.inner.variance, rng).map_err(Error::MathError)?;
        b -= &a_s;
        b += p;

        // It is now safe to enable variable time computations.
        let variable_time = crate::VariableTime::new(crate::PublicData::assert_public());
        a.allow_variable_time_computations(variable_time);
        b.allow_variable_time_computations(variable_time);

        Ok(Ciphertext {
            par: self.par.clone(),
            seed: Some(seed),
            c: vec![b, a],
            level,
        })
    }
}

// Also clears partially decoded coefficients when protobuf parsing fails.
impl Drop for SecretKeyProto {
    fn drop(&mut self) {
        self.coeffs.zeroize();
    }
}

impl From<&SecretKey> for SecretKeyProto {
    fn from(sk: &SecretKey) -> Self {
        Self {
            coeffs: sk.coeffs.to_vec(),
        }
    }
}

impl SecretKey {
    /// Export unencrypted secret material in the existing protobuf format.
    /// The returned bytes and temporary coefficient copies are zeroized on
    /// drop.
    #[must_use]
    pub fn export_secret_bytes(&self) -> Zeroizing<Vec<u8>> {
        Zeroizing::new(SecretKeyProto::from(self).encode_to_vec())
    }
}

impl SecretKey {
    /// Import validated protobuf bytes, binding contextual values to the
    /// supplied parameters.
    pub fn from_bytes(bytes: &[u8], par: &Parameters) -> Result<Self> {
        let mut proto: SecretKeyProto = Message::decode(bytes).map_err(|_| {
            Error::SerializationError(SerializationError::Decode {
                object: crate::SerializedObject::SecretKey,
            })
        })?;

        if proto.coeffs.len() != par.degree() {
            return Err(Error::SerializationError(
                SerializationError::InvalidSecretKeyCoefficientCount {
                    actual: proto.coeffs.len(),
                    expected: par.degree(),
                },
            ));
        }

        Ok(Self {
            par: par.clone(),
            coeffs: std::mem::take(&mut proto.coeffs).into_boxed_slice(),
        })
    }
}

impl SecretKey {
    /// Encrypt a plaintext with compatible parameters using caller-owned
    /// cryptographic randomness.
    pub fn encrypt<R: RngCore + CryptoRng>(
        &self,
        pt: &Plaintext,
        rng: &mut R,
    ) -> Result<Ciphertext> {
        pt.validate_for(&self.par)?;
        let m = Zeroizing::new(pt.to_poly());
        self.encrypt_poly(m.as_ref(), rng)
    }
}

impl SecretKey {
    /// Decrypt a ciphertext with compatible parameters at its current level.
    pub fn decrypt(&self, ct: &Ciphertext) -> Result<Plaintext> {
        ct.validate_for(&self.par)?;
        // Let's create a secret key with the ciphertext context
        let s = Zeroizing::new(
            Poly::<PowerBasis>::try_convert_from(self.coeffs.as_ref(), ct.c[0].ctx())?.into_ntt(),
        );
        let mut si = s.clone();

        let mut c = Zeroizing::new(ct.c[0].clone());
        c.disallow_variable_time_computations();

        // Compute the phase c0 + c1*s + c2*s^2 + ... where the secret power
        // s^k is computed on-the-fly
        for i in 1..ct.len() {
            let mut cis = Zeroizing::new(ct.c[i].clone());
            cis.disallow_variable_time_computations();
            *cis.as_mut() *= si.as_ref();
            *c.as_mut() += &cis;
            if i + 1 < ct.len() {
                *si.as_mut() *= s.as_ref();
            }
        }
        let ctx_lvl = self.par.context_level_at(ct.level)?;
        let ctx = c.ctx().clone();
        let c_inner = std::mem::replace(c.as_mut(), Poly::<Ntt>::zero(&ctx));
        let c_pb = Zeroizing::new(c_inner.into_power_basis());
        let d = Zeroizing::new(c_pb.as_ref().scale(&ctx_lvl.cipher_plain_context.scaler)?);

        let poly = match self.par.inner.plaintext.small() {
            Some(plaintext_modulus) if **plaintext_modulus < self.par.inner.moduli[0] => {
                let mut v = Vec::<u64>::try_from(d.as_ref())?;
                v.iter_mut().for_each(|vi| *vi += **plaintext_modulus);
                let mut w = v[..self.par.degree()].to_vec();

                let q = Modulus::new(self.par.inner.moduli[0]).map_err(Error::MathError)?;
                q.reduce_vec(&mut w);
                plaintext_modulus.reduce_vec(&mut w);
                Poly::<PowerBasis>::try_convert_from(w.as_slice(), ct.c[0].ctx())?.into_ntt()
            }
            Some(_) | None => {
                // A single residue cannot recover values modulo t when t is
                // larger than q0, even if t itself fits in a machine word.
                let v: Vec<BigUint> = Vec::<BigUint>::from(d.as_ref())
                    .into_iter()
                    .map(|vi| vi + self.par.plaintext_modulus())
                    .collect_vec();

                let mut w = v[..self.par.degree()].to_vec();
                let q_poly = d.as_ref().ctx().modulus();
                w.iter_mut().for_each(|wi| *wi %= q_poly);

                self.par.inner.plaintext.reduce_vec(&mut w);
                Poly::<PowerBasis>::try_convert_from(w.as_slice(), ct.c[0].ctx())?.into_ntt()
            }
        };

        let pt = Plaintext {
            par: self.par.clone(),

            poly_ntt: poly,
        };

        Ok(pt)
    }
}

#[cfg(test)]
mod tests {
    use super::SecretKey;
    use crate::bfv::{Encoding, Plaintext, parameters::Parameters};
    use crate::proto::bfv::SecretKey as SecretKeyProto;

    use prost::Message;
    use rand::rng;
    use std::error::Error;

    #[test]
    fn keygen() {
        let mut rng = rng();
        let params = Parameters::test_parameters(1, 16);
        let sk = SecretKey::generate(&params, &mut rng);
        assert_eq!(sk.par, params);

        sk.coeffs.iter().for_each(|ci| {
            // Check that this is a small polynomial
            assert!((*ci).abs() <= 2 * sk.par.inner.variance as i64)
        })
    }

    #[test]
    fn decrypt_word_plaintext_larger_than_first_prime() -> Result<(), Box<dyn Error>> {
        use crate::bfv::{Ciphertext, ParametersBuilder, PublicKey};

        let t = 1u64 << 40;
        let params = ParametersBuilder::new()
            .degree(16)
            .plaintext_modulus(t)
            .ciphertext_modulus_bits([30, 30, 30, 30])
            .build()?;
        let mut rng = rng();
        let sk = SecretKey::generate(&params, &mut rng);
        let pk = PublicKey::from_secret_key(&sk, &mut rng);
        let mut values = vec![0u64, 1, 12345, t / 2, t - 1];
        values.resize(params.degree(), 0);
        // Retain enough ciphertext modulus for a meaningful noise budget.
        for level in [0, 1] {
            let encoding = Encoding::Polynomial;
            let pt = Plaintext::encode_at_level(&params, &values, encoding, level)?;
            for ct in [sk.encrypt(&pt, &mut rng)?, pk.encrypt(&pt, &mut rng)?] {
                let ct: Ciphertext = ct;
                assert_eq!(sk.decrypt(&ct)?.decode(encoding)?, values);
                let doubled = ct.add(&ct).unwrap();
                let expected = values.iter().map(|v| (v * 2) % t).collect::<Vec<_>>();
                assert_eq!(sk.decrypt(&doubled)?.decode(encoding)?, expected);
            }
        }
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
                    let q = fhe_math::zq::Modulus::new(params.plaintext_modulus_u64().unwrap())
                        .unwrap();

                    let pt = Plaintext::encode_at_level(
                        &params,
                        &q.random_vec(params.degree(), &mut rng),
                        Encoding::Polynomial,
                        level,
                    )?;
                    let ct = sk.encrypt(&pt, &mut rng)?;
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
    fn encrypt_decrypt_reject_invalid_inputs() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let params = Parameters::test_parameters(1, 16);
        let other_params = Parameters::test_parameters(1, 32);
        let sk = SecretKey::generate(&params, &mut rng);
        let other_pt = Plaintext::encode(&other_params, &[1u64][..], Encoding::Polynomial)?;
        let encrypted: crate::Result<crate::bfv::Ciphertext> = sk.encrypt(&other_pt, &mut rng);

        assert!(encrypted.is_err());
        assert!(matches!(
            sk.decrypt(&crate::bfv::Ciphertext::invalid_empty(&params)),
            Err(crate::Error::Ciphertext(_))
        ));
        Ok(())
    }

    #[test]
    fn measure_noise_within_modulus_bits() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let params = Parameters::test_parameters(1, 16);
        let sk = SecretKey::generate(&params, &mut rng);
        let q = fhe_math::zq::Modulus::new(params.plaintext_modulus_u64().unwrap()).unwrap();

        let pt = Plaintext::encode(
            &params,
            &q.random_vec(params.degree(), &mut rng),
            Encoding::Polynomial,
        )?;
        let ct = sk.encrypt(&pt, &mut rng)?;
        let noise = sk.measure_noise_vartime(
            &ct,
            crate::SecretDependentDiagnostics::acknowledge_leakage(),
        )?;

        let modulus_bits = ct.c[0].ctx().modulus().bits() as usize;
        assert!(noise <= modulus_bits);

        Ok(())
    }

    #[test]
    fn serialize_roundtrip() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let params = Parameters::test_parameters(2, 16);
        let sk = SecretKey::generate(&params, &mut rng);

        let bytes = sk.export_secret_bytes();
        let decoded = SecretKey::from_bytes(&bytes, &params)?;

        assert_eq!(decoded, sk);
        Ok(())
    }

    #[test]
    fn deserialize_invalid_length() {
        let params = Parameters::test_parameters(1, 16);
        let mut proto = SecretKeyProto {
            coeffs: vec![0; params.degree()],
        };
        proto.coeffs.pop();

        let bytes = proto.encode_to_vec();
        let err = SecretKey::from_bytes(&bytes, &params).unwrap_err();

        assert!(matches!(
            err,
            crate::Error::SerializationError(
                crate::SerializationError::InvalidSecretKeyCoefficientCount { .. }
            )
        ));
    }

    #[test]
    fn encryption_uses_only_the_supplied_rng() -> Result<(), Box<dyn Error>> {
        use rand::{Rng, SeedableRng};
        use rand_chacha::ChaCha8Rng;
        let par = Parameters::test_parameters(2, 16);
        let sk = SecretKey::generate(&par, &mut rng());
        let pt = Plaintext::encode(&par, &[42u64], Encoding::Polynomial)?;
        let mut first = ChaCha8Rng::seed_from_u64(123);
        let mut second = ChaCha8Rng::seed_from_u64(123);
        let a: crate::bfv::Ciphertext = sk.encrypt(&pt, &mut first)?;
        let b: crate::bfv::Ciphertext = sk.encrypt(&pt, &mut second)?;
        assert_eq!(a.to_bytes(), b.to_bytes());
        assert_eq!(first.next_u64(), second.next_u64());
        let next: crate::bfv::Ciphertext = sk.encrypt(&pt, &mut first)?;
        assert_ne!(a.to_bytes(), next.to_bytes());
        assert_eq!(sk.decrypt(&a)?, sk.decrypt(&next)?);
        Ok(())
    }
}
