//! Relinearization keys for the BFV encryption scheme

use super::key_switching_key::KeySwitchingKey;
use crate::bfv::{Ciphertext, Parameters, SecretKey, wire::FromProto};
use crate::proto::bfv::{
    KeySwitchingKey as KeySwitchingKeyProto, RelinearizationKey as RelinearizationKeyProto,
};
use crate::{Error, Result, error::SerializationError};
use fhe_math::rq::{Ntt, Poly, PowerBasis, switcher::Switcher};

use prost::Message;
use rand::{CryptoRng, Rng as RngCore};
use zeroize::Zeroizing;

/// Relinearization key for the BFV encryption scheme.
/// A relinearization key is a special type of key switching key,
/// which switch from `s^2` to `s` where `s` is the secret key.
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct RelinearizationKey {
    pub(crate) ksk: KeySwitchingKey,
}

impl RelinearizationKey {
    /// Generate a [`RelinearizationKey`] from a [`SecretKey`].
    pub fn new<R: RngCore + CryptoRng>(sk: &SecretKey, rng: &mut R) -> Result<Self> {
        Self::new_leveled_internal(sk, 0, 0, rng)
    }

    /// Configure distinct ciphertext and evaluation-key levels before
    /// generation.
    #[must_use]
    pub fn builder(sk: &SecretKey) -> RelinearizationKeyBuilder<'_> {
        RelinearizationKeyBuilder {
            sk,
            ciphertext_level: 0,
            key_level: 0,
        }
    }

    fn new_leveled_internal<R: RngCore + CryptoRng>(
        sk: &SecretKey,
        ciphertext_level: usize,
        key_level: usize,
        rng: &mut R,
    ) -> Result<Self> {
        let ctx_ciphertext = sk.par.context_at_level(ciphertext_level)?;
        if key_level > ciphertext_level {
            return Err(Error::InvalidLevel {
                level: key_level,
                min_level: 0,
                max_level: ciphertext_level,
            });
        }
        let ctx_relin_key = sk.par.context_at_level(key_level)?;

        if ctx_relin_key.moduli().len() == 1 {
            return Err(crate::error::EvaluationKeyError::KeySwitchingNotSupported.into());
        }

        let s = Zeroizing::new(
            Poly::<PowerBasis>::from_signed_coefficients(sk.coeffs.as_ref(), ctx_ciphertext)?
                .into_ntt(),
        );
        let s2 = Zeroizing::new((s.as_ref() * s.as_ref()).into_power_basis());
        let switcher_up = Switcher::new(ctx_ciphertext, ctx_relin_key)?;
        let s2_switched_up = Zeroizing::new(s2.switch(&switcher_up)?);
        let ksk = KeySwitchingKey::new(sk, &s2_switched_up, ciphertext_level, key_level, rng)?;
        Ok(Self { ksk })
    }

    /// Relinearizes the supplied `(c0, c1, c2)` ciphertext in place, reducing
    /// it to two components.
    pub fn relinearize(&self, ct: &mut Ciphertext) -> Result<()> {
        ct.validate_for(&self.ksk.par)?;
        if ct.len() != 3 {
            Err(crate::error::CiphertextError::InvalidPolynomialCount {
                operation: crate::error::CiphertextOperation::Relinearization,
                actual: ct.len(),
                expected: 3,
            }
            .into())
        } else if ct.level != self.ksk.ciphertext_level {
            Err(Error::InvalidLevel {
                level: ct.level,
                min_level: self.ksk.ciphertext_level,
                max_level: self.ksk.ciphertext_level,
            })
        } else {
            let c2 = ct.c[2].clone().into_power_basis();
            let (mut c0, mut c1) = self.relinearizes_poly(&c2)?;

            if c0.ctx() != ct.c[0].ctx() {
                let mut c0_pb = c0.into_power_basis();
                let mut c1_pb = c1.into_power_basis();
                c0_pb.switch_down_to(ct.c[0].ctx())?;
                c1_pb.switch_down_to(ct.c[1].ctx())?;
                c0 = c0_pb.into_ntt();
                c1 = c1_pb.into_ntt();
            }

            ct.c[0] += &c0;
            ct.c[1] += &c1;
            ct.truncate(2);
            Ok(())
        }
    }

    /// Relinearize using polynomials.
    pub(crate) fn relinearizes_poly(
        &self,
        c2: &Poly<PowerBasis>,
    ) -> Result<(Poly<Ntt>, Poly<Ntt>)> {
        self.ksk.key_switch(c2)
    }
}

/// Consuming configuration for a relinearization key, borrowing its secret.
/// Both levels default to zero. Validation happens before randomness is used.
#[derive(Debug, Clone)]
pub struct RelinearizationKeyBuilder<'key> {
    sk: &'key SecretKey,
    ciphertext_level: usize,
    key_level: usize,
}

impl RelinearizationKeyBuilder<'_> {
    /// Level of the ciphertexts this key will relinearize.
    #[must_use]
    pub fn ciphertext_level(mut self, level: usize) -> Self {
        self.ciphertext_level = level;
        self
    }

    /// Level at which the key-switching material is stored.
    /// Must not exceed the ciphertext level; a one-modulus key is unsupported.
    #[must_use]
    pub fn key_level(mut self, level: usize) -> Self {
        self.key_level = level;
        self
    }

    /// Validate the levels and generate the key using caller-owned randomness.
    pub fn build<R: RngCore + CryptoRng>(self, rng: &mut R) -> Result<RelinearizationKey> {
        RelinearizationKey::new_leveled_internal(
            self.sk,
            self.ciphertext_level,
            self.key_level,
            rng,
        )
    }
}

impl From<&RelinearizationKey> for RelinearizationKeyProto {
    fn from(value: &RelinearizationKey) -> Self {
        RelinearizationKeyProto {
            ksk: Some(KeySwitchingKeyProto::from(&value.ksk)),
        }
    }
}

impl FromProto<&RelinearizationKeyProto> for RelinearizationKey {
    fn from_proto(
        value: &RelinearizationKeyProto,
        par: &Parameters,
        limits: &crate::DecodeLimits,
    ) -> Result<Self> {
        if let Some(ksk) = &value.ksk {
            Ok(RelinearizationKey {
                ksk: KeySwitchingKey::from_proto(ksk, par, limits)?,
            })
        } else {
            Err(Error::SerializationError(
                SerializationError::MissingField {
                    field: crate::error::SerializedField::RelinearizationKeySwitchingKey,
                },
            ))
        }
    }
}

impl RelinearizationKey {
    /// Serialize in the existing protobuf wire format.
    #[must_use]
    pub fn to_bytes(&self) -> Vec<u8> {
        RelinearizationKeyProto::from(self).encode_to_vec()
    }
}

impl RelinearizationKey {
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
            crate::error::SerializedObject::RelinearizationKey,
            Some(par),
            limits,
        )?;
        let rk = Message::decode(bytes).map_err(|source| {
            Error::SerializationError(SerializationError::Decode {
                object: crate::error::SerializedObject::RelinearizationKey,
                source,
            })
        })?;
        RelinearizationKey::from_proto(&rk, par, limits)
    }
}

#[cfg(test)]
mod tests {
    use super::RelinearizationKey;
    use crate::bfv::{Ciphertext, Encoding, Parameters, SecretKey, wire::FromProto};
    use crate::proto::bfv::RelinearizationKey as RelinearizationKeyProto;
    use fhe_math::rq::{Ntt, Poly, PowerBasis};

    use rand::rng;
    use std::error::Error;

    #[test]
    fn relinearization_discards_a_serialized_last_component_seed() -> Result<(), Box<dyn Error>> {
        let params = Parameters::test_parameters(2, 16);
        let mut rng = rng();
        let sk = SecretKey::generate(&params, &mut rng);
        let rk = RelinearizationKey::new(&sk, &mut rng)?;
        let ctx = params.context_at_level(0)?;
        let seed = [23; 32];
        let mut ct = Ciphertext::from_components(
            vec![
                Poly::zero(ctx),
                Poly::zero(ctx),
                Poly::random_from_seed(ctx, seed),
            ],
            &params,
        )?;
        ct.seed = Some(seed);
        let mut restored = Ciphertext::from_bytes(&ct.to_bytes(), &params)?;
        assert!(restored.seed.is_some());
        rk.relinearize(&mut restored)?;
        assert!(restored.seed.is_none());
        let round_trip = Ciphertext::from_bytes(&restored.to_bytes(), &params)?;
        assert_eq!(restored, round_trip);
        assert_eq!(restored.to_bytes(), round_trip.to_bytes());
        Ok(())
    }

    #[test]
    fn relinearization() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for params in [Parameters::test_parameters(6, 16)] {
            for _ in 0..100 {
                let sk = SecretKey::generate(&params, &mut rng);
                let rk = RelinearizationKey::new(&sk, &mut rng)?;

                let ctx = params.context_at_level(0)?;
                let s = Poly::<PowerBasis>::from_signed_coefficients(sk.coeffs.as_ref(), ctx)
                    .map_err(crate::Error::MathError)?
                    .into_ntt();
                let s2 = &s * &s;

                // Let's generate manually an "extended" ciphertext (c0 = e - c1 * s - c2 * s^2,
                // c1, c2) encrypting 0.
                let c2 = Poly::<Ntt>::random(ctx, &mut rng);
                let c1 = Poly::<Ntt>::random(ctx, &mut rng);
                let mut c0 = Poly::<PowerBasis>::small(ctx, 16, &mut rng)?.into_ntt();
                c0 -= &(&c1 * &s);
                c0 -= &(&c2 * &s2);
                let mut ct =
                    Ciphertext::from_components(vec![c0.clone(), c1.clone(), c2.clone()], &params)?;

                // Relinearize the extended ciphertext!
                rk.relinearize(&mut ct)?;
                assert_eq!(ct.len(), 2);

                // Check that the relinearization by polynomials works the same way
                let c2_pb = c2.clone().into_power_basis();
                let (c0r, c1r) = rk.relinearizes_poly(&c2_pb)?;
                let mut c0r_pb = c0r.into_power_basis();
                c0r_pb.switch_down_to(c0.ctx())?;
                let mut c1r_pb = c1r.into_power_basis();
                c1r_pb.switch_down_to(c1.ctx())?;
                let c0r = c0r_pb.into_ntt();
                let c1r = c1r_pb.into_ntt();
                assert_eq!(
                    ct,
                    Ciphertext::from_components(vec![&c0 + &c0r, &c1 + &c1r], &params)?
                );

                // Print the noise and decrypt
                println!(
                    "Noise: {}",
                    sk.measure_noise_vartime(
                        &ct,
                        crate::SecretDependentDiagnostics::acknowledge_leakage()
                    )?
                );
                let pt = sk.decrypt(&ct)?;
                let w = pt.decode(Encoding::Polynomial)?;
                assert_eq!(w, &[0u64; 16]);
            }
        }
        Ok(())
    }

    #[test]
    fn relinearization_leveled() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for params in [Parameters::test_parameters(5, 16)] {
            for ciphertext_level in 0..params.max_level() {
                for key_level in 0..=ciphertext_level {
                    for _ in 0..10 {
                        let sk = SecretKey::generate(&params, &mut rng);
                        let rk = RelinearizationKey::builder(&sk)
                            .ciphertext_level(ciphertext_level)
                            .key_level(key_level)
                            .build(&mut rng)?;

                        let ctx = params.context_at_level(ciphertext_level)?;
                        let s =
                            Poly::<PowerBasis>::from_signed_coefficients(sk.coeffs.as_ref(), ctx)
                                .map_err(crate::Error::MathError)?
                                .into_ntt();
                        let s2 = &s * &s;
                        // Let's generate manually an "extended" ciphertext (c0 = e - c1 * s - c2 *
                        // s^2, c1, c2) encrypting 0.
                        let c2 = Poly::<Ntt>::random(ctx, &mut rng);
                        let c1 = Poly::<Ntt>::random(ctx, &mut rng);
                        let mut c0 = Poly::<PowerBasis>::small(ctx, 16, &mut rng)?.into_ntt();
                        c0 -= &(&c1 * &s);
                        c0 -= &(&c2 * &s2);
                        let mut ct = Ciphertext::from_components(
                            vec![c0.clone(), c1.clone(), c2.clone()],
                            &params,
                        )?;

                        // Relinearize the extended ciphertext!
                        rk.relinearize(&mut ct)?;
                        assert_eq!(ct.len(), 2);

                        // Check that the relinearization by polynomials works the same way
                        let c2_pb = c2.clone().into_power_basis();
                        let (c0r, c1r) = rk.relinearizes_poly(&c2_pb)?;
                        let mut c0r_pb = c0r.into_power_basis();
                        c0r_pb.switch_down_to(c0.ctx())?;
                        let mut c1r_pb = c1r.into_power_basis();
                        c1r_pb.switch_down_to(c1.ctx())?;
                        let c0r = c0r_pb.into_ntt();
                        let c1r = c1r_pb.into_ntt();
                        assert_eq!(
                            ct,
                            Ciphertext::from_components(vec![&c0 + &c0r, &c1 + &c1r], &params)?
                        );

                        // Print the noise and decrypt
                        println!(
                            "Noise: {}",
                            sk.measure_noise_vartime(
                                &ct,
                                crate::SecretDependentDiagnostics::acknowledge_leakage()
                            )?
                        );
                        let pt = sk.decrypt(&ct)?;
                        let w = pt.decode(Encoding::Polynomial)?;
                        assert_eq!(w, &[0u64; 16]);
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn proto_conversion() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for params in [
            Parameters::test_parameters(6, 16),
            Parameters::test_parameters(3, 16),
        ] {
            let sk = SecretKey::generate(&params, &mut rng);
            let rk = RelinearizationKey::new(&sk, &mut rng)?;
            let proto = RelinearizationKeyProto::from(&rk);
            assert_eq!(
                rk,
                RelinearizationKey::from_proto(&proto, &params, &crate::DecodeLimits::default())?
            );
        }
        Ok(())
    }
}
