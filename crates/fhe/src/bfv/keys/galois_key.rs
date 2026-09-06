//! Galois keys for the BFV encryption scheme

use super::key_switching_key::KeySwitchingKey;
use crate::bfv::{BfvParameters, Ciphertext, SecretKey, traits::TryConvertFrom};
use crate::proto::bfv::{GaloisKey as GaloisKeyProto, KeySwitchingKey as KeySwitchingKeyProto};
use crate::{Error, Result, SerializationError};
use fhe_math::rq::{
    Ntt, Poly, PowerBasis, SubstitutionExponent, switcher::Switcher,
    traits::TryConvertFrom as TryConvertFromPoly,
};
use rand::{CryptoRng, Rng as RngCore};
use std::sync::Arc;
use zeroize::{Zeroize, Zeroizing};

/// Galois key for the BFV encryption scheme.
/// A Galois key is a special type of key switching key,
/// which switch from `s(x^i)` to `s(x)` where `s(x)` is the secret key.
#[derive(Debug, PartialEq, Eq)]
pub struct GaloisKey {
    pub(crate) element: SubstitutionExponent,
    pub(crate) ksk: KeySwitchingKey,
}

impl GaloisKey {
    /// Generate a [`GaloisKey`] from a [`SecretKey`].
    pub fn new<R: RngCore + CryptoRng>(
        sk: &SecretKey,
        exponent: usize,
        ciphertext_level: usize,
        galois_key_level: usize,
        rng: &mut R,
    ) -> Result<Self> {
        let ctx_galois_key = sk.par.context_at_level(galois_key_level)?;
        let ctx_ciphertext = sk.par.context_at_level(ciphertext_level)?;

        let ciphertext_exponent =
            SubstitutionExponent::new(ctx_ciphertext, exponent).map_err(Error::MathError)?;

        let switcher_up = Switcher::new(ctx_ciphertext, ctx_galois_key)?;
        let s = Zeroizing::new(Poly::<PowerBasis>::try_convert_from(
            sk.coeffs.as_ref(),
            ctx_ciphertext,
            false,
        )?);
        let s_sub = Zeroizing::new(s.substitute(&ciphertext_exponent)?);
        let s_sub_switched_up = Zeroizing::new(s_sub.switch(&switcher_up)?);

        let ksk = KeySwitchingKey::new(
            sk,
            &s_sub_switched_up,
            ciphertext_level,
            galois_key_level,
            rng,
        )?;

        Ok(Self {
            element: ciphertext_exponent,
            ksk,
        })
    }

    /// Relinearize a [`Ciphertext`] using the [`GaloisKey`]
    pub fn relinearize(&self, ct: &Ciphertext) -> Result<Ciphertext> {
        self.validate_ciphertext(ct)?;

        let c2 = ct[1].substitute(&self.element)?;
        let (mut c0, mut c1) = self.ksk.key_switch_ntt(c2)?;

        if c0.ctx() != ct[0].ctx() {
            c0.switch_down_to(ct[0].ctx())?;
            c1.switch_down_to(ct[1].ctx())?;
        }

        c0 += &ct[0].substitute(&self.element)?;

        Ok(Ciphertext {
            par: ct.par.clone(),
            seed: None,
            c: vec![c0, c1],
            level: self.ksk.ciphertext_level,
        })
    }

    /// Relinearize a [`Ciphertext`] writing the result into `out`.
    pub fn relinearize_into(&self, ct: &Ciphertext, out: &mut Ciphertext) -> Result<()> {
        self.validate_ciphertext(ct)?;

        if out.len() != 2 || out[0].ctx() != ct[0].ctx() || out[1].ctx() != ct[1].ctx() {
            out.c = vec![
                Poly::<Ntt>::zero(ct[0].ctx()),
                Poly::<Ntt>::zero(ct[1].ctx()),
            ];
        }
        out.par = ct.par.clone();
        out.seed = None;
        out.level = self.ksk.ciphertext_level;

        let (out0_slice, out1_slice) = out.split_at_mut(1);
        let out0 = &mut out0_slice[0];
        let out1 = &mut out1_slice[0];

        out0.zeroize();
        out1.zeroize();

        let c2 = ct[1].substitute(&self.element)?;
        self.ksk.key_switch_ntt_assign(c2, out0, out1)?;

        if out0.ctx() != ct[0].ctx() {
            out0.switch_down_to(ct[0].ctx())?;
            out1.switch_down_to(ct[1].ctx())?;
        }

        *out0 += &ct[0].substitute(&self.element)?;
        Ok(())
    }

    fn validate_ciphertext(&self, ct: &Ciphertext) -> Result<()> {
        ct.validate_for(&self.ksk.par)?;
        if ct.len() != 2 {
            return Err(crate::CiphertextError::InvalidPolynomialCount {
                operation: crate::CiphertextOperation::Galois,
                actual: ct.len(),
                expected: 2,
            }
            .into());
        }
        if ct.level != self.ksk.ciphertext_level {
            return Err(Error::InvalidLevel {
                level: ct.level,
                min_level: self.ksk.ciphertext_level,
                max_level: self.ksk.ciphertext_level,
            });
        }
        Ok(())
    }
}

impl From<&GaloisKey> for GaloisKeyProto {
    fn from(value: &GaloisKey) -> Self {
        GaloisKeyProto {
            exponent: value.element.exponent as u32,
            ksk: Some(KeySwitchingKeyProto::from(&value.ksk)),
        }
    }
}

impl TryConvertFrom<&GaloisKeyProto> for GaloisKey {
    fn try_convert_from(value: &GaloisKeyProto, par: &Arc<BfvParameters>) -> Result<Self> {
        if let Some(ksk) = &value.ksk {
            let ksk = KeySwitchingKey::try_convert_from(ksk, par)?;

            let ctx = par.context_at_level(ksk.ciphertext_level)?;
            let element = SubstitutionExponent::new(ctx, value.exponent as usize)
                .map_err(Error::MathError)?;

            Ok(GaloisKey { element, ksk })
        } else {
            Err(Error::SerializationError(
                SerializationError::MissingField {
                    field: crate::SerializedField::GaloisKeySwitchingKey,
                },
            ))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::GaloisKey;
    use crate::bfv::{
        BfvParameters, Ciphertext, Encoding, Plaintext, SecretKey, traits::TryConvertFrom,
    };
    use crate::proto::bfv::GaloisKey as GaloisKeyProto;
    use fhe_traits::{FheDecoder, FheDecrypter, FheEncoder, FheEncrypter};
    use rand::rng;
    use std::error::Error;

    #[test]
    fn relinearization() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for params in [
            BfvParameters::default_arc(6, 16),
            BfvParameters::default_arc(3, 16),
        ] {
            for _ in 0..30 {
                let sk = SecretKey::random(&params, &mut rng);
                let v = fhe_math::zq::Modulus::new(params.plaintext())
                    .unwrap()
                    .random_vec(params.degree(), &mut rng);
                let row_size = params.degree() >> 1;

                let pt = Plaintext::try_encode(&v, Encoding::simd(), &params)?;
                let ct = sk.try_encrypt(&pt, &mut rng)?;

                for i in 1..2 * params.degree() {
                    if i & 1 == 0 {
                        assert!(GaloisKey::new(&sk, i, 0, 0, &mut rng).is_err())
                    } else {
                        let gk = GaloisKey::new(&sk, i, 0, 0, &mut rng)?;
                        let ct2 = gk.relinearize(&ct)?;
                        println!("Noise: {}", unsafe { sk.measure_noise(&ct2)? });

                        if i == 3 {
                            let pt = sk.try_decrypt(&ct2)?;

                            // The expected result is rotated one on the left
                            let mut expected = vec![0u64; params.degree()];
                            expected[..row_size - 1].copy_from_slice(&v[1..row_size]);
                            expected[row_size - 1] = v[0];
                            expected[row_size..2 * row_size - 1]
                                .copy_from_slice(&v[row_size + 1..]);
                            expected[2 * row_size - 1] = v[row_size];
                            assert_eq!(&Vec::<u64>::try_decode(&pt, Encoding::simd())?, &expected)
                        } else if i == params.degree() * 2 - 1 {
                            let pt = sk.try_decrypt(&ct2)?;

                            // The expected result has its rows swapped
                            let mut expected = vec![0u64; params.degree()];
                            expected[..row_size].copy_from_slice(&v[row_size..]);
                            expected[row_size..].copy_from_slice(&v[..row_size]);
                            assert_eq!(&Vec::<u64>::try_decode(&pt, Encoding::simd())?, &expected)
                        }
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn relinearization_into() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for params in [
            BfvParameters::default_arc(6, 16),
            BfvParameters::default_arc(3, 16),
        ] {
            let sk = SecretKey::random(&params, &mut rng);
            let pt = Plaintext::try_encode(&[1u64, 2, 3, 4][..], Encoding::simd(), &params)?;
            let ct = sk.try_encrypt(&pt, &mut rng)?;
            let gk = GaloisKey::new(&sk, 3, 0, 0, &mut rng)?;

            let ct_expected = gk.relinearize(&ct)?;

            let mut out = Ciphertext::zero(&ct.par);
            gk.relinearize_into(&ct, &mut out)?;

            assert_eq!(ct_expected, out);
        }
        Ok(())
    }

    #[test]
    fn leveled_relinearization_matches_power_basis_switching() -> Result<(), Box<dyn Error>> {
        use rand::SeedableRng;
        use rand_chacha::ChaCha8Rng;

        let params = BfvParameters::default_arc(4, 16);
        let mut rng = ChaCha8Rng::seed_from_u64(0x6a1015);
        let sk = SecretKey::random(&params, &mut rng);
        for level in 0..=params.max_level() {
            let pt = Plaintext::try_encode(
                &[1u64, 2, 3, 4][..],
                Encoding::simd_at_level(level),
                &params,
            )?;
            let ct: Ciphertext = sk.try_encrypt(&pt, &mut rng)?;
            for key_level in 0..=level {
                for exponent in [3, params.degree() + 1] {
                    let gk = GaloisKey::new(&sk, exponent, level, key_level, &mut rng)?;
                    let mut out = Ciphertext::zero(&params);
                    // Reuse the output across permission changes and check
                    // both allocating and in-place entry points against the
                    // previous transform/switch/transform implementation.
                    for public in [true, false, true] {
                        let mut input = ct.clone();
                        if !public {
                            input[1].disallow_variable_time_computations();
                        }
                        let c2 = input[1].substitute(&gk.element)?.into_power_basis();
                        let (c0, c1) = gk.ksk.key_switch(&c2)?;
                        let mut c0 = c0.into_power_basis();
                        let mut c1 = c1.into_power_basis();
                        c0.switch_down_to(input[0].ctx())?;
                        c1.switch_down_to(input[1].ctx())?;
                        let mut c0 = c0.into_ntt();
                        c0 += &input[0].substitute(&gk.element)?;
                        let expected = Ciphertext::new(vec![c0, c1.into_ntt()], &params)?;
                        assert_eq!(gk.relinearize(&input)?, expected);
                        gk.relinearize_into(&input, &mut out)?;
                        assert_eq!(out, expected);
                        assert!(
                            out.iter()
                                .all(|p| p.allows_variable_time_computations() == public)
                        );
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn relinearization_rejects_invalid_ciphertexts() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let params = BfvParameters::default_arc(3, 16);
        let sk = SecretKey::random(&params, &mut rng);
        let gk = GaloisKey::new(&sk, 3, 0, 0, &mut rng)?;
        let invalid = Ciphertext::zero(&params);

        assert!(matches!(
            gk.relinearize(&invalid),
            Err(crate::Error::Ciphertext(_))
        ));
        let mut out = Ciphertext::zero(&params);
        assert!(matches!(
            gk.relinearize_into(&invalid, &mut out),
            Err(crate::Error::Ciphertext(_))
        ));
        Ok(())
    }

    #[test]
    fn proto_conversion() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for params in [
            BfvParameters::default_arc(6, 16),
            BfvParameters::default_arc(4, 16),
        ] {
            let sk = SecretKey::random(&params, &mut rng);
            let gk = GaloisKey::new(&sk, 9, 0, 0, &mut rng)?;
            let proto = GaloisKeyProto::from(&gk);
            assert_eq!(gk, GaloisKey::try_convert_from(&proto, &params)?);
        }
        Ok(())
    }
}
