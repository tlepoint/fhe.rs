//! Galois keys for the BFV encryption scheme

use super::key_switching_key::KeySwitchingKey;
use crate::bfv::{traits::TryConvertFrom, BfvParameters, Ciphertext, SecretKey};
use crate::proto::bfv::{GaloisKey as GaloisKeyProto, KeySwitchingKey as KeySwitchingKeyProto};
use crate::{Error, Result};
use fhe_math::rq::{
    switcher::Switcher, traits::TryConvertFrom as TryConvertFromPoly, Poly, Representation,
    SubstitutionExponent,
};
use rand::{CryptoRng, RngCore};
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
        let s = Zeroizing::new(Poly::try_convert_from(
            sk.coeffs.as_ref(),
            ctx_ciphertext,
            false,
            Representation::PowerBasis,
        )?);
        let s_sub = Zeroizing::new(s.substitute(&ciphertext_exponent)?);
        let mut s_sub_switched_up = Zeroizing::new(s_sub.switch(&switcher_up)?);
        s_sub_switched_up.change_representation(Representation::PowerBasis);

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
        // assert_eq!(ct.par, self.ksk.par);
        assert_eq!(ct.len(), 2);

        let mut c2 = ct[1].substitute(&self.element)?;
        c2.change_representation(Representation::PowerBasis);
        let (mut c0, mut c1) = self.ksk.key_switch(&c2)?;

        if c0.ctx() != ct[0].ctx() {
            c0.change_representation(Representation::PowerBasis);
            c1.change_representation(Representation::PowerBasis);
            c0.switch_down_to(ct[0].ctx())?;
            c1.switch_down_to(ct[1].ctx())?;
            c0.change_representation(Representation::Ntt);
            c1.change_representation(Representation::Ntt);
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
        assert_eq!(ct.len(), 2);

        if out.len() != 2 || out[0].ctx() != ct[0].ctx() || out[1].ctx() != ct[1].ctx() {
            out.c = vec![
                Poly::zero(ct[0].ctx(), Representation::Ntt),
                Poly::zero(ct[1].ctx(), Representation::Ntt),
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

        let mut c2 = ct[1].substitute(&self.element)?;
        c2.change_representation(Representation::PowerBasis);
        self.ksk.key_switch_assign(&c2, out0, out1)?;

        if out0.ctx() != ct[0].ctx() {
            out0.change_representation(Representation::PowerBasis);
            out1.change_representation(Representation::PowerBasis);
            out0.switch_down_to(ct[0].ctx())?;
            out1.switch_down_to(ct[1].ctx())?;
            out0.change_representation(Representation::Ntt);
            out1.change_representation(Representation::Ntt);
        }

        *out0 += &ct[0].substitute(&self.element)?;
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
        if value.ksk.is_some() {
            let ksk = KeySwitchingKey::try_convert_from(value.ksk.as_ref().unwrap(), par)?;

            let ctx = par.context_at_level(ksk.ciphertext_level)?;
            let element = SubstitutionExponent::new(ctx, value.exponent as usize)
                .map_err(Error::MathError)?;

            Ok(GaloisKey { element, ksk })
        } else {
            Err(Error::DefaultError("Invalid serialization".to_string()))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::GaloisKey;
    use crate::bfv::{
        traits::TryConvertFrom, BfvParameters, Ciphertext, Encoding, Plaintext, SecretKey,
    };
    use crate::proto::bfv::GaloisKey as GaloisKeyProto;
    use fhe_traits::{FheDecoder, FheDecrypter, FheEncoder, FheEncrypter};
    use rand::thread_rng;
    use std::error::Error;

    #[test]
    fn relinearization() -> Result<(), Box<dyn Error>> {
        let mut rng = thread_rng();
        for params in [
            BfvParameters::default_arc(6, 16),
            BfvParameters::default_arc(3, 16),
        ] {
            for _ in 0..30 {
                let sk = SecretKey::random(&params, &mut rng);
                let v = params.plaintext.random_vec(params.degree(), &mut rng);
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
        let mut rng = thread_rng();
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
    fn proto_conversion() -> Result<(), Box<dyn Error>> {
        let mut rng = thread_rng();
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
