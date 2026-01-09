//! Relinearization keys for the BFV encryption scheme

use std::sync::Arc;

use super::key_switching_key::KeySwitchingKey;
use crate::bfv::{BfvParameters, Ciphertext, SecretKey, traits::TryConvertFrom};
use crate::proto::bfv::{
    KeySwitchingKey as KeySwitchingKeyProto, RelinearizationKey as RelinearizationKeyProto,
};
use crate::{Error, Result};
use fhe_math::rq::{
    Poly, Representation, switcher::Switcher, traits::TryConvertFrom as TryConvertFromPoly,
};
use fhe_traits::{DeserializeParametrized, FheParametrized, Serialize};
use prost::Message;
use rand::{CryptoRng, RngCore};
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

    /// Generate a [`RelinearizationKey`] from a [`SecretKey`].
    pub fn new_leveled<R: RngCore + CryptoRng>(
        sk: &SecretKey,
        ciphertext_level: usize,
        key_level: usize,
        rng: &mut R,
    ) -> Result<Self> {
        Self::new_leveled_internal(sk, ciphertext_level, key_level, rng)
    }

    fn new_leveled_internal<R: RngCore + CryptoRng>(
        sk: &SecretKey,
        ciphertext_level: usize,
        key_level: usize,
        rng: &mut R,
    ) -> Result<Self> {
        let ctx_relin_key = sk.par.context_at_level(key_level)?;
        let ctx_ciphertext = sk.par.context_at_level(ciphertext_level)?;

        if ctx_relin_key.moduli().len() == 1 {
            return Err(Error::DefaultError(
                "These parameters do not support key switching".to_string(),
            ));
        }

        let mut s = Zeroizing::new(Poly::try_convert_from(
            sk.coeffs.as_ref(),
            ctx_ciphertext,
            false,
            Representation::PowerBasis,
        )?);
        s.change_representation(Representation::Ntt);
        let mut s2 = Zeroizing::new(s.as_ref() * s.as_ref());
        s2.change_representation(Representation::PowerBasis);
        let switcher_up = Switcher::new(ctx_ciphertext, ctx_relin_key)?;
        let s2_switched_up = Zeroizing::new(s2.switch(&switcher_up)?);
        let ksk = KeySwitchingKey::new(sk, &s2_switched_up, ciphertext_level, key_level, rng)?;
        Ok(Self { ksk })
    }

    /// Relinearizes the supplied `(c0, c1, c2)` ciphertext in place, reducing
    /// it to two components.
    pub fn relinearizes(&self, ct: &mut Ciphertext) -> Result<()> {
        if ct.len() != 3 {
            Err(Error::DefaultError(
                "Only supports relinearization of ciphertext with 3 parts".to_string(),
            ))
        } else if ct.level != self.ksk.ciphertext_level {
            Err(Error::DefaultError(
                "Ciphertext has incorrect level".to_string(),
            ))
        } else {
            let mut c2 = ct[2].clone();
            c2.change_representation(Representation::PowerBasis);

            #[allow(unused_mut)]
            let (mut c0, mut c1) = self.relinearizes_poly(&c2)?;

            if c0.ctx() != ct[0].ctx() {
                c0.change_representation(Representation::PowerBasis);
                c1.change_representation(Representation::PowerBasis);
                c0.switch_down_to(ct[0].ctx())?;
                c1.switch_down_to(ct[1].ctx())?;
                c0.change_representation(Representation::Ntt);
                c1.change_representation(Representation::Ntt);
            }

            ct[0] += &c0;
            ct[1] += &c1;
            ct.truncate(2);
            Ok(())
        }
    }

    /// Relinearize using polynomials.
    pub(crate) fn relinearizes_poly(&self, c2: &Poly) -> Result<(Poly, Poly)> {
        self.ksk.key_switch(c2)
    }
}

impl From<&RelinearizationKey> for RelinearizationKeyProto {
    fn from(value: &RelinearizationKey) -> Self {
        RelinearizationKeyProto {
            ksk: Some(KeySwitchingKeyProto::from(&value.ksk)),
        }
    }
}

impl TryConvertFrom<&RelinearizationKeyProto> for RelinearizationKey {
    fn try_convert_from(value: &RelinearizationKeyProto, par: &Arc<BfvParameters>) -> Result<Self> {
        if value.ksk.is_some() {
            Ok(RelinearizationKey {
                ksk: KeySwitchingKey::try_convert_from(value.ksk.as_ref().unwrap(), par)?,
            })
        } else {
            Err(Error::DefaultError("Invalid serialization".to_string()))
        }
    }
}

impl Serialize for RelinearizationKey {
    fn to_bytes(&self) -> Vec<u8> {
        RelinearizationKeyProto::from(self).encode_to_vec()
    }
}

impl FheParametrized for RelinearizationKey {
    type Parameters = BfvParameters;
}

impl DeserializeParametrized for RelinearizationKey {
    type Error = Error;

    fn from_bytes(bytes: &[u8], par: &Arc<Self::Parameters>) -> Result<Self> {
        let rk = Message::decode(bytes);
        if let Ok(rk) = rk {
            RelinearizationKey::try_convert_from(&rk, par)
        } else {
            Err(Error::DefaultError("Invalid serialization".to_string()))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::RelinearizationKey;
    use crate::bfv::{BfvParameters, Ciphertext, Encoding, SecretKey, traits::TryConvertFrom};
    use crate::proto::bfv::RelinearizationKey as RelinearizationKeyProto;
    use fhe_math::rq::{Poly, Representation, traits::TryConvertFrom as TryConvertFromPoly};
    use fhe_traits::{FheDecoder, FheDecrypter};
    use rand::rng;
    use std::error::Error;

    #[test]
    fn relinearization() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for params in [BfvParameters::default_arc(6, 16)] {
            for _ in 0..100 {
                let sk = SecretKey::random(&params, &mut rng);
                let rk = RelinearizationKey::new(&sk, &mut rng)?;

                let ctx = params.context_at_level(0)?;
                let mut s = Poly::try_convert_from(
                    sk.coeffs.as_ref(),
                    ctx,
                    false,
                    Representation::PowerBasis,
                )
                .map_err(crate::Error::MathError)?;
                s.change_representation(Representation::Ntt);
                let s2 = &s * &s;

                // Let's generate manually an "extended" ciphertext (c0 = e - c1 * s - c2 * s^2,
                // c1, c2) encrypting 0.
                let mut c2 = Poly::random(ctx, Representation::Ntt, &mut rng);
                let c1 = Poly::random(ctx, Representation::Ntt, &mut rng);
                let mut c0 = Poly::small(ctx, Representation::PowerBasis, 16, &mut rng)?;
                c0.change_representation(Representation::Ntt);
                c0 -= &(&c1 * &s);
                c0 -= &(&c2 * &s2);
                let mut ct = Ciphertext::new(vec![c0.clone(), c1.clone(), c2.clone()], &params)?;

                // Relinearize the extended ciphertext!
                rk.relinearizes(&mut ct)?;
                assert_eq!(ct.len(), 2);

                // Check that the relinearization by polynomials works the same way
                c2.change_representation(Representation::PowerBasis);
                let (mut c0r, mut c1r) = rk.relinearizes_poly(&c2)?;
                c0r.change_representation(Representation::PowerBasis);
                c0r.switch_down_to(c0.ctx())?;
                c1r.change_representation(Representation::PowerBasis);
                c1r.switch_down_to(c1.ctx())?;
                c0r.change_representation(Representation::Ntt);
                c1r.change_representation(Representation::Ntt);
                assert_eq!(ct, Ciphertext::new(vec![&c0 + &c0r, &c1 + &c1r], &params)?);

                // Print the noise and decrypt
                println!("Noise: {}", unsafe { sk.measure_noise(&ct)? });
                let pt = sk.try_decrypt(&ct)?;
                let w = Vec::<u64>::try_decode(&pt, Encoding::poly())?;
                assert_eq!(w, &[0u64; 16]);
            }
        }
        Ok(())
    }

    #[test]
    fn relinearization_leveled() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for params in [BfvParameters::default_arc(5, 16)] {
            for ciphertext_level in 0..params.max_level() {
                for key_level in 0..=ciphertext_level {
                    for _ in 0..10 {
                        let sk = SecretKey::random(&params, &mut rng);
                        let rk = RelinearizationKey::new_leveled(
                            &sk,
                            ciphertext_level,
                            key_level,
                            &mut rng,
                        )?;

                        let ctx = params.context_at_level(ciphertext_level)?;
                        let mut s = Poly::try_convert_from(
                            sk.coeffs.as_ref(),
                            ctx,
                            false,
                            Representation::PowerBasis,
                        )
                        .map_err(crate::Error::MathError)?;
                        s.change_representation(Representation::Ntt);
                        let s2 = &s * &s;
                        // Let's generate manually an "extended" ciphertext (c0 = e - c1 * s - c2 *
                        // s^2, c1, c2) encrypting 0.
                        let mut c2 = Poly::random(ctx, Representation::Ntt, &mut rng);
                        let c1 = Poly::random(ctx, Representation::Ntt, &mut rng);
                        let mut c0 = Poly::small(ctx, Representation::PowerBasis, 16, &mut rng)?;
                        c0.change_representation(Representation::Ntt);
                        c0 -= &(&c1 * &s);
                        c0 -= &(&c2 * &s2);
                        let mut ct =
                            Ciphertext::new(vec![c0.clone(), c1.clone(), c2.clone()], &params)?;

                        // Relinearize the extended ciphertext!
                        rk.relinearizes(&mut ct)?;
                        assert_eq!(ct.len(), 2);

                        // Check that the relinearization by polynomials works the same way
                        c2.change_representation(Representation::PowerBasis);
                        let (mut c0r, mut c1r) = rk.relinearizes_poly(&c2)?;
                        c0r.change_representation(Representation::PowerBasis);
                        c0r.switch_down_to(c0.ctx())?;
                        c1r.change_representation(Representation::PowerBasis);
                        c1r.switch_down_to(c1.ctx())?;
                        c0r.change_representation(Representation::Ntt);
                        c1r.change_representation(Representation::Ntt);
                        assert_eq!(ct, Ciphertext::new(vec![&c0 + &c0r, &c1 + &c1r], &params)?);

                        // Print the noise and decrypt
                        println!("Noise: {}", unsafe { sk.measure_noise(&ct)? });
                        let pt = sk.try_decrypt(&ct)?;
                        let w = Vec::<u64>::try_decode(&pt, Encoding::poly())?;
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
            BfvParameters::default_arc(6, 16),
            BfvParameters::default_arc(3, 16),
        ] {
            let sk = SecretKey::random(&params, &mut rng);
            let rk = RelinearizationKey::new(&sk, &mut rng)?;
            let proto = RelinearizationKeyProto::from(&rk);
            assert_eq!(rk, RelinearizationKey::try_convert_from(&proto, &params)?);
        }
        Ok(())
    }
}
