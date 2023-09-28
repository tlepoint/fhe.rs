use std::sync::Arc;

use fhe_math::rq::traits::TryConvertFrom;
use fhe_math::rq::{Poly, Representation};
use fhe_math::zq::Modulus;
use itertools::Itertools;
use rand::{CryptoRng, RngCore};
use zeroize::Zeroizing;

use crate::bfv::{BfvParameters, Ciphertext, Plaintext, SecretKey};
use crate::mbfv::Aggregate;
use crate::{Error, Result};

/// Each party uses the `SecretKeySwitchShare` to generate their share of the new ciphertext and
/// participate in the "Protocol 3: KeySwitch" protocol detailed in Multiparty BFV (p7).
///
/// Note: it appears the MBFV paper assumes the output key is split into the same number of parties
/// as the input key?
pub struct SecretKeySwitchShare {
    pub(crate) par: Arc<BfvParameters>,
    /// The original input ciphertext
    // Probably doesn't need to be Arc in real usage but w/e
    pub(crate) ct: Arc<Ciphertext>,
    pub(crate) h_share: Poly,
}

impl SecretKeySwitchShare {
    /// Participate in a new KeySwitch protocol
    ///
    /// 1. *Private input*: BFV input secret key share
    /// 2. *Private input*: BFV output secret key share
    /// 3. *Public input*: Ciphertext
    /// 4. *Public input*: TODO: variance of the ciphertext noise
    pub fn new<R: RngCore + CryptoRng>(
        sk_input_share: &SecretKey,
        sk_output_share: &SecretKey,
        ct: Arc<Ciphertext>,
        rng: &mut R,
    ) -> Result<Self> {
        if sk_input_share.par != sk_output_share.par || sk_output_share.par != ct.par {
            return Err(Error::DefaultError(
                "Incompatible BFV parameters".to_string(),
            ));
        }
        let par = sk_input_share.par.clone();
        let mut s_in = Zeroizing::new(Poly::try_convert_from(
            sk_input_share.coeffs.as_ref(),
            ct.c[0].ctx(),
            false,
            Representation::PowerBasis,
        )?);
        s_in.change_representation(Representation::Ntt);
        let mut s_out = Zeroizing::new(Poly::try_convert_from(
            sk_output_share.coeffs.as_ref(),
            ct.c[0].ctx(),
            false,
            Representation::PowerBasis,
        )?);
        s_out.change_representation(Representation::Ntt);

        // Sample error
        // TODO this should be exponential in ciphertext noise!
        let e = Zeroizing::new(Poly::small(
            ct.c[0].ctx(),
            Representation::Ntt,
            par.variance,
            rng,
        )?);

        // Create h_i share
        // TODO look at SecretKey::try_decrypt, probably need the `for i in 1..ct.c.len()` loop here.
        // although I do think this is correct for len == 2, so for now:
        assert_eq!(ct.c.len(), 2);
        let mut h_share = s_in.as_ref() - s_out.as_ref();
        h_share.disallow_variable_time_computations();
        h_share *= &ct.c[1];
        h_share += e.as_ref();

        Ok(Self { par, ct, h_share })
    }
}

impl Aggregate for SecretKeySwitchShare {
    type Output = Ciphertext;

    fn aggregate<I>(shares: I) -> Result<Self::Output>
    where
        I: IntoIterator<Item = Self>,
    {
        let mut shares = shares.into_iter();
        let share = shares.next().ok_or(Error::TooFewValues(0, 1))?;
        let mut h = share.h_share;
        for sh in shares {
            h += &sh.h_share;
        }

        let c0 = &share.ct.c[0] + &h;
        let c1 = share.ct.c[1].clone();

        Ciphertext::new(vec![c0, c1], &share.par)
    }
}

/// Each party uses the `DecryptionShare` to generate their share of the decrypted ciphertext.
///
/// This is the same thing as the `SecretKeySwitchShare` protocol, just with an output key of zero.
pub struct DecryptionShare {
    pub(crate) sks_share: SecretKeySwitchShare,
}

impl DecryptionShare {
    /// Participate in a new Decryption protocol (i.e. KeySwitch protocol with a zero output key).
    ///
    /// 1. *Private input*: BFV input secret key share
    /// 3. *Public input*: Ciphertext
    /// 4. *Public input*: TODO: variance of the ciphertext noise
    pub fn new<R: RngCore + CryptoRng>(
        sk_input_share: &SecretKey,
        ct: &Arc<Ciphertext>,
        rng: &mut R,
    ) -> Result<Self> {
        let par = &sk_input_share.par;
        let zero = SecretKey::new(vec![0; par.degree()], par);
        let sks_share = SecretKeySwitchShare::new(sk_input_share, &zero, ct.clone(), rng)?;
        Ok(DecryptionShare { sks_share })
    }
}

impl Aggregate for DecryptionShare {
    type Output = Plaintext;

    fn aggregate<I>(shares: I) -> Result<Self::Output>
    where
        I: IntoIterator<Item = Self>,
    {
        let sks_shares = shares.into_iter().map(|s| s.sks_share);
        let ct = SecretKeySwitchShare::aggregate(sks_shares)?;
        let par = ct.par;

        // Note: during SKS, c[1]*sk has already been added to c[0].
        let mut c = Zeroizing::new(ct.c[0].clone());
        c.disallow_variable_time_computations();
        c.change_representation(Representation::PowerBasis);

        // The true decryption part is done during SKS; all that is left is to scale
        let d = Zeroizing::new(c.scale(&par.scalers[ct.level])?);
        let v = Zeroizing::new(
            Vec::<u64>::from(d.as_ref())
                .iter_mut()
                .map(|vi| *vi + par.plaintext.modulus())
                .collect_vec(),
        );
        let mut w = v[..par.degree()].to_vec();
        let q = Modulus::new(par.moduli[0]).map_err(Error::MathError)?;
        q.reduce_vec(&mut w);
        par.plaintext.reduce_vec(&mut w);

        let mut poly =
            Poly::try_convert_from(&w, ct.c[0].ctx(), false, Representation::PowerBasis)?;
        poly.change_representation(Representation::Ntt);

        let pt = Plaintext {
            par: par.clone(),
            value: w.into_boxed_slice(),
            encoding: None,
            poly_ntt: poly,
            level: ct.level,
        };

        Ok(pt)
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use fhe_math::rq::{Poly, Representation};
    use fhe_traits::{FheEncoder, FheEncrypter};
    use rand::thread_rng;

    use crate::{
        bfv::{BfvParameters, Encoding, Plaintext, SecretKey},
        mbfv::{protocols::PublicKeyShare, Aggregate},
    };

    use super::*;

    const NUM_PARTIES: usize = 11;

    struct Party {
        sk_share: SecretKey,
        pk_share: PublicKeyShare,
    }

    #[test]
    fn encrypt_decrypt() {
        let mut rng = thread_rng();
        for par in [
            BfvParameters::default_arc(1, 8),
            BfvParameters::default_arc(6, 8),
        ] {
            for level in 0..=par.max_level() {
                for _ in 0..20 {
                    let crp =
                        Poly::random(par.ctx_at_level(0).unwrap(), Representation::Ntt, &mut rng);

                    let mut parties: Vec<Party> = vec![];

                    // Parties collectively generate public key
                    for _ in 0..NUM_PARTIES {
                        let sk_share = SecretKey::random(&par, &mut rng);
                        let pk_share =
                            PublicKeyShare::new(&sk_share, crp.clone(), &mut rng).unwrap();
                        parties.push(Party { sk_share, pk_share })
                    }
                    let public_key =
                        PublicKeyShare::aggregate(parties.iter().map(|p| p.pk_share.clone()))
                            .unwrap();

                    // Use it to encrypt a random polynomial
                    let pt1 = Plaintext::try_encode(
                        &par.plaintext.random_vec(par.degree(), &mut rng),
                        Encoding::poly_at_level(level),
                        &par,
                    )
                    .unwrap();
                    let ct = Arc::new(public_key.try_encrypt(&pt1, &mut rng).unwrap());

                    // Parties perform a collective decryption
                    let decryption_shares = parties
                        .iter()
                        .map(|p| DecryptionShare::new(&p.sk_share, &ct, &mut rng).unwrap());
                    let pt2 = DecryptionShare::aggregate(decryption_shares).unwrap();

                    assert_eq!(pt1, pt2);
                }
            }
        }
    }

    #[test]
    fn encrypt_keyswitch_decrypt() {
        let mut rng = thread_rng();
        for par in [
            BfvParameters::default_arc(1, 8),
            BfvParameters::default_arc(6, 8),
        ] {
            for level in 0..=par.max_level() {
                for _ in 0..20 {
                    let crp =
                        Poly::random(par.ctx_at_level(0).unwrap(), Representation::Ntt, &mut rng);

                    // Parties collectively generate public key
                    let mut parties: Vec<Party> = vec![];
                    for _ in 0..NUM_PARTIES {
                        let sk_share = SecretKey::random(&par, &mut rng);
                        let pk_share =
                            PublicKeyShare::new(&sk_share, crp.clone(), &mut rng).unwrap();
                        parties.push(Party { sk_share, pk_share })
                    }

                    let public_key =
                        PublicKeyShare::aggregate(parties.iter().map(|p| p.pk_share.clone()))
                            .unwrap();

                    // Use it to encrypt a random polynomial ct1
                    let pt1 = Plaintext::try_encode(
                        &[1u64],
                        // &par.plaintext.random_vec(par.degree(), &mut rng),
                        Encoding::poly_at_level(level),
                        &par,
                    )
                    .unwrap();
                    let ct1 = Arc::new(public_key.try_encrypt(&pt1, &mut rng).unwrap());

                    // Key switch ct1 to a different set of parties
                    let mut out_parties = Vec::new();
                    for _ in 0..NUM_PARTIES {
                        let sk_share = SecretKey::random(&par, &mut rng);
                        let pk_share =
                            PublicKeyShare::new(&sk_share, crp.clone(), &mut rng).unwrap();
                        out_parties.push(Party { sk_share, pk_share })
                    }
                    let skss: Vec<SecretKeySwitchShare> = parties
                        .iter()
                        .zip(out_parties.iter())
                        .map(|(ip, op)| {
                            SecretKeySwitchShare::new(
                                &ip.sk_share,
                                &op.sk_share,
                                ct1.clone(),
                                &mut rng,
                            )
                            .unwrap()
                        })
                        .collect();
                    let ct2 = Arc::new(SecretKeySwitchShare::aggregate(skss).unwrap());

                    // The second set of parties then does a collective decryption
                    let decryption_shares = out_parties
                        .iter()
                        .map(|p| DecryptionShare::new(&p.sk_share, &ct2, &mut rng).unwrap());

                    let pt2 = DecryptionShare::aggregate(decryption_shares).unwrap();
                    assert_eq!(pt1, pt2);
                }
            }
        }
    }
}
