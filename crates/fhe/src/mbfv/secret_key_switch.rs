use std::sync::Arc;

use fhe_math::rq::{Poly, Representation, traits::TryConvertFrom};
use itertools::Itertools;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use rand::{CryptoRng, RngCore};
use zeroize::Zeroizing;

use crate::bfv::{BfvParameters, Ciphertext, Plaintext, PlaintextValues, SecretKey};
use crate::{Error, Result};

use super::Aggregate;

/// A party's share in the secret key switch protocol.
///
/// Each party uses the `SecretKeySwitchShare` to generate their share of the
/// new ciphertext and participate in the "Protocol 3: KeySwitch" protocol
/// detailed in [Multiparty BFV](https://eprint.iacr.org/2020/304.pdf) (p7). Use the [`Aggregate`] impl to combine the
/// shares into a [`Ciphertext`].
///
/// Note: this protocol assumes the output key is split into the same number of
/// parties as the input key, and is likely only useful for niche scenarios.
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
    /// 3. *Public input*: Input ciphertext to keyswitch
    // 4. *Public input*: TODO: variance of the ciphertext noise
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
        // Note: M-BFV implementation only supports ciphertext of length 2
        if ct.len() != 2 {
            return Err(Error::TooManyValues {
                actual: ct.len(),
                limit: 2,
            });
        }

        let par = sk_input_share.par.clone();
        let mut s_in = Zeroizing::new(Poly::try_convert_from(
            sk_input_share.coeffs.as_ref(),
            ct[0].ctx(),
            false,
            Representation::PowerBasis,
        )?);
        s_in.change_representation(Representation::Ntt);
        let mut s_out = Zeroizing::new(Poly::try_convert_from(
            sk_output_share.coeffs.as_ref(),
            ct[0].ctx(),
            false,
            Representation::PowerBasis,
        )?);
        s_out.change_representation(Representation::Ntt);

        // Sample error
        // TODO this should be exponential in ciphertext noise!
        let e = Zeroizing::new(Poly::small(
            ct[0].ctx(),
            Representation::Ntt,
            par.variance,
            rng,
        )?);

        // Create h_i share
        let mut h_share = s_in.as_ref() - s_out.as_ref();
        h_share.disallow_variable_time_computations();
        h_share *= &ct[1];
        h_share += e.as_ref();

        Ok(Self { par, ct, h_share })
    }
}

impl Aggregate<SecretKeySwitchShare> for Ciphertext {
    fn from_shares<T>(iter: T) -> Result<Self>
    where
        T: IntoIterator<Item = SecretKeySwitchShare>,
    {
        let mut shares = iter.into_iter();
        let share = shares.next().ok_or(Error::TooFewValues {
            actual: 0,
            minimum: 1,
        })?;
        let mut h = share.h_share;
        for sh in shares {
            h += &sh.h_share;
        }

        let c0 = &share.ct[0] + &h;
        let c1 = share.ct[1].clone();

        Ciphertext::new(vec![c0, c1], &share.par)
    }
}

/// A party's share in the decryption protocol.
///
/// Each party uses the `DecryptionShare` to generate their share of the
/// plaintext output. Note that this is a special case of the "Protocol 3:
/// KeySwitch" protocol detailed in [Multiparty BFV](https://eprint.iacr.org/2020/304.pdf) (p7), using an output key of zero. Use the
/// [`Aggregate`] impl to combine the shares into a [`Plaintext`].
pub struct DecryptionShare {
    pub(crate) sks_share: SecretKeySwitchShare,
}

impl DecryptionShare {
    /// Participate in a new Decryption protocol.
    ///
    /// 1. *Private input*: BFV input secret key share
    /// 3. *Public input*: Ciphertext to decrypt
    // 4. *Public input*: TODO: variance of the ciphertext noise
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

impl Aggregate<DecryptionShare> for Plaintext {
    fn from_shares<T>(iter: T) -> Result<Self>
    where
        T: IntoIterator<Item = DecryptionShare>,
    {
        let sks_shares = iter.into_iter().map(|s| s.sks_share);
        let ct = Ciphertext::from_shares(sks_shares)?;

        // Note: during SKS, c[1]*sk has already been added to c[0].
        let mut c = Zeroizing::new(ct[0].clone());
        c.disallow_variable_time_computations();
        c.change_representation(Representation::PowerBasis);

        // The true decryption part is done during SKS; all that is left is to scale
        let ctx_lvl = ct.par.context_level_at(ct.level)?;
        let d = Zeroizing::new(c.scale(&ctx_lvl.cipher_plain_context.scaler)?);

        let v: Vec<BigUint> = Vec::<BigUint>::from(d.as_ref())
            .into_iter()
            .map(|vi| vi + ct.par.plaintext_big())
            .collect_vec();

        let mut w = v[..ct.par.degree()].to_vec();
        let q_poly = d.as_ref().ctx().modulus();
        w.iter_mut().for_each(|wi| *wi %= q_poly);

        ct.par.plaintext.reduce_vec(&mut w);

        let mut poly =
            Poly::try_convert_from(w.as_slice(), ct[0].ctx(), false, Representation::PowerBasis)?;
        poly.change_representation(Representation::Ntt);

        let value = match ct.par.plaintext {
            crate::bfv::PlaintextModulus::Small(_) => PlaintextValues::Small(
                w.iter()
                    .map(|x| x.to_u64().unwrap())
                    .collect::<Vec<_>>()
                    .into_boxed_slice(),
            ),
            crate::bfv::PlaintextModulus::Large(_) => PlaintextValues::Large(w.into_boxed_slice()),
        };

        let pt = Plaintext {
            par: ct.par.clone(),
            value,
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

    use fhe_traits::{FheDecoder, FheEncoder, FheEncrypter};
    use rand::rng;

    use crate::{
        bfv::{BfvParameters, Encoding, Plaintext, PublicKey, SecretKey},
        mbfv::{
            Aggregate, AggregateIter, CommonRandomPoly, DecryptionShare, PublicKeyShare,
            SecretKeySwitchShare,
        },
    };

    const NUM_PARTIES: usize = 11;

    struct Party {
        sk_share: SecretKey,
        pk_share: PublicKeyShare,
    }

    #[test]
    fn encrypt_decrypt() {
        let mut rng = rng();
        for par in [
            BfvParameters::default_arc(1, 16),
            BfvParameters::default_arc(6, 32),
        ] {
            for level in 0..=par.max_level() {
                for _ in 0..20 {
                    let crp = CommonRandomPoly::new(&par, &mut rng).unwrap();

                    let mut parties: Vec<Party> = vec![];

                    // Parties collectively generate public key
                    for _ in 0..NUM_PARTIES {
                        let sk_share = SecretKey::random(&par, &mut rng);
                        let pk_share =
                            PublicKeyShare::new(&sk_share, crp.clone(), &mut rng).unwrap();
                        parties.push(Party { sk_share, pk_share })
                    }
                    let public_key: PublicKey = parties
                        .iter()
                        .map(|p| p.pk_share.clone())
                        .aggregate()
                        .unwrap();

                    // Use it to encrypt a random polynomial
                    let q = fhe_math::zq::Modulus::new(par.plaintext()).unwrap();
                    let pt1 = Plaintext::try_encode(
                        &q.random_vec(par.degree(), &mut rng),
                        Encoding::poly_at_level(level),
                        &par,
                    )
                    .unwrap();
                    let ct = Arc::new(public_key.try_encrypt(&pt1, &mut rng).unwrap());

                    // Parties perform a collective decryption
                    let decryption_shares = parties
                        .iter()
                        .map(|p| DecryptionShare::new(&p.sk_share, &ct, &mut rng));
                    let pt2 = Plaintext::from_shares(decryption_shares).unwrap();

                    assert_eq!(pt1, pt2);
                }
            }
        }
    }

    #[test]
    fn encrypt_keyswitch_decrypt() {
        let mut rng = rng();
        for par in [
            BfvParameters::default_arc(1, 16),
            BfvParameters::default_arc(6, 32),
        ] {
            for level in 0..=par.max_level() {
                for _ in 0..20 {
                    let crp = CommonRandomPoly::new(&par, &mut rng).unwrap();

                    // Parties collectively generate public key
                    let mut parties: Vec<Party> = vec![];
                    for _ in 0..NUM_PARTIES {
                        let sk_share = SecretKey::random(&par, &mut rng);
                        let pk_share =
                            PublicKeyShare::new(&sk_share, crp.clone(), &mut rng).unwrap();
                        parties.push(Party { sk_share, pk_share })
                    }

                    let public_key =
                        PublicKey::from_shares(parties.iter().map(|p| p.pk_share.clone())).unwrap();

                    // Use it to encrypt a random polynomial ct1
                    let q = fhe_math::zq::Modulus::new(par.plaintext()).unwrap();
                    let pt1 = Plaintext::try_encode(
                        &q.random_vec(par.degree(), &mut rng),
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
                    let ct2 = parties
                        .iter()
                        .zip(out_parties.iter())
                        .map(|(ip, op)| {
                            SecretKeySwitchShare::new(
                                &ip.sk_share,
                                &op.sk_share,
                                ct1.clone(),
                                &mut rng,
                            )
                        })
                        .aggregate()
                        .unwrap();
                    let ct2 = Arc::new(ct2);

                    // The second set of parties then does a collective decryption
                    let pt2 = out_parties
                        .iter()
                        .map(|p| DecryptionShare::new(&p.sk_share, &ct2, &mut rng))
                        .aggregate()
                        .unwrap();

                    assert_eq!(pt1, pt2);
                }
            }
        }
    }

    #[test]
    fn collective_keys_enable_homomorphic_addition() {
        let mut rng = rng();
        for par in [
            BfvParameters::default_arc(1, 16),
            BfvParameters::default_arc(6, 32),
        ] {
            for level in 0..=par.max_level() {
                for _ in 0..20 {
                    let crp = CommonRandomPoly::new(&par, &mut rng).unwrap();

                    let mut parties: Vec<Party> = vec![];

                    // Parties collectively generate public key
                    for _ in 0..NUM_PARTIES {
                        let sk_share = SecretKey::random(&par, &mut rng);
                        let pk_share =
                            PublicKeyShare::new(&sk_share, crp.clone(), &mut rng).unwrap();
                        parties.push(Party { sk_share, pk_share })
                    }
                    let public_key: PublicKey = parties
                        .iter()
                        .map(|p| p.pk_share.clone())
                        .aggregate()
                        .unwrap();

                    // Parties encrypt two plaintexts
                    let q = fhe_math::zq::Modulus::new(par.plaintext()).unwrap();
                    let a = q.random_vec(par.degree(), &mut rng);
                    let b = q.random_vec(par.degree(), &mut rng);
                    let mut expected = a.clone();
                    q.add_vec(&mut expected, &b);

                    let pt_a =
                        Plaintext::try_encode(&a, Encoding::poly_at_level(level), &par).unwrap();
                    let pt_b =
                        Plaintext::try_encode(&b, Encoding::poly_at_level(level), &par).unwrap();
                    let ct_a = public_key.try_encrypt(&pt_a, &mut rng).unwrap();
                    let ct_b = public_key.try_encrypt(&pt_b, &mut rng).unwrap();

                    // and add them together
                    let ct = Arc::new(&ct_a + &ct_b);

                    // Parties perform a collective decryption
                    let pt = parties
                        .iter()
                        .map(|p| DecryptionShare::new(&p.sk_share, &ct, &mut rng))
                        .aggregate()
                        .unwrap();

                    assert_eq!(
                        Vec::<u64>::try_decode(&pt, Encoding::poly_at_level(level)).unwrap(),
                        expected
                    );
                }
            }
        }
    }
}
