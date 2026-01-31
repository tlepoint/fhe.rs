use std::marker::PhantomData;
use std::sync::Arc;

use crate::Error;
use crate::bfv::{BfvParameters, KeySwitchingKey, RelinearizationKey, SecretKey};
use crate::errors::Result;
use fhe_math::rns::RnsContext;
use fhe_math::rq::{Poly, Representation, traits::TryConvertFrom};
use itertools::izip;
use rand::{CryptoRng, RngCore};
use zeroize::Zeroizing;

use super::round::{R1, R1Aggregated, R2, Round};
use super::{Aggregate, CommonRandomPoly};

/// A party's share in the relinearization key generation protocol.
/// Use the [`RelinKeyGenerator`] to create these shares.
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct RelinKeyShare<R: Round = R1> {
    pub(crate) par: Arc<BfvParameters>,
    pub(crate) h0: Box<[Poly]>,
    pub(crate) h1: Box<[Poly]>,
    last_round: Option<Arc<RelinKeyShare<R1Aggregated>>>,
    _phantom_data: PhantomData<R>,
}

/// A builder for creating relinearization key generation shares per party.
///
/// Each party uses the `RelinKeyGenerator` to generate their shares and
/// participate in the "Protocol 2: RelinKeyGen" protocol detailed in
/// [Multiparty BFV](https://eprint.iacr.org/2020/304.pdf) (p6). The shares need to be aggregated between
/// rounds:
///
/// ```rust
/// use std::sync::Arc;
/// use fhe::bfv::{BfvParametersBuilder, RelinearizationKey, SecretKey};
/// use fhe::mbfv::{Aggregate, CommonRandomPoly, RelinKeyGenerator, RelinKeyShare, round::*};
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let parameters = BfvParametersBuilder::new()
///         .set_degree(4096)
///         .set_moduli(&[0xffffee001, 0xffffc4001, 0x1ffffe0001])
///         .set_plaintext_modulus(1 << 10)
///         .build_arc()?;
///
/// // Party perspective
/// let mut rng = rand::rng();
/// let sk_share = SecretKey::random(&parameters, &mut rng);
/// let crp = CommonRandomPoly::new_vec(&parameters, &mut rng)?;
/// let rlk_generator = RelinKeyGenerator::new(&sk_share, &crp, &mut rng)?;
/// let rlk_r1_share = rlk_generator.round_1(&mut rng)?;
///
/// // Aggregator perspective
/// let r1_shares = vec![rlk_r1_share]; // all party shares go here
/// let rlk_r1_aggregated = RelinKeyShare::<R1Aggregated>::from_shares(r1_shares)?;
///
/// // Party perspective
/// let rlk_r2_share = rlk_generator.round_2(&Arc::new(rlk_r1_aggregated), &mut rng)?;
///
/// // Aggregator perspective
/// let r2_shares = vec![rlk_r2_share]; // all party shares go here
/// let rlk = RelinearizationKey::from_shares(r2_shares)?;
/// # Ok(())
/// # }
/// ```
pub struct RelinKeyGenerator<'a, 'b> {
    sk_share: &'a SecretKey,
    crp: &'b [CommonRandomPoly],
    u: Zeroizing<Poly>,
}

impl<'a, 'b> RelinKeyGenerator<'a, 'b> {
    /// Create a new relin key generator for a given party.
    ///
    /// 1. *Private input*: BFV secret key share
    /// 2. *Public input*: common random polynomial vector
    pub fn new<R: RngCore + CryptoRng>(
        sk_share: &'a SecretKey,
        crp: &'b [CommonRandomPoly],
        rng: &mut R,
    ) -> Result<Self> {
        let par = sk_share.par.clone();
        let ctx = par.context_at_level(0)?;
        if ctx.moduli().len() == 1 {
            Err(Error::DefaultError(
                "These parameters do not support key switching".to_string(),
            ))
        } else if crp.len() != ctx.moduli().len() {
            Err(Error::DefaultError(
                "The size of the CRP polynomial vector must equal the number of ciphertext moduli."
                    .to_string(),
            ))
        } else {
            let u = Zeroizing::new(Poly::small(ctx, Representation::Ntt, par.variance, rng)?);
            Ok(Self { sk_share, crp, u })
        }
    }

    /// Generate share for round 1
    pub fn round_1<R: RngCore + CryptoRng>(&self, rng: &mut R) -> Result<RelinKeyShare<R1>> {
        <RelinKeyShare<R1>>::new(self.sk_share, self.crp, &self.u, rng)
    }

    /// Generate share for round 2
    pub fn round_2<R: RngCore + CryptoRng>(
        &self,
        r1: &Arc<RelinKeyShare<R1Aggregated>>,
        rng: &mut R,
    ) -> Result<RelinKeyShare<R2>> {
        <RelinKeyShare<R2>>::new(self.sk_share, &self.u, r1, rng)
    }
}

impl RelinKeyShare<R1> {
    fn new<R: RngCore + CryptoRng>(
        sk_share: &SecretKey,
        crp: &[CommonRandomPoly],
        u: &Zeroizing<Poly>,
        rng: &mut R,
    ) -> Result<Self> {
        let par = sk_share.par.clone();

        if crp.len() != par.context_at_level(0)?.moduli().len() {
            Err(Error::DefaultError(
                "The size of the CRP polynomial vector must equal the number of ciphertext moduli."
                    .to_string(),
            ))
        } else {
            let h0 = Self::generate_h0(sk_share, crp, u, rng)?;
            let h1 = Self::generate_h1(sk_share, crp, rng)?;
            Ok(Self {
                par,
                h0,
                h1,
                last_round: None,
                _phantom_data: PhantomData,
            })
        }
    }

    fn generate_h0<R: RngCore + CryptoRng>(
        sk_share: &SecretKey,
        crp: &[CommonRandomPoly],
        u: &Zeroizing<Poly>,
        rng: &mut R,
    ) -> Result<Box<[Poly]>> {
        let par = sk_share.par.clone();
        let ctx = par.context_at_level(0)?;

        let s = Zeroizing::new(Poly::try_convert_from(
            sk_share.coeffs.as_ref(),
            ctx,
            false,
            Representation::PowerBasis,
        )?);
        let rns = RnsContext::new(&sk_share.par.moduli[..crp.len()])?;
        let h0 = crp
            .iter()
            .enumerate()
            .map(|(i, a)| {
                let w = rns.get_garner(i).unwrap();
                let mut w_s = Zeroizing::new(w * s.as_ref());
                w_s.change_representation(Representation::Ntt);

                let e = Zeroizing::new(Poly::small(ctx, Representation::Ntt, par.variance, rng)?);

                let mut h = -a.poly.clone();
                h.disallow_variable_time_computations();
                h.change_representation(Representation::Ntt);
                h *= u.as_ref();
                h += w_s.as_ref();
                h += e.as_ref();
                Ok(h)
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(h0.into_boxed_slice())
    }

    fn generate_h1<R: RngCore + CryptoRng>(
        sk_share: &SecretKey,
        crp: &[CommonRandomPoly],
        rng: &mut R,
    ) -> Result<Box<[Poly]>> {
        let par = sk_share.par.clone();
        let ctx = par.context_at_level(0)?;
        let mut s = Zeroizing::new(Poly::try_convert_from(
            sk_share.coeffs.as_ref(),
            ctx,
            false,
            Representation::PowerBasis,
        )?);
        s.change_representation(Representation::Ntt);

        let h1 = crp
            .iter()
            .map(|a| {
                let mut h = a.poly.clone();
                h.disallow_variable_time_computations();
                h.change_representation(Representation::Ntt);
                let e = Zeroizing::new(Poly::small(ctx, Representation::Ntt, par.variance, rng)?);
                h *= s.as_ref();
                h += e.as_ref();
                Ok(h)
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(h1.into_boxed_slice())
    }
}

impl Aggregate<RelinKeyShare<R1>> for RelinKeyShare<R1Aggregated> {
    fn from_shares<T>(iter: T) -> Result<Self>
    where
        T: IntoIterator<Item = RelinKeyShare<R1>>,
    {
        let mut shares = iter.into_iter();
        let share = shares.next().ok_or(Error::TooFewValues {
            actual: 0,
            minimum: 1,
        })?;
        let mut h0 = share.h0;
        let mut h1 = share.h1;
        for sh in shares {
            izip!(h0.iter_mut(), sh.h0.iter()).for_each(|(h0i, sh_h0i)| *h0i += sh_h0i);
            izip!(h1.iter_mut(), sh.h1.iter()).for_each(|(h1i, sh_h1i)| *h1i += sh_h1i);
        }

        Ok(RelinKeyShare {
            par: share.par,
            h0,
            h1,
            last_round: None,
            _phantom_data: PhantomData,
        })
    }
}

impl RelinKeyShare<R2> {
    fn new<R: RngCore + CryptoRng>(
        sk_share: &SecretKey,
        u: &Zeroizing<Poly>,
        r1: &Arc<RelinKeyShare<R1Aggregated>>,
        rng: &mut R,
    ) -> Result<Self> {
        let par = sk_share.par.clone();
        let h0 = Self::generate_h0(sk_share, &r1.h0, rng)?;
        let h1 = Self::generate_h1(sk_share, u, &r1.h1, rng)?;
        Ok(Self {
            par,
            h0,
            h1,
            last_round: Some(Arc::clone(r1)),
            _phantom_data: PhantomData,
        })
    }

    fn generate_h0<R: RngCore + CryptoRng>(
        sk_share: &SecretKey,
        r1_h0: &[Poly],
        rng: &mut R,
    ) -> Result<Box<[Poly]>> {
        let par = sk_share.par.clone();
        let ctx = par.context_at_level(0)?;

        let mut s = Zeroizing::new(Poly::try_convert_from(
            sk_share.coeffs.as_ref(),
            ctx,
            false,
            Representation::PowerBasis,
        )?);
        s.change_representation(Representation::Ntt);
        let h0 = r1_h0
            .iter()
            .map(|h| {
                let e = Zeroizing::new(Poly::small(ctx, Representation::Ntt, par.variance, rng)?);

                let mut h_prime = h.clone();
                h_prime.disallow_variable_time_computations();
                h_prime.change_representation(Representation::Ntt);
                h_prime *= s.as_ref();

                h_prime += e.as_ref();
                Ok(h_prime)
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(h0.into_boxed_slice())
    }

    fn generate_h1<R: RngCore + CryptoRng>(
        sk_share: &SecretKey,
        u: &Zeroizing<Poly>,
        r1_h1: &[Poly],
        rng: &mut R,
    ) -> Result<Box<[Poly]>> {
        let par = sk_share.par.clone();
        let ctx = par.context_at_level(0)?;
        let mut s = Zeroizing::new(Poly::try_convert_from(
            sk_share.coeffs.as_ref(),
            ctx,
            false,
            Representation::PowerBasis,
        )?);
        s.change_representation(Representation::Ntt);

        let u_s = Zeroizing::new(u.as_ref() - s.as_ref());

        let h1 = r1_h1
            .iter()
            .map(|h| {
                let mut h_prime = h.clone();
                h_prime.disallow_variable_time_computations();
                h_prime.change_representation(Representation::Ntt);
                let e = Zeroizing::new(Poly::small(ctx, Representation::Ntt, par.variance, rng)?);
                h_prime *= u_s.as_ref();
                h_prime += e.as_ref();
                Ok(h_prime)
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(h1.into_boxed_slice())
    }
}

impl Aggregate<RelinKeyShare<R2>> for RelinearizationKey {
    fn from_shares<T>(iter: T) -> Result<Self>
    where
        T: IntoIterator<Item = RelinKeyShare<R2>>,
    {
        let mut shares = iter.into_iter();
        let share = shares.next().ok_or(Error::TooFewValues {
            actual: 0,
            minimum: 1,
        })?;
        let par = share.par.clone();
        let ctx = par.context_at_level(0)?.clone();
        let r1 = share.last_round.ok_or(Error::DefaultError(
            "Shares from round 2 should include a copy for the round 1 aggregation.".to_string(),
        ))?;

        let mut h0 = share.h0;
        let mut h1 = share.h1;
        for sh in shares {
            izip!(h0.iter_mut(), h1.iter_mut(), sh.h0.iter(), sh.h1.iter()).for_each(
                |(h0, h1, h0i, h1i)| {
                    *h0 += h0i;
                    *h1 += h1i;
                },
            );
        }

        let mut c0 = h0;
        izip!(c0.iter_mut(), h1.iter()).for_each(|(c0, h1)| {
            *c0 += h1;
            c0.change_representation(Representation::NttShoup);
        });

        let mut c1 = r1.h1.clone();
        c1.iter_mut().for_each(|c1| {
            c1.change_representation(Representation::NttShoup);
        });

        let ksk = KeySwitchingKey {
            par,
            c0,
            c1,
            seed: None,
            ciphertext_level: 0,
            ctx_ciphertext: ctx.clone(),
            ksk_level: 0,
            ctx_ksk: ctx.clone(),
            log_base: 0,
        };
        Ok(RelinearizationKey { ksk })
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use fhe_traits::{FheDecoder, FheEncoder, FheEncrypter};
    use rand::rng;

    use crate::{
        bfv::{BfvParameters, Encoding, Multiplicator, Plaintext, PublicKey, SecretKey},
        mbfv::{
            Aggregate as _, AggregateIter, CommonRandomPoly, DecryptionShare, PublicKeyShare,
            RelinKeyGenerator,
        },
    };

    const NUM_PARTIES: usize = 5;

    #[test]
    fn relinearization_works() {
        let mut rng = rng();
        for par in [
            BfvParameters::default_arc(3, 16),
            BfvParameters::default_arc(6, 32),
        ] {
            // Just support level 0 for now.
            let level = 0;
            for _ in 0..10 {
                let crp = CommonRandomPoly::new_vec(&par, &mut rng).unwrap();

                let mut party_sks: Vec<SecretKey> = vec![];
                let mut party_pks: Vec<PublicKeyShare> = vec![];
                let mut party_rlks: Vec<RelinKeyGenerator> = vec![];

                // Parties undergo round 1
                for _ in 0..NUM_PARTIES {
                    let sk_share = SecretKey::random(&par, &mut rng);
                    party_sks.push(sk_share);
                }
                let crp_pk = CommonRandomPoly::new(&par, &mut rng).unwrap();
                (0..NUM_PARTIES).for_each(|i| {
                    let pk_share =
                        PublicKeyShare::new(&party_sks[i], crp_pk.clone(), &mut rng).unwrap();
                    let rlk_generator =
                        RelinKeyGenerator::new(&party_sks[i], &crp, &mut rng).unwrap();
                    party_pks.push(pk_share);
                    party_rlks.push(rlk_generator);
                });

                // Aggregate pk shares into public key
                let public_key = PublicKey::from_shares(party_pks).unwrap();

                // Aggregate rlk r1 shares
                let rlk_r1 = Arc::new(
                    party_rlks
                        .iter()
                        .map(|g| g.round_1(&mut rng))
                        .aggregate()
                        .unwrap(),
                );

                // Aggregate rlk r2 shares into relin key
                let rlk = party_rlks
                    .iter()
                    .map(|g| g.round_2(&rlk_r1, &mut rng))
                    .aggregate()
                    .unwrap();

                // Create a couple random encrypted polynomials
                let v1 = fhe_math::zq::Modulus::new(par.plaintext())
                    .unwrap()
                    .random_vec(par.degree(), &mut rng);
                let v2 = fhe_math::zq::Modulus::new(par.plaintext())
                    .unwrap()
                    .random_vec(par.degree(), &mut rng);
                let pt1 = Plaintext::try_encode(&v1, Encoding::simd_at_level(level), &par).unwrap();
                let pt2 = Plaintext::try_encode(&v2, Encoding::simd_at_level(level), &par).unwrap();
                let ct1 = public_key.try_encrypt(&pt1, &mut rng).unwrap();
                let ct2 = public_key.try_encrypt(&pt2, &mut rng).unwrap();

                // Multiply them
                let mut multiplicator = Multiplicator::default(&rlk).unwrap();
                if par.moduli().len() > 1 {
                    multiplicator.enable_mod_switching().unwrap();
                }
                let ct = Arc::new(multiplicator.multiply(&ct1, &ct2).unwrap());
                assert_eq!(ct.len(), 2);

                // Parties perform a collective decryption
                let pt = party_sks
                    .iter()
                    .map(|s| DecryptionShare::new(s, &ct, &mut rng))
                    .aggregate()
                    .unwrap();

                let mut expected = v1.clone();
                fhe_math::zq::Modulus::new(par.plaintext())
                    .unwrap()
                    .mul_vec(&mut expected, &v2);
                assert_eq!(
                    Vec::<u64>::try_decode(&pt, Encoding::simd_at_level(pt.level)).unwrap(),
                    expected
                );
            }
        }
    }
}
