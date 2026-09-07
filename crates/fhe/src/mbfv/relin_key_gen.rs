use std::marker::PhantomData;
use std::sync::Arc;

use crate::bfv::{KeySwitchingKey, Parameters, SecretKey, evaluation::RelinearizationKey};
use crate::errors::Result;
use fhe_math::rns::RnsContext;
use fhe_math::rq::{Ntt, NttShoup, Poly, PowerBasis};
use itertools::izip;
use rand::{CryptoRng, Rng as RngCore};
use zeroize::Zeroizing;

use super::round::{R1, R1Aggregated, R2, Round};
use super::{Aggregate, CommonRandomPoly};

/// A party's share in the relinearization key generation protocol.
/// Use the [`RelinKeyGenerator`] to create these shares.
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct RelinKeyShare<R: Round = R1> {
    pub(crate) par: Parameters,
    pub(crate) h0: Box<[Poly<Ntt>]>,
    pub(crate) h1: Box<[Poly<Ntt>]>,
    last_round: R::RelinDependency,
}

/// One party's state for a single execution of relinearization key generation.
/// Round transitions consume the state, preventing reuse of its ephemeral
/// secret. Round markers and parameter checks do not authenticate participants
/// or bind shares to an application session; callers must establish that
/// separately.
///
/// ```compile_fail
/// use fhe::mbfv::{RelinKeyGenerator, RelinKeyShare, round::R1Aggregated};
/// use std::sync::Arc;
/// fn skip_first_round(g: RelinKeyGenerator<'_, '_>, r1: &Arc<RelinKeyShare<R1Aggregated>>) {
///     g.round_2(r1, &mut rand::rng());
/// }
/// ```
///
/// Each party uses the `RelinKeyGenerator` to generate their shares and
/// participate in the "Protocol 2: RelinKeyGen" protocol detailed in
/// [Multiparty BFV](https://eprint.iacr.org/2020/304.pdf) p6. The shares need to be aggregated between
/// rounds:
///
/// ```rust
/// use std::sync::Arc;
/// use fhe::bfv::{ParametersBuilder, evaluation::RelinearizationKey, SecretKey};
/// use fhe::mbfv::{Aggregate, CommonRandomPoly, RelinKeyGenerator, RelinKeyShare, round::*};
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let parameters = ParametersBuilder::new()
///         .degree(4096)
///         .ciphertext_moduli(&[0xffffee001, 0xffffc4001, 0x1ffffe0001])
///         .plaintext_modulus(1_u64 << 10)
///         .build()?;
///
/// // Party perspective
/// let mut rng = rand::rng();
/// let sk_share = SecretKey::generate(&parameters, &mut rng);
/// let crp = CommonRandomPoly::new_vec(&parameters, &mut rng)?;
/// let rlk_generator = RelinKeyGenerator::new(&sk_share, &crp, &mut rng)?;
/// let (rlk_r1_share, rlk_generator) = rlk_generator.round_1(&mut rng)?;
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
pub struct RelinKeyGenerator<'a, 'b, R: Round = R1> {
    sk_share: &'a SecretKey,
    crp: &'b [CommonRandomPoly],
    u: Zeroizing<Poly<Ntt>>,
    round: PhantomData<R>,
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
            Err(crate::EvaluationKeyError::KeySwitchingNotSupported.into())
        } else if crp.len() != ctx.moduli().len() {
            Err(crate::MultipartyError::InvalidCommonRandomPolynomialCount {
                actual: crp.len(),
                expected: ctx.moduli().len(),
            }
            .into())
        } else {
            if crp.iter().any(|a| a.poly.ctx() != ctx) {
                return Err(crate::MultipartyError::IncompatibleShares.into());
            }
            let u = Zeroizing::new(Poly::<Ntt>::small(ctx, par.inner.variance, rng)?);
            Ok(Self {
                sk_share,
                crp,
                u,
                round: PhantomData,
            })
        }
    }

    /// Generate round one once, returning its share and the state needed for
    /// round two. Each generator owns a fresh ephemeral secret for one
    /// execution.
    pub fn round_1<R: RngCore + CryptoRng>(
        self,
        rng: &mut R,
    ) -> Result<(RelinKeyShare<R1>, RelinKeyGenerator<'a, 'b, R2>)> {
        let share = <RelinKeyShare<R1>>::new(self.sk_share, self.crp, &self.u, rng)?;
        Ok((
            share,
            RelinKeyGenerator {
                sk_share: self.sk_share,
                crp: self.crp,
                u: self.u,
                round: PhantomData,
            },
        ))
    }
}

impl RelinKeyGenerator<'_, '_, R2> {
    /// Generate round two once and clear the ephemeral secret when this
    /// consumed state is dropped, including on an error. Clone the
    /// resulting share for retransmission; create a new generator for a new
    /// protocol execution.
    pub fn round_2<R: RngCore + CryptoRng>(
        self,
        r1: &Arc<RelinKeyShare<R1Aggregated>>,
        rng: &mut R,
    ) -> Result<RelinKeyShare<R2>> {
        <RelinKeyShare<R2>>::new(self.sk_share, &self.u, r1, rng)
    }
}

impl<R: Round> RelinKeyShare<R> {
    fn validate_for(&self, par: &Parameters) -> Result<()> {
        let ctx = par.context_at_level(0)?;
        if !Parameters::compatible(&self.par, par)
            || self.h0.len() != ctx.moduli().len()
            || self.h1.len() != ctx.moduli().len()
            || self
                .h0
                .iter()
                .chain(self.h1.iter())
                .any(|p| p.ctx() != ctx || !p.is_canonical())
        {
            return Err(crate::MultipartyError::IncompatibleShares.into());
        }
        Ok(())
    }
}

impl RelinKeyShare<R1> {
    fn new<R: RngCore + CryptoRng>(
        sk_share: &SecretKey,
        crp: &[CommonRandomPoly],
        u: &Zeroizing<Poly<Ntt>>,
        rng: &mut R,
    ) -> Result<Self> {
        let par = sk_share.par.clone();

        let expected_crp_count = par.context_at_level(0)?.moduli().len();
        if crp.len() != expected_crp_count {
            Err(crate::MultipartyError::InvalidCommonRandomPolynomialCount {
                actual: crp.len(),
                expected: expected_crp_count,
            }
            .into())
        } else {
            let h0 = Self::generate_h0(sk_share, crp, u, rng)?;
            let h1 = Self::generate_h1(sk_share, crp, rng)?;
            Ok(Self {
                par,
                h0,
                h1,
                last_round: (),
            })
        }
    }

    fn generate_h0<R: RngCore + CryptoRng>(
        sk_share: &SecretKey,
        crp: &[CommonRandomPoly],
        u: &Zeroizing<Poly<Ntt>>,
        rng: &mut R,
    ) -> Result<Box<[Poly<Ntt>]>> {
        let par = sk_share.par.clone();
        let ctx = par.context_at_level(0)?;

        let s = Zeroizing::new(
            Poly::<PowerBasis>::from_signed_coefficients(sk_share.coeffs.as_ref(), ctx)?.into_ntt(),
        );
        let rns = RnsContext::new(&sk_share.par.inner.moduli[..crp.len()])?;
        let h0 = crp
            .iter()
            .enumerate()
            .map(|(i, a)| {
                let w = rns.get_garner(i).unwrap();
                let w_s = Zeroizing::new(w * s.as_ref());

                let e = Zeroizing::new(Poly::<Ntt>::small(ctx, par.inner.variance, rng)?);

                let mut h = -a.poly.clone();
                h.disallow_variable_time_computations();
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
    ) -> Result<Box<[Poly<Ntt>]>> {
        let par = sk_share.par.clone();
        let ctx = par.context_at_level(0)?;
        let s = Zeroizing::new(
            Poly::<PowerBasis>::from_signed_coefficients(sk_share.coeffs.as_ref(), ctx)?.into_ntt(),
        );

        let h1 = crp
            .iter()
            .map(|a| {
                let mut h = a.poly.clone();
                h.disallow_variable_time_computations();
                let e = Zeroizing::new(Poly::<Ntt>::small(ctx, par.inner.variance, rng)?);
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
        let share = shares.next().ok_or(crate::MultipartyError::NoShares)?;
        share.validate_for(&share.par)?;
        let mut h0 = share.h0;
        let mut h1 = share.h1;
        for sh in shares {
            sh.validate_for(&share.par)?;
            izip!(h0.iter_mut(), sh.h0.iter()).for_each(|(h0i, sh_h0i)| *h0i += sh_h0i);
            izip!(h1.iter_mut(), sh.h1.iter()).for_each(|(h1i, sh_h1i)| *h1i += sh_h1i);
        }

        Ok(RelinKeyShare {
            par: share.par,
            h0,
            h1,
            last_round: (),
        })
    }
}

impl RelinKeyShare<R2> {
    fn new<R: RngCore + CryptoRng>(
        sk_share: &SecretKey,
        u: &Zeroizing<Poly<Ntt>>,
        r1: &Arc<RelinKeyShare<R1Aggregated>>,
        rng: &mut R,
    ) -> Result<Self> {
        let par = sk_share.par.clone();
        r1.validate_for(&par)?;
        let h0 = Self::generate_h0(sk_share, &r1.h0, rng)?;
        let h1 = Self::generate_h1(sk_share, u, &r1.h1, rng)?;
        Ok(Self {
            par,
            h0,
            h1,
            last_round: Arc::clone(r1),
        })
    }

    fn generate_h0<R: RngCore + CryptoRng>(
        sk_share: &SecretKey,
        r1_h0: &[Poly<Ntt>],
        rng: &mut R,
    ) -> Result<Box<[Poly<Ntt>]>> {
        let par = sk_share.par.clone();
        let ctx = par.context_at_level(0)?;

        let s = Zeroizing::new(
            Poly::<PowerBasis>::from_signed_coefficients(sk_share.coeffs.as_ref(), ctx)?.into_ntt(),
        );
        let h0 = r1_h0
            .iter()
            .map(|h| {
                let e = Zeroizing::new(Poly::<Ntt>::small(ctx, par.inner.variance, rng)?);

                let mut h_prime = h.clone();
                h_prime.disallow_variable_time_computations();
                h_prime *= s.as_ref();

                h_prime += e.as_ref();
                Ok(h_prime)
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(h0.into_boxed_slice())
    }

    fn generate_h1<R: RngCore + CryptoRng>(
        sk_share: &SecretKey,
        u: &Zeroizing<Poly<Ntt>>,
        r1_h1: &[Poly<Ntt>],
        rng: &mut R,
    ) -> Result<Box<[Poly<Ntt>]>> {
        let par = sk_share.par.clone();
        let ctx = par.context_at_level(0)?;
        let s = Zeroizing::new(
            Poly::<PowerBasis>::from_signed_coefficients(sk_share.coeffs.as_ref(), ctx)?.into_ntt(),
        );

        let u_s = Zeroizing::new(u.as_ref() - s.as_ref());

        let h1 = r1_h1
            .iter()
            .map(|h| {
                let mut h_prime = h.clone();
                h_prime.disallow_variable_time_computations();
                let e = Zeroizing::new(Poly::<Ntt>::small(ctx, par.inner.variance, rng)?);
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
        let share = shares.next().ok_or(crate::MultipartyError::NoShares)?;
        share.validate_for(&share.par)?;
        let par = share.par.clone();
        let ctx = par.context_at_level(0)?.clone();
        let r1 = share.last_round;
        r1.validate_for(&par)?;

        let mut h0 = share.h0;
        let mut h1 = share.h1;
        for sh in shares {
            sh.validate_for(&par)?;
            if !Arc::ptr_eq(&r1, &sh.last_round) && r1 != sh.last_round {
                return Err(crate::MultipartyError::RoundOneAggregationMismatch.into());
            }
            izip!(h0.iter_mut(), h1.iter_mut(), sh.h0.iter(), sh.h1.iter()).for_each(
                |(h0, h1, h0i, h1i)| {
                    *h0 += h0i;
                    *h1 += h1i;
                },
            );
        }

        let mut c0 = Vec::from(h0);
        izip!(c0.iter_mut(), h1.iter()).for_each(|(c0, h1)| *c0 += h1);
        let c0 = c0
            .into_iter()
            .map(Poly::<Ntt>::into_ntt_shoup)
            .collect::<Vec<Poly<NttShoup>>>()
            .into_boxed_slice();

        let c1 = r1
            .h1
            .iter()
            .cloned()
            .map(Poly::<Ntt>::into_ntt_shoup)
            .collect::<Vec<Poly<NttShoup>>>()
            .into_boxed_slice();

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

    use rand::rng;

    use crate::{
        bfv::{
            Encoding, Parameters, Plaintext, PublicKey, SecretKey,
            evaluation::{MultiplicationPlan, RelinearizationKey},
        },
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
            Parameters::test_parameters(3, 16),
            Parameters::test_parameters(6, 32),
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
                    let sk_share = SecretKey::generate(&par, &mut rng);
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

                let (r1_shares, round_two): (Vec<_>, Vec<_>) = party_rlks
                    .into_iter()
                    .map(|g| g.round_1(&mut rng).unwrap())
                    .unzip();
                let rlk_r1 = Arc::new(
                    super::RelinKeyShare::<super::R1Aggregated>::from_shares(r1_shares).unwrap(),
                );
                let rlk: RelinearizationKey = round_two
                    .into_iter()
                    .map(|g| g.round_2(&rlk_r1, &mut rng))
                    .aggregate()
                    .unwrap();

                // Create a couple random encrypted polynomials
                let v1 = fhe_math::zq::Modulus::new(par.plaintext_modulus_u64().unwrap())
                    .unwrap()
                    .random_vec(par.degree(), &mut rng);
                let v2 = fhe_math::zq::Modulus::new(par.plaintext_modulus_u64().unwrap())
                    .unwrap()
                    .random_vec(par.degree(), &mut rng);
                let pt1 = Plaintext::encode_at_level(&par, &v1, Encoding::Simd, level).unwrap();
                let pt2 = Plaintext::encode_at_level(&par, &v2, Encoding::Simd, level).unwrap();
                let ct1 = public_key.encrypt(&pt1, &mut rng).unwrap();
                let ct2 = public_key.encrypt(&pt2, &mut rng).unwrap();

                // Multiply them
                let multiplicator = MultiplicationPlan::builder(&rlk.ksk.par)
                    .relinearization(&rlk)
                    .modulus_switching(par.max_level() > 0)
                    .build()
                    .unwrap();
                let ct = Arc::new(multiplicator.multiply(&ct1, &ct2).unwrap());
                assert_eq!(ct.len(), 2);

                // Parties perform a collective decryption
                let pt: Plaintext = party_sks
                    .iter()
                    .map(|s| DecryptionShare::new(s, &ct, &mut rng))
                    .aggregate()
                    .unwrap();

                let mut expected = v1.clone();
                fhe_math::zq::Modulus::new(par.plaintext_modulus_u64().unwrap())
                    .unwrap()
                    .mul_vec(&mut expected, &v2);
                assert_eq!(pt.decode(Encoding::Simd).unwrap(), expected);
            }
        }
    }
}

#[cfg(test)]
mod state_tests {
    use super::*;
    use rand::{RngExt, SeedableRng};
    use rand_chacha::ChaCha8Rng;

    #[test]
    fn rounds_reject_foreign_parameters_shapes_and_aggregations() -> Result<()> {
        let par = Parameters::test_parameters(2, 16);
        let foreign = Parameters::builder()
            .degree(8)
            .plaintext_modulus(17_u64)
            .ciphertext_modulus_bits([50, 50])
            .build()?;
        let mut rng = ChaCha8Rng::seed_from_u64(123);
        let sk = SecretKey::generate(&par, &mut rng);
        let foreign_sk = SecretKey::generate(&foreign, &mut rng);
        let crp = CommonRandomPoly::new_vec(&par, &mut rng)?;
        let foreign_crp = CommonRandomPoly::new_vec(&foreign, &mut rng)?;
        assert!(RelinKeyGenerator::new(&sk, &foreign_crp, &mut rng).is_err());
        let (one, next) = RelinKeyGenerator::new(&sk, &crp, &mut rng)?.round_1(&mut rng)?;
        let (other, _) =
            RelinKeyGenerator::new(&foreign_sk, &foreign_crp, &mut rng)?.round_1(&mut rng)?;
        assert!(RelinKeyShare::<R1Aggregated>::from_shares([one.clone(), other.clone()]).is_err());
        let mut short = one.clone();
        short.h0 = Box::new([]);
        assert!(RelinKeyShare::<R1Aggregated>::from_shares([one.clone(), short]).is_err());
        let wrong_r1 = Arc::new(RelinKeyShare::<R1Aggregated>::from_shares([other])?);
        let mut saved_rng = rng.clone();
        assert!(next.round_2(&wrong_r1, &mut rng).is_err());
        assert_eq!(rng.random::<u64>(), saved_rng.random::<u64>());

        let (a, next_a) = RelinKeyGenerator::new(&sk, &crp, &mut rng)?.round_1(&mut rng)?;
        let (b, next_b) = RelinKeyGenerator::new(&sk, &crp, &mut rng)?.round_1(&mut rng)?;
        let r1_a = Arc::new(RelinKeyShare::<R1Aggregated>::from_shares([a])?);
        let r1_b = Arc::new(RelinKeyShare::<R1Aggregated>::from_shares([b])?);
        let two_a = next_a.round_2(&r1_a, &mut rng)?;
        let two_b = next_b.round_2(&r1_b, &mut rng)?;
        assert!(matches!(
            RelinearizationKey::from_shares([two_a.clone(), two_b]),
            Err(crate::Error::Multiparty(
                crate::MultipartyError::RoundOneAggregationMismatch
            ))
        ));
        // An equivalent independently allocated aggregate is accepted.
        let mut two_equivalent = two_a.clone();
        two_equivalent.last_round = Arc::new((*r1_a).clone());
        RelinearizationKey::from_shares([two_a, two_equivalent])?;
        Ok(())
    }
}
