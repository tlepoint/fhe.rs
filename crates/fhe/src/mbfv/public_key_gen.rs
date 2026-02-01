use std::sync::Arc;

use crate::Error;
use crate::bfv::{BfvParameters, Ciphertext, PublicKey, SecretKey};
use crate::errors::Result;
use fhe_math::rq::{Poly, Representation, traits::TryConvertFrom};
use rand::{CryptoRng, RngCore};
use zeroize::Zeroizing;

use super::{Aggregate, CommonRandomPoly};

/// A party's share in public key generation protocol.
///
/// Each party uses the `PublicKeyShare` to generate their share of the public key and participate in the in the "Protocol 1: EncKeyGen", as detailed in [Multiparty BFV](https://eprint.iacr.org/2020/304.pdf) (p6). Use the [`Aggregate`] impl to combine the shares into a [`PublicKey`].
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct PublicKeyShare {
    pub(crate) par: Arc<BfvParameters>,
    pub(crate) crp: CommonRandomPoly,
    pub(crate) p0_share: Poly,
}

impl PublicKeyShare {
    /// Participate in a new EncKeyGen protocol.
    ///
    /// 1. *Private input*: BFV secret key share
    /// 2. *Public input*: common random polynomial
    //
    // Implementation note: This is largely the same approach taken by fhe.rs, a
    // symmetric encryption of zero, the difference being that the crp is used
    // instead of a random poly. Might be possible to just pass a valid seed to
    // each party and basically take the SecretKey::try_encrypt implementation,
    // but with the hardcoded seed.
    pub fn new<R: RngCore + CryptoRng>(
        sk_share: &SecretKey,
        crp: CommonRandomPoly,
        rng: &mut R,
    ) -> Result<Self> {
        let par = sk_share.par.clone();
        let ctx = par.context_at_level(0)?;

        // Convert secret key to usable polynomial
        let mut s = Zeroizing::new(Poly::try_convert_from(
            sk_share.coeffs.as_ref(),
            ctx,
            false,
            Representation::PowerBasis,
        )?);
        s.change_representation(Representation::Ntt);

        // Sample error
        let e = Zeroizing::new(Poly::small(ctx, Representation::Ntt, par.variance, rng)?);
        // Create p0_i share
        let mut p0_share = -crp.poly.clone();
        p0_share.disallow_variable_time_computations();
        p0_share.change_representation(Representation::Ntt);
        p0_share *= s.as_ref();
        p0_share += e.as_ref();
        unsafe { p0_share.allow_variable_time_computations() }

        Ok(Self { par, crp, p0_share })
    }
}

impl Aggregate<PublicKeyShare> for PublicKey {
    fn from_shares<T>(iter: T) -> Result<Self>
    where
        T: IntoIterator<Item = PublicKeyShare>,
    {
        let mut shares = iter.into_iter();
        let share = shares.next().ok_or(Error::TooFewValues {
            actual: 0,
            minimum: 1,
        })?;
        let mut p0 = share.p0_share;
        for sh in shares {
            p0 += &sh.p0_share;
        }

        Ok(PublicKey {
            c: Ciphertext::new(vec![p0, share.crp.poly], &share.par)?,
            par: share.par,
        })
    }
}

#[cfg(test)]
mod tests {
    use fhe_traits::{FheEncoder, FheEncrypter};
    use rand::rng;

    use crate::{
        bfv::{BfvParameters, Encoding, Plaintext, PublicKey, SecretKey},
        mbfv::{Aggregate as _, CommonRandomPoly},
    };

    use super::PublicKeyShare;

    const NUM_PARTIES: usize = 11;

    #[test]
    // This just makes sure the public key creation is successful, and arbitrary
    // encryptions complete without error. See a full encrypt->decrypt test in
    // `secret_key_switch`.
    fn protocol_creates_valid_pk() {
        let mut rng = rng();
        for par in [
            BfvParameters::default_arc(1, 16),
            BfvParameters::default_arc(6, 32),
        ] {
            for level in 0..=par.max_level() {
                for _ in 0..20 {
                    let crp = CommonRandomPoly::new(&par, &mut rng).unwrap();

                    let mut pk_shares: Vec<PublicKeyShare> = vec![];

                    // Parties collectively generate public key
                    for _ in 0..NUM_PARTIES {
                        let sk_share = SecretKey::random(&par, &mut rng);
                        let pk_share =
                            PublicKeyShare::new(&sk_share, crp.clone(), &mut rng).unwrap();
                        pk_shares.push(pk_share);
                    }
                    let public_key = PublicKey::from_shares(pk_shares).unwrap();

                    // Use it to encrypt a random polynomial
                    let pt = Plaintext::try_encode(
                        &fhe_math::zq::Modulus::new(par.plaintext())
                            .unwrap()
                            .random_vec(par.degree(), &mut rng),
                        Encoding::poly_at_level(level),
                        &par,
                    )
                    .unwrap();
                    let _ct = public_key.try_encrypt(&pt, &mut rng).unwrap();
                }
            }
        }
    }
}
