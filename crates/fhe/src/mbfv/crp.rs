use std::sync::Arc;

use crate::Result;
use crate::bfv::BfvParameters;
use fhe_math::rq::{Ntt, Poly};
use rand::{CryptoRng, RngCore};

/// A polynomial sampled from a random _common reference string_.
// TODO CRS->CRP implementation. For now just a random polynomial.
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct CommonRandomPoly {
    pub(crate) poly: Poly<Ntt>,
}

impl CommonRandomPoly {
    /// Generate a new random CRP.
    pub fn new<R: RngCore + CryptoRng>(par: &Arc<BfvParameters>, rng: &mut R) -> Result<Self> {
        Self::new_leveled(par, 0, rng)
    }

    /// Generate a new random CRP vector.
    ///
    /// The size of the vector is equal to the number of ciphertext moduli, as
    /// required for the relinearization key generation protocol.
    pub fn new_vec<R: RngCore + CryptoRng>(
        par: &Arc<BfvParameters>,
        rng: &mut R,
    ) -> Result<Vec<Self>> {
        (0..par.moduli().len())
            .map(|_| Self::new(par, rng))
            .collect()
    }

    /// Generate a new random leveled CRP.
    pub fn new_leveled<R: RngCore + CryptoRng>(
        par: &Arc<BfvParameters>,
        level: usize,
        rng: &mut R,
    ) -> Result<Self> {
        let ctx = par.context_at_level(level)?;
        let poly = Poly::<Ntt>::random(ctx, rng);
        Ok(Self { poly })
    }
}
