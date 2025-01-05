use std::sync::Arc;

use crate::bfv::BfvParameters;
use crate::Result;
use fhe_math::rq::Poly;
use rand::{CryptoRng, RngCore};

/// A polynomial sampled from a random _common reference string_.
// TODO CRS->CRP implementation. For now just a random polynomial.
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct CommonRandomPoly {
    pub(crate) poly: Poly,
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
        let ctx = par.ctx_at_level(level)?;
        let poly = Poly::random(ctx, fhe_math::rq::Representation::Ntt, rng);
        Ok(Self { poly })
    }
}
