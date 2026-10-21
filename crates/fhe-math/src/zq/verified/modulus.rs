//! Validate a modulus and construct the arithmetic metadata used by the
//! kernels.

use super::MODULUS_LIMIT;
#[cfg(verus_keep_ghost)]
use super::lazy_reduce::{base, leading_zero_bounds, split_u128, valid_barrett, valid_opt};
use super::supports_opt;
#[cfg(verus_keep_ghost)]
use vstd::arithmetic::div_mod::{
    lemma_fundamental_div_mod, lemma_fundamental_div_mod_converse_div, lemma_mod_bound,
};
use vstd::prelude::*;
#[cfg(verus_keep_ghost)]
use vstd::std_specs::bits::u64_leading_zeros;

verus! {

/// The arithmetic portion of a `Modulus`, independent of RNG and CPU dispatch.
pub(crate) struct ModulusParameters {
    pub p: u64,
    pub barrett_hi: u64,
    pub barrett_lo: u64,
    pub leading_zeros: u32,
    pub supports_opt: bool,
}

impl ModulusParameters {
    /// Exactly the metadata conditions consumed by the reduction kernels.
    pub(crate) open spec fn valid(&self) -> bool {
        &&& valid_barrett(self.p, self.barrett_hi, self.barrett_lo)
        &&& self.leading_zeros == u64_leading_zeros(self.p)
        &&& 2 <= self.leading_zeros <= 62
        &&& self.supports_opt == valid_opt(self.p, self.leading_zeros)
    }
}

/// Incrementing the dividend increments the quotient only at a remainder of p-1.
#[verifier::nonlinear]
proof fn reciprocal_bounds(p: int, q: int, r: int)
    requires
        2 <= p,
        0 <= q,
        0 <= r < p,
        u128::MAX as int == q * p + r,
    ensures
        q < u128::MAX as int,
        ({ let mu = q + if r == p - 1 { 1int } else { 0int };
            mu * p <= base() * base() < (mu + 1) * p }),
{
}

/// Compute floor(2^128 / p) without representing 2^128 in a machine integer.
fn barrett_reciprocal(p: u64) -> (mu: u128)
    requires 2 <= p < MODULUS_LIMIT,
    ensures
        mu as int == (base() * base()) / p as int,
        mu * p <= base() * base() < (mu + 1) * p,
{
    let divisor = p as u128;
    let q = u128::MAX / divisor;
    let r = u128::MAX % divisor;
    proof {
        lemma_fundamental_div_mod(u128::MAX as int, p as int);
        lemma_mod_bound(u128::MAX as int, p as int);
        assert(u128::MAX as int == q as int * p as int + r as int) by (nonlinear_arith)
            requires u128::MAX as int == p as int * q as int + r as int;
        reciprocal_bounds(p as int, q as int, r as int);
    }
    let mu = q + if r == divisor - 1 { 1 } else { 0 };
    proof {
        let remainder = base() * base() - mu as int * p as int;
        assert(0 <= remainder < p) by (nonlinear_arith)
            requires mu * p <= base() * base() < (mu + 1) * p,
                remainder == base() * base() - mu as int * p as int;
        lemma_fundamental_div_mod_converse_div(
            base() * base(), p as int, mu as int, remainder,
        );
    }
    mu
}

/// Accept exactly the supported moduli and construct valid reduction metadata.
pub(crate) fn new_modulus(p: u64) -> (result: Option<ModulusParameters>)
    ensures
        result.is_some() <==> 2 <= p < MODULUS_LIMIT,
        result matches Some(parameters) ==> parameters.p == p && parameters.valid(),
{
    if !(2..MODULUS_LIMIT).contains(&p) {
        None
    } else {
        let mu = barrett_reciprocal(p);
        let barrett_hi = (mu >> 64) as u64;
        let barrett_lo = mu as u64;
        let leading_zeros = p.leading_zeros();
        let supports_opt = supports_opt(p);
        proof {
            split_u128(mu);
            leading_zero_bounds(p);
        }
        Some(ModulusParameters { p, barrett_hi, barrett_lo, leading_zeros, supports_opt })
    }
}

}
