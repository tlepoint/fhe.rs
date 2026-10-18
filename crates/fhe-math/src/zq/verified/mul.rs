//! Modular multiplication, Shoup precomputation, and in-place vector kernels.

#[cfg(verus_keep_ghost)]
use super::MODULUS_LIMIT;
#[cfg(verus_keep_ghost)]
use super::lazy_reduce::{base, split_u128, valid_barrett, valid_opt};
use super::{
    reduce_opt_u128, reduce_opt_u128_vt, reduce_u128, reduce_u128_vt, reduce1, reduce1_vt,
};
#[cfg(verus_keep_ghost)]
use vstd::arithmetic::div_mod::{lemma_fundamental_div_mod, lemma_mod_multiples_vanish};
use vstd::prelude::*;

verus! {

/// Multiply canonical residues using the full-width Barrett reducer.
#[inline]
#[verifier::nonlinear]
pub(crate) const fn mul_mod(a: u64, b: u64, p: u64, hi: u64, lo: u64) -> (r: u64)
    requires valid_barrett(p, hi, lo), a < p, b < p,
    ensures r < p, r as int == (a as int * b as int) % p as int,
{
    reduce_u128((a as u128) * (b as u128), p, hi, lo)
}

/// Multiply canonical residues using variable-time final reduction.
///
/// # Safety
/// Timing may reveal information about the operands.
#[inline]
#[verifier::nonlinear]
pub(crate) const unsafe fn mul_mod_vt(a: u64, b: u64, p: u64, hi: u64, lo: u64) -> (r: u64)
    requires valid_barrett(p, hi, lo), a < p, b < p,
    ensures r < p, r as int == (a as int * b as int) % p as int,
{
    unsafe { reduce_u128_vt((a as u128) * (b as u128), p, hi, lo) }
}

/// Multiply canonical residues using the optimized Barrett reducer.
#[inline]
#[verifier::nonlinear]
pub(crate) const fn mul_opt(a: u64, b: u64, p: u64, hi: u64, lo: u64, s: u32) -> (r: u64)
    requires valid_barrett(p, hi, lo), valid_opt(p, s), a < p, b < p,
    ensures r < p, r as int == (a as int * b as int) % p as int,
{
    reduce_opt_u128((a as u128) * (b as u128), p, hi, lo, s)
}

/// Multiply canonical residues using optimized variable-time reduction.
///
/// # Safety
/// Timing may reveal information about the operands.
#[inline]
#[verifier::nonlinear]
pub(crate) const unsafe fn mul_opt_vt(a: u64, b: u64, p: u64, hi: u64, lo: u64, s: u32) -> (r: u64)
    requires valid_barrett(p, hi, lo), valid_opt(p, s), a < p, b < p,
    ensures r < p, r as int == (a as int * b as int) % p as int,
{
    unsafe { reduce_opt_u128_vt((a as u128) * (b as u128), p, hi, lo, s) }
}

/// The precomputed word is floor(b * 2^64 / p).
pub open spec fn valid_shoup(b: u64, p: u64, b_shoup: u64) -> bool {
    &&& 2 <= p < MODULUS_LIMIT
    &&& b < p
    &&& b_shoup * p <= b * base() < (b_shoup + 1) * p
}

proof fn shoup_shift(a: u64)
    ensures ((a as u128) << 64) as int == a * base(),
{
    assert((a as u128) << 64 == (a as u128) * 0x1_0000_0000_0000_0000u128)
        by (bit_vector);
}

/// Compute the exact Shoup word without truncating the quotient.
#[inline]
#[verifier::nonlinear]
pub(crate) const fn shoup(a: u64, p: u64) -> (r: u64)
    requires 2 <= p < MODULUS_LIMIT, a < p,
    ensures valid_shoup(a, p, r), r as int == a * base() / (p as int),
{
    proof {
        shoup_shift(a);
        lemma_fundamental_div_mod(a * base(), p as int);
        assert(0 <= a * base() / (p as int) < base());
    }
    (((a as u128) << 64) / (p as u128)) as u64
}

/// The Shoup quotient underestimates a*b/p by at most one, even for a >= p.
#[verifier::nonlinear]
proof fn shoup_remainder(a: int, b: int, p: int, s: int, q: int, radix: int)
    requires
        0 <= a < radix,
        0 <= b < p,
        0 < p < radix,
        0 <= s,
        0 <= q,
        s * p <= b * radix < (s + 1) * p,
        q * radix <= a * s < (q + 1) * radix,
    ensures
        0 <= a * b - q * p < 2 * p,
        (a * b - q * p) % p == (a * b) % p,
{
    assert(q * p <= a * b);
    assert(a * b < (q + 2) * p);
    lemma_mod_multiples_vanish(-q, a * b, p);
}

/// Multiply any u64 by a canonical residue, returning a residue below 2*p.
#[inline]
#[verifier::nonlinear]
pub(crate) const fn lazy_mul_shoup(a: u64, b: u64, b_shoup: u64, p: u64) -> (r: u64)
    requires valid_shoup(b, p, b_shoup),
    ensures r < 2 * p, r as int % p as int == (a as int * b as int) % p as int,
{
    let q = ((a as u128) * (b_shoup as u128)) >> 64;
    proof {
        split_u128((a * b_shoup) as u128);
        shoup_remainder(a as int, b as int, p as int, b_shoup as int, q as int, base());
    }
    ((a as u128) * (b as u128) - q * (p as u128)) as u64
}

/// Canonical Shoup multiplication; the first operand may be any u64.
#[inline]
pub(crate) const fn mul_shoup(a: u64, b: u64, b_shoup: u64, p: u64) -> (r: u64)
    requires valid_shoup(b, p, b_shoup),
    ensures r < p, r as int == (a as int * b as int) % p as int,
{
    reduce1(lazy_mul_shoup(a, b, b_shoup, p), p)
}

/// Canonical Shoup multiplication using variable-time final reduction.
///
/// # Safety
/// Timing may reveal information about the operands.
#[inline]
pub(crate) const unsafe fn mul_shoup_vt(a: u64, b: u64, b_shoup: u64, p: u64) -> (r: u64)
    requires valid_shoup(b, p, b_shoup),
    ensures r < p, r as int == (a as int * b as int) % p as int,
{
    unsafe { reduce1_vt(lazy_mul_shoup(a, b, b_shoup, p), p) }
}

/// Multiply vectors in place; the const parameter selects the Barrett reducer.
#[inline]
pub(crate) fn mul_vec<const OPTIMIZED: bool>(
    a: &mut [u64], b: &[u64], p: u64, hi: u64, lo: u64, s: u32,
)
    requires
        valid_barrett(p, hi, lo),
        OPTIMIZED ==> valid_opt(p, s),
        old(a)@.len() == b@.len(),
        forall|j: int| 0 <= j < old(a)@.len() ==> #[trigger] old(a)@[j] < p,
        forall|j: int| 0 <= j < b@.len() ==> #[trigger] b@[j] < p,
    ensures
        final(a)@.len() == old(a)@.len(),
        forall|j: int| 0 <= j < final(a)@.len() ==>
            #[trigger] final(a)@[j] < p
            && final(a)@[j] as int == (old(a)@[j] as int * b@[j] as int) % p as int,
{
    let mut i = 0;
    while i < a.len() && i < b.len()
        invariant
            valid_barrett(p, hi, lo),
            OPTIMIZED ==> valid_opt(p, s),
            a@.len() == old(a)@.len(),
            a@.len() == b@.len(),
            0 <= i <= a.len(),
            forall|k: int| 0 <= k < old(a)@.len() ==> #[trigger] old(a)@[k] < p,
            forall|k: int| 0 <= k < b@.len() ==> #[trigger] b@[k] < p,
            forall|k: int| 0 <= k < i ==>
                #[trigger] a@[k] < p
                && a@[k] as int == (old(a)@[k] as int * b@[k] as int) % p as int,
            forall|k: int| i <= k < a.len() ==> #[trigger] a@[k] == old(a)@[k],
        decreases a.len() - i,
    {
        a[i] = if OPTIMIZED {
            mul_opt(a[i], b[i], p, hi, lo, s)
        } else {
            mul_mod(a[i], b[i], p, hi, lo)
        };
        i += 1;
    }
}

/// Multiply vectors in place; the const parameter selects the Barrett reducer.
///
/// # Safety
/// Timing may reveal information about the entries.
#[inline]
pub(crate) unsafe fn mul_vec_vt<const OPTIMIZED: bool>(
    a: &mut [u64], b: &[u64], p: u64, hi: u64, lo: u64, s: u32,
)
    requires
        valid_barrett(p, hi, lo),
        OPTIMIZED ==> valid_opt(p, s),
        old(a)@.len() == b@.len(),
        forall|j: int| 0 <= j < old(a)@.len() ==> #[trigger] old(a)@[j] < p,
        forall|j: int| 0 <= j < b@.len() ==> #[trigger] b@[j] < p,
    ensures
        final(a)@.len() == old(a)@.len(),
        forall|j: int| 0 <= j < final(a)@.len() ==>
            #[trigger] final(a)@[j] < p
            && final(a)@[j] as int == (old(a)@[j] as int * b@[j] as int) % p as int,
{
    let mut i = 0;
    while i < a.len() && i < b.len()
        invariant
            valid_barrett(p, hi, lo),
            OPTIMIZED ==> valid_opt(p, s),
            a@.len() == old(a)@.len(),
            a@.len() == b@.len(),
            0 <= i <= a.len(),
            forall|k: int| 0 <= k < old(a)@.len() ==> #[trigger] old(a)@[k] < p,
            forall|k: int| 0 <= k < b@.len() ==> #[trigger] b@[k] < p,
            forall|k: int| 0 <= k < i ==>
                #[trigger] a@[k] < p
                && a@[k] as int == (old(a)@[k] as int * b@[k] as int) % p as int,
            forall|k: int| i <= k < a.len() ==> #[trigger] a@[k] == old(a)@[k],
        decreases a.len() - i,
    {
        a[i] = if OPTIMIZED {
            unsafe { mul_opt_vt(a[i], b[i], p, hi, lo, s) }
        } else {
            unsafe { mul_mod_vt(a[i], b[i], p, hi, lo) }
        };
        i += 1;
    }
}

/// Shoup multiplication in place; entries of a may be any u64.
#[inline]
pub(crate) fn mul_shoup_vec(a: &mut [u64], b: &[u64], b_shoup: &[u64], p: u64)
    requires
        2 <= p < MODULUS_LIMIT,
        old(a)@.len() == b@.len(),
        b@.len() == b_shoup@.len(),
        forall|j: int| 0 <= j < b@.len() ==> valid_shoup(b@[j], p, #[trigger] b_shoup@[j]),
    ensures
        final(a)@.len() == old(a)@.len(),
        forall|j: int| 0 <= j < final(a)@.len() ==>
            #[trigger] final(a)@[j] < p
            && final(a)@[j] as int == (old(a)@[j] as int * b@[j] as int) % p as int,
{
    let mut i = 0;
    while i < a.len() && i < b.len()
        invariant
            2 <= p < MODULUS_LIMIT,
            a@.len() == old(a)@.len(),
            a@.len() == b@.len(),
            0 <= i <= a.len(),
            b@.len() == b_shoup@.len(),
            forall|k: int| 0 <= k < b@.len() ==> valid_shoup(b@[k], p, #[trigger] b_shoup@[k]),
            forall|k: int| 0 <= k < i ==>
                #[trigger] a@[k] < p
                && a@[k] as int == (old(a)@[k] as int * b@[k] as int) % p as int,
            forall|k: int| i <= k < a.len() ==> #[trigger] a@[k] == old(a)@[k],
        decreases a.len() - i,
    {
        a[i] = mul_shoup(a[i], b[i], b_shoup[i], p);
        i += 1;
    }
}

/// Shoup multiplication in place; entries of a may be any u64.
///
/// # Safety
/// Timing may reveal information about the entries.
#[inline]
pub(crate) unsafe fn mul_shoup_vec_vt(a: &mut [u64], b: &[u64], b_shoup: &[u64], p: u64)
    requires
        2 <= p < MODULUS_LIMIT,
        old(a)@.len() == b@.len(),
        b@.len() == b_shoup@.len(),
        forall|j: int| 0 <= j < b@.len() ==> valid_shoup(b@[j], p, #[trigger] b_shoup@[j]),
    ensures
        final(a)@.len() == old(a)@.len(),
        forall|j: int| 0 <= j < final(a)@.len() ==>
            #[trigger] final(a)@[j] < p
            && final(a)@[j] as int == (old(a)@[j] as int * b@[j] as int) % p as int,
{
    let mut i = 0;
    while i < a.len() && i < b.len()
        invariant
            2 <= p < MODULUS_LIMIT,
            a@.len() == old(a)@.len(),
            a@.len() == b@.len(),
            0 <= i <= a.len(),
            b@.len() == b_shoup@.len(),
            forall|k: int| 0 <= k < b@.len() ==> valid_shoup(b@[k], p, #[trigger] b_shoup@[k]),
            forall|k: int| 0 <= k < i ==>
                #[trigger] a@[k] < p
                && a@[k] as int == (old(a)@[k] as int * b@[k] as int) % p as int,
            forall|k: int| i <= k < a.len() ==> #[trigger] a@[k] == old(a)@[k],
        decreases a.len() - i,
    {
        a[i] = unsafe { mul_shoup_vt(a[i], b[i], b_shoup[i], p) };
        i += 1;
    }
}

/// Multiply arbitrary u64 entries by one canonical scalar, precomputed once.
#[inline]
pub(crate) fn scalar_mul_vec(a: &mut [u64], b: u64, p: u64)
    requires
        2 <= p < MODULUS_LIMIT,
        b < p,
    ensures
        final(a)@.len() == old(a)@.len(),
        forall|j: int| 0 <= j < final(a)@.len() ==>
            #[trigger] final(a)@[j] < p
            && final(a)@[j] as int == (old(a)@[j] as int * b as int) % p as int,
{
    let b_shoup = shoup(b, p);
    let mut i = 0;
    while i < a.len()
        invariant
            2 <= p < MODULUS_LIMIT,
            a@.len() == old(a)@.len(),
            valid_shoup(b, p, b_shoup),
            0 <= i <= a.len(),
            forall|k: int| 0 <= k < i ==>
                #[trigger] a@[k] < p
                && a@[k] as int == (old(a)@[k] as int * b as int) % p as int,
            forall|k: int| i <= k < a.len() ==> #[trigger] a@[k] == old(a)@[k],
        decreases a.len() - i,
    {
        a[i] = mul_shoup(a[i], b, b_shoup, p);
        i += 1;
    }
}

/// Multiply arbitrary u64 entries by one canonical scalar, precomputed once.
///
/// # Safety
/// Timing may reveal information about the entries.
#[inline]
pub(crate) unsafe fn scalar_mul_vec_vt(a: &mut [u64], b: u64, p: u64)
    requires
        2 <= p < MODULUS_LIMIT,
        b < p,
    ensures
        final(a)@.len() == old(a)@.len(),
        forall|j: int| 0 <= j < final(a)@.len() ==>
            #[trigger] final(a)@[j] < p
            && final(a)@[j] as int == (old(a)@[j] as int * b as int) % p as int,
{
    let b_shoup = shoup(b, p);
    let mut i = 0;
    while i < a.len()
        invariant
            2 <= p < MODULUS_LIMIT,
            a@.len() == old(a)@.len(),
            valid_shoup(b, p, b_shoup),
            0 <= i <= a.len(),
            forall|k: int| 0 <= k < i ==>
                #[trigger] a@[k] < p
                && a@[k] as int == (old(a)@[k] as int * b as int) % p as int,
            forall|k: int| i <= k < a.len() ==> #[trigger] a@[k] == old(a)@[k],
        decreases a.len() - i,
    {
        a[i] = unsafe { mul_shoup_vt(a[i], b, b_shoup, p) };
        i += 1;
    }
}

/// Precompute Shoup words in input order.
#[inline]
pub(crate) fn shoup_vec(a: &[u64], p: u64) -> (r: Vec<u64>)
    requires
        2 <= p < MODULUS_LIMIT,
        forall|j: int| 0 <= j < a@.len() ==> #[trigger] a@[j] < p,
    ensures
        r@.len() == a@.len(),
        forall|j: int| 0 <= j < r@.len() ==>
            valid_shoup(a@[j], p, #[trigger] r@[j])
            && r@[j] as int == a@[j] * base() / (p as int),
{
    let mut r = Vec::with_capacity(a.len());
    let mut i = 0;
    while i < a.len()
        invariant
            2 <= p < MODULUS_LIMIT,
            forall|j: int| 0 <= j < a@.len() ==> #[trigger] a@[j] < p,
            0 <= i <= a.len(),
            r@.len() == i,
            forall|j: int| 0 <= j < i ==>
                valid_shoup(a@[j], p, #[trigger] r@[j])
                && r@[j] as int == a@[j] * base() / (p as int),
        decreases a.len() - i,
    {
        r.push(shoup(a[i], p));
        i += 1;
    }
    r
}

/// Fill an equally sized output slice with Shoup words.
#[inline]
pub(crate) fn shoup_into(a: &[u64], out: &mut [u64], p: u64)
    requires
        2 <= p < MODULUS_LIMIT,
        a@.len() == old(out)@.len(),
        forall|j: int| 0 <= j < a@.len() ==> #[trigger] a@[j] < p,
    ensures
        final(out)@.len() == old(out)@.len(),
        forall|j: int| 0 <= j < final(out)@.len() ==>
            valid_shoup(a@[j], p, #[trigger] final(out)@[j])
            && final(out)@[j] as int == a@[j] * base() / (p as int),
{
    let mut i = 0;
    while i < a.len()
        invariant
            2 <= p < MODULUS_LIMIT,
            forall|j: int| 0 <= j < a@.len() ==> #[trigger] a@[j] < p,
            a@.len() == out@.len(),
            out@.len() == old(out)@.len(),
            0 <= i <= a.len(),
            forall|j: int| 0 <= j < i ==>
                valid_shoup(a@[j], p, #[trigger] out@[j])
                && out@[j] as int == a@[j] * base() / (p as int),
            forall|j: int| i <= j < out.len() ==> #[trigger] out@[j] == old(out)@[j],
        decreases a.len() - i,
    {
        out[i] = shoup(a[i], p);
        i += 1;
    }
}

}
