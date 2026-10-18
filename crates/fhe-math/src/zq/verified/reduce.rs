//! One-subtraction kernels and canonical reductions composed from verified
//! lazy reductions.

#[cfg(verus_keep_ghost)]
use super::lazy_reduce::{base, valid_barrett, valid_opt};
use super::{
    const_time_cond_select, lazy_reduce, lazy_reduce_opt, lazy_reduce_opt_u128, lazy_reduce_u128,
};
#[cfg(verus_keep_ghost)]
use vstd::arithmetic::div_mod::{lemma_mod_mod, lemma_mod_multiples_vanish};
use vstd::prelude::*;

verus! {

/// Reduce a value below `2*p` to its canonical residue.
///
/// The wider modulus bound includes the NTT's calls with twice a 62-bit modulus.
#[inline]
pub(crate) const fn reduce1(x: u64, p: u64) -> (result: u64)
    requires
        0 < p < 0x8000_0000_0000_0000u64,
        x < 2 * p,
    ensures
        result < p,
        result as int == (x as int) % (p as int),
{
    let result = const_time_cond_select(x, x.wrapping_sub(p), x < p);

    proof {
        // The selection contract supplies the two cases; one subtraction suffices.
        assert((x as int) % (p as int) == if x < p {
            x as int
        } else {
            x as int - p as int
        }) by (nonlinear_arith)
            requires
                0 < p,
                x < 2 * p,
        ;
    }

    result
}

/// Reduce a value below `2*p` using the target's variable-time implementation.
///
/// # Safety
/// Timing may reveal information about `x`.
#[cfg(any(target_os = "macos", target_feature = "avx2"))]
#[inline]
pub(crate) const unsafe fn reduce1_vt(x: u64, p: u64) -> (result: u64)
    requires
        0 < p < 0x8000_0000_0000_0000u64,
        x < 2 * p,
    ensures
        result < p,
        result as int == (x as int) % (p as int),
{
    proof {
        assert((x as int) % (p as int) == if x >= p {
            x as int - p as int
        } else {
            x as int
        }) by (nonlinear_arith)
            requires
                0 < p,
                x < 2 * p,
        ;
    }

    if x >= p { x - p } else { x }
}

/// Use the verified selector-based reduction on targets without the branch path.
///
/// # Safety
/// Callers must permit variable timing, as on the other targets.
#[cfg(all(not(target_os = "macos"), not(target_feature = "avx2")))]
#[inline]
pub(crate) const unsafe fn reduce1_vt(x: u64, p: u64) -> (result: u64)
    requires
        0 < p < 0x8000_0000_0000_0000u64,
        x < 2 * p,
    ensures
        result < p,
        result as int == (x as int) % (p as int),
{
    reduce1(x, p)
}

/// Canonical u64 reduction.
#[inline]
pub(crate) const fn reduce(a: u64, p: u64, hi: u64, lo: u64) -> (r: u64)
    requires valid_barrett(p, hi, lo),
    ensures r < p, r as int == a as int % p as int,
{
    reduce1(lazy_reduce(a, p, hi, lo), p)
}

/// Canonical u64 reduction using the variable-time final subtraction.
///
/// # Safety
/// Timing may reveal information about the input.
#[inline]
pub(crate) const unsafe fn reduce_vt(a: u64, p: u64, hi: u64, lo: u64) -> (r: u64)
    requires valid_barrett(p, hi, lo),
    ensures r < p, r as int == a as int % p as int,
{
    unsafe { reduce1_vt(lazy_reduce(a, p, hi, lo), p) }
}

/// Canonical u128 reduction.
#[inline]
pub(crate) const fn reduce_u128(a: u128, p: u64, hi: u64, lo: u64) -> (r: u64)
    requires valid_barrett(p, hi, lo),
    ensures r < p, r as int == a as int % p as int,
{
    reduce1(lazy_reduce_u128(a, p, hi, lo), p)
}

/// Canonical u128 reduction using the variable-time final subtraction.
///
/// # Safety
/// Timing may reveal information about the input.
#[inline]
pub(crate) const unsafe fn reduce_u128_vt(a: u128, p: u64, hi: u64, lo: u64) -> (r: u64)
    requires valid_barrett(p, hi, lo),
    ensures r < p, r as int == a as int % p as int,
{
    unsafe { reduce1_vt(lazy_reduce_u128(a, p, hi, lo), p) }
}

/// Canonical u64 reduction.
#[inline]
pub(crate) const fn reduce_opt(a: u64, p: u64, s: u32) -> (r: u64)
    requires valid_opt(p, s),
    ensures r < p, r as int == a as int % p as int,
{
    reduce1(lazy_reduce_opt(a, p, s), p)
}

/// Canonical u64 reduction using the variable-time final subtraction.
///
/// # Safety
/// Timing may reveal information about the input.
#[inline]
pub(crate) const unsafe fn reduce_opt_vt(a: u64, p: u64, s: u32) -> (r: u64)
    requires valid_opt(p, s),
    ensures r < p, r as int == a as int % p as int,
{
    unsafe { reduce1_vt(lazy_reduce_opt(a, p, s), p) }
}

/// Canonical u128 reduction.
#[inline]
pub(crate) const fn reduce_opt_u128(a: u128, p: u64, hi: u64, lo: u64, s: u32) -> (r: u64)
    requires valid_barrett(p, hi, lo),
        valid_opt(p, s),
        a < p * p,
    ensures r < p, r as int == a as int % p as int,
{
    reduce1(lazy_reduce_opt_u128(a, p, hi, lo, s), p)
}

/// Canonical u128 reduction using the variable-time final subtraction.
///
/// # Safety
/// Timing may reveal information about the input.
#[inline]
pub(crate) const unsafe fn reduce_opt_u128_vt(a: u128, p: u64, hi: u64, lo: u64, s: u32) -> (r: u64)
    requires valid_barrett(p, hi, lo),
        valid_opt(p, s),
        a < p * p,
    ensures r < p, r as int == a as int % p as int,
{
    unsafe { reduce1_vt(lazy_reduce_opt_u128(a, p, hi, lo, s), p) }
}

/// Lift a signed input by a positive multiple of p without losing any bits.
#[inline]
const fn lift_i64(a: i64, p: u64) -> (r: u128)
    requires 2 <= p < 0x4000_0000_0000_0000u64,
    ensures
        r as int == p * base() + a,
        r as int % p as int == a as int % p as int,
{
    proof {
        assert((p as i128) << 64 == (p as i128) * 0x1_0000_0000_0000_0000i128)
            by (bit_vector)
            requires p < 0x4000_0000_0000_0000u64;
        lemma_mod_multiples_vanish(base(), a as int, p as int);
    }
    (((p as i128) << 64) + (a as i128)) as u128
}

/// Reduce any signed input to its nonnegative Euclidean remainder.
#[inline]
pub(crate) const fn reduce_i64(a: i64, p: u64, hi: u64, lo: u64) -> (r: u64)
    requires valid_barrett(p, hi, lo),
    ensures r < p, r as int == a as int % p as int,
{
    reduce_u128(lift_i64(a, p), p, hi, lo)
}

/// Reduce any signed input to its nonnegative Euclidean remainder.
///
/// # Safety
/// Timing may reveal information about the input.
#[inline]
pub(crate) const unsafe fn reduce_i64_vt(a: i64, p: u64, hi: u64, lo: u64) -> (r: u64)
    requires valid_barrett(p, hi, lo),
    ensures r < p, r as int == a as int % p as int,
{
    unsafe { reduce_u128_vt(lift_i64(a, p), p, hi, lo) }
}

/// Reduce an NTT coefficient below four times the modulus.
#[inline]
pub(crate) const fn reduce3(a: u64, p: u64, p_twice: u64) -> (r: u64)
    requires
        0 < p < 0x4000_0000_0000_0000u64,
        p_twice == 2 * p,
        a < 4 * p,
    ensures r < p, r as int == a as int % p as int,
{
    let y = reduce1(a, p_twice);
    proof { lemma_mod_mod(a as int, p as int, 2); }
    reduce1(y, p)
}

/// Reduce an NTT coefficient below four times the modulus.
///
/// # Safety
/// Timing may reveal information about the input.
#[inline]
pub(crate) const unsafe fn reduce3_vt(a: u64, p: u64, p_twice: u64) -> (r: u64)
    requires
        0 < p < 0x4000_0000_0000_0000u64,
        p_twice == 2 * p,
        a < 4 * p,
    ensures r < p, r as int == a as int % p as int,
{
    let y = unsafe { reduce1_vt(a, p_twice) };
    proof { lemma_mod_mod(a as int, p as int, 2); }
    unsafe { reduce1_vt(y, p) }
}

/// Reduce every entry in place, preserving the slice length.
#[inline]
pub(crate) fn reduce_vec(a: &mut [u64], p: u64, hi: u64, lo: u64)
    requires valid_barrett(p, hi, lo),
    ensures
        final(a)@.len() == old(a)@.len(),
        forall|i: int| 0 <= i < final(a)@.len() ==>
            #[trigger] final(a)@[i] < p
            && final(a)@[i] as int == old(a)@[i] as int % p as int,
{
    let mut i = 0;
    while i < a.len()
        invariant
            valid_barrett(p, hi, lo),
            0 <= i <= a.len(),
            a@.len() == old(a)@.len(),
            forall|j: int| 0 <= j < i ==>
                #[trigger] a@[j] < p
                && a@[j] as int == old(a)@[j] as int % p as int,
            forall|j: int| i <= j < a.len() ==> #[trigger] a@[j] == old(a)@[j],
        decreases a.len() - i,
    {
        a[i] = reduce(a[i], p, hi, lo);
        i += 1;
    }
}

/// Reduce every entry in place, preserving the slice length.
///
/// # Safety
/// Timing may reveal information about the entries.
#[inline]
pub(crate) unsafe fn reduce_vec_vt(a: &mut [u64], p: u64, hi: u64, lo: u64)
    requires valid_barrett(p, hi, lo),
    ensures
        final(a)@.len() == old(a)@.len(),
        forall|i: int| 0 <= i < final(a)@.len() ==>
            #[trigger] final(a)@[i] < p
            && final(a)@[i] as int == old(a)@[i] as int % p as int,
{
    let mut i = 0;
    while i < a.len()
        invariant
            valid_barrett(p, hi, lo),
            0 <= i <= a.len(),
            a@.len() == old(a)@.len(),
            forall|j: int| 0 <= j < i ==>
                #[trigger] a@[j] < p
                && a@[j] as int == old(a)@[j] as int % p as int,
            forall|j: int| i <= j < a.len() ==> #[trigger] a@[j] == old(a)@[j],
        decreases a.len() - i,
    {
        a[i] = unsafe { reduce_vt(a[i], p, hi, lo) };
        i += 1;
    }
}

/// Allocate canonical residues in input order, preserving the input and length.
#[inline]
pub(crate) fn reduce_vec_new(a: &[u64], p: u64, hi: u64, lo: u64) -> (r: Vec<u64>)
    requires valid_barrett(p, hi, lo),
    ensures
        r@.len() == a@.len(),
        forall|i: int| 0 <= i < r@.len() ==>
            #[trigger] r@[i] < p && r@[i] as int == a@[i] as int % p as int,
{
    let mut r: Vec<u64> = Vec::with_capacity(a.len());
    let mut i = 0;
    while i < a.len()
        invariant
            valid_barrett(p, hi, lo),
            0 <= i <= a.len(),
            r@.len() == i,
            forall|j: int| 0 <= j < i ==>
                #[trigger] r@[j] < p && r@[j] as int == a@[j] as int % p as int,
        decreases a.len() - i,
    {
        let x = reduce(a[i], p, hi, lo);
        r.push(x);
        i += 1;
    }
    r
}

/// Allocate canonical residues in input order, preserving the input and length.
///
/// # Safety
/// Timing may reveal information about the entries.
#[inline]
pub(crate) unsafe fn reduce_vec_new_vt(a: &[u64], p: u64, hi: u64, lo: u64) -> (r: Vec<u64>)
    requires valid_barrett(p, hi, lo),
    ensures
        r@.len() == a@.len(),
        forall|i: int| 0 <= i < r@.len() ==>
            #[trigger] r@[i] < p && r@[i] as int == a@[i] as int % p as int,
{
    let mut r: Vec<u64> = Vec::with_capacity(a.len());
    let mut i = 0;
    while i < a.len()
        invariant
            valid_barrett(p, hi, lo),
            0 <= i <= a.len(),
            r@.len() == i,
            forall|j: int| 0 <= j < i ==>
                #[trigger] r@[j] < p && r@[j] as int == a@[j] as int % p as int,
        decreases a.len() - i,
    {
        let x = unsafe { reduce_vt(a[i], p, hi, lo) };
        r.push(x);
        i += 1;
    }
    r
}

/// Allocate canonical residues in input order, preserving the input and length.
#[inline]
pub(crate) fn reduce_vec_i64(a: &[i64], p: u64, hi: u64, lo: u64) -> (r: Vec<u64>)
    requires valid_barrett(p, hi, lo),
    ensures
        r@.len() == a@.len(),
        forall|i: int| 0 <= i < r@.len() ==>
            #[trigger] r@[i] < p && r@[i] as int == a@[i] as int % p as int,
{
    let mut r: Vec<u64> = Vec::with_capacity(a.len());
    let mut i = 0;
    while i < a.len()
        invariant
            valid_barrett(p, hi, lo),
            0 <= i <= a.len(),
            r@.len() == i,
            forall|j: int| 0 <= j < i ==>
                #[trigger] r@[j] < p && r@[j] as int == a@[j] as int % p as int,
        decreases a.len() - i,
    {
        let x = reduce_i64(a[i], p, hi, lo);
        r.push(x);
        i += 1;
    }
    r
}

/// Allocate canonical residues in input order, preserving the input and length.
///
/// # Safety
/// Timing may reveal information about the entries.
#[inline]
pub(crate) unsafe fn reduce_vec_i64_vt(a: &[i64], p: u64, hi: u64, lo: u64) -> (r: Vec<u64>)
    requires valid_barrett(p, hi, lo),
    ensures
        r@.len() == a@.len(),
        forall|i: int| 0 <= i < r@.len() ==>
            #[trigger] r@[i] < p && r@[i] as int == a@[i] as int % p as int,
{
    let mut r: Vec<u64> = Vec::with_capacity(a.len());
    let mut i = 0;
    while i < a.len()
        invariant
            valid_barrett(p, hi, lo),
            0 <= i <= a.len(),
            r@.len() == i,
            forall|j: int| 0 <= j < i ==>
                #[trigger] r@[j] < p && r@[j] as int == a@[j] as int % p as int,
        decreases a.len() - i,
    {
        let x = unsafe { reduce_i64_vt(a[i], p, hi, lo) };
        r.push(x);
        i += 1;
    }
    r
}

}
