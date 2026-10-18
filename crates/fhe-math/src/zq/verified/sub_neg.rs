//! Modular subtraction and negation, including their in-place vector loops.

use super::{reduce1, reduce1_vt};
#[cfg(verus_keep_ghost)]
use vstd::arithmetic::div_mod::lemma_mod_multiples_vanish;
use vstd::prelude::*;

verus! {

/// Subtract canonical residues.
#[inline]
pub(crate) const fn sub_mod(a: u64, b: u64, p: u64) -> (r: u64)
    requires
        2 <= p < 0x4000_0000_0000_0000u64,
        a < p,
        b < p,
    ensures
        r < p,
        r as int == (a as int - b as int) % p as int,
{
    proof { lemma_mod_multiples_vanish(1, a as int - b as int, p as int); }
    reduce1(a + p - b, p)
}

/// Subtract canonical residues using the variable-time reducer.
///
/// # Safety
/// Timing may reveal information about the operands.
#[inline]
pub(crate) const unsafe fn sub_mod_vt(a: u64, b: u64, p: u64) -> (r: u64)
    requires
        2 <= p < 0x4000_0000_0000_0000u64,
        a < p,
        b < p,
    ensures
        r < p,
        r as int == (a as int - b as int) % p as int,
{
    proof { lemma_mod_multiples_vanish(1, a as int - b as int, p as int); }
    unsafe { reduce1_vt(a + p - b, p) }
}

/// Negate a canonical residue.
#[inline]
pub(crate) const fn neg_mod(a: u64, p: u64) -> (r: u64)
    requires
        2 <= p < 0x4000_0000_0000_0000u64,
        a < p,
    ensures
        r < p,
        r as int == (-(a as int)) % p as int,
{
    proof { lemma_mod_multiples_vanish(1, -(a as int), p as int); }
    reduce1(p - a, p)
}

/// Negate a canonical residue using the variable-time reducer.
///
/// # Safety
/// Timing may reveal information about the operands.
#[inline]
pub(crate) const unsafe fn neg_mod_vt(a: u64, p: u64) -> (r: u64)
    requires
        2 <= p < 0x4000_0000_0000_0000u64,
        a < p,
    ensures
        r < p,
        r as int == (-(a as int)) % p as int,
{
    proof { lemma_mod_multiples_vanish(1, -(a as int), p as int); }
    unsafe { reduce1_vt(p - a, p) }
}

/// Subtract vectors in place, preserving length and the unprocessed suffix.
#[inline]
pub(crate) fn sub_vec(a: &mut [u64], b: &[u64], p: u64)
    requires
        2 <= p < 0x4000_0000_0000_0000u64,
        old(a)@.len() == b@.len(),
        forall|j: int| 0 <= j < old(a)@.len() ==> #[trigger] old(a)@[j] < p,
        forall|j: int| 0 <= j < b@.len() ==> #[trigger] b@[j] < p,
    ensures
        final(a)@.len() == old(a)@.len(),
        forall|j: int| 0 <= j < final(a)@.len() ==>
            #[trigger] final(a)@[j] < p
            && final(a)@[j] as int == (old(a)@[j] as int - b@[j] as int) % p as int,
{
    let mut i = 0;
    while i < a.len() && i < b.len()
        invariant
            2 <= p < 0x4000_0000_0000_0000u64,
            a@.len() == old(a)@.len(),
            a@.len() == b@.len(),
            0 <= i <= a.len(),
            forall|k: int| 0 <= k < old(a)@.len() ==> #[trigger] old(a)@[k] < p,
            forall|k: int| 0 <= k < b@.len() ==> #[trigger] b@[k] < p,
            forall|k: int| 0 <= k < i ==>
                #[trigger] a@[k] < p
                && a@[k] as int == (old(a)@[k] as int - b@[k] as int) % p as int,
            forall|k: int| i <= k < a.len() ==> #[trigger] a@[k] == old(a)@[k],
        decreases a.len() - i,
    {
        a[i] = sub_mod(a[i], b[i], p);
        i += 1;
    }
}

/// Subtract complete blocks of 16 and then the remaining entries.
///
/// # Safety
/// Timing may reveal information about the entries.
#[inline]
pub(crate) unsafe fn sub_vec_vt(a: &mut [u64], b: &[u64], p: u64)
    requires
        2 <= p < 0x4000_0000_0000_0000u64,
        old(a)@.len() == b@.len(),
        forall|j: int| 0 <= j < old(a)@.len() ==> #[trigger] old(a)@[j] < p,
        forall|j: int| 0 <= j < b@.len() ==> #[trigger] b@[j] < p,
    ensures
        final(a)@.len() == old(a)@.len(),
        forall|j: int| 0 <= j < final(a)@.len() ==>
            #[trigger] final(a)@[j] < p
            && final(a)@[j] as int == (old(a)@[j] as int - b@[j] as int) % p as int,
{
    let mut i = 0;
    // Keep the original block size and tail structure.
    while a.len() - i >= 16
        invariant
            2 <= p < 0x4000_0000_0000_0000u64,
            a@.len() == old(a)@.len(),
            a@.len() == b@.len(),
            0 <= i <= a.len(),
            forall|k: int| 0 <= k < old(a)@.len() ==> #[trigger] old(a)@[k] < p,
            forall|k: int| 0 <= k < b@.len() ==> #[trigger] b@[k] < p,
            forall|k: int| 0 <= k < i ==>
                #[trigger] a@[k] < p
                && a@[k] as int == (old(a)@[k] as int - b@[k] as int) % p as int,
            forall|k: int| i <= k < a.len() ==> #[trigger] a@[k] == old(a)@[k],
        decreases a.len() - i,
    {
        let end = i + 16;
        let mut j = i;
        while j < end
            invariant
                2 <= p < 0x4000_0000_0000_0000u64,
                a@.len() == old(a)@.len(),
                a@.len() == b@.len(),
                0 <= j <= a.len(),
                forall|k: int| 0 <= k < old(a)@.len() ==> #[trigger] old(a)@[k] < p,
                forall|k: int| 0 <= k < b@.len() ==> #[trigger] b@[k] < p,
                forall|k: int| 0 <= k < j ==>
                    #[trigger] a@[k] < p
                    && a@[k] as int == (old(a)@[k] as int - b@[k] as int) % p as int,
                forall|k: int| j <= k < a.len() ==> #[trigger] a@[k] == old(a)@[k],
                i <= j <= end <= a.len(),
                end == i + 16,
            decreases end - j,
        {
            a[j] = unsafe { sub_mod_vt(a[j], b[j], p) };
            j += 1;
        }
        i = end;
    }

    while i < a.len()
        invariant
            2 <= p < 0x4000_0000_0000_0000u64,
            a@.len() == old(a)@.len(),
            a@.len() == b@.len(),
            0 <= i <= a.len(),
            forall|k: int| 0 <= k < old(a)@.len() ==> #[trigger] old(a)@[k] < p,
            forall|k: int| 0 <= k < b@.len() ==> #[trigger] b@[k] < p,
            forall|k: int| 0 <= k < i ==>
                #[trigger] a@[k] < p
                && a@[k] as int == (old(a)@[k] as int - b@[k] as int) % p as int,
            forall|k: int| i <= k < a.len() ==> #[trigger] a@[k] == old(a)@[k],
        decreases a.len() - i,
    {
        a[i] = unsafe { sub_mod_vt(a[i], b[i], p) };
        i += 1;
    }
}

/// Negate each entry in place, preserving the slice length.
#[inline]
pub(crate) fn neg_vec(a: &mut [u64], p: u64)
    requires
        2 <= p < 0x4000_0000_0000_0000u64,
        forall|j: int| 0 <= j < old(a)@.len() ==> #[trigger] old(a)@[j] < p,
    ensures
        final(a)@.len() == old(a)@.len(),
        forall|j: int| 0 <= j < final(a)@.len() ==>
            #[trigger] final(a)@[j] < p
            && final(a)@[j] as int == -(old(a)@[j] as int) % p as int,
{
    let mut i = 0;
    while i < a.len()
        invariant
            2 <= p < 0x4000_0000_0000_0000u64,
            a@.len() == old(a)@.len(),
            0 <= i <= a.len(),
            forall|j: int| 0 <= j < old(a)@.len() ==> #[trigger] old(a)@[j] < p,
            forall|j: int| 0 <= j < i ==>
                #[trigger] a@[j] < p
                && a@[j] as int == -(old(a)@[j] as int) % p as int,
            forall|j: int| i <= j < a.len() ==> #[trigger] a@[j] == old(a)@[j],
        decreases a.len() - i,
    {
        a[i] = neg_mod(a[i], p);
        i += 1;
    }
}

/// Negate each entry in place, preserving the slice length.
///
/// # Safety
/// Timing may reveal information about the entries.
#[inline]
pub(crate) unsafe fn neg_vec_vt(a: &mut [u64], p: u64)
    requires
        2 <= p < 0x4000_0000_0000_0000u64,
        forall|j: int| 0 <= j < old(a)@.len() ==> #[trigger] old(a)@[j] < p,
    ensures
        final(a)@.len() == old(a)@.len(),
        forall|j: int| 0 <= j < final(a)@.len() ==>
            #[trigger] final(a)@[j] < p
            && final(a)@[j] as int == -(old(a)@[j] as int) % p as int,
{
    let mut i = 0;
    while i < a.len()
        invariant
            2 <= p < 0x4000_0000_0000_0000u64,
            a@.len() == old(a)@.len(),
            0 <= i <= a.len(),
            forall|j: int| 0 <= j < old(a)@.len() ==> #[trigger] old(a)@[j] < p,
            forall|j: int| 0 <= j < i ==>
                #[trigger] a@[j] < p
                && a@[j] as int == -(old(a)@[j] as int) % p as int,
            forall|j: int| i <= j < a.len() ==> #[trigger] a@[j] == old(a)@[j],
        decreases a.len() - i,
    {
        a[i] = unsafe { neg_mod_vt(a[i], p) };
        i += 1;
    }
}

}
