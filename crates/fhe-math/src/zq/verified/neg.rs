//! Modular negation, including the in-place vector loops.

#[cfg(verus_keep_ghost)]
use super::MODULUS_LIMIT;
use super::{reduce1, reduce1_vt};
#[cfg(verus_keep_ghost)]
use vstd::arithmetic::div_mod::lemma_mod_multiples_vanish;
use vstd::prelude::*;

verus! {

/// Negate a canonical residue.
#[inline]
pub(crate) const fn neg_mod(a: u64, p: u64) -> (r: u64)
    requires
        2 <= p < MODULUS_LIMIT,
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
        2 <= p < MODULUS_LIMIT,
        a < p,
    ensures
        r < p,
        r as int == (-(a as int)) % p as int,
{
    proof { lemma_mod_multiples_vanish(1, -(a as int), p as int); }
    unsafe { reduce1_vt(p - a, p) }
}

/// Negate each entry in place, preserving the slice length.
#[inline]
pub(crate) fn neg_vec(a: &mut [u64], p: u64)
    requires
        2 <= p < MODULUS_LIMIT,
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
            2 <= p < MODULUS_LIMIT,
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
        2 <= p < MODULUS_LIMIT,
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
            2 <= p < MODULUS_LIMIT,
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
