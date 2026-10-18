//! Modular addition, including the in-place vector loops.

#[cfg(verus_keep_ghost)]
use super::MODULUS_LIMIT;
use super::{reduce1, reduce1_vt};
use vstd::prelude::*;

verus! {

/// Add canonical residues for a modulus accepted by `Modulus::new`.
///
/// The contract proves the result and arithmetic safety under these input bounds.
/// It does not prove constant-time machine code or validate callers' inputs.
#[inline]
pub(crate) const fn add_mod(a: u64, b: u64, p: u64) -> (result: u64)
    requires
        2 <= p < MODULUS_LIMIT,
        a < p,
        b < p,
    ensures
        result < p,
        result as int == (a as int + b as int) % (p as int),
{
    reduce1(a + b, p)
}

/// Add canonical residues using the verified variable-time reducer.
///
/// # Safety
/// Timing may reveal information about the operands.
#[inline]
pub(crate) const unsafe fn add_mod_vt(a: u64, b: u64, p: u64) -> (result: u64)
    requires
        2 <= p < MODULUS_LIMIT,
        a < p,
        b < p,
    ensures
        result < p,
        result as int == (a as int + b as int) % (p as int),
{
    unsafe { reduce1_vt(a + b, p) }
}

/// Add vectors in place, preserving length and the unprocessed suffix.
#[inline]
pub(crate) fn add_vec(a: &mut [u64], b: &[u64], p: u64)
    requires
        2 <= p < MODULUS_LIMIT,
        old(a)@.len() == b@.len(),
        forall|j: int| 0 <= j < old(a)@.len() ==> #[trigger] old(a)@[j] < p,
        forall|j: int| 0 <= j < b@.len() ==> #[trigger] b@[j] < p,
    ensures
        final(a)@.len() == old(a)@.len(),
        forall|j: int| 0 <= j < final(a)@.len() ==>
            #[trigger] final(a)@[j] < p
            && final(a)@[j] as int == (old(a)@[j] as int + b@[j] as int) % p as int,
{
    let mut i = 0;
    while i < a.len() && i < b.len()
        invariant
            2 <= p < MODULUS_LIMIT,
            a@.len() == old(a)@.len(),
            a@.len() == b@.len(),
            0 <= i <= a.len(),
            forall|k: int| 0 <= k < old(a)@.len() ==> #[trigger] old(a)@[k] < p,
            forall|k: int| 0 <= k < b@.len() ==> #[trigger] b@[k] < p,
            forall|k: int| 0 <= k < i ==>
                #[trigger] a@[k] < p
                && a@[k] as int == (old(a)@[k] as int + b@[k] as int) % p as int,
            forall|k: int| i <= k < a.len() ==> #[trigger] a@[k] == old(a)@[k],
        decreases a.len() - i,
    {
        a[i] = add_mod(a[i], b[i], p);
        i += 1;
    }
}

/// Add complete blocks of 16 and then the remaining entries.
///
/// # Safety
/// Timing may reveal information about the entries.
#[inline]
pub(crate) unsafe fn add_vec_vt(a: &mut [u64], b: &[u64], p: u64)
    requires
        2 <= p < MODULUS_LIMIT,
        old(a)@.len() == b@.len(),
        forall|j: int| 0 <= j < old(a)@.len() ==> #[trigger] old(a)@[j] < p,
        forall|j: int| 0 <= j < b@.len() ==> #[trigger] b@[j] < p,
    ensures
        final(a)@.len() == old(a)@.len(),
        forall|j: int| 0 <= j < final(a)@.len() ==>
            #[trigger] final(a)@[j] < p
            && final(a)@[j] as int == (old(a)@[j] as int + b@[j] as int) % p as int,
{
    let mut i = 0;
    // Keep the original block size and tail structure.
    while a.len() - i >= 16
        invariant
            2 <= p < MODULUS_LIMIT,
            a@.len() == old(a)@.len(),
            a@.len() == b@.len(),
            0 <= i <= a.len(),
            forall|k: int| 0 <= k < old(a)@.len() ==> #[trigger] old(a)@[k] < p,
            forall|k: int| 0 <= k < b@.len() ==> #[trigger] b@[k] < p,
            forall|k: int| 0 <= k < i ==>
                #[trigger] a@[k] < p
                && a@[k] as int == (old(a)@[k] as int + b@[k] as int) % p as int,
            forall|k: int| i <= k < a.len() ==> #[trigger] a@[k] == old(a)@[k],
        decreases a.len() - i,
    {
        let end = i + 16;
        let mut j = i;
        while j < end
            invariant
                2 <= p < MODULUS_LIMIT,
                a@.len() == old(a)@.len(),
                a@.len() == b@.len(),
                0 <= j <= a.len(),
                forall|k: int| 0 <= k < old(a)@.len() ==> #[trigger] old(a)@[k] < p,
                forall|k: int| 0 <= k < b@.len() ==> #[trigger] b@[k] < p,
                forall|k: int| 0 <= k < j ==>
                    #[trigger] a@[k] < p
                    && a@[k] as int == (old(a)@[k] as int + b@[k] as int) % p as int,
                forall|k: int| j <= k < a.len() ==> #[trigger] a@[k] == old(a)@[k],
                i <= j <= end <= a.len(),
                end == i + 16,
            decreases end - j,
        {
            a[j] = unsafe { add_mod_vt(a[j], b[j], p) };
            j += 1;
        }
        i = end;
    }

    while i < a.len()
        invariant
            2 <= p < MODULUS_LIMIT,
            a@.len() == old(a)@.len(),
            a@.len() == b@.len(),
            0 <= i <= a.len(),
            forall|k: int| 0 <= k < old(a)@.len() ==> #[trigger] old(a)@[k] < p,
            forall|k: int| 0 <= k < b@.len() ==> #[trigger] b@[k] < p,
            forall|k: int| 0 <= k < i ==>
                #[trigger] a@[k] < p
                && a@[k] as int == (old(a)@[k] as int + b@[k] as int) % p as int,
            forall|k: int| i <= k < a.len() ==> #[trigger] a@[k] == old(a)@[k],
        decreases a.len() - i,
    {
        a[i] = unsafe { add_mod_vt(a[i], b[i], p) };
        i += 1;
    }
}

}
