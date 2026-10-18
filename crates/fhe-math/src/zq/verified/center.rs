//! Center canonical residues, preserving their value modulo p.

#[cfg(verus_keep_ghost)]
use super::MODULUS_LIMIT;
use super::const_time_cond_select;
#[cfg(verus_keep_ghost)]
use vstd::arithmetic::div_mod::{lemma_mod_multiples_vanish, lemma_small_mod};
use vstd::prelude::*;

verus! {

/// Reinterpreting a signed word as unsigned and back preserves its bits.
proof fn signed_round_trip(x: i64)
    ensures (x as u64) as i64 == x,
{
    assert((x as u64) as i64 == x) by (bit_vector);
}

/// Return the unique representative in the mathematical interval [-p/2, p/2).
#[inline]
pub(crate) const fn center(a: u64, p: u64) -> (r: i64)
    requires
        2 <= p < MODULUS_LIMIT,
        a < p,
    ensures
        -(p as int) <= 2 * (r as int) < p as int,
        r as int % p as int == a as int,
        r as int == if 2 * (a as int) < p as int { a as int } else { a as int - p as int },
{
    let threshold = (p + 1) >> 1;
    let cond = a >= threshold;
    let on_true = (a as i64).wrapping_sub(p as i64) as u64;
    let on_false = a;

    let result = const_time_cond_select(on_true, on_false, cond) as i64;
    proof {
        assert(threshold == ((p + 1) as u64) / 2u64) by (bit_vector)
            requires p < MODULUS_LIMIT,
                threshold == ((p + 1) as u64) >> 1;
        let difference = (a as i64).wrapping_sub(p as i64);
        assert(difference as int == a as int - p as int);
        signed_round_trip(difference);
        assert(result as int == if cond { a as int - p as int } else { a as int });
        lemma_small_mod(a as nat, p as nat);
        if cond {
            lemma_mod_multiples_vanish(-1, a as int, p as int);
            assert(result as int % p as int == a as int);
        }
    }
    result
}

/// Allocate centered representatives in input order, preserving the length.
#[inline]
pub(crate) fn center_vec(a: &[u64], p: u64) -> (r: Vec<i64>)
    requires
        2 <= p < MODULUS_LIMIT,
        forall|j: int| 0 <= j < a@.len() ==> #[trigger] a@[j] < p,
    ensures
        r@.len() == a@.len(),
        forall|j: int| 0 <= j < r@.len() ==>
            -(p as int) <= 2 * (#[trigger] r@[j] as int) < p as int
            && r@[j] as int % p as int == a@[j] as int
            && r@[j] as int == if 2 * (a@[j] as int) < p as int {
                a@[j] as int
            } else {
                a@[j] as int - p as int
            },
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
                -(p as int) <= 2 * (#[trigger] r@[j] as int) < p as int
                && r@[j] as int % p as int == a@[j] as int
                && r@[j] as int == if 2 * (a@[j] as int) < p as int {
                    a@[j] as int
                } else {
                    a@[j] as int - p as int
                },
        decreases a.len() - i,
    {
        r.push(center(a[i], p));
        i += 1;
    }
    r
}

}
