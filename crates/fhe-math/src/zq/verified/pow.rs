//! Modular exponentiation by left-to-right square-and-multiply.

#[cfg(verus_keep_ghost)]
use super::lazy_reduce::{leading_zero_bounds, valid_barrett};
use super::mul_mod;
#[cfg(verus_keep_ghost)]
use vstd::arithmetic::div_mod::{lemma_mul_mod_noop_general, lemma_small_mod};
#[cfg(verus_keep_ghost)]
use vstd::arithmetic::power::{lemma_pow_adds, lemma_pow0, lemma_pow1, pow as int_pow};
#[cfg(verus_keep_ghost)]
use vstd::arithmetic::power2::{lemma_pow2_pos, lemma_pow2_unfold, pow2};
#[cfg(verus_keep_ghost)]
use vstd::bits::lemma_u64_shr_is_div;
use vstd::prelude::*;
#[cfg(verus_keep_ghost)]
use vstd::std_specs::bits::u64_leading_zeros;

verus! {

/// The highest set bit is the initial one-bit exponent prefix.
#[verifier::nonlinear]
proof fn highest_bit(n: u64)
    requires 2 <= n < 0x4000_0000_0000_0000u64,
    ensures
        2 <= u64_leading_zeros(n) <= 62,
        n >> ((63 - u64_leading_zeros(n)) as u64) == 1,
{
    leading_zero_bounds(n);
    let k = (63 - u64_leading_zeros(n)) as nat;
    lemma_pow2_pos(k);
    lemma_pow2_unfold((k + 1) as nat);
    lemma_u64_shr_is_div(n, k as u64);
    assert(pow2(k) <= n < 2 * pow2(k));
    assert(n as nat / pow2(k) == 1);
}

/// Appending one binary digit doubles the prefix and adds that digit.
proof fn prefix_step(n: u64, i: u64)
    requires i < 63,
    ensures
        (n >> i) as int == 2 * ((n >> ((i + 1) as u64)) as int)
            + (((n >> i) & 1) as int),
        (n >> i) & 1 == 0 || (n >> i) & 1 == 1,
{
    assert((n >> i) as int == 2 * ((n >> ((i + 1) as u64)) as int)
        + (((n >> i) & 1) as int)) by (bit_vector)
        requires i < 63;
    assert((n >> i) & 1 == 0 || (n >> i) & 1 == 1) by (bit_vector);
}

/// Multiplication adds the exponents represented by canonical residues.
proof fn multiply_powers(a: int, e1: nat, e2: nat, p: int, r1: int, r2: int)
    requires
        0 < p,
        r1 == int_pow(a, e1) % p,
        r2 == int_pow(a, e2) % p,
    ensures (r1 * r2) % p == int_pow(a, e1 + e2) % p,
{
    lemma_pow_adds(a, e1, e2);
    lemma_mul_mod_noop_general(int_pow(a, e1), int_pow(a, e2), p);
}

/// Compute a^n modulo p with the original variable-time bit traversal.
#[inline]
pub(crate) fn pow_mod(a: u64, n: u64, p: u64, hi: u64, lo: u64) -> (result: u64)
    requires valid_barrett(p, hi, lo), a < p, n < p,
    ensures result < p, result as int == int_pow(a as int, n as nat) % p as int,
{
    proof {
        lemma_pow0(a as int);
        lemma_pow1(a as int);
        lemma_small_mod(1, p as nat);
        lemma_small_mod(a as nat, p as nat);
    }
    if n == 0 {
        1
    } else if n == 1 {
        a
    } else {
        proof { highest_bit(n); }
        let mut r = a;
        let mut i = (62 - n.leading_zeros()) as isize;
        while i >= 0
            invariant
                valid_barrett(p, hi, lo),
                a < p,
                2 <= n < p,
                -1 <= i <= 60,
                r < p,
                r as int == int_pow(a as int, (n >> ((i + 1) as u64)) as nat) % p as int,
            decreases i + 1,
        {
            let ghost prefix = (n >> ((i + 1) as u64)) as nat;
            proof { multiply_powers(a as int, prefix, prefix, p as int, r as int, r as int); }
            r = mul_mod(r, r, p, hi, lo);
            if (n >> i) & 1 == 1 {
                proof {
                    lemma_pow1(a as int);
                    lemma_small_mod(a as nat, p as nat);
                    multiply_powers(a as int, prefix + prefix, 1, p as int, r as int, a as int);
                }
                r = mul_mod(r, a, p, hi, lo);
            }
            proof {
                prefix_step(n, i as u64);
                assert(r as int == int_pow(a as int, (n >> i) as nat) % p as int);
            }
            i -= 1;
        }
        proof { assert(n >> 0 == n) by (bit_vector); }
        r
    }
}

}
