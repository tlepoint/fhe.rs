//! Barrett reduction kernels, eligibility checks, and their arithmetic lemmas.

#[cfg(verus_keep_ghost)]
use super::MODULUS_LIMIT;
#[cfg(verus_keep_ghost)]
use vstd::arithmetic::div_mod::*;
#[cfg(verus_keep_ghost)]
use vstd::arithmetic::power2::*;
#[cfg(verus_keep_ghost)]
use vstd::bits::*;
use vstd::prelude::*;
#[cfg(verus_keep_ghost)]
use vstd::std_specs::bits::{axiom_u64_leading_zeros, u64_leading_zeros};

verus! {

pub open spec fn base() -> int { 0x1_0000_0000_0000_0000 }

pub open spec fn valid_barrett(p: u64, hi: u64, lo: u64) -> bool {
    &&& 2 <= p < MODULUS_LIMIT
    &&& (hi * base() + lo) * p <= base() * base()
    &&& base() * base() < (hi * base() + lo + 1) * p
}

/// The strict integer form of the NFL optimization condition.
pub open spec fn optimization_condition(p: u64, s: u32) -> bool {
    let t = pow2(s as nat) as int;
    (t * t * t + 1) * base() < t * t * t * (t + 1) * p
}

/// The optimization condition together with the modulus and bit-length bounds
/// required by the lazy-reduction kernels.
pub open spec fn valid_opt(p: u64, s: u32) -> bool {
    let k = pow2((64 - s) as nat) as int;
    &&& 2 <= p < MODULUS_LIMIT
    &&& 2 <= s <= 62
    &&& k / 2 <= p < k
    &&& optimization_condition(p, s)
}

/// Connect the standard-library leading-zero specification to integer bounds.
#[verifier::nonlinear]
pub(super) proof fn leading_zero_bounds(p: u64)
    ensures
        0 <= u64_leading_zeros(p) <= 64,
        p == 0 <==> u64_leading_zeros(p) == 64,
        p > 0 ==> pow2((63 - u64_leading_zeros(p)) as nat) <= p
            < pow2((64 - u64_leading_zeros(p)) as nat),
        2 <= p < MODULUS_LIMIT ==> 2 <= u64_leading_zeros(p) <= 62,
    decreases p,
{
    axiom_u64_leading_zeros(p);
    lemma2_to64();
    lemma2_to64_rest();
    if p != 0 {
        reveal(u64_leading_zeros);
        if p > 1 {
            leading_zero_bounds((p / 2) as u64);
            let s = u64_leading_zeros(p);
            lemma_pow2_unfold((63 - s) as nat);
            lemma_pow2_unfold((64 - s) as nat);
            assert(2 * (p / 2) <= p < 2 * (p / 2 + 1));
        } else {
            axiom_u64_leading_zeros(0);
            assert(u64_leading_zeros(p) == 63);
        }
    }
}

spec fn correction(s: u32) -> int {
    if 3 * s < 64 { pow2((64 - 3 * s) as nat) as int } else { 0 }
}

/// `correction(s)` is exactly floor(B / 2^(3*s)), without a large integer.
#[verifier::nonlinear]
proof fn correction_bounds(s: u32)
    requires 1 <= s <= 63,
    ensures
        0 <= correction(s) < base(),
        correction(s) * pow2(s as nat) * pow2(s as nat) * pow2(s as nat) <= base()
            < (correction(s) + 1) * pow2(s as nat) * pow2(s as nat) * pow2(s as nat),
{
    lemma2_to64();
    lemma_pow2_pos(s as nat);
    lemma_pow2_adds(s as nat, s as nat);
    lemma_pow2_adds((2 * s) as nat, s as nat);
    if 3 * s < 64 {
        lemma_pow2_adds((3 * s) as nat, (64 - 3 * s) as nat);
        lemma_pow2_strictly_increases((64 - 3 * s) as nat, 64);
    } else {
        assert(64 < 3 * s);
        lemma_pow2_strictly_increases(64, (3 * s) as nat);
    }
}

/// Returns whether the modulus supports optimized multiplication and reduction.
///
/// Tests Equation (1) of <https://hal.archives-ouvertes.fr/hal-01242273/document>.
/// Returns false for zero and for values with the top bit set. The condition
/// does not require primality or restrict the input to `Modulus`' 62-bit range.
#[must_use]
#[verifier::nonlinear]
pub fn supports_opt(p: u64) -> (result: bool)
    ensures
        result == optimization_condition(p, u64_leading_zeros(p) as u32),
        2 <= p < MODULUS_LIMIT ==>
            result == valid_opt(p, u64_leading_zeros(p) as u32),
{
    let s = p.leading_zeros();
    proof {
        leading_zero_bounds(p);
        lemma_pow2_pos(s as nat);
        lemma2_to64();
        if 2 <= p < MODULUS_LIMIT {
            lemma_pow2_unfold((64 - s) as nat);
        }
    }
    if p == 0 || s == 0 {
        return false;
    }

    proof {
        lemma_pow2_strictly_increases(s as nat, 64);
        shl_u128(1, s);
        correction_bounds(s);
    }
    let scaled = ((1u128 << s) + 1) * (p as u128);
    let threshold = if 3 * s >= 64 {
        0u128
    } else {
        proof {
            lemma_pow2_strictly_increases((64 - 3 * s) as nat, 64);
            shl_u128(1, (64 - 3 * s) as u32);
        }
        1u128 << (64 - 3 * s)
    };

    // Divide the NFL inequality by 2^(3*s), retaining its strict comparison:
    // (2^s + 1)*p > B + floor(B / 2^(3*s)).
    scaled > 0x1_0000_0000_0000_0000u128 + threshold
}

#[verifier::nonlinear]
proof fn optimized_reciprocal(p: int, hi: int, lo: int, b: int, t: int)
    requires
        2 <= p < b,
        4 <= t,
        p * t < b,
        0 <= hi < b,
        0 <= lo < b,
        (hi * b + lo) * p <= b * b < (hi * b + lo + 1) * p,
        (t * t * t + 1) * b < t * t * t * (t + 1) * p,
    ensures
        hi == t,
        b < (t + 1) * p,
        p * p + b * lo < b * b,
{
    assert(b < (t + 1) * p);
    assert(hi == t);
    assert(p * p * p * t * t * t < b * b * b);
    assert(b * b * (t + 1) * p > b * b * b + p * p * p);
    assert(p * p + b * lo < b * b);
}

#[verifier::nonlinear]
proof fn optimized_remainder(a: int, p: int, mu: int, q: int, lo: int, al: int, b: int)
    requires
        0 <= a < p * p,
        0 < p < b,
        0 <= mu,
        0 <= q,
        0 <= lo < b,
        0 <= al < b,
        mu * p <= b * b < (mu + 1) * p,
        q * b * b <= a * mu - al * lo < (q + 1) * b * b,
        p * p + b * lo < b * b,
    ensures
        0 <= a - q * p < 2 * p,
        (a - q * p) % p == a % p,
{
    assert(q * p <= a);
    assert(a * (b * b - mu * p) + al * lo * p < p * b * b);
    assert(a < (q + 2) * p);
    lemma_mod_multiples_vanish(-q, a, p);
}

proof fn shift_parameters(s: u32)
    requires 2 <= s <= 62,
    ensures
        4 <= pow2(s as nat),
        0 < pow2((64 - s) as nat),
        pow2(s as nat) * pow2((64 - s) as nat) == base(),
{
    lemma2_to64();
    if s > 2 { lemma_pow2_strictly_increases(2, s as nat); }
    lemma_pow2_pos((64 - s) as nat);
    lemma_pow2_adds(s as nat, (64 - s) as nat);
}

/// One step of the connection between machine shifts and integer multiplication.
proof fn shl_u128_bits(x: u128, s: u32, prev: u128)
    requires
        0 < s < 128,
        prev == x << ((s - 1) as u32),
        prev < 0x8000_0000_0000_0000_0000_0000_0000_0000u128,
    ensures
        (x << s) as int == 2 * prev,
{
    assert(x << s == prev << 1) by (bit_vector)
        requires 0 < s < 128, prev == x << ((s - 1) as u32);
    assert(prev << 1 == 2 * prev) by (bit_vector)
        requires prev < 0x8000_0000_0000_0000_0000_0000_0000_0000u128;
}

proof fn shl_u128_zero(x: u128)
    ensures x << 0 == x,
{
    assert(x << 0 == x) by (bit_vector);
}

/// Relate a non-overflowing 128-bit left shift to mathematical multiplication.
#[verifier::nonlinear]
proof fn shl_u128(x: u128, s: u32)
    requires
        s < 128,
        x * pow2(s as nat) < base() * base(),
    ensures
        (x << s) as int == x * pow2(s as nat),
    decreases s,
{
    if s == 0 {
        lemma2_to64();
        shl_u128_zero(x);
    } else {
        lemma_pow2_unfold(s as nat);
        lemma_pow2_pos((s - 1) as nat);
        shl_u128(x, (s - 1) as u32);
        let prev = x << ((s - 1) as u32);
        shl_u128_bits(x, s, prev);
    }
}

/// The reciprocal's two limbs sum to at most the radix. This also bounds
/// the intermediate sum in the full-width implementation.
#[verifier::nonlinear]
proof fn reciprocal_limbs(p: int, hi: int, lo: int, b: int)
    requires
        2 <= p < b,
        0 <= hi < b,
        0 <= lo < b,
        (hi * b + lo) * p <= b * b,
        b * b < (hi * b + lo + 1) * p,
    ensures
        hi * p <= b,
        b < (hi + 1) * p,
        hi + lo <= b,
{
    assert(hi * p <= b);
    assert(b < (hi + 1) * p);
    assert(lo * p <= b * (b - hi * p));
    assert(hi + lo <= b);
}

/// A floor reciprocal underestimates the quotient by at most one.
#[verifier::nonlinear]
proof fn barrett_remainder(a: int, p: int, mu: int, q: int, scale: int)
    requires
        0 <= a < scale,
        0 < p,
        0 <= mu,
        0 <= q,
        mu * p <= scale < (mu + 1) * p,
        q * scale <= a * mu < (q + 1) * scale,
    ensures
        0 <= a - q * p < 2 * p,
        (a - q * p) % p == a % p,
{
    assert(q * p <= a);
    assert(a < (q + 2) * p);
    lemma_mod_multiples_vanish(-q, a, p);
}

pub(super) proof fn split_u128(x: u128)
    ensures
        (x >> 64) as int == x as int / base(),
        (x as u64) as int == x as int % base(),
        x as int == (x >> 64) as int * base() + (x as u64) as int,
{
    assert(x >> 64 == x / 0x1_0000_0000_0000_0000u128) by (bit_vector);
    assert((x as u64) as u128 == x % 0x1_0000_0000_0000_0000u128) by (bit_vector);
    lemma_fundamental_div_mod(x as int, base());
}

/// Reduce any 128-bit input to a congruent value below twice the modulus.
#[inline]
#[verifier::nonlinear]
pub(crate) const fn lazy_reduce_u128(a: u128, p: u64, hi: u64, lo: u64) -> (r: u64)
    requires valid_barrett(p, hi, lo),
    ensures
        r < 2 * p,
        r as int % p as int == a as int % p as int,
{
    let a_lo = a as u64;
    let a_hi = (a >> 64) as u64;
    let p_lo_lo = ((a_lo as u128) * (lo as u128)) >> 64;
    let p_hi_lo = (a_hi as u128) * (lo as u128);
    let p_lo_hi = (a_lo as u128) * (hi as u128);

    proof {
        split_u128(a);
        split_u128((a_lo * lo) as u128);
        reciprocal_limbs(p as int, hi as int, lo as int, base());
        assert(p_lo_lo <= lo as u128);
        assert(p_lo_hi + p_hi_lo + p_lo_lo < base() * base());
    }

    let middle = p_lo_hi + p_hi_lo + p_lo_lo;
    proof { split_u128(middle); }
    let q = (middle >> 64) + (a_hi as u128) * (hi as u128);

    proof {
        let mu = hi * base() + lo;
        assert(q * base() * base() <= a * mu < (q + 1) * base() * base());
        barrett_remainder(a as int, p as int, mu, q as int, base() * base());
    }

    (a - q * (p as u128)) as u64
}

/// The single-limb specialization of Barrett reduction.
#[inline]
#[verifier::nonlinear]
pub(crate) const fn lazy_reduce(a: u64, p: u64, hi: u64, lo: u64) -> (r: u64)
    requires valid_barrett(p, hi, lo),
    ensures
        r < 2 * p,
        r as int % p as int == a as int % p as int,
{
    let p_lo_lo = ((a as u128) * (lo as u128)) >> 64;
    let p_lo_hi = (a as u128) * (hi as u128);
    proof {
        split_u128((a * lo) as u128);
        reciprocal_limbs(p as int, hi as int, lo as int, base());
        assert(p_lo_lo <= lo as u128);
    }

    let middle = p_lo_hi + p_lo_lo;
    proof { split_u128(middle); }
    let q = middle >> 64;

    proof {
        let mu = hi * base() + lo;
        assert(q * base() * base() <= a * mu < (q + 1) * base() * base());
        barrett_remainder(a as int, p as int, mu, q as int, base() * base());
    }

    (a as u128 - q * (p as u128)) as u64
}

/// Optimized reduction for inputs below `p^2` and an eligible modulus.
/// `_hi` identifies the full reciprocal in the proof; it is unused at runtime.
#[inline]
#[verifier::nonlinear]
pub(crate) const fn lazy_reduce_opt_u128(a: u128, p: u64, _hi: u64, lo: u64, s: u32) -> (r: u64)
    requires
        valid_barrett(p, _hi, lo),
        valid_opt(p, s),
        a < p * p,
    ensures
        r < 2 * p,
        r as int % p as int == a as int % p as int,
{
    proof {
        shift_parameters(s);
        optimized_reciprocal(p as int, _hi as int, lo as int, base(), pow2(s as nat) as int);
        shl_u128(a, s);
        split_u128(a);
    }

    let middle = ((lo as u128) * (a >> 64)) + (a << s);
    proof { split_u128(middle); }
    let q = middle >> 64;

    proof {
        let mu = _hi * base() + lo;
        assert(q * base() * base() <= a * mu - (a as u64) * lo
            < (q + 1) * base() * base());
        optimized_remainder(a as int, p as int, mu, q as int, lo as int, (a as u64) as int, base());
    }

    (a - q * (p as u128)) as u64
}

/// Optimized single-limb reduction for an eligible modulus.
#[inline]
#[verifier::nonlinear]
pub(crate) const fn lazy_reduce_opt(a: u64, p: u64, s: u32) -> (r: u64)
    requires valid_opt(p, s),
    ensures
        r < 2 * p,
        r as int % p as int == a as int % p as int,
{
    proof {
        shift_parameters(s);
        lemma_u64_shr_is_div(a, (64 - s) as u64);
        let t = pow2(s as nat) as int;
        assert(base() < (t + 1) * p);
    }

    let q = a >> (64 - s);
    proof {
        let k = pow2((64 - s) as nat) as int;
        let t = pow2(s as nat) as int;
        assert(q * k <= a < (q + 1) * k);
        assert(q < t);
        assert(0 <= a - q * p < 2 * p);
        lemma_mod_multiples_vanish(-(q as int), a as int, p as int);
    }

    ((a as u128) - (q as u128) * (p as u128)) as u64
}

/// Reduce every entry in place, preserving its residue and the slice length.
/// The const parameter keeps selection of the implementation outside the loop.
#[inline]
pub(crate) fn lazy_reduce_vec<const OPTIMIZED: bool>(
    a: &mut [u64], p: u64, hi: u64, lo: u64, s: u32,
)
    requires
        valid_barrett(p, hi, lo),
        OPTIMIZED ==> valid_opt(p, s),
    ensures
        final(a)@.len() == old(a)@.len(),
        forall|i: int| 0 <= i < final(a)@.len() ==>
            #[trigger] final(a)@[i] < 2 * p
            && final(a)@[i] as int % p as int == old(a)@[i] as int % p as int,
{
    let mut i = 0;
    while i < a.len()
        invariant
            valid_barrett(p, hi, lo),
            OPTIMIZED ==> valid_opt(p, s),
            0 <= i <= a.len(),
            a@.len() == old(a)@.len(),
            forall|j: int| 0 <= j < i ==>
                #[trigger] a@[j] < 2 * p
                && a@[j] as int % p as int == old(a)@[j] as int % p as int,
            forall|j: int| i <= j < a.len() ==> #[trigger] a@[j] == old(a)@[j],
        decreases a.len() - i,
    {
        let x = a[i];
        a[i] = if OPTIMIZED {
            lazy_reduce_opt(x, p, s)
        } else {
            lazy_reduce(x, p, hi, lo)
        };
        i += 1;
    }
}

}
