//! Prime search and optimization eligibility for the modular arithmetic layer.

#![cfg_attr(verus_keep_ghost, allow(unused_imports))]

use super::backend::{Backend, Runtime};
#[cfg(verus_keep_ghost)]
use super::verified::inversion::prime;
#[cfg(verus_keep_ghost)]
use vstd::arithmetic::div_mod::{
    lemma_fundamental_div_mod, lemma_fundamental_div_mod_converse_mod, lemma_mod_bound,
    lemma_mod_sub_multiples_vanish, lemma_small_mod,
};
#[cfg(verus_keep_ghost)]
use vstd::arithmetic::mul::{
    lemma_mul_inequality, lemma_mul_is_commutative, lemma_mul_strict_inequality_converse,
};
#[cfg(verus_keep_ghost)]
use vstd::arithmetic::power2::{lemma_pow2_strictly_increases, lemma2_to64, pow2};
#[cfg(verus_keep_ghost)]
use vstd::bits::lemma_u64_shl_is_mul;
use vstd::prelude::*;

#[doc(inline)]
pub use super::verified::supports_opt;

verus! {

/// Supported search parameters; a small upper bound may describe an empty range.
pub open spec fn valid_prime_search(num_bits: usize, modulo: u64, upper_bound: u64) -> bool {
    &&& 10 <= num_bits <= 62
    &&& 2 <= modulo
    &&& upper_bound <= pow2(num_bits as nat)
}

/// A prime in the requested bit-width interval and residue class.
pub open spec fn eligible_prime(p: u64, num_bits: usize, modulo: u64, upper_bound: u64) -> bool {
    &&& valid_prime_search(num_bits, modulo, upper_bound)
    &&& pow2((num_bits - 1) as nat) <= p < upper_bound
    &&& p % modulo == 1
    &&& prime(p as int)
}

/// Return the greatest eligible prime, or None exactly when none exists.
pub open spec fn prime_search_result(
    result: Option<u64>, num_bits: usize, modulo: u64, upper_bound: u64,
) -> bool {
    match result {
        Some(p) => eligible_prime(p, num_bits, modulo, upper_bound)
            && forall|q: u64| #[trigger] eligible_prime(q, num_bits, modulo, upper_bound) ==> q <= p,
        None => forall|q: u64| !#[trigger] eligible_prime(q, num_bits, modulo, upper_bound),
    }
}

/// Two distinct integers in the same residue class are at least one step apart.
proof fn residue_spacing(a: int, b: int, modulo: int)
    requires 1 < modulo, a < b, a % modulo == 1, b % modulo == 1,
    ensures a + modulo <= b,
{
    lemma_fundamental_div_mod(a, modulo);
    lemma_fundamental_div_mod(b, modulo);
    let qa = a / modulo;
    let qb = b / modulo;
    lemma_mul_is_commutative(qa, modulo);
    lemma_mul_is_commutative(qb, modulo);
    lemma_mul_strict_inequality_converse(qa, qb, modulo);
    lemma_mul_inequality(qa + 1, qb, modulo);
    assert((qa + 1) * modulo == qa * modulo + modulo) by (nonlinear_arith);
}

/// Generate the greatest `num_bits`-bit prime congruent to 1 modulo `modulo`
/// and strictly smaller than `upper_bound`.
///
/// Returns `None` if no such prime exists, or if `num_bits` is outside
/// `10..=62`, `modulo < 2`, or `upper_bound > 2^num_bits`. Empty ranges,
/// including `upper_bound == 0`, also return `None`.
#[must_use]
pub fn generate_prime(num_bits: usize, modulo: u64, upper_bound: u64) -> (result: Option<u64>)
    ensures prime_search_result(result, num_bits, modulo, upper_bound),
{
    generate_prime_with_backend::<Runtime>(num_bits, modulo, upper_bound)
}

// Only primality detection uses the existing trusted runtime contract. The
// public wrapper, validation, alignment, and entire search loop are verified.
fn generate_prime_with_backend<B: Backend>(
    num_bits: usize, modulo: u64, upper_bound: u64,
) -> (result: Option<u64>)
    ensures prime_search_result(result, num_bits, modulo, upper_bound),
{
    if !(10..=62).contains(&num_bits) || modulo < 2 {
        return None;
    }
    proof {
        lemma2_to64();
        lemma_pow2_strictly_increases(num_bits as nat, 64);
        lemma_pow2_strictly_increases((num_bits - 1) as nat, num_bits as nat);
        lemma_pow2_strictly_increases(0, (num_bits - 1) as nat);
        lemma_u64_shl_is_mul(1, num_bits as u64);
        lemma_u64_shl_is_mul(1, (num_bits - 1) as u64);
    }
    let limit = 1u64 << num_bits;
    let lower = 1u64 << (num_bits - 1);
    if upper_bound > limit || upper_bound <= lower {
        return None;
    }

    // Round down directly to residue 1 instead of scanning every intervening
    // integer. The initial candidate can be below the requested bit width.
    let start = upper_bound - 1;
    let gap = (start - 1) % modulo;
    proof {
        lemma_mod_bound(start as int - 1, modulo as int);
        lemma_fundamental_div_mod(start as int - 1, modulo as int);
        assert(gap <= start - 1) by (nonlinear_arith)
            requires start >= 1, modulo >= 2,
                gap == (start - 1) % (modulo as int);
    }
    let mut candidate = start - gap;
    proof {
        lemma_mul_is_commutative((start as int - 1) / modulo as int, modulo as int);
        lemma_fundamental_div_mod_converse_mod(
            candidate as int, modulo as int, (start as int - 1) / modulo as int, 1,
        );
        assert forall|q: u64| #[trigger] eligible_prime(q, num_bits, modulo, upper_bound)
            implies q <= candidate by {
            if candidate < q {
                residue_spacing(candidate as int, q as int, modulo as int);
            }
        }
    }

    while candidate >= lower
        invariant
            valid_prime_search(num_bits, modulo, upper_bound),
            lower == pow2((num_bits - 1) as nat), 2 <= lower,
            candidate < upper_bound,
            candidate % modulo == 1,
            forall|q: u64| #[trigger] eligible_prime(q, num_bits, modulo, upper_bound) ==> q <= candidate,
        decreases candidate,
    {
        if B::is_prime(candidate) {
            return Some(candidate);
        }
        let ghost previous = candidate;
        proof {
            if candidate < modulo {
                lemma_small_mod(candidate as nat, modulo as nat);
            }
            lemma_mod_sub_multiples_vanish(candidate as int, modulo as int);
        }
        candidate -= modulo;
        proof {
            assert forall|q: u64| #[trigger] eligible_prime(q, num_bits, modulo, upper_bound)
                implies q <= candidate by {
                if candidate < q {
                    if q < previous {
                        residue_spacing(q as int, previous as int, modulo as int);
                    }
                }
            }
        }
    }
    None
}

}

#[cfg(test)]
mod tests {
    use num_bigint::BigUint;
    use num_traits::ToPrimitive;
    use proptest::prelude::*;

    use super::{generate_prime, supports_opt};

    // Evaluate the original NFL inequality with unbounded arithmetic, including
    // zero and 64-bit inputs. This does not use the verified implementation's
    // division-based simplification.
    fn supports_opt_reference(p: u64) -> bool {
        let t = BigUint::from(1u64) << p.leading_zeros() as usize;
        let cube = &t * &t * &t;
        ((&cube + 1u64) << 64) < cube * (t + 1u64) * p
    }

    #[test]
    fn supports_opt_small_inputs() {
        for p in 0..=u16::MAX as u64 {
            assert_eq!(supports_opt(p), supports_opt_reference(p), "p={p}");
        }
    }

    #[test]
    fn supports_opt_boundaries() {
        assert!(!supports_opt(0));
        assert!(!supports_opt(1));
        assert!(!supports_opt(1 << 63));
        assert!(!supports_opt(u64::MAX));
        assert!(supports_opt((1 << 62) - 1));
        assert!(supports_opt((1 << 63) - 1));

        for bit in 0..64 {
            let power = 1u64 << bit;
            for p in [power - 1, power, power + 1] {
                assert_eq!(supports_opt(p), supports_opt_reference(p), "p={p}");
            }
        }

        // The strict inequality is an equality at this 63-bit input (s = 1).
        let equality = 3u64 << 61;
        assert!(!supports_opt(equality - 1));
        assert!(!supports_opt(equality));
        assert!(supports_opt(equality + 1));

        // Locate the first admissible integer independently for every possible
        // nonzero leading-zero count, and test its neighbors within that width.
        for bits in 1usize..=63 {
            let t = BigUint::from(1u64) << (64 - bits);
            let cube = &t * &t * &t;
            let first = (((&cube + 1u64) << 64usize) / (&cube * (&t + 1u64)))
                .to_u64()
                .unwrap()
                + 1;
            let min = 1u64 << (bits - 1);
            let max = (1u64 << bits) - 1;
            for p in [first - 1, first, first + 1] {
                if (min..=max).contains(&p) {
                    assert_eq!(supports_opt(p), p >= first, "p={p}, bits={bits}");
                    assert_eq!(supports_opt(p), supports_opt_reference(p), "p={p}");
                }
            }
        }
    }

    proptest! {
        #[test]
        fn supports_opt_matches_biguint(p in any::<u64>()) {
            prop_assert_eq!(supports_opt(p), supports_opt_reference(p));
        }
    }

    // Verifies that the same moduli as in the NFLlib library are generated.
    // <https://github.com/quarkslab/NFLlib/blob/master/include/nfl/params.hpp>
    #[test]
    fn nfl_62bit_primes() {
        let mut generated = vec![];
        let mut upper_bound = u64::MAX >> 2;
        while generated.len() != 20 {
            let p = generate_prime(62, 2 * 1048576, upper_bound);
            assert!(p.is_some());
            upper_bound = p.unwrap();
            generated.push(upper_bound);
        }
        assert_eq!(
            generated,
            vec![
                4611686018326724609,
                4611686018309947393,
                4611686018282684417,
                4611686018257518593,
                4611686018232352769,
                4611686018171535361,
                4611686018106523649,
                4611686018058289153,
                4611686018051997697,
                4611686017974403073,
                4611686017812922369,
                4611686017781465089,
                4611686017773076481,
                4611686017678704641,
                4611686017666121729,
                4611686017647247361,
                4611686017590624257,
                4611686017554972673,
                4611686017529806849,
                4611686017517223937
            ]
        )
    }

    #[test]
    fn invalid_searches_and_empty_intervals() {
        for bits in [0, 1, 9, 63, 64, usize::MAX] {
            assert_eq!(generate_prime(bits, 2, u64::MAX), None);
        }
        for bits in 10..=62 {
            let limit = 1u64 << bits;
            let lower = limit / 2;
            for modulo in [0, 1] {
                assert_eq!(generate_prime(bits, modulo, limit), None);
            }
            for bound in [0, 1, lower - 1, lower, limit + 1, u64::MAX] {
                assert_eq!(generate_prime(bits, 2, bound), None);
            }
            for modulo in [limit, 1 << 63, u64::MAX] {
                assert_eq!(generate_prime(bits, modulo, limit), None);
            }
        }
    }

    fn trial_prime(n: u64) -> bool {
        if n < 2 {
            return false;
        }
        let mut divisor = 2;
        while divisor <= n / divisor {
            if n.is_multiple_of(divisor) {
                return false;
            }
            divisor += 1;
        }
        true
    }

    #[test]
    fn prime_search_matches_exhaustive_trial_division() {
        // Independent ascending oracle: as the exclusive bound increases,
        // remember the latest eligible prime, including gaps with no solution.
        for bits in [10, 11] {
            let lower = 1u64 << (bits - 1);
            let limit = 1u64 << bits;
            let primes: Vec<_> = (lower..limit).filter(|p| trial_prime(*p)).collect();
            for modulo in (0..=64).chain([127, 256, 511, 512, 513, 1024, u64::MAX]) {
                let mut expected = None;
                for bound in 0..=limit {
                    if bound > lower {
                        let candidate = bound - 1;
                        if modulo >= 2 && candidate % modulo == 1 && primes.contains(&candidate) {
                            expected = Some(candidate);
                        }
                    }
                    assert_eq!(
                        generate_prime(bits, modulo, bound),
                        expected,
                        "bits={bits}, modulo={modulo}, bound={bound}"
                    );
                }
            }
        }
    }

    #[test]
    fn prime_search_alignment_and_strict_bounds() {
        assert_eq!(generate_prime(10, 2, 521), None);
        assert_eq!(generate_prime(10, 2, 522), Some(521));
        assert_eq!(generate_prime(10, 500, 1024), None);
        assert_eq!(generate_prime(10, 510, 1024), Some(1021));
        assert_eq!(generate_prime(10, 1020, 1024), Some(1021));
        assert_eq!(generate_prime(10, 1021, 1024), None);

        // One eligible candidate at a large distance from the upper bound.
        // Arithmetic alignment avoids a linear scan across that distance.
        let prime = 4611686018326724609;
        assert_eq!(generate_prime(62, prime - 1, 1 << 62), Some(prime));
        assert_eq!(generate_prime(62, prime - 1, prime), None);

        // The only candidate in each interval is composite (divisible by 3).
        assert_eq!(generate_prime(62, (1 << 62) - 2, 1 << 62), None);
        assert_eq!(generate_prime(62, 1 << 61, 1 << 62), None);
    }

    #[test]
    fn modulo_too_large() {
        assert!(generate_prime(10, 2048, 1 << 10).is_none());
    }

    #[test]
    fn not_found() {
        // 1033 is the smallest 11-bit prime congruent to 1 modulo 16, so looking for a
        // smaller one should fail.
        assert!(generate_prime(11, 16, 1033).is_none());
    }
}
