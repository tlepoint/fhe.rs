//! NTT eligibility and exact multiplicative order for power-of-two roots.

#![cfg_attr(verus_keep_ghost, allow(unused_imports))]

use super::native::NttOperator;
#[cfg(verus_keep_ghost)]
use crate::zq::verified::prime;
use crate::zq::{
    Modulus,
    backend::{Backend, Runtime},
};
#[cfg(verus_keep_ghost)]
use vstd::arithmetic::div_mod::{
    lemma_fundamental_div_mod, lemma_fundamental_div_mod_converse_mod, lemma_mod_bound,
    lemma_mul_mod_noop_general, lemma_small_mod,
};
#[cfg(verus_keep_ghost)]
use vstd::arithmetic::power::{
    lemma_pow_adds, lemma_pow_mod_noop, lemma_pow_multiplies, lemma1_pow, pow,
};
#[cfg(verus_keep_ghost)]
use vstd::arithmetic::power2::is_pow2;
use vstd::prelude::*;

verus! {

/// Parameters for a negacyclic transform with a primitive 2n-th root.
pub(crate) open spec fn ntt_supported(p: u64, n: usize) -> bool {
    &&& n >= 8
    &&& is_pow2(n as int)
    &&& prime(p as int)
    &&& 2 * n < p
    &&& p as int % (2 * n as int) == 1
}

/// Exact multiplicative order, rather than just an n-th root of unity.
pub(crate) open spec fn primitive_root(a: int, n: nat, p: int) -> bool {
    &&& p > 1
    &&& n > 0
    &&& pow(a, n) % p == 1
    &&& forall|k: nat| 0 < k < n ==> #[trigger] (pow(a, k) % p) != 1
}

/// Connect the usual single-set-bit test with mathematical powers of two.
proof fn power_of_two_bits(n: u64)
    ensures (n != 0 && n & n.wrapping_sub(1) == 0) == is_pow2(n as int),
    decreases n,
{
    reveal(is_pow2);
    if n > 1 {
        power_of_two_bits(n / 2);
        assert((n & n.wrapping_sub(1) == 0)
            == (n % 2 == 0 && (n / 2) & (n / 2).wrapping_sub(1) == 0))
            by (bit_vector) requires n > 1;
    } else {
        assert(n == 1 ==> n & n.wrapping_sub(1) == 0) by (bit_vector);
    }
}

#[inline]
fn power_of_two(n: u64) -> (result: bool)
    ensures result == is_pow2(n as int),
{
    proof { power_of_two_bits(n); }
    n != 0 && n & (n - 1) == 0
}

/// Whether a prime modulus supports a negacyclic NTT of size n.
/// Invalid sizes and moduli return false, including sizes whose double would
/// overflow u64. No multiplication or narrowing happens before validation.
pub(crate) fn supports_ntt(p: u64, n: usize) -> (result: bool)
    ensures result == ntt_supported(p, n),
{
    supports_ntt_with_backend::<Runtime>(p, n)
}

fn supports_ntt_with_backend<B: Backend>(p: u64, n: usize) -> (result: bool)
    ensures result == ntt_supported(p, n),
{
    if p < 2 || n < 8 || n as u128 > ((p - 1) / 2) as u128 {
        return false;
    }
    let size = n as u64;
    if !power_of_two(size) {
        return false;
    }
    p % (2 * size) == 1 && B::is_prime(p)
}

/// A multiple of a period is still a period.
proof fn period_multiple(a: int, p: int, d: nat, q: nat)
    requires p > 1, pow(a, d) % p == 1,
    ensures pow(a, d * q) % p == 1,
{
    lemma_pow_multiplies(a, d, q);
    lemma_pow_mod_noop(pow(a, d), q, p);
    lemma1_pow(q);
    lemma_small_mod(1, p as nat);
}

/// Reduce an exponent modulo a known positive period.
proof fn period_remainder(a: int, p: int, n: nat, d: nat)
    requires p > 1, d > 0, pow(a, d) % p == 1,
    ensures pow(a, n % d) % p == pow(a, n) % p,
{
    lemma_fundamental_div_mod(n as int, d as int);
    period_multiple(a, p, d, n / d);
    lemma_pow_adds(a, d * (n / d), n % d);
    lemma_mul_mod_noop_general(pow(a, d * (n / d)), pow(a, n % d), p);
    lemma_mod_bound(pow(a, n % d), p);
    lemma_small_mod((pow(a, n % d) % p) as nat, p as nat);
}

/// Euclid's algorithm produces a common period dividing both input periods.
proof fn common_period(a: int, p: int, n: nat, k: nat) -> (d: nat)
    requires p > 1, n > 0, k > 0, pow(a, n) % p == 1, pow(a, k) % p == 1,
    ensures d > 0, d <= k, n % d == 0, k % d == 0, pow(a, d) % p == 1,
    decreases k,
{
    lemma_mod_bound(n as int, k as int);
    if n % k == 0 {
        k
    } else {
        period_remainder(a, p, n, k);
        let d = common_period(a, p, k, n % k);
        lemma_fundamental_div_mod(n as int, k as int);
        lemma_mul_mod_noop_general(k as int, (n / k) as int, d as int);
        vstd::arithmetic::div_mod::lemma_add_mod_noop(
            (k * (n / k)) as int, (n % k) as int, d as int,
        );
        assert((k * (n / k)) % d == 0) by (nonlinear_arith)
            requires k % d == 0, d > 0,
                (k * (n / k)) % d == ((k % d) * ((n / k) % d)) % d;
        assert(n as int == k * (n / k) + n % k);
        assert(n % d == 0);
        d
    }
}

/// Every factor greater than one of a power of two is even.
proof fn pow2_factor_even(n: int, d: int, q: int)
    requires is_pow2(n), d > 0, q > 1, n == d * q,
    ensures q % 2 == 0,
    decreases n,
{
    reveal(is_pow2);
    assert(n >= q) by (nonlinear_arith) requires d > 0, q > 1, n == d * q;
    if q % 2 != 0 {
        lemma_mod_bound(d, 2);
        lemma_mod_bound(q, 2);
        lemma_mul_mod_noop_general(d, q, 2);
        assert(((d % 2) * (q % 2)) % 2 == d % 2) by (nonlinear_arith)
            requires 0 <= d % 2 < 2, q % 2 == 1;
        assert(d % 2 == 0);
        lemma_fundamental_div_mod(d, 2);
        lemma_fundamental_div_mod(n, 2);
        assert(n / 2 == (d / 2) * q) by (nonlinear_arith)
            requires n == d * q, n == 2 * (n / 2), d == 2 * (d / 2);
        pow2_factor_even(n / 2, d / 2, q);
    }
}

/// Any proper positive divisor of a power of two divides its half.
proof fn proper_divisor_divides_half(n: nat, d: nat)
    requires is_pow2(n as int), 0 < d < n, n % d == 0,
    ensures (n / 2) % d == 0,
{
    lemma_fundamental_div_mod(n as int, d as int);
    let q = n / d;
    assert(q > 1) by (nonlinear_arith) requires n == d * q, 0 < d < n;
    pow2_factor_even(n as int, d as int, q as int);
    lemma_fundamental_div_mod(q as int, 2);
    assert(n / 2 == d * (q / 2)) by (nonlinear_arith)
        requires n == d * q, q == 2 * (q / 2);
    lemma_fundamental_div_mod_converse_mod((n / 2) as int, d as int, (q / 2) as int, 0);
}

/// For powers of two, the half-exponent check excludes every smaller order.
proof fn power_of_two_root_criterion(a: int, n: nat, p: int)
    requires p > 1, n >= 2, is_pow2(n as int),
    ensures primitive_root(a, n, p)
        == (pow(a, n) % p == 1 && pow(a, n / 2) % p != 1),
{
    if pow(a, n) % p == 1 && pow(a, n / 2) % p != 1 {
        assert forall|k: nat| 0 < k < n implies #[trigger] (pow(a, k) % p) != 1 by {
            if pow(a, k) % p == 1 {
                let d = common_period(a, p, n, k);
                proper_divisor_divides_half(n, d);
                lemma_fundamental_div_mod((n / 2) as int, d as int);
                period_multiple(a, p, d, (n / 2) / d);
            }
        }
    }
}

impl NttOperator {
    /// Test exact order n modulo p, for power-of-two orders.
    ///
    /// Debug builds reject noncanonical a, n < 2, non-power-of-two n, or n >= p.
    /// Neither primality nor NTT eligibility is needed for this criterion.
    pub(super) fn is_primitive_root(a: u64, n: usize, p: &Modulus) -> (result: bool)
        requires a < p.value(), 2 <= n < p.value(), is_pow2(n as int),
        ensures result == primitive_root(a as int, n as nat, p.value() as int),
    {
        let modulus = p.as_u64();
        debug_assert!(a < modulus);
        debug_assert!(2 <= n && (n as u128) < modulus as u128);
        debug_assert!(power_of_two(n as u64));
        proof { power_of_two_root_criterion(a as int, n as nat, modulus as int); }
        p.pow(a, n as u64) == 1 && p.pow(a, (n / 2) as u64) != 1
    }
}

}

#[cfg(test)]
mod tests {
    use super::*;

    fn trial_prime(p: u64) -> bool {
        p >= 2
            && (2..)
                .take_while(|d| d * d <= p)
                .all(|d| !p.is_multiple_of(d))
    }

    #[test]
    fn supports_ntt_exhausts_small_parameters() {
        for p in 0..=2048 {
            let prime = trial_prime(p);
            for n in 0usize..=128 {
                let expected = prime
                    && n >= 8
                    && n.is_power_of_two()
                    && (p as u128 - 1).is_multiple_of(2 * n as u128);
                assert_eq!(supports_ntt(p, n), expected, "p={p}, n={n}");
            }
        }
    }

    #[test]
    fn supports_ntt_covers_size_and_modulus_boundaries() {
        let mut sizes = vec![0, 1, 7, 8, 9, usize::MAX];
        for bit in 1..usize::BITS {
            let n = 1usize << bit;
            sizes.extend([n - 1, n, n + 1]);
        }
        for p in [
            0,
            1,
            2,
            17,
            97,
            65537,
            4611686018326724609,
            u64::MAX - 58,
            u64::MAX,
        ] {
            let prime = fhe_util::is_prime(p);
            for &n in &sizes {
                // Wide arithmetic independently checks the congruence without
                // the implementation's upper-bound guard or narrowing cast.
                let expected = prime
                    && n >= 8
                    && n.is_power_of_two()
                    && (p as u128 - 1).is_multiple_of(2 * n as u128);
                assert_eq!(supports_ntt(p, n), expected, "p={p}, n={n}");
            }
        }
    }

    #[test]
    fn ntt_constructor_rejects_invalid_sizes() {
        let p = Modulus::new(17).unwrap();
        for n in [0, 1, 2, 4, 7, 9, 15, 16, usize::MAX, 1 << (usize::BITS - 1)] {
            assert!(NttOperator::new(&p, n).is_none());
            assert!(super::super::NttOperator::new(&p, n).is_none());
        }
        assert!(NttOperator::new(&p, 8).is_some());
    }

    // Find the first positive exponent returning one using repeated modular
    // multiplication, independently of the two exponentiation checks.
    fn exact_order_reference(a: u64, n: usize, p: u64) -> bool {
        let mut r = 1;
        for exponent in 1..=n {
            r = ((r as u128 * a as u128) % p as u128) as u64;
            if r == 1 {
                return exponent == n;
            }
        }
        false
    }

    #[test]
    fn primitive_root_exhausts_small_orders_and_moduli() {
        for p in 2..=128 {
            let modulus = Modulus::new(p).unwrap();
            for n in [2, 4, 8, 16, 32, 64] {
                if n >= p as usize {
                    continue;
                }
                for a in 0..p {
                    assert_eq!(
                        NttOperator::is_primitive_root(a, n, &modulus),
                        exact_order_reference(a, n, p),
                        "a={a}, n={n}, p={p}",
                    );
                }
            }
        }
    }

    #[test]
    fn primitive_root_covers_large_moduli_and_orders() {
        for p in [
            65537,
            2013265921,
            4611686018326724609,
            crate::zq::verified::MODULUS_LIMIT - 1,
        ] {
            let modulus = Modulus::new(p).unwrap();
            for n in [2, 4, 8, 16, 128, 1024] {
                for a in [0, 1, 2, p / 2, p - 1] {
                    assert_eq!(
                        NttOperator::is_primitive_root(a, n, &modulus),
                        exact_order_reference(a, n, p),
                        "a={a}, n={n}, p={p}",
                    );
                }
            }
            for bit in 1..usize::BITS {
                let n = 1usize << bit;
                if n as u128 >= p as u128 {
                    continue;
                }
                assert!(!NttOperator::is_primitive_root(0, n, &modulus));
                assert!(!NttOperator::is_primitive_root(1, n, &modulus));
                assert_eq!(NttOperator::is_primitive_root(p - 1, n, &modulus), n == 2);
            }
        }
    }

    #[test]
    #[cfg(debug_assertions)]
    fn primitive_root_rejects_invalid_parameters() {
        let p = Modulus::new(17).unwrap();
        for (a, n) in [
            (17, 8),
            (u64::MAX, 8),
            (1, 0),
            (1, 1),
            (1, 3),
            (1, 17),
            (1, usize::MAX),
        ] {
            assert!(std::panic::catch_unwind(|| NttOperator::is_primitive_root(a, n, &p)).is_err());
        }
    }
}
