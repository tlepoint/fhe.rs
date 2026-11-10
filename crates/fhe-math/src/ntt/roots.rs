//! Existence of primitive 2n-th roots of unity and the search that finds one.

#![cfg_attr(verus_keep_ghost, allow(unused_imports))]

use super::native::NttOperator;
use super::predicates::supports_ntt;
#[cfg(verus_keep_ghost)]
use super::predicates::{ntt_supported, power_of_two_root_criterion, primitive_root};
use crate::zq::Modulus;
#[cfg(verus_keep_ghost)]
use crate::zq::verified::{cancel, fermat, prime};
use rand::{RngExt, SeedableRng};
use rand_chacha::ChaCha8Rng;
#[cfg(verus_keep_ghost)]
use vstd::arithmetic::div_mod::{
    lemma_fundamental_div_mod, lemma_fundamental_div_mod_converse, lemma_mod_add_multiples_vanish,
    lemma_mul_mod_noop_left, lemma_small_mod, lemma_sub_mod_noop,
};
#[cfg(verus_keep_ghost)]
use vstd::arithmetic::power::{
    lemma_pow_adds, lemma_pow_mod_noop, lemma_pow_multiplies, lemma_pow0, lemma_pow1, lemma1_pow,
    pow,
};
#[cfg(verus_keep_ghost)]
use vstd::arithmetic::power2::is_pow2;
use vstd::prelude::*;

/// The candidate stream of the original search: ChaCha8 seeded with zero.
/// Retaining it keeps the selected roots, and hence the transform, unchanged.
pub(super) struct Seeded(ChaCha8Rng);

verus! {

/// Number of seeded pseudorandom candidates tried before the exhaustive scan.
const SEEDED_ATTEMPTS: usize = 100;

/// Horner evaluation of a polynomial, leading coefficient first.
spec fn eval(c: Seq<int>, x: int) -> int
    decreases c.len(),
{
    if c.len() == 0 { 0 } else { eval(c.drop_last(), x) * x + c.last() }
}

/// Synthetic division by X - r: the quotient's coefficients are the partial
/// Horner values at r.
spec fn quotient(c: Seq<int>, r: int) -> Seq<int> {
    Seq::new((c.len() - 1) as nat, |k: int| eval(c.take(k + 1), r))
}

/// The factor theorem over the integers: c(x) = (x - r) q(x) + c(r).
proof fn factor(c: Seq<int>, r: int, x: int)
    requires c.len() >= 1,
    ensures eval(c, x) == (x - r) * eval(quotient(c, r), x) + eval(c, r),
    decreases c.len(),
{
    let q = quotient(c, r);
    let d = c.drop_last();
    if c.len() == 1 {
        assert(eval(d, x) == 0 && eval(d, r) == 0 && eval(q, x) == 0);
    } else {
        factor(d, r, x);
        assert(q.drop_last() =~= quotient(d, r)) by {
            assert forall|k: int| 0 <= k < d.len() - 1 implies #[trigger] c.take(k + 1) == d.take(k + 1) by {
                assert(c.take(k + 1) =~= d.take(k + 1));
            }
        }
        assert(c.take(c.len() - 1) =~= d);
        assert(q.last() == eval(d, r));
        let e = eval(quotient(d, r), x);
        assert((eval(d, x)) * x + c.last()
            == (x - r) * (e * x + eval(d, r)) + (eval(d, r) * r + c.last())) by (nonlinear_arith)
            requires eval(d, x) == (x - r) * e + eval(d, r);
    }
}

/// A nonzero factor of absolute value below a prime p can be cancelled from a
/// product divisible by p.
proof fn zero_product(a: int, b: int, p: int)
    requires prime(p), -p < a < p, a != 0, (a * b) % p == 0,
    ensures b % p == 0,
{
    let reduced = a % p;
    if a > 0 {
        lemma_small_mod(a as nat, p as nat);
    } else {
        lemma_mod_add_multiples_vanish(a, p);
        lemma_small_mod((p + a) as nat, p as nat);
    }
    lemma_mul_mod_noop_left(a, b, p);
    assert(reduced * 0 == 0);
    cancel(reduced, b, 0, p);
    lemma_small_mod(0, p as nat);
}

/// Lagrange's bound on polynomial roots modulo a prime: a polynomial with
/// more distinct roots than its degree has a leading coefficient divisible by p.
proof fn leading_coefficient_vanishes(c: Seq<int>, roots: Seq<int>, p: int)
    requires
        prime(p),
        c.len() >= 1,
        roots.len() == c.len(),
        forall|i: int| 0 <= i < roots.len() ==> 0 <= #[trigger] roots[i] < p,
        forall|i: int, j: int| 0 <= i < j < roots.len() ==> roots[i] != roots[j],
        forall|i: int| 0 <= i < roots.len() ==> eval(c, #[trigger] roots[i]) % p == 0,
    ensures c[0] % p == 0,
    decreases c.len(),
{
    if c.len() == 1 {
        assert(eval(c.drop_last(), roots[0]) == 0);
        assert(eval(c, roots[0]) == c[0]);
    } else {
        let r = roots.last();
        let q = quotient(c, r);
        let rest = roots.drop_last();
        assert forall|i: int| 0 <= i < rest.len() implies eval(q, #[trigger] rest[i]) % p == 0 by {
            let x = roots[i];
            factor(c, r, x);
            lemma_sub_mod_noop(eval(c, x), eval(c, r), p);
            lemma_small_mod(0, p as nat);
            zero_product(x - r, eval(q, x), p);
        }
        let lead = c.take(1);
        assert(eval(lead.drop_last(), r) == 0 && lead.last() == c[0]);
        assert(eval(lead, r) == c[0]);
        leading_coefficient_vanishes(q, rest, p);
    }
}

/// The coefficients of X^m - 1, and of its leading monomial alone.
spec fn power_minus_one(m: nat) -> Seq<int> {
    Seq::new(m + 1, |k: int| if k == 0 { 1 } else if k == m { -1 } else { 0 })
}

spec fn monomial(m: nat) -> Seq<int> {
    Seq::new(m + 1, |k: int| if k == 0 { 1int } else { 0int })
}

proof fn eval_monomial(m: nat, x: int)
    ensures eval(monomial(m), x) == pow(x, m),
    decreases m,
{
    let c = monomial(m);
    if m == 0 {
        lemma_pow0(x);
        assert(eval(c.drop_last(), x) == 0);
    } else {
        assert(c.drop_last() =~= monomial((m - 1) as nat));
        eval_monomial((m - 1) as nat, x);
        lemma_pow_adds(x, (m - 1) as nat, 1);
        lemma_pow1(x);
    }
}

proof fn eval_power_minus_one(m: nat, x: int)
    requires m >= 1,
    ensures eval(power_minus_one(m), x) == pow(x, m) - 1,
{
    assert(power_minus_one(m).drop_last() =~= monomial((m - 1) as nat));
    eval_monomial((m - 1) as nat, x);
    lemma_pow_adds(x, (m - 1) as nat, 1);
    lemma_pow1(x);
}

/// X^m - 1 has at most m roots modulo a prime, so one of 1, ..., m + 1 is not
/// an m-th root of unity.
proof fn not_a_root_of_unity(p: int, m: nat) -> (g: int)
    requires prime(p), 1 <= m, m + 1 < p,
    ensures 1 <= g < p, pow(g, m) % p != 1,
{
    if forall|g: int| 1 <= g <= m + 1 ==> #[trigger] pow(g, m) % p == 1 {
        let c = power_minus_one(m);
        let roots = Seq::new(m + 1, |i: int| i + 1);
        assert forall|i: int| 0 <= i < roots.len() implies eval(c, #[trigger] roots[i]) % p == 0 by {
            let g = roots[i];
            eval_power_minus_one(m, g);
            assert(pow(g, m) % p == 1);
            lemma_small_mod(1, p as nat);
            lemma_small_mod(0, p as nat);
            lemma_sub_mod_noop(pow(g, m), 1, p);
        }
        leading_coefficient_vanishes(c, roots, p);
        lemma_small_mod(1, p as nat);
        assert(false);
    }
    choose|g: int| 1 <= g <= m + 1 && pow(g, m) % p != 1
}

/// Exponentiating a canonical power: (g^e mod p)^k = g^(e k) modulo p.
proof fn power_of_power(g: int, e: nat, k: nat, p: int)
    requires p > 1,
    ensures pow(pow(g, e) % p, k) % p == pow(g, e * k) % p,
{
    lemma_pow_mod_noop(pow(g, e), k, p);
    lemma_pow_multiplies(g, e, k);
}

/// Arithmetic facts about the cofactor lambda = (p - 1) / (2n).
proof fn search_parameters(p: u64, n: usize)
    requires ntt_supported(p, n),
    ensures
        is_pow2(2 * n as int),
        p - 1 == 2 * n * ((p - 1) / (2 * n)),
        1 <= (p - 1) / (2 * n) < p,
{
    let order = 2 * n as int;
    reveal(is_pow2);
    assert(order / 2 == n as int && order % 2 == 0);
    lemma_fundamental_div_mod(p as int, order);
    lemma_fundamental_div_mod_converse(p as int - 1, order, p as int / order, 0);
    let lambda = (p as int - 1) / order;
    assert(lambda >= 1 && lambda < p) by (nonlinear_arith)
        requires p - 1 == order * lambda, 0 < order < p;
}

/// The candidate g^lambda has order dividing 2n by Fermat's little theorem.
/// It has exact order 2n unless g is a ((p - 1) / 2)-th root of unity.
proof fn rejected_candidate(g: int, n: nat, lambda: nat, p: int)
    requires
        prime(p),
        0 < g < p,
        is_pow2(2 * n as int),
        n >= 1,
        p - 1 == 2 * n * lambda,
        !primitive_root(pow(g, lambda) % p, 2 * n, p),
    ensures pow(g, n * lambda) % p == 1,
{
    let root = pow(g, lambda) % p;
    fermat(g, p);
    power_of_power(g, lambda, 2 * n, p);
    power_of_power(g, lambda, n, p);
    assert(lambda * (2 * n) == p - 1 && lambda * n == n * lambda) by (nonlinear_arith)
        requires p - 1 == 2 * n * lambda;
    power_of_two_root_criterion(root, 2 * n, p);
    assert((2 * n) / 2 == n);
}

/// A stream of candidate residues for the first phase of the root search.
///
/// The contract deliberately promises nothing about the drawn values: the
/// search validates every candidate, so its result does not depend on them.
pub(super) trait Candidates: Sized {
    /// Start the stream.
    fn start() -> Self;

    /// Draw the next candidate for the modulus p.
    fn draw(&mut self, p: u64) -> u64
        requires 2 <= p;
}

#[cfg(verus_keep_ghost)]
#[verifier::external_type_specification]
#[verifier::external_body]
pub(super) struct SeededSpec(Seeded);

#[verifier::external]
impl Candidates for Seeded {
    fn start() -> Self { Self(SeedableRng::seed_from_u64(0)) }

    fn draw(&mut self, p: u64) -> u64 { self.0.random_range(0..p) }
}

impl NttOperator {
    /// Returns a primitive 2n-th root of unity modulo p.
    ///
    /// The seeded pseudorandom candidates are tried first; an exhaustive scan,
    /// which provably succeeds, follows if none of them qualifies.
    pub(super) fn primitive_root(n: usize, p: &Modulus) -> (root: u64)
        requires ntt_supported(p.value(), n),
        ensures
            root < p.value(),
            primitive_root(root as int, 2 * n as nat, p.value() as int),
    {
        primitive_root_from::<Seeded>(n, p)
    }
}

/// Search for a primitive 2n-th root of unity among the powers g^lambda,
/// where lambda = (p - 1) / (2n).
pub(super) fn primitive_root_from<C: Candidates>(n: usize, p: &Modulus) -> (root: u64)
    requires ntt_supported(p.value(), n),
    ensures
        root < p.value(),
        primitive_root(root as int, 2 * n as nat, p.value() as int),
{
    let modulus = p.as_u64();
    debug_assert!(supports_ntt(modulus, n));

    let order = 2 * (n as u64);
    let lambda = (modulus - 1) / order;
    proof { search_parameters(modulus, n); }

    let mut candidates = C::start();
    for _attempt in 0..SEEDED_ATTEMPTS
        invariant
            modulus == p.value(),
            order == 2 * n,
            2 <= order < modulus,
            is_pow2(order as int),
            lambda < modulus,
    {
        let candidate = candidates.draw(modulus);
        // Nothing is assumed about the stream, so check the range here.
        if candidate < modulus {
            let root = p.pow(candidate, lambda);
            if NttOperator::is_primitive_root(root, order, p) {
                return root;
            }
        }
    }

    // Every g below p is tried in turn. A rejected g is a ((p - 1) / 2)-th
    // root of unity, and some g is not, so the scan stops at or before it.
    let ghost half = n as nat * lambda as nat;
    let ghost witness: int;
    proof {
        assert(1 <= half && 2 * half == modulus - 1) by (nonlinear_arith)
            requires modulus - 1 == 2 * n * lambda, half == n * lambda, n >= 1, lambda >= 1;
        witness = not_a_root_of_unity(modulus as int, half);
        lemma1_pow(half);
        lemma_small_mod(1, modulus as nat);
    }
    let mut generator: u64 = 2;
    loop
        invariant
            ntt_supported(modulus, n),
            modulus == p.value(),
            order == 2 * n,
            2 <= order < modulus,
            is_pow2(order as int),
            lambda < modulus,
            modulus - 1 == 2 * n * lambda,
            half == n as nat * lambda as nat,
            2 <= generator <= witness < modulus,
            pow(witness, half) % (modulus as int) != 1,
        decreases witness - generator,
    {
        let root = p.pow(generator, lambda);
        if NttOperator::is_primitive_root(root, order, p) {
            return root;
        }
        proof {
            rejected_candidate(generator as int, n as nat, lambda as nat, modulus as int);
        }
        generator += 1;
    }
}

}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::zq::primes::generate_prime;

    // Independent square-and-multiply over u128.
    fn pow_reference(a: u64, mut e: u64, p: u64) -> u64 {
        let (mut r, mut b) = (1u128, u128::from(a));
        while e > 0 {
            if e & 1 == 1 {
                r = r * b % u128::from(p);
            }
            b = b * b % u128::from(p);
            e >>= 1;
        }
        r as u64
    }

    // Modulo a prime, -1 is the only element of order two, so a has order 2n
    // exactly when a^n = -1.
    fn has_order_2n(a: u64, n: usize, p: u64) -> bool {
        pow_reference(a, n as u64, p) == p - 1
    }

    // The search as it was before verification, with an independent order test.
    fn original_search(n: usize, p: u64) -> Option<u64> {
        let lambda = (p - 1) / (2 * n as u64);
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        (0..100)
            .map(|_| pow_reference(rng.random_range(0..p), lambda, p))
            .find(|&root| has_order_2n(root, n, p))
    }

    fn parameters() -> Vec<(u64, usize)> {
        let mut parameters = vec![(17, 8), (97, 8), (97, 16), (193, 32), (1153, 64)];
        for log_n in 3..=15 {
            let n = 1usize << log_n;
            parameters.extend([(65537, n), (2013265921, n), (4611686018326724609, n)]);
            for bits in [20, 31, 45, 58, 62] {
                let p = generate_prime(bits, 2 * n as u64, 1 << bits).unwrap();
                parameters.push((p, n));
            }
        }
        parameters
    }

    #[test]
    fn seeded_search_selects_the_original_roots() {
        for (p, n) in parameters() {
            let modulus = Modulus::new(p).unwrap();
            let root = NttOperator::primitive_root(n, &modulus);
            assert!(has_order_2n(root, n, p), "p={p}, n={n}");
            assert_eq!(Some(root), original_search(n, p), "p={p}, n={n}");
        }
    }

    // Candidate streams on which the first phase never succeeds: zero and one
    // are not primitive roots, and out-of-range draws are skipped.
    struct Fixed<const VALUE: u64>;
    impl<const VALUE: u64> Candidates for Fixed<VALUE> {
        fn start() -> Self {
            Self
        }

        fn draw(&mut self, _p: u64) -> u64 {
            VALUE
        }
    }

    // The scan returns g^lambda for the least g >= 2 that is not a
    // ((p - 1) / 2)-th root of unity.
    fn scan_reference(n: usize, p: u64) -> u64 {
        let g = (2..p)
            .find(|&g| pow_reference(g, (p - 1) / 2, p) != 1)
            .unwrap();
        pow_reference(g, (p - 1) / (2 * n as u64), p)
    }

    #[test]
    fn exhaustive_scan_follows_rejected_candidates() {
        for (p, n) in parameters() {
            let modulus = Modulus::new(p).unwrap();
            let expected = scan_reference(n, p);
            assert!(has_order_2n(expected, n, p), "p={p}, n={n}");
            assert_eq!(primitive_root_from::<Fixed<0>>(n, &modulus), expected);
            assert_eq!(primitive_root_from::<Fixed<1>>(n, &modulus), expected);
            assert_eq!(
                primitive_root_from::<Fixed<{ u64::MAX }>>(n, &modulus),
                expected
            );
        }
    }

    #[test]
    fn every_small_supported_modulus_has_a_root() {
        let mut supported = 0;
        for p in 17..=20000 {
            for n in [8, 16, 32, 64, 128, 256, 512, 1024] {
                if !supports_ntt(p, n) {
                    continue;
                }
                supported += 1;
                let modulus = Modulus::new(p).unwrap();
                for root in [
                    NttOperator::primitive_root(n, &modulus),
                    primitive_root_from::<Fixed<0>>(n, &modulus),
                ] {
                    assert!(root < p && has_order_2n(root, n, p), "p={p}, n={n}");
                    // No smaller positive exponent returns one.
                    let mut power = 1;
                    for exponent in 1..2 * n as u64 {
                        power = power * root % p;
                        assert_ne!(power, 1, "p={p}, n={n}, exponent={exponent}");
                    }
                }
            }
        }
        assert!(supported > 400);
    }

    #[test]
    #[cfg(debug_assertions)]
    fn search_rejects_unsupported_parameters() {
        for (p, n) in [(17, 16), (19, 8), (65537, 12), (65535, 8)] {
            let modulus = Modulus::new(p).unwrap();
            assert!(std::panic::catch_unwind(|| NttOperator::primitive_root(n, &modulus)).is_err());
        }
    }
}
