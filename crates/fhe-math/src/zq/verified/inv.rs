//! Modular inversion by Fermat exponentiation, with a proved number-theoretic
//! foundation and an explicit primality-check boundary.

#[cfg(verus_keep_ghost)]
use super::lazy_reduce::valid_barrett;
use super::pow_mod;
#[cfg(verus_keep_ghost)]
use vstd::arithmetic::div_mod::{
    lemma_fundamental_div_mod, lemma_mod_adds, lemma_mod_bound, lemma_mod_multiples_vanish,
    lemma_mul_mod_noop_general, lemma_small_mod,
};
#[cfg(verus_keep_ghost)]
use vstd::arithmetic::power::{lemma_pow_adds, lemma_pow0, lemma_pow1, pow as int_pow};
use vstd::prelude::*;
#[cfg(verus_keep_ghost)]
use vstd::set_lib::lemma_subset_equality;

verus! {

broadcast use vstd::arithmetic::mul::group_mul_basics;

/// A prime has no positive proper divisor other than one.
pub open spec fn prime(p: int) -> bool {
    p > 1 && forall|d: int| 1 < d < p ==> #[trigger] (p % d) != 0
}

/// Euclid's recursion constructs a positive common divisor and Bezout witnesses.
proof fn bezout(x: int, y: int) -> (r: (int, int, int))
    requires x > 0, y >= 0,
    ensures
        r.0 > 0,
        x % r.0 == 0,
        y % r.0 == 0,
        r.0 == r.1 * x + r.2 * y,
    decreases y,
{
    if y == 0 {
        (x, 1, 0)
    } else {
        lemma_mod_bound(x, y);
        lemma_fundamental_div_mod(x, y);
        let (d, u, v) = bezout(y, x % y);
        let q = x / y;
        lemma_mul_mod_noop_general(y, q, d);
        assert((y % d) * q == 0) by (nonlinear_arith) requires y % d == 0;
        assert((y * q) % d == 0);
        lemma_mod_adds(y * q, x % y, d);
        assert(x == y * q + x % y);
        assert(x % d == 0);
        assert(d == v * x + (u - v * q) * y) by (nonlinear_arith)
            requires d == u * y + v * (x % y), x == y * q + x % y;
        (d, v, u - v * q)
    }
}

/// Every nonzero canonical residue has a Bezout inverse modulo a prime.
proof fn inverse_witness(a: int, p: int) -> (k: int)
    requires prime(p), 0 < a < p,
    ensures (k * a) % p == 1,
{
    let (d, u, v) = bezout(p, a);
    if d > a {
        assert(a % d == a) by (nonlinear_arith) requires 0 < a < d;
    }
    assert(d <= a);
    assert(d == 1);
    lemma_mod_multiples_vanish(u, v * a, p);
    lemma_small_mod(1, p as nat);
    assert((u * p + v * a) % p == 1);
    v
}

/// Cancelling a nonzero factor modulo a prime follows from its Bezout inverse.
proof fn cancel(a: int, b: int, c: int, p: int)
    requires prime(p), 0 < a < p, (a * b) % p == (a * c) % p,
    ensures b % p == c % p,
{
    let k = inverse_witness(a, p);
    lemma_mul_mod_noop_general(k, a * b, p);
    lemma_mul_mod_noop_general(k, a * c, p);
    lemma_mul_mod_noop_general(k * a, b, p);
    lemma_mul_mod_noop_general(k * a, c, p);
    assert(k * (a * b) == (k * a) * b) by (nonlinear_arith);
    assert(k * (a * c) == (k * a) * c) by (nonlinear_arith);
}

/// Multiplication of canonical residues is associative.
proof fn mod_associative(a: int, b: int, c: int, p: int)
    requires p > 0,
    ensures ((a * b) % p * c) % p == (a * ((b * c) % p)) % p,
{
    lemma_mul_mod_noop_general(a * b, c, p);
    lemma_mul_mod_noop_general(a, b * c, p);
    assert((a * b) * c == a * (b * c)) by (nonlinear_arith);
}

/// Products are reduced at each step; the empty product is one.
spec fn product(s: Seq<int>, p: int) -> int
    decreases s.len(),
{
    if s.len() == 0 { 1 } else { (product(s.drop_last(), p) * s.last()) % p }
}

/// Removing a factor changes the product by precisely that factor.
proof fn product_remove(s: Seq<int>, i: int, p: int)
    requires p > 1, 0 <= i < s.len(),
    ensures product(s, p) == (product(s.remove(i), p) * s[i]) % p,
    decreases s.len(),
{
    let n = s.len() as int - 1;
    if i == n {
        assert(s.remove(i) =~= s.drop_last());
    } else {
        let t = s.drop_last();
        product_remove(t, i, p);
        assert(s.remove(i).drop_last() =~= t.remove(i));
        assert(s.remove(i).last() == s.last());
        let r = product(t.remove(i), p);
        mod_associative(r, s[i], s.last(), p);
        mod_associative(r, s.last(), s[i], p);
        assert(s[i] * s.last() == s.last() * s[i]) by (nonlinear_arith);
    }
}

/// Removing one member of a duplicate-free sequence removes it from its set.
proof fn remove_member(s: Seq<int>, i: int)
    requires s.no_duplicates(), 0 <= i < s.len(),
    ensures
        s.remove(i).no_duplicates(),
        s.remove(i).to_set() == s.to_set().remove(s[i]),
{
    let t = s.remove(i);
    s.remove_ensures(i);
    assert forall|j: int, k: int| 0 <= j < t.len() && 0 <= k < t.len() && j != k
        implies #[trigger] t[j] != #[trigger] t[k] by {
        let sj = if j < i { j } else { j + 1 };
        let sk = if k < i { k } else { k + 1 };
        assert(s[sj] != s[sk]);
    }
    assert(t.to_set() =~= s.to_set().remove(s[i])) by {
        assert forall|x: int| #[trigger] t.to_set().contains(x)
            <==> s.to_set().remove(s[i]).contains(x) by {
            if t.to_set().contains(x) {
                let j = choose|j: int| 0 <= j < t.len() && t[j] == x;
                let sj = if j < i { j } else { j + 1 };
                assert(s[sj] == x && s[sj] != s[i]);
            } else if s.to_set().remove(s[i]).contains(x) {
                let j = choose|j: int| 0 <= j < s.len() && s[j] == x;
                let tj = if j < i { j } else { j - 1 };
                assert(t[tj] == x);
            }
        }
    }
}

/// Reordering distinct factors preserves their modular product.
proof fn product_permutation(s: Seq<int>, t: Seq<int>, p: int)
    requires p > 1, s.no_duplicates(), t.no_duplicates(), s.to_set() == t.to_set(),
    ensures product(s, p) == product(t, p),
    decreases s.len(),
{
    if s.len() == 0 {
        if t.len() != 0 { assert(s.to_set().contains(t[0])); }
    } else {
        let last = s.len() as int - 1;
        assert(t.to_set().contains(s[last]));
        let i = choose|i: int| 0 <= i < t.len() && t[i] == s[last];
        remove_member(s, last);
        remove_member(t, i);
        product_remove(t, i, p);
        assert(s.remove(last) =~= s.drop_last());
        product_permutation(s.drop_last(), t.remove(i), p);
    }
}

/// Scaling each of n factors multiplies their product by a^n modulo p.
proof fn product_scaled(s: Seq<int>, a: int, p: int)
    requires p > 1,
    ensures product(s.map_values(|x: int| (a * x) % p), p)
        == (int_pow(a, s.len()) * product(s, p)) % p,
    decreases s.len(),
{
    let t = s.map_values(|x: int| (a * x) % p);
    if s.len() == 0 {
        lemma_pow0(a);
        lemma_small_mod(1, p as nat);
    } else {
        let rest = s.drop_last();
        let scaled_rest = rest.map_values(|x: int| (a * x) % p);
        assert(t.drop_last() =~= scaled_rest);
        product_scaled(rest, a, p);
        lemma_pow_adds(a, rest.len(), 1);
        lemma_pow1(a);
        let power = int_pow(a, rest.len());
        let r = product(rest, p);
        let last = s.last();
        lemma_mul_mod_noop_general(power * r, a * last, p);
        lemma_mul_mod_noop_general(power * a, r * last, p);
        assert((power * r) * (a * last) == (power * a) * (r * last)) by (nonlinear_arith);
    }
}

/// A product of nonzero canonical residues is itself nonzero modulo a prime.
proof fn product_nonzero(s: Seq<int>, p: int)
    requires prime(p), forall|i: int| 0 <= i < s.len() ==> 0 < #[trigger] s[i] < p,
    ensures 0 < product(s, p) < p,
    decreases s.len(),
{
    if s.len() != 0 {
        product_nonzero(s.drop_last(), p);
        let r = product(s.drop_last(), p);
        let a = s.last();
        lemma_mod_bound(r * a, p);
        if (r * a) % p == 0 {
            lemma_small_mod(0, p as nat);
            cancel(r, a, 0, p);
            lemma_small_mod(a as nat, p as nat);
        }
    }
}

/// Fermat's little theorem from the ordinary divisor-based definition of prime.
proof fn fermat(a: int, p: int)
    requires prime(p), 0 < a < p,
    ensures int_pow(a, (p - 1) as nat) % p == 1,
{
    let s = Seq::new((p - 1) as nat, |i: int| i + 1);
    let t = s.map_values(|x: int| (a * x) % p);
    assert(s.no_duplicates());
    assert forall|i: int| 0 <= i < t.len() implies 0 < #[trigger] t[i] < p by {
        lemma_mod_bound(a * s[i], p);
        if t[i] == 0 {
            lemma_small_mod(0, p as nat);
            cancel(a, s[i], 0, p);
            lemma_small_mod(s[i] as nat, p as nat);
        }
    }
    assert forall|i: int, j: int| 0 <= i < t.len() && 0 <= j < t.len() && i != j
        implies #[trigger] t[i] != #[trigger] t[j] by {
        if t[i] == t[j] {
            cancel(a, s[i], s[j], p);
            lemma_small_mod(s[i] as nat, p as nat);
            lemma_small_mod(s[j] as nat, p as nat);
        }
    }
    assert(t.to_set().subset_of(s.to_set())) by {
        assert forall|x: int| #[trigger] t.to_set().contains(x)
            implies s.to_set().contains(x) by {
            let i = choose|i: int| 0 <= i < t.len() && t[i] == x;
            assert(s[x - 1] == x);
        }
    }
    s.unique_seq_to_set();
    t.unique_seq_to_set();
    lemma_subset_equality(t.to_set(), s.to_set());
    product_permutation(s, t, p);
    product_scaled(s, a, p);
    product_nonzero(s, p);
    let r = product(s, p);
    let power = int_pow(a, (p - 1) as nat);
    lemma_small_mod(r as nat, p as nat);
    assert(r * power == power * r) by (nonlinear_arith);
    cancel(r, power, 1, p);
    lemma_small_mod(1, p as nat);
}

/// Invert a canonical residue using the original Fermat exponentiation.
///
/// `p_is_prime` comes from the external primality check. Its equivalence to
/// mathematical primality is a precondition, not a theorem proved here.
#[inline]
pub(crate) fn inv_mod(a: u64, p: u64, hi: u64, lo: u64, p_is_prime: bool) -> (result: Option<u64>)
    requires
        valid_barrett(p, hi, lo),
        p_is_prime == prime(p as int),
        p_is_prime ==> a < p,
    ensures
        result.is_none() <==> !p_is_prime || a == 0,
        result matches Some(r) ==> 0 < r < p
            && (a as int * r as int) % p as int == 1
            && r as int == int_pow(a as int, (p - 2) as nat) % p as int,
{
    if !p_is_prime || a == 0 {
        None
    } else {
        let r = pow_mod(a, p - 2, p, hi, lo);
        proof {
            fermat(a as int, p as int);
            lemma_pow_adds(a as int, (p - 2) as nat, 1);
            lemma_pow1(a as int);
            let power = int_pow(a as int, (p - 2) as nat);
            lemma_mul_mod_noop_general(a as int, power, p as int);
            assert(a as int * power == power * a as int) by (nonlinear_arith);
            lemma_small_mod(0, p as nat);
        }
        Some(r)
    }
}

}
