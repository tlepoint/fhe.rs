//! Verified construction of the native NTT operator: the root of unity, both
//! twiddle tables in bit-reversed order, their Shoup words, and the inverse of
//! the size. The transform loops are a separate proof obligation.

#![cfg_attr(verus_keep_ghost, allow(unused_imports))]

use super::NttOperator;
#[cfg(verus_keep_ghost)]
use crate::ntt::predicates::{ntt_supported, primitive_root};
use crate::ntt::supports_ntt;
use crate::zq::Modulus;
#[cfg(verus_keep_ghost)]
use crate::zq::verified::valid_shoup;
#[cfg(verus_keep_ghost)]
use vstd::arithmetic::div_mod::{lemma_mul_mod_noop_left, lemma_small_mod};
#[cfg(verus_keep_ghost)]
use vstd::arithmetic::power::{lemma_pow_adds, lemma_pow0, lemma_pow1, lemma0_pow, pow};
#[cfg(verus_keep_ghost)]
use vstd::arithmetic::power2::is_pow2;
use vstd::prelude::*;

verus! {

/// Reversal of the log2(size) low bits of i, for a power-of-two size: the
/// lowest bit of i becomes the highest bit of the result.
pub(crate) open spec fn bit_reverse(i: nat, size: nat) -> nat
    decreases size,
{
    if size <= 1 {
        0
    } else {
        (if i % 2 == 1 { size / 2 } else { 0 }) + bit_reverse(i / 2, size / 2)
    }
}

proof fn bit_reverse_zero(size: nat)
    ensures bit_reverse(0, size) == 0,
    decreases size,
{
    if size > 1 {
        bit_reverse_zero(size / 2);
    }
}

/// The highest bit is the reversal of the lowest one.
proof fn bit_reverse_half(size: nat)
    requires is_pow2(size as int), size >= 2,
    ensures bit_reverse(size / 2, size) == 1,
    decreases size,
{
    reveal(is_pow2);
    if size == 2 {
        assert(bit_reverse(0, 1) == 0);
    } else {
        assert(is_pow2(size as int / 2));
        bit_reverse_half(size / 2);
    }
}

/// For a power-of-two size, the reversal is again below the size.
proof fn bit_reverse_bound(i: nat, size: nat)
    requires is_pow2(size as int),
    ensures bit_reverse(i, size) < size,
    decreases size,
{
    reveal(is_pow2);
    if size > 1 {
        bit_reverse_bound(i / 2, size / 2);
    }
}

/// An index below half the size has a zero top bit, so reversing it at the
/// full width doubles its reversal at one bit less.
proof fn bit_reverse_doubles(i: nat, size: nat)
    requires is_pow2(size as int), size >= 2, i < size / 2,
    ensures bit_reverse(i, size) == 2 * bit_reverse(i, size / 2),
    decreases size,
{
    reveal(is_pow2);
    if size > 2 {
        assert(is_pow2(size as int / 2));
        bit_reverse_doubles(i / 2, size / 2);
    }
}

/// The bit-reversal permutation of 0..size, filled in index order: the
/// reversal of i is the halved reversal of i / 2, with the lowest bit of i
/// placed on top.
fn bit_reversal(size: usize) -> (r: Vec<usize>)
    requires is_pow2(size as int),
    ensures
        r@.len() == size,
        forall|i: int| 0 <= i < size ==>
            #[trigger] r@[i] == bit_reverse(i as nat, size as nat) && r@[i] < size,
{
    proof {
        reveal(is_pow2);
        bit_reverse_zero(size as nat);
    }
    let half = size / 2;
    let mut table = Vec::with_capacity(size);
    table.push(0);
    while table.len() < size
        invariant
            is_pow2(size as int),
            half == size / 2,
            1 <= table@.len() <= size,
            forall|k: int| 0 <= k < table@.len() ==>
                #[trigger] table@[k] == bit_reverse(k as nat, size as nat) && table@[k] < size,
        decreases size - table@.len(),
    {
        let i = table.len();
        proof {
            reveal(is_pow2);
            assert(size % 2 == 0 && is_pow2(half as int));
            bit_reverse_doubles((i / 2) as nat, size as nat);
            bit_reverse_bound((i / 2) as nat, half as nat);
        }
        let low = table[i / 2] / 2;
        table.push(if i % 2 == 1 { low + half } else { low });
    }
    table
}

/// Consecutive powers: entry j is base^(j + offset) modulo p.
pub(crate) open spec fn power_table(table: Seq<u64>, base: int, offset: nat, p: u64) -> bool {
    forall|j: int| 0 <= j < table.len() ==>
        #[trigger] table[j] as int == pow(base, j as nat + offset) % (p as int)
}

/// A twiddle table: entry i is root^(bit_reverse(i) + offset) modulo p, in
/// canonical form and paired with its exact Shoup word.
pub(crate) open spec fn twiddle_table(
    table: Seq<u64>, shoup: Seq<u64>, root: int, offset: nat, p: u64,
) -> bool {
    &&& table.len() == shoup.len()
    &&& forall|i: int| 0 <= i < table.len() ==>
        #[trigger] table[i] as int
            == pow(root, bit_reverse(i as nat, table.len()) + offset) % (p as int)
        && valid_shoup(table[i], p, shoup[i])
}

/// Multiplying a canonical power by its base yields the next power.
proof fn next_power(base: int, e: nat, p: int)
    requires p > 0,
    ensures ((pow(base, e) % p) * base) % p == pow(base, e + 1) % p,
{
    lemma_mul_mod_noop_left(pow(base, e), base, p);
    lemma_pow_adds(base, e, 1);
    lemma_pow1(base);
}

/// The first `size` powers of `base`, starting from `first = base^offset`.
fn consecutive_powers(
    p: &Modulus, first: u64, base: u64, size: usize, Ghost(offset): Ghost<nat>,
) -> (r: Vec<u64>)
    requires
        base < p.value(),
        first < p.value(),
        first as int == pow(base as int, offset) % (p.value() as int),
    ensures
        r@.len() == size,
        power_table(r@, base as int, offset, p.value()),
        forall|j: int| 0 <= j < r@.len() ==> #[trigger] r@[j] < p.value(),
{
    let mut table = Vec::with_capacity(size);
    let mut power = first;
    while table.len() < size
        invariant
            table@.len() <= size,
            base < p.value(),
            power < p.value(),
            power as int == pow(base as int, table@.len() + offset) % (p.value() as int),
            power_table(table@, base as int, offset, p.value()),
            forall|j: int| 0 <= j < table@.len() ==> #[trigger] table@[j] < p.value(),
        decreases size - table@.len(),
    {
        proof { next_power(base as int, table@.len() + offset, p.value() as int); }
        table.push(power);
        power = p.mul(power, base);
    }
    table
}

/// Reorder two tables of equal power-of-two length into bit-reversed order.
fn bit_reversed(a: &[u64], b: &[u64]) -> (r: (Vec<u64>, Vec<u64>))
    requires is_pow2(a@.len() as int), a@.len() == b@.len(),
    ensures
        r.0@.len() == a@.len(),
        r.1@.len() == a@.len(),
        forall|i: int| 0 <= i < a@.len() ==> #[trigger] bit_reverse(i as nat, a@.len()) < a@.len(),
        forall|i: int| 0 <= i < a@.len() ==>
            #[trigger] r.0@[i] == a@[bit_reverse(i as nat, a@.len()) as int],
        forall|i: int| 0 <= i < a@.len() ==>
            #[trigger] r.1@[i] == b@[bit_reverse(i as nat, a@.len()) as int],
{
    let size = a.len();
    let reversal = bit_reversal(size);
    let mut left = Vec::with_capacity(size);
    let mut right = Vec::with_capacity(size);
    let mut i = 0;
    while i < size
        invariant
            size == a@.len(),
            size == b@.len(),
            reversal@.len() == size,
            forall|k: int| 0 <= k < size ==>
                #[trigger] reversal@[k] == bit_reverse(k as nat, size as nat) && reversal@[k] < size,
            i <= size,
            left@.len() == i,
            right@.len() == i,
            forall|k: int| 0 <= k < i ==> #[trigger] bit_reverse(k as nat, size as nat) < size,
            forall|k: int| 0 <= k < i ==>
                #[trigger] left@[k] == a@[bit_reverse(k as nat, size as nat) as int],
            forall|k: int| 0 <= k < i ==>
                #[trigger] right@[k] == b@[bit_reverse(k as nat, size as nat) as int],
        decreases size - i,
    {
        let j = reversal[i];
        left.push(a[j]);
        right.push(b[j]);
        i += 1;
    }
    (left, right)
}

impl NttOperator {
    /// Create an NTT operator given a modulus for a specific size.
    ///
    /// Returns None if the size is not a power of two at least 8, or the
    /// modulus does not support the NTT for this size. Every supported pair
    /// yields an operator.
    #[must_use]
    pub fn new(p: &Modulus, size: usize) -> (result: Option<Self>)
        ensures
            result.is_some() == ntt_supported(p.value(), size),
            result matches Some(op) ==> op.valid() && op.spec_modulus() == p.value()
                && op.spec_size() == size,
    {
        let modulus = p.as_u64();
        if !supports_ntt(modulus, size) {
            return None;
        }
        let size_inv = p.inv(size as u64)?;

        let omega = Self::primitive_root(size, p);
        proof {
            // Zero is not a root of unity, so omega is invertible.
            lemma0_pow(2 * size as nat);
            lemma_small_mod(0, modulus as nat);
        }
        let omega_inv = p.inv(omega)?;

        proof {
            lemma_pow0(omega as int);
            lemma_pow1(omega as int);
            lemma_pow1(omega_inv as int);
            lemma_small_mod(1, modulus as nat);
            lemma_small_mod(omega as nat, modulus as nat);
            lemma_small_mod(omega_inv as nat, modulus as nat);
        }
        let powers = consecutive_powers(p, 1, omega, size, Ghost(0));
        let powers_inv = consecutive_powers(p, omega_inv, omega_inv, size, Ghost(1));
        let (omegas, zetas_inv) = bit_reversed(powers.as_slice(), powers_inv.as_slice());

        let omegas_shoup = p.shoup_vec(omegas.as_slice());
        let zetas_inv_shoup = p.shoup_vec(zetas_inv.as_slice());

        let operator = Self {
            p: p.clone(),
            p_twice: modulus * 2,
            size,
            omegas,
            omegas_shoup,
            zetas_inv,
            zetas_inv_shoup,
            size_inv,
            size_inv_shoup: p.shoup(size_inv),
        };
        proof {
            bit_reverse_half(size as nat);
            bit_reverse_zero(size as nat);
            assert(operator.root() == omega);
            assert(operator.root_inv() == omega_inv);
        }
        Some(operator)
    }
}

}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::zq::primes::generate_prime;
    use itertools::Itertools;
    use std::iter::successors;

    #[test]
    fn bit_reversal_matches_the_word_reversal() {
        assert_eq!(bit_reversal(1), [0]);
        for log_size in 1..=20 {
            let size = 1usize << log_size;
            let expected = (0..size)
                .map(|i| i.reverse_bits() >> (size.leading_zeros() + 1))
                .collect_vec();
            assert_eq!(bit_reversal(size), expected, "size={size}");
        }
    }

    fn parameters() -> Vec<(u64, usize)> {
        let mut parameters = vec![(17, 8), (97, 8), (97, 16), (193, 32), (1153, 64)];
        for log_n in 3..=14 {
            let n = 1usize << log_n;
            parameters.extend([(65537, n), (2013265921, n), (4611686018326724609, n)]);
            for bits in [20, 45, 62] {
                let p = generate_prime(bits, 2 * n as u64, 1 << bits).unwrap();
                parameters.push((p, n));
            }
        }
        parameters
    }

    // The tables as they were built before verification, from the same root.
    #[test]
    fn constructor_matches_the_original_tables() {
        for (modulus, size) in parameters() {
            let p = Modulus::new(modulus).unwrap();
            let op = NttOperator::new(&p, size).unwrap();

            let omega = NttOperator::primitive_root(size, &p);
            let omega_inv = p.inv(omega).unwrap();
            let powers = successors(Some(1u64), |n| Some(p.mul(*n, omega)))
                .take(size)
                .collect_vec();
            let powers_inv = successors(Some(omega_inv), |n| Some(p.mul(*n, omega_inv)))
                .take(size)
                .collect_vec();
            let (omegas, zetas_inv): (Vec<u64>, Vec<u64>) = (0..size)
                .map(|i| {
                    let j = i.reverse_bits() >> (size.leading_zeros() + 1);
                    (powers[j], powers_inv[j])
                })
                .unzip();
            let size_inv = p.inv(size as u64).unwrap();

            assert_eq!(op.p, p);
            assert_eq!(op.p_twice, 2 * modulus);
            assert_eq!(op.size, size);
            assert_eq!(op.omegas_shoup, p.shoup_vec(&omegas));
            assert_eq!(op.zetas_inv_shoup, p.shoup_vec(&zetas_inv));
            assert_eq!(op.omegas, omegas);
            assert_eq!(op.zetas_inv, zetas_inv);
            assert_eq!(op.size_inv, size_inv);
            assert_eq!(op.size_inv_shoup, p.shoup(size_inv));
        }
    }

    fn mul_reference(a: u64, b: u64, p: u64) -> u64 {
        (u128::from(a) * u128::from(b) % u128::from(p)) as u64
    }

    fn shoup_reference(a: u64, p: u64) -> u64 {
        ((u128::from(a) << 64) / u128::from(p)) as u64
    }

    // Check the `valid` predicate directly, with u128 arithmetic and a
    // bit-at-a-time reversal as independent oracles.
    #[test]
    fn constructor_establishes_its_contract() {
        for (p, n) in parameters() {
            let op = NttOperator::new(&Modulus::new(p).unwrap(), n).unwrap();
            assert_eq!((op.p.as_u64(), op.size, op.p_twice), (p, n, 2 * p));
            for table in [
                &op.omegas,
                &op.omegas_shoup,
                &op.zetas_inv,
                &op.zetas_inv_shoup,
            ] {
                assert_eq!(table.len(), n);
            }

            let root = op.omegas[n / 2];
            let root_inv = op.zetas_inv[0];
            assert_eq!(mul_reference(root, root_inv, p), 1);
            assert_eq!(mul_reference(n as u64, op.size_inv, p), 1);
            assert_eq!(op.size_inv_shoup, shoup_reference(op.size_inv, p));

            // powers[j] = root^j, and inverse_powers[j] = root_inv^(j + 1).
            let mut powers = vec![1; 2 * n + 1];
            let mut inverse_powers = vec![root_inv; n];
            for j in 1..=2 * n {
                powers[j] = mul_reference(powers[j - 1], root, p);
            }
            for j in 1..n {
                inverse_powers[j] = mul_reference(inverse_powers[j - 1], root_inv, p);
            }
            // Exact order 2n: no smaller positive power is one.
            assert_eq!(powers[2 * n], 1);
            assert!(powers[1..2 * n].iter().all(|&power| power != 1));

            let bits = n.trailing_zeros();
            for i in 0..n {
                let reversed = (0..bits)
                    .filter(|bit| (i >> bit) & 1 == 1)
                    .map(|bit| 1 << (bits - 1 - bit))
                    .sum::<usize>();
                assert_eq!(op.omegas[i], powers[reversed], "p={p}, n={n}, i={i}");
                assert_eq!(op.zetas_inv[i], inverse_powers[reversed]);
                assert_eq!(op.omegas_shoup[i], shoup_reference(op.omegas[i], p));
                assert_eq!(op.zetas_inv_shoup[i], shoup_reference(op.zetas_inv[i], p));
            }
        }
    }

    #[test]
    fn constructor_succeeds_exactly_for_supported_parameters() {
        for p in 2..=4000 {
            let modulus = Modulus::new(p).unwrap();
            for n in [0, 1, 4, 7, 8, 9, 16, 24, 32, 64, 128, 256, 512] {
                assert_eq!(
                    NttOperator::new(&modulus, n).is_some(),
                    supports_ntt(p, n),
                    "p={p}, n={n}"
                );
            }
        }
    }
}
