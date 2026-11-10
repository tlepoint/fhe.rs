//! Verification regressions using only the concrete public Modulus API.
//!
//! These functions are proof clients, not alternate implementations. They do
//! not inspect fields or assume arithmetic metadata predicates.
#![allow(dead_code)]
#![allow(unused_imports)] // Cargo Verus also checks an erased Rust view.

use super::primes::{generate_prime, prime_search_result};
use super::{Modulus, verified::MODULUS_LIMIT};
use vstd::arithmetic::power2::{lemma_pow2_strictly_increases, lemma2_to64, lemma2_to64_rest};
use vstd::prelude::*;

verus! {

fn scalar_client(p: u64, a: u64, b: u64) -> (r: u64)
    requires 2 <= p < MODULUS_LIMIT, a < p, b < p,
    ensures r < p, r as int == (a + b) % (p as int),
{
    let modulus = Modulus::new(p).unwrap();
    let copy = modulus.clone();
    let value = *copy;
    assert(value == p);
    let equal = copy == modulus;
    assert(equal);
    copy.add(a, b)
}

fn vector_client(p: u64, values: &mut [u64])
    requires 2 <= p < MODULUS_LIMIT,
    ensures final(values)@.len() == old(values)@.len(),
        forall|i: int| 0 <= i < final(values)@.len() ==>
            #[trigger] final(values)@[i] < p
            && final(values)@[i] as int == old(values)@[i] as int % (p as int),
{
    let modulus = Modulus::new(p).unwrap();
    modulus.reduce_vec(values);
}

fn optimized_client(p: u64, a: u64, b: u64) -> (r: u64)
    requires 2 <= p < MODULUS_LIMIT, a < p, b < p,
    ensures r < p, r as int == (a as int * b as int) % (p as int),
{
    let modulus = Modulus::new(p).unwrap();
    if super::primes::supports_opt(p) {
        modulus.mul_opt(a, b)
    } else {
        modulus.mul(a, b)
    }
}

fn prime_search_client(num_bits: usize, modulo: u64, upper_bound: u64) -> (r: Option<Modulus>)
    ensures match r {
        Some(m) => prime_search_result(Some(m.value()), num_bits, modulo, upper_bound),
        None => prime_search_result(None, num_bits, modulo, upper_bound),
    },
{
    match generate_prime(num_bits, modulo, upper_bound) {
        None => None,
        Some(p) => {
            proof {
                lemma2_to64();
                lemma2_to64_rest();
                lemma_pow2_strictly_increases(0, (num_bits - 1) as nat);
                if num_bits < 62 {
                    lemma_pow2_strictly_increases(num_bits as nat, 62);
                }
            }
            Some(Modulus::new(p).unwrap())
        },
    }
}

}
