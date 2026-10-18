//! Fill vectors with canonical samples from an abstract residue sampler.
//!
//! The production adapter to `rand::Uniform` lives outside the verified module.
//! Its sampling contract is a trusted boundary; this proves the consuming
//! loops, not the distribution's uniformity or the RNG's cryptographic
//! properties.

use vstd::prelude::*;

verus! {

/// Contract for one draw from a fixed modulus' distribution.
///
/// `history` records returned samples, not raw RNG words or rejection attempts.
/// Both specifications erase in ordinary Rust builds. Implementations outside
/// Verus must uphold this contract for the vector proofs to apply.
pub(crate) trait ResidueSampler {
    spec fn modulus(&self) -> u64;
    spec fn history(&self) -> Seq<u64>;

    fn sample(&mut self) -> (r: u64)
        ensures
            final(self).modulus() == old(self).modulus(),
            r < old(self).modulus(),
            final(self).history() == old(self).history().push(r);
}

/// Allocate exactly `size` samples, in draw order.
#[inline]
pub(crate) fn random_vec<S: ResidueSampler>(size: usize, sampler: &mut S) -> (r: Vec<u64>)
    requires
        2 <= old(sampler).modulus() < 0x4000_0000_0000_0000u64,
    ensures
        r@.len() == size,
        final(sampler).modulus() == old(sampler).modulus(),
        forall|j: int| 0 <= j < r@.len() ==> #[trigger] r@[j] < final(sampler).modulus(),
        final(sampler).history() == old(sampler).history() + r@,
{
    let ghost initial_history = sampler.history();
    let ghost p = sampler.modulus();
    let mut r: Vec<u64> = Vec::with_capacity(size);
    let mut i = 0;
    while i < size
        invariant
            0 <= i <= size,
            r@.len() == i,
            sampler.modulus() == p,
            forall|j: int| 0 <= j < i ==> #[trigger] r@[j] < p,
            sampler.history() == initial_history + r@,
        decreases size - i,
    {
        let value = sampler.sample();
        r.push(value);
        proof { assert(sampler.history() =~= initial_history + r@); }
        i += 1;
    }
    r
}

/// Replace every entry with one sample, in draw order, preserving the length.
#[inline]
pub(crate) fn random_into<S: ResidueSampler>(out: &mut [u64], sampler: &mut S)
    requires
        2 <= old(sampler).modulus() < 0x4000_0000_0000_0000u64,
    ensures
        final(out)@.len() == old(out)@.len(),
        final(sampler).modulus() == old(sampler).modulus(),
        forall|j: int| 0 <= j < final(out)@.len() ==>
            #[trigger] final(out)@[j] < final(sampler).modulus(),
        final(sampler).history() == old(sampler).history() + final(out)@,
{
    let ghost initial_history = sampler.history();
    let ghost p = sampler.modulus();
    let mut i = 0;
    while i < out.len()
        invariant
            0 <= i <= out.len(),
            out@.len() == old(out)@.len(),
            sampler.modulus() == p,
            forall|j: int| 0 <= j < i ==> #[trigger] out@[j] < p,
            forall|j: int| i <= j < out@.len() ==> #[trigger] out@[j] == old(out)@[j],
            sampler.history() == initial_history + out@.subrange(0, i as int),
        decreases out.len() - i,
    {
        let value = sampler.sample();
        out[i] = value;
        proof {
            assert(sampler.history() =~= initial_history + out@.subrange(0, i as int + 1));
        }
        i += 1;
    }
    proof { assert(out@.subrange(0, i as int) =~= out@); }
}

}
