//! Contracts for the existing external CPU, RNG, and primality dependencies.
//!
//! The public Modulus alias fixes the verified implementation to Runtime. The
//! implementation of this trait is the explicit trusted boundary; arithmetic,
//! field initialization, cloning, and public wrapper bodies are checked.

#![cfg_attr(verus_keep_ghost, allow(unused_imports))]

#[cfg(verus_keep_ghost)]
use super::verified::{MODULUS_LIMIT, ResidueSampler, inversion::prime};
use rand::{CryptoRng, Rng as RngCore, RngExt, distr::Uniform};
use vstd::prelude::*;

verus! {

#[cfg(verus_keep_ghost)]
#[verifier::external_trait_specification]
pub trait ErrorSpec: core::fmt::Debug + core::fmt::Display {
    type ExternalTraitSpecificationFor: core::error::Error;
}

#[cfg(verus_keep_ghost)]
#[verifier::external_trait_specification]
pub trait TryRngSpec {
    type ExternalTraitSpecificationFor: rand::TryRng;
    type Error: core::error::Error;
}
#[cfg(verus_keep_ghost)]
#[verifier::external_trait_specification]
pub trait TryCryptoRngSpec: rand::TryRng {
    type ExternalTraitSpecificationFor: rand::TryCryptoRng;
}
#[cfg(verus_keep_ghost)]
#[verifier::external_trait_specification]
pub trait RngSpec: rand::TryRng<Error = core::convert::Infallible> {
    type ExternalTraitSpecificationFor: RngCore;
}
#[cfg(verus_keep_ghost)]
#[verifier::external_trait_specification]
pub trait CryptoRngSpec: RngCore + rand::TryCryptoRng<Error = core::convert::Infallible> {
    type ExternalTraitSpecificationFor: CryptoRng;
}

/// The contracts required of the external runtime services.
pub trait Backend: Sized {
    /// The upper bound of this backend's uniform distribution.
    spec fn modulus(&self) -> u64;

    /// Construct a distribution and a CPU dispatcher.
    fn new(p: u64) -> (r: Self)
        requires 2 <= p < MODULUS_LIMIT,
        ensures r.modulus() == p;

    /// Preserve the backend's logical state when cloning a modulus.
    fn duplicate(&self) -> (r: Self)
        ensures r == *self;

    /// Dispatch a callback exactly according to its functional contract.
    fn dispatch<T, F: FnOnce() -> T>(&self, f: F) -> (r: T)
        requires f.requires(()),
        ensures f.ensures((), r);

    /// Dispatch an in-place operation, retaining the caller's mutable borrow.
    fn dispatch_mut<'a, F: FnOnce(&'a mut [u64]) -> &'a mut [u64]>(
        &self, a: &'a mut [u64], f: F,
    ) -> (r: &'a mut [u64])
        requires f.requires((old(a),)),
        ensures f.ensures((old(a),), r), final(a)@ == final(r)@;

    /// Draw a residue from the distribution constructed for this modulus.
    fn sample<R: RngCore + CryptoRng>(&self, rng: &mut R) -> (r: u64)
        ensures r < self.modulus();

    /// Connect the external primality detector to mathematical primality.
    fn is_prime(p: u64) -> (r: bool)
        ensures r == prime(p as int);
}

}

/// The production backend; its operations retain the existing library calls.
#[derive(Debug, Clone)]
pub struct Runtime {
    distribution: Uniform<u64>,
    arch: pulp::Arch,
}

verus! {

#[cfg(verus_keep_ghost)]
#[verifier::external_type_specification]
#[verifier::external_body]
pub struct RuntimeSpec(Runtime);

#[verifier::external]
impl Backend for Runtime {
    // Required by Cargo's ghost-enabled Rust check. Verus ignores this external
    // implementation and reasons about the Backend trait contracts instead.
    spec fn modulus(&self) -> u64 { arbitrary() }

    fn new(p: u64) -> Self {
        Self { distribution: Uniform::new(0, p).unwrap(), arch: pulp::Arch::new() }
    }
    fn duplicate(&self) -> Self { self.clone() }

    #[inline(always)]
    fn dispatch<T, F: FnOnce() -> T>(&self, f: F) -> T { self.arch.dispatch(f) }

    #[inline(always)]
    fn dispatch_mut<'a, F: FnOnce(&'a mut [u64]) -> &'a mut [u64]>(
        &self, a: &'a mut [u64], f: F,
    ) -> &'a mut [u64] { self.arch.dispatch(|| f(a)) }

    #[inline]
    fn sample<R: RngCore + CryptoRng>(&self, rng: &mut R) -> u64 {
        rng.sample(self.distribution)
    }
    fn is_prime(p: u64) -> bool { fhe_util::is_prime(p) }
}

}

verus! {

/// A checked bridge between one backend draw and the verified sampling loops.
pub(super) struct Sampler<'a, B: Backend, R> {
    backend: &'a B,
    rng: &'a mut R,
    #[allow(dead_code)] // Read by specifications; erased in ordinary builds.
    history: Ghost<Seq<u64>>,
}

impl<'a, B: Backend, R: RngCore + CryptoRng> Sampler<'a, B, R> {
    pub(super) fn new(backend: &'a B, rng: &'a mut R) -> (r: Self)
        ensures r.modulus() == backend.modulus(), r.history() == Seq::<u64>::empty(),
    { Self { backend, rng, history: Ghost(Seq::empty()) } }
}

impl<B: Backend, R: RngCore + CryptoRng> super::verified::ResidueSampler for Sampler<'_, B, R> {
    closed spec fn modulus(&self) -> u64 { self.backend.modulus() }
    closed spec fn history(&self) -> Seq<u64> { self.history@ }

    fn sample(&mut self) -> (r: u64) {
        let r = self.backend.sample(self.rng);
        proof { self.history = Ghost(self.history@.push(r)); }
        r
    }
}

}
