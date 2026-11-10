//! The public modulus implementation, verified against explicit runtime
//! contracts.

#![cfg_attr(verus_keep_ghost, allow(unused_imports))]

#[cfg(verus_keep_ghost)]
use super::verified::{
    inversion::prime,
    lazy_reduce::{base, valid_barrett, valid_opt},
    multiplication::{shoup_quotient, valid_shoup},
    serialization,
};
use super::{backend::Backend, verified};
use crate::error::{Error, Result};
#[cfg(verus_keep_ghost)]
use fhe_util::packing;
use rand::{CryptoRng, Rng as RngCore};
#[cfg(verus_keep_ghost)]
use vstd::arithmetic::power::pow as int_pow;
use vstd::prelude::*;
#[cfg(verus_keep_ghost)]
use vstd::std_specs::bits::u64_leading_zeros;

verus! {

#[cfg(verus_keep_ghost)]
#[verifier::external_type_specification]
pub struct ErrorSpec(Error);
#[cfg(verus_keep_ghost)]
#[verifier::external_type_specification]
#[verifier::external_body]
pub struct DecodeLimitErrorSpec(fhe_util::DecodeLimitError);
#[cfg(verus_keep_ghost)]
#[verifier::external_type_specification]
#[verifier::external_body]
pub struct PolynomialSerializationErrorSpec(crate::error::PolynomialSerializationError);
#[cfg(verus_keep_ghost)]
#[verifier::external_type_specification]
#[verifier::external_body]
pub struct RepresentationSpec(crate::rq::Representation);

/// Arithmetic state and the external runtime services it depends on.
#[derive(Debug)]
pub struct Modulus<B: Backend> {
    p: u64,
    barrett_hi: u64,
    barrett_lo: u64,
    leading_zeros: u32,
    supports_opt: bool,
    backend: B,
}

impl<B: Backend> Modulus<B> {
    /// The represented modulus, available to downstream specifications.
    pub closed spec fn value(&self) -> u64 { self.p }

    /// Whether the modulus admits the optimized reduction kernels.
    pub closed spec fn optimized(&self) -> bool { self.supports_opt }

    /// Initialized arithmetic metadata and the matching sampling distribution.
    #[verifier::type_invariant]
    pub closed spec fn invariant(&self) -> bool {
        &&& valid_barrett(self.p, self.barrett_hi, self.barrett_lo)
        &&& self.leading_zeros == u64_leading_zeros(self.p)
        &&& 2 <= self.leading_zeros <= 62
        &&& self.supports_opt == valid_opt(self.p, self.leading_zeros)
        &&& self.backend.modulus() == self.p
    }
}

impl<B: Backend> Clone for Modulus<B> {
    fn clone(&self) -> (r: Self)
        ensures r == *self,
    {
        proof {
            use_type_invariant(self);
        }
        Self {
            p: self.p,
            barrett_hi: self.barrett_hi,
            barrett_lo: self.barrett_lo,
            leading_zeros: self.leading_zeros,
            supports_opt: self.supports_opt,
            backend: self.backend.duplicate(),
        }
    }
}

impl<B: Backend> Eq for Modulus<B> {}
#[cfg(verus_keep_ghost)]
impl<B: Backend> vstd::std_specs::cmp::PartialEqSpecImpl for Modulus<B> {
    open spec fn obeys_eq_spec() -> bool { true }
    open spec fn eq_spec(&self, other: &Self) -> bool { self.value() == other.value() }
}
impl<B: Backend> PartialEq for Modulus<B> {
    fn eq(&self, other: &Self) -> (r: bool)
        ensures r == (self.value() == other.value()),
    { self.p == other.p }
}
impl<B: Backend> std::ops::Deref for Modulus<B> {
    type Target = u64;
    fn deref(&self) -> (r: &u64)
        ensures *r == self.value(), 2 <= *r < verified::MODULUS_LIMIT,
    {
        proof {
            use_type_invariant(self);
        }
        &self.p
    }
}

impl<B: Backend> Modulus<B> {
    /// Read the modulus in const contexts without exposing mutable metadata.
    pub(crate) const fn as_u64(&self) -> (r: u64)
        ensures r == self.value(), 2 <= r < verified::MODULUS_LIMIT,
    {
        proof {
            use_type_invariant(self);
        }
        self.p
    }

    fn canonical(&self, a: &[u64]) -> (r: bool)
        ensures r == (forall|i: int| 0 <= i < a@.len() ==> #[trigger] a@[i] < self.value()),
    {
        let mut i = 0;
        while i < a.len()
            invariant 0 <= i <= a.len(),
                forall|j: int| 0 <= j < i ==> #[trigger] a@[j] < self.value(),
            decreases a.len() - i,
        {
            if a[i] >= self.p { return false; }
            i += 1;
        }
        true
    }

    fn shoup_words_match(&self, b: &[u64], words: &[u64]) -> (r: bool)
        requires forall|i: int| 0 <= i < b@.len() ==> #[trigger] b@[i] < self.value(),
        ensures r == (b@.len() == words@.len() && forall|i: int| 0 <= i < b@.len() ==>
            valid_shoup(b@[i], self.value(), #[trigger] words@[i])),
    {
        proof {
            use_type_invariant(self);
        }
        if b.len() != words.len() { return false; }
        let mut i = 0;
        while i < b.len()
            invariant 0 <= i <= b.len(), b.len() == words.len(), self.invariant(),
                forall|j: int| 0 <= j < b@.len() ==> #[trigger] b@[j] < self.value(),
                forall|j: int| 0 <= j < i ==> valid_shoup(b@[j], self.value(), #[trigger] words@[j]),
            decreases b.len() - i,
        {
            let expected = self.shoup(b[i]);
            if expected != words[i] {
                proof {
                    if valid_shoup(b@[i as int], self.value(), words@[i as int]) {
                        shoup_quotient(b@[i as int], self.p, words@[i as int]);
                    }
                }
                return false;
            }
            i += 1;
        }
        true
    }

    /// Create a modulus from an integer of at most 62 bits.
    pub fn new(p: u64) -> (r: Result<Self>)
        ensures
            r.is_ok() <==> 2 <= p < verified::MODULUS_LIMIT,
            r matches Ok(m) ==> m.value() == p
                && m.optimized() == valid_opt(p, u64_leading_zeros(p) as u32),
            r matches Err(e) ==> e == Error::InvalidModulus(p),
    {
        match verified::new_modulus(p) {
            None => Err(Error::InvalidModulus(p)),
            Some(parameters) => Ok(Self {
                p: parameters.p,
                barrett_hi: parameters.barrett_hi,
                barrett_lo: parameters.barrett_lo,
                leading_zeros: parameters.leading_zeros,
                supports_opt: parameters.supports_opt,
                backend: B::new(p),
            }),
        }
    }

    /// Performs the modular addition of a and b in constant time.
    /// Aborts if a >= p or b >= p in debug mode.
    #[must_use]
    pub const fn add(&self, a: u64, b: u64) -> (r: u64)
        requires a < self.value(), b < self.value(),
        ensures r < self.value(), r as int == (a + b) % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(a < self.p && b < self.p);
        verified::add_mod(a, b, self.p)
    }

    /// Performs the modular addition of a and b in variable time.
    /// Aborts if a >= p or b >= p in debug mode.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being added.
    #[must_use]
    pub const unsafe fn add_vt(&self, a: u64, b: u64) -> (r: u64)
        requires a < self.value(), b < self.value(),
        ensures r < self.value(), r as int == (a + b) % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(a < self.p && b < self.p);
        unsafe { verified::add_mod_vt(a, b, self.p) }
    }

    /// Performs the modular subtraction of a and b in constant time.
    /// Aborts if a >= p or b >= p in debug mode.
    #[must_use]
    pub const fn sub(&self, a: u64, b: u64) -> (r: u64)
        requires a < self.value(), b < self.value(),
        ensures r < self.value(), r as int == (a as int - b as int) % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(a < self.p && b < self.p);
        verified::sub_mod(a, b, self.p)
    }

    /// Performs the modular multiplication of a and b in constant time.
    /// Aborts if a >= p or b >= p in debug mode.
    #[must_use]
    pub const fn mul(&self, a: u64, b: u64) -> (r: u64)
        requires a < self.value(), b < self.value(),
        ensures r < self.value(), r as int == (a as int * b as int) % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(a < self.p && b < self.p);
        verified::mul_mod(a, b, self.p, self.barrett_hi, self.barrett_lo)
    }

    /// Performs the modular multiplication of a and b in variable time.
    /// Aborts if a >= p or b >= p in debug mode.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being multiplied.
    // Production vector loops call the verified scalar kernel directly.
    #[cfg(test)]
    pub(super) const unsafe fn mul_vt(&self, a: u64, b: u64) -> (r: u64)
        requires a < self.value(), b < self.value(),
        ensures r < self.value(), r as int == (a as int * b as int) % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(a < self.p && b < self.p);
        unsafe { verified::mul_mod_vt(a, b, self.p, self.barrett_hi, self.barrett_lo) }
    }

    /// Optimized modular multiplication of a and b in constant time.
    ///
    /// Aborts if a >= p or b >= p in debug mode.
    #[must_use]
    pub const fn mul_opt(&self, a: u64, b: u64) -> (r: u64)
        requires a < self.value(), b < self.value(), self.optimized(),
        ensures r < self.value(), r as int == (a as int * b as int) % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(self.supports_opt);
        debug_assert!(a < self.p && b < self.p);

        verified::mul_opt(
            a,
            b,
            self.p,
            self.barrett_hi,
            self.barrett_lo,
            self.leading_zeros,
        )
    }

    /// Optimized modular multiplication of a and b in variable time.
    /// Aborts if a >= p or b >= p in debug mode.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being multiplied.
    #[cfg(test)]
    pub(super) const unsafe fn mul_opt_vt(&self, a: u64, b: u64) -> (r: u64)
        requires a < self.value(), b < self.value(), self.optimized(),
        ensures r < self.value(), r as int == (a as int * b as int) % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(self.supports_opt);
        debug_assert!(a < self.p && b < self.p);

        unsafe {
            verified::mul_opt_vt(
                a,
                b,
                self.p,
                self.barrett_hi,
                self.barrett_lo,
                self.leading_zeros,
            )
        }
    }

    /// Modular negation in constant time.
    ///
    /// Aborts if a >= p in debug mode.
    #[must_use]
    pub const fn neg(&self, a: u64) -> (r: u64)
        requires a < self.value(),
        ensures r < self.value(), r as int == (-(a as int)) % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(a < self.p);
        verified::neg_mod(a, self.p)
    }

    /// Modular negation in variable time.
    /// Aborts if a >= p in debug mode.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the value being negated.
    // Production vector loops call the verified scalar kernel directly.
    #[cfg(test)]
    pub(super) const unsafe fn neg_vt(&self, a: u64) -> (r: u64)
        requires a < self.value(),
        ensures r < self.value(), r as int == (-(a as int)) % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(a < self.p);
        unsafe { verified::neg_mod_vt(a, self.p) }
    }

    /// Compute the Shoup representation of a.
    ///
    /// Aborts if a >= p in debug mode.
    #[must_use]
    pub const fn shoup(&self, a: u64) -> (r: u64)
        requires a < self.value(),
        ensures valid_shoup(a, self.value(), r), r as int == a * base() / (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(a < self.p);
        verified::shoup(a, self.p)
    }

    /// Shoup multiplication of a and b in constant time.
    /// The first operand a may be any u64; b must be canonical.
    ///
    /// Aborts if b >= p or b_shoup != shoup(b) in debug mode.
    #[must_use]
    pub const fn mul_shoup(&self, a: u64, b: u64, b_shoup: u64) -> (r: u64)
        requires valid_shoup(b, self.value(), b_shoup),
        ensures r < self.value(), r as int == (a as int * b as int) % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(b < self.p);
        proof {
            shoup_quotient(b, self.p, b_shoup);
        }
        debug_assert!(b_shoup == self.shoup(b));
        verified::mul_shoup(a, b, b_shoup, self.p)
    }

    /// Shoup multiplication of a and b in variable time.
    /// The first operand a may be any u64; b must be canonical.
    /// Aborts if b >= p or b_shoup != shoup(b) in debug mode.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being multiplied.
    #[cfg(test)]
    pub(super) const unsafe fn mul_shoup_vt(&self, a: u64, b: u64, b_shoup: u64) -> (r: u64)
        requires valid_shoup(b, self.value(), b_shoup),
        ensures r < self.value(), r as int == (a as int * b as int) % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(b < self.p);
        proof {
            shoup_quotient(b, self.p, b_shoup);
        }
        debug_assert!(b_shoup == self.shoup(b));
        unsafe { verified::mul_shoup_vt(a, b, b_shoup, self.p) }
    }

    /// Lazy Shoup multiplication of a and b in constant time.
    /// The output is in the interval [0, 2 * p).
    /// The first operand a may be any u64; b must be canonical.
    ///
    /// Aborts if b >= p or b_shoup != shoup(b) in debug mode.
    #[must_use]
    pub const fn lazy_mul_shoup(&self, a: u64, b: u64, b_shoup: u64) -> (r: u64)
        requires valid_shoup(b, self.value(), b_shoup),
        ensures
            r < 2 * self.value(),
            r as int % (self.value() as int) == (a as int * b as int) % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(b < self.p);
        proof {
            shoup_quotient(b, self.p, b_shoup);
        }
        debug_assert!(b_shoup == self.shoup(b));

        let r = verified::lazy_mul_shoup(a, b, b_shoup, self.p);
        debug_assert!(r < 2 * self.p);

        r
    }

    /// Modular addition of vectors in place in constant time.
    ///
    /// Aborts if a and b differ in size, and if any of their values is >= p in
    /// debug mode.
    pub fn add_vec(&self, a: &mut [u64], b: &[u64])
        requires
            old(a)@.len() == b@.len(),
            forall|i: int| 0 <= i < old(a)@.len() ==> #[trigger] old(a)@[i] < self.value(),
            forall|i: int| 0 <= i < b@.len() ==> #[trigger] b@[i] < self.value(),
        ensures
            final(a)@.len() == old(a)@.len(),
            forall|i: int| 0 <= i < final(a)@.len() ==>
                #[trigger] final(a)@[i] < self.value()
                && final(a)@[i] as int == (old(a)@[i] as int + b@[i] as int) % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(a.len() == b.len());
        debug_assert!(self.canonical(a));
        debug_assert!(self.canonical(b));
        let ghost before = a@;
        let _ = self.backend.dispatch_mut(a, |a: &mut [u64]| -> (result: &mut [u64])
            requires
                old(a)@.len() == b@.len(),
                forall|i: int| 0 <= i < old(a)@.len() ==> #[trigger] old(a)@[i] < self.value(),
                forall|i: int| 0 <= i < b@.len() ==> #[trigger] b@[i] < self.value(),
                old(a)@ == before,
            ensures
                result@.len() == before.len(),
                forall|i: int| 0 <= i < result@.len() ==>
                    #[trigger] result@[i] < self.value()
                    && result@[i] as int == (before[i] as int + b@[i] as int) % (self.value() as int),
                final(a)@ == final(result)@,
        {
            proof {
                use_type_invariant(self);
            }
            verified::add_vec(a, b, self.p);
            a
        });
    }

    /// Modular addition of vectors in place in variable time.
    /// Aborts if a and b differ in size, and if any of their values is >= p in
    /// debug mode.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being added.
    pub unsafe fn add_vec_vt(&self, a: &mut [u64], b: &[u64])
        requires
            old(a)@.len() == b@.len(),
            forall|i: int| 0 <= i < old(a)@.len() ==> #[trigger] old(a)@[i] < self.value(),
            forall|i: int| 0 <= i < b@.len() ==> #[trigger] b@[i] < self.value(),
        ensures
            final(a)@.len() == old(a)@.len(),
            forall|i: int| 0 <= i < final(a)@.len() ==>
                #[trigger] final(a)@[i] < self.value()
                && final(a)@[i] as int == (old(a)@[i] as int + b@[i] as int) % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(a.len() == b.len());
        debug_assert!(self.canonical(a));
        debug_assert!(self.canonical(b));
        let ghost before = a@;
        let _ = self.backend.dispatch_mut(a, |a: &mut [u64]| -> (result: &mut [u64])
            requires
                old(a)@.len() == b@.len(),
                forall|i: int| 0 <= i < old(a)@.len() ==> #[trigger] old(a)@[i] < self.value(),
                forall|i: int| 0 <= i < b@.len() ==> #[trigger] b@[i] < self.value(),
                old(a)@ == before,
            ensures
                result@.len() == before.len(),
                forall|i: int| 0 <= i < result@.len() ==>
                    #[trigger] result@[i] < self.value()
                    && result@[i] as int == (before[i] as int + b@[i] as int) % (self.value() as int),
                final(a)@ == final(result)@,
        {
            proof {
                use_type_invariant(self);
            }
            unsafe { verified::add_vec_vt(a, b, self.p) };
            a
        });
    }

    /// Modular subtraction of vectors in place in constant time.
    ///
    /// Aborts if a and b differ in size, and if any of their values is >= p in
    /// debug mode.
    pub fn sub_vec(&self, a: &mut [u64], b: &[u64])
        requires
            old(a)@.len() == b@.len(),
            forall|i: int| 0 <= i < old(a)@.len() ==> #[trigger] old(a)@[i] < self.value(),
            forall|i: int| 0 <= i < b@.len() ==> #[trigger] b@[i] < self.value(),
        ensures
            final(a)@.len() == old(a)@.len(),
            forall|i: int| 0 <= i < final(a)@.len() ==>
                #[trigger] final(a)@[i] < self.value()
                && final(a)@[i] as int == (old(a)@[i] as int - b@[i] as int) % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(a.len() == b.len());
        debug_assert!(self.canonical(a));
        debug_assert!(self.canonical(b));
        let ghost before = a@;
        let _ = self.backend.dispatch_mut(a, |a: &mut [u64]| -> (result: &mut [u64])
            requires
                old(a)@.len() == b@.len(),
                forall|i: int| 0 <= i < old(a)@.len() ==> #[trigger] old(a)@[i] < self.value(),
                forall|i: int| 0 <= i < b@.len() ==> #[trigger] b@[i] < self.value(),
                old(a)@ == before,
            ensures
                result@.len() == before.len(),
                forall|i: int| 0 <= i < result@.len() ==>
                    #[trigger] result@[i] < self.value()
                    && result@[i] as int == (before[i] as int - b@[i] as int) % (self.value() as int),
                final(a)@ == final(result)@,
        {
            proof {
                use_type_invariant(self);
            }
            verified::sub_vec(a, b, self.p);
            a
        });
    }

    /// Modular subtraction of vectors in place in variable time.
    /// Aborts if a and b differ in size, and if any of their values is >= p in
    /// debug mode.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being subtracted.
    pub unsafe fn sub_vec_vt(&self, a: &mut [u64], b: &[u64])
        requires
            old(a)@.len() == b@.len(),
            forall|i: int| 0 <= i < old(a)@.len() ==> #[trigger] old(a)@[i] < self.value(),
            forall|i: int| 0 <= i < b@.len() ==> #[trigger] b@[i] < self.value(),
        ensures
            final(a)@.len() == old(a)@.len(),
            forall|i: int| 0 <= i < final(a)@.len() ==>
                #[trigger] final(a)@[i] < self.value()
                && final(a)@[i] as int == (old(a)@[i] as int - b@[i] as int) % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(a.len() == b.len());
        debug_assert!(self.canonical(a));
        debug_assert!(self.canonical(b));
        let ghost before = a@;
        let _ = self.backend.dispatch_mut(a, |a: &mut [u64]| -> (result: &mut [u64])
            requires
                old(a)@.len() == b@.len(),
                forall|i: int| 0 <= i < old(a)@.len() ==> #[trigger] old(a)@[i] < self.value(),
                forall|i: int| 0 <= i < b@.len() ==> #[trigger] b@[i] < self.value(),
                old(a)@ == before,
            ensures
                result@.len() == before.len(),
                forall|i: int| 0 <= i < result@.len() ==>
                    #[trigger] result@[i] < self.value()
                    && result@[i] as int == (before[i] as int - b@[i] as int) % (self.value() as int),
                final(a)@ == final(result)@,
        {
            proof {
                use_type_invariant(self);
            }
            unsafe { verified::sub_vec_vt(a, b, self.p) };
            a
        });
    }

    /// Modular multiplication of vectors in place in constant time.
    ///
    /// Aborts if a and b differ in size, and if any of their values is >= p in
    /// debug mode.
    pub fn mul_vec(&self, a: &mut [u64], b: &[u64])
        requires
            old(a)@.len() == b@.len(),
            forall|i: int| 0 <= i < old(a)@.len() ==> #[trigger] old(a)@[i] < self.value(),
            forall|i: int| 0 <= i < b@.len() ==> #[trigger] b@[i] < self.value(),
        ensures
            final(a)@.len() == old(a)@.len(),
            forall|i: int| 0 <= i < final(a)@.len() ==>
                #[trigger] final(a)@[i] < self.value()
                && final(a)@[i] as int == (old(a)@[i] as int * b@[i] as int) % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(a.len() == b.len());
        debug_assert!(self.canonical(a));
        debug_assert!(self.canonical(b));
        let ghost before = a@;
        let _ = self.backend.dispatch_mut(a, |a: &mut [u64]| -> (result: &mut [u64])
            requires
                old(a)@.len() == b@.len(),
                forall|i: int| 0 <= i < old(a)@.len() ==> #[trigger] old(a)@[i] < self.value(),
                forall|i: int| 0 <= i < b@.len() ==> #[trigger] b@[i] < self.value(),
                old(a)@ == before,
            ensures
                result@.len() == before.len(),
                forall|i: int| 0 <= i < result@.len() ==>
                    #[trigger] result@[i] < self.value()
                    && result@[i] as int == (before[i] as int * b@[i] as int) % (self.value() as int),
                final(a)@ == final(result)@,
        {
            proof {
                use_type_invariant(self);
            }
            if self.supports_opt {
                verified::mul_vec::<true>(a, b, self.p, self.barrett_hi, self.barrett_lo, self.leading_zeros);
            } else {
                verified::mul_vec::<false>(a, b, self.p, self.barrett_hi, self.barrett_lo, self.leading_zeros);
            }
            a
        });
    }

    /// Modular scalar multiplication of vectors in place in constant time.
    /// Entries of a may be any u64; b must be canonical.
    ///
    /// Aborts if b >= p in debug mode.
    pub fn scalar_mul_vec(&self, a: &mut [u64], b: u64)
        requires b < self.value(),
        ensures
            final(a)@.len() == old(a)@.len(),
            forall|i: int| 0 <= i < final(a)@.len() ==>
                #[trigger] final(a)@[i] < self.value()
                && final(a)@[i] as int == (old(a)@[i] as int * b as int) % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(b < self.p);
        let ghost before = a@;
        let _ = self.backend.dispatch_mut(a, |a: &mut [u64]| -> (result: &mut [u64])
            requires b < self.value(), old(a)@ == before,
            ensures
                result@.len() == before.len(),
                forall|i: int| 0 <= i < result@.len() ==>
                    #[trigger] result@[i] < self.value()
                    && result@[i] as int == (before[i] as int * b as int) % (self.value() as int),
                final(a)@ == final(result)@,
        {
            proof {
                use_type_invariant(self);
            }
            verified::scalar_mul_vec(a, b, self.p);
            a
        });
    }

    /// Modular scalar multiplication of vectors in place in variable time.
    /// Entries of a may be any u64; b must be canonical.
    /// Aborts if b >= p in debug mode.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being multiplied.
    pub unsafe fn scalar_mul_vec_vt(&self, a: &mut [u64], b: u64)
        requires b < self.value(),
        ensures
            final(a)@.len() == old(a)@.len(),
            forall|i: int| 0 <= i < final(a)@.len() ==>
                #[trigger] final(a)@[i] < self.value()
                && final(a)@[i] as int == (old(a)@[i] as int * b as int) % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(b < self.p);
        let ghost before = a@;
        let _ = self.backend.dispatch_mut(a, |a: &mut [u64]| -> (result: &mut [u64])
            requires b < self.value(), old(a)@ == before,
            ensures
                result@.len() == before.len(),
                forall|i: int| 0 <= i < result@.len() ==>
                    #[trigger] result@[i] < self.value()
                    && result@[i] as int == (before[i] as int * b as int) % (self.value() as int),
                final(a)@ == final(result)@,
        {
            proof {
                use_type_invariant(self);
            }
            unsafe { verified::scalar_mul_vec_vt(a, b, self.p) };
            a
        });
    }

    /// Modular multiplication of vectors in place in variable time.
    /// Aborts if a and b differ in size, and if any of their values is >= p in
    /// debug mode.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being multiplied.
    pub unsafe fn mul_vec_vt(&self, a: &mut [u64], b: &[u64])
        requires
            old(a)@.len() == b@.len(),
            forall|i: int| 0 <= i < old(a)@.len() ==> #[trigger] old(a)@[i] < self.value(),
            forall|i: int| 0 <= i < b@.len() ==> #[trigger] b@[i] < self.value(),
        ensures
            final(a)@.len() == old(a)@.len(),
            forall|i: int| 0 <= i < final(a)@.len() ==>
                #[trigger] final(a)@[i] < self.value()
                && final(a)@[i] as int == (old(a)@[i] as int * b@[i] as int) % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(a.len() == b.len());
        debug_assert!(self.canonical(a));
        debug_assert!(self.canonical(b));
        let ghost before = a@;
        let _ = self.backend.dispatch_mut(a, |a: &mut [u64]| -> (result: &mut [u64])
            requires
                old(a)@.len() == b@.len(),
                forall|i: int| 0 <= i < old(a)@.len() ==> #[trigger] old(a)@[i] < self.value(),
                forall|i: int| 0 <= i < b@.len() ==> #[trigger] b@[i] < self.value(),
                old(a)@ == before,
            ensures
                result@.len() == before.len(),
                forall|i: int| 0 <= i < result@.len() ==>
                    #[trigger] result@[i] < self.value()
                    && result@[i] as int == (before[i] as int * b@[i] as int) % (self.value() as int),
                final(a)@ == final(result)@,
        {
            proof {
                use_type_invariant(self);
            }
            if self.supports_opt {
                unsafe { verified::mul_vec_vt::<true>(a, b, self.p, self.barrett_hi, self.barrett_lo, self.leading_zeros) };
            } else {
                unsafe { verified::mul_vec_vt::<false>(a, b, self.p, self.barrett_hi, self.barrett_lo, self.leading_zeros) };
            }
            a
        });
    }

    /// Compute the Shoup representation of a vector.
    ///
    /// Aborts if any of the values of the vector is >= p in debug mode.
    #[must_use]
    pub fn shoup_vec(&self, a: &[u64]) -> (r: Vec<u64>)
        requires forall|i: int| 0 <= i < a@.len() ==> #[trigger] a@[i] < self.value(),
        ensures
            r@.len() == a@.len(),
            forall|i: int| 0 <= i < r@.len() ==>
                valid_shoup(a@[i], self.value(), #[trigger] r@[i])
                && r@[i] as int == a@[i] * base() / (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(self.canonical(a));
        self.backend.dispatch(|| -> (r: Vec<u64>)
            ensures
                r@.len() == a@.len(),
                forall|i: int| 0 <= i < r@.len() ==>
                    valid_shoup(a@[i], self.value(), #[trigger] r@[i])
                    && r@[i] as int == a@[i] * base() / (self.value() as int),
        {
            proof {
                use_type_invariant(self);
            }
            verified::shoup_vec(a, self.p)
        })
    }

    pub(crate) fn shoup_into(&self, a: &[u64], out: &mut [u64])
        requires
            a@.len() == old(out)@.len(),
            forall|i: int| 0 <= i < a@.len() ==> #[trigger] a@[i] < self.value(),
        ensures
            final(out)@.len() == old(out)@.len(),
            forall|i: int| 0 <= i < final(out)@.len() ==>
                valid_shoup(a@[i], self.value(), #[trigger] final(out)@[i])
                && final(out)@[i] as int == a@[i] * base() / (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        assert!(a.len() == out.len());
        debug_assert!(self.canonical(a));
        let ghost before = out@;
        let _ = self.backend.dispatch_mut(out, |out: &mut [u64]| -> (result: &mut [u64])
            requires
                a@.len() == old(out)@.len(),
                forall|i: int| 0 <= i < a@.len() ==> #[trigger] a@[i] < self.value(),
                old(out)@ == before,
            ensures
                result@.len() == before.len(),
                forall|i: int| 0 <= i < result@.len() ==>
                    valid_shoup(a@[i], self.value(), #[trigger] result@[i])
                    && result@[i] as int == a@[i] * base() / (self.value() as int),
                final(out)@ == final(result)@,
        {
            proof {
                use_type_invariant(self);
            }
            verified::shoup_into(a, out, self.p);
            out
        });
    }

    /// Shoup modular multiplication of vectors in place in constant time.
    /// Entries of a may be any u64; b must contain canonical residues.
    ///
    /// Aborts in debug mode if lengths differ, an entry of b is >= p, or a
    /// precomputed word differs from shoup(b).
    pub fn mul_shoup_vec(&self, a: &mut [u64], b: &[u64], b_shoup: &[u64])
        requires
            old(a)@.len() == b@.len(),
            b@.len() == b_shoup@.len(),
            forall|i: int| 0 <= i < b@.len() ==> valid_shoup(b@[i], self.value(), #[trigger] b_shoup@[i]),
        ensures
            final(a)@.len() == old(a)@.len(),
            forall|i: int| 0 <= i < final(a)@.len() ==>
                #[trigger] final(a)@[i] < self.value()
                && final(a)@[i] as int == (old(a)@[i] as int * b@[i] as int) % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(a.len() == b.len());
        debug_assert!(a.len() == b_shoup.len());
        proof {
            assert forall|i: int| 0 <= i < b@.len() implies #[trigger] b@[i] < self.value() by {
                assert(valid_shoup(b@[i], self.value(), b_shoup@[i]));
            }
        }
        debug_assert!(self.canonical(b));
        debug_assert!(self.shoup_words_match(b, b_shoup));
        let ghost before = a@;
        let _ = self.backend.dispatch_mut(a, |a: &mut [u64]| -> (result: &mut [u64])
            requires
                old(a)@.len() == b@.len(),
                b@.len() == b_shoup@.len(),
                forall|i: int| 0 <= i < b@.len() ==> valid_shoup(b@[i], self.value(), #[trigger] b_shoup@[i]),
                old(a)@ == before,
            ensures
                result@.len() == before.len(),
                forall|i: int| 0 <= i < result@.len() ==>
                    #[trigger] result@[i] < self.value()
                    && result@[i] as int == (before[i] as int * b@[i] as int) % (self.value() as int),
                final(a)@ == final(result)@,
        {
            proof {
                use_type_invariant(self);
            }
            verified::mul_shoup_vec(a, b, b_shoup, self.p);
            a
        });
    }

    /// Shoup modular multiplication of vectors in place in variable time.
    /// Entries of a may be any u64; b must contain canonical residues.
    /// Aborts in debug mode if lengths differ, an entry of b is >= p, or a
    /// precomputed word differs from shoup(b).
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being multiplied.
    pub unsafe fn mul_shoup_vec_vt(&self, a: &mut [u64], b: &[u64], b_shoup: &[u64])
        requires
            old(a)@.len() == b@.len(),
            b@.len() == b_shoup@.len(),
            forall|i: int| 0 <= i < b@.len() ==> valid_shoup(b@[i], self.value(), #[trigger] b_shoup@[i]),
        ensures
            final(a)@.len() == old(a)@.len(),
            forall|i: int| 0 <= i < final(a)@.len() ==>
                #[trigger] final(a)@[i] < self.value()
                && final(a)@[i] as int == (old(a)@[i] as int * b@[i] as int) % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(a.len() == b.len());
        debug_assert!(a.len() == b_shoup.len());
        proof {
            assert forall|i: int| 0 <= i < b@.len() implies #[trigger] b@[i] < self.value() by {
                assert(valid_shoup(b@[i], self.value(), b_shoup@[i]));
            }
        }
        debug_assert!(self.canonical(b));
        debug_assert!(self.shoup_words_match(b, b_shoup));
        let ghost before = a@;
        let _ = self.backend.dispatch_mut(a, |a: &mut [u64]| -> (result: &mut [u64])
            requires
                old(a)@.len() == b@.len(),
                b@.len() == b_shoup@.len(),
                forall|i: int| 0 <= i < b@.len() ==> valid_shoup(b@[i], self.value(), #[trigger] b_shoup@[i]),
                old(a)@ == before,
            ensures
                result@.len() == before.len(),
                forall|i: int| 0 <= i < result@.len() ==>
                    #[trigger] result@[i] < self.value()
                    && result@[i] as int == (before[i] as int * b@[i] as int) % (self.value() as int),
                final(a)@ == final(result)@,
        {
            proof {
                use_type_invariant(self);
            }
            unsafe { verified::mul_shoup_vec_vt(a, b, b_shoup, self.p) };
            a
        });
    }

    /// Reduce a vector in place in constant time.
    pub fn reduce_vec(&self, a: &mut [u64])
        ensures
            final(a)@.len() == old(a)@.len(),
            forall|i: int| 0 <= i < final(a)@.len() ==>
                #[trigger] final(a)@[i] < self.value()
                && final(a)@[i] as int == (old(a)@[i] as int) % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        let ghost before = a@;
        let _ = self.backend.dispatch_mut(a, |a: &mut [u64]| -> (result: &mut [u64])
            requires old(a)@ == before,
            ensures
                result@.len() == before.len(),
                forall|i: int| 0 <= i < result@.len() ==>
                    #[trigger] result@[i] < self.value()
                    && result@[i] as int == (before[i] as int) % (self.value() as int),
                final(a)@ == final(result)@,
        {
            proof {
                use_type_invariant(self);
            }
            verified::reduce_vec(a, self.p, self.barrett_hi, self.barrett_lo);
            a
        });
    }

    /// Center a value modulo p as i64 in constant time.
    ///
    /// The output is congruent to a modulo p and satisfies `-p <= 2*result <
    /// p`. For even p, the midpoint p/2 maps to -p/2.
    /// Aborts if `a >= p` in debug mode.
    #[must_use]
    pub const fn center(&self, a: u64) -> (r: i64)
        requires a < self.value(),
        ensures
            -(self.value() as int) <= 2 * (r as int) < self.value() as int,
            r as int % (self.value() as int) == a as int,
            r as int == if 2 * (a as int) < self.value() as int { a as int } else { a as int - self.value() as int },
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(a < self.p);
        verified::center(a, self.p)
    }

    /// Return centered representatives in input order, in constant time.
    /// Each output satisfies the range and congruence guarantees of
    /// [`Self::center`]. Aborts if any input is >= p in debug mode.
    #[must_use]
    pub fn center_vec(&self, a: &[u64]) -> (r: Vec<i64>)
        requires forall|i: int| 0 <= i < a@.len() ==> #[trigger] a@[i] < self.value(),
        ensures
            r@.len() == a@.len(),
            forall|i: int| 0 <= i < r@.len() ==>
                -(self.value() as int) <= 2 * (#[trigger] r@[i] as int) < self.value() as int
                && r@[i] as int % (self.value() as int) == a@[i] as int
                && r@[i] as int == if 2 * (a@[i] as int) < self.value() as int { a@[i] as int } else { a@[i] as int - self.value() as int },
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(self.canonical(a));
        self.backend.dispatch(|| -> (r: Vec<i64>)
            ensures
                r@.len() == a@.len(),
                forall|i: int| 0 <= i < r@.len() ==>
                    -(self.value() as int) <= 2 * (#[trigger] r@[i] as int) < self.value() as int
                    && r@[i] as int % (self.value() as int) == a@[i] as int
                    && r@[i] as int == if 2 * (a@[i] as int) < self.value() as int { a@[i] as int } else { a@[i] as int - self.value() as int },
        {
            proof {
                use_type_invariant(self);
            }
            verified::center_vec(a, self.p)
        })
    }

    /// Reduce a vector in place in variable time.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being reduced.
    pub unsafe fn reduce_vec_vt(&self, a: &mut [u64])
        ensures
            final(a)@.len() == old(a)@.len(),
            forall|i: int| 0 <= i < final(a)@.len() ==>
                #[trigger] final(a)@[i] < self.value()
                && final(a)@[i] as int == (old(a)@[i] as int) % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        let ghost before = a@;
        let _ = self.backend.dispatch_mut(a, |a: &mut [u64]| -> (result: &mut [u64])
            requires old(a)@ == before,
            ensures
                result@.len() == before.len(),
                forall|i: int| 0 <= i < result@.len() ==>
                    #[trigger] result@[i] < self.value()
                    && result@[i] as int == (before[i] as int) % (self.value() as int),
                final(a)@ == final(result)@,
        {
            proof {
                use_type_invariant(self);
            }
            unsafe { verified::reduce_vec_vt(a, self.p, self.barrett_hi, self.barrett_lo) };
            a
        });
    }

    /// Modular reduction of a i64 in constant time.
    // Production vector loops call the verified scalar kernel directly.
    #[cfg(test)]
    pub(super) const fn reduce_i64(&self, a: i64) -> (r: u64)
        ensures r < self.value(), r as int == a as int % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        let r = verified::reduce_i64(a, self.p, self.barrett_hi, self.barrett_lo);
        debug_assert!(r < self.p && r as i128 == (a as i128).rem_euclid(self.p as i128));
        r
    }

    /// Modular reduction of a i64 in variable time.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being reduced.
    #[cfg(test)]
    pub(super) const unsafe fn reduce_i64_vt(&self, a: i64) -> (r: u64)
        ensures r < self.value(), r as int == a as int % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        let r = unsafe { verified::reduce_i64_vt(a, self.p, self.barrett_hi, self.barrett_lo) };
        debug_assert!(r < self.p && r as i128 == (a as i128).rem_euclid(self.p as i128));
        r
    }

    /// Return a new vector of canonical residues of signed inputs in constant
    /// time.
    #[must_use]
    pub fn reduce_vec_i64(&self, a: &[i64]) -> (r: Vec<u64>)
        ensures
            r@.len() == a@.len(),
            forall|i: int| 0 <= i < r@.len() ==>
                #[trigger] r@[i] < self.value()
                && r@[i] as int == (a@[i] as int) % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        self.backend.dispatch(|| -> (r: Vec<u64>)
            ensures
                r@.len() == a@.len(),
                forall|i: int| 0 <= i < r@.len() ==>
                    #[trigger] r@[i] < self.value()
                    && r@[i] as int == (a@[i] as int) % (self.value() as int),
        {
            proof {
                use_type_invariant(self);
            }
            verified::reduce_vec_i64(a, self.p, self.barrett_hi, self.barrett_lo)
        })
    }

    /// Return a new vector of canonical residues of signed inputs in variable
    /// time.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being reduced.
    #[must_use]
    pub unsafe fn reduce_vec_i64_vt(&self, a: &[i64]) -> (r: Vec<u64>)
        ensures
            r@.len() == a@.len(),
            forall|i: int| 0 <= i < r@.len() ==>
                #[trigger] r@[i] < self.value()
                && r@[i] as int == (a@[i] as int) % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        self.backend.dispatch(|| -> (r: Vec<u64>)
            ensures
                r@.len() == a@.len(),
                forall|i: int| 0 <= i < r@.len() ==>
                    #[trigger] r@[i] < self.value()
                    && r@[i] as int == (a@[i] as int) % (self.value() as int),
        {
            proof {
                use_type_invariant(self);
            }
            unsafe { verified::reduce_vec_i64_vt(a, self.p, self.barrett_hi, self.barrett_lo) }
        })
    }

    /// Reduce a vector in constant time.
    #[must_use]
    pub fn reduce_vec_new(&self, a: &[u64]) -> (r: Vec<u64>)
        ensures
            r@.len() == a@.len(),
            forall|i: int| 0 <= i < r@.len() ==>
                #[trigger] r@[i] < self.value()
                && r@[i] as int == (a@[i] as int) % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        self.backend.dispatch(|| -> (r: Vec<u64>)
            ensures
                r@.len() == a@.len(),
                forall|i: int| 0 <= i < r@.len() ==>
                    #[trigger] r@[i] < self.value()
                    && r@[i] as int == (a@[i] as int) % (self.value() as int),
        {
            proof {
                use_type_invariant(self);
            }
            verified::reduce_vec_new(a, self.p, self.barrett_hi, self.barrett_lo)
        })
    }

    /// Reduce a vector in variable time.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being reduced.
    #[must_use]
    pub unsafe fn reduce_vec_new_vt(&self, a: &[u64]) -> (r: Vec<u64>)
        ensures
            r@.len() == a@.len(),
            forall|i: int| 0 <= i < r@.len() ==>
                #[trigger] r@[i] < self.value()
                && r@[i] as int == (a@[i] as int) % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        self.backend.dispatch(|| -> (r: Vec<u64>)
            ensures
                r@.len() == a@.len(),
                forall|i: int| 0 <= i < r@.len() ==>
                    #[trigger] r@[i] < self.value()
                    && r@[i] as int == (a@[i] as int) % (self.value() as int),
        {
            proof {
                use_type_invariant(self);
            }
            unsafe { verified::reduce_vec_new_vt(a, self.p, self.barrett_hi, self.barrett_lo) }
        })
    }

    /// Modular negation of a vector in place in constant time.
    ///
    /// Aborts if any of the values in the vector is >= p in debug mode.
    pub fn neg_vec(&self, a: &mut [u64])
        requires forall|i: int| 0 <= i < old(a)@.len() ==> #[trigger] old(a)@[i] < self.value(),
        ensures
            final(a)@.len() == old(a)@.len(),
            forall|i: int| 0 <= i < final(a)@.len() ==>
                #[trigger] final(a)@[i] < self.value()
                && final(a)@[i] as int == (-(old(a)@[i] as int)) % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(self.canonical(a));
        let ghost before = a@;
        let _ = self.backend.dispatch_mut(a, |a: &mut [u64]| -> (result: &mut [u64])
            requires
                forall|i: int| 0 <= i < old(a)@.len() ==> #[trigger] old(a)@[i] < self.value(),
                old(a)@ == before,
            ensures
                result@.len() == before.len(),
                forall|i: int| 0 <= i < result@.len() ==>
                    #[trigger] result@[i] < self.value()
                    && result@[i] as int == (-(before[i] as int)) % (self.value() as int),
                final(a)@ == final(result)@,
        {
            proof {
                use_type_invariant(self);
            }
            verified::neg_vec(a, self.p);
            a
        });
    }

    /// Modular negation of a vector in place in variable time.
    /// Aborts if any of the values in the vector is >= p in debug mode.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being negated.
    pub unsafe fn neg_vec_vt(&self, a: &mut [u64])
        requires forall|i: int| 0 <= i < old(a)@.len() ==> #[trigger] old(a)@[i] < self.value(),
        ensures
            final(a)@.len() == old(a)@.len(),
            forall|i: int| 0 <= i < final(a)@.len() ==>
                #[trigger] final(a)@[i] < self.value()
                && final(a)@[i] as int == (-(old(a)@[i] as int)) % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(self.canonical(a));
        let ghost before = a@;
        let _ = self.backend.dispatch_mut(a, |a: &mut [u64]| -> (result: &mut [u64])
            requires
                forall|i: int| 0 <= i < old(a)@.len() ==> #[trigger] old(a)@[i] < self.value(),
                old(a)@ == before,
            ensures
                result@.len() == before.len(),
                forall|i: int| 0 <= i < result@.len() ==>
                    #[trigger] result@[i] < self.value()
                    && result@[i] as int == (-(before[i] as int)) % (self.value() as int),
                final(a)@ == final(result)@,
        {
            proof {
                use_type_invariant(self);
            }
            unsafe { verified::neg_vec_vt(a, self.p) };
            a
        });
    }

    /// Modular exponentiation in variable time.
    ///
    /// Returns the canonical residue of a^n modulo p. Exponent zero returns 1,
    /// including when a is zero.
    /// Aborts if a >= p or n >= p in debug mode.
    #[must_use]
    pub fn pow(&self, a: u64, n: u64) -> (r: u64)
        requires a < self.value(), n < self.value(),
        ensures r < self.value(), r as int == int_pow(a as int, n as nat) % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(a < self.p && n < self.p);

        verified::pow_mod(a, n, self.p, self.barrett_hi, self.barrett_lo)
    }

    /// Modular inversion in variable time.
    ///
    /// Returns the canonical multiplicative inverse of a, or None if p is not
    /// prime or a = 0. Composite moduli are rejected even when a is invertible.
    /// Aborts if p is prime and a >= p in debug mode.
    #[must_use]
    pub fn inv(&self, a: u64) -> (r: std::option::Option<u64>)
        requires prime(self.value() as int) ==> a < self.value(),
        ensures
            r.is_some() <==> prime(self.value() as int) && a != 0,
            r matches Some(x) ==> 0 < x < self.value() && (a as int * x as int) % (self.value() as int) == 1,
    {
        proof {
            use_type_invariant(self);
        }
        let p_is_prime = B::is_prime(self.p);
        debug_assert!(!p_is_prime || a < self.p);
        let result = verified::inv_mod(a, self.p, self.barrett_hi, self.barrett_lo, p_is_prime);
        if let Some(r) = result {
            debug_assert!(self.mul(a, r) == 1);
        }
        result
    }

    /// Modular reduction of a u128 in constant time.
    #[must_use]
    pub const fn reduce_u128(&self, a: u128) -> (r: u64)
        ensures r < self.value(), r as int == a as int % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        let r = verified::reduce_u128(a, self.p, self.barrett_hi, self.barrett_lo);
        debug_assert!(r < self.p && r as u128 == a % self.p as u128);
        r
    }

    /// Modular reduction of a u128 in variable time.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the value being reduced.
    #[must_use]
    pub const unsafe fn reduce_u128_vt(&self, a: u128) -> (r: u64)
        ensures r < self.value(), r as int == a as int % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        let r = unsafe { verified::reduce_u128_vt(a, self.p, self.barrett_hi, self.barrett_lo) };
        debug_assert!(r < self.p && r as u128 == a % self.p as u128);
        r
    }

    /// Modular reduction of a u64 in constant time.
    #[must_use]
    pub const fn reduce(&self, a: u64) -> (r: u64)
        ensures r < self.value(), r as int == a as int % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        let r = verified::reduce(a, self.p, self.barrett_hi, self.barrett_lo);
        debug_assert!(r < self.p && r == a % self.p);
        r
    }

    /// Modular reduction of a u64 in variable time.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the value being reduced.
    #[must_use]
    pub const unsafe fn reduce_vt(&self, a: u64) -> (r: u64)
        ensures r < self.value(), r as int == a as int % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        let r = unsafe { verified::reduce_vt(a, self.p, self.barrett_hi, self.barrett_lo) };
        debug_assert!(r < self.p && r == a % self.p);
        r
    }

    /// Optimized modular reduction of a u128 in constant time.
    /// Requires an eligible modulus and `a < p^2`; checked in debug mode.
    #[must_use]
    pub const fn reduce_opt_u128(&self, a: u128) -> (r: u64)
        requires self.optimized(), a < self.value() * self.value(),
        ensures r < self.value(), r as int == a as int % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(self.supports_opt);
        proof {
            assert((self.p as int) * (self.p as int) <= u128::MAX) by (nonlinear_arith)
                requires 2 <= self.p < verified::MODULUS_LIMIT;
        }
        debug_assert!(a < (self.p as u128) * (self.p as u128));
        let r = verified::reduce_opt_u128(
            a,
            self.p,
            self.barrett_hi,
            self.barrett_lo,
            self.leading_zeros,
        );
        debug_assert!(r < self.p && r as u128 == a % self.p as u128);
        r
    }

    /// Optimized modular reduction of a u128 in variable time.
    /// Requires an eligible modulus and `a < p^2`; checked in debug mode.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the value being reduced.
    // Production multiplication calls the verified reducer directly.
    #[cfg(test)]
    pub(crate) const unsafe fn reduce_opt_u128_vt(&self, a: u128) -> (r: u64)
        requires self.optimized(), a < self.value() * self.value(),
        ensures r < self.value(), r as int == a as int % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(self.supports_opt);
        proof {
            assert((self.p as int) * (self.p as int) <= u128::MAX) by (nonlinear_arith)
                requires 2 <= self.p < verified::MODULUS_LIMIT;
        }
        debug_assert!(a < (self.p as u128) * (self.p as u128));
        let r = unsafe {
            verified::reduce_opt_u128_vt(
                a,
                self.p,
                self.barrett_hi,
                self.barrett_lo,
                self.leading_zeros,
            )
        };
        debug_assert!(r < self.p && r as u128 == a % self.p as u128);
        r
    }

    /// Optimized modular reduction of a u64 in constant time.
    /// Requires a modulus supporting optimized reduction; checked in debug
    /// mode.
    #[must_use]
    pub const fn reduce_opt(&self, a: u64) -> (r: u64)
        requires self.optimized(),
        ensures r < self.value(), r as int == a as int % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(self.supports_opt);
        let r = verified::reduce_opt(a, self.p, self.leading_zeros);
        debug_assert!(r < self.p && r == a % self.p);
        r
    }

    /// Optimized modular reduction of a u64 in variable time.
    /// Requires a modulus supporting optimized reduction; checked in debug
    /// mode.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the value being reduced.
    #[must_use]
    pub const unsafe fn reduce_opt_vt(&self, a: u64) -> (r: u64)
        requires self.optimized(),
        ensures r < self.value(), r as int == a as int % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(self.supports_opt);
        let r = unsafe { verified::reduce_opt_vt(a, self.p, self.leading_zeros) };
        debug_assert!(r < self.p && r == a % self.p);
        r
    }

    /// Return x mod p in constant time.
    /// Aborts if x >= 2 * p in debug mode.
    pub(crate) const fn reduce1(x: u64, p: u64) -> (r: u64)
        requires 0 < p < 0x8000_0000_0000_0000u64, x < 2 * p,
        ensures r < p, r as int == x as int % p as int,
    {
        proof {
            assert(p >> 63 == 0) by (bit_vector)
                requires p < 0x8000_0000_0000_0000u64;
        }
        debug_assert!(p >> 63 == 0);
        debug_assert!(x < 2 * p);

        let r = verified::reduce1(x, p);

        debug_assert!(r == x % p);

        r
    }

    /// Return x mod p in variable time.
    /// Aborts if x >= 2 * p in debug mode.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the value being reduced.
    #[inline]
    pub(crate) const unsafe fn reduce1_vt(x: u64, p: u64) -> (r: u64)
        requires 0 < p < 0x8000_0000_0000_0000u64, x < 2 * p,
        ensures r < p, r as int == x as int % p as int,
    {
        proof {
            assert(p >> 63 == 0) by (bit_vector)
                requires p < 0x8000_0000_0000_0000u64;
        }
        debug_assert!(p >> 63 == 0);
        debug_assert!(x < 2 * p);

        let r = unsafe { verified::reduce1_vt(x, p) };
        debug_assert!(r == x % p);

        r
    }

    /// Lazy modular reduction of a in constant time.
    /// The output is in the interval [0, 2 * p).
    #[must_use]
    pub const fn lazy_reduce_u128(&self, a: u128) -> (r: u64)
        ensures r < 2 * self.value(), r as int % (self.value() as int) == a as int % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        let r = verified::lazy_reduce_u128(a, self.p, self.barrett_hi, self.barrett_lo);

        debug_assert!((r as u128) < 2 * (self.p as u128));
        debug_assert!(r % self.p == (a % (self.p as u128)) as u64);

        r
    }

    /// Lazy modular reduction of a in constant time.
    /// The output is in the interval [0, 2 * p).
    #[must_use]
    pub const fn lazy_reduce(&self, a: u64) -> (r: u64)
        ensures r < 2 * self.value(), r as int % (self.value() as int) == a as int % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        let r = verified::lazy_reduce(a, self.p, self.barrett_hi, self.barrett_lo);

        debug_assert!((r as u128) < 2 * (self.p as u128));
        debug_assert!(r % self.p == a % self.p);

        r
    }

    /// Lazy optimized modular reduction of a in constant time.
    /// The output is in the interval [0, 2 * p).
    ///
    /// Requires a modulus supporting optimized reduction and an input < p ^ 2.
    /// Aborts if either condition is violated in debug mode.
    #[must_use]
    pub const fn lazy_reduce_opt_u128(&self, a: u128) -> (r: u64)
        requires self.optimized(), a < self.value() * self.value(),
        ensures r < 2 * self.value(), r as int % (self.value() as int) == a as int % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(self.supports_opt);
        proof {
            assert((self.p as int) * (self.p as int) <= u128::MAX) by (nonlinear_arith)
                requires 2 <= self.p < verified::MODULUS_LIMIT;
        }
        debug_assert!(a < (self.p as u128) * (self.p as u128));

        let r = verified::lazy_reduce_opt_u128(
            a,
            self.p,
            self.barrett_hi,
            self.barrett_lo,
            self.leading_zeros,
        );

        debug_assert!((r as u128) < 2 * (self.p as u128));
        debug_assert!(r % self.p == (a % (self.p as u128)) as u64);

        r
    }

    /// Lazy optimized modular reduction of a in constant time.
    /// The output is in the interval [0, 2 * p).
    /// Requires a modulus supporting optimized reduction.
    // Production canonical/vector kernels compose with this helper in Verus.
    #[cfg(test)]
    pub(super) const fn lazy_reduce_opt(&self, a: u64) -> (r: u64)
        requires self.optimized(),
        ensures r < 2 * self.value(), r as int % (self.value() as int) == a as int % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        debug_assert!(self.supports_opt);
        let r = verified::lazy_reduce_opt(a, self.p, self.leading_zeros);

        debug_assert!((r as u128) < 2 * (self.p as u128));
        debug_assert!(r % self.p == a % self.p);

        r
    }

    /// Lazy modular reduction of a vector in constant time.
    /// The output coefficients are in the interval [0, 2 * p).
    pub fn lazy_reduce_vec(&self, a: &mut [u64])
        ensures
            final(a)@.len() == old(a)@.len(),
            forall|i: int| 0 <= i < final(a)@.len() ==>
                #[trigger] final(a)@[i] < 2 * self.value()
                && final(a)@[i] as int % (self.value() as int) == old(a)@[i] as int % (self.value() as int),
    {
        proof {
            use_type_invariant(self);
        }
        if self.supports_opt {
            verified::lazy_reduce_vec::<true>(a, self.p, self.barrett_hi, self.barrett_lo, self.leading_zeros);
        } else {
            verified::lazy_reduce_vec::<false>(a, self.p, self.barrett_hi, self.barrett_lo, self.leading_zeros);
        }
    }

    /// Returns `size` uniformly sampled residues in `[0, p)`.
    pub fn random_vec<R: RngCore + CryptoRng>(&self, size: usize, rng: &mut R) -> (r: Vec<u64>)
        ensures r@.len() == size,
            forall|i: int| 0 <= i < r@.len() ==> #[trigger] r@[i] < self.value(),
    {
        proof {
            use_type_invariant(self);
        }
        let mut sampler = super::backend::Sampler::new(&self.backend, rng);
        verified::random_vec(size, &mut sampler)
    }

    /// Fills `out` with uniformly sampled residues in `[0, p)`.
    pub(crate) fn random_into<R: RngCore + CryptoRng>(&self, out: &mut [u64], rng: &mut R)
        ensures final(out)@.len() == old(out)@.len(),
            forall|i: int| 0 <= i < final(out)@.len() ==> #[trigger] final(out)@[i] < self.value(),
    {
        proof {
            use_type_invariant(self);
        }
        let mut sampler = super::backend::Sampler::new(&self.backend, rng);
        verified::random_into(out, &mut sampler);
    }

    pub(crate) fn serialize_into(&self, a: &[u64], out: &mut Vec<u8>)
        requires
            packing::valid_words(a@, serialization::width(self.value())),
            a.len() * serialization::width(self.value()) <= usize::MAX,
            old(out).len() + (a.len() * serialization::width(self.value()) + 7) / 8 <= usize::MAX,
        ensures
            final(out)@.take(old(out).len() as int) == old(out)@,
            packing::packed(a@, serialization::width(self.value()), final(out)@.skip(old(out).len() as int), 8),
    {
        proof {
            use_type_invariant(self);
        }
        verified::serialize_into(a, self.p, out);
    }

    /// Length of the serialization of a vector of size `size`.
    ///
    /// Panics if the size is not a multiple of 8 or the byte count overflows
    /// `usize`.
    #[must_use]
    #[expect(
        clippy::expect_used,
        reason = "the documented API panics for invalid lengths"
    )]
    pub const fn serialization_length(&self, size: usize) -> (r: usize)
        requires size % 8 == 0, size * serialization::width(self.value()) / 8 <= usize::MAX,
        ensures r == size * serialization::width(self.value()) / 8,
    {
        proof {
            use_type_invariant(self);
        }
        assert!(size.is_multiple_of(8));
        verified::serialization_length(self.p, size).expect("serialization length overflows usize")
    }

    /// Pack words at this modulus' bit width, rejecting truncation and padding
    /// the final byte with zeros. This does not validate that every word is <
    /// p.
    pub fn serialize_vec(
        &self,
        a: &[u64],
    ) -> (r: std::result::Result<Vec<u8>, fhe_util::TranscodeError>)
        ensures
            match r {
                Ok(b) => packing::pack_error(a@, serialization::width(self.value()), 0).is_none()
                    && packing::packed(a@, serialization::width(self.value()), b@, 8),
                Err(e) => packing::pack_error(a@, serialization::width(self.value()), 0) == Some(e),
            },
    {
        proof {
            use_type_invariant(self);
        }
        verified::serialize_vec(a, self.p)
    }

    /// Decode every bit at this modulus' width, zero-extending a final partial
    /// word. This does not reduce or validate residues against p.
    pub fn deserialize_vec(
        &self,
        b: &[u8],
    ) -> (r: std::result::Result<Vec<u64>, fhe_util::TranscodeError>)
        ensures
            match r {
                Ok(a) => packing::unpack_error(b@, serialization::width(self.value())).is_none()
                    && packing::packed(b@, 8, a@, serialization::width(self.value())),
                Err(e) => packing::unpack_error(b@, serialization::width(self.value())) == Some(e),
            },
    {
        proof {
            use_type_invariant(self);
        }
        verified::deserialize_vec(b, self.p)
    }
}

}

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
