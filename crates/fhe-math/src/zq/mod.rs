#![warn(missing_docs, unused_imports)]

//! Ring operations for moduli up to 62 bits.

// Expect indexing/slicing in this performance-critical polynomial arithmetic module.
#![expect(
    clippy::indexing_slicing,
    reason = "performance or example code relies on validated indices"
)]

pub mod primes;
pub(crate) mod verified;

use std::ops::Deref;

use crate::error::{Error, Result};
use fhe_util::is_prime;
use pulp::Arch;
use rand::{CryptoRng, Rng as RngCore, RngExt, distr::Uniform};
#[cfg(test)]
use verified::const_time_cond_select;

/// Structure encapsulating an integer modulus up to 62 bits.
#[derive(Debug, Clone)]
pub struct Modulus {
    pub(crate) p: u64,
    barrett_hi: u64,
    barrett_lo: u64,
    leading_zeros: u32,
    pub(crate) supports_opt: bool,
    distribution: Uniform<u64>,
    arch: Arch,
}

// Trusted sampling boundary: Modulus::new constructs Uniform(0, p), and each
// call returns one sample from that same distribution. The verified loops use
// the ResidueSampler contract; rand and this adapter are outside the proof.
struct UniformSampler<'a, R> {
    distribution: Uniform<u64>,
    rng: &'a mut R,
}

// This adapter remains the existing trusted sampling boundary. Give Cargo's
// ghost-enabled Rust build the trait's ghost methods without claiming to
// verify rand or this implementation.
vstd::prelude::verus! {
    #[verifier::external]
    impl<R: RngCore + CryptoRng> verified::ResidueSampler for UniformSampler<'_, R> {
        uninterp spec fn modulus(&self) -> u64;
        uninterp spec fn history(&self) -> vstd::seq::Seq<u64>;

        #[inline]
        fn sample(&mut self) -> u64 {
            self.rng.sample(self.distribution)
        }
    }
}

// We need to declare Eq manually because of the `Uniform` member.
impl Eq for Modulus {}

impl PartialEq for Modulus {
    fn eq(&self, other: &Self) -> bool {
        let Self {
            p,
            barrett_hi: _,
            barrett_lo: _,
            leading_zeros: _,
            supports_opt: _,
            distribution: _,
            arch: _,
        } = self;
        let Self {
            p: other_p,
            barrett_hi: _,
            barrett_lo: _,
            leading_zeros: _,
            supports_opt: _,
            distribution: _,
            arch: _,
        } = other;

        // All other fields are deterministically derived from p, so we only compare p.
        // The destructuring ensures the compiler will warn us if new fields are added.
        p == other_p
    }
}

// Override the dereference to return the underlying modulus.
impl Deref for Modulus {
    type Target = u64;

    fn deref(&self) -> &Self::Target {
        &self.p
    }
}

impl Modulus {
    /// Create a modulus from an integer of at most 62 bits.
    pub fn new(p: u64) -> Result<Self> {
        let verified::ModulusParameters {
            p,
            barrett_hi,
            barrett_lo,
            leading_zeros,
            supports_opt,
        } = verified::new_modulus(p).ok_or(Error::InvalidModulus(p))?;
        Ok(Self {
            p,
            barrett_hi,
            barrett_lo,
            leading_zeros,
            supports_opt,
            distribution: Uniform::new(0, p).unwrap(),
            arch: Arch::new(),
        })
    }

    /// Performs the modular addition of a and b in constant time.
    /// Aborts if a >= p or b >= p in debug mode.
    #[must_use]
    pub const fn add(&self, a: u64, b: u64) -> u64 {
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
    pub const unsafe fn add_vt(&self, a: u64, b: u64) -> u64 {
        debug_assert!(a < self.p && b < self.p);
        unsafe { verified::add_mod_vt(a, b, self.p) }
    }

    /// Performs the modular subtraction of a and b in constant time.
    /// Aborts if a >= p or b >= p in debug mode.
    #[must_use]
    pub const fn sub(&self, a: u64, b: u64) -> u64 {
        debug_assert!(a < self.p && b < self.p);
        verified::sub_mod(a, b, self.p)
    }

    /// Performs the modular multiplication of a and b in constant time.
    /// Aborts if a >= p or b >= p in debug mode.
    #[must_use]
    pub const fn mul(&self, a: u64, b: u64) -> u64 {
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
    const unsafe fn mul_vt(&self, a: u64, b: u64) -> u64 {
        debug_assert!(a < self.p && b < self.p);
        unsafe { verified::mul_mod_vt(a, b, self.p, self.barrett_hi, self.barrett_lo) }
    }

    /// Optimized modular multiplication of a and b in constant time.
    ///
    /// Aborts if a >= p or b >= p in debug mode.
    #[must_use]
    pub const fn mul_opt(&self, a: u64, b: u64) -> u64 {
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
    const unsafe fn mul_opt_vt(&self, a: u64, b: u64) -> u64 {
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
    pub const fn neg(&self, a: u64) -> u64 {
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
    const unsafe fn neg_vt(&self, a: u64) -> u64 {
        debug_assert!(a < self.p);
        unsafe { verified::neg_mod_vt(a, self.p) }
    }

    /// Compute the Shoup representation of a.
    ///
    /// Aborts if a >= p in debug mode.
    #[must_use]
    pub const fn shoup(&self, a: u64) -> u64 {
        debug_assert!(a < self.p);
        verified::shoup(a, self.p)
    }

    /// Shoup multiplication of a and b in constant time.
    /// The first operand a may be any u64; b must be canonical.
    ///
    /// Aborts if b >= p or b_shoup != shoup(b) in debug mode.
    #[must_use]
    pub const fn mul_shoup(&self, a: u64, b: u64, b_shoup: u64) -> u64 {
        debug_assert!(b < self.p);
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
    const unsafe fn mul_shoup_vt(&self, a: u64, b: u64, b_shoup: u64) -> u64 {
        debug_assert!(b < self.p);
        debug_assert!(b_shoup == self.shoup(b));
        unsafe { verified::mul_shoup_vt(a, b, b_shoup, self.p) }
    }

    /// Lazy Shoup multiplication of a and b in constant time.
    /// The output is in the interval [0, 2 * p).
    /// The first operand a may be any u64; b must be canonical.
    ///
    /// Aborts if b >= p or b_shoup != shoup(b) in debug mode.
    #[must_use]
    pub const fn lazy_mul_shoup(&self, a: u64, b: u64, b_shoup: u64) -> u64 {
        debug_assert!(b < self.p);
        debug_assert!(b_shoup == self.shoup(b));

        let r = verified::lazy_mul_shoup(a, b, b_shoup, self.p);
        debug_assert!(r < 2 * self.p);

        r
    }

    /// Modular addition of vectors in place in constant time.
    ///
    /// Aborts if a and b differ in size, and if any of their values is >= p in
    /// debug mode.
    pub fn add_vec(&self, a: &mut [u64], b: &[u64]) {
        debug_assert_eq!(a.len(), b.len());
        debug_assert!(a.iter().all(|ai| *ai < self.p));
        debug_assert!(b.iter().all(|bi| *bi < self.p));
        self.arch.dispatch(|| verified::add_vec(a, b, self.p))
    }

    /// Modular addition of vectors in place in variable time.
    /// Aborts if a and b differ in size, and if any of their values is >= p in
    /// debug mode.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being added.
    pub unsafe fn add_vec_vt(&self, a: &mut [u64], b: &[u64]) {
        debug_assert_eq!(a.len(), b.len());
        debug_assert!(a.iter().all(|ai| *ai < self.p));
        debug_assert!(b.iter().all(|bi| *bi < self.p));
        self.arch
            .dispatch(|| unsafe { verified::add_vec_vt(a, b, self.p) })
    }

    /// Modular subtraction of vectors in place in constant time.
    ///
    /// Aborts if a and b differ in size, and if any of their values is >= p in
    /// debug mode.
    pub fn sub_vec(&self, a: &mut [u64], b: &[u64]) {
        debug_assert_eq!(a.len(), b.len());
        debug_assert!(a.iter().all(|ai| *ai < self.p));
        debug_assert!(b.iter().all(|bi| *bi < self.p));
        self.arch.dispatch(|| verified::sub_vec(a, b, self.p))
    }

    /// Modular subtraction of vectors in place in variable time.
    /// Aborts if a and b differ in size, and if any of their values is >= p in
    /// debug mode.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being subtracted.
    pub unsafe fn sub_vec_vt(&self, a: &mut [u64], b: &[u64]) {
        debug_assert_eq!(a.len(), b.len());
        debug_assert!(a.iter().all(|ai| *ai < self.p));
        debug_assert!(b.iter().all(|bi| *bi < self.p));
        self.arch
            .dispatch(|| unsafe { verified::sub_vec_vt(a, b, self.p) })
    }

    /// Modular multiplication of vectors in place in constant time.
    ///
    /// Aborts if a and b differ in size, and if any of their values is >= p in
    /// debug mode.
    pub fn mul_vec(&self, a: &mut [u64], b: &[u64]) {
        debug_assert_eq!(a.len(), b.len());
        debug_assert!(a.iter().all(|ai| *ai < self.p));
        debug_assert!(b.iter().all(|bi| *bi < self.p));
        if self.supports_opt {
            self.arch.dispatch(|| {
                verified::mul_vec::<true>(
                    a,
                    b,
                    self.p,
                    self.barrett_hi,
                    self.barrett_lo,
                    self.leading_zeros,
                )
            })
        } else {
            self.arch.dispatch(|| {
                verified::mul_vec::<false>(
                    a,
                    b,
                    self.p,
                    self.barrett_hi,
                    self.barrett_lo,
                    self.leading_zeros,
                )
            })
        }
    }

    /// Modular scalar multiplication of vectors in place in constant time.
    /// Entries of a may be any u64; b must be canonical.
    ///
    /// Aborts if b >= p in debug mode.
    pub fn scalar_mul_vec(&self, a: &mut [u64], b: u64) {
        debug_assert!(b < self.p);
        self.arch
            .dispatch(|| verified::scalar_mul_vec(a, b, self.p))
    }

    /// Modular scalar multiplication of vectors in place in variable time.
    /// Entries of a may be any u64; b must be canonical.
    /// Aborts if b >= p in debug mode.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being multiplied.
    pub unsafe fn scalar_mul_vec_vt(&self, a: &mut [u64], b: u64) {
        debug_assert!(b < self.p);
        self.arch
            .dispatch(|| unsafe { verified::scalar_mul_vec_vt(a, b, self.p) })
    }

    /// Modular multiplication of vectors in place in variable time.
    /// Aborts if a and b differ in size, and if any of their values is >= p in
    /// debug mode.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being multiplied.
    pub unsafe fn mul_vec_vt(&self, a: &mut [u64], b: &[u64]) {
        debug_assert_eq!(a.len(), b.len());
        debug_assert!(a.iter().all(|ai| *ai < self.p));
        debug_assert!(b.iter().all(|bi| *bi < self.p));
        if self.supports_opt {
            self.arch.dispatch(|| unsafe {
                verified::mul_vec_vt::<true>(
                    a,
                    b,
                    self.p,
                    self.barrett_hi,
                    self.barrett_lo,
                    self.leading_zeros,
                )
            })
        } else {
            self.arch.dispatch(|| unsafe {
                verified::mul_vec_vt::<false>(
                    a,
                    b,
                    self.p,
                    self.barrett_hi,
                    self.barrett_lo,
                    self.leading_zeros,
                )
            })
        }
    }

    /// Compute the Shoup representation of a vector.
    ///
    /// Aborts if any of the values of the vector is >= p in debug mode.
    #[must_use]
    pub fn shoup_vec(&self, a: &[u64]) -> Vec<u64> {
        debug_assert!(a.iter().all(|ai| *ai < self.p));
        self.arch.dispatch(|| verified::shoup_vec(a, self.p))
    }

    pub(crate) fn shoup_into(&self, a: &[u64], out: &mut [u64]) {
        assert_eq!(a.len(), out.len());
        debug_assert!(a.iter().all(|ai| *ai < self.p));
        self.arch.dispatch(|| verified::shoup_into(a, out, self.p));
    }

    /// Shoup modular multiplication of vectors in place in constant time.
    /// Entries of a may be any u64; b must contain canonical residues.
    ///
    /// Aborts in debug mode if lengths differ, an entry of b is >= p, or a
    /// precomputed word differs from shoup(b).
    pub fn mul_shoup_vec(&self, a: &mut [u64], b: &[u64], b_shoup: &[u64]) {
        debug_assert_eq!(a.len(), b.len());
        debug_assert_eq!(a.len(), b_shoup.len());
        debug_assert_eq!(&b_shoup, &self.shoup_vec(b));
        self.arch
            .dispatch(|| verified::mul_shoup_vec(a, b, b_shoup, self.p))
    }

    /// Shoup modular multiplication of vectors in place in variable time.
    /// Entries of a may be any u64; b must contain canonical residues.
    /// Aborts in debug mode if lengths differ, an entry of b is >= p, or a
    /// precomputed word differs from shoup(b).
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being multiplied.
    pub unsafe fn mul_shoup_vec_vt(&self, a: &mut [u64], b: &[u64], b_shoup: &[u64]) {
        debug_assert_eq!(a.len(), b.len());
        debug_assert_eq!(a.len(), b_shoup.len());
        debug_assert_eq!(&b_shoup, &self.shoup_vec(b));
        self.arch
            .dispatch(|| unsafe { verified::mul_shoup_vec_vt(a, b, b_shoup, self.p) })
    }

    /// Reduce a vector in place in constant time.
    pub fn reduce_vec(&self, a: &mut [u64]) {
        self.arch
            .dispatch(|| verified::reduce_vec(a, self.p, self.barrett_hi, self.barrett_lo))
    }

    /// Center a value modulo p as i64 in constant time.
    ///
    /// The output is congruent to a modulo p and satisfies `-p <= 2*result <
    /// p`. For even p, the midpoint p/2 maps to -p/2.
    /// Aborts if `a >= p` in debug mode.
    #[must_use]
    pub const fn center(&self, a: u64) -> i64 {
        debug_assert!(a < self.p);
        verified::center(a, self.p)
    }

    /// Return centered representatives in input order, in constant time.
    /// Each output satisfies the range and congruence guarantees of
    /// [`Self::center`]. Aborts if any input is >= p in debug mode.
    #[must_use]
    pub fn center_vec(&self, a: &[u64]) -> Vec<i64> {
        debug_assert!(a.iter().all(|ai| *ai < self.p));
        self.arch.dispatch(|| verified::center_vec(a, self.p))
    }

    /// Reduce a vector in place in variable time.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being reduced.
    pub unsafe fn reduce_vec_vt(&self, a: &mut [u64]) {
        self.arch.dispatch(|| unsafe {
            verified::reduce_vec_vt(a, self.p, self.barrett_hi, self.barrett_lo)
        })
    }

    /// Modular reduction of a i64 in constant time.
    // Production vector loops call the verified scalar kernel directly.
    #[cfg(test)]
    const fn reduce_i64(&self, a: i64) -> u64 {
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
    const unsafe fn reduce_i64_vt(&self, a: i64) -> u64 {
        let r = unsafe { verified::reduce_i64_vt(a, self.p, self.barrett_hi, self.barrett_lo) };
        debug_assert!(r < self.p && r as i128 == (a as i128).rem_euclid(self.p as i128));
        r
    }

    /// Return a new vector of canonical residues of signed inputs in constant
    /// time.
    #[must_use]
    pub fn reduce_vec_i64(&self, a: &[i64]) -> Vec<u64> {
        self.arch
            .dispatch(|| verified::reduce_vec_i64(a, self.p, self.barrett_hi, self.barrett_lo))
    }

    /// Return a new vector of canonical residues of signed inputs in variable
    /// time.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being reduced.
    #[must_use]
    pub unsafe fn reduce_vec_i64_vt(&self, a: &[i64]) -> Vec<u64> {
        self.arch.dispatch(|| unsafe {
            verified::reduce_vec_i64_vt(a, self.p, self.barrett_hi, self.barrett_lo)
        })
    }

    /// Reduce a vector in constant time.
    #[must_use]
    pub fn reduce_vec_new(&self, a: &[u64]) -> Vec<u64> {
        self.arch
            .dispatch(|| verified::reduce_vec_new(a, self.p, self.barrett_hi, self.barrett_lo))
    }

    /// Reduce a vector in variable time.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being reduced.
    #[must_use]
    pub unsafe fn reduce_vec_new_vt(&self, a: &[u64]) -> Vec<u64> {
        self.arch.dispatch(|| unsafe {
            verified::reduce_vec_new_vt(a, self.p, self.barrett_hi, self.barrett_lo)
        })
    }

    /// Modular negation of a vector in place in constant time.
    ///
    /// Aborts if any of the values in the vector is >= p in debug mode.
    pub fn neg_vec(&self, a: &mut [u64]) {
        debug_assert!(a.iter().all(|ai| *ai < self.p));
        self.arch.dispatch(|| verified::neg_vec(a, self.p))
    }

    /// Modular negation of a vector in place in variable time.
    /// Aborts if any of the values in the vector is >= p in debug mode.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being negated.
    pub unsafe fn neg_vec_vt(&self, a: &mut [u64]) {
        debug_assert!(a.iter().all(|ai| *ai < self.p));
        self.arch
            .dispatch(|| unsafe { verified::neg_vec_vt(a, self.p) })
    }

    /// Modular exponentiation in variable time.
    ///
    /// Returns the canonical residue of a^n modulo p. Exponent zero returns 1,
    /// including when a is zero.
    /// Aborts if a >= p or n >= p in debug mode.
    #[must_use]
    pub fn pow(&self, a: u64, n: u64) -> u64 {
        debug_assert!(a < self.p && n < self.p);

        verified::pow_mod(a, n, self.p, self.barrett_hi, self.barrett_lo)
    }

    /// Modular inversion in variable time.
    ///
    /// Returns the canonical multiplicative inverse of a, or None if p is not
    /// prime or a = 0. Composite moduli are rejected even when a is invertible.
    /// Aborts if p is prime and a >= p in debug mode.
    #[must_use]
    pub fn inv(&self, a: u64) -> std::option::Option<u64> {
        let p_is_prime = is_prime(self.p);
        debug_assert!(!p_is_prime || a < self.p);
        let result = verified::inv_mod(a, self.p, self.barrett_hi, self.barrett_lo, p_is_prime);
        if let Some(r) = result {
            debug_assert_eq!(self.mul(a, r), 1);
        }
        result
    }

    /// Modular reduction of a u128 in constant time.
    #[must_use]
    pub const fn reduce_u128(&self, a: u128) -> u64 {
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
    pub const unsafe fn reduce_u128_vt(&self, a: u128) -> u64 {
        let r = unsafe { verified::reduce_u128_vt(a, self.p, self.barrett_hi, self.barrett_lo) };
        debug_assert!(r < self.p && r as u128 == a % self.p as u128);
        r
    }

    /// Modular reduction of a u64 in constant time.
    #[must_use]
    pub const fn reduce(&self, a: u64) -> u64 {
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
    pub const unsafe fn reduce_vt(&self, a: u64) -> u64 {
        let r = unsafe { verified::reduce_vt(a, self.p, self.barrett_hi, self.barrett_lo) };
        debug_assert!(r < self.p && r == a % self.p);
        r
    }

    /// Optimized modular reduction of a u128 in constant time.
    /// Requires an eligible modulus and `a < p^2`; checked in debug mode.
    #[must_use]
    pub const fn reduce_opt_u128(&self, a: u128) -> u64 {
        debug_assert!(self.supports_opt);
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
    pub(crate) const unsafe fn reduce_opt_u128_vt(&self, a: u128) -> u64 {
        debug_assert!(self.supports_opt);
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
    pub const fn reduce_opt(&self, a: u64) -> u64 {
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
    pub const unsafe fn reduce_opt_vt(&self, a: u64) -> u64 {
        debug_assert!(self.supports_opt);
        let r = unsafe { verified::reduce_opt_vt(a, self.p, self.leading_zeros) };
        debug_assert!(r < self.p && r == a % self.p);
        r
    }

    /// Return x mod p in constant time.
    /// Aborts if x >= 2 * p in debug mode.
    pub(crate) const fn reduce1(x: u64, p: u64) -> u64 {
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
    pub(crate) const unsafe fn reduce1_vt(x: u64, p: u64) -> u64 {
        debug_assert!(p >> 63 == 0);
        debug_assert!(x < 2 * p);

        let r = unsafe { verified::reduce1_vt(x, p) };
        debug_assert!(r == x % p);

        r
    }

    /// Lazy modular reduction of a in constant time.
    /// The output is in the interval [0, 2 * p).
    #[must_use]
    pub const fn lazy_reduce_u128(&self, a: u128) -> u64 {
        let r = verified::lazy_reduce_u128(a, self.p, self.barrett_hi, self.barrett_lo);

        debug_assert!((r as u128) < 2 * (self.p as u128));
        debug_assert!(r % self.p == (a % (self.p as u128)) as u64);

        r
    }

    /// Lazy modular reduction of a in constant time.
    /// The output is in the interval [0, 2 * p).
    #[must_use]
    pub const fn lazy_reduce(&self, a: u64) -> u64 {
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
    pub const fn lazy_reduce_opt_u128(&self, a: u128) -> u64 {
        debug_assert!(self.supports_opt);
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
    const fn lazy_reduce_opt(&self, a: u64) -> u64 {
        debug_assert!(self.supports_opt);
        let r = verified::lazy_reduce_opt(a, self.p, self.leading_zeros);

        debug_assert!((r as u128) < 2 * (self.p as u128));
        debug_assert!(r % self.p == a % self.p);

        r
    }

    /// Lazy modular reduction of a vector in constant time.
    /// The output coefficients are in the interval [0, 2 * p).
    pub fn lazy_reduce_vec(&self, a: &mut [u64]) {
        if self.supports_opt {
            verified::lazy_reduce_vec::<true>(
                a,
                self.p,
                self.barrett_hi,
                self.barrett_lo,
                self.leading_zeros,
            );
        } else {
            verified::lazy_reduce_vec::<false>(
                a,
                self.p,
                self.barrett_hi,
                self.barrett_lo,
                self.leading_zeros,
            );
        }
    }

    /// Returns `size` uniformly sampled residues in `[0, p)`.
    pub fn random_vec<R: RngCore + CryptoRng>(&self, size: usize, rng: &mut R) -> Vec<u64> {
        verified::random_vec(
            size,
            &mut UniformSampler {
                distribution: self.distribution,
                rng,
            },
        )
    }

    /// Fills `out` with uniformly sampled residues in `[0, p)`.
    pub(crate) fn random_into<R: RngCore + CryptoRng>(&self, out: &mut [u64], rng: &mut R) {
        verified::random_into(
            out,
            &mut UniformSampler {
                distribution: self.distribution,
                rng,
            },
        );
    }

    pub(crate) fn serialize_into(&self, a: &[u64], out: &mut Vec<u8>) {
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
    pub const fn serialization_length(&self, size: usize) -> usize {
        assert!(size.is_multiple_of(8));
        verified::serialization_length(self.p, size).expect("serialization length overflows usize")
    }

    /// Pack words at this modulus' bit width, rejecting truncation and padding
    /// the final byte with zeros. This does not validate that every word is <
    /// p.
    pub fn serialize_vec(
        &self,
        a: &[u64],
    ) -> std::result::Result<Vec<u8>, fhe_util::TranscodeError> {
        verified::serialize_vec(a, self.p)
    }

    /// Decode every bit at this modulus' width, zero-extending a final partial
    /// word. This does not reduce or validate residues against p.
    pub fn deserialize_vec(
        &self,
        b: &[u8],
    ) -> std::result::Result<Vec<u64>, fhe_util::TranscodeError> {
        verified::deserialize_vec(b, self.p)
    }
}

#[cfg(test)]
mod tests {
    use super::Modulus;
    use itertools::{Itertools, izip};
    use proptest::collection::vec as prop_vec;
    use proptest::prelude::{BoxedStrategy, Just, Strategy, any};
    use rand::{RngExt, SeedableRng, distr::Uniform};
    use rand_chacha::ChaCha8Rng;

    // Independent width oracle: find the smallest capacity containing [0, p).
    fn reference_serialization_width(p: u64) -> usize {
        let mut capacity = 1u128;
        let mut width = 0;
        while capacity < p as u128 {
            capacity *= 2;
            width += 1;
        }
        width
    }

    fn reference_serialization(words: &[u64], width: usize) -> Vec<u8> {
        let mut bytes = Vec::new();
        let mut byte = 0;
        let mut offset = 0;
        for word in words {
            for bit in 0..width {
                byte |= (((word >> bit) & 1) as u8) << offset;
                offset += 1;
                if offset == 8 {
                    bytes.push(byte);
                    byte = 0;
                    offset = 0;
                }
            }
        }
        if offset != 0 {
            bytes.push(byte);
        }
        bytes
    }

    fn reference_deserialization(bytes: &[u8], width: usize) -> Vec<u64> {
        let mut words = Vec::new();
        let mut word = 0;
        let mut offset = 0;
        for byte in bytes {
            for bit in 0..8 {
                word |= (((byte >> bit) & 1) as u64) << offset;
                offset += 1;
                if offset == width {
                    words.push(word);
                    word = 0;
                    offset = 0;
                }
            }
        }
        if offset != 0 {
            words.push(word);
        }
        words
    }

    #[test]
    fn serialization_covers_modulus_widths_and_partial_bytes() {
        for p in reduction_boundary_moduli() {
            let modulus = Modulus::new(p).unwrap();
            let width = reference_serialization_width(p);
            assert_eq!(modulus.serialization_length(8), width);
            for len in 0usize..=17 {
                let words: Vec<_> = (0..len)
                    .map(|i| match i % 3 {
                        0 => 0,
                        1 => 1,
                        _ => p - 1,
                    })
                    .collect();
                let bytes = modulus.serialize_vec(&words).unwrap();
                assert_eq!(
                    bytes,
                    reference_serialization(&words, width),
                    "p={p}, len={len}"
                );
                assert_eq!(bytes.len(), (len * width).div_ceil(8));
                let decoded = modulus.deserialize_vec(&bytes).unwrap();
                assert_eq!(decoded, reference_deserialization(&bytes, width));
                assert_eq!(&decoded[..len], words);
                assert!(decoded[len..].iter().all(|x| *x == 0));
                if len.is_multiple_of(8) {
                    assert_eq!(bytes.len(), modulus.serialization_length(len));
                    assert_eq!(decoded, words);
                }
                let mut appended = vec![0xa5, 0x5a];
                modulus.serialize_into(&words, &mut appended);
                assert_eq!(&appended[..2], [0xa5, 0x5a]);
                assert_eq!(&appended[2..], bytes);
            }
            for bytes in [&[][..], &[0xff][..], &[0xff, 0x80, 0, 0x55][..]] {
                assert_eq!(
                    modulus.deserialize_vec(bytes).unwrap(),
                    reference_deserialization(bytes, width)
                );
            }
        }
    }

    #[test]
    fn serialization_preserves_width_checks_and_stream_semantics() {
        use fhe_util::TranscodeError;
        use std::panic::{AssertUnwindSafe, catch_unwind};

        let modulus = Modulus::new(5).unwrap();
        // Fitting words remain accepted even when they are not canonical residues.
        assert_eq!(modulus.serialize_vec(&[5, 7]).unwrap(), [61]);
        assert_eq!(modulus.deserialize_vec(&[61]).unwrap(), [5, 7, 0]);
        assert_eq!(modulus.deserialize_vec(&[0xff]).unwrap(), [7, 7, 3]);
        assert_eq!(
            modulus.serialize_vec(&[8]),
            Err(TranscodeError::ValueTooWide)
        );
        assert_eq!(
            modulus.serialize_vec(&[u64::MAX]),
            Err(TranscodeError::ValueTooWide)
        );

        let mut bytes = vec![0xa5];
        modulus.serialize_into(&[1], &mut bytes);
        modulus.serialize_into(&[2], &mut bytes);
        assert_eq!(bytes, [0xa5, 1, 2]);
        let pointer = bytes.as_ptr();
        let capacity = bytes.capacity();
        assert!(
            catch_unwind(AssertUnwindSafe(
                || modulus.serialize_into(&[0, 8], &mut bytes)
            ))
            .is_err()
        );
        assert_eq!(bytes, [0xa5, 1, 2]);
        assert_eq!(bytes.as_ptr(), pointer);
        assert_eq!(bytes.capacity(), capacity);
    }

    #[test]
    fn serialization_length_checks_machine_boundaries_and_remains_const() {
        const LARGE: Option<usize> = super::verified::serialization_length(3, usize::MAX - 7);
        const INVALID: Option<usize> = super::verified::serialization_length(3, 7);
        assert_eq!(LARGE, Some((usize::MAX - 7) / 4));
        assert_eq!(INVALID, None);

        for p in reduction_boundary_moduli() {
            let modulus = Modulus::new(p).unwrap();
            let width = reference_serialization_width(p);
            let max_size = (usize::MAX / width).min(usize::MAX / 8) * 8;
            for size in [0, 8, 16, max_size] {
                let expected = size as u128 * width as u128 / 8;
                assert_eq!(modulus.serialization_length(size) as u128, expected);
            }
            if max_size <= usize::MAX - 8 {
                let invalid = max_size + 8;
                assert_eq!(super::verified::serialization_length(p, invalid), None);
                // This check is required in release builds too.
                assert!(
                    std::panic::catch_unwind(|| modulus.serialization_length(invalid)).is_err()
                );
            }
            for invalid in [1, 7, 9, usize::MAX] {
                assert_eq!(super::verified::serialization_length(p, invalid), None);
                assert!(
                    std::panic::catch_unwind(|| modulus.serialization_length(invalid)).is_err()
                );
            }
        }
    }

    fn check_constructor_metadata(p: u64) {
        use num_bigint::BigUint;
        use num_traits::ToPrimitive;

        let modulus = Modulus::new(p).unwrap();
        let expected_reciprocal = ((BigUint::from(1u64) << 128usize) / p).to_u128().unwrap();
        let actual_reciprocal = ((modulus.barrett_hi as u128) << 64) | modulus.barrett_lo as u128;
        assert_eq!(modulus.p, p);
        assert_eq!(actual_reciprocal, expected_reciprocal, "p={p}");
        assert_eq!(modulus.leading_zeros, p.leading_zeros(), "p={p}");

        // Evaluate the original NFL inequality independently of supports_opt.
        let t = BigUint::from(1u64) << p.leading_zeros() as usize;
        let cube = &t * &t * &t;
        let expected_opt = ((&cube + 1u64) << 64usize) < &cube * (&t + 1u64) * p;
        assert_eq!(modulus.supports_opt, expected_opt, "p={p}");
    }

    #[test]
    fn constructor_exhausts_small_moduli() {
        for p in 2..=u16::MAX as u64 {
            check_constructor_metadata(p);
        }
    }

    #[test]
    fn constructor_checks_boundaries_and_optimization_thresholds() {
        // Keep an independent numeric check of the shared, exclusive limit.
        assert_eq!(super::verified::MODULUS_LIMIT, 1u64 << 62);
        check_constructor_metadata(super::verified::MODULUS_LIMIT - 1);
        // Includes every supported power of two (where the reciprocal needs
        // its +1 correction), adjacent values, and optimization thresholds.
        for p in reduction_boundary_moduli() {
            check_constructor_metadata(p);
        }
        for p in [
            0,
            1,
            super::verified::MODULUS_LIMIT,
            super::verified::MODULUS_LIMIT + 1,
            (1 << 63) - 1,
            1 << 63,
            u64::MAX,
        ] {
            assert_eq!(Modulus::new(p), Err(crate::Error::InvalidModulus(p)));
        }
    }

    #[test]
    fn conditional_selection_handles_word_patterns_and_remains_const() {
        const SELECTED: u64 = super::const_time_cond_select(u64::MAX, 0, true);
        assert_eq!(SELECTED, u64::MAX);

        let words = [
            0,
            1,
            (1 << 63) - 1,
            1 << 63,
            u64::MAX,
            0x5555_5555_5555_5555,
            0xaaaa_aaaa_aaaa_aaaa,
        ];
        for on_true in words {
            for on_false in words {
                for cond in [false, true] {
                    assert_eq!(
                        super::const_time_cond_select(on_true, on_false, cond),
                        if cond { on_true } else { on_false },
                    );
                }
            }
        }
    }

    #[test]
    fn reduce1_exhausts_small_moduli() {
        for p in 1..=64 {
            for x in 0..2 * p {
                for result in [Modulus::reduce1(x, p), unsafe { Modulus::reduce1_vt(x, p) }] {
                    assert!(result < p);
                    assert_eq!(result, x % p);
                }
            }
        }
    }

    #[test]
    fn reduce1_handles_wide_moduli_and_remains_const() {
        const REDUCED: u64 = super::verified::reduce1(u64::MAX - 2, (1 << 63) - 1);
        const REDUCED_VT: u64 = unsafe { super::verified::reduce1_vt(u64::MAX - 2, (1 << 63) - 1) };
        assert_eq!(REDUCED, (1 << 63) - 2);
        assert_eq!(REDUCED_VT, (1 << 63) - 2);

        // Include the wider helper domain used by the NTT with p = 2*q.
        for p in [1, 2, 3, (1 << 62) - 1, 1 << 62, (1 << 63) - 1] {
            for x in [0, p - 1, p, 2 * p - 1] {
                for result in [Modulus::reduce1(x, p), unsafe { Modulus::reduce1_vt(x, p) }] {
                    assert!(result < p);
                    assert_eq!(result, x % p);
                }
            }
        }
    }

    #[test]
    fn add_exhausts_small_moduli() {
        // Include even and composite moduli: addition does not require a prime.
        for p in 2..=64 {
            let modulus = Modulus::new(p).unwrap();
            for a in 0..p {
                for b in 0..p {
                    for result in [modulus.add(a, b), unsafe { modulus.add_vt(a, b) }] {
                        assert!(result < p);
                        assert_eq!(result as u128, (a as u128 + b as u128) % p as u128);
                    }
                }
            }
        }
    }

    #[test]
    fn add_handles_boundaries_and_remains_const() {
        const SUM: u64 = super::verified::add_mod(4, 6, 7);
        const SUM_VT: u64 = unsafe { super::verified::add_mod_vt(4, 6, 7) };
        assert_eq!(SUM, 3);
        assert_eq!(SUM_VT, 3);

        for p in [
            2,
            3,
            4,
            (1 << 61) - 1,
            1 << 61,
            (1 << 62) - 2,
            (1 << 62) - 1,
        ] {
            let modulus = Modulus::new(p).unwrap();
            for (a, b) in [
                (0, 0),
                (0, p - 1),
                (p - 1, 0),
                (p - 1, 1),
                (p - 1, p - 1),
                (p / 2, p - p / 2),
            ] {
                for result in [modulus.add(a, b), unsafe { modulus.add_vt(a, b) }] {
                    assert!(result < p);
                    assert_eq!(result as u128, (a as u128 + b as u128) % p as u128);
                }
            }
        }
    }

    // Include powers of two, maximal moduli of each width, and both sides of
    // the strict NFL inequality used to select the optimized implementation.
    fn reduction_boundary_moduli() -> Vec<u64> {
        use num_bigint::BigUint;
        use num_traits::ToPrimitive;

        let mut moduli: Vec<u64> = (2..=64).collect();
        for bits in 2usize..=62 {
            let upper = 1u64 << bits;
            moduli.extend([upper / 2, upper / 2 + 1, upper - 1]);
            let t = BigUint::from(1u64) << (64 - bits);
            let cube = &t * &t * &t;
            let threshold = (((&cube + 1u64) << 64usize) / (&cube * (&t + 1u64)))
                .to_u64()
                .unwrap();
            for p in [threshold, threshold + 1] {
                if (2..upper).contains(&p) {
                    moduli.push(p);
                }
            }
        }
        moduli.sort_unstable();
        moduli.dedup();
        moduli
    }

    #[test]
    fn lazy_reducers_cover_boundaries_and_optimization_thresholds() {
        let mut optimized_count = 0;
        let mut general_count = 0;
        for p in reduction_boundary_moduli() {
            let modulus = Modulus::new(p).unwrap();
            let square = (p as u128) * (p as u128);
            let words = [0, 1, p - 1, p, 2 * p - 1, 2 * p, u64::MAX];
            for a in words {
                let r = modulus.lazy_reduce(a);
                assert!(r < 2 * p, "p={p}, a={a}");
                assert_eq!(r % p, a % p, "p={p}, a={a}");
                if modulus.supports_opt {
                    let r = modulus.lazy_reduce_opt(a);
                    assert!(r < 2 * p, "p={p}, a={a}");
                    assert_eq!(r % p, a % p, "p={p}, a={a}");
                }
            }
            for a in [
                0,
                1,
                p as u128 - 1,
                p as u128,
                2 * p as u128 - 1,
                (1u128 << 64) - 1,
                1u128 << 64,
                (1u128 << 64) + 1,
                square / 2,
                square - p as u128,
                square - 1,
                square,
                u128::MAX - p as u128,
                u128::MAX,
            ] {
                let r = modulus.lazy_reduce_u128(a);
                assert!(r < 2 * p, "p={p}, a={a}");
                assert_eq!(r as u128 % p as u128, a % p as u128, "p={p}, a={a}");
                if modulus.supports_opt && a < square {
                    let r = modulus.lazy_reduce_opt_u128(a);
                    assert!(r < 2 * p, "p={p}, a={a}");
                    assert_eq!(r as u128 % p as u128, a % p as u128, "p={p}, a={a}");
                }
            }
            for len in [0, 1, 3, 16, 17, 65] {
                let original: Vec<_> = words.into_iter().cycle().take(len).collect();
                let mut reduced = original.clone();
                modulus.lazy_reduce_vec(&mut reduced);
                assert_eq!(reduced.len(), original.len());
                for (a, r) in original.into_iter().zip(reduced) {
                    assert!(r < 2 * p, "p={p}, a={a}");
                    assert_eq!(r % p, a % p, "p={p}, a={a}");
                }
            }
            if modulus.supports_opt {
                optimized_count += 1;
            } else {
                general_count += 1;
            }
        }
        assert!(optimized_count > 0 && general_count > 0);
    }

    #[test]
    fn lazy_reducers_exhaust_small_inputs() {
        for p in 2..=64 {
            let modulus = Modulus::new(p).unwrap();
            for a in 0..4 * p {
                for r in [modulus.lazy_reduce(a), modulus.lazy_reduce_u128(a as u128)] {
                    assert!(r < 2 * p);
                    assert_eq!(r % p, a % p);
                }
            }
        }
    }

    #[test]
    fn lazy_scalar_kernels_remain_const() {
        const P: u64 = (1 << 62) - 1;
        const MU: u128 = u128::MAX / P as u128;
        const HI: u64 = (MU >> 64) as u64;
        const LO: u64 = MU as u64;
        const A: u128 = P as u128 * P as u128 - 1;
        const GENERAL: u64 = super::verified::lazy_reduce(u64::MAX, P, HI, LO);
        const WIDE: u64 = super::verified::lazy_reduce_u128(u128::MAX, P, HI, LO);
        const OPT: u64 = super::verified::lazy_reduce_opt(u64::MAX, P, P.leading_zeros());
        const OPT_WIDE: u64 =
            super::verified::lazy_reduce_opt_u128(A, P, HI, LO, P.leading_zeros());
        for (a, r) in [
            (u64::MAX as u128, GENERAL),
            (u128::MAX, WIDE),
            (u64::MAX as u128, OPT),
            (A, OPT_WIDE),
        ] {
            assert!(r < 2 * P);
            assert_eq!(r as u128 % P as u128, a % P as u128);
        }
    }

    #[test]
    #[cfg(debug_assertions)]
    fn lazy_optimized_reducers_reject_invalid_inputs() {
        let unsupported = Modulus::new(17).unwrap();
        assert!(std::panic::catch_unwind(|| unsupported.lazy_reduce_opt(0)).is_err());
        assert!(std::panic::catch_unwind(|| unsupported.lazy_reduce_opt_u128(0)).is_err());
        let modulus = Modulus::new((1 << 62) - 1).unwrap();
        let square = (*modulus as u128) * (*modulus as u128);
        assert!(std::panic::catch_unwind(|| modulus.lazy_reduce_opt_u128(square)).is_err());
    }

    #[test]
    fn canonical_reducers_cover_boundaries() {
        for p in reduction_boundary_moduli() {
            let modulus = Modulus::new(p).unwrap();
            let words = [
                0,
                1,
                p - 1,
                p,
                p + 1,
                2 * p - 1,
                2 * p,
                4 * p - 1,
                1 << 63,
                u64::MAX - 1,
                u64::MAX,
            ];
            for a in words {
                assert_eq!(modulus.reduce(a), a % p);
                assert_eq!(unsafe { modulus.reduce_vt(a) }, a % p);
                if modulus.supports_opt {
                    assert_eq!(modulus.reduce_opt(a), a % p);
                    assert_eq!(unsafe { modulus.reduce_opt_vt(a) }, a % p);
                }
            }

            let wide_p = p as u128;
            let square = wide_p * wide_p;
            for a in [
                0,
                1,
                wide_p - 1,
                wide_p,
                wide_p + 1,
                square - 1,
                square,
                square + 1,
                (1u128 << 64) - 1,
                1 << 64,
                1 << 127,
                u128::MAX,
            ] {
                let expected = (a % wide_p) as u64;
                assert_eq!(modulus.reduce_u128(a), expected);
                assert_eq!(unsafe { modulus.reduce_u128_vt(a) }, expected);
                if modulus.supports_opt && a < square {
                    assert_eq!(modulus.reduce_opt_u128(a), expected);
                    assert_eq!(unsafe { modulus.reduce_opt_u128_vt(a) }, expected);
                }
            }

            let signed_p = p as i64;
            let signed = [
                i64::MIN,
                i64::MIN + 1,
                -signed_p - 1,
                -signed_p,
                -signed_p + 1,
                -1,
                0,
                1,
                signed_p - 1,
                signed_p,
                signed_p + 1,
                i64::MAX - 1,
                i64::MAX,
            ];
            for a in signed {
                let expected = (a as i128).rem_euclid(p as i128) as u64;
                assert_eq!(modulus.reduce_i64(a), expected);
                assert_eq!(unsafe { modulus.reduce_i64_vt(a) }, expected);
            }

            for len in [0, 1, 2, 7, 8, 9, 31, 32, 33, 63, 64, 65, 129] {
                let input: Vec<_> = words.into_iter().cycle().take(len).collect();
                let expected: Vec<_> = input.iter().map(|a| a % p).collect();
                // Slice sentinels also check that only the requested region changes.
                let mut buffer = vec![u64::MAX; len + 2];
                buffer[1..=len].copy_from_slice(&input);
                modulus.reduce_vec(&mut buffer[1..=len]);
                assert_eq!(&buffer[1..=len], expected);
                assert_eq!((buffer[0], buffer[len + 1]), (u64::MAX, u64::MAX));
                buffer[1..=len].copy_from_slice(&input);
                unsafe {
                    modulus.reduce_vec_vt(&mut buffer[1..=len]);
                }
                assert_eq!(&buffer[1..=len], expected);
                assert_eq!((buffer[0], buffer[len + 1]), (u64::MAX, u64::MAX));
                assert_eq!(modulus.reduce_vec_new(&input), expected);
                assert_eq!(unsafe { modulus.reduce_vec_new_vt(&input) }, expected);

                let input: Vec<_> = signed.into_iter().cycle().take(len).collect();
                let expected: Vec<_> = input
                    .iter()
                    .map(|a| (*a as i128).rem_euclid(p as i128) as u64)
                    .collect();
                assert_eq!(modulus.reduce_vec_i64(&input), expected);
                assert_eq!(unsafe { modulus.reduce_vec_i64_vt(&input) }, expected);
            }
        }
    }

    #[test]
    fn canonical_reducers_exhaust_small_inputs() {
        for p in 2..=64 {
            let modulus = Modulus::new(p).unwrap();
            for a in 0..4 * p {
                let expected = a % p;
                assert_eq!(modulus.reduce(a), expected);
                assert_eq!(unsafe { modulus.reduce_vt(a) }, expected);
                assert_eq!(modulus.reduce_u128(a as u128), expected);
                assert_eq!(unsafe { modulus.reduce_u128_vt(a as u128) }, expected);
            }
            for a in -(4 * p as i64)..=4 * p as i64 {
                let expected = (a as i128).rem_euclid(p as i128) as u64;
                assert_eq!(modulus.reduce_i64(a), expected);
                assert_eq!(unsafe { modulus.reduce_i64_vt(a) }, expected);
            }
        }
    }

    #[test]
    fn canonical_scalar_kernels_remain_const() {
        use super::verified;
        const P: u64 = (1 << 62) - 1;
        const MU: u128 = u128::MAX / P as u128;
        const HI: u64 = (MU >> 64) as u64;
        const LO: u64 = MU as u64;
        const A: u128 = P as u128 * P as u128 - 1;
        const RESULTS: [u64; 12] = [
            verified::reduce(u64::MAX, P, HI, LO),
            unsafe { verified::reduce_vt(u64::MAX, P, HI, LO) },
            verified::reduce_u128(u128::MAX, P, HI, LO),
            unsafe { verified::reduce_u128_vt(u128::MAX, P, HI, LO) },
            verified::reduce_opt(u64::MAX, P, P.leading_zeros()),
            unsafe { verified::reduce_opt_vt(u64::MAX, P, P.leading_zeros()) },
            verified::reduce_opt_u128(A, P, HI, LO, P.leading_zeros()),
            unsafe { verified::reduce_opt_u128_vt(A, P, HI, LO, P.leading_zeros()) },
            verified::reduce_i64(i64::MIN, P, HI, LO),
            unsafe { verified::reduce_i64_vt(i64::MIN, P, HI, LO) },
            verified::reduce3(4 * P - 1, P, 2 * P),
            unsafe { verified::reduce3_vt(4 * P - 1, P, 2 * P) },
        ];
        let expected = [
            u64::MAX % P,
            (u128::MAX % P as u128) as u64,
            u64::MAX % P,
            (A % P as u128) as u64,
            (i64::MIN as i128).rem_euclid(P as i128) as u64,
            P - 1,
        ];
        for (&pair, expected) in RESULTS.as_chunks::<2>().0.iter().zip(expected) {
            assert_eq!(pair, [expected, expected]);
        }
    }

    #[test]
    #[cfg(debug_assertions)]
    fn canonical_optimized_reducers_reject_invalid_inputs() {
        let unsupported = Modulus::new(17).unwrap();
        assert!(std::panic::catch_unwind(|| unsupported.reduce_opt(0)).is_err());
        assert!(std::panic::catch_unwind(|| unsafe { unsupported.reduce_opt_vt(0) }).is_err());
        assert!(std::panic::catch_unwind(|| unsupported.reduce_opt_u128(0)).is_err());
        assert!(std::panic::catch_unwind(|| unsafe { unsupported.reduce_opt_u128_vt(0) }).is_err());
        let modulus = Modulus::new((1 << 62) - 1).unwrap();
        let square = (*modulus as u128) * (*modulus as u128);
        assert!(std::panic::catch_unwind(|| modulus.reduce_opt_u128(square)).is_err());
        assert!(
            std::panic::catch_unwind(|| unsafe { modulus.reduce_opt_u128_vt(square) }).is_err()
        );
    }

    #[test]
    fn addition_vectors_cover_blocks_tails_and_slices() {
        for p in [
            2,
            3,
            4,
            17,
            65536,
            65537,
            1 << 61,
            (1 << 62) - 2,
            (1 << 62) - 1,
        ] {
            let modulus = Modulus::new(p).unwrap();
            // Include sums below, equal to, and above p, including 2*p - 2.
            let left = [0, 0, 1, p - 1, p - 1, p / 2, p - 2, p - 1];
            let right = [0, p - 1, p - 1, 0, p - 1, p - p / 2, 1, 1];
            for len in [0, 1, 2, 15, 16, 17, 31, 32, 33, 47, 48, 49, 63, 64, 65, 129] {
                let original: Vec<_> = left.into_iter().cycle().take(len).collect();
                let rhs: Vec<_> = right.into_iter().cycle().take(len).collect();
                let expected: Vec<_> = original
                    .iter()
                    .zip(&rhs)
                    .map(|(&a, &b)| ((a as u128 + b as u128) % p as u128) as u64)
                    .collect();
                for variable_time in [false, true] {
                    let mut buffer = vec![u64::MAX; len + 2];
                    buffer[1..=len].copy_from_slice(&original);
                    if variable_time {
                        unsafe {
                            modulus.add_vec_vt(&mut buffer[1..=len], &rhs);
                        }
                    } else {
                        modulus.add_vec(&mut buffer[1..=len], &rhs);
                    }
                    assert_eq!(&buffer[1..=len], expected);
                    assert_eq!((buffer[0], buffer[len + 1]), (u64::MAX, u64::MAX));
                }
            }
        }
    }

    #[test]
    #[cfg(debug_assertions)]
    fn addition_vectors_reject_invalid_inputs() {
        let modulus = Modulus::new(17).unwrap();
        for invalid in [17, 18, u64::MAX] {
            for index in [0, 15, 16] {
                let mut invalid_values = [0; 17];
                invalid_values[index] = invalid;
                for variable_time in [false, true] {
                    assert!(
                        std::panic::catch_unwind(|| {
                            let mut a = invalid_values;
                            if variable_time {
                                unsafe {
                                    modulus.add_vec_vt(&mut a, &[0; 17]);
                                }
                            } else {
                                modulus.add_vec(&mut a, &[0; 17]);
                            }
                        })
                        .is_err()
                    );
                    assert!(
                        std::panic::catch_unwind(|| {
                            let mut a = [0; 17];
                            if variable_time {
                                unsafe {
                                    modulus.add_vec_vt(&mut a, &invalid_values);
                                }
                            } else {
                                modulus.add_vec(&mut a, &invalid_values);
                            }
                        })
                        .is_err()
                    );
                }
            }
        }
        for (left, right) in [(0, 1), (1, 0), (16, 17), (17, 16)] {
            assert!(
                std::panic::catch_unwind(|| {
                    modulus.add_vec(&mut vec![0; left], &vec![0; right]);
                })
                .is_err()
            );
            assert!(
                std::panic::catch_unwind(|| unsafe {
                    modulus.add_vec_vt(&mut vec![0; left], &vec![0; right]);
                })
                .is_err()
            );
        }
    }

    #[test]
    fn subtraction_negation_exhaust_small_moduli() {
        for p in 2..=64 {
            let modulus = Modulus::new(p).unwrap();
            for a in 0..p {
                let expected = (-(a as i128)).rem_euclid(p as i128) as u64;
                assert_eq!(modulus.neg(a), expected);
                assert_eq!(unsafe { modulus.neg_vt(a) }, expected);
                assert_eq!(modulus.neg(modulus.neg(a)), a);
                for b in 0..p {
                    let expected = (a as i128 - b as i128).rem_euclid(p as i128) as u64;
                    assert_eq!(modulus.sub(a, b), expected);
                    assert_eq!(unsafe { super::verified::sub_mod_vt(a, b, p) }, expected);
                }
            }
        }
    }

    #[test]
    fn subtraction_negation_boundaries_and_const_kernels() {
        use super::verified;
        const P: u64 = (1 << 62) - 1;
        const RESULTS: [u64; 8] = [
            verified::sub_mod(0, P - 1, P),
            unsafe { verified::sub_mod_vt(0, P - 1, P) },
            verified::sub_mod(P - 1, P - 1, P),
            unsafe { verified::sub_mod_vt(P - 1, P - 1, P) },
            verified::neg_mod(0, P),
            unsafe { verified::neg_mod_vt(0, P) },
            verified::neg_mod(P - 1, P),
            unsafe { verified::neg_mod_vt(P - 1, P) },
        ];
        assert_eq!(RESULTS, [1, 1, 0, 0, 0, 0, 1, 1]);

        for p in reduction_boundary_moduli() {
            let modulus = Modulus::new(p).unwrap();
            let values = [0, 1, p / 2, p - 2, p - 1];
            for a in values {
                let expected = (-(a as i128)).rem_euclid(p as i128) as u64;
                assert_eq!(modulus.neg(a), expected);
                assert_eq!(unsafe { modulus.neg_vt(a) }, expected);
                for b in values {
                    let expected = (a as i128 - b as i128).rem_euclid(p as i128) as u64;
                    assert_eq!(modulus.sub(a, b), expected);
                    assert_eq!(unsafe { verified::sub_mod_vt(a, b, p) }, expected);
                }
            }
        }
    }

    #[test]
    fn subtraction_negation_vectors_cover_blocks_tails_and_slices() {
        for p in [
            2,
            3,
            4,
            17,
            65536,
            65537,
            1 << 61,
            (1 << 62) - 2,
            (1 << 62) - 1,
        ] {
            let modulus = Modulus::new(p).unwrap();
            let values = [0, 1, p - 1, p / 2, p - 2];
            for len in [0, 1, 2, 15, 16, 17, 31, 32, 33, 47, 48, 49, 63, 64, 65, 129] {
                let original: Vec<_> = values.into_iter().cycle().take(len).collect();
                let rhs: Vec<_> = values.into_iter().rev().cycle().take(len).collect();
                let differences: Vec<_> = original
                    .iter()
                    .zip(&rhs)
                    .map(|(&a, &b)| (a as i128 - b as i128).rem_euclid(p as i128) as u64)
                    .collect();
                let negatives: Vec<_> = original
                    .iter()
                    .map(|&a| (-(a as i128)).rem_euclid(p as i128) as u64)
                    .collect();

                for variable_time in [false, true] {
                    let mut buffer = vec![u64::MAX; len + 2];
                    buffer[1..=len].copy_from_slice(&original);
                    if variable_time {
                        unsafe {
                            modulus.sub_vec_vt(&mut buffer[1..=len], &rhs);
                        }
                    } else {
                        modulus.sub_vec(&mut buffer[1..=len], &rhs);
                    }
                    assert_eq!(&buffer[1..=len], differences);
                    assert_eq!((buffer[0], buffer[len + 1]), (u64::MAX, u64::MAX));

                    buffer[1..=len].copy_from_slice(&original);
                    if variable_time {
                        unsafe {
                            modulus.neg_vec_vt(&mut buffer[1..=len]);
                        }
                    } else {
                        modulus.neg_vec(&mut buffer[1..=len]);
                    }
                    assert_eq!(&buffer[1..=len], negatives);
                    assert_eq!((buffer[0], buffer[len + 1]), (u64::MAX, u64::MAX));
                }
            }
        }
    }

    #[test]
    #[cfg(debug_assertions)]
    fn subtraction_negation_reject_invalid_inputs() {
        let modulus = Modulus::new(17).unwrap();
        for invalid in [17, 18, u64::MAX] {
            assert!(std::panic::catch_unwind(|| modulus.sub(invalid, 0)).is_err());
            assert!(std::panic::catch_unwind(|| modulus.sub(0, invalid)).is_err());
            assert!(std::panic::catch_unwind(|| modulus.neg(invalid)).is_err());
            assert!(std::panic::catch_unwind(|| unsafe { modulus.neg_vt(invalid) }).is_err());
            for index in [0, 15, 16] {
                let mut invalid_values = [0; 17];
                invalid_values[index] = invalid;
                for variable_time in [false, true] {
                    assert!(
                        std::panic::catch_unwind(|| {
                            let mut a = invalid_values;
                            if variable_time {
                                unsafe {
                                    modulus.sub_vec_vt(&mut a, &[0; 17]);
                                }
                            } else {
                                modulus.sub_vec(&mut a, &[0; 17]);
                            }
                        })
                        .is_err()
                    );
                    assert!(
                        std::panic::catch_unwind(|| {
                            let mut a = [0; 17];
                            if variable_time {
                                unsafe {
                                    modulus.sub_vec_vt(&mut a, &invalid_values);
                                }
                            } else {
                                modulus.sub_vec(&mut a, &invalid_values);
                            }
                        })
                        .is_err()
                    );
                    assert!(
                        std::panic::catch_unwind(|| {
                            let mut a = invalid_values;
                            if variable_time {
                                unsafe {
                                    modulus.neg_vec_vt(&mut a);
                                }
                            } else {
                                modulus.neg_vec(&mut a);
                            }
                        })
                        .is_err()
                    );
                }
            }
        }
        for (left, right) in [(0, 1), (1, 0), (16, 17), (17, 16)] {
            assert!(
                std::panic::catch_unwind(|| {
                    modulus.sub_vec(&mut vec![0; left], &vec![0; right]);
                })
                .is_err()
            );
            assert!(
                std::panic::catch_unwind(|| unsafe {
                    modulus.sub_vec_vt(&mut vec![0; left], &vec![0; right]);
                })
                .is_err()
            );
        }
    }

    fn check_multiplication(modulus: &Modulus, a: u64, b: u64) {
        let p = **modulus;
        let expected = ((a as u128 * b as u128) % p as u128) as u64;
        let b_shoup = (((b as u128) << 64) / p as u128) as u64;
        assert_eq!(modulus.shoup(b), b_shoup);
        assert_eq!(modulus.mul_shoup(a, b, b_shoup), expected);
        assert_eq!(unsafe { modulus.mul_shoup_vt(a, b, b_shoup) }, expected);
        let lazy = modulus.lazy_mul_shoup(a, b, b_shoup);
        assert!(lazy < 2 * p);
        assert_eq!(lazy % p, expected);
        if a < p {
            assert_eq!(modulus.mul(a, b), expected);
            assert_eq!(unsafe { modulus.mul_vt(a, b) }, expected);
            if modulus.supports_opt {
                assert_eq!(modulus.mul_opt(a, b), expected);
                assert_eq!(unsafe { modulus.mul_opt_vt(a, b) }, expected);
            }
        }
    }

    #[test]
    fn multiplication_exhausts_small_moduli() {
        for p in 2..=64 {
            let modulus = Modulus::new(p).unwrap();
            for a in 0..2 * p {
                for b in 0..p {
                    check_multiplication(&modulus, a, b);
                }
            }
        }
    }

    #[test]
    fn multiplication_boundaries_and_const_kernels() {
        use super::verified;
        const P: u64 = (1 << 62) - 1;
        const MU: u128 = u128::MAX / P as u128;
        const HI: u64 = (MU >> 64) as u64;
        const LO: u64 = MU as u64;
        const B: u64 = P - 1;
        const B_SHOUP: u64 = verified::shoup(B, P);
        const RESULTS: [u64; 6] = [
            verified::mul_mod(B, B, P, HI, LO),
            unsafe { verified::mul_mod_vt(B, B, P, HI, LO) },
            verified::mul_opt(B, B, P, HI, LO, P.leading_zeros()),
            unsafe { verified::mul_opt_vt(B, B, P, HI, LO, P.leading_zeros()) },
            verified::mul_shoup(B, B, B_SHOUP, P),
            unsafe { verified::mul_shoup_vt(B, B, B_SHOUP, P) },
        ];
        const LAZY: u64 = verified::lazy_mul_shoup(u64::MAX, B, B_SHOUP, P);
        assert_eq!(RESULTS, [1; 6]);
        assert_eq!(B_SHOUP as u128, ((B as u128) << 64) / P as u128);
        let lazy = LAZY;
        assert!(lazy < 2 * P);
        assert_eq!(
            LAZY as u128 % P as u128,
            (u64::MAX as u128 * B as u128) % P as u128
        );

        for p in reduction_boundary_moduli() {
            let modulus = Modulus::new(p).unwrap();
            for a in [0, 1, p / 2, p - 1, p, 2 * p - 1, 1 << 63, u64::MAX] {
                for b in [0, 1, p / 2, p - 2, p - 1] {
                    check_multiplication(&modulus, a, b);
                }
            }
        }
    }

    #[test]
    fn multiplication_vectors_cover_boundaries_and_slices() {
        for p in reduction_boundary_moduli() {
            let modulus = Modulus::new(p).unwrap();
            let words = [0, 1, p - 1, p, 2 * p - 1, 1 << 63, u64::MAX];
            let rhs_words = [p - 1, p - 1, p - 1, 1, 0, p / 2, p - 2];
            for len in [0, 1, 2, 7, 8, 9, 15, 16, 17, 31, 32, 33, 65, 129] {
                let wide: Vec<_> = words.into_iter().cycle().take(len).collect();
                let canonical: Vec<_> = wide.iter().map(|a| a % p).collect();
                let rhs: Vec<_> = rhs_words.into_iter().cycle().take(len).collect();
                let shoup: Vec<_> = rhs
                    .iter()
                    .map(|&b| (((b as u128) << 64) / p as u128) as u64)
                    .collect();
                assert_eq!(modulus.shoup_vec(&rhs), shoup);
                let mut shoup_buffer = vec![u64::MAX; len + 2];
                modulus.shoup_into(&rhs, &mut shoup_buffer[1..=len]);
                assert_eq!(&shoup_buffer[1..=len], shoup);
                assert_eq!(
                    (shoup_buffer[0], shoup_buffer[len + 1]),
                    (u64::MAX, u64::MAX)
                );

                for variable_time in [false, true] {
                    for use_shoup in [false, true] {
                        let original = if use_shoup { &wide } else { &canonical };
                        let expected: Vec<_> = original
                            .iter()
                            .zip(&rhs)
                            .map(|(&a, &b)| ((a as u128 * b as u128) % p as u128) as u64)
                            .collect();
                        let mut buffer = vec![u64::MAX; len + 2];
                        buffer[1..=len].copy_from_slice(original);
                        match (use_shoup, variable_time) {
                            (false, false) => modulus.mul_vec(&mut buffer[1..=len], &rhs),
                            (false, true) => unsafe {
                                modulus.mul_vec_vt(&mut buffer[1..=len], &rhs)
                            },
                            (true, false) => {
                                modulus.mul_shoup_vec(&mut buffer[1..=len], &rhs, &shoup)
                            }
                            (true, true) => unsafe {
                                modulus.mul_shoup_vec_vt(&mut buffer[1..=len], &rhs, &shoup)
                            },
                        }
                        assert_eq!(&buffer[1..=len], expected);
                        assert_eq!((buffer[0], buffer[len + 1]), (u64::MAX, u64::MAX));
                    }
                    for b in [0, 1, p - 1] {
                        let expected: Vec<_> = wide
                            .iter()
                            .map(|&a| ((a as u128 * b as u128) % p as u128) as u64)
                            .collect();
                        let mut buffer = vec![u64::MAX; len + 2];
                        buffer[1..=len].copy_from_slice(&wide);
                        if variable_time {
                            unsafe {
                                modulus.scalar_mul_vec_vt(&mut buffer[1..=len], b);
                            }
                        } else {
                            modulus.scalar_mul_vec(&mut buffer[1..=len], b);
                        }
                        assert_eq!(&buffer[1..=len], expected);
                        assert_eq!((buffer[0], buffer[len + 1]), (u64::MAX, u64::MAX));
                    }
                }
            }
        }
    }

    #[test]
    #[cfg(debug_assertions)]
    fn multiplication_rejects_invalid_inputs() {
        for p in [17, (1 << 62) - 1] {
            let modulus = Modulus::new(p).unwrap();
            for invalid in [p, p + 1, u64::MAX] {
                assert!(
                    std::panic::catch_unwind(|| unsafe { modulus.mul_vt(invalid, 0) }).is_err()
                );
                assert!(
                    std::panic::catch_unwind(|| unsafe { modulus.mul_vt(0, invalid) }).is_err()
                );
                if modulus.supports_opt {
                    assert!(
                        std::panic::catch_unwind(|| unsafe { modulus.mul_opt_vt(invalid, 0) })
                            .is_err()
                    );
                    assert!(
                        std::panic::catch_unwind(|| unsafe { modulus.mul_opt_vt(0, invalid) })
                            .is_err()
                    );
                }
                assert!(
                    std::panic::catch_unwind(|| modulus.lazy_mul_shoup(0, invalid, 0)).is_err()
                );
                assert!(
                    std::panic::catch_unwind(|| unsafe { modulus.mul_shoup_vt(0, invalid, 0) })
                        .is_err()
                );
                assert!(
                    std::panic::catch_unwind(|| modulus.scalar_mul_vec(&mut [], invalid)).is_err()
                );
                assert!(
                    std::panic::catch_unwind(|| unsafe {
                        modulus.scalar_mul_vec_vt(&mut [], invalid)
                    })
                    .is_err()
                );
                for index in [0, 15, 16] {
                    let mut bad = [0; 17];
                    bad[index] = invalid;
                    assert!(
                        std::panic::catch_unwind(|| modulus.mul_vec(&mut { bad }, &[0; 17]))
                            .is_err()
                    );
                    assert!(
                        std::panic::catch_unwind(|| modulus.mul_vec(&mut [0; 17], &bad)).is_err()
                    );
                    assert!(
                        std::panic::catch_unwind(|| unsafe {
                            modulus.mul_vec_vt(&mut { bad }, &[0; 17])
                        })
                        .is_err()
                    );
                    assert!(
                        std::panic::catch_unwind(|| unsafe {
                            modulus.mul_vec_vt(&mut [0; 17], &bad)
                        })
                        .is_err()
                    );
                    assert!(std::panic::catch_unwind(|| modulus.shoup_vec(&bad)).is_err());
                    assert!(
                        std::panic::catch_unwind(|| modulus.shoup_into(&bad, &mut [0; 17]))
                            .is_err()
                    );
                    assert!(
                        std::panic::catch_unwind(|| modulus.mul_shoup_vec(
                            &mut [0; 17],
                            &bad,
                            &[0; 17]
                        ))
                        .is_err()
                    );
                    assert!(
                        std::panic::catch_unwind(|| unsafe {
                            modulus.mul_shoup_vec_vt(&mut [0; 17], &bad, &[0; 17])
                        })
                        .is_err()
                    );
                }
            }
            let b = p - 1;
            let wrong = modulus.shoup(b) ^ 1;
            assert!(std::panic::catch_unwind(|| modulus.lazy_mul_shoup(0, b, wrong)).is_err());
            assert!(std::panic::catch_unwind(|| modulus.mul_shoup(0, b, wrong)).is_err());
            assert!(
                std::panic::catch_unwind(|| unsafe { modulus.mul_shoup_vt(0, b, wrong) }).is_err()
            );
            assert!(
                std::panic::catch_unwind(|| modulus.mul_shoup_vec(&mut [0], &[b], &[wrong]))
                    .is_err()
            );
            assert!(
                std::panic::catch_unwind(|| unsafe {
                    modulus.mul_shoup_vec_vt(&mut [0], &[b], &[wrong])
                })
                .is_err()
            );
            for (left, right) in [(0, 1), (1, 0), (16, 17), (17, 16)] {
                assert!(
                    std::panic::catch_unwind(
                        || modulus.mul_vec(&mut vec![0; left], &vec![0; right])
                    )
                    .is_err()
                );
                assert!(
                    std::panic::catch_unwind(|| unsafe {
                        modulus.mul_vec_vt(&mut vec![0; left], &vec![0; right])
                    })
                    .is_err()
                );
                assert!(
                    std::panic::catch_unwind(
                        || modulus.shoup_into(&vec![0; left], &mut vec![0; right])
                    )
                    .is_err()
                );
                assert!(
                    std::panic::catch_unwind(|| modulus.mul_shoup_vec(
                        &mut vec![0; left],
                        &vec![0; right],
                        &vec![0; right]
                    ))
                    .is_err()
                );
                assert!(
                    std::panic::catch_unwind(|| unsafe {
                        modulus.mul_shoup_vec_vt(
                            &mut vec![0; left],
                            &vec![0; right],
                            &vec![0; right],
                        )
                    })
                    .is_err()
                );
                assert!(
                    std::panic::catch_unwind(|| modulus.mul_shoup_vec(
                        &mut vec![0; left],
                        &vec![0; left],
                        &vec![0; right]
                    ))
                    .is_err()
                );
                assert!(
                    std::panic::catch_unwind(|| unsafe {
                        modulus.mul_shoup_vec_vt(
                            &mut vec![0; left],
                            &vec![0; left],
                            &vec![0; right],
                        )
                    })
                    .is_err()
                );
            }
        }
        let unsupported = Modulus::new(17).unwrap();
        assert!(std::panic::catch_unwind(|| unsupported.mul_opt(0, 0)).is_err());
        assert!(std::panic::catch_unwind(|| unsafe { unsupported.mul_opt_vt(0, 0) }).is_err());
    }

    #[test]
    fn random_sampling_preserves_uniform_stream_and_slice_boundaries() {
        let mut saw_rejection = false;
        for p in reduction_boundary_moduli() {
            let modulus = Modulus::new(p).unwrap();
            for size in [0, 1, 2, 15, 16, 17, 63, 64, 65, 257] {
                for seed in [0, 0x5a17] {
                    let mut reference_rng = ChaCha8Rng::seed_from_u64(seed);
                    // Start partway through a block, including an odd word position.
                    reference_rng.set_word_pos(13);
                    let mut vector_rng = reference_rng.clone();
                    let mut slice_rng = reference_rng.clone();
                    // Preserve the original iterator implementation as the oracle.
                    let expected: Vec<u64> = (&mut reference_rng)
                        .sample_iter(Uniform::new(0, p).unwrap())
                        .take(size)
                        .collect();
                    let actual = modulus.random_vec(size, &mut vector_rng);
                    assert_eq!(actual.len(), size);
                    assert!(actual.iter().all(|&value| value < p));
                    assert_eq!(actual, expected, "p={p}, size={size}, seed={seed}");

                    let mut guarded = vec![u64::MAX; size + 4];
                    modulus.random_into(&mut guarded[2..size + 2], &mut slice_rng);
                    assert_eq!(&guarded[2..size + 2], expected);
                    assert_eq!(&guarded[..2], &[u64::MAX; 2]);
                    assert_eq!(&guarded[size + 2..], &[u64::MAX; 2]);
                    assert_eq!(vector_rng.get_word_pos(), reference_rng.get_word_pos());
                    assert_eq!(slice_rng.get_word_pos(), reference_rng.get_word_pos());
                    if size == 0 {
                        assert_eq!(reference_rng.get_word_pos(), 13);
                    }
                    // A u64 draw consumes two ChaCha words; extra words show
                    // that these compatibility cases also exercise rejection.
                    saw_rejection |= reference_rng.get_word_pos() > 13 + 2 * size as u128;
                }
            }
        }
        assert!(saw_rejection);
    }

    #[test]
    fn random_sampling_chunks_match_one_call() {
        let chunks = [0, 1, 0, 15, 16, 17, 255, 0, 256, 0];
        let size = chunks.iter().sum();
        for p in [2, 3, 16, 1153, (1 << 61) + 1, (1 << 62) - 1] {
            let modulus = Modulus::new(p).unwrap();
            let mut whole_rng = ChaCha8Rng::seed_from_u64(42);
            let mut chunk_rng = whole_rng.clone();
            let expected = modulus.random_vec(size, &mut whole_rng);
            let mut actual = vec![u64::MAX; size];
            let mut offset = 0;
            for (i, size) in chunks.into_iter().enumerate() {
                let chunk = &mut actual[offset..offset + size];
                if i % 2 == 0 {
                    chunk.copy_from_slice(&modulus.random_vec(size, &mut chunk_rng));
                } else {
                    modulus.random_into(chunk, &mut chunk_rng);
                }
                offset += size;
            }
            assert_eq!(actual, expected, "p={p}");
            assert_eq!(chunk_rng.get_word_pos(), whole_rng.get_word_pos());
        }
    }

    // Utility functions for the proptests.

    fn valid_moduli() -> impl Strategy<Value = Modulus> {
        any::<u64>().prop_filter_map("filter invalid moduli", |p| Modulus::new(p).ok())
    }

    fn valid_moduli_opt() -> impl Strategy<Value = Modulus> {
        valid_moduli().prop_filter("filter moduli not supporting opt", |p| p.supports_opt)
    }

    fn vecs() -> BoxedStrategy<(Vec<u64>, Vec<u64>)> {
        prop_vec(any::<u64>(), 1..100)
            .prop_flat_map(|vec| {
                let len = vec.len();
                (Just(vec), prop_vec(any::<u64>(), len))
            })
            .boxed()
    }

    proptest! {
        #[test]
        fn constructor(p: u64) {
            if (2..1 << 62).contains(&p) {
                check_constructor_metadata(p);
            } else {
                prop_assert_eq!(Modulus::new(p), Err(crate::Error::InvalidModulus(p)));
            }
            // Also exercise a valid 62-bit input on every iteration.
            check_constructor_metadata((p >> 2).max(2));
        }

        #[test]
        fn neg(p in valid_moduli(), mut a: u64) {
            a = p.reduce(a);
            let expected = (-(a as i128)).rem_euclid(*p as i128) as u64;
            prop_assert_eq!(p.neg(a), expected);
            unsafe { prop_assert_eq!(p.neg_vt(a), expected) }

            #[cfg(debug_assertions)]
            {
                prop_assert!(std::panic::catch_unwind(|| p.neg(*p)).is_err());
                prop_assert!(std::panic::catch_unwind(|| p.neg(*p + 1)).is_err());
            }
        }

        #[test]
        fn add(p in valid_moduli(), mut a: u64, mut b: u64) {
            a = p.reduce(a);
            b = p.reduce(b);
            prop_assert_eq!(p.add(a, b), (a + b) % *p);
            unsafe { prop_assert_eq!(p.add_vt(a, b), (a + b) % *p) }

            #[cfg(debug_assertions)]
            {
                prop_assert!(std::panic::catch_unwind(|| p.add(*p, a)).is_err());
                prop_assert!(std::panic::catch_unwind(|| p.add(a, *p)).is_err());
                prop_assert!(std::panic::catch_unwind(|| p.add(*p + 1, a)).is_err());
                prop_assert!(std::panic::catch_unwind(|| p.add(a, *p + 1)).is_err());
                unsafe {
                    prop_assert!(std::panic::catch_unwind(|| p.add_vt(*p, a)).is_err());
                    prop_assert!(std::panic::catch_unwind(|| p.add_vt(a, *p)).is_err());
                    prop_assert!(std::panic::catch_unwind(|| p.add_vt(*p + 1, a)).is_err());
                    prop_assert!(std::panic::catch_unwind(|| p.add_vt(a, *p + 1)).is_err());
                }
            }
        }

        #[test]
        fn sub(p in valid_moduli(), mut a: u64, mut b: u64) {
            a = p.reduce(a);
            b = p.reduce(b);
            let expected = (a as i128 - b as i128).rem_euclid(*p as i128) as u64;
            prop_assert_eq!(p.sub(a, b), expected);
            prop_assert_eq!(unsafe { super::verified::sub_mod_vt(a, b, *p) }, expected);

            #[cfg(debug_assertions)]
            {
                prop_assert!(std::panic::catch_unwind(|| p.sub(*p, a)).is_err());
                prop_assert!(std::panic::catch_unwind(|| p.sub(a, *p)).is_err());
                prop_assert!(std::panic::catch_unwind(|| p.sub(*p + 1, a)).is_err());
                prop_assert!(std::panic::catch_unwind(|| p.sub(a, *p + 1)).is_err());
            }
        }

        #[test]
        fn mul(p in valid_moduli(), mut a: u64, mut b: u64) {
            a = p.reduce(a);
            b = p.reduce(b);
            prop_assert_eq!(p.mul(a, b) as u128, ((a as u128) * (b as u128)) % (*p as u128));
            unsafe { prop_assert_eq!(p.mul_vt(a, b) as u128, ((a as u128) * (b as u128)) % (*p as u128)) }

            #[cfg(debug_assertions)]
            {
                prop_assert!(std::panic::catch_unwind(|| p.mul(*p, a)).is_err());
                prop_assert!(std::panic::catch_unwind(|| p.mul(a, *p)).is_err());
                prop_assert!(std::panic::catch_unwind(|| p.mul(*p + 1, a)).is_err());
                prop_assert!(std::panic::catch_unwind(|| p.mul(a, *p + 1)).is_err());
            }
        }

        #[test]
        fn mul_shoup(p in valid_moduli(), a: u64, mut b: u64) {
            b = p.reduce(b);

            // Compare the precomputation and lazy result to independent arithmetic.
            let b_shoup = p.shoup(b);
            prop_assert_eq!(b_shoup as u128, ((b as u128) << 64) / *p as u128);
            let lazy = p.lazy_mul_shoup(a, b, b_shoup);
            prop_assert!(lazy < 2 * *p);
            prop_assert_eq!(lazy as u128 % *p as u128, (a as u128 * b as u128) % *p as u128);

            #[cfg(debug_assertions)]
            {
                prop_assert!(std::panic::catch_unwind(|| p.shoup(*p)).is_err());
                prop_assert!(std::panic::catch_unwind(|| p.shoup(*p + 1)).is_err());
            }

            // Check that the multiplication yields the expected result
            prop_assert_eq!(p.mul_shoup(a, b, b_shoup) as u128, ((a as u128) * (b as u128)) % (*p as u128));
            unsafe { prop_assert_eq!(p.mul_shoup_vt(a, b, b_shoup) as u128, ((a as u128) * (b as u128)) % (*p as u128)) }

            // Check that the multiplication with incorrect b_shoup panics in debug mode
            #[cfg(debug_assertions)]
            {
                prop_assert!(std::panic::catch_unwind(|| p.mul_shoup(a, *p, b_shoup)).is_err());
                prop_assert!(std::panic::catch_unwind(|| p.mul_shoup(a, b, b_shoup ^ 1)).is_err());
            }
        }

        #[test]
        fn reduce(p in valid_moduli(), a: u64) {
            prop_assert_eq!(p.reduce(a), a % *p);
            unsafe { prop_assert_eq!(p.reduce_vt(a), a % *p) }
            if p.supports_opt {
                prop_assert_eq!(p.reduce_opt(a), a % *p);
                unsafe { prop_assert_eq!(p.reduce_opt_vt(a), a % *p) }
            }
        }

        #[test]
        fn lazy_reduce(p in valid_moduli(), a: u64) {
            prop_assert!(p.lazy_reduce(a) < 2 * *p);
            prop_assert_eq!(p.lazy_reduce(a) % *p, a % *p);
        }

        #[test]
        fn lazy_reduce_u128(p in valid_moduli(), a: u128) {
            let r = p.lazy_reduce_u128(a);
            prop_assert!(r < 2 * *p);
            prop_assert_eq!(r as u128 % *p as u128, a % *p as u128);
        }

        #[test]
        fn lazy_reduce_optimized(p in valid_moduli_opt(), a: u64, wide: u128) {
            let r = p.lazy_reduce_opt(a);
            prop_assert!(r < 2 * *p);
            prop_assert_eq!(r % *p, a % *p);
            let a = wide % ((*p as u128) * (*p as u128));
            let r = p.lazy_reduce_opt_u128(a);
            prop_assert!(r < 2 * *p);
            prop_assert_eq!(r as u128 % *p as u128, a % *p as u128);
        }

        #[test]
        fn reduce_i64(p in valid_moduli(), a: i64) {
            let b = (a as i128).rem_euclid(*p as i128) as u64;
            prop_assert_eq!(p.reduce_i64(a), b);
            unsafe { prop_assert_eq!(p.reduce_i64_vt(a), b) }
        }

        #[test]
        fn reduce_u128(p in valid_moduli(), mut a: u128) {
            prop_assert_eq!(p.reduce_u128(a) as u128, a % (*p as u128));
            unsafe { prop_assert_eq!(p.reduce_u128_vt(a) as u128, a % (*p as u128)) }
            if p.supports_opt {
                let p_square = (*p as u128) * (*p as u128);
                a %= p_square;
                prop_assert_eq!(p.reduce_opt_u128(a) as u128, a % (*p as u128));
                unsafe { prop_assert_eq!(p.reduce_opt_u128_vt(a) as u128, a % (*p as u128)) }
            }
        }

        #[test]
        fn add_vec(p in valid_moduli(), (mut a, mut b) in vecs()) {
            p.reduce_vec(&mut a);
            p.reduce_vec(&mut b);
            let c = a.clone();
            let expected: Vec<_> = c.iter().zip(&b)
                .map(|(&a, &b)| ((a as u128 + b as u128) % *p as u128) as u64)
                .collect();
            p.add_vec(&mut a, &b);
            prop_assert_eq!(&a, &expected);
            a.clone_from(&c);
            unsafe { p.add_vec_vt(&mut a, &b) }
            prop_assert_eq!(a, expected);
        }

        #[test]
        fn sub_vec(p in valid_moduli(), (mut a, mut b) in vecs()) {
            p.reduce_vec(&mut a);
            p.reduce_vec(&mut b);
            let c = a.clone();
            p.sub_vec(&mut a, &b);
            prop_assert_eq!(a.clone(), izip!(b.iter(), c.iter()).map(|(bi, ci)| (*ci as i128 - *bi as i128).rem_euclid(*p as i128) as u64).collect_vec());
            a.clone_from(&c);
            unsafe { p.sub_vec_vt(&mut a, &b) }
            prop_assert_eq!(a, izip!(b.iter(), c.iter()).map(|(bi, ci)| (*ci as i128 - *bi as i128).rem_euclid(*p as i128) as u64).collect_vec());
        }

        #[test]
        fn mul_vec(p in valid_moduli(), (mut a, mut b) in vecs()) {
            p.reduce_vec(&mut a);
            p.reduce_vec(&mut b);
            let c = a.clone();
            p.mul_vec(&mut a, &b);
            prop_assert_eq!(a.clone(), izip!(b.iter(), c.iter()).map(|(bi, ci)| ((*ci as u128 * *bi as u128) % *p as u128) as u64).collect_vec());
            a.clone_from(&c);
            unsafe { p.mul_vec_vt(&mut a, &b); }
            prop_assert_eq!(a, izip!(b.iter(), c.iter()).map(|(bi, ci)| ((*ci as u128 * *bi as u128) % *p as u128) as u64).collect_vec());
        }

        #[test]
        fn scalar_mul_vec(p in valid_moduli(), mut a: Vec<u64>, mut b: u64) {
            b = p.reduce(b);
            let c = a.clone();

            p.scalar_mul_vec(&mut a, b);
            prop_assert_eq!(a.clone(), c.iter().map(|ci| ((*ci as u128 * b as u128) % *p as u128) as u64).collect_vec());

            a.clone_from(&c);
            unsafe { p.scalar_mul_vec_vt(&mut a, b) }
            prop_assert_eq!(a, c.iter().map(|ci| ((*ci as u128 * b as u128) % *p as u128) as u64).collect_vec());
        }

        #[test]
        fn mul_shoup_vec(p in valid_moduli(), (mut a, mut b) in vecs()) {
            p.reduce_vec(&mut b);
            let b_shoup = p.shoup_vec(&b);
            let c = a.clone();
            p.mul_shoup_vec(&mut a, &b, &b_shoup);
            prop_assert_eq!(a.clone(), izip!(b.iter(), c.iter()).map(|(bi, ci)| ((*ci as u128 * *bi as u128) % *p as u128) as u64).collect_vec());
            a.clone_from(&c);
            unsafe { p.mul_shoup_vec_vt(&mut a, &b, &b_shoup) }
            prop_assert_eq!(a, izip!(b.iter(), c.iter()).map(|(bi, ci)| ((*ci as u128 * *bi as u128) % *p as u128) as u64).collect_vec());
        }

        #[test]
        fn reduce_vec(p in valid_moduli(), a: Vec<u64>) {
            let mut b = a.clone();
            p.reduce_vec(&mut b);
            prop_assert_eq!(b.clone(), a.iter().map(|ai| *ai % *p).collect_vec());

            b.clone_from(&a);
            unsafe { p.reduce_vec_vt(&mut b) }
            prop_assert_eq!(b, a.iter().map(|ai| *ai % *p).collect_vec());
        }

        #[test]
        fn lazy_reduce_vec(p in valid_moduli(), a: Vec<u64>) {
            let mut b = a.clone();
            p.lazy_reduce_vec(&mut b);
            prop_assert!(b.iter().all(|bi| *bi < 2 * *p));
            prop_assert!(izip!(a, b).all(|(ai, bi)| bi % *p == ai % *p));
        }

        #[test]
        fn reduce_vec_new(p in valid_moduli(), a: Vec<u64>) {
            let b = p.reduce_vec_new(&a);
            prop_assert_eq!(b, a.iter().map(|ai| *ai % *p).collect_vec());
            prop_assert_eq!(p.reduce_vec_new(&a), unsafe { p.reduce_vec_new_vt(&a) });
        }

        #[test]
        fn reduce_vec_i64(p in valid_moduli(), a: Vec<i64>) {
            let b = p.reduce_vec_i64(&a);
            prop_assert_eq!(b, a.iter().map(|ai| (*ai as i128).rem_euclid(*p as i128) as u64).collect_vec());
            let b = unsafe { p.reduce_vec_i64_vt(&a) };
            prop_assert_eq!(b, a.iter().map(|ai| (*ai as i128).rem_euclid(*p as i128) as u64).collect_vec());
        }

        #[test]
        fn neg_vec(p in valid_moduli(), mut a: Vec<u64>) {
            p.reduce_vec(&mut a);
            let mut b = a.clone();
            p.neg_vec(&mut b);
            prop_assert_eq!(b.clone(), a.iter().map(|ai| (-(*ai as i128)).rem_euclid(*p as i128) as u64).collect_vec());
            b.clone_from(&a);
            unsafe { p.neg_vec_vt(&mut b); }
            prop_assert_eq!(b, a.iter().map(|ai| (-(*ai as i128)).rem_euclid(*p as i128) as u64).collect_vec());
        }

        #[test]
        fn random_vec(p in valid_moduli(), size in 0..1000usize, seed: u64) {
            let mut actual_rng = ChaCha8Rng::seed_from_u64(seed);
            let mut reference_rng = actual_rng.clone();
            let expected: Vec<u64> = (&mut reference_rng)
                .sample_iter(Uniform::new(0, *p).unwrap()).take(size).collect();
            let actual = p.random_vec(size, &mut actual_rng);
            prop_assert_eq!(actual.len(), size);
            prop_assert!(actual.iter().all(|&value| value < *p));
            prop_assert_eq!(actual, expected);
            prop_assert_eq!(actual_rng.get_word_pos(), reference_rng.get_word_pos());
        }

        #[test]
        fn random_into(p in valid_moduli(), mut out: Vec<u64>, seed: u64) {
            let mut actual_rng = ChaCha8Rng::seed_from_u64(seed);
            let mut reference_rng = actual_rng.clone();
            let expected: Vec<u64> = (&mut reference_rng)
                .sample_iter(Uniform::new(0, *p).unwrap()).take(out.len()).collect();
            p.random_into(&mut out, &mut actual_rng);
            prop_assert!(out.iter().all(|&value| value < *p));
            prop_assert_eq!(out, expected);
            prop_assert_eq!(actual_rng.get_word_pos(), reference_rng.get_word_pos());
        }

        #[test]
        fn serialize(p in valid_moduli(), mut a in prop_vec(any::<u64>(), 8)) {
            p.reduce_vec(&mut a);
            let b = p.serialize_vec(&a).unwrap();
            let c = p.deserialize_vec(&b).unwrap();
            prop_assert_eq!(a, c);
        }

        #[test]
        fn center(p in valid_moduli(), a: u64) {
            let a = a % *p;
            let b = p.center(a);
            let half = (*p / 2) as i128;
            prop_assert_eq!(b as i128, (a as i128 + half).rem_euclid(*p as i128) - half);
            prop_assert!(-(*p as i128) <= 2 * b as i128 && 2 * (b as i128) < *p as i128);
            prop_assert_eq!((b as i128).rem_euclid(*p as i128), a as i128);
        }

        #[test]
        fn center_vec(p in valid_moduli(), mut a: Vec<u64>) {
            a.iter_mut().for_each(|value| *value %= *p);
            let half = (*p / 2) as i128;
            let expected: Vec<_> = a.iter()
                .map(|&value| ((value as i128 + half).rem_euclid(*p as i128) - half) as i64)
                .collect();
            prop_assert_eq!(p.center_vec(&a), expected);
        }

        #[test]
        fn mul_opt(p in valid_moduli_opt(), mut a: u64, mut b: u64) {
            a = p.reduce(a);
            b = p.reduce(b);

            prop_assert_eq!(p.mul_opt(a, b) as u128, ((a as u128) * (b as u128)) % (*p as u128));
            unsafe { prop_assert_eq!(p.mul_opt_vt(a, b) as u128, ((a as u128) * (b as u128)) % (*p as u128)) }

            #[cfg(debug_assertions)]
            {
                prop_assert!(std::panic::catch_unwind(|| p.mul_opt(*p, a)).is_err());
                prop_assert!(std::panic::catch_unwind(|| p.mul_opt(a, *p)).is_err());
                prop_assert!(std::panic::catch_unwind(|| p.mul_opt(*p + 1, a)).is_err());
                prop_assert!(std::panic::catch_unwind(|| p.mul_opt(a, *p + 1)).is_err());
            }
        }

        #[test]
        fn pow(p in valid_moduli(), mut a: u64, mut b: u64) {
            use num_bigint::BigUint;
            a %= *p;
            b %= *p;

            prop_assert_eq!(p.pow(a, 0), 1);
            prop_assert_eq!(p.pow(a, 1), a);
            if *p > 2 {
                prop_assert_eq!(p.pow(a, 2) as u128, (a as u128 * a as u128) % *p as u128);
            }

            let result = p.pow(a, b);
            prop_assert!(result < *p);
            prop_assert_eq!(BigUint::from(result),
                BigUint::from(a).modpow(&BigUint::from(b), &BigUint::from(*p)));

            #[cfg(debug_assertions)]
            {
                prop_assert!(std::panic::catch_unwind(|| p.pow(*p, 1)).is_err());
                prop_assert!(std::panic::catch_unwind(|| p.pow(*p << 1, 1)).is_err());
                prop_assert!(std::panic::catch_unwind(|| p.pow(0, *p)).is_err());
                prop_assert!(std::panic::catch_unwind(|| p.pow(0, *p << 1)).is_err());
            }
        }

        #[test]
        fn inv_matches_extended_euclid(
            p in proptest::sample::select(INVERSION_PRIMES), a: u64,
        ) {
            let modulus = Modulus::new(p).unwrap();
            let a = a % p;
            let actual = modulus.inv(a);
            prop_assert_eq!(actual, fhe_util::inverse(a, p));
            if let Some(r) = actual {
                prop_assert!(0 < r && r < p);
                prop_assert_eq!((a as u128 * r as u128) % p as u128, 1);
            }
        }
    }

    #[test]
    fn pow_exhausts_small_moduli() {
        for p in 2..=64 {
            let modulus = Modulus::new(p).unwrap();
            for a in 0..p {
                let mut expected = 1;
                for n in 0..p {
                    assert_eq!(modulus.pow(a, n), expected, "p={p}, a={a}, n={n}");
                    expected = ((expected as u128 * a as u128) % p as u128) as u64;
                }
            }
        }
    }

    #[test]
    fn pow_matches_bigint_at_bit_boundaries() {
        use num_bigint::BigUint;

        for p in reduction_boundary_moduli() {
            let modulus = Modulus::new(p).unwrap();
            let mut exponents = vec![
                0,
                1,
                2,
                p / 2 - 1,
                p / 2,
                p - 2,
                p - 1,
                0x1555_5555_5555_5555 % p,
                0x2aaa_aaaa_aaaa_aaaa % p,
            ];
            // At the widest modulus, cover each possible leading-bit position
            // and the transition on either side of every power of two.
            if p == (1 << 62) - 1 {
                for bit in 0..62 {
                    let power = 1u64 << bit;
                    exponents.extend([power - 1, power, power + 1]);
                }
            }
            exponents.retain(|&n| n < p);
            exponents.sort_unstable();
            exponents.dedup();
            let big_p = BigUint::from(p);
            for a in [0, 1, 2 % p, p / 2, p - 1] {
                let big_a = BigUint::from(a);
                for &n in &exponents {
                    let result = modulus.pow(a, n);
                    assert!(result < p);
                    assert_eq!(
                        BigUint::from(result),
                        big_a.modpow(&BigUint::from(n), &big_p),
                        "p={p}, a={a}, n={n}"
                    );
                }
            }
        }
    }

    #[test]
    #[cfg(debug_assertions)]
    fn pow_rejects_invalid_base_and_exponent() {
        for p in [2, 3, 16, 17, (1 << 62) - 1] {
            let modulus = Modulus::new(p).unwrap();
            for invalid in [p, p + 1, u64::MAX] {
                for valid in [0, 1, p - 1] {
                    assert!(std::panic::catch_unwind(|| modulus.pow(invalid, valid)).is_err());
                    assert!(std::panic::catch_unwind(|| modulus.pow(valid, invalid)).is_err());
                }
            }
        }
    }

    const INVERSION_PRIMES: &[u64] = &[
        2,
        3,
        17,
        1987,
        65537,
        1_000_000_007,
        2_147_483_647,
        2_305_843_009_213_693_951,
        4_611_686_018_326_724_609,
    ];

    #[test]
    fn inv_exhausts_small_moduli() {
        for p in 2..=256 {
            let modulus = Modulus::new(p).unwrap();
            // Independent trial division and exhaustive inverse search.
            let prime = (2..p).all(|divisor| p % divisor != 0);
            for a in 0..p {
                let expected = if prime {
                    (1..p).find(|&r| a * r % p == 1)
                } else {
                    None
                };
                assert_eq!(modulus.inv(a), expected, "p={p}, a={a}");
            }
        }
    }

    #[test]
    fn inv_matches_extended_euclid_at_boundaries() {
        for &p in INVERSION_PRIMES {
            assert!(fhe_util::is_prime(p));
            let modulus = Modulus::new(p).unwrap();
            let mut operands = vec![0, 1, p / 2, p - 2, p - 1];
            for bit in 1..62 {
                for a in [(1u64 << bit) - 1, 1u64 << bit, (1u64 << bit) + 1] {
                    if a < p {
                        operands.push(a);
                    }
                }
            }
            for a in operands {
                let actual = modulus.inv(a);
                assert_eq!(actual, fhe_util::inverse(a, p), "p={p}, a={a}");
                if let Some(r) = actual {
                    assert!(0 < r && r < p);
                    assert_eq!((a as u128 * r as u128) % p as u128, 1);
                }
            }
        }
    }

    #[test]
    fn inv_rejects_composites_including_invertible_operands() {
        for p in [
            4,
            6,
            9,
            15,
            561,
            1105,
            1729,
            1 << 61,
            (1 << 62) - 2,
            (1 << 62) - 1,
        ] {
            let modulus = Modulus::new(p).unwrap();
            assert_eq!(fhe_util::inverse(1, p), Some(1));
            for a in [0, 1, 2, p / 2, p - 1, p, p + 1, u64::MAX] {
                assert_eq!(modulus.inv(a), None, "p={p}, a={a}");
            }
        }
    }

    #[test]
    #[cfg(debug_assertions)]
    fn inv_rejects_noncanonical_operands_for_primes() {
        for &p in INVERSION_PRIMES {
            let modulus = Modulus::new(p).unwrap();
            for a in [p, p + 1, 2 * p, u64::MAX] {
                assert!(std::panic::catch_unwind(|| modulus.inv(a)).is_err());
            }
        }
    }

    #[test]
    fn centered_boundary_is_correct_for_odd_and_even_moduli() {
        for p in (2u64..=64).chain([1024, 1153]) {
            let modulus = Modulus::new(p).unwrap();
            let half = (p / 2) as i128;
            let input: Vec<_> = (0..p).collect();
            let expected: Vec<_> = input
                .iter()
                .map(|&value| ((value as i128 + half).rem_euclid(p as i128) - half) as i64)
                .collect();
            assert_eq!(modulus.center_vec(&input), expected);
            for value in 0..p {
                let centered = modulus.center(value);
                assert!(2 * centered >= -(p as i64));
                assert!(2 * centered < p as i64);
                assert_eq!((centered as i128).rem_euclid(p as i128), value as i128);
                assert_eq!(centered, expected[value as usize]);
            }
        }
        assert_eq!(Modulus::new(1153).unwrap().center(576), 576);
    }

    #[test]
    fn center_boundaries_vectors_and_const_kernel() {
        use super::verified;
        const P: u64 = (1 << 62) - 1;
        const RESULTS: [i64; 6] = [
            verified::center(0, P),
            verified::center(P / 2, P),
            verified::center(P / 2 + 1, P),
            verified::center(P - 1, P),
            verified::center((P - 1) / 2, P - 1),
            verified::center(1, 2),
        ];
        assert_eq!(
            RESULTS,
            [
                0,
                (P / 2) as i64,
                -((P / 2) as i64),
                -1,
                -(((P - 1) / 2) as i64),
                -1
            ]
        );

        for p in reduction_boundary_moduli() {
            let modulus = Modulus::new(p).unwrap();
            let half = (p / 2) as i128;
            let threshold = p.div_ceil(2);
            let values: Vec<_> = [0, 1, threshold - 1, threshold, threshold + 1, p - 1]
                .into_iter()
                .filter(|&a| a < p)
                .collect();
            for &a in &values {
                let expected = (a as i128 + half).rem_euclid(p as i128) - half;
                let centered = modulus.center(a) as i128;
                assert_eq!(centered, expected);
                assert!(-(p as i128) <= 2 * centered && 2 * centered < p as i128);
                assert_eq!(centered.rem_euclid(p as i128), a as i128);
            }
            for len in [0, 1, 2, 7, 8, 9, 15, 16, 17, 31, 32, 33, 65, 129] {
                let input: Vec<_> = values.iter().copied().cycle().take(len).collect();
                let expected: Vec<_> = input
                    .iter()
                    .map(|&a| ((a as i128 + half).rem_euclid(p as i128) - half) as i64)
                    .collect();
                let mut buffer = vec![u64::MAX; len + 2];
                buffer[1..=len].copy_from_slice(&input);
                assert_eq!(modulus.center_vec(&buffer[1..=len]), expected);
                assert_eq!(&buffer[1..=len], input);
                assert_eq!((buffer[0], buffer[len + 1]), (u64::MAX, u64::MAX));
            }
        }
    }

    #[test]
    #[cfg(debug_assertions)]
    fn center_rejects_noncanonical_inputs() {
        for p in [2, 17, (1 << 62) - 1] {
            let modulus = Modulus::new(p).unwrap();
            for invalid in [p, p + 1, u64::MAX] {
                assert!(std::panic::catch_unwind(|| modulus.center(invalid)).is_err());
                for index in [0, 7, 16] {
                    let mut input = [0; 17];
                    input[index] = invalid;
                    assert!(std::panic::catch_unwind(|| modulus.center_vec(&input)).is_err());
                }
            }
        }
    }
}
