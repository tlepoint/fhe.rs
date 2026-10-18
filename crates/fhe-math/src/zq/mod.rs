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
use fhe_util::{is_prime, transcode_from_bytes, transcode_to_bytes, transcode_to_bytes_into};
use itertools::{Itertools, izip};
use num_bigint::BigUint;
use num_traits::cast::ToPrimitive;
use pulp::Arch;
use rand::{CryptoRng, Rng as RngCore, RngExt, distr::Uniform};
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
        if p < 2 || (p >> 62) != 0 {
            Err(Error::InvalidModulus(p))
        } else {
            let barrett = ((BigUint::from(1u64) << 128usize) / p).to_u128().unwrap(); // 2^128 / p
            Ok(Self {
                p,
                barrett_hi: (barrett >> 64) as u64,
                barrett_lo: barrett as u64,
                leading_zeros: p.leading_zeros(),
                supports_opt: primes::supports_opt(p),
                distribution: Uniform::new(0, p).unwrap(),
                arch: Arch::new(),
            })
        }
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
        self.reduce_u128((a as u128) * (b as u128))
    }

    /// Performs the modular multiplication of a and b in constant time.
    /// Aborts if a >= p or b >= p in debug mode.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being multiplied.
    const unsafe fn mul_vt(&self, a: u64, b: u64) -> u64 {
        debug_assert!(a < self.p && b < self.p);
        unsafe { Self::reduce1_vt(self.lazy_reduce_u128((a as u128) * (b as u128)), self.p) }
    }

    /// Optimized modular multiplication of a and b in constant time.
    ///
    /// Aborts if a >= p or b >= p in debug mode.
    #[must_use]
    pub const fn mul_opt(&self, a: u64, b: u64) -> u64 {
        debug_assert!(self.supports_opt);
        debug_assert!(a < self.p && b < self.p);

        self.reduce_opt_u128((a as u128) * (b as u128))
    }

    /// Optimized modular multiplication of a and b in variable time.
    /// Aborts if a >= p or b >= p in debug mode.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being multiplied.
    const unsafe fn mul_opt_vt(&self, a: u64, b: u64) -> u64 {
        debug_assert!(self.supports_opt);
        debug_assert!(a < self.p && b < self.p);

        unsafe { self.reduce_opt_u128_vt((a as u128) * (b as u128)) }
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

        (((a as u128) << 64) / (self.p as u128)) as u64
    }

    /// Shoup multiplication of a and b in constant time.
    ///
    /// Aborts if b >= p or b_shoup != shoup(b) in debug mode.
    #[must_use]
    pub const fn mul_shoup(&self, a: u64, b: u64, b_shoup: u64) -> u64 {
        Self::reduce1(self.lazy_mul_shoup(a, b, b_shoup), self.p)
    }

    /// Shoup multiplication of a and b in variable time.
    /// Aborts if b >= p or b_shoup != shoup(b) in debug mode.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being multiplied.
    const unsafe fn mul_shoup_vt(&self, a: u64, b: u64, b_shoup: u64) -> u64 {
        unsafe { Self::reduce1_vt(self.lazy_mul_shoup(a, b, b_shoup), self.p) }
    }

    /// Lazy Shoup multiplication of a and b in constant time.
    /// The output is in the interval [0, 2 * p).
    ///
    /// Aborts if b >= p or b_shoup != shoup(b) in debug mode.
    #[must_use]
    pub const fn lazy_mul_shoup(&self, a: u64, b: u64, b_shoup: u64) -> u64 {
        debug_assert!(b < self.p);
        debug_assert!(b_shoup == self.shoup(b));

        let q = ((a as u128) * (b_shoup as u128)) >> 64;
        let r = ((a as u128) * (b as u128) - q * (self.p as u128)) as u64;

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

        if self.supports_opt {
            self.arch.dispatch(|| {
                izip!(a.iter_mut(), b.iter()).for_each(|(ai, bi)| *ai = self.mul_opt(*ai, *bi))
            })
        } else {
            self.arch.dispatch(|| {
                izip!(a.iter_mut(), b.iter()).for_each(|(ai, bi)| *ai = self.mul(*ai, *bi))
            })
        }
    }

    /// Modular scalar multiplication of vectors in place in constant time.
    ///
    /// Aborts if any of the values in a is >= p in debug mode.
    pub fn scalar_mul_vec(&self, a: &mut [u64], b: u64) {
        let b_shoup = self.shoup(b);
        self.arch.dispatch(|| {
            a.iter_mut()
                .for_each(|ai| *ai = self.mul_shoup(*ai, b, b_shoup))
        })
    }

    /// Modular scalar multiplication of vectors in place in variable time.
    /// Aborts if any of the values in a is >= p in debug mode.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being multiplied.
    pub unsafe fn scalar_mul_vec_vt(&self, a: &mut [u64], b: u64) {
        let b_shoup = self.shoup(b);
        self.arch.dispatch(|| {
            a.iter_mut()
                .for_each(|ai| *ai = unsafe { self.mul_shoup_vt(*ai, b, b_shoup) })
        })
    }

    /// Modular multiplication of vectors in place in variable time.
    /// Aborts if a and b differ in size, and if any of their values is >= p in
    /// debug mode.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being subtracted.
    pub unsafe fn mul_vec_vt(&self, a: &mut [u64], b: &[u64]) {
        debug_assert_eq!(a.len(), b.len());

        if self.supports_opt {
            self.arch.dispatch(|| {
                izip!(a.iter_mut(), b.iter())
                    .for_each(|(ai, bi)| *ai = unsafe { self.mul_opt_vt(*ai, *bi) })
            })
        } else {
            self.arch.dispatch(|| {
                izip!(a.iter_mut(), b.iter())
                    .for_each(|(ai, bi)| *ai = unsafe { self.mul_vt(*ai, *bi) })
            })
        }
    }

    /// Compute the Shoup representation of a vector.
    ///
    /// Aborts if any of the values of the vector is >= p in debug mode.
    #[must_use]
    pub fn shoup_vec(&self, a: &[u64]) -> Vec<u64> {
        self.arch
            .dispatch(|| a.iter().map(|ai| self.shoup(*ai)).collect_vec())
    }

    pub(crate) fn shoup_into(&self, a: &[u64], out: &mut [u64]) {
        assert_eq!(a.len(), out.len());
        self.arch.dispatch(|| {
            for (ai, output) in a.iter().zip(out.iter_mut()) {
                *output = self.shoup(*ai);
            }
        });
    }

    /// Shoup modular multiplication of vectors in place in constant time.
    ///
    /// Aborts if a and b differ in size, and if any of their values is >= p in
    /// debug mode.
    pub fn mul_shoup_vec(&self, a: &mut [u64], b: &[u64], b_shoup: &[u64]) {
        debug_assert_eq!(a.len(), b.len());
        debug_assert_eq!(a.len(), b_shoup.len());
        debug_assert_eq!(&b_shoup, &self.shoup_vec(b));

        self.arch.dispatch(|| {
            izip!(a.iter_mut(), b.iter(), b_shoup.iter())
                .for_each(|(ai, bi, bi_shoup)| *ai = self.mul_shoup(*ai, *bi, *bi_shoup))
        })
    }

    /// Shoup modular multiplication of vectors in place in variable time.
    /// Aborts if a and b differ in size, and if any of their values is >= p in
    /// debug mode.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being multiplied.
    pub unsafe fn mul_shoup_vec_vt(&self, a: &mut [u64], b: &[u64], b_shoup: &[u64]) {
        debug_assert_eq!(a.len(), b.len());
        debug_assert_eq!(a.len(), b_shoup.len());
        debug_assert_eq!(&b_shoup, &self.shoup_vec(b));

        self.arch.dispatch(|| {
            izip!(a.iter_mut(), b.iter(), b_shoup.iter()).for_each(|(ai, bi, bi_shoup)| {
                *ai = unsafe { self.mul_shoup_vt(*ai, *bi, *bi_shoup) }
            })
        })
    }

    /// Reduce a vector in place in constant time.
    pub fn reduce_vec(&self, a: &mut [u64]) {
        self.arch
            .dispatch(|| verified::reduce_vec(a, self.p, self.barrett_hi, self.barrett_lo))
    }

    /// Center a value modulo p as i64 in constant time.
    ///
    /// The output is in the interval `[-p/2, p/2)`.
    /// Aborts if `a >= p` in debug mode.
    #[must_use]
    pub const fn center(&self, a: u64) -> i64 {
        debug_assert!(a < self.p);

        let threshold = (self.p + 1) >> 1;
        let cond = a >= threshold;
        let on_true = (a as i64).wrapping_sub(self.p as i64) as u64;
        let on_false = a;

        const_time_cond_select(on_true, on_false, cond) as i64
    }

    /// Center a vector in constant time.
    #[must_use]
    pub fn center_vec(&self, a: &[u64]) -> Vec<i64> {
        self.arch
            .dispatch(|| a.iter().map(|ai| self.center(*ai)).collect_vec())
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
    /// Aborts if a >= p or n >= p in debug mode.
    #[must_use]
    pub fn pow(&self, a: u64, n: u64) -> u64 {
        debug_assert!(a < self.p && n < self.p);

        if n == 0 {
            1
        } else if n == 1 {
            a
        } else {
            let mut r = a;
            let mut i = (62 - n.leading_zeros()) as isize;
            while i >= 0 {
                r = self.mul(r, r);
                if (n >> i) & 1 == 1 {
                    r = self.mul(r, a);
                }
                i -= 1;
            }
            r
        }
    }

    /// Modular inversion in variable time.
    ///
    /// Returns None if p is not prime or a = 0.
    /// Aborts if a >= p in debug mode.
    #[must_use]
    pub fn inv(&self, a: u64) -> std::option::Option<u64> {
        if !is_prime(self.p) || a == 0 {
            None
        } else {
            let r = self.pow(a, self.p - 2);
            debug_assert_eq!(self.mul(a, r), 1);
            Some(r)
        }
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

    /// Returns a random vector.
    pub fn random_vec<R: RngCore + CryptoRng>(&self, size: usize, rng: &mut R) -> Vec<u64> {
        rng.sample_iter(self.distribution).take(size).collect_vec()
    }

    pub(crate) fn random_into<R: RngCore + CryptoRng>(&self, out: &mut [u64], rng: &mut R) {
        for (output, sample) in out.iter_mut().zip(rng.sample_iter(self.distribution)) {
            *output = sample;
        }
    }

    pub(crate) fn serialize_into(&self, a: &[u64], out: &mut Vec<u8>) {
        let p_nbits = 64 - (self.p - 1).leading_zeros() as usize;
        transcode_to_bytes_into(a, p_nbits, out).unwrap();
    }

    /// Length of the serialization of a vector of size `size`.
    ///
    /// Panics if the size is not a multiple of 8.
    #[must_use]
    pub const fn serialization_length(&self, size: usize) -> usize {
        assert!(size.is_multiple_of(8));
        let p_nbits = 64 - (self.p - 1).leading_zeros() as usize;
        p_nbits * size / 8
    }

    /// Pack words at this modulus' bit width, rejecting truncation and padding
    /// the final byte with zeros. This does not validate that every word is <
    /// p.
    pub fn serialize_vec(
        &self,
        a: &[u64],
    ) -> std::result::Result<Vec<u8>, fhe_util::TranscodeError> {
        let p_nbits = 64 - (self.p - 1).leading_zeros() as usize;
        transcode_to_bytes(a, p_nbits)
    }

    /// Decode every bit at this modulus' width, zero-extending a final partial
    /// word. This does not reduce or validate residues against p.
    pub fn deserialize_vec(
        &self,
        b: &[u8],
    ) -> std::result::Result<Vec<u64>, fhe_util::TranscodeError> {
        let p_nbits = 64 - (self.p - 1).leading_zeros() as usize;
        transcode_from_bytes(b, p_nbits)
    }
}

#[cfg(test)]
mod tests {
    use super::Modulus;
    use itertools::{Itertools, izip};
    use proptest::collection::vec as prop_vec;
    use proptest::prelude::{BoxedStrategy, Just, Strategy, any};
    use rand::{Rng as RngCore, rng};

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
        const REDUCED: u64 = Modulus::reduce1(u64::MAX - 2, (1 << 63) - 1);
        const REDUCED_VT: u64 = unsafe { Modulus::reduce1_vt(u64::MAX - 2, (1 << 63) - 1) };
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
            // 63 and 64-bit integers do not work.
            prop_assert!(Modulus::new(p | (1u64 << 62)).is_err());
            prop_assert!(Modulus::new(p | (1u64 << 63)).is_err());

            // p = 0 & 1 do not work.
            prop_assert!(Modulus::new(0u64).is_err());
            prop_assert!(Modulus::new(1u64).is_err());

            // Otherwise, all moduli should work.
            prop_assume!(p >> 2 >= 2);
            let q = Modulus::new(p >> 2);
            prop_assert!(q.is_ok());
            prop_assert_eq!(*q.unwrap(), p >> 2);
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
        fn mul_shoup(p in valid_moduli(), mut a: u64, mut b: u64) {
            a = p.reduce(a);
            b = p.reduce(b);

            // Compute shoup representation
            let b_shoup = p.shoup(b);

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
                prop_assume!(a != b);
                prop_assert!(std::panic::catch_unwind(|| p.mul_shoup(a, a, b_shoup)).is_err());
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
            prop_assert_eq!(a.clone(), izip!(b.iter(), c.iter()).map(|(bi, ci)| p.mul(*ci, *bi)).collect_vec());
            a.clone_from(&c);
            unsafe { p.mul_vec_vt(&mut a, &b); }
            prop_assert_eq!(a, izip!(b.iter(), c.iter()).map(|(bi, ci)| p.mul(*ci, *bi)).collect_vec());
        }

        #[test]
        fn scalar_mul_vec(p in valid_moduli(), mut a: Vec<u64>, mut b: u64) {
            p.reduce_vec(&mut a);
            b = p.reduce(b);
            let c = a.clone();

            p.scalar_mul_vec(&mut a, b);
            prop_assert_eq!(a.clone(), c.iter().map(|ci| p.mul(*ci, b)).collect_vec());

            a.clone_from(&c);
            unsafe { p.scalar_mul_vec_vt(&mut a, b) }
            prop_assert_eq!(a, c.iter().map(|ci| p.mul(*ci, b)).collect_vec());
        }

        #[test]
        fn mul_shoup_vec(p in valid_moduli(), (mut a, mut b) in vecs()) {
            p.reduce_vec(&mut a);
            p.reduce_vec(&mut b);
            let b_shoup = p.shoup_vec(&b);
            let c = a.clone();
            p.mul_shoup_vec(&mut a, &b, &b_shoup);
            prop_assert_eq!(a.clone(), izip!(b.iter(), c.iter()).map(|(bi, ci)| p.mul(*ci, *bi)).collect_vec());
            a.clone_from(&c);
            unsafe { p.mul_shoup_vec_vt(&mut a, &b, &b_shoup) }
            prop_assert_eq!(a, izip!(b.iter(), c.iter()).map(|(bi, ci)| p.mul(*ci, *bi)).collect_vec());
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
        fn random_vec(p in valid_moduli(), size in 1..1000usize) {
            let mut rng = rng();

            let v = p.random_vec(size, &mut rng);
            prop_assert_eq!(v.len(), size);

            let w = p.random_vec(size, &mut rng);
            prop_assert_eq!(w.len(), size);

            if (*p).leading_zeros() <= 30 {
                prop_assert_ne!(v, w); // This will hold with probability at least 2^(-30)
            }
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
            let a = p.reduce(a);
            let b = p.center(a);
            if a >= (*p + 1) >> 1 {
                prop_assert_eq!(b, (a as i64) - (*p as i64));
            } else {
                prop_assert_eq!(b, a as i64);
            }
            prop_assert_eq!(p.reduce_i64(b), a);
        }

        #[test]
        fn center_vec(p in valid_moduli(), a: Vec<u64>) {
             let mut a = a.clone();
             p.reduce_vec(&mut a);
             let b = p.center_vec(&a);
             prop_assert_eq!(b.len(), a.len());
             for (ai, bi) in izip!(a.iter(), b.iter()) {
                 prop_assert_eq!(p.center(*ai), *bi);
             }
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
            a = p.reduce(a);
            b = p.reduce(b);

            prop_assert_eq!(p.pow(a, 0), 1);
            prop_assert_eq!(p.pow(a, 1), a);
            if *p > 2 {
                prop_assert_eq!(p.pow(a, 2), p.mul(a, a));
            }

            let b_small = b % 1000;
            let mut r = 1;
            for _ in 0..b_small {
                r = p.mul(r, a);
            }
            prop_assert_eq!(p.pow(a, b_small), r);

            #[cfg(debug_assertions)]
            {
                prop_assert!(std::panic::catch_unwind(|| p.pow(*p, 1)).is_err());
                prop_assert!(std::panic::catch_unwind(|| p.pow(*p << 1, 1)).is_err());
                prop_assert!(std::panic::catch_unwind(|| p.pow(0, *p)).is_err());
                prop_assert!(std::panic::catch_unwind(|| p.pow(0, *p << 1)).is_err());
            }
        }
    }

    // TODO: Make a proptest.
    #[test]
    fn inv() {
        let ntests = 100;
        let mut rng = rand::rng();

        for p in [2u64, 3, 17, 1987, 4611686018326724609] {
            let q = Modulus::new(p).unwrap();

            assert!(q.inv(0).is_none());
            assert_eq!(q.inv(1).unwrap(), 1);
            assert_eq!(q.inv(p - 1).unwrap(), p - 1);

            #[cfg(debug_assertions)]
            {
                assert!(std::panic::catch_unwind(|| q.inv(p)).is_err());
                assert!(std::panic::catch_unwind(|| q.inv(p << 1)).is_err());
            }

            for _ in 0..ntests {
                let a = rng.next_u64() % p;
                let b = q.inv(a);

                if a == 0 {
                    assert!(b.is_none())
                } else {
                    assert!(b.is_some());
                    assert_eq!(q.mul(a, b.unwrap()), 1)
                }
            }
        }
    }

    #[test]
    fn centered_boundary_is_correct_for_odd_and_even_moduli() {
        for p in [2u64, 3, 17, 1024, 1153] {
            let modulus = Modulus::new(p).unwrap();
            for value in 0..p {
                let centered = modulus.center(value);
                assert!(2 * centered >= -(p as i64));
                assert!(2 * centered < p as i64);
                assert_eq!(modulus.reduce_i64(centered), value);
            }
        }
        assert_eq!(Modulus::new(1153).unwrap().center(576), 576);
    }
}
