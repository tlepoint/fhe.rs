#![warn(missing_docs, unused_imports)]

//! Ring operations for moduli up to 62 bits.

pub mod primes;

use std::ops::Deref;

use crate::errors::{Error, Result};
use fhe_util::{is_prime, transcode_from_bytes, transcode_to_bytes};
use itertools::{izip, Itertools};
use num_bigint::BigUint;
use num_traits::cast::ToPrimitive;
use pulp::Arch;
use rand::{distr::Uniform, CryptoRng, Rng, RngCore};

/// cond ? on_true : on_false
const fn const_time_cond_select(on_true: u64, on_false: u64, cond: bool) -> u64 {
    let mask = -(cond as i64) as u64;
    let diff = on_true ^ on_false;
    (diff & mask) ^ on_false
}

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
        Self::reduce1(a + b, self.p)
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
        Self::reduce1_vt(a + b, self.p)
    }

    /// Performs the modular subtraction of a and b in constant time.
    /// Aborts if a >= p or b >= p in debug mode.
    #[must_use]
    pub const fn sub(&self, a: u64, b: u64) -> u64 {
        debug_assert!(a < self.p && b < self.p);
        Self::reduce1(a + self.p - b, self.p)
    }

    /// Performs the modular subtraction of a and b in constant time.
    /// Aborts if a >= p or b >= p in debug mode.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being subtracted.
    const unsafe fn sub_vt(&self, a: u64, b: u64) -> u64 {
        debug_assert!(a < self.p && b < self.p);
        Self::reduce1_vt(a + self.p - b, self.p)
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
        Self::reduce1_vt(self.lazy_reduce_u128((a as u128) * (b as u128)), self.p)
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

        self.reduce_opt_u128_vt((a as u128) * (b as u128))
    }

    /// Modular negation in constant time.
    ///
    /// Aborts if a >= p in debug mode.
    #[must_use]
    pub const fn neg(&self, a: u64) -> u64 {
        debug_assert!(a < self.p);
        Self::reduce1(self.p - a, self.p)
    }

    /// Modular negation in variable time.
    /// Aborts if a >= p in debug mode.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the value being negated.
    const unsafe fn neg_vt(&self, a: u64) -> u64 {
        debug_assert!(a < self.p);
        Self::reduce1_vt(self.p - a, self.p)
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
        Self::reduce1_vt(self.lazy_mul_shoup(a, b, b_shoup), self.p)
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
        self.arch.dispatch(|| {
            izip!(a.iter_mut(), b.iter()).for_each(|(ai, bi)| *ai = self.add(*ai, *bi))
        })
    }

    /// Modular addition of vectors in place in variable time.
    /// Aborts if a and b differ in size, and if any of their values is >= p in
    /// debug mode.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being added.
    pub unsafe fn add_vec_vt(&self, a: &mut [u64], b: &[u64]) {
        let n = a.len();
        debug_assert_eq!(n, b.len());

        let p = self.p;
        macro_rules! add_at {
            ($idx:expr) => {
                *a.get_unchecked_mut($idx) =
                    Self::reduce1_vt(*a.get_unchecked_mut($idx) + *b.get_unchecked($idx), p);
            };
        }

        if n % 16 == 0 {
            self.arch.dispatch(|| {
                for i in 0..n / 16 {
                    add_at!(16 * i);
                    add_at!(16 * i + 1);
                    add_at!(16 * i + 2);
                    add_at!(16 * i + 3);
                    add_at!(16 * i + 4);
                    add_at!(16 * i + 5);
                    add_at!(16 * i + 6);
                    add_at!(16 * i + 7);
                    add_at!(16 * i + 8);
                    add_at!(16 * i + 9);
                    add_at!(16 * i + 10);
                    add_at!(16 * i + 11);
                    add_at!(16 * i + 12);
                    add_at!(16 * i + 13);
                    add_at!(16 * i + 14);
                    add_at!(16 * i + 15);
                }
            })
        } else {
            self.arch.dispatch(|| {
                izip!(a.iter_mut(), b.iter()).for_each(|(ai, bi)| *ai = self.add_vt(*ai, *bi))
            })
        }
    }

    /// Modular subtraction of vectors in place in constant time.
    ///
    /// Aborts if a and b differ in size, and if any of their values is >= p in
    /// debug mode.
    pub fn sub_vec(&self, a: &mut [u64], b: &[u64]) {
        debug_assert_eq!(a.len(), b.len());
        self.arch.dispatch(|| {
            izip!(a.iter_mut(), b.iter()).for_each(|(ai, bi)| *ai = self.sub(*ai, *bi))
        })
    }

    /// Modular subtraction of vectors in place in variable time.
    /// Aborts if a and b differ in size, and if any of their values is >= p in
    /// debug mode.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being subtracted.
    pub unsafe fn sub_vec_vt(&self, a: &mut [u64], b: &[u64]) {
        let n = a.len();
        debug_assert_eq!(n, b.len());

        let p = self.p;
        macro_rules! sub_at {
            ($idx:expr) => {
                *a.get_unchecked_mut($idx) =
                    Self::reduce1_vt(p + *a.get_unchecked_mut($idx) - *b.get_unchecked($idx), p);
            };
        }

        if n % 16 == 0 {
            self.arch.dispatch(|| {
                for i in 0..n / 16 {
                    sub_at!(16 * i);
                    sub_at!(16 * i + 1);
                    sub_at!(16 * i + 2);
                    sub_at!(16 * i + 3);
                    sub_at!(16 * i + 4);
                    sub_at!(16 * i + 5);
                    sub_at!(16 * i + 6);
                    sub_at!(16 * i + 7);
                    sub_at!(16 * i + 8);
                    sub_at!(16 * i + 9);
                    sub_at!(16 * i + 10);
                    sub_at!(16 * i + 11);
                    sub_at!(16 * i + 12);
                    sub_at!(16 * i + 13);
                    sub_at!(16 * i + 14);
                    sub_at!(16 * i + 15);
                }
            })
        } else {
            self.arch.dispatch(|| {
                izip!(a.iter_mut(), b.iter()).for_each(|(ai, bi)| *ai = self.sub_vt(*ai, *bi))
            })
        }
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
                .for_each(|ai| *ai = self.mul_shoup_vt(*ai, b, b_shoup))
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
                izip!(a.iter_mut(), b.iter()).for_each(|(ai, bi)| *ai = self.mul_opt_vt(*ai, *bi))
            })
        } else {
            self.arch.dispatch(|| {
                izip!(a.iter_mut(), b.iter()).for_each(|(ai, bi)| *ai = self.mul_vt(*ai, *bi))
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
            izip!(a.iter_mut(), b.iter(), b_shoup.iter())
                .for_each(|(ai, bi, bi_shoup)| *ai = self.mul_shoup_vt(*ai, *bi, *bi_shoup))
        })
    }

    /// Reduce a vector in place in constant time.
    pub fn reduce_vec(&self, a: &mut [u64]) {
        self.arch
            .dispatch(|| a.iter_mut().for_each(|ai| *ai = self.reduce(*ai)))
    }

    /// Center a value modulo p as i64 in variable time.
    /// TODO: To test and to make constant time?
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the value being centered.
    const unsafe fn center_vt(&self, a: u64) -> i64 {
        debug_assert!(a < self.p);

        if a >= self.p >> 1 {
            (a as i64) - (self.p as i64)
        } else {
            a as i64
        }
    }

    /// Center a vector in variable time.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being centered.
    #[must_use]
    pub unsafe fn center_vec_vt(&self, a: &[u64]) -> Vec<i64> {
        self.arch
            .dispatch(|| a.iter().map(|ai| self.center_vt(*ai)).collect_vec())
    }

    /// Reduce a vector in place in variable time.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being reduced.
    pub unsafe fn reduce_vec_vt(&self, a: &mut [u64]) {
        self.arch
            .dispatch(|| a.iter_mut().for_each(|ai| *ai = self.reduce_vt(*ai)))
    }

    /// Modular reduction of a i64 in constant time.
    const fn reduce_i64(&self, a: i64) -> u64 {
        self.reduce_u128((((self.p as i128) << 64) + (a as i128)) as u128)
    }

    /// Modular reduction of a i64 in variable time.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being reduced.
    const unsafe fn reduce_i64_vt(&self, a: i64) -> u64 {
        self.reduce_u128_vt((((self.p as i128) << 64) + (a as i128)) as u128)
    }

    /// Reduce a vector in place in constant time.
    #[must_use]
    pub fn reduce_vec_i64(&self, a: &[i64]) -> Vec<u64> {
        self.arch
            .dispatch(|| a.iter().map(|ai| self.reduce_i64(*ai)).collect_vec())
    }

    /// Reduce a vector in place in variable time.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being reduced.
    #[must_use]
    pub unsafe fn reduce_vec_i64_vt(&self, a: &[i64]) -> Vec<u64> {
        self.arch
            .dispatch(|| a.iter().map(|ai| self.reduce_i64_vt(*ai)).collect())
    }

    /// Reduce a vector in constant time.
    #[must_use]
    pub fn reduce_vec_new(&self, a: &[u64]) -> Vec<u64> {
        self.arch
            .dispatch(|| a.iter().map(|ai| self.reduce(*ai)).collect())
    }

    /// Reduce a vector in variable time.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being reduced.
    #[must_use]
    pub unsafe fn reduce_vec_new_vt(&self, a: &[u64]) -> Vec<u64> {
        self.arch
            .dispatch(|| a.iter().map(|bi| self.reduce_vt(*bi)).collect())
    }

    /// Modular negation of a vector in place in constant time.
    ///
    /// Aborts if any of the values in the vector is >= p in debug mode.
    pub fn neg_vec(&self, a: &mut [u64]) {
        self.arch
            .dispatch(|| a.iter_mut().for_each(|ai| *ai = self.neg(*ai)))
    }

    /// Modular negation of a vector in place in variable time.
    /// Aborts if any of the values in the vector is >= p in debug mode.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the values being negated.
    pub unsafe fn neg_vec_vt(&self, a: &mut [u64]) {
        self.arch
            .dispatch(|| a.iter_mut().for_each(|ai| *ai = self.neg_vt(*ai)))
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
        Self::reduce1(self.lazy_reduce_u128(a), self.p)
    }

    /// Modular reduction of a u128 in variable time.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the value being reduced.
    #[must_use]
    pub const unsafe fn reduce_u128_vt(&self, a: u128) -> u64 {
        Self::reduce1_vt(self.lazy_reduce_u128(a), self.p)
    }

    /// Modular reduction of a u64 in constant time.
    #[must_use]
    pub const fn reduce(&self, a: u64) -> u64 {
        Self::reduce1(self.lazy_reduce(a), self.p)
    }

    /// Modular reduction of a u64 in variable time.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the value being reduced.
    #[must_use]
    pub const unsafe fn reduce_vt(&self, a: u64) -> u64 {
        Self::reduce1_vt(self.lazy_reduce(a), self.p)
    }

    /// Optimized modular reduction of a u128 in constant time.
    #[must_use]
    pub const fn reduce_opt_u128(&self, a: u128) -> u64 {
        debug_assert!(self.supports_opt);
        Self::reduce1(self.lazy_reduce_opt_u128(a), self.p)
    }

    /// Optimized modular reduction of a u128 in constant time.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the value being reduced.
    pub(crate) const unsafe fn reduce_opt_u128_vt(&self, a: u128) -> u64 {
        debug_assert!(self.supports_opt);
        Self::reduce1_vt(self.lazy_reduce_opt_u128(a), self.p)
    }

    /// Optimized modular reduction of a u64 in constant time.
    #[must_use]
    pub const fn reduce_opt(&self, a: u64) -> u64 {
        Self::reduce1(self.lazy_reduce_opt(a), self.p)
    }

    /// Optimized modular reduction of a u64 in variable time.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the value being reduced.
    #[must_use]
    pub const unsafe fn reduce_opt_vt(&self, a: u64) -> u64 {
        Self::reduce1_vt(self.lazy_reduce_opt(a), self.p)
    }

    /// Return x mod p in constant time.
    /// Aborts if x >= 2 * p in debug mode.
    pub(crate) const fn reduce1(x: u64, p: u64) -> u64 {
        debug_assert!(p >> 63 == 0);
        debug_assert!(x < 2 * p);

        let r = const_time_cond_select(x, x.wrapping_sub(p), x < p);

        debug_assert!(r == x % p);

        r
    }

    /// Return x mod p in variable time.
    /// Aborts if x >= 2 * p in debug mode.
    ///
    /// # Safety
    /// This function is not constant time and its timing may reveal information
    /// about the value being reduced.
    #[cfg(any(target_os = "macos", target_feature = "avx2"))]
    pub(crate) const unsafe fn reduce1_vt(x: u64, p: u64) -> u64 {
        debug_assert!(p >> 63 == 0);
        debug_assert!(x < 2 * p);

        if x >= p {
            x - p
        } else {
            x
        }
    }

    #[cfg(all(not(target_os = "macos"), not(target_feature = "avx2")))]
    #[inline]
    pub(crate) const unsafe fn reduce1_vt(x: u64, p: u64) -> u64 {
        Self::reduce1(x, p)
    }

    /// Lazy modular reduction of a in constant time.
    /// The output is in the interval [0, 2 * p).
    #[must_use]
    pub const fn lazy_reduce_u128(&self, a: u128) -> u64 {
        let a_lo = a as u64;
        let a_hi = (a >> 64) as u64;
        let p_lo_lo = ((a_lo as u128) * (self.barrett_lo as u128)) >> 64;
        let p_hi_lo = (a_hi as u128) * (self.barrett_lo as u128);
        let p_lo_hi = (a_lo as u128) * (self.barrett_hi as u128);

        let q = ((p_lo_hi + p_hi_lo + p_lo_lo) >> 64) + (a_hi as u128) * (self.barrett_hi as u128);
        let r = (a - q * (self.p as u128)) as u64;

        debug_assert!((r as u128) < 2 * (self.p as u128));
        debug_assert!(r % self.p == (a % (self.p as u128)) as u64);

        r
    }

    /// Lazy modular reduction of a in constant time.
    /// The output is in the interval [0, 2 * p).
    #[must_use]
    pub const fn lazy_reduce(&self, a: u64) -> u64 {
        let p_lo_lo = ((a as u128) * (self.barrett_lo as u128)) >> 64;
        let p_lo_hi = (a as u128) * (self.barrett_hi as u128);

        let q = (p_lo_hi + p_lo_lo) >> 64;
        let r = (a as u128 - q * (self.p as u128)) as u64;

        debug_assert!((r as u128) < 2 * (self.p as u128));
        debug_assert!(r % self.p == a % self.p);

        r
    }

    /// Lazy optimized modular reduction of a in constant time.
    /// The output is in the interval [0, 2 * p).
    ///
    /// Aborts if the input is >= p ^ 2 in debug mode.
    #[must_use]
    pub const fn lazy_reduce_opt_u128(&self, a: u128) -> u64 {
        debug_assert!(a < (self.p as u128) * (self.p as u128));

        let q = (((self.barrett_lo as u128) * (a >> 64)) + (a << self.leading_zeros)) >> 64;
        let r = (a - q * (self.p as u128)) as u64;

        debug_assert!((r as u128) < 2 * (self.p as u128));
        debug_assert!(r % self.p == (a % (self.p as u128)) as u64);

        r
    }

    /// Lazy optimized modular reduction of a in constant time.
    /// The output is in the interval [0, 2 * p).
    const fn lazy_reduce_opt(&self, a: u64) -> u64 {
        let q = a >> (64 - self.leading_zeros);
        let r = ((a as u128) - (q as u128) * (self.p as u128)) as u64;

        debug_assert!((r as u128) < 2 * (self.p as u128));
        debug_assert!(r % self.p == a % self.p);

        r
    }

    /// Lazy modular reduction of a vector in constant time.
    /// The output coefficients are in the interval [0, 2 * p).
    pub fn lazy_reduce_vec(&self, a: &mut [u64]) {
        if self.supports_opt {
            a.iter_mut().for_each(|ai| *ai = self.lazy_reduce_opt(*ai))
        } else {
            a.iter_mut().for_each(|ai| *ai = self.lazy_reduce(*ai))
        }
    }

    /// Returns a random vector.
    pub fn random_vec<R: RngCore + CryptoRng>(&self, size: usize, rng: &mut R) -> Vec<u64> {
        rng.sample_iter(self.distribution).take(size).collect_vec()
    }

    /// Length of the serialization of a vector of size `size`.
    ///
    /// Panics if the size is not a multiple of 8.
    #[must_use]
    pub const fn serialization_length(&self, size: usize) -> usize {
        assert!(size % 8 == 0);
        let p_nbits = 64 - (self.p - 1).leading_zeros() as usize;
        p_nbits * size / 8
    }

    /// Serialize a vector of elements of length a multiple of 8.
    ///
    /// Panics if the length of the vector is not a multiple of 8.
    #[must_use]
    pub fn serialize_vec(&self, a: &[u64]) -> Vec<u8> {
        let p_nbits = 64 - (self.p - 1).leading_zeros() as usize;
        transcode_to_bytes(a, p_nbits)
    }

    /// Deserialize a vector of bytes into a vector of elements mod p.
    #[must_use]
    pub fn deserialize_vec(&self, b: &[u8]) -> Vec<u64> {
        let p_nbits = 64 - (self.p - 1).leading_zeros() as usize;
        transcode_from_bytes(b, p_nbits)
    }
}

#[cfg(test)]
mod tests {
    use super::{primes, Modulus};
    use itertools::{izip, Itertools};
    use proptest::collection::vec as prop_vec;
    use proptest::prelude::{any, BoxedStrategy, Just, Strategy};
    use rand::{rng, RngCore};

    // Utility functions for the proptests.

    fn valid_moduli() -> impl Strategy<Value = Modulus> {
        any::<u64>().prop_filter_map("filter invalid moduli", |p| Modulus::new(p).ok())
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
            prop_assert_eq!(p.neg(a), (*p - a) % *p);
            unsafe { prop_assert_eq!(p.neg_vt(a), (*p - a) % *p) }

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
            }
        }

        #[test]
        fn sub(p in valid_moduli(), mut a: u64, mut b: u64) {
            a = p.reduce(a);
            b = p.reduce(b);
            prop_assert_eq!(p.sub(a, b), (a + *p - b) % *p);
            unsafe { prop_assert_eq!(p.sub_vt(a, b), (a + *p - b) % *p) }

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
            prop_assert_eq!(p.lazy_reduce(a) % *p, p.reduce(a));
        }

        #[test]
        fn reduce_i64(p in valid_moduli(), a: i64) {
            let b = if a < 0 { p.neg(p.reduce(-a as u64)) } else { p.reduce(a as u64) };
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
            p.add_vec(&mut a, &b);
            prop_assert_eq!(a.clone(), izip!(b.iter(), c.iter()).map(|(bi, ci)| p.add(*bi, *ci)).collect_vec());
            a.clone_from(&c);
            unsafe { p.add_vec_vt(&mut a, &b) }
            prop_assert_eq!(a, izip!(b.iter(), c.iter()).map(|(bi, ci)| p.add(*bi, *ci)).collect_vec());
        }

        #[test]
        fn sub_vec(p in valid_moduli(), (mut a, mut b) in vecs()) {
            p.reduce_vec(&mut a);
            p.reduce_vec(&mut b);
            let c = a.clone();
            p.sub_vec(&mut a, &b);
            prop_assert_eq!(a.clone(), izip!(b.iter(), c.iter()).map(|(bi, ci)| p.sub(*ci, *bi)).collect_vec());
            a.clone_from(&c);
            unsafe { p.sub_vec_vt(&mut a, &b) }
            prop_assert_eq!(a, izip!(b.iter(), c.iter()).map(|(bi, ci)| p.sub(*ci, *bi)).collect_vec());
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
            prop_assert_eq!(b.clone(), a.iter().map(|ai| p.reduce(*ai)).collect_vec());

            b.clone_from(&a);
            unsafe { p.reduce_vec_vt(&mut b) }
            prop_assert_eq!(b, a.iter().map(|ai| p.reduce(*ai)).collect_vec());
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
            prop_assert_eq!(b, a.iter().map(|ai| p.reduce(*ai)).collect_vec());
            prop_assert_eq!(p.reduce_vec_new(&a), unsafe { p.reduce_vec_new_vt(&a) });
        }

        #[test]
        fn reduce_vec_i64(p in valid_moduli(), a: Vec<i64>) {
            let b = p.reduce_vec_i64(&a);
            prop_assert_eq!(b, a.iter().map(|ai| p.reduce_i64(*ai)).collect_vec());
            let b = unsafe { p.reduce_vec_i64_vt(&a) };
            prop_assert_eq!(b, a.iter().map(|ai| p.reduce_i64(*ai)).collect_vec());
        }

        #[test]
        fn neg_vec(p in valid_moduli(), mut a: Vec<u64>) {
            p.reduce_vec(&mut a);
            let mut b = a.clone();
            p.neg_vec(&mut b);
            prop_assert_eq!(b.clone(), a.iter().map(|ai| p.neg(*ai)).collect_vec());
            b.clone_from(&a);
            unsafe { p.neg_vec_vt(&mut b); }
            prop_assert_eq!(b, a.iter().map(|ai| p.neg(*ai)).collect_vec());
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
            let b = p.serialize_vec(&a);
            let c = p.deserialize_vec(&b);
            prop_assert_eq!(a, c);
        }
    }

    // TODO: Make a proptest.
    #[test]
    fn mul_opt() {
        let ntests = 100;
        let mut rng = rand::rng();

        #[allow(clippy::single_element_loop)]
        for p in [4611686018326724609] {
            let q = Modulus::new(p).unwrap();
            assert!(primes::supports_opt(p));

            assert_eq!(q.mul_opt(0, 1), 0);
            assert_eq!(q.mul_opt(1, 1), 1);
            assert_eq!(q.mul_opt(2 % p, 3 % p), 6 % p);
            assert_eq!(q.mul_opt(p - 1, 1), p - 1);
            assert_eq!(q.mul_opt(p - 1, 2 % p), p - 2);

            #[cfg(debug_assertions)]
            {
                assert!(std::panic::catch_unwind(|| q.mul_opt(p, 1)).is_err());
                assert!(std::panic::catch_unwind(|| q.mul_opt(p << 1, 1)).is_err());
                assert!(std::panic::catch_unwind(|| q.mul_opt(0, p)).is_err());
                assert!(std::panic::catch_unwind(|| q.mul_opt(0, p << 1)).is_err());
            }

            for _ in 0..ntests {
                let a = rng.next_u64() % p;
                let b = rng.next_u64() % p;
                assert_eq!(
                    q.mul_opt(a, b),
                    (((a as u128) * (b as u128)) % (p as u128)) as u64
                );
            }
        }
    }

    // TODO: Make a proptest.
    #[test]
    fn pow() {
        let ntests = 10;
        let mut rng = rand::rng();

        for p in [2u64, 3, 17, 1987, 4611686018326724609] {
            let q = Modulus::new(p).unwrap();

            assert_eq!(q.pow(p - 1, 0), 1);
            assert_eq!(q.pow(p - 1, 1), p - 1);
            assert_eq!(q.pow(p - 1, 2 % p), 1);
            assert_eq!(q.pow(1, p - 2), 1);
            assert_eq!(q.pow(1, p - 1), 1);

            #[cfg(debug_assertions)]
            {
                assert!(std::panic::catch_unwind(|| q.pow(p, 1)).is_err());
                assert!(std::panic::catch_unwind(|| q.pow(p << 1, 1)).is_err());
                assert!(std::panic::catch_unwind(|| q.pow(0, p)).is_err());
                assert!(std::panic::catch_unwind(|| q.pow(0, p << 1)).is_err());
            }

            for _ in 0..ntests {
                let a = rng.next_u64() % p;
                let b = (rng.next_u64() % p) % 1000;
                let mut c = b;
                let mut r = 1;
                while c > 0 {
                    r = q.mul(r, a);
                    c -= 1;
                }
                assert_eq!(q.pow(a, b), r);
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
}
