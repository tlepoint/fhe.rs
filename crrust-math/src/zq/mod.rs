#![warn(missing_docs, unused_imports)]

//! Ring operations for moduli up to 63 bits.

use num_bigint::BigUint;
use num_traits::cast::ToPrimitive;

/// Structure holding a modulus up to 63 bits.
#[derive(Debug, Clone, PartialEq)]
pub struct Modulus {
    p: u64,
    p_twice: u64,
    barrett: u128,
    is_prime: bool,
}

impl Modulus {
    /// Create a modulus from a prime number of at most 63 bits.
    pub fn new(p: u64) -> std::option::Option<Self> {
        if p < 2 || (p >> 63) != 0 {
            None
        } else {
            let p_twice = 2 * p;
            let barrett = ((BigUint::from(1u64) << 128u64) / p).to_u128().unwrap(); // 2^128 / p

            Some(Self {
                p,
                p_twice,
                barrett,
                is_prime: Self::is_prime(p),
            })
        }
    }

    /// Returns whether the modulus p supports the Number Theoretic Transform of size 2^n.
    /// Aborts if n >= 63.
    pub fn supports_ntt(&self, n: usize) -> bool {
        assert!(n <= 62);

        self.p % (1 << (n + 1)) == 1
    }

    /// Modular addition of a and b in variable time.
    ///
    /// Aborts if a >= p or b >= p in debug mode.
    pub fn add(&self, a: u64, b: u64) -> u64 {
        debug_assert!(a < self.p && b < self.p);

        let r = a + b;
        if r >= self.p {
            r - self.p
        } else {
            r
        }
    }

    /// Modular subtraction of a and b in variable time.
    ///
    /// Aborts if a >= p or b >= p in debug mode.
    pub fn sub(&self, a: u64, b: u64) -> u64 {
        debug_assert!(a < self.p && b < self.p);

        let r = a + self.p - b;
        if r >= self.p {
            r - self.p
        } else {
            r
        }
    }

    /// Modular multiplication of a and b in variable time.
    ///
    /// Aborts if a >= p or b >= p in debug mode.
    pub fn mul(&self, a: u64, b: u64) -> u64 {
        debug_assert!(a < self.p && b < self.p);

        self.reduce_u128((a as u128) * (b as u128))
    }

    /// Modular exponentiation in variable time.
    ///
    /// Aborts if a >= p or n >= p in debug mode.
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
    /// Returns None if p is not prime.
    /// Aborts if a >= p in debug mode.
    pub fn inv(&self, a: u64) -> std::option::Option<u64> {
        if !self.is_prime || a == 0 {
            None
        } else {
            let r = self.pow(a, self.p - 2);
            debug_assert_eq!(self.mul(a, r), 1);
            Some(r)
        }
    }

    /// Modular reduction of a u128 in variable time.
    fn reduce_u128(&self, a: u128) -> u64 {
        let r = self.lazy_reduce_u128(a);
        if r >= self.p {
            r - self.p
        } else {
            r
        }
    }

    /// Lazy modular reduction of a in variable time. The output is in the interval [0, 2 * p).
    fn lazy_reduce_u128(&self, a: u128) -> u64 {
        let a_lo = a as u64;
        let a_hi = (a >> 64) as u64;
        let p_lo_lo = ((a_lo as u128) * ((self.barrett as u64) as u128)) >> 64;
        let p_lo_hi = (a_lo as u128) * (self.barrett >> 64);
        let p_hi_lo = (a_hi as u128) * ((self.barrett as u64) as u128);
        let mut q = (p_lo_hi + p_hi_lo + p_lo_lo) >> 64;
        q += (a_hi as u128) * (self.barrett >> 64);

        let r = (a - q * (self.p as u128)) as u64;
        debug_assert!(r < 2 * self.p);

        r
    }

    fn is_prime(p: u64) -> bool {
        // TODO: To implement
        !(p < 2 || (p != 2 && p & 1 == 0))
    }
}

#[cfg(test)]
mod tests {
    use super::Modulus;
    use rand::RngCore;
    use std::panic::UnwindSafe;

    // Redefine catch_unwind to silence the panic.
    pub fn catch_unwind<F, R>(f: F) -> std::thread::Result<R>
    where
        F: FnOnce() -> R + UnwindSafe,
    {
        let prev_hook = std::panic::take_hook();
        std::panic::set_hook(Box::new(|_| {}));
        let r = std::panic::catch_unwind(f);
        std::panic::set_hook(prev_hook);
        r
    }

    #[test]
    fn test_constructor() {
        // Initialize using a prime within the correct bound.
        for p in [2u64, 3, 17, 1987, 4611686018326724609] {
            assert!(Modulus::new(p).is_some())
        }

        // Initialize using out-of-bound numbers.
        for p in [0u64, 1, 9223372036854775837 /* = next_prime(2**63) */] {
            assert!(Modulus::new(p).is_none())
        }
    }

    #[test]
    fn test_add() {
        let ntests = 100;
        let mut rng = rand::thread_rng();

        for p in [2u64, 3, 17, 1987, 4611686018326724609] {
            let q = Modulus::new(p).unwrap();

            assert_eq!(q.add(0, 1), 1);
            assert_eq!(q.add(1, 1), 2 % p);
            assert_eq!(q.add(p - 1, 1), 0);
            assert_eq!(q.add(p - 1, 2 % p), 1);

            assert!(catch_unwind(|| q.add(p, 1)).is_err());
            assert!(catch_unwind(|| q.add(p << 1, 1)).is_err());
            assert!(catch_unwind(|| q.add(0, p)).is_err());
            assert!(catch_unwind(|| q.add(0, p << 1)).is_err());

            for _ in 0..ntests {
                let a = rng.next_u64() % p;
                let b = rng.next_u64() % p;
                assert_eq!(q.add(a, b), (a + b) % p);
            }
        }
    }

    #[test]
    fn test_sub() {
        let ntests = 100;
        let mut rng = rand::thread_rng();

        for p in [2u64, 3, 17, 1987, 4611686018326724609] {
            let q = Modulus::new(p).unwrap();

            assert_eq!(q.sub(0, 1), p - 1);
            assert_eq!(q.sub(1, 1), 0);
            assert_eq!(q.sub(p - 1, 1), p - 2);
            assert_eq!(q.sub(0, p - 1), 1);

            assert!(catch_unwind(|| q.sub(p, 1)).is_err());
            assert!(catch_unwind(|| q.sub(p << 1, 1)).is_err());
            assert!(catch_unwind(|| q.sub(0, p)).is_err());
            assert!(catch_unwind(|| q.sub(0, p << 1)).is_err());

            for _ in 0..ntests {
                let a = rng.next_u64() % p;
                let b = rng.next_u64() % p;
                assert_eq!(q.sub(a, b), (a + p - b) % p);
            }
        }
    }

    #[test]
    fn test_mul() {
        let ntests = 100;
        let mut rng = rand::thread_rng();

        for p in [2u64, 3, 17, 1987, 4611686018326724609] {
            let q = Modulus::new(p).unwrap();

            assert_eq!(q.mul(0, 1), 0);
            assert_eq!(q.mul(1, 1), 1);
            assert_eq!(q.mul(2 % p, 3 % p), 6 % p);
            assert_eq!(q.mul(p - 1, 1), p - 1);
            assert_eq!(q.mul(p - 1, 2 % p), p - 2);

            assert!(catch_unwind(|| q.mul(p, 1)).is_err());
            assert!(catch_unwind(|| q.mul(p << 1, 1)).is_err());
            assert!(catch_unwind(|| q.mul(0, p)).is_err());
            assert!(catch_unwind(|| q.mul(0, p << 1)).is_err());

            for _ in 0..ntests {
                let a = rng.next_u64() % p;
                let b = rng.next_u64() % p;
                let ab = (a as u128) * (b as u128);
                assert_eq!(q.mul(a, b), (ab % (p as u128)) as u64);
            }
        }
    }

    #[test]
    fn test_pow() {
        let ntests = 10;
        let mut rng = rand::thread_rng();

        for p in [2u64, 3, 17, 1987, 4611686018326724609] {
            let q = Modulus::new(p).unwrap();

            assert_eq!(q.pow(p - 1, 0), 1);
            assert_eq!(q.pow(p - 1, 1), p - 1);
            assert_eq!(q.pow(p - 1, 2 % p), 1);
            assert_eq!(q.pow(1, p - 2), 1);
            assert_eq!(q.pow(1, p - 1), 1);

            assert!(catch_unwind(|| q.pow(p, 1)).is_err());
            assert!(catch_unwind(|| q.pow(p << 1, 1)).is_err());
            assert!(catch_unwind(|| q.pow(0, p)).is_err());
            assert!(catch_unwind(|| q.pow(0, p << 1)).is_err());

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

    #[test]
    fn test_inv() {
        let ntests = 100;
        let mut rng = rand::thread_rng();

        for p in [2u64, 3, 17, 1987, 4611686018326724609] {
            let q = Modulus::new(p).unwrap();

            assert!(q.inv(0).is_none());
            assert_eq!(q.inv(1).unwrap(), 1);
            assert_eq!(q.inv(p - 1).unwrap(), p - 1);

            assert!(catch_unwind(|| q.inv(p)).is_err());
            assert!(catch_unwind(|| q.inv(p << 1)).is_err());

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
