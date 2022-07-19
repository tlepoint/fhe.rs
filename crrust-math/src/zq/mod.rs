#![warn(missing_docs, unused_imports)]

//! Ring operations for primes up to 62 bits.

use num_bigint::BigUint;
use num_traits::cast::ToPrimitive;

/// Structure holding a prime modulus up to 62 bits.
#[derive(Debug, Clone, PartialEq)]
pub struct Modulus {
    p: u64,
    p_twice: u64,
    barrett: u128,
}

impl Modulus {
    /// Create a modulus from a prime number of at most 63 bits.
    pub fn new(p: u64) -> std::option::Option<Self> {
        if !Self::is_prime(p) || (p >> 63) != 0 {
            None
        } else {
            let p_twice = 2 * p;
            let barrett = ((BigUint::from(1u64) << 128u64) / p).to_u128().unwrap(); // 2^128 / p

            Some(Self {
                p,
                p_twice,
                barrett,
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
    /// Aborts if a or b >= p in debug mode.
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
    /// Aborts if a or b >= p in debug mode.
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
    /// Aborts if a or b >= p in debug mode.
    pub fn mul(&self, a: u64, b: u64) -> u64 {
        debug_assert!(a < self.p && b < self.p);

        self.reduce_u128((a as u128) * (b as u128))
    }

    /// Modular reduction of a in variable time.
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
        for p in [0u64, 1, 9223372036854775837] {
            assert!(Modulus::new(p).is_none())
        }

        // TODO: Verifies that it fails when p is not prime.
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
}
