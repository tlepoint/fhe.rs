//! Number-Theoretic Transform in ZZ_q.

use fhe_util::is_prime;

mod native;

#[cfg(any(feature = "concrete-ntt", feature = "concrete-ntt-nightly"))]
mod concrete;

#[cfg(any(feature = "concrete-ntt", feature = "concrete-ntt-nightly"))]
pub use concrete::NttOperator;
#[cfg(not(any(feature = "concrete-ntt", feature = "concrete-ntt-nightly")))]
pub use native::NttOperator;

/// Returns whether a modulus p is prime and supports the Number Theoretic
/// Transform of size n.
///
/// Aborts if n is not a power of 2 that is >= 8.
pub(crate) fn supports_ntt(p: u64, n: usize) -> bool {
    assert!(n >= 8 && n.is_power_of_two());

    p % ((n as u64) << 1) == 1 && is_prime(p)
}

#[cfg(test)]
mod tests {
    use rand::thread_rng;

    use super::{supports_ntt, NttOperator};
    use crate::zq::Modulus;

    #[test]
    fn constructor() {
        for size in [32, 1024] {
            for p in [1153, 4611686018326724609] {
                let q = Modulus::new(p).unwrap();
                let supports_ntt = supports_ntt(p, size);

                let op = NttOperator::new(&q, size);

                if supports_ntt {
                    assert!(op.is_some());
                } else {
                    assert!(op.is_none());
                }
            }
        }
    }

    #[test]
    fn bijection() {
        let ntests = 100;
        let mut rng = thread_rng();

        for size in [32, 1024] {
            for p in [1153, 4611686018326724609] {
                let q = Modulus::new(p).unwrap();

                if supports_ntt(p, size) {
                    let op = NttOperator::new(&q, size).unwrap();

                    for _ in 0..ntests {
                        let mut a = q.random_vec(size, &mut rng);
                        let a_clone = a.clone();
                        let mut b = a.clone();

                        op.forward(&mut a);
                        assert_ne!(a, a_clone);

                        unsafe { op.forward_vt(b.as_mut_ptr()) }
                        assert_eq!(a, b);

                        op.backward(&mut a);
                        assert_eq!(a, a_clone);

                        unsafe { op.backward_vt(b.as_mut_ptr()) }
                        assert_eq!(a, b);
                    }
                }
            }
        }
    }

    #[test]
    fn forward_lazy() {
        let ntests = 100;
        let mut rng = thread_rng();

        for size in [32, 1024] {
            for p in [1153, 4611686018326724609] {
                let q = Modulus::new(p).unwrap();

                if supports_ntt(p, size) {
                    let op = NttOperator::new(&q, size).unwrap();

                    for _ in 0..ntests {
                        let mut a = q.random_vec(size, &mut rng);
                        let mut a_lazy = a.clone();

                        op.forward(&mut a);

                        unsafe {
                            op.forward_vt_lazy(a_lazy.as_mut_ptr());
                            q.reduce_vec(&mut a_lazy);
                        }

                        assert_eq!(a, a_lazy);
                    }
                }
            }
        }
    }
}
