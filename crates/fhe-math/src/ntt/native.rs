use crate::zq::Modulus;
use itertools::Itertools;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use std::iter::successors;

/// Number-Theoretic Transform operator.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NttOperator {
    p: Modulus,
    p_twice: u64,
    size: usize,
    omegas: Box<[u64]>,
    omegas_shoup: Box<[u64]>,
    zetas_inv: Box<[u64]>,
    zetas_inv_shoup: Box<[u64]>,
    size_inv: u64,
    size_inv_shoup: u64,
}

impl NttOperator {
    /// Create an NTT operator given a modulus for a specific size.
    ///
    /// Aborts if the size is not a power of 2 that is >= 8 in debug mode.
    /// Returns None if the modulus does not support the NTT for this specific
    /// size.
    pub fn new(p: &Modulus, size: usize) -> Option<Self> {
        if !super::supports_ntt(p.p, size) {
            None
        } else {
            let size_inv = p.inv(size as u64)?;

            let omega = Self::primitive_root(size, p);
            let omega_inv = p.inv(omega)?;

            let powers = successors(Some(1u64), |n| Some(p.mul(*n, omega)))
                .take(size)
                .collect_vec();
            let powers_inv = successors(Some(omega_inv), |n| Some(p.mul(*n, omega_inv)))
                .take(size)
                .collect_vec();

            let mut omegas = Vec::with_capacity(size);
            let mut zetas_inv = Vec::with_capacity(size);
            for i in 0..size {
                let j = i.reverse_bits() >> (size.leading_zeros() + 1);
                omegas.push(powers[j]);
                zetas_inv.push(powers_inv[j]);
            }

            let omegas_shoup = p.shoup_vec(&omegas);
            let zetas_inv_shoup = p.shoup_vec(&zetas_inv);

            Some(Self {
                p: p.clone(),
                p_twice: p.p * 2,
                size,
                omegas: omegas.into_boxed_slice(),
                omegas_shoup: omegas_shoup.into_boxed_slice(),
                zetas_inv: zetas_inv.into_boxed_slice(),
                zetas_inv_shoup: zetas_inv_shoup.into_boxed_slice(),
                size_inv,
                size_inv_shoup: p.shoup(size_inv),
            })
        }
    }

    /// Compute the forward NTT in place.
    /// Aborts if a is not of the size handled by the operator.
    pub fn forward(&self, a: &mut [u64]) {
        debug_assert_eq!(a.len(), self.size);

        let mut l = self.size >> 1;
        let mut k = 1;
        while l > 0 {
            for chunk in a.chunks_exact_mut(2 * l) {
                let omega = self.omegas[k];
                let omega_shoup = self.omegas_shoup[k];
                k += 1;

                let (left, right) = chunk.split_at_mut(l);
                if l == 1 {
                    // The last level should reduce the output
                    self.butterfly(&mut left[0], &mut right[0], omega, omega_shoup);
                    left[0] = self.reduce3(left[0]);
                    right[0] = self.reduce3(right[0]);
                } else {
                    for (x, y) in left.iter_mut().zip(right.iter_mut()) {
                        self.butterfly(x, y, omega, omega_shoup);
                    }
                }
            }
            l >>= 1;
        }
    }

    /// Compute the backward NTT in place.
    /// Aborts if a is not of the size handled by the operator.
    pub fn backward(&self, a: &mut [u64]) {
        debug_assert_eq!(a.len(), self.size);

        let mut k = 0;
        let mut l = 1;

        while l < self.size {
            for chunk in a.chunks_exact_mut(2 * l) {
                let zeta_inv = self.zetas_inv[k];
                let zeta_inv_shoup = self.zetas_inv_shoup[k];
                k += 1;

                let (left, right) = chunk.split_at_mut(l);
                if l == 1 {
                    self.inv_butterfly(&mut left[0], &mut right[0], zeta_inv, zeta_inv_shoup);
                } else {
                    for (x, y) in left.iter_mut().zip(right.iter_mut()) {
                        self.inv_butterfly(x, y, zeta_inv, zeta_inv_shoup);
                    }
                }
            }
            l <<= 1;
        }

        a.iter_mut()
            .for_each(|ai| *ai = self.p.mul_shoup(*ai, self.size_inv, self.size_inv_shoup));
    }

    /// Compute the forward NTT in place in variable time in a lazily fashion.
    /// This means that the output coefficients may be up to 4 times the
    /// modulus.
    ///
    /// # Safety
    /// This function assumes that a_ptr points to at least `size` elements.
    /// This function is not constant time and its timing may reveal information
    /// about the value being reduced.
    pub(crate) unsafe fn forward_vt_lazy(&self, a_ptr: *mut u64) {
        let mut l = self.size >> 1;
        let mut m = 1;
        let mut k = 1;
        while l > 0 {
            for i in 0..m {
                let omega = *self.omegas.get_unchecked(k);
                let omega_shoup = *self.omegas_shoup.get_unchecked(k);
                k += 1;

                let s = 2 * i * l;
                match l {
                    1 => {
                        self.butterfly_vt(
                            &mut *a_ptr.add(s),
                            &mut *a_ptr.add(s + l),
                            omega,
                            omega_shoup,
                        );
                    }
                    _ => {
                        for j in s..(s + l) {
                            self.butterfly_vt(
                                &mut *a_ptr.add(j),
                                &mut *a_ptr.add(j + l),
                                omega,
                                omega_shoup,
                            );
                        }
                    }
                }
            }
            l >>= 1;
            m <<= 1;
        }
    }

    /// Compute the forward NTT in place in variable time.
    ///
    /// # Safety
    /// This function assumes that a_ptr points to at least `size` elements.
    /// This function is not constant time and its timing may reveal information
    /// about the value being reduced.
    pub unsafe fn forward_vt(&self, a_ptr: *mut u64) {
        self.forward_vt_lazy(a_ptr);
        for i in 0..self.size {
            *a_ptr.add(i) = self.reduce3_vt(*a_ptr.add(i))
        }
    }

    /// Compute the backward NTT in place in variable time.
    ///
    /// # Safety
    /// This function assumes that a_ptr points to at least `size` elements.
    /// This function is not constant time and its timing may reveal information
    /// about the value being reduced.
    pub unsafe fn backward_vt(&self, a_ptr: *mut u64) {
        let mut k = 0;
        let mut m = self.size >> 1;
        let mut l = 1;
        while m > 0 {
            for i in 0..m {
                let s = 2 * i * l;
                let zeta_inv = *self.zetas_inv.get_unchecked(k);
                let zeta_inv_shoup = *self.zetas_inv_shoup.get_unchecked(k);
                k += 1;
                match l {
                    1 => {
                        self.inv_butterfly_vt(
                            &mut *a_ptr.add(s),
                            &mut *a_ptr.add(s + l),
                            zeta_inv,
                            zeta_inv_shoup,
                        );
                    }
                    _ => {
                        for j in s..(s + l) {
                            self.inv_butterfly_vt(
                                &mut *a_ptr.add(j),
                                &mut *a_ptr.add(j + l),
                                zeta_inv,
                                zeta_inv_shoup,
                            );
                        }
                    }
                }
            }
            l <<= 1;
            m >>= 1;
        }

        for i in 0..self.size as isize {
            *a_ptr.offset(i) =
                self.p
                    .mul_shoup(*a_ptr.offset(i), self.size_inv, self.size_inv_shoup)
        }
    }

    /// Reduce a modulo p.
    ///
    /// Aborts if a >= 4 * p.
    const fn reduce3(&self, a: u64) -> u64 {
        debug_assert!(a < 4 * self.p.p);

        let y = Modulus::reduce1(a, self.p_twice);
        Modulus::reduce1(y, self.p.p)
    }

    /// Reduce a modulo p in variable time.
    ///
    /// Aborts if a >= 4 * p.
    const unsafe fn reduce3_vt(&self, a: u64) -> u64 {
        debug_assert!(a < 4 * self.p.p);

        let y = Modulus::reduce1_vt(a, self.p_twice);
        Modulus::reduce1_vt(y, self.p.p)
    }

    /// NTT Butterfly.
    fn butterfly(&self, x: &mut u64, y: &mut u64, w: u64, w_shoup: u64) {
        debug_assert!(*x < 4 * self.p.p);
        debug_assert!(*y < 4 * self.p.p);
        debug_assert!(w < self.p.p);
        debug_assert_eq!(self.p.shoup(w), w_shoup);

        *x = Modulus::reduce1(*x, self.p_twice);
        let t = self.p.lazy_mul_shoup(*y, w, w_shoup);
        *y = *x + self.p_twice - t;
        *x += t;

        debug_assert!(*x < 4 * self.p.p);
        debug_assert!(*y < 4 * self.p.p);
    }

    /// NTT Butterfly in variable time.
    unsafe fn butterfly_vt(&self, x: &mut u64, y: &mut u64, w: u64, w_shoup: u64) {
        debug_assert!(*x < 4 * self.p.p);
        debug_assert!(*y < 4 * self.p.p);
        debug_assert!(w < self.p.p);
        debug_assert_eq!(self.p.shoup(w), w_shoup);

        *x = Modulus::reduce1_vt(*x, self.p_twice);
        let t = self.p.lazy_mul_shoup(*y, w, w_shoup);
        *y = *x + self.p_twice - t;
        *x += t;

        debug_assert!(*x < 4 * self.p.p);
        debug_assert!(*y < 4 * self.p.p);
    }

    /// Inverse NTT butterfly.
    fn inv_butterfly(&self, x: &mut u64, y: &mut u64, z: u64, z_shoup: u64) {
        debug_assert!(*x < self.p_twice);
        debug_assert!(*y < self.p_twice);
        debug_assert!(z < self.p.p);
        debug_assert_eq!(self.p.shoup(z), z_shoup);

        let t = *x;
        *x = Modulus::reduce1(*y + t, self.p_twice);
        *y = self.p.lazy_mul_shoup(self.p_twice + t - *y, z, z_shoup);

        debug_assert!(*x < self.p_twice);
        debug_assert!(*y < self.p_twice);
    }

    /// Inverse NTT butterfly in variable time
    unsafe fn inv_butterfly_vt(&self, x: &mut u64, y: &mut u64, z: u64, z_shoup: u64) {
        debug_assert!(*x < self.p_twice);
        debug_assert!(*y < self.p_twice);
        debug_assert!(z < self.p.p);
        debug_assert_eq!(self.p.shoup(z), z_shoup);

        let t = *x;
        *x = Modulus::reduce1_vt(*y + t, self.p_twice);
        *y = self.p.lazy_mul_shoup(self.p_twice + t - *y, z, z_shoup);

        debug_assert!(*x < self.p_twice);
        debug_assert!(*y < self.p_twice);
    }

    /// Returns a 2n-th primitive root modulo p.
    ///
    /// Aborts if p is not prime or n is not a power of 2 that is >= 8.
    fn primitive_root(n: usize, p: &Modulus) -> u64 {
        debug_assert!(super::supports_ntt(p.p, n));

        let lambda = (p.p - 1) / (2 * n as u64);

        let mut rng: ChaCha8Rng = SeedableRng::seed_from_u64(0);
        for _ in 0..100 {
            let mut root = rng.gen_range(0..p.p);
            root = p.pow(root, lambda);
            if Self::is_primitive_root(root, 2 * n, p) {
                return root;
            }
        }

        debug_assert!(false, "Couldn't find primitive root");
        0
    }

    /// Returns whether a is a n-th primitive root of unity.
    ///
    /// Aborts if a >= p in debug mode.
    fn is_primitive_root(a: u64, n: usize, p: &Modulus) -> bool {
        debug_assert!(a < p.p);
        debug_assert!(super::supports_ntt(p.p, n >> 1)); // TODO: This is not exactly the right condition here.

        // A primitive root of unity is such that x^n = 1 mod p, and x^(n/p) != 1 mod p
        // for all prime p dividing n.
        (p.pow(a, n as u64) == 1) && (p.pow(a, (n / 2) as u64) != 1)
    }
}
