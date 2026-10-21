// Expect indexing in this performance-critical NTT implementation.
// The unsafe blocks already indicate this is performance-sensitive code.
#![expect(
    clippy::indexing_slicing,
    reason = "performance or example code relies on validated indices"
)]

use crate::zq::Modulus;
use itertools::Itertools;
use rand::{RngExt, SeedableRng};
use rand_chacha::ChaCha8Rng;
use std::iter::successors;

mod butterflies;

vstd::prelude::verus! {

/// Number-Theoretic Transform operator.
// The fields are visible to Verus; construction and derived cloning are not
// verified yet. Butterfly contracts explicitly require the p_twice relation.
#[verifier::allow(autoderive_clone_without_spec)]
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

}

impl NttOperator {
    /// Number of coefficients transformed by this operator.
    #[must_use]
    #[cfg(not(feature = "tfhe-ntt"))]
    pub const fn size(&self) -> usize {
        self.size
    }

    /// Prime modulus of the transform.
    #[must_use]
    #[cfg(not(feature = "tfhe-ntt"))]
    pub const fn modulus(&self) -> u64 {
        self.p.as_u64()
    }

    /// Create an NTT operator given a modulus for a specific size.
    ///
    /// Returns None if the size is not a power of two at least 8, or the
    /// modulus does not support the NTT for this size.
    #[must_use]
    pub fn new(p: &Modulus, size: usize) -> Option<Self> {
        if !super::supports_ntt(p.as_u64(), size) {
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

            let (omegas, zetas_inv): (Vec<u64>, Vec<u64>) = (0..size)
                .map(|i| {
                    let j = i.reverse_bits() >> (size.leading_zeros() + 1);
                    (powers[j], powers_inv[j])
                })
                .unzip();

            let omegas_shoup = p.shoup_vec(&omegas);
            let zetas_inv_shoup = p.shoup_vec(&zetas_inv);

            Some(Self {
                p: p.clone(),
                p_twice: p.as_u64() * 2,
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
        assert_eq!(a.len(), self.size, "NTT input length mismatch");

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
        assert_eq!(a.len(), self.size, "NTT input length mismatch");

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
        unsafe { self.forward_vt_impl::<false>(a_ptr) };
    }

    // CANONICAL is selected at compile time; the lazy path retains its range.
    unsafe fn forward_vt_impl<const CANONICAL: bool>(&self, a_ptr: *mut u64) {
        let a = unsafe { std::slice::from_raw_parts_mut(a_ptr, self.size) };

        let mut l = self.size >> 1;
        let mut m = 1;
        let mut k = 1;
        while l > 1 {
            for i in 0..m {
                let omega = unsafe { *self.omegas.get_unchecked(k) };
                let omega_shoup = unsafe { *self.omegas_shoup.get_unchecked(k) };
                k += 1;
                let s = 2 * i * l;
                for j in s..(s + l) {
                    // SAFETY: j and j + l are distinct and in bounds because
                    // j + l < s + 2 * l <= 2 * m * l = size.
                    let [x, y] = unsafe { a.get_disjoint_unchecked_mut([j, j + l]) };
                    unsafe { self.butterfly_vt(x, y, omega, omega_shoup) };
                }
            }
            l >>= 1;
            m <<= 1;
        }

        // The final stage has adjacent outputs and one twiddle per pair.
        // Butterfly outputs are < 4p < 2^64, so the same two reductions as
        // the former standalone pass yield canonical values without overflow.
        for (([x, y], &omega), &omega_shoup) in a
            .as_chunks_mut::<2>()
            .0
            .iter_mut()
            .zip(self.omegas[k..].iter())
            .zip(self.omegas_shoup[k..].iter())
        {
            unsafe { self.butterfly_vt(x, y, omega, omega_shoup) };
            if CANONICAL {
                *x = unsafe { self.reduce3_vt(*x) };
                *y = unsafe { self.reduce3_vt(*y) };
            }
        }
    }

    /// Compute the forward NTT in place in variable time.
    ///
    /// # Safety
    /// This function assumes that a_ptr points to at least `size` elements.
    /// This function is not constant time and its timing may reveal information
    /// about the value being reduced.
    pub unsafe fn forward_vt(&self, a_ptr: *mut u64) {
        unsafe { self.forward_vt_impl::<true>(a_ptr) };
    }

    /// Compute the backward NTT in place in variable time.
    ///
    /// # Safety
    /// This function assumes that a_ptr points to at least `size` elements.
    /// This function is not constant time and its timing may reveal information
    /// about the value being reduced.
    pub unsafe fn backward_vt(&self, a_ptr: *mut u64) {
        let a = unsafe { std::slice::from_raw_parts_mut(a_ptr, self.size) };

        let mut k = 0;
        let mut m = self.size >> 1;
        let mut l = 1;
        while m > 0 {
            for i in 0..m {
                let s = 2 * i * l;
                let zeta_inv = unsafe { *self.zetas_inv.get_unchecked(k) };
                let zeta_inv_shoup = unsafe { *self.zetas_inv_shoup.get_unchecked(k) };
                k += 1;
                match l {
                    1 => {
                        // SAFETY: s and s + l are distinct (l > 0) and in-bounds
                        // (s + l < 2 * m * l = size)
                        let [x, y] = unsafe { a.get_disjoint_unchecked_mut([s, s + l]) };
                        unsafe { self.inv_butterfly_vt(x, y, zeta_inv, zeta_inv_shoup) };
                    }
                    _ => {
                        for j in s..(s + l) {
                            // SAFETY: j and j + l are distinct (l > 0) and in-bounds
                            // (j + l < s + 2 * l <= 2 * m * l = size)
                            let [x, y] = unsafe { a.get_disjoint_unchecked_mut([j, j + l]) };
                            unsafe { self.inv_butterfly_vt(x, y, zeta_inv, zeta_inv_shoup) };
                        }
                    }
                }
            }
            l <<= 1;
            m >>= 1;
        }

        for ai in a.iter_mut() {
            *ai = self.p.mul_shoup(*ai, self.size_inv, self.size_inv_shoup);
        }
    }

    /// Reduce a modulo p.
    ///
    /// Aborts if a >= 4 * p.
    const fn reduce3(&self, a: u64) -> u64 {
        debug_assert!(a < 4 * self.p.as_u64());

        let r = crate::zq::verified::reduce3(a, self.p.as_u64(), self.p_twice);
        debug_assert!(r < self.p.as_u64() && r == a % self.p.as_u64());
        r
    }

    /// Reduce a modulo p in variable time.
    ///
    /// Aborts if a >= 4 * p.
    const unsafe fn reduce3_vt(&self, a: u64) -> u64 {
        debug_assert!(a < 4 * self.p.as_u64());

        let r = unsafe { crate::zq::verified::reduce3_vt(a, self.p.as_u64(), self.p_twice) };
        debug_assert!(r < self.p.as_u64() && r == a % self.p.as_u64());
        r
    }

    /// Returns a 2n-th primitive root modulo p.
    ///
    /// Aborts if p is not prime or n is not a power of 2 that is >= 8.
    fn primitive_root(n: usize, p: &Modulus) -> u64 {
        debug_assert!(super::supports_ntt(p.as_u64(), n));

        let lambda = (p.as_u64() - 1) / (2 * n as u64);

        let mut rng: ChaCha8Rng = SeedableRng::seed_from_u64(0);
        for _ in 0..100 {
            let mut root = rng.random_range(0..p.as_u64());
            root = p.pow(root, lambda);
            if Self::is_primitive_root(root, 2 * n, p) {
                return root;
            }
        }

        debug_assert!(false, "Couldn't find primitive root");
        0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Butterfly arithmetic needs no roots of unity or prime modulus. Empty
    // transform tables let these tests cover every supported modulus shape.
    fn butterfly_operator(p: u64) -> NttOperator {
        NttOperator {
            p: Modulus::new(p).unwrap(),
            p_twice: 2 * p,
            size: 0,
            omegas: Box::default(),
            omegas_shoup: Box::default(),
            zetas_inv: Box::default(),
            zetas_inv_shoup: Box::default(),
            size_inv: 0,
            size_inv_shoup: 0,
        }
    }

    fn check_butterfly(op: &NttOperator, a: u64, b: u64, w: u64, inverse: bool) {
        let p = op.p.as_u64();
        let shoup = (((w as u128) << 64) / p as u128) as u64;
        let (mut x, mut y) = (a, b);
        let (mut x_vt, mut y_vt) = (a, b);
        let (expected_x, expected_y, bound) = if inverse {
            op.inv_butterfly(&mut x, &mut y, w, shoup);
            unsafe { op.inv_butterfly_vt(&mut x_vt, &mut y_vt, w, shoup) };
            (
                a as i128 + b as i128,
                (a as i128 - b as i128) * w as i128,
                2 * p,
            )
        } else {
            op.butterfly(&mut x, &mut y, w, shoup);
            unsafe { op.butterfly_vt(&mut x_vt, &mut y_vt, w, shoup) };
            (
                a as i128 + b as i128 * w as i128,
                a as i128 - b as i128 * w as i128,
                4 * p,
            )
        };
        assert!(x < bound && y < bound);
        assert_eq!(x % p, expected_x.rem_euclid(p as i128) as u64);
        assert_eq!(y % p, expected_y.rem_euclid(p as i128) as u64);
        assert_eq!((x, y), (x_vt, y_vt));
    }

    #[test]
    fn butterflies_exhaust_small_moduli() {
        for p in 2..=13 {
            let op = butterfly_operator(p);
            for a in 0..4 * p {
                for b in 0..4 * p {
                    for w in 0..p {
                        check_butterfly(&op, a, b, w, false);
                        if a < 2 * p && b < 2 * p {
                            check_butterfly(&op, a, b, w, true);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn butterflies_cover_lazy_boundaries_and_full_width_moduli() {
        let limit = crate::zq::verified::MODULUS_LIMIT;
        for p in [2, 3, 4, 17, 1 << 31, (1 << 61) + 1, limit - 2, limit - 1] {
            let op = butterfly_operator(p);
            let operands = [
                0,
                1,
                p - 1,
                p,
                p + 1,
                2 * p - 1,
                2 * p,
                2 * p + 1,
                3 * p - 1,
                3 * p,
                3 * p + 1,
                4 * p - 1,
            ];
            for a in operands {
                for b in operands {
                    for w in [0, 1, p / 2, p - 1] {
                        check_butterfly(&op, a, b, w, false);
                        if a < 2 * p && b < 2 * p {
                            check_butterfly(&op, a, b, w, true);
                        }
                    }
                }
            }
        }

        let mut rng = ChaCha8Rng::seed_from_u64(0);
        for _ in 0..1000 {
            let p = rng.random_range(2..limit);
            let op = butterfly_operator(p);
            let w = rng.random_range(0..p);
            for inverse in [false, true] {
                let bound = if inverse { 2 * p } else { 4 * p };
                let a = rng.random_range(0..bound);
                let b = rng.random_range(0..bound);
                check_butterfly(&op, a, b, w, inverse);
            }
        }
    }

    #[test]
    #[cfg(debug_assertions)]
    fn butterflies_reject_invalid_operands_and_shoup_words() {
        type Butterfly = unsafe fn(&NttOperator, &mut u64, &mut u64, u64, u64);
        let p = 17;
        let op = butterfly_operator(p);
        for (butterfly, bound) in [
            (NttOperator::butterfly as Butterfly, 4 * p),
            (NttOperator::butterfly_vt as Butterfly, 4 * p),
            (NttOperator::inv_butterfly as Butterfly, 2 * p),
            (NttOperator::inv_butterfly_vt as Butterfly, 2 * p),
        ] {
            for (a, b, w, shoup) in [
                (bound, 0, 0, 0),
                (0, bound, 0, 0),
                (0, 0, p, 0),
                (0, 0, 1, 0),
            ] {
                assert!(
                    std::panic::catch_unwind(|| {
                        let (mut x, mut y) = (a, b);
                        unsafe { butterfly(&op, &mut x, &mut y, w, shoup) };
                    })
                    .is_err()
                );
            }
        }
    }

    #[test]
    fn reduce3_covers_small_inputs_and_boundaries() {
        for p in [17, 97, 65537, 2013265921, 4611686018326724609] {
            let op = NttOperator::new(&Modulus::new(p).unwrap(), 8).unwrap();
            let mut inputs = vec![
                0,
                1,
                p - 1,
                p,
                p + 1,
                2 * p - 1,
                2 * p,
                2 * p + 1,
                3 * p - 1,
                3 * p,
                3 * p + 1,
                4 * p - 1,
            ];
            if p < 100 {
                inputs.extend(0..4 * p);
            }
            for a in inputs {
                assert_eq!(op.reduce3(a), a % p);
                assert_eq!(unsafe { op.reduce3_vt(a) }, a % p);
            }
            #[cfg(debug_assertions)]
            {
                assert!(std::panic::catch_unwind(|| op.reduce3(4 * p)).is_err());
                assert!(std::panic::catch_unwind(|| unsafe { op.reduce3_vt(4 * p) }).is_err());
            }
        }
    }

    // Test the native implementation even when the TFHE backend is enabled.
    #[test]
    fn final_pass_ranges_and_round_trips() {
        let mut rng = ChaCha8Rng::seed_from_u64(10);
        for n in [8, 16, 32, 128, 1024, 8192, 16384] {
            for modulus in [65537, 2013265921, 4611686018326724609] {
                let p = Modulus::new(modulus).unwrap();
                let op = NttOperator::new(&p, n).unwrap();
                let boundaries = [
                    0,
                    1,
                    modulus - 1,
                    modulus,
                    2 * modulus - 1,
                    2 * modulus,
                    4 * modulus - 1,
                ];
                let inputs = [
                    vec![0; n],
                    vec![modulus - 1; n],
                    vec![4 * modulus - 1; n],
                    (0..n).map(|i| boundaries[i % boundaries.len()]).collect(),
                    p.random_vec(n, &mut rng),
                ];
                for input in inputs {
                    let canonical: Vec<_> = input.iter().map(|x| x % modulus).collect();
                    let mut expected = canonical.clone();
                    op.forward(&mut expected);
                    let mut ct = input.clone();
                    let mut vt = input.clone();
                    let mut lazy = input;
                    op.forward(&mut ct);
                    unsafe {
                        op.forward_vt(vt.as_mut_ptr());
                        op.forward_vt_lazy(lazy.as_mut_ptr());
                    }
                    assert_eq!(ct, expected);
                    assert_eq!(vt, expected);
                    assert!(vt.iter().all(|&x| x < modulus));
                    assert!(lazy.iter().all(|&x| x < 4 * modulus));
                    assert!(lazy.iter().zip(&expected).all(|(&x, &y)| x % modulus == y));
                    // Inverse inputs can be lazy in [0, 2p), but not [0, 4p).
                    for x in &mut vt {
                        *x += modulus;
                    }
                    op.backward(&mut ct);
                    unsafe {
                        op.backward_vt(vt.as_mut_ptr());
                    }
                    assert_eq!(ct, canonical);
                    assert_eq!(vt, canonical);
                }
                let mut inverse = vec![2 * modulus - 1; n];
                let mut inverse_vt = inverse.clone();
                op.backward(&mut inverse);
                unsafe {
                    op.backward_vt(inverse_vt.as_mut_ptr());
                }
                assert_eq!(inverse, inverse_vt);
                assert!(inverse.iter().all(|&x| x < modulus));
                op.forward(&mut inverse);
                assert_eq!(inverse, vec![modulus - 1; n]);
            }
        }
    }

    #[test]
    fn final_pass_matches_reference_ring_product() {
        let mut rng = ChaCha8Rng::seed_from_u64(11);
        for n in [8, 32, 128, 1024, 8192] {
            for modulus in [65537, 2013265921, 4611686018326724609] {
                let p = Modulus::new(modulus).unwrap();
                let op = NttOperator::new(&p, n).unwrap();
                let a = p.random_vec(n, &mut rng);
                let mut b = p.random_vec(n, &mut rng);
                // Dense schoolbook products at small sizes, sparse at large
                // sizes to check wraparound without quadratic test cost.
                if n > 128 {
                    b.fill(0);
                    b[0] = modulus - 1;
                    b[n / 2] = modulus - 1;
                    b[n - 1] = modulus - 1;
                }
                let mut expected = vec![0u64; n];
                for (j, &y) in b.iter().enumerate().filter(|(_, y)| **y != 0) {
                    for (i, &x) in a.iter().enumerate() {
                        let product = (u128::from(x) * u128::from(y)) % u128::from(modulus);
                        let k = (i + j) % n;
                        let term = if i + j >= n {
                            u128::from(modulus) - product
                        } else {
                            product
                        };
                        expected[k] =
                            ((u128::from(expected[k]) + term) % u128::from(modulus)) as u64;
                    }
                }
                for variable_time in [false, true] {
                    let mut lhs = a.clone();
                    let mut rhs = b.clone();
                    if variable_time {
                        unsafe {
                            op.forward_vt(lhs.as_mut_ptr());
                            op.forward_vt(rhs.as_mut_ptr());
                        }
                    } else {
                        op.forward(&mut lhs);
                        op.forward(&mut rhs);
                    }
                    for (x, y) in lhs.iter_mut().zip(rhs) {
                        *x = p.mul(*x, y);
                    }
                    if variable_time {
                        unsafe {
                            op.backward_vt(lhs.as_mut_ptr());
                        }
                    } else {
                        op.backward(&mut lhs);
                    }
                    assert_eq!(lhs, expected, "degree={n}, modulus={modulus}");
                }
            }
        }
    }
}
