//! Number-Theoretic Transform in ZZ_q.

use super::Modulus;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use std::rc::Rc;

/// Number-Theoretic Transform operator.
#[derive(Debug, Clone, PartialEq)]
pub struct NttOperator {
	p: Rc<Modulus>,
	p_twice: u64,
	size: usize,
	omegas: Vec<u64>,
	omegas_shoup: Vec<u64>,
	zetas_inv: Vec<u64>,
	zetas_inv_shoup: Vec<u64>,
	size_inv: u64,
	size_inv_shoup: u64,
}

impl NttOperator {
	/// Create an NTT operator given a modulus for a specific size.
	///
	/// Aborts if the size is not a power of 2 that is >= 8 in debug mode.
	/// Returns None if the modulus does not support the NTT for this specific size.
	pub fn new(p: &Rc<Modulus>, size: usize) -> std::option::Option<Self> {
		if !p.supports_ntt(size) {
			None
		} else {
			let omega = Self::primitive_root(size, p);
			let omega_inv = p.inv(omega).unwrap();

			let mut omegas = vec![];
			let mut zetas_inv = vec![];
			for i in 0..size {
				let j = i.reverse_bits() >> (size.leading_zeros() + 1);
				let w = p.pow(omega, j as u64);
				let z = p.pow(omega_inv, (j + 1) as u64);
				omegas.push(w);
				zetas_inv.push(z);
			}

			let size_inv = p.inv(size as u64).unwrap();
			let omegas_shoup = p.shoup_vec(&omegas);
			let zetas_inv_shoup = p.shoup_vec(&zetas_inv);
			Some(Self {
				p: p.clone(),
				p_twice: p.p * 2,
				size,
				omegas,
				omegas_shoup,
				zetas_inv,
				zetas_inv_shoup,
				size_inv,
				size_inv_shoup: p.shoup(size_inv),
			})
		}
	}

	/// Compute the forward NTT in place.
	///
	/// Aborts if a is not of the size handled by the operator.
	pub fn forward(&self, a: &mut [u64]) {
		debug_assert!(a.len() == self.size);

		let n = self.size;
		let a_ptr = a.as_mut_ptr();

		let mut l = n >> 1;
		let mut m = 1;
		let mut k = 1;
		while l > 0 {
			for i in 0..m {
				unsafe {
					let omega = *self.omegas.get_unchecked(k);
					let omega_shoup = *self.omegas_shoup.get_unchecked(k);
					k += 1;

					let s = 2 * i * l;
					match l {
						1 => {
							let uj = a_ptr.add(s);
							let ujl = a_ptr.add(s + l);
							self.butterfly(&mut *uj, &mut *ujl, omega, omega_shoup);
							*uj = self.reduce3(*uj);
							*ujl = self.reduce3(*ujl);
						}
						_ => {
							for j in s..(s + l) {
								self.butterfly(
									&mut *a_ptr.add(j),
									&mut *a_ptr.add(j + l),
									omega,
									omega_shoup,
								);
							}
						}
					}
				}
			}
			l >>= 1;
			m <<= 1;
		}
	}

	/// Compute the backward NTT in place.
	///
	/// Aborts if a is not of the size handled by the operator.
	pub fn backward(&self, a: &mut [u64]) {
		debug_assert!(a.len() == self.size);

		let a_ptr = a.as_mut_ptr();

		let mut k = 0;
		let mut m = self.size >> 1;
		let mut l = 1;
		while m > 0 {
			for i in 0..m {
				let s = 2 * i * l;
				unsafe {
					let zeta_inv = *self.zetas_inv.get_unchecked(k);
					let zeta_inv_shoup = *self.zetas_inv_shoup.get_unchecked(k);
					k += 1;
					match l {
						1 => {
							self.inv_butterfly(
								&mut *a_ptr.add(s),
								&mut *a_ptr.add(s + l),
								zeta_inv,
								zeta_inv_shoup,
							);
						}
						2 => {
							self.inv_butterfly(
								&mut *a_ptr.add(s),
								&mut *a_ptr.add(s + 2),
								zeta_inv,
								zeta_inv_shoup,
							);
							self.inv_butterfly(
								&mut *a_ptr.add(s + 1),
								&mut *a_ptr.add(s + 3),
								zeta_inv,
								zeta_inv_shoup,
							);
						}
						_ => {
							for j in s..(s + l) {
								self.inv_butterfly(
									&mut *a_ptr.add(j),
									&mut *a_ptr.add(j + l),
									zeta_inv,
									zeta_inv_shoup,
								);
							}
						}
					}
				}
			}
			l <<= 1;
			m >>= 1;
		}

		a.iter_mut()
			.for_each(|ai| *ai = self.p.mul_shoup(*ai, self.size_inv, self.size_inv_shoup));
	}

	/// Reduce a modulo p.
	///
	/// Aborts if a >= 4 * p.
	fn reduce3(&self, a: u64) -> u64 {
		debug_assert!(a < 4 * self.p.p);

		let y = Modulus::reduce1(a, 2 * self.p.p);
		Modulus::reduce1(y, self.p.p)
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

	/// Inverse NTT butterfly.
	fn inv_butterfly(&self, x: &mut u64, y: &mut u64, z: u64, z_shoup: u64) {
		debug_assert!(*x < 2 * self.p.p);
		debug_assert!(*y < 2 * self.p.p);
		debug_assert!(z < self.p.p);
		debug_assert_eq!(self.p.shoup(z), z_shoup);

		let t = *x;
		*x = Modulus::reduce1(*y + t, self.p_twice);
		*y = self.p.lazy_mul_shoup(self.p_twice + t - *y, z, z_shoup);

		debug_assert!(*x < 2 * self.p.p);
		debug_assert!(*y < 2 * self.p.p);
	}

	/// Returns a 2n-th primitive root modulo p.
	///
	/// Aborts if p is not prime or n is not a power of 2 that is >= 8.
	fn primitive_root(n: usize, p: &Rc<Modulus>) -> u64 {
		debug_assert!(p.supports_ntt(n));

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
	fn is_primitive_root(a: u64, n: usize, p: &Rc<Modulus>) -> bool {
		debug_assert!(a < p.p);
		debug_assert!(p.supports_ntt(n));

		// A primitive root of unity is such that x^n = 1 mod p, and x^(n/p) != 1 mod p
		// for all prime p dividing n.
		(p.pow(a, n as u64) == 1) && (p.pow(a, (n / 2) as u64) != 1)
	}
}

#[cfg(test)]
mod tests {
	use super::NttOperator;
	use crate::zq::Modulus;
	use rand::RngCore;
	use std::rc::Rc;

	#[test]
	fn test_constructor() {
		for size in [8, 1024] {
			for p in [1153, 4611686018326724609] {
				let q = Modulus::new(p).unwrap();
				let supports_ntt = q.supports_ntt(size);

				let op = NttOperator::new(&Rc::new(q), size);

				if supports_ntt {
					assert!(op.is_some());
				} else {
					assert!(op.is_none());
				}
			}
		}
	}

	fn random_vector(size: usize, p: u64) -> Vec<u64> {
		let mut rng = rand::thread_rng();
		let mut v = vec![];
		for _ in 0..size {
			v.push(rng.next_u64() % p)
		}
		v
	}

	#[test]
	fn test_bijection() {
		let ntests = 100;

		for size in [8, 1024] {
			for p in [1153, 4611686018326724609] {
				let q = Modulus::new(p).unwrap();

				if q.supports_ntt(size) {
					let op = NttOperator::new(&Rc::new(q), size).unwrap();

					for _ in 0..ntests {
						let mut a = random_vector(size, p);
						let a_clone = a.clone();

						op.forward(&mut a);
						assert_ne!(a, a_clone);

						op.backward(&mut a);
						assert_eq!(a, a_clone);
					}
				}
			}
		}
	}
}
