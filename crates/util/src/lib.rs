#![crate_name = "util"]
#![crate_type = "lib"]
#![warn(missing_docs, unused_imports)]
#![feature(is_some_with)]
#![feature(test)]

//! Utilities for the fhe.rs library.

use num_bigint::{prime::probably_prime, BigUint};
use rand::{thread_rng, RngCore};
use std::panic::UnwindSafe;

/// Define catch_unwind to silence the panic in unit tests.
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

/// Returns whether the modulus p is prime; this function is 100% accurate.
pub fn is_prime(p: u64) -> bool {
	probably_prime(&BigUint::from(p), 0)
}

/// Sample a vector of independent centered binomial distributions of a given variance.
/// Currently, supports variance up to 16.
pub fn sample_vec_cbd(vector_size: usize, variance: usize) -> Result<Vec<i64>, &'static str> {
	if !(1..=16).contains(&variance) {
		return Err("The variance should be between 1 and 16");
	}

	let mut out = Vec::with_capacity(vector_size);

	let number_bits = 4 * variance;
	let mask_add = ((u64::MAX >> (64 - number_bits)) >> (2 * variance)) as u128;
	let mask_sub = (mask_add << (2 * variance)) as u128;

	let mut current_pool = 0u128;
	let mut current_pool_nbits = 0;
	let mut rng = thread_rng();

	for _ in 0..vector_size {
		if current_pool_nbits < number_bits {
			current_pool |= (rng.next_u64() as u128) << current_pool_nbits;
			current_pool_nbits += 64;
		}
		debug_assert!(current_pool_nbits >= number_bits);
		out.push(
			((current_pool & mask_add).count_ones() as i64)
				- ((current_pool & mask_sub).count_ones() as i64),
		);
		current_pool >>= number_bits;
		current_pool_nbits -= number_bits;
	}

	Ok(out)
}

#[cfg(test)]
mod tests {
	extern crate test;
	use super::{is_prime, sample_vec_cbd};
	use itertools::Itertools;

	#[test]
	fn test_prime() {
		assert!(is_prime(2));
		assert!(is_prime(3));
		assert!(is_prime(5));
		assert!(is_prime(7));
		assert!(is_prime(4611686018326724609));

		assert!(!is_prime(0));
		assert!(!is_prime(1));
		assert!(!is_prime(4));
		assert!(!is_prime(6));
		assert!(!is_prime(8));
		assert!(!is_prime(9));
		assert!(!is_prime(4611686018326724607));
	}

	#[test]
	fn test_sample_cbd() {
		assert!(sample_vec_cbd(10, 0).is_err());
		assert!(sample_vec_cbd(10, 17).is_err());

		for var in 1..=16 {
			for size in 0..=100 {
				let v = sample_vec_cbd(size, var).unwrap();
				assert_eq!(v.len(), size);
			}

			// Verifies that the min, max are in absolute value smaller than 2 * var
			let v = sample_vec_cbd(100000, var).unwrap();
			let w = test::stats::Summary::new(&v.iter().map(|vi| *vi as f64).collect_vec());
			assert!(w.max <= (2.0 * var as f64));
			assert!(w.min >= (-2.0 * var as f64));

			// Verifies that the variance is correct. We could probably refine the bound
			// but for now, we will just check that the rounded value is equal to the variance.
			assert!(w.var.round() == (var as f64));
		}
	}
}
