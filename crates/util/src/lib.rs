#![crate_name = "util"]
#![crate_type = "lib"]
#![warn(missing_docs, unused_imports)]
#![feature(is_some_with)]
#![feature(int_roundings)]
#![feature(test)]

//! Utilities for the fhe.rs library.

use num_bigint::{prime::probably_prime, BigUint, ModInverse};
use num_traits::cast::ToPrimitive;
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

/// Sample a vector of independent centered binomial distributions of a given
/// variance. Currently, supports variance up to 16.
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

/// Transcodes a vector of u64 of nbits numbers into a vector of bytes
/// TODO: To test
pub fn transcode_forward(a: &[u64], nbits: usize) -> Vec<u8> {
	assert!(nbits <= 64);
	let mask = (u64::MAX >> (64 - nbits)) as u128;
	let nbytes = (a.len() * nbits).div_ceil(8);
	let mut out = Vec::with_capacity(nbytes);

	let mut current_index = 0;
	let mut current_value = 0u128;
	let mut current_value_nbits = 0;
	while current_index < a.len() {
		if current_value_nbits < 8 {
			debug_assert!(64 - a[current_index].leading_zeros() <= nbits as u32);
			current_value |= ((a[current_index] as u128) & mask) << current_value_nbits;
			current_value_nbits += nbits;
			current_index += 1;
		}
		while current_value_nbits >= 8 {
			out.push(current_value as u8);
			current_value >>= 8;
			current_value_nbits -= 8;
		}
	}
	if current_value_nbits > 0 {
		assert!(current_value_nbits < 8);
		assert_eq!(out.len(), nbytes - 1);
		out.push(current_value as u8)
	} else {
		assert_eq!(out.len(), nbytes);
		assert_eq!(current_value, 0);
	}
	out
}

/// Transcodes a vector of u8 into a vector of u64 of nbits numbers
/// TODO: To test
pub fn transcode_backward(b: &[u8], nbits: usize) -> Vec<u64> {
	assert!(nbits <= 64);
	let mask = (u64::MAX >> (64 - nbits)) as u128;

	let nelements = (b.len() * 8).div_ceil(nbits);
	let mut out = Vec::with_capacity(nelements);

	let mut current_value = 0u128;
	let mut current_value_nbits = 0;
	let mut current_index = 0;
	while current_index < b.len() {
		if current_value_nbits < nbits {
			current_value |= (b[current_index] as u128) << current_value_nbits;
			current_value_nbits += 8;
			current_index += 1;
		}
		while current_value_nbits >= nbits {
			out.push((current_value & mask) as u64);
			current_value >>= nbits;
			current_value_nbits -= nbits;
		}
	}
	if current_value_nbits > 0 {
		assert_eq!(out.len(), nelements - 1);
		out.push(current_value as u64);
	} else {
		assert_eq!(out.len(), nelements);
		assert_eq!(current_value, 0);
	}
	out
}

/// Computes the inverse
/// TODO: To test
pub fn inverse(a: u64, p: u64) -> Option<u64> {
	let p = BigUint::from(p);
	let a = BigUint::from(a);
	a.mod_inverse(p)?.to_u64()
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
			// but for now, we will just check that the rounded value is equal to the
			// variance.
			assert!(w.var.round() == (var as f64));
		}
	}
}
