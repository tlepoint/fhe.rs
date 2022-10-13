#![crate_name = "fhe_util"]
#![crate_type = "lib"]
#![warn(missing_docs, unused_imports)]

//! Utilities for the fhe.rs library.

#[cfg(test)]
#[macro_use]
extern crate proptest;

mod u256;
use rand::{CryptoRng, RngCore};
pub use u256::U256;

use num_bigint_dig::{prime::probably_prime, BigUint, ModInverse};
use num_traits::{cast::ToPrimitive, PrimInt};
use std::{mem::size_of, panic::UnwindSafe};

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
/// variance. Returns an error if the variance is strictly larger than 16.
pub fn sample_vec_cbd<R: RngCore + CryptoRng>(
	vector_size: usize,
	variance: usize,
	rng: &mut R,
) -> Result<Vec<i64>, &'static str> {
	if !(1..=16).contains(&variance) {
		return Err("The variance should be between 1 and 16");
	}

	let mut out = Vec::with_capacity(vector_size);

	let number_bits = 4 * variance;
	let mask_add = ((u64::MAX >> (64 - number_bits)) >> (2 * variance)) as u128;
	let mask_sub = (mask_add << (2 * variance)) as u128;

	let mut current_pool = 0u128;
	let mut current_pool_nbits = 0;

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

/// Transcodes a vector of u64 of `nbits`-bit numbers into a vector of bytes.
pub fn transcode_to_bytes(a: &[u64], nbits: usize) -> Vec<u8> {
	assert!(nbits <= 64);
	assert!(nbits > 0);

	let mask = (u64::MAX >> (64 - nbits)) as u128;
	let nbytes = div_ceil(a.len() * nbits, 8);
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

/// Transcodes a vector of u8 into a vector of u64 of `nbits`-bit numbers.
pub fn transcode_from_bytes(b: &[u8], nbits: usize) -> Vec<u64> {
	assert!(nbits <= 64);
	assert!(nbits > 0);
	let mask = (u64::MAX >> (64 - nbits)) as u128;

	let nelements = div_ceil(b.len() * 8, nbits);
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

/// Transcodes a vector of u64 of `input_nbits`-bit numbers into a vector of u64
/// of `output_nbits`-bit numbers.
pub fn transcode_bidirectional(a: &[u64], input_nbits: usize, output_nbits: usize) -> Vec<u64> {
	assert!(input_nbits <= 64);
	assert!(output_nbits <= 64);
	assert!(input_nbits > 0);
	assert!(output_nbits > 0);
	let input_mask = (u64::MAX >> (64 - input_nbits)) as u128;
	let output_mask = (u64::MAX >> (64 - output_nbits)) as u128;
	let output_size = div_ceil(a.len() * input_nbits, output_nbits);
	let mut out = Vec::with_capacity(output_size);

	let mut current_index = 0;
	let mut current_value = 0u128;
	let mut current_value_nbits = 0;
	while current_index < a.len() {
		if current_value_nbits < output_nbits {
			debug_assert!(64 - a[current_index].leading_zeros() <= input_nbits as u32);
			current_value |= ((a[current_index] as u128) & input_mask) << current_value_nbits;
			current_value_nbits += input_nbits;
			current_index += 1;
		}
		while current_value_nbits >= output_nbits {
			out.push((current_value & output_mask) as u64);
			current_value >>= output_nbits;
			current_value_nbits -= output_nbits;
		}
	}
	if current_value_nbits > 0 {
		assert!(current_value_nbits < output_nbits);
		assert_eq!(out.len(), output_size - 1);
		out.push(current_value as u64)
	} else {
		assert_eq!(out.len(), output_size);
		assert_eq!(current_value, 0);
	}
	out
}

/// Computes the modular multiplicative inverse of `a` modulo `p`. Returns
/// `None` if `a` is not invertible modulo `p`.
pub fn inverse(a: u64, p: u64) -> Option<u64> {
	let p = BigUint::from(p);
	let a = BigUint::from(a);
	a.mod_inverse(p)?.to_u64()
}

/// Returns the number of bits b such that 2^b <= value
/// to simulate the `.ilog2()` function from <https://github.com/rust-lang/rust/issues/70887>.
/// Panics when `value` is 0.
pub fn ilog2<T: PrimInt>(value: T) -> usize {
	assert!(value > T::zero());
	// For this, we compute sizeof(T) - 1 - value.leading_zeros(). Indeed, when 2^b
	// <= value < 2^(b+1), then value.leading_zeros() = sizeof(T) - (b + 1).
	size_of::<T>() * 8 - 1 - value.leading_zeros() as usize
}

/// Returns the ceil of a divided by b, to simulate the
/// `.div_ceil()` function from <https://github.com/rust-lang/rust/issues/88581>.
/// Panics when `b` is 0.
pub fn div_ceil<T: PrimInt>(a: T, b: T) -> T {
	assert!(b > T::zero());
	(a + b - T::one()) / b
}

/// Compute the sample variance of a list of values.
/// Panics if the length of value is < 2.
pub fn variance<T: PrimInt>(values: &[T]) -> f64 {
	assert!(values.len() > 1);
	let mean = values.iter().fold(0f64, |acc, i| acc + i.to_f64().unwrap()) / (values.len() as f64);
	values.iter().fold(0f64, |acc, i| {
		acc + (i.to_f64().unwrap() - mean) * (i.to_f64().unwrap() - mean)
	}) / ((values.len() as f64) - 1.0)
}

#[cfg(test)]
mod tests {
	use itertools::Itertools;
	use rand::{thread_rng, RngCore};

	use crate::{div_ceil, ilog2, variance};

	use super::{
		inverse, is_prime, sample_vec_cbd, transcode_bidirectional, transcode_from_bytes,
		transcode_to_bytes,
	};

	#[test]
	fn prime() {
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
	fn ilog2_is_correct() {
		assert_eq!(ilog2(1), 0);
		assert_eq!(ilog2(2), 1);
		assert_eq!(ilog2(3), 1);
		assert_eq!(ilog2(4), 2);
		for i in 2..=110 {
			assert_eq!(ilog2(1u128 << i), i);
			assert_eq!(ilog2((1u128 << i) + 1), i);
			assert_eq!(ilog2((1u128 << (i + 1)) - 1), i);
		}
	}

	#[test]
	fn div_ceil_is_correct() {
		for _ in 0..100 {
			let a = (thread_rng().next_u32() >> 1) as usize;
			assert_eq!(div_ceil(a, 1), a);
			assert_eq!(div_ceil(a, 2), (a >> 1) + (a & 1));
			assert_eq!(div_ceil(a, 8), (a + 7) / 8);
		}
	}

	#[test]
	fn sample_cbd() {
		assert!(sample_vec_cbd(10, 0, &mut thread_rng()).is_err());
		assert!(sample_vec_cbd(10, 17, &mut thread_rng()).is_err());

		for var in 1..=16 {
			for size in 0..=100 {
				let v = sample_vec_cbd(size, var, &mut thread_rng()).unwrap();
				assert_eq!(v.len(), size);
			}

			// Verifies that the min, max are in absolute value smaller than 2 * var
			let v = sample_vec_cbd(100000, var, &mut thread_rng()).unwrap();
			assert!(v.iter().map(|vi| vi.abs()).max().unwrap() <= 2 * var as i64);

			// Verifies that the variance is correct. We could probably refine the bound
			// but for now, we will just check that the rounded value is equal to the
			// variance.
			assert!(variance(&v).round() == (var as f64));
		}
	}

	#[test]
	fn transcode_self_consistency() {
		let mut rng = thread_rng();

		for size in 1..=100 {
			let input = (0..size).map(|_| rng.next_u64()).collect_vec();
			for input_nbits in 1..63 {
				let masked_input = input
					.iter()
					.map(|i| (*i) & (u64::MAX >> (64 - input_nbits)))
					.collect_vec();
				let bytes = transcode_to_bytes(&masked_input, input_nbits);
				let bytes_as_u64 = transcode_bidirectional(&masked_input, input_nbits, 8);
				assert_eq!(bytes, bytes_as_u64.iter().map(|e| *e as u8).collect_vec());

				let input_from_bytes = transcode_from_bytes(&bytes, input_nbits);
				assert!(input_from_bytes.len() >= masked_input.len());
				assert_eq!(input_from_bytes[..masked_input.len()], masked_input);

				let input_from_u64 = transcode_bidirectional(&bytes_as_u64, 8, input_nbits);
				assert!(input_from_u64.len() >= masked_input.len());
				assert_eq!(input_from_u64[..masked_input.len()], masked_input);

				for output_nbits in 1..63 {
					let output = transcode_bidirectional(&masked_input, input_nbits, output_nbits);
					let input_from_output =
						transcode_bidirectional(&output, output_nbits, input_nbits);
					assert!(input_from_output.len() >= masked_input.len());
					assert_eq!(input_from_output[..masked_input.len()], masked_input);
				}
			}
		}
	}

	#[test]
	fn inv_kats() {
		// KATs for inversion generated in Sage using the following code.
		/*
		sage: for p in range(2, 1000, 7):
		....:     for a in range(1, 30, 3):
		....:         if gcd(a, p) == 1:
		....:             i = ZZ(a)^(-1) % p
		....:             print("assert_eq!(inverse({}, {}), Some({}));".format(a, p, i))
		....:         else:
		....:             print("assert!(inverse({}, {}).is_none());".format(a, p))
		 */
		assert_eq!(inverse(1, 2), Some(1));
		assert!(inverse(4, 2).is_none());
		assert_eq!(inverse(7, 2), Some(1));
		assert!(inverse(10, 2).is_none());
		assert_eq!(inverse(13, 2), Some(1));
		assert!(inverse(16, 2).is_none());
		assert_eq!(inverse(19, 2), Some(1));
		assert!(inverse(22, 2).is_none());
		assert_eq!(inverse(25, 2), Some(1));
		assert!(inverse(28, 2).is_none());
		assert_eq!(inverse(1, 9), Some(1));
		assert_eq!(inverse(4, 9), Some(7));
		assert_eq!(inverse(7, 9), Some(4));
		assert_eq!(inverse(10, 9), Some(1));
		assert_eq!(inverse(13, 9), Some(7));
		assert_eq!(inverse(16, 9), Some(4));
		assert_eq!(inverse(19, 9), Some(1));
		assert_eq!(inverse(22, 9), Some(7));
		assert_eq!(inverse(25, 9), Some(4));
		assert_eq!(inverse(28, 9), Some(1));
		assert_eq!(inverse(1, 16), Some(1));
		assert!(inverse(4, 16).is_none());
		assert_eq!(inverse(7, 16), Some(7));
		assert!(inverse(10, 16).is_none());
		assert_eq!(inverse(13, 16), Some(5));
		assert!(inverse(16, 16).is_none());
		assert_eq!(inverse(19, 16), Some(11));
		assert!(inverse(22, 16).is_none());
		assert_eq!(inverse(25, 16), Some(9));
		assert!(inverse(28, 16).is_none());
		assert_eq!(inverse(1, 23), Some(1));
		assert_eq!(inverse(4, 23), Some(6));
		assert_eq!(inverse(7, 23), Some(10));
		assert_eq!(inverse(10, 23), Some(7));
		assert_eq!(inverse(13, 23), Some(16));
		assert_eq!(inverse(16, 23), Some(13));
		assert_eq!(inverse(19, 23), Some(17));
		assert_eq!(inverse(22, 23), Some(22));
		assert_eq!(inverse(25, 23), Some(12));
		assert_eq!(inverse(28, 23), Some(14));
		assert_eq!(inverse(1, 30), Some(1));
		assert!(inverse(4, 30).is_none());
		assert_eq!(inverse(7, 30), Some(13));
		assert!(inverse(10, 30).is_none());
		assert_eq!(inverse(13, 30), Some(7));
		assert!(inverse(16, 30).is_none());
		assert_eq!(inverse(19, 30), Some(19));
		assert!(inverse(22, 30).is_none());
		assert!(inverse(25, 30).is_none());
		assert!(inverse(28, 30).is_none());
		assert_eq!(inverse(1, 37), Some(1));
		assert_eq!(inverse(4, 37), Some(28));
		assert_eq!(inverse(7, 37), Some(16));
		assert_eq!(inverse(10, 37), Some(26));
		assert_eq!(inverse(13, 37), Some(20));
		assert_eq!(inverse(16, 37), Some(7));
		assert_eq!(inverse(19, 37), Some(2));
		assert_eq!(inverse(22, 37), Some(32));
		assert_eq!(inverse(25, 37), Some(3));
		assert_eq!(inverse(28, 37), Some(4));
		assert_eq!(inverse(1, 44), Some(1));
		assert!(inverse(4, 44).is_none());
		assert_eq!(inverse(7, 44), Some(19));
		assert!(inverse(10, 44).is_none());
		assert_eq!(inverse(13, 44), Some(17));
		assert!(inverse(16, 44).is_none());
		assert_eq!(inverse(19, 44), Some(7));
		assert!(inverse(22, 44).is_none());
		assert_eq!(inverse(25, 44), Some(37));
		assert!(inverse(28, 44).is_none());
		assert_eq!(inverse(1, 51), Some(1));
		assert_eq!(inverse(4, 51), Some(13));
		assert_eq!(inverse(7, 51), Some(22));
		assert_eq!(inverse(10, 51), Some(46));
		assert_eq!(inverse(13, 51), Some(4));
		assert_eq!(inverse(16, 51), Some(16));
		assert_eq!(inverse(19, 51), Some(43));
		assert_eq!(inverse(22, 51), Some(7));
		assert_eq!(inverse(25, 51), Some(49));
		assert_eq!(inverse(28, 51), Some(31));
		assert_eq!(inverse(1, 58), Some(1));
		assert!(inverse(4, 58).is_none());
		assert_eq!(inverse(7, 58), Some(25));
		assert!(inverse(10, 58).is_none());
		assert_eq!(inverse(13, 58), Some(9));
		assert!(inverse(16, 58).is_none());
		assert_eq!(inverse(19, 58), Some(55));
		assert!(inverse(22, 58).is_none());
		assert_eq!(inverse(25, 58), Some(7));
		assert!(inverse(28, 58).is_none());
		assert_eq!(inverse(1, 65), Some(1));
		assert_eq!(inverse(4, 65), Some(49));
		assert_eq!(inverse(7, 65), Some(28));
		assert!(inverse(10, 65).is_none());
		assert!(inverse(13, 65).is_none());
		assert_eq!(inverse(16, 65), Some(61));
		assert_eq!(inverse(19, 65), Some(24));
		assert_eq!(inverse(22, 65), Some(3));
		assert!(inverse(25, 65).is_none());
		assert_eq!(inverse(28, 65), Some(7));
		assert_eq!(inverse(1, 72), Some(1));
		assert!(inverse(4, 72).is_none());
		assert_eq!(inverse(7, 72), Some(31));
		assert!(inverse(10, 72).is_none());
		assert_eq!(inverse(13, 72), Some(61));
		assert!(inverse(16, 72).is_none());
		assert_eq!(inverse(19, 72), Some(19));
		assert!(inverse(22, 72).is_none());
		assert_eq!(inverse(25, 72), Some(49));
		assert!(inverse(28, 72).is_none());
		assert_eq!(inverse(1, 79), Some(1));
		assert_eq!(inverse(4, 79), Some(20));
		assert_eq!(inverse(7, 79), Some(34));
		assert_eq!(inverse(10, 79), Some(8));
		assert_eq!(inverse(13, 79), Some(73));
		assert_eq!(inverse(16, 79), Some(5));
		assert_eq!(inverse(19, 79), Some(25));
		assert_eq!(inverse(22, 79), Some(18));
		assert_eq!(inverse(25, 79), Some(19));
		assert_eq!(inverse(28, 79), Some(48));
		assert_eq!(inverse(1, 86), Some(1));
		assert!(inverse(4, 86).is_none());
		assert_eq!(inverse(7, 86), Some(37));
		assert!(inverse(10, 86).is_none());
		assert_eq!(inverse(13, 86), Some(53));
		assert!(inverse(16, 86).is_none());
		assert_eq!(inverse(19, 86), Some(77));
		assert!(inverse(22, 86).is_none());
		assert_eq!(inverse(25, 86), Some(31));
		assert!(inverse(28, 86).is_none());
		assert_eq!(inverse(1, 93), Some(1));
		assert_eq!(inverse(4, 93), Some(70));
		assert_eq!(inverse(7, 93), Some(40));
		assert_eq!(inverse(10, 93), Some(28));
		assert_eq!(inverse(13, 93), Some(43));
		assert_eq!(inverse(16, 93), Some(64));
		assert_eq!(inverse(19, 93), Some(49));
		assert_eq!(inverse(22, 93), Some(55));
		assert_eq!(inverse(25, 93), Some(67));
		assert_eq!(inverse(28, 93), Some(10));
		assert_eq!(inverse(1, 100), Some(1));
		assert!(inverse(4, 100).is_none());
		assert_eq!(inverse(7, 100), Some(43));
		assert!(inverse(10, 100).is_none());
		assert_eq!(inverse(13, 100), Some(77));
		assert!(inverse(16, 100).is_none());
		assert_eq!(inverse(19, 100), Some(79));
		assert!(inverse(22, 100).is_none());
		assert!(inverse(25, 100).is_none());
		assert!(inverse(28, 100).is_none());
		assert_eq!(inverse(1, 107), Some(1));
		assert_eq!(inverse(4, 107), Some(27));
		assert_eq!(inverse(7, 107), Some(46));
		assert_eq!(inverse(10, 107), Some(75));
		assert_eq!(inverse(13, 107), Some(33));
		assert_eq!(inverse(16, 107), Some(87));
		assert_eq!(inverse(19, 107), Some(62));
		assert_eq!(inverse(22, 107), Some(73));
		assert_eq!(inverse(25, 107), Some(30));
		assert_eq!(inverse(28, 107), Some(65));
		assert_eq!(inverse(1, 114), Some(1));
		assert!(inverse(4, 114).is_none());
		assert_eq!(inverse(7, 114), Some(49));
		assert!(inverse(10, 114).is_none());
		assert_eq!(inverse(13, 114), Some(79));
		assert!(inverse(16, 114).is_none());
		assert!(inverse(19, 114).is_none());
		assert!(inverse(22, 114).is_none());
		assert_eq!(inverse(25, 114), Some(73));
		assert!(inverse(28, 114).is_none());
		assert_eq!(inverse(1, 121), Some(1));
		assert_eq!(inverse(4, 121), Some(91));
		assert_eq!(inverse(7, 121), Some(52));
		assert_eq!(inverse(10, 121), Some(109));
		assert_eq!(inverse(13, 121), Some(28));
		assert_eq!(inverse(16, 121), Some(53));
		assert_eq!(inverse(19, 121), Some(51));
		assert!(inverse(22, 121).is_none());
		assert_eq!(inverse(25, 121), Some(92));
		assert_eq!(inverse(28, 121), Some(13));
		assert_eq!(inverse(1, 128), Some(1));
		assert!(inverse(4, 128).is_none());
		assert_eq!(inverse(7, 128), Some(55));
		assert!(inverse(10, 128).is_none());
		assert_eq!(inverse(13, 128), Some(69));
		assert!(inverse(16, 128).is_none());
		assert_eq!(inverse(19, 128), Some(27));
		assert!(inverse(22, 128).is_none());
		assert_eq!(inverse(25, 128), Some(41));
		assert!(inverse(28, 128).is_none());
		assert_eq!(inverse(1, 135), Some(1));
		assert_eq!(inverse(4, 135), Some(34));
		assert_eq!(inverse(7, 135), Some(58));
		assert!(inverse(10, 135).is_none());
		assert_eq!(inverse(13, 135), Some(52));
		assert_eq!(inverse(16, 135), Some(76));
		assert_eq!(inverse(19, 135), Some(64));
		assert_eq!(inverse(22, 135), Some(43));
		assert!(inverse(25, 135).is_none());
		assert_eq!(inverse(28, 135), Some(82));
		assert_eq!(inverse(1, 142), Some(1));
		assert!(inverse(4, 142).is_none());
		assert_eq!(inverse(7, 142), Some(61));
		assert!(inverse(10, 142).is_none());
		assert_eq!(inverse(13, 142), Some(11));
		assert!(inverse(16, 142).is_none());
		assert_eq!(inverse(19, 142), Some(15));
		assert!(inverse(22, 142).is_none());
		assert_eq!(inverse(25, 142), Some(125));
		assert!(inverse(28, 142).is_none());
		assert_eq!(inverse(1, 149), Some(1));
		assert_eq!(inverse(4, 149), Some(112));
		assert_eq!(inverse(7, 149), Some(64));
		assert_eq!(inverse(10, 149), Some(15));
		assert_eq!(inverse(13, 149), Some(23));
		assert_eq!(inverse(16, 149), Some(28));
		assert_eq!(inverse(19, 149), Some(102));
		assert_eq!(inverse(22, 149), Some(61));
		assert_eq!(inverse(25, 149), Some(6));
		assert_eq!(inverse(28, 149), Some(16));
		assert_eq!(inverse(1, 156), Some(1));
		assert!(inverse(4, 156).is_none());
		assert_eq!(inverse(7, 156), Some(67));
		assert!(inverse(10, 156).is_none());
		assert!(inverse(13, 156).is_none());
		assert!(inverse(16, 156).is_none());
		assert_eq!(inverse(19, 156), Some(115));
		assert!(inverse(22, 156).is_none());
		assert_eq!(inverse(25, 156), Some(25));
		assert!(inverse(28, 156).is_none());
		assert_eq!(inverse(1, 163), Some(1));
		assert_eq!(inverse(4, 163), Some(41));
		assert_eq!(inverse(7, 163), Some(70));
		assert_eq!(inverse(10, 163), Some(49));
		assert_eq!(inverse(13, 163), Some(138));
		assert_eq!(inverse(16, 163), Some(51));
		assert_eq!(inverse(19, 163), Some(103));
		assert_eq!(inverse(22, 163), Some(126));
		assert_eq!(inverse(25, 163), Some(150));
		assert_eq!(inverse(28, 163), Some(99));
		assert_eq!(inverse(1, 170), Some(1));
		assert!(inverse(4, 170).is_none());
		assert_eq!(inverse(7, 170), Some(73));
		assert!(inverse(10, 170).is_none());
		assert_eq!(inverse(13, 170), Some(157));
		assert!(inverse(16, 170).is_none());
		assert_eq!(inverse(19, 170), Some(9));
		assert!(inverse(22, 170).is_none());
		assert!(inverse(25, 170).is_none());
		assert!(inverse(28, 170).is_none());
		assert_eq!(inverse(1, 177), Some(1));
		assert_eq!(inverse(4, 177), Some(133));
		assert_eq!(inverse(7, 177), Some(76));
		assert_eq!(inverse(10, 177), Some(124));
		assert_eq!(inverse(13, 177), Some(109));
		assert_eq!(inverse(16, 177), Some(166));
		assert_eq!(inverse(19, 177), Some(28));
		assert_eq!(inverse(22, 177), Some(169));
		assert_eq!(inverse(25, 177), Some(85));
		assert_eq!(inverse(28, 177), Some(19));
		assert_eq!(inverse(1, 184), Some(1));
		assert!(inverse(4, 184).is_none());
		assert_eq!(inverse(7, 184), Some(79));
		assert!(inverse(10, 184).is_none());
		assert_eq!(inverse(13, 184), Some(85));
		assert!(inverse(16, 184).is_none());
		assert_eq!(inverse(19, 184), Some(155));
		assert!(inverse(22, 184).is_none());
		assert_eq!(inverse(25, 184), Some(81));
		assert!(inverse(28, 184).is_none());
		assert_eq!(inverse(1, 191), Some(1));
		assert_eq!(inverse(4, 191), Some(48));
		assert_eq!(inverse(7, 191), Some(82));
		assert_eq!(inverse(10, 191), Some(172));
		assert_eq!(inverse(13, 191), Some(147));
		assert_eq!(inverse(16, 191), Some(12));
		assert_eq!(inverse(19, 191), Some(181));
		assert_eq!(inverse(22, 191), Some(165));
		assert_eq!(inverse(25, 191), Some(107));
		assert_eq!(inverse(28, 191), Some(116));
		assert_eq!(inverse(1, 198), Some(1));
		assert!(inverse(4, 198).is_none());
		assert_eq!(inverse(7, 198), Some(85));
		assert!(inverse(10, 198).is_none());
		assert_eq!(inverse(13, 198), Some(61));
		assert!(inverse(16, 198).is_none());
		assert_eq!(inverse(19, 198), Some(73));
		assert!(inverse(22, 198).is_none());
		assert_eq!(inverse(25, 198), Some(103));
		assert!(inverse(28, 198).is_none());
		assert_eq!(inverse(1, 205), Some(1));
		assert_eq!(inverse(4, 205), Some(154));
		assert_eq!(inverse(7, 205), Some(88));
		assert!(inverse(10, 205).is_none());
		assert_eq!(inverse(13, 205), Some(142));
		assert_eq!(inverse(16, 205), Some(141));
		assert_eq!(inverse(19, 205), Some(54));
		assert_eq!(inverse(22, 205), Some(28));
		assert!(inverse(25, 205).is_none());
		assert_eq!(inverse(28, 205), Some(22));
		assert_eq!(inverse(1, 212), Some(1));
		assert!(inverse(4, 212).is_none());
		assert_eq!(inverse(7, 212), Some(91));
		assert!(inverse(10, 212).is_none());
		assert_eq!(inverse(13, 212), Some(49));
		assert!(inverse(16, 212).is_none());
		assert_eq!(inverse(19, 212), Some(67));
		assert!(inverse(22, 212).is_none());
		assert_eq!(inverse(25, 212), Some(17));
		assert!(inverse(28, 212).is_none());
		assert_eq!(inverse(1, 219), Some(1));
		assert_eq!(inverse(4, 219), Some(55));
		assert_eq!(inverse(7, 219), Some(94));
		assert_eq!(inverse(10, 219), Some(22));
		assert_eq!(inverse(13, 219), Some(118));
		assert_eq!(inverse(16, 219), Some(178));
		assert_eq!(inverse(19, 219), Some(196));
		assert_eq!(inverse(22, 219), Some(10));
		assert_eq!(inverse(25, 219), Some(184));
		assert_eq!(inverse(28, 219), Some(133));
		assert_eq!(inverse(1, 226), Some(1));
		assert!(inverse(4, 226).is_none());
		assert_eq!(inverse(7, 226), Some(97));
		assert!(inverse(10, 226).is_none());
		assert_eq!(inverse(13, 226), Some(87));
		assert!(inverse(16, 226).is_none());
		assert_eq!(inverse(19, 226), Some(119));
		assert!(inverse(22, 226).is_none());
		assert_eq!(inverse(25, 226), Some(217));
		assert!(inverse(28, 226).is_none());
		assert_eq!(inverse(1, 233), Some(1));
		assert_eq!(inverse(4, 233), Some(175));
		assert_eq!(inverse(7, 233), Some(100));
		assert_eq!(inverse(10, 233), Some(70));
		assert_eq!(inverse(13, 233), Some(18));
		assert_eq!(inverse(16, 233), Some(102));
		assert_eq!(inverse(19, 233), Some(184));
		assert_eq!(inverse(22, 233), Some(53));
		assert_eq!(inverse(25, 233), Some(28));
		assert_eq!(inverse(28, 233), Some(25));
		assert_eq!(inverse(1, 240), Some(1));
		assert!(inverse(4, 240).is_none());
		assert_eq!(inverse(7, 240), Some(103));
		assert!(inverse(10, 240).is_none());
		assert_eq!(inverse(13, 240), Some(37));
		assert!(inverse(16, 240).is_none());
		assert_eq!(inverse(19, 240), Some(139));
		assert!(inverse(22, 240).is_none());
		assert!(inverse(25, 240).is_none());
		assert!(inverse(28, 240).is_none());
		assert_eq!(inverse(1, 247), Some(1));
		assert_eq!(inverse(4, 247), Some(62));
		assert_eq!(inverse(7, 247), Some(106));
		assert_eq!(inverse(10, 247), Some(173));
		assert!(inverse(13, 247).is_none());
		assert_eq!(inverse(16, 247), Some(139));
		assert!(inverse(19, 247).is_none());
		assert_eq!(inverse(22, 247), Some(146));
		assert_eq!(inverse(25, 247), Some(168));
		assert_eq!(inverse(28, 247), Some(150));
		assert_eq!(inverse(1, 254), Some(1));
		assert!(inverse(4, 254).is_none());
		assert_eq!(inverse(7, 254), Some(109));
		assert!(inverse(10, 254).is_none());
		assert_eq!(inverse(13, 254), Some(215));
		assert!(inverse(16, 254).is_none());
		assert_eq!(inverse(19, 254), Some(107));
		assert!(inverse(22, 254).is_none());
		assert_eq!(inverse(25, 254), Some(61));
		assert!(inverse(28, 254).is_none());
		assert_eq!(inverse(1, 261), Some(1));
		assert_eq!(inverse(4, 261), Some(196));
		assert_eq!(inverse(7, 261), Some(112));
		assert_eq!(inverse(10, 261), Some(235));
		assert_eq!(inverse(13, 261), Some(241));
		assert_eq!(inverse(16, 261), Some(49));
		assert_eq!(inverse(19, 261), Some(55));
		assert_eq!(inverse(22, 261), Some(178));
		assert_eq!(inverse(25, 261), Some(94));
		assert_eq!(inverse(28, 261), Some(28));
		assert_eq!(inverse(1, 268), Some(1));
		assert!(inverse(4, 268).is_none());
		assert_eq!(inverse(7, 268), Some(115));
		assert!(inverse(10, 268).is_none());
		assert_eq!(inverse(13, 268), Some(165));
		assert!(inverse(16, 268).is_none());
		assert_eq!(inverse(19, 268), Some(127));
		assert!(inverse(22, 268).is_none());
		assert_eq!(inverse(25, 268), Some(193));
		assert!(inverse(28, 268).is_none());
		assert_eq!(inverse(1, 275), Some(1));
		assert_eq!(inverse(4, 275), Some(69));
		assert_eq!(inverse(7, 275), Some(118));
		assert!(inverse(10, 275).is_none());
		assert_eq!(inverse(13, 275), Some(127));
		assert_eq!(inverse(16, 275), Some(86));
		assert_eq!(inverse(19, 275), Some(29));
		assert!(inverse(22, 275).is_none());
		assert!(inverse(25, 275).is_none());
		assert_eq!(inverse(28, 275), Some(167));
		assert_eq!(inverse(1, 282), Some(1));
		assert!(inverse(4, 282).is_none());
		assert_eq!(inverse(7, 282), Some(121));
		assert!(inverse(10, 282).is_none());
		assert_eq!(inverse(13, 282), Some(217));
		assert!(inverse(16, 282).is_none());
		assert_eq!(inverse(19, 282), Some(193));
		assert!(inverse(22, 282).is_none());
		assert_eq!(inverse(25, 282), Some(79));
		assert!(inverse(28, 282).is_none());
		assert_eq!(inverse(1, 289), Some(1));
		assert_eq!(inverse(4, 289), Some(217));
		assert_eq!(inverse(7, 289), Some(124));
		assert_eq!(inverse(10, 289), Some(29));
		assert_eq!(inverse(13, 289), Some(89));
		assert_eq!(inverse(16, 289), Some(271));
		assert_eq!(inverse(19, 289), Some(213));
		assert_eq!(inverse(22, 289), Some(92));
		assert_eq!(inverse(25, 289), Some(185));
		assert_eq!(inverse(28, 289), Some(31));
		assert_eq!(inverse(1, 296), Some(1));
		assert!(inverse(4, 296).is_none());
		assert_eq!(inverse(7, 296), Some(127));
		assert!(inverse(10, 296).is_none());
		assert_eq!(inverse(13, 296), Some(205));
		assert!(inverse(16, 296).is_none());
		assert_eq!(inverse(19, 296), Some(187));
		assert!(inverse(22, 296).is_none());
		assert_eq!(inverse(25, 296), Some(225));
		assert!(inverse(28, 296).is_none());
		assert_eq!(inverse(1, 303), Some(1));
		assert_eq!(inverse(4, 303), Some(76));
		assert_eq!(inverse(7, 303), Some(130));
		assert_eq!(inverse(10, 303), Some(91));
		assert_eq!(inverse(13, 303), Some(70));
		assert_eq!(inverse(16, 303), Some(19));
		assert_eq!(inverse(19, 303), Some(16));
		assert_eq!(inverse(22, 303), Some(124));
		assert_eq!(inverse(25, 303), Some(97));
		assert_eq!(inverse(28, 303), Some(184));
		assert_eq!(inverse(1, 310), Some(1));
		assert!(inverse(4, 310).is_none());
		assert_eq!(inverse(7, 310), Some(133));
		assert!(inverse(10, 310).is_none());
		assert_eq!(inverse(13, 310), Some(167));
		assert!(inverse(16, 310).is_none());
		assert_eq!(inverse(19, 310), Some(49));
		assert!(inverse(22, 310).is_none());
		assert!(inverse(25, 310).is_none());
		assert!(inverse(28, 310).is_none());
		assert_eq!(inverse(1, 317), Some(1));
		assert_eq!(inverse(4, 317), Some(238));
		assert_eq!(inverse(7, 317), Some(136));
		assert_eq!(inverse(10, 317), Some(222));
		assert_eq!(inverse(13, 317), Some(122));
		assert_eq!(inverse(16, 317), Some(218));
		assert_eq!(inverse(19, 317), Some(267));
		assert_eq!(inverse(22, 317), Some(245));
		assert_eq!(inverse(25, 317), Some(279));
		assert_eq!(inverse(28, 317), Some(34));
		assert_eq!(inverse(1, 324), Some(1));
		assert!(inverse(4, 324).is_none());
		assert_eq!(inverse(7, 324), Some(139));
		assert!(inverse(10, 324).is_none());
		assert_eq!(inverse(13, 324), Some(25));
		assert!(inverse(16, 324).is_none());
		assert_eq!(inverse(19, 324), Some(307));
		assert!(inverse(22, 324).is_none());
		assert_eq!(inverse(25, 324), Some(13));
		assert!(inverse(28, 324).is_none());
		assert_eq!(inverse(1, 331), Some(1));
		assert_eq!(inverse(4, 331), Some(83));
		assert_eq!(inverse(7, 331), Some(142));
		assert_eq!(inverse(10, 331), Some(298));
		assert_eq!(inverse(13, 331), Some(51));
		assert_eq!(inverse(16, 331), Some(269));
		assert_eq!(inverse(19, 331), Some(122));
		assert_eq!(inverse(22, 331), Some(316));
		assert_eq!(inverse(25, 331), Some(53));
		assert_eq!(inverse(28, 331), Some(201));
		assert_eq!(inverse(1, 338), Some(1));
		assert!(inverse(4, 338).is_none());
		assert_eq!(inverse(7, 338), Some(145));
		assert!(inverse(10, 338).is_none());
		assert!(inverse(13, 338).is_none());
		assert!(inverse(16, 338).is_none());
		assert_eq!(inverse(19, 338), Some(89));
		assert!(inverse(22, 338).is_none());
		assert_eq!(inverse(25, 338), Some(311));
		assert!(inverse(28, 338).is_none());
		assert_eq!(inverse(1, 345), Some(1));
		assert_eq!(inverse(4, 345), Some(259));
		assert_eq!(inverse(7, 345), Some(148));
		assert!(inverse(10, 345).is_none());
		assert_eq!(inverse(13, 345), Some(292));
		assert_eq!(inverse(16, 345), Some(151));
		assert_eq!(inverse(19, 345), Some(109));
		assert_eq!(inverse(22, 345), Some(298));
		assert!(inverse(25, 345).is_none());
		assert_eq!(inverse(28, 345), Some(37));
		assert_eq!(inverse(1, 352), Some(1));
		assert!(inverse(4, 352).is_none());
		assert_eq!(inverse(7, 352), Some(151));
		assert!(inverse(10, 352).is_none());
		assert_eq!(inverse(13, 352), Some(325));
		assert!(inverse(16, 352).is_none());
		assert_eq!(inverse(19, 352), Some(315));
		assert!(inverse(22, 352).is_none());
		assert_eq!(inverse(25, 352), Some(169));
		assert!(inverse(28, 352).is_none());
		assert_eq!(inverse(1, 359), Some(1));
		assert_eq!(inverse(4, 359), Some(90));
		assert_eq!(inverse(7, 359), Some(154));
		assert_eq!(inverse(10, 359), Some(36));
		assert_eq!(inverse(13, 359), Some(221));
		assert_eq!(inverse(16, 359), Some(202));
		assert_eq!(inverse(19, 359), Some(189));
		assert_eq!(inverse(22, 359), Some(49));
		assert_eq!(inverse(25, 359), Some(158));
		assert_eq!(inverse(28, 359), Some(218));
		assert_eq!(inverse(1, 366), Some(1));
		assert!(inverse(4, 366).is_none());
		assert_eq!(inverse(7, 366), Some(157));
		assert!(inverse(10, 366).is_none());
		assert_eq!(inverse(13, 366), Some(169));
		assert!(inverse(16, 366).is_none());
		assert_eq!(inverse(19, 366), Some(289));
		assert!(inverse(22, 366).is_none());
		assert_eq!(inverse(25, 366), Some(205));
		assert!(inverse(28, 366).is_none());
		assert_eq!(inverse(1, 373), Some(1));
		assert_eq!(inverse(4, 373), Some(280));
		assert_eq!(inverse(7, 373), Some(160));
		assert_eq!(inverse(10, 373), Some(112));
		assert_eq!(inverse(13, 373), Some(287));
		assert_eq!(inverse(16, 373), Some(70));
		assert_eq!(inverse(19, 373), Some(216));
		assert_eq!(inverse(22, 373), Some(17));
		assert_eq!(inverse(25, 373), Some(194));
		assert_eq!(inverse(28, 373), Some(40));
		assert_eq!(inverse(1, 380), Some(1));
		assert!(inverse(4, 380).is_none());
		assert_eq!(inverse(7, 380), Some(163));
		assert!(inverse(10, 380).is_none());
		assert_eq!(inverse(13, 380), Some(117));
		assert!(inverse(16, 380).is_none());
		assert!(inverse(19, 380).is_none());
		assert!(inverse(22, 380).is_none());
		assert!(inverse(25, 380).is_none());
		assert!(inverse(28, 380).is_none());
		assert_eq!(inverse(1, 387), Some(1));
		assert_eq!(inverse(4, 387), Some(97));
		assert_eq!(inverse(7, 387), Some(166));
		assert_eq!(inverse(10, 387), Some(271));
		assert_eq!(inverse(13, 387), Some(268));
		assert_eq!(inverse(16, 387), Some(121));
		assert_eq!(inverse(19, 387), Some(163));
		assert_eq!(inverse(22, 387), Some(88));
		assert_eq!(inverse(25, 387), Some(31));
		assert_eq!(inverse(28, 387), Some(235));
		assert_eq!(inverse(1, 394), Some(1));
		assert!(inverse(4, 394).is_none());
		assert_eq!(inverse(7, 394), Some(169));
		assert!(inverse(10, 394).is_none());
		assert_eq!(inverse(13, 394), Some(91));
		assert!(inverse(16, 394).is_none());
		assert_eq!(inverse(19, 394), Some(83));
		assert!(inverse(22, 394).is_none());
		assert_eq!(inverse(25, 394), Some(331));
		assert!(inverse(28, 394).is_none());
		assert_eq!(inverse(1, 401), Some(1));
		assert_eq!(inverse(4, 401), Some(301));
		assert_eq!(inverse(7, 401), Some(172));
		assert_eq!(inverse(10, 401), Some(361));
		assert_eq!(inverse(13, 401), Some(216));
		assert_eq!(inverse(16, 401), Some(376));
		assert_eq!(inverse(19, 401), Some(190));
		assert_eq!(inverse(22, 401), Some(237));
		assert_eq!(inverse(25, 401), Some(385));
		assert_eq!(inverse(28, 401), Some(43));
		assert_eq!(inverse(1, 408), Some(1));
		assert!(inverse(4, 408).is_none());
		assert_eq!(inverse(7, 408), Some(175));
		assert!(inverse(10, 408).is_none());
		assert_eq!(inverse(13, 408), Some(157));
		assert!(inverse(16, 408).is_none());
		assert_eq!(inverse(19, 408), Some(43));
		assert!(inverse(22, 408).is_none());
		assert_eq!(inverse(25, 408), Some(49));
		assert!(inverse(28, 408).is_none());
		assert_eq!(inverse(1, 415), Some(1));
		assert_eq!(inverse(4, 415), Some(104));
		assert_eq!(inverse(7, 415), Some(178));
		assert!(inverse(10, 415).is_none());
		assert_eq!(inverse(13, 415), Some(32));
		assert_eq!(inverse(16, 415), Some(26));
		assert_eq!(inverse(19, 415), Some(284));
		assert_eq!(inverse(22, 415), Some(283));
		assert!(inverse(25, 415).is_none());
		assert_eq!(inverse(28, 415), Some(252));
		assert_eq!(inverse(1, 422), Some(1));
		assert!(inverse(4, 422).is_none());
		assert_eq!(inverse(7, 422), Some(181));
		assert!(inverse(10, 422).is_none());
		assert_eq!(inverse(13, 422), Some(65));
		assert!(inverse(16, 422).is_none());
		assert_eq!(inverse(19, 422), Some(311));
		assert!(inverse(22, 422).is_none());
		assert_eq!(inverse(25, 422), Some(287));
		assert!(inverse(28, 422).is_none());
		assert_eq!(inverse(1, 429), Some(1));
		assert_eq!(inverse(4, 429), Some(322));
		assert_eq!(inverse(7, 429), Some(184));
		assert_eq!(inverse(10, 429), Some(43));
		assert!(inverse(13, 429).is_none());
		assert_eq!(inverse(16, 429), Some(295));
		assert_eq!(inverse(19, 429), Some(271));
		assert!(inverse(22, 429).is_none());
		assert_eq!(inverse(25, 429), Some(103));
		assert_eq!(inverse(28, 429), Some(46));
		assert_eq!(inverse(1, 436), Some(1));
		assert!(inverse(4, 436).is_none());
		assert_eq!(inverse(7, 436), Some(187));
		assert!(inverse(10, 436).is_none());
		assert_eq!(inverse(13, 436), Some(369));
		assert!(inverse(16, 436).is_none());
		assert_eq!(inverse(19, 436), Some(23));
		assert!(inverse(22, 436).is_none());
		assert_eq!(inverse(25, 436), Some(157));
		assert!(inverse(28, 436).is_none());
		assert_eq!(inverse(1, 443), Some(1));
		assert_eq!(inverse(4, 443), Some(111));
		assert_eq!(inverse(7, 443), Some(190));
		assert_eq!(inverse(10, 443), Some(133));
		assert_eq!(inverse(13, 443), Some(409));
		assert_eq!(inverse(16, 443), Some(360));
		assert_eq!(inverse(19, 443), Some(70));
		assert_eq!(inverse(22, 443), Some(141));
		assert_eq!(inverse(25, 443), Some(319));
		assert_eq!(inverse(28, 443), Some(269));
		assert_eq!(inverse(1, 450), Some(1));
		assert!(inverse(4, 450).is_none());
		assert_eq!(inverse(7, 450), Some(193));
		assert!(inverse(10, 450).is_none());
		assert_eq!(inverse(13, 450), Some(277));
		assert!(inverse(16, 450).is_none());
		assert_eq!(inverse(19, 450), Some(379));
		assert!(inverse(22, 450).is_none());
		assert!(inverse(25, 450).is_none());
		assert!(inverse(28, 450).is_none());
		assert_eq!(inverse(1, 457), Some(1));
		assert_eq!(inverse(4, 457), Some(343));
		assert_eq!(inverse(7, 457), Some(196));
		assert_eq!(inverse(10, 457), Some(320));
		assert_eq!(inverse(13, 457), Some(211));
		assert_eq!(inverse(16, 457), Some(200));
		assert_eq!(inverse(19, 457), Some(433));
		assert_eq!(inverse(22, 457), Some(187));
		assert_eq!(inverse(25, 457), Some(128));
		assert_eq!(inverse(28, 457), Some(49));
		assert_eq!(inverse(1, 464), Some(1));
		assert!(inverse(4, 464).is_none());
		assert_eq!(inverse(7, 464), Some(199));
		assert!(inverse(10, 464).is_none());
		assert_eq!(inverse(13, 464), Some(357));
		assert!(inverse(16, 464).is_none());
		assert_eq!(inverse(19, 464), Some(171));
		assert!(inverse(22, 464).is_none());
		assert_eq!(inverse(25, 464), Some(297));
		assert!(inverse(28, 464).is_none());
		assert_eq!(inverse(1, 471), Some(1));
		assert_eq!(inverse(4, 471), Some(118));
		assert_eq!(inverse(7, 471), Some(202));
		assert_eq!(inverse(10, 471), Some(424));
		assert_eq!(inverse(13, 471), Some(145));
		assert_eq!(inverse(16, 471), Some(265));
		assert_eq!(inverse(19, 471), Some(124));
		assert_eq!(inverse(22, 471), Some(364));
		assert_eq!(inverse(25, 471), Some(358));
		assert_eq!(inverse(28, 471), Some(286));
		assert_eq!(inverse(1, 478), Some(1));
		assert!(inverse(4, 478).is_none());
		assert_eq!(inverse(7, 478), Some(205));
		assert!(inverse(10, 478).is_none());
		assert_eq!(inverse(13, 478), Some(331));
		assert!(inverse(16, 478).is_none());
		assert_eq!(inverse(19, 478), Some(151));
		assert!(inverse(22, 478).is_none());
		assert_eq!(inverse(25, 478), Some(153));
		assert!(inverse(28, 478).is_none());
		assert_eq!(inverse(1, 485), Some(1));
		assert_eq!(inverse(4, 485), Some(364));
		assert_eq!(inverse(7, 485), Some(208));
		assert!(inverse(10, 485).is_none());
		assert_eq!(inverse(13, 485), Some(112));
		assert_eq!(inverse(16, 485), Some(91));
		assert_eq!(inverse(19, 485), Some(434));
		assert_eq!(inverse(22, 485), Some(463));
		assert!(inverse(25, 485).is_none());
		assert_eq!(inverse(28, 485), Some(52));
		assert_eq!(inverse(1, 492), Some(1));
		assert!(inverse(4, 492).is_none());
		assert_eq!(inverse(7, 492), Some(211));
		assert!(inverse(10, 492).is_none());
		assert_eq!(inverse(13, 492), Some(265));
		assert!(inverse(16, 492).is_none());
		assert_eq!(inverse(19, 492), Some(259));
		assert!(inverse(22, 492).is_none());
		assert_eq!(inverse(25, 492), Some(433));
		assert!(inverse(28, 492).is_none());
		assert_eq!(inverse(1, 499), Some(1));
		assert_eq!(inverse(4, 499), Some(125));
		assert_eq!(inverse(7, 499), Some(214));
		assert_eq!(inverse(10, 499), Some(50));
		assert_eq!(inverse(13, 499), Some(192));
		assert_eq!(inverse(16, 499), Some(156));
		assert_eq!(inverse(19, 499), Some(394));
		assert_eq!(inverse(22, 499), Some(431));
		assert_eq!(inverse(25, 499), Some(20));
		assert_eq!(inverse(28, 499), Some(303));
		assert_eq!(inverse(1, 506), Some(1));
		assert!(inverse(4, 506).is_none());
		assert_eq!(inverse(7, 506), Some(217));
		assert!(inverse(10, 506).is_none());
		assert_eq!(inverse(13, 506), Some(39));
		assert!(inverse(16, 506).is_none());
		assert_eq!(inverse(19, 506), Some(293));
		assert!(inverse(22, 506).is_none());
		assert_eq!(inverse(25, 506), Some(81));
		assert!(inverse(28, 506).is_none());
		assert_eq!(inverse(1, 513), Some(1));
		assert_eq!(inverse(4, 513), Some(385));
		assert_eq!(inverse(7, 513), Some(220));
		assert_eq!(inverse(10, 513), Some(154));
		assert_eq!(inverse(13, 513), Some(79));
		assert_eq!(inverse(16, 513), Some(481));
		assert!(inverse(19, 513).is_none());
		assert_eq!(inverse(22, 513), Some(70));
		assert_eq!(inverse(25, 513), Some(472));
		assert_eq!(inverse(28, 513), Some(55));
		assert_eq!(inverse(1, 520), Some(1));
		assert!(inverse(4, 520).is_none());
		assert_eq!(inverse(7, 520), Some(223));
		assert!(inverse(10, 520).is_none());
		assert!(inverse(13, 520).is_none());
		assert!(inverse(16, 520).is_none());
		assert_eq!(inverse(19, 520), Some(219));
		assert!(inverse(22, 520).is_none());
		assert!(inverse(25, 520).is_none());
		assert!(inverse(28, 520).is_none());
		assert_eq!(inverse(1, 527), Some(1));
		assert_eq!(inverse(4, 527), Some(132));
		assert_eq!(inverse(7, 527), Some(226));
		assert_eq!(inverse(10, 527), Some(369));
		assert_eq!(inverse(13, 527), Some(446));
		assert_eq!(inverse(16, 527), Some(33));
		assert_eq!(inverse(19, 527), Some(111));
		assert_eq!(inverse(22, 527), Some(24));
		assert_eq!(inverse(25, 527), Some(253));
		assert_eq!(inverse(28, 527), Some(320));
		assert_eq!(inverse(1, 534), Some(1));
		assert!(inverse(4, 534).is_none());
		assert_eq!(inverse(7, 534), Some(229));
		assert!(inverse(10, 534).is_none());
		assert_eq!(inverse(13, 534), Some(493));
		assert!(inverse(16, 534).is_none());
		assert_eq!(inverse(19, 534), Some(253));
		assert!(inverse(22, 534).is_none());
		assert_eq!(inverse(25, 534), Some(235));
		assert!(inverse(28, 534).is_none());
		assert_eq!(inverse(1, 541), Some(1));
		assert_eq!(inverse(4, 541), Some(406));
		assert_eq!(inverse(7, 541), Some(232));
		assert_eq!(inverse(10, 541), Some(487));
		assert_eq!(inverse(13, 541), Some(333));
		assert_eq!(inverse(16, 541), Some(372));
		assert_eq!(inverse(19, 541), Some(57));
		assert_eq!(inverse(22, 541), Some(123));
		assert_eq!(inverse(25, 541), Some(303));
		assert_eq!(inverse(28, 541), Some(58));
		assert_eq!(inverse(1, 548), Some(1));
		assert!(inverse(4, 548).is_none());
		assert_eq!(inverse(7, 548), Some(235));
		assert!(inverse(10, 548).is_none());
		assert_eq!(inverse(13, 548), Some(253));
		assert!(inverse(16, 548).is_none());
		assert_eq!(inverse(19, 548), Some(375));
		assert!(inverse(22, 548).is_none());
		assert_eq!(inverse(25, 548), Some(285));
		assert!(inverse(28, 548).is_none());
		assert_eq!(inverse(1, 555), Some(1));
		assert_eq!(inverse(4, 555), Some(139));
		assert_eq!(inverse(7, 555), Some(238));
		assert!(inverse(10, 555).is_none());
		assert_eq!(inverse(13, 555), Some(427));
		assert_eq!(inverse(16, 555), Some(451));
		assert_eq!(inverse(19, 555), Some(409));
		assert_eq!(inverse(22, 555), Some(328));
		assert!(inverse(25, 555).is_none());
		assert_eq!(inverse(28, 555), Some(337));
		assert_eq!(inverse(1, 562), Some(1));
		assert!(inverse(4, 562).is_none());
		assert_eq!(inverse(7, 562), Some(241));
		assert!(inverse(10, 562).is_none());
		assert_eq!(inverse(13, 562), Some(173));
		assert!(inverse(16, 562).is_none());
		assert_eq!(inverse(19, 562), Some(355));
		assert!(inverse(22, 562).is_none());
		assert_eq!(inverse(25, 562), Some(45));
		assert!(inverse(28, 562).is_none());
		assert_eq!(inverse(1, 569), Some(1));
		assert_eq!(inverse(4, 569), Some(427));
		assert_eq!(inverse(7, 569), Some(244));
		assert_eq!(inverse(10, 569), Some(57));
		assert_eq!(inverse(13, 569), Some(394));
		assert_eq!(inverse(16, 569), Some(249));
		assert_eq!(inverse(19, 569), Some(30));
		assert_eq!(inverse(22, 569), Some(388));
		assert_eq!(inverse(25, 569), Some(478));
		assert_eq!(inverse(28, 569), Some(61));
		assert_eq!(inverse(1, 576), Some(1));
		assert!(inverse(4, 576).is_none());
		assert_eq!(inverse(7, 576), Some(247));
		assert!(inverse(10, 576).is_none());
		assert_eq!(inverse(13, 576), Some(133));
		assert!(inverse(16, 576).is_none());
		assert_eq!(inverse(19, 576), Some(91));
		assert!(inverse(22, 576).is_none());
		assert_eq!(inverse(25, 576), Some(553));
		assert!(inverse(28, 576).is_none());
		assert_eq!(inverse(1, 583), Some(1));
		assert_eq!(inverse(4, 583), Some(146));
		assert_eq!(inverse(7, 583), Some(250));
		assert_eq!(inverse(10, 583), Some(175));
		assert_eq!(inverse(13, 583), Some(314));
		assert_eq!(inverse(16, 583), Some(328));
		assert_eq!(inverse(19, 583), Some(491));
		assert!(inverse(22, 583).is_none());
		assert_eq!(inverse(25, 583), Some(70));
		assert_eq!(inverse(28, 583), Some(354));
		assert_eq!(inverse(1, 590), Some(1));
		assert!(inverse(4, 590).is_none());
		assert_eq!(inverse(7, 590), Some(253));
		assert!(inverse(10, 590).is_none());
		assert_eq!(inverse(13, 590), Some(227));
		assert!(inverse(16, 590).is_none());
		assert_eq!(inverse(19, 590), Some(559));
		assert!(inverse(22, 590).is_none());
		assert!(inverse(25, 590).is_none());
		assert!(inverse(28, 590).is_none());
		assert_eq!(inverse(1, 597), Some(1));
		assert_eq!(inverse(4, 597), Some(448));
		assert_eq!(inverse(7, 597), Some(256));
		assert_eq!(inverse(10, 597), Some(418));
		assert_eq!(inverse(13, 597), Some(46));
		assert_eq!(inverse(16, 597), Some(112));
		assert_eq!(inverse(19, 597), Some(220));
		assert_eq!(inverse(22, 597), Some(190));
		assert_eq!(inverse(25, 597), Some(406));
		assert_eq!(inverse(28, 597), Some(64));
		assert_eq!(inverse(1, 604), Some(1));
		assert!(inverse(4, 604).is_none());
		assert_eq!(inverse(7, 604), Some(259));
		assert!(inverse(10, 604).is_none());
		assert_eq!(inverse(13, 604), Some(93));
		assert!(inverse(16, 604).is_none());
		assert_eq!(inverse(19, 604), Some(159));
		assert!(inverse(22, 604).is_none());
		assert_eq!(inverse(25, 604), Some(145));
		assert!(inverse(28, 604).is_none());
		assert_eq!(inverse(1, 611), Some(1));
		assert_eq!(inverse(4, 611), Some(153));
		assert_eq!(inverse(7, 611), Some(262));
		assert_eq!(inverse(10, 611), Some(550));
		assert!(inverse(13, 611).is_none());
		assert_eq!(inverse(16, 611), Some(191));
		assert_eq!(inverse(19, 611), Some(193));
		assert_eq!(inverse(22, 611), Some(250));
		assert_eq!(inverse(25, 611), Some(220));
		assert_eq!(inverse(28, 611), Some(371));
		assert_eq!(inverse(1, 618), Some(1));
		assert!(inverse(4, 618).is_none());
		assert_eq!(inverse(7, 618), Some(265));
		assert!(inverse(10, 618).is_none());
		assert_eq!(inverse(13, 618), Some(523));
		assert!(inverse(16, 618).is_none());
		assert_eq!(inverse(19, 618), Some(553));
		assert!(inverse(22, 618).is_none());
		assert_eq!(inverse(25, 618), Some(445));
		assert!(inverse(28, 618).is_none());
		assert_eq!(inverse(1, 625), Some(1));
		assert_eq!(inverse(4, 625), Some(469));
		assert_eq!(inverse(7, 625), Some(268));
		assert!(inverse(10, 625).is_none());
		assert_eq!(inverse(13, 625), Some(577));
		assert_eq!(inverse(16, 625), Some(586));
		assert_eq!(inverse(19, 625), Some(329));
		assert_eq!(inverse(22, 625), Some(483));
		assert!(inverse(25, 625).is_none());
		assert_eq!(inverse(28, 625), Some(67));
		assert_eq!(inverse(1, 632), Some(1));
		assert!(inverse(4, 632).is_none());
		assert_eq!(inverse(7, 632), Some(271));
		assert!(inverse(10, 632).is_none());
		assert_eq!(inverse(13, 632), Some(389));
		assert!(inverse(16, 632).is_none());
		assert_eq!(inverse(19, 632), Some(499));
		assert!(inverse(22, 632).is_none());
		assert_eq!(inverse(25, 632), Some(177));
		assert!(inverse(28, 632).is_none());
		assert_eq!(inverse(1, 639), Some(1));
		assert_eq!(inverse(4, 639), Some(160));
		assert_eq!(inverse(7, 639), Some(274));
		assert_eq!(inverse(10, 639), Some(64));
		assert_eq!(inverse(13, 639), Some(295));
		assert_eq!(inverse(16, 639), Some(40));
		assert_eq!(inverse(19, 639), Some(370));
		assert_eq!(inverse(22, 639), Some(610));
		assert_eq!(inverse(25, 639), Some(409));
		assert_eq!(inverse(28, 639), Some(388));
		assert_eq!(inverse(1, 646), Some(1));
		assert!(inverse(4, 646).is_none());
		assert_eq!(inverse(7, 646), Some(277));
		assert!(inverse(10, 646).is_none());
		assert_eq!(inverse(13, 646), Some(497));
		assert!(inverse(16, 646).is_none());
		assert!(inverse(19, 646).is_none());
		assert!(inverse(22, 646).is_none());
		assert_eq!(inverse(25, 646), Some(491));
		assert!(inverse(28, 646).is_none());
		assert_eq!(inverse(1, 653), Some(1));
		assert_eq!(inverse(4, 653), Some(490));
		assert_eq!(inverse(7, 653), Some(280));
		assert_eq!(inverse(10, 653), Some(196));
		assert_eq!(inverse(13, 653), Some(201));
		assert_eq!(inverse(16, 653), Some(449));
		assert_eq!(inverse(19, 653), Some(275));
		assert_eq!(inverse(22, 653), Some(564));
		assert_eq!(inverse(25, 653), Some(209));
		assert_eq!(inverse(28, 653), Some(70));
		assert_eq!(inverse(1, 660), Some(1));
		assert!(inverse(4, 660).is_none());
		assert_eq!(inverse(7, 660), Some(283));
		assert!(inverse(10, 660).is_none());
		assert_eq!(inverse(13, 660), Some(457));
		assert!(inverse(16, 660).is_none());
		assert_eq!(inverse(19, 660), Some(139));
		assert!(inverse(22, 660).is_none());
		assert!(inverse(25, 660).is_none());
		assert!(inverse(28, 660).is_none());
		assert_eq!(inverse(1, 667), Some(1));
		assert_eq!(inverse(4, 667), Some(167));
		assert_eq!(inverse(7, 667), Some(286));
		assert_eq!(inverse(10, 667), Some(467));
		assert_eq!(inverse(13, 667), Some(154));
		assert_eq!(inverse(16, 667), Some(542));
		assert_eq!(inverse(19, 667), Some(316));
		assert_eq!(inverse(22, 667), Some(91));
		assert_eq!(inverse(25, 667), Some(587));
		assert_eq!(inverse(28, 667), Some(405));
		assert_eq!(inverse(1, 674), Some(1));
		assert!(inverse(4, 674).is_none());
		assert_eq!(inverse(7, 674), Some(289));
		assert!(inverse(10, 674).is_none());
		assert_eq!(inverse(13, 674), Some(363));
		assert!(inverse(16, 674).is_none());
		assert_eq!(inverse(19, 674), Some(71));
		assert!(inverse(22, 674).is_none());
		assert_eq!(inverse(25, 674), Some(27));
		assert!(inverse(28, 674).is_none());
		assert_eq!(inverse(1, 681), Some(1));
		assert_eq!(inverse(4, 681), Some(511));
		assert_eq!(inverse(7, 681), Some(292));
		assert_eq!(inverse(10, 681), Some(613));
		assert_eq!(inverse(13, 681), Some(262));
		assert_eq!(inverse(16, 681), Some(298));
		assert_eq!(inverse(19, 681), Some(466));
		assert_eq!(inverse(22, 681), Some(31));
		assert_eq!(inverse(25, 681), Some(109));
		assert_eq!(inverse(28, 681), Some(73));
		assert_eq!(inverse(1, 688), Some(1));
		assert!(inverse(4, 688).is_none());
		assert_eq!(inverse(7, 688), Some(295));
		assert!(inverse(10, 688).is_none());
		assert_eq!(inverse(13, 688), Some(53));
		assert!(inverse(16, 688).is_none());
		assert_eq!(inverse(19, 688), Some(507));
		assert!(inverse(22, 688).is_none());
		assert_eq!(inverse(25, 688), Some(633));
		assert!(inverse(28, 688).is_none());
		assert_eq!(inverse(1, 695), Some(1));
		assert_eq!(inverse(4, 695), Some(174));
		assert_eq!(inverse(7, 695), Some(298));
		assert!(inverse(10, 695).is_none());
		assert_eq!(inverse(13, 695), Some(107));
		assert_eq!(inverse(16, 695), Some(391));
		assert_eq!(inverse(19, 695), Some(439));
		assert_eq!(inverse(22, 695), Some(158));
		assert!(inverse(25, 695).is_none());
		assert_eq!(inverse(28, 695), Some(422));
		assert_eq!(inverse(1, 702), Some(1));
		assert!(inverse(4, 702).is_none());
		assert_eq!(inverse(7, 702), Some(301));
		assert!(inverse(10, 702).is_none());
		assert!(inverse(13, 702).is_none());
		assert!(inverse(16, 702).is_none());
		assert_eq!(inverse(19, 702), Some(37));
		assert!(inverse(22, 702).is_none());
		assert_eq!(inverse(25, 702), Some(337));
		assert!(inverse(28, 702).is_none());
		assert_eq!(inverse(1, 709), Some(1));
		assert_eq!(inverse(4, 709), Some(532));
		assert_eq!(inverse(7, 709), Some(304));
		assert_eq!(inverse(10, 709), Some(71));
		assert_eq!(inverse(13, 709), Some(600));
		assert_eq!(inverse(16, 709), Some(133));
		assert_eq!(inverse(19, 709), Some(112));
		assert_eq!(inverse(22, 709), Some(419));
		assert_eq!(inverse(25, 709), Some(312));
		assert_eq!(inverse(28, 709), Some(76));
		assert_eq!(inverse(1, 716), Some(1));
		assert!(inverse(4, 716).is_none());
		assert_eq!(inverse(7, 716), Some(307));
		assert!(inverse(10, 716).is_none());
		assert_eq!(inverse(13, 716), Some(661));
		assert!(inverse(16, 716).is_none());
		assert_eq!(inverse(19, 716), Some(603));
		assert!(inverse(22, 716).is_none());
		assert_eq!(inverse(25, 716), Some(401));
		assert!(inverse(28, 716).is_none());
		assert_eq!(inverse(1, 723), Some(1));
		assert_eq!(inverse(4, 723), Some(181));
		assert_eq!(inverse(7, 723), Some(310));
		assert_eq!(inverse(10, 723), Some(217));
		assert_eq!(inverse(13, 723), Some(445));
		assert_eq!(inverse(16, 723), Some(226));
		assert_eq!(inverse(19, 723), Some(685));
		assert_eq!(inverse(22, 723), Some(493));
		assert_eq!(inverse(25, 723), Some(376));
		assert_eq!(inverse(28, 723), Some(439));
		assert_eq!(inverse(1, 730), Some(1));
		assert!(inverse(4, 730).is_none());
		assert_eq!(inverse(7, 730), Some(313));
		assert!(inverse(10, 730).is_none());
		assert_eq!(inverse(13, 730), Some(337));
		assert!(inverse(16, 730).is_none());
		assert_eq!(inverse(19, 730), Some(269));
		assert!(inverse(22, 730).is_none());
		assert!(inverse(25, 730).is_none());
		assert!(inverse(28, 730).is_none());
		assert_eq!(inverse(1, 737), Some(1));
		assert_eq!(inverse(4, 737), Some(553));
		assert_eq!(inverse(7, 737), Some(316));
		assert_eq!(inverse(10, 737), Some(516));
		assert_eq!(inverse(13, 737), Some(567));
		assert_eq!(inverse(16, 737), Some(691));
		assert_eq!(inverse(19, 737), Some(194));
		assert!(inverse(22, 737).is_none());
		assert_eq!(inverse(25, 737), Some(59));
		assert_eq!(inverse(28, 737), Some(79));
		assert_eq!(inverse(1, 744), Some(1));
		assert!(inverse(4, 744).is_none());
		assert_eq!(inverse(7, 744), Some(319));
		assert!(inverse(10, 744).is_none());
		assert_eq!(inverse(13, 744), Some(229));
		assert!(inverse(16, 744).is_none());
		assert_eq!(inverse(19, 744), Some(235));
		assert!(inverse(22, 744).is_none());
		assert_eq!(inverse(25, 744), Some(625));
		assert!(inverse(28, 744).is_none());
		assert_eq!(inverse(1, 751), Some(1));
		assert_eq!(inverse(4, 751), Some(188));
		assert_eq!(inverse(7, 751), Some(322));
		assert_eq!(inverse(10, 751), Some(676));
		assert_eq!(inverse(13, 751), Some(520));
		assert_eq!(inverse(16, 751), Some(47));
		assert_eq!(inverse(19, 751), Some(672));
		assert_eq!(inverse(22, 751), Some(239));
		assert_eq!(inverse(25, 751), Some(721));
		assert_eq!(inverse(28, 751), Some(456));
		assert_eq!(inverse(1, 758), Some(1));
		assert!(inverse(4, 758).is_none());
		assert_eq!(inverse(7, 758), Some(325));
		assert!(inverse(10, 758).is_none());
		assert_eq!(inverse(13, 758), Some(175));
		assert!(inverse(16, 758).is_none());
		assert_eq!(inverse(19, 758), Some(399));
		assert!(inverse(22, 758).is_none());
		assert_eq!(inverse(25, 758), Some(91));
		assert!(inverse(28, 758).is_none());
		assert_eq!(inverse(1, 765), Some(1));
		assert_eq!(inverse(4, 765), Some(574));
		assert_eq!(inverse(7, 765), Some(328));
		assert!(inverse(10, 765).is_none());
		assert_eq!(inverse(13, 765), Some(412));
		assert_eq!(inverse(16, 765), Some(526));
		assert_eq!(inverse(19, 765), Some(604));
		assert_eq!(inverse(22, 765), Some(313));
		assert!(inverse(25, 765).is_none());
		assert_eq!(inverse(28, 765), Some(82));
		assert_eq!(inverse(1, 772), Some(1));
		assert!(inverse(4, 772).is_none());
		assert_eq!(inverse(7, 772), Some(331));
		assert!(inverse(10, 772).is_none());
		assert_eq!(inverse(13, 772), Some(297));
		assert!(inverse(16, 772).is_none());
		assert_eq!(inverse(19, 772), Some(447));
		assert!(inverse(22, 772).is_none());
		assert_eq!(inverse(25, 772), Some(525));
		assert!(inverse(28, 772).is_none());
		assert_eq!(inverse(1, 779), Some(1));
		assert_eq!(inverse(4, 779), Some(195));
		assert_eq!(inverse(7, 779), Some(334));
		assert_eq!(inverse(10, 779), Some(78));
		assert_eq!(inverse(13, 779), Some(60));
		assert_eq!(inverse(16, 779), Some(633));
		assert!(inverse(19, 779).is_none());
		assert_eq!(inverse(22, 779), Some(602));
		assert_eq!(inverse(25, 779), Some(187));
		assert_eq!(inverse(28, 779), Some(473));
		assert_eq!(inverse(1, 786), Some(1));
		assert!(inverse(4, 786).is_none());
		assert_eq!(inverse(7, 786), Some(337));
		assert!(inverse(10, 786).is_none());
		assert_eq!(inverse(13, 786), Some(121));
		assert!(inverse(16, 786).is_none());
		assert_eq!(inverse(19, 786), Some(331));
		assert!(inverse(22, 786).is_none());
		assert_eq!(inverse(25, 786), Some(283));
		assert!(inverse(28, 786).is_none());
		assert_eq!(inverse(1, 793), Some(1));
		assert_eq!(inverse(4, 793), Some(595));
		assert_eq!(inverse(7, 793), Some(340));
		assert_eq!(inverse(10, 793), Some(238));
		assert!(inverse(13, 793).is_none());
		assert_eq!(inverse(16, 793), Some(347));
		assert_eq!(inverse(19, 793), Some(167));
		assert_eq!(inverse(22, 793), Some(757));
		assert_eq!(inverse(25, 793), Some(571));
		assert_eq!(inverse(28, 793), Some(85));
		assert_eq!(inverse(1, 800), Some(1));
		assert!(inverse(4, 800).is_none());
		assert_eq!(inverse(7, 800), Some(343));
		assert!(inverse(10, 800).is_none());
		assert_eq!(inverse(13, 800), Some(677));
		assert!(inverse(16, 800).is_none());
		assert_eq!(inverse(19, 800), Some(379));
		assert!(inverse(22, 800).is_none());
		assert!(inverse(25, 800).is_none());
		assert!(inverse(28, 800).is_none());
		assert_eq!(inverse(1, 807), Some(1));
		assert_eq!(inverse(4, 807), Some(202));
		assert_eq!(inverse(7, 807), Some(346));
		assert_eq!(inverse(10, 807), Some(565));
		assert_eq!(inverse(13, 807), Some(745));
		assert_eq!(inverse(16, 807), Some(454));
		assert_eq!(inverse(19, 807), Some(85));
		assert_eq!(inverse(22, 807), Some(697));
		assert_eq!(inverse(25, 807), Some(226));
		assert_eq!(inverse(28, 807), Some(490));
		assert_eq!(inverse(1, 814), Some(1));
		assert!(inverse(4, 814).is_none());
		assert_eq!(inverse(7, 814), Some(349));
		assert!(inverse(10, 814).is_none());
		assert_eq!(inverse(13, 814), Some(501));
		assert!(inverse(16, 814).is_none());
		assert_eq!(inverse(19, 814), Some(557));
		assert!(inverse(22, 814).is_none());
		assert_eq!(inverse(25, 814), Some(521));
		assert!(inverse(28, 814).is_none());
		assert_eq!(inverse(1, 821), Some(1));
		assert_eq!(inverse(4, 821), Some(616));
		assert_eq!(inverse(7, 821), Some(352));
		assert_eq!(inverse(10, 821), Some(739));
		assert_eq!(inverse(13, 821), Some(379));
		assert_eq!(inverse(16, 821), Some(154));
		assert_eq!(inverse(19, 821), Some(605));
		assert_eq!(inverse(22, 821), Some(112));
		assert_eq!(inverse(25, 821), Some(624));
		assert_eq!(inverse(28, 821), Some(88));
		assert_eq!(inverse(1, 828), Some(1));
		assert!(inverse(4, 828).is_none());
		assert_eq!(inverse(7, 828), Some(355));
		assert!(inverse(10, 828).is_none());
		assert_eq!(inverse(13, 828), Some(637));
		assert!(inverse(16, 828).is_none());
		assert_eq!(inverse(19, 828), Some(523));
		assert!(inverse(22, 828).is_none());
		assert_eq!(inverse(25, 828), Some(265));
		assert!(inverse(28, 828).is_none());
		assert_eq!(inverse(1, 835), Some(1));
		assert_eq!(inverse(4, 835), Some(209));
		assert_eq!(inverse(7, 835), Some(358));
		assert!(inverse(10, 835).is_none());
		assert_eq!(inverse(13, 835), Some(257));
		assert_eq!(inverse(16, 835), Some(261));
		assert_eq!(inverse(19, 835), Some(44));
		assert_eq!(inverse(22, 835), Some(38));
		assert!(inverse(25, 835).is_none());
		assert_eq!(inverse(28, 835), Some(507));
		assert_eq!(inverse(1, 842), Some(1));
		assert!(inverse(4, 842).is_none());
		assert_eq!(inverse(7, 842), Some(361));
		assert!(inverse(10, 842).is_none());
		assert_eq!(inverse(13, 842), Some(583));
		assert!(inverse(16, 842).is_none());
		assert_eq!(inverse(19, 842), Some(133));
		assert!(inverse(22, 842).is_none());
		assert_eq!(inverse(25, 842), Some(741));
		assert!(inverse(28, 842).is_none());
		assert_eq!(inverse(1, 849), Some(1));
		assert_eq!(inverse(4, 849), Some(637));
		assert_eq!(inverse(7, 849), Some(364));
		assert_eq!(inverse(10, 849), Some(85));
		assert_eq!(inverse(13, 849), Some(196));
		assert_eq!(inverse(16, 849), Some(796));
		assert_eq!(inverse(19, 849), Some(715));
		assert_eq!(inverse(22, 849), Some(193));
		assert_eq!(inverse(25, 849), Some(34));
		assert_eq!(inverse(28, 849), Some(91));
		assert_eq!(inverse(1, 856), Some(1));
		assert!(inverse(4, 856).is_none());
		assert_eq!(inverse(7, 856), Some(367));
		assert!(inverse(10, 856).is_none());
		assert_eq!(inverse(13, 856), Some(461));
		assert!(inverse(16, 856).is_none());
		assert_eq!(inverse(19, 856), Some(811));
		assert!(inverse(22, 856).is_none());
		assert_eq!(inverse(25, 856), Some(137));
		assert!(inverse(28, 856).is_none());
		assert_eq!(inverse(1, 863), Some(1));
		assert_eq!(inverse(4, 863), Some(216));
		assert_eq!(inverse(7, 863), Some(370));
		assert_eq!(inverse(10, 863), Some(259));
		assert_eq!(inverse(13, 863), Some(332));
		assert_eq!(inverse(16, 863), Some(54));
		assert_eq!(inverse(19, 863), Some(318));
		assert_eq!(inverse(22, 863), Some(510));
		assert_eq!(inverse(25, 863), Some(794));
		assert_eq!(inverse(28, 863), Some(524));
		assert_eq!(inverse(1, 870), Some(1));
		assert!(inverse(4, 870).is_none());
		assert_eq!(inverse(7, 870), Some(373));
		assert!(inverse(10, 870).is_none());
		assert_eq!(inverse(13, 870), Some(67));
		assert!(inverse(16, 870).is_none());
		assert_eq!(inverse(19, 870), Some(229));
		assert!(inverse(22, 870).is_none());
		assert!(inverse(25, 870).is_none());
		assert!(inverse(28, 870).is_none());
		assert_eq!(inverse(1, 877), Some(1));
		assert_eq!(inverse(4, 877), Some(658));
		assert_eq!(inverse(7, 877), Some(376));
		assert_eq!(inverse(10, 877), Some(614));
		assert_eq!(inverse(13, 877), Some(135));
		assert_eq!(inverse(16, 877), Some(603));
		assert_eq!(inverse(19, 877), Some(277));
		assert_eq!(inverse(22, 877), Some(598));
		assert_eq!(inverse(25, 877), Some(421));
		assert_eq!(inverse(28, 877), Some(94));
		assert_eq!(inverse(1, 884), Some(1));
		assert!(inverse(4, 884).is_none());
		assert_eq!(inverse(7, 884), Some(379));
		assert!(inverse(10, 884).is_none());
		assert!(inverse(13, 884).is_none());
		assert!(inverse(16, 884).is_none());
		assert_eq!(inverse(19, 884), Some(791));
		assert!(inverse(22, 884).is_none());
		assert_eq!(inverse(25, 884), Some(389));
		assert!(inverse(28, 884).is_none());
		assert_eq!(inverse(1, 891), Some(1));
		assert_eq!(inverse(4, 891), Some(223));
		assert_eq!(inverse(7, 891), Some(382));
		assert_eq!(inverse(10, 891), Some(802));
		assert_eq!(inverse(13, 891), Some(754));
		assert_eq!(inverse(16, 891), Some(724));
		assert_eq!(inverse(19, 891), Some(469));
		assert!(inverse(22, 891).is_none());
		assert_eq!(inverse(25, 891), Some(499));
		assert_eq!(inverse(28, 891), Some(541));
		assert_eq!(inverse(1, 898), Some(1));
		assert!(inverse(4, 898).is_none());
		assert_eq!(inverse(7, 898), Some(385));
		assert!(inverse(10, 898).is_none());
		assert_eq!(inverse(13, 898), Some(829));
		assert!(inverse(16, 898).is_none());
		assert_eq!(inverse(19, 898), Some(709));
		assert!(inverse(22, 898).is_none());
		assert_eq!(inverse(25, 898), Some(467));
		assert!(inverse(28, 898).is_none());
		assert_eq!(inverse(1, 905), Some(1));
		assert_eq!(inverse(4, 905), Some(679));
		assert_eq!(inverse(7, 905), Some(388));
		assert!(inverse(10, 905).is_none());
		assert_eq!(inverse(13, 905), Some(557));
		assert_eq!(inverse(16, 905), Some(396));
		assert_eq!(inverse(19, 905), Some(524));
		assert_eq!(inverse(22, 905), Some(288));
		assert!(inverse(25, 905).is_none());
		assert_eq!(inverse(28, 905), Some(97));
		assert_eq!(inverse(1, 912), Some(1));
		assert!(inverse(4, 912).is_none());
		assert_eq!(inverse(7, 912), Some(391));
		assert!(inverse(10, 912).is_none());
		assert_eq!(inverse(13, 912), Some(421));
		assert!(inverse(16, 912).is_none());
		assert!(inverse(19, 912).is_none());
		assert!(inverse(22, 912).is_none());
		assert_eq!(inverse(25, 912), Some(73));
		assert!(inverse(28, 912).is_none());
		assert_eq!(inverse(1, 919), Some(1));
		assert_eq!(inverse(4, 919), Some(230));
		assert_eq!(inverse(7, 919), Some(394));
		assert_eq!(inverse(10, 919), Some(92));
		assert_eq!(inverse(13, 919), Some(707));
		assert_eq!(inverse(16, 919), Some(517));
		assert_eq!(inverse(19, 919), Some(387));
		assert_eq!(inverse(22, 919), Some(376));
		assert_eq!(inverse(25, 919), Some(772));
		assert_eq!(inverse(28, 919), Some(558));
		assert_eq!(inverse(1, 926), Some(1));
		assert!(inverse(4, 926).is_none());
		assert_eq!(inverse(7, 926), Some(397));
		assert!(inverse(10, 926).is_none());
		assert_eq!(inverse(13, 926), Some(285));
		assert!(inverse(16, 926).is_none());
		assert_eq!(inverse(19, 926), Some(195));
		assert!(inverse(22, 926).is_none());
		assert_eq!(inverse(25, 926), Some(889));
		assert!(inverse(28, 926).is_none());
		assert_eq!(inverse(1, 933), Some(1));
		assert_eq!(inverse(4, 933), Some(700));
		assert_eq!(inverse(7, 933), Some(400));
		assert_eq!(inverse(10, 933), Some(280));
		assert_eq!(inverse(13, 933), Some(646));
		assert_eq!(inverse(16, 933), Some(175));
		assert_eq!(inverse(19, 933), Some(442));
		assert_eq!(inverse(22, 933), Some(721));
		assert_eq!(inverse(25, 933), Some(112));
		assert_eq!(inverse(28, 933), Some(100));
		assert_eq!(inverse(1, 940), Some(1));
		assert!(inverse(4, 940).is_none());
		assert_eq!(inverse(7, 940), Some(403));
		assert!(inverse(10, 940).is_none());
		assert_eq!(inverse(13, 940), Some(217));
		assert!(inverse(16, 940).is_none());
		assert_eq!(inverse(19, 940), Some(99));
		assert!(inverse(22, 940).is_none());
		assert!(inverse(25, 940).is_none());
		assert!(inverse(28, 940).is_none());
		assert_eq!(inverse(1, 947), Some(1));
		assert_eq!(inverse(4, 947), Some(237));
		assert_eq!(inverse(7, 947), Some(406));
		assert_eq!(inverse(10, 947), Some(663));
		assert_eq!(inverse(13, 947), Some(510));
		assert_eq!(inverse(16, 947), Some(296));
		assert_eq!(inverse(19, 947), Some(648));
		assert_eq!(inverse(22, 947), Some(904));
		assert_eq!(inverse(25, 947), Some(644));
		assert_eq!(inverse(28, 947), Some(575));
		assert_eq!(inverse(1, 954), Some(1));
		assert!(inverse(4, 954).is_none());
		assert_eq!(inverse(7, 954), Some(409));
		assert!(inverse(10, 954).is_none());
		assert_eq!(inverse(13, 954), Some(367));
		assert!(inverse(16, 954).is_none());
		assert_eq!(inverse(19, 954), Some(703));
		assert!(inverse(22, 954).is_none());
		assert_eq!(inverse(25, 954), Some(229));
		assert!(inverse(28, 954).is_none());
		assert_eq!(inverse(1, 961), Some(1));
		assert_eq!(inverse(4, 961), Some(721));
		assert_eq!(inverse(7, 961), Some(412));
		assert_eq!(inverse(10, 961), Some(865));
		assert_eq!(inverse(13, 961), Some(74));
		assert_eq!(inverse(16, 961), Some(901));
		assert_eq!(inverse(19, 961), Some(607));
		assert_eq!(inverse(22, 961), Some(830));
		assert_eq!(inverse(25, 961), Some(346));
		assert_eq!(inverse(28, 961), Some(103));
		assert_eq!(inverse(1, 968), Some(1));
		assert!(inverse(4, 968).is_none());
		assert_eq!(inverse(7, 968), Some(415));
		assert!(inverse(10, 968).is_none());
		assert_eq!(inverse(13, 968), Some(149));
		assert!(inverse(16, 968).is_none());
		assert_eq!(inverse(19, 968), Some(51));
		assert!(inverse(22, 968).is_none());
		assert_eq!(inverse(25, 968), Some(697));
		assert!(inverse(28, 968).is_none());
		assert_eq!(inverse(1, 975), Some(1));
		assert_eq!(inverse(4, 975), Some(244));
		assert_eq!(inverse(7, 975), Some(418));
		assert!(inverse(10, 975).is_none());
		assert!(inverse(13, 975).is_none());
		assert_eq!(inverse(16, 975), Some(61));
		assert_eq!(inverse(19, 975), Some(154));
		assert_eq!(inverse(22, 975), Some(133));
		assert!(inverse(25, 975).is_none());
		assert_eq!(inverse(28, 975), Some(592));
		assert_eq!(inverse(1, 982), Some(1));
		assert!(inverse(4, 982).is_none());
		assert_eq!(inverse(7, 982), Some(421));
		assert!(inverse(10, 982).is_none());
		assert_eq!(inverse(13, 982), Some(831));
		assert!(inverse(16, 982).is_none());
		assert_eq!(inverse(19, 982), Some(827));
		assert!(inverse(22, 982).is_none());
		assert_eq!(inverse(25, 982), Some(275));
		assert!(inverse(28, 982).is_none());
		assert_eq!(inverse(1, 989), Some(1));
		assert_eq!(inverse(4, 989), Some(742));
		assert_eq!(inverse(7, 989), Some(424));
		assert_eq!(inverse(10, 989), Some(99));
		assert_eq!(inverse(13, 989), Some(913));
		assert_eq!(inverse(16, 989), Some(680));
		assert_eq!(inverse(19, 989), Some(937));
		assert_eq!(inverse(22, 989), Some(45));
		assert_eq!(inverse(25, 989), Some(633));
		assert_eq!(inverse(28, 989), Some(106));
		assert_eq!(inverse(1, 996), Some(1));
		assert!(inverse(4, 996).is_none());
		assert_eq!(inverse(7, 996), Some(427));
		assert!(inverse(10, 996).is_none());
		assert_eq!(inverse(13, 996), Some(613));
		assert!(inverse(16, 996).is_none());
		assert_eq!(inverse(19, 996), Some(367));
		assert!(inverse(22, 996).is_none());
		assert_eq!(inverse(25, 996), Some(757));
		assert!(inverse(28, 996).is_none());
	}
}
