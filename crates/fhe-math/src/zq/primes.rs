//! Optimized primes generated as in the NFLlib library.

use fhe_util::is_prime;
use num_bigint::BigUint;

/// Returns whether the modulus supports optimized multiplication and reduction.
/// These optimized operations are possible when the modulus verifies
/// Equation (1) of <https://hal.archives-ouvertes.fr/hal-01242273/document>.
pub fn supports_opt(p: u64) -> bool {
	if p.leading_zeros() < 1 {
		return false;
	}

	// Let's multiply the inequality by (2^s0+1)*2^(3s0):
	// we want to output true when
	//    (2^(3s0)+1) * 2^64 < 2^(3s0) * (2^s0+1) * p
	let mut middle = BigUint::from(1u64) << (3 * p.leading_zeros() as usize);
	let left_side = (&middle + 1u64) << 64;
	middle *= (1u64 << p.leading_zeros()) + 1;
	middle *= p;

	left_side < middle
}

/// Generate a `num_bits`-bit prime, congruent to 1 mod `modulo`, strictly
/// smaller than `upper_bound`. Note that `num_bits` must belong to (10..=62),
/// and upper_bound must be <= 1 << num_bits.
pub fn generate_prime(num_bits: usize, modulo: u64, upper_bound: u64) -> Option<u64> {
	if !(10..=62).contains(&num_bits) {
		None
	} else {
		debug_assert!(
			(1u64 << num_bits) >= upper_bound,
			"upper_bound larger than number of bits"
		);

		let leading_zeros = (64 - num_bits) as u32;

		let mut tentative_prime = upper_bound - 1;
		while tentative_prime % modulo != 1 && tentative_prime.leading_zeros() == leading_zeros {
			tentative_prime -= 1
		}

		while tentative_prime.leading_zeros() == leading_zeros
			&& !is_prime(tentative_prime)
			&& tentative_prime >= modulo
		{
			tentative_prime -= modulo
		}

		if tentative_prime.leading_zeros() == leading_zeros && is_prime(tentative_prime) {
			Some(tentative_prime)
		} else {
			None
		}
	}
}

#[cfg(test)]
mod tests {
	use super::generate_prime;
	use fhe_util::catch_unwind;

	// Verifies that the same moduli as in the NFLlib library are generated.
	// <https://github.com/quarkslab/NFLlib/blob/master/include/nfl/params.hpp>
	#[test]
	fn nfl_62bit_primes() {
		let mut generated = vec![];
		let mut upper_bound = u64::MAX >> 2;
		while generated.len() != 20 {
			let p = generate_prime(62, 2 * 1048576, upper_bound);
			assert!(p.is_some());
			upper_bound = p.unwrap();
			generated.push(upper_bound);
		}
		assert_eq!(
			generated,
			vec![
				4611686018326724609,
				4611686018309947393,
				4611686018282684417,
				4611686018257518593,
				4611686018232352769,
				4611686018171535361,
				4611686018106523649,
				4611686018058289153,
				4611686018051997697,
				4611686017974403073,
				4611686017812922369,
				4611686017781465089,
				4611686017773076481,
				4611686017678704641,
				4611686017666121729,
				4611686017647247361,
				4611686017590624257,
				4611686017554972673,
				4611686017529806849,
				4611686017517223937
			]
		)
	}

	#[test]
	fn upper_bound() {
		debug_assert!(catch_unwind(|| generate_prime(62, 2 * 1048576, (1 << 62) + 1)).is_err());
	}

	#[test]
	fn modulo_too_large() {
		assert!(generate_prime(10, 2048, 1 << 10).is_none());
	}

	#[test]
	fn not_found() {
		// 1033 is the smallest 11-bit prime congruent to 1 modulo 16, so looking for a
		// smaller one should fail.
		assert!(generate_prime(11, 16, 1033).is_none());
	}
}
