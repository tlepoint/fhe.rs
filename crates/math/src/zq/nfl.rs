//! Optimized primes generated as in the NFLlib library.

use num_bigint::BigUint;
use util::is_prime;

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

/// Generate an optimized `num_bits`-bit prime, congruent to 1 mod `modulo`, smaller than `upper_bound`.
/// Note that `num_bits` must belong to (30..=62), and upper_bound must be <= 1 << num_bits.
pub fn generate_opt_prime(num_bits: usize, modulo: u64, upper_bound: u64) -> Option<u64> {
	if !(30..=62).contains(&num_bits) {
		None
	} else {
		debug_assert!(
			(1u64 << num_bits) >= upper_bound,
			"upper_bound larger than number of bits"
		);

		let mut tentative_prime = upper_bound - 1;
		while tentative_prime % modulo != 1
			&& tentative_prime.leading_zeros() == ((64 - num_bits) as u32)
		{
			tentative_prime -= 1
		}
		if tentative_prime.leading_zeros() != ((64 - num_bits) as u32) {
			return None;
		}

		while tentative_prime.leading_zeros() == (64 - num_bits) as u32
			&& !is_prime(tentative_prime)
		{
			tentative_prime -= modulo
		}

		if tentative_prime.leading_zeros() == (64 - num_bits) as u32
			&& supports_opt(tentative_prime)
		{
			Some(tentative_prime)
		} else {
			None
		}
	}
}

#[cfg(test)]
mod tests {
	use super::generate_opt_prime;

	// Verifies that the same moduli as in the NFLlib library are generated.
	// <https://github.com/quarkslab/NFLlib/blob/master/include/nfl/params.hpp>
	#[test]
	fn test_nfl_62bit_primes() {
		let mut generated = vec![];
		let mut upper_bound = u64::MAX >> 2;
		while generated.len() != 20 {
			let p = generate_opt_prime(62, 2 * 1048576, upper_bound);
			assert!(p.is_some());
			upper_bound = p.unwrap();
			generated.push(p.unwrap());
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
}
