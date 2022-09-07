//! Utility functions for the examples

use fhe::bfv::{BfvParameters, Encoding, Plaintext};
use fhe_traits::FheEncoder;
use fhe_util::transcode_from_bytes;
use std::{cmp::min, fmt, sync::Arc, time::Duration};

// Utility for displaying duration
pub struct DisplayDuration(pub Duration);

impl fmt::Display for DisplayDuration {
	fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
		let duration_us = (self.0.as_secs_f64() * 1000000.0).round();
		if duration_us >= 1000.0 {
			write!(f, "{} ms", (duration_us / 100.0).round() / 10.0)
		} else {
			write!(f, "{} us", duration_us)
		}
	}
}

// Utility functions for PIR
pub fn generate_database(database_size: usize, elements_size: usize) -> Vec<Vec<u8>> {
	assert!(elements_size >= 4);
	let mut database = vec![vec![0u8; elements_size]; database_size];
	for (i, element) in database.iter_mut().enumerate() {
		element[..min(4, elements_size)]
			.copy_from_slice(&(i as u32).to_le_bytes()[..min(4, elements_size)]);
	}
	database
}

pub fn number_elements_per_plaintext(
	degree: usize,
	plaintext_nbits: usize,
	elements_size: usize,
) -> usize {
	(plaintext_nbits * degree) / (elements_size * 8)
}

pub fn encode_database(
	database: &Vec<Vec<u8>>,
	par: Arc<BfvParameters>,
	level: usize,
) -> (Vec<Plaintext>, (usize, usize)) {
	let elements_size = database[0].len();
	// We want to capture the number of bits b such that 2^b <= par.plaintext()
	// to simulate the `.ilog2()` function.
	// For this, we compute 63 - par.plaintext().leading_zeros(). Indeed, when 2^b
	// <= par.plaintext() < 2^(b+1), then par.plaintext().leading_zeros() = 64 - (b
	// - 1).
	let plaintext_nbits = 63 - par.plaintext().leading_zeros() as usize;
	let number_elements_per_plaintext =
		number_elements_per_plaintext(par.degree(), plaintext_nbits, elements_size);
	let number_rows =
		(database.len() + number_elements_per_plaintext - 1) / number_elements_per_plaintext;
	println!("number_rows = {}", number_rows);
	println!(
		"number_elements_per_plaintext = {}",
		number_elements_per_plaintext
	);
	let dimension_1 = (number_rows as f64).sqrt().ceil() as usize;
	let dimension_2 = (number_rows + dimension_1 - 1) / dimension_1;
	println!("dimensions = {} {}", dimension_1, dimension_2);
	println!("dimension = {}", dimension_1 * dimension_2);
	let mut preprocessed_database = vec![
		Plaintext::zero(Encoding::poly_at_level(level), &par)
			.unwrap();
		dimension_1 * dimension_2
	];
	(0..number_rows).for_each(|i| {
		let mut serialized_plaintext = vec![0u8; number_elements_per_plaintext * elements_size];
		for j in 0..number_elements_per_plaintext {
			if let Some(pt) = database.get(j + i * number_elements_per_plaintext) {
				serialized_plaintext[j * elements_size..(j + 1) * elements_size].copy_from_slice(pt)
			}
		}
		let pt_values = transcode_from_bytes(&serialized_plaintext, plaintext_nbits);
		preprocessed_database[i] =
			Plaintext::try_encode(&pt_values as &[u64], Encoding::poly_at_level(level), &par)
				.unwrap();
	});
	(preprocessed_database, (dimension_1, dimension_2))
}

#[allow(dead_code)]
fn main() {}
