#![feature(int_log)]
#![feature(int_roundings)]

//! Utility functions for the examples

use fhers::bfv::{BfvParameters, Encoding, Plaintext};
use fhers_traits::FheEncoder;
use ndarray::Array2;
use std::{fmt, sync::Arc, time::Duration};
use util::transcode_backward;

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

// Utility macros for timing
#[macro_export]
macro_rules! timeit {
	($name:expr, $code:expr) => {{
		timeit_n!($name, 1, $code)
	}};
}

#[macro_export]
macro_rules! timeit_n {
	($name:expr, $loops:expr, $code:expr) => {{
		use utilities::DisplayDuration;
		let start = std::time::Instant::now();
		for i in 1..$loops {
			let r = $code;
		}
		let r = $code;
		println!(
			"⏱  {}: {}",
			$name,
			DisplayDuration(start.elapsed() / $loops)
		);
		r
	}};
}

// Utility functions for PIR
pub fn generate_database(database_size: usize, elements_size: usize) -> Vec<Vec<u8>> {
	assert!(elements_size >= 4);
	let mut database = vec![vec![0u8; elements_size]; database_size];
	for (i, element) in database.iter_mut().enumerate() {
		element[..4].copy_from_slice(&(i as u32).to_le_bytes());
	}
	database
}

pub fn number_elements_per_plaintext(par: Arc<BfvParameters>, elements_size: usize) -> usize {
	(par.plaintext().ilog2() as usize * par.degree()) / (elements_size * 8)
}

pub fn encode_database(
	database: &Vec<Vec<u8>>,
	par: Arc<BfvParameters>,
	level: usize,
) -> Array2<Plaintext> {
	let elements_size = database[0].len();
	let number_elements_per_plaintext = number_elements_per_plaintext(par.clone(), elements_size);
	let number_rows = database.len().div_ceil(number_elements_per_plaintext);
	println!("number_rows = {}", number_rows);
	println!(
		"number_elements_per_plaintext = {}",
		number_elements_per_plaintext
	);
	let dimension_1 = 1 << ((63 - (number_rows as u64).leading_zeros()).div_ceil(2));
	let dimension_2 = number_rows.div_ceil(dimension_1);
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
		let pt_values = transcode_backward(&serialized_plaintext, par.plaintext().ilog2() as usize);
		preprocessed_database[i] =
			Plaintext::try_encode(&pt_values as &[u64], Encoding::poly_at_level(level), &par)
				.unwrap();
	});
	Array2::from_shape_vec((dimension_1, dimension_2), preprocessed_database).unwrap()
}
