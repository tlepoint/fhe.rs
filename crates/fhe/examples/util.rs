// Expect indexing in examples for simplicity
#![expect(
    clippy::indexing_slicing,
    reason = "performance or example code relies on validated indices"
)]

//! Utility functions for the examples

// Example utilities are shared across multiple binaries, so some items are unused per-target.
#![allow(dead_code, unused_imports, unused_macros)]

use fhe::bfv;
use fhe_traits::FheEncoder;
use fhe_util::transcode_from_bytes;
use std::{cmp::min, fmt, sync::Arc, time::Duration};

/// Macros to time code and display a human-readable duration.
pub mod timeit {
    macro_rules! timeit_n {
        ($name:expr, $loops:expr, $code:expr) => {{
            use util::DisplayDuration;
            let start = std::time::Instant::now();
            let r = $code;
            for _ in 1..$loops {
                let _ = $code;
            }
            println!(
                "⏱  {}: {}",
                $name,
                DisplayDuration(start.elapsed() / $loops)
            );
            r
        }};
    }

    macro_rules! timeit {
        ($name:expr, $code:expr) => {{
            use util::DisplayDuration;
            let start = std::time::Instant::now();
            let r = $code;
            println!("⏱  {}: {}", $name, DisplayDuration(start.elapsed()));
            r
        }};
    }

    pub(crate) use timeit;
    pub(crate) use timeit_n;
}

/// Utility struct for displaying human-readable duration of the form "10.5 ms",
/// "350 μs", or "27 ns".
pub struct DisplayDuration(pub Duration);

impl fmt::Display for DisplayDuration {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let duration_ns = self.0.as_nanos();
        if duration_ns < 1_000_u128 {
            write!(f, "{duration_ns} ns")
        } else if duration_ns < 1_000_000_u128 {
            write!(f, "{} μs", (duration_ns + 500) / 1_000)
        } else {
            let duration_ms_times_10 = (duration_ns + 50_000) / (100_000);
            write!(f, "{} ms", (duration_ms_times_10 as f64) / 10.0)
        }
    }
}

// Utility functions for Private Information Retrieval.

/// Generate a database of elements of the form [i || 0...0] where i is the 4B
/// little endian encoding of the index. When the element size is less than 4B,
/// the encoding is truncated.
#[must_use]
pub fn generate_database(database_size: usize, elements_size: usize) -> Vec<Vec<u8>> {
    assert!(database_size > 0 && elements_size > 0);
    let mut database = vec![vec![0u8; elements_size]; database_size];
    for (i, element) in database.iter_mut().enumerate() {
        element[..min(4, elements_size)]
            .copy_from_slice(&(i as u32).to_le_bytes()[..min(4, elements_size)]);
    }
    database
}

/// Compute the number of elements per plaintext given a size configuration.
#[must_use]
pub fn number_elements_per_plaintext(
    degree: usize,
    plaintext_nbits: usize,
    elements_size: usize,
) -> usize {
    (plaintext_nbits * degree) / (elements_size * 8)
}

/// Encode a database into BFV plaintexts, returning the encoded rows and
/// layout.
#[must_use]
pub fn encode_database(
    database: &[Vec<u8>],
    par: Arc<bfv::BfvParameters>,
    level: usize,
) -> (Vec<bfv::Plaintext>, (usize, usize)) {
    assert!(!database.is_empty());

    let elements_size = database[0].len();
    let plaintext_nbits = par.plaintext().ilog2() as usize;
    let number_elements_per_plaintext =
        number_elements_per_plaintext(par.degree(), plaintext_nbits, elements_size);
    let number_rows = database.len().div_ceil(number_elements_per_plaintext);
    println!("number_rows = {number_rows}");
    println!("number_elements_per_plaintext = {number_elements_per_plaintext}");
    let dimension_1 = (number_rows as f64).sqrt().ceil() as usize;
    let dimension_2 = number_rows.div_ceil(dimension_1);
    println!("dimensions = {dimension_1} {dimension_2}");
    println!("dimension = {}", dimension_1 * dimension_2);

    let mut preprocessed_database =
        vec![
            bfv::Plaintext::zero(bfv::Encoding::poly_at_level(level), &par).unwrap();
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
            bfv::Plaintext::try_encode(&pt_values, bfv::Encoding::poly_at_level(level), &par)
                .unwrap();
    });
    (preprocessed_database, (dimension_1, dimension_2))
}

fn main() {}
