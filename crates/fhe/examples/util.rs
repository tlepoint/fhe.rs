// Expect indexing in examples for simplicity
#![expect(
    clippy::indexing_slicing,
    reason = "performance or example code relies on validated indices"
)]

//! Utility functions for the examples

// Example utilities are shared across multiple binaries, so some items are unused per-target.
#![allow(dead_code, unused_imports, unused_macros)]

use fhe::bfv;
use fhe_traits::{FheEncoder, FheEncoderVariableTime};
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

/// Matrix layout for the two stages of a PIR response.
#[derive(Clone, Copy)]
pub enum DatabaseLayout {
    /// Balance the two dimensions.
    Square,
    /// Minimize columns without increasing query-expansion rounds. This reduces
    /// ciphertext multiplications in MulPIR and intermediate ciphertext folding
    /// and modulus switching in SealPIR.
    FewerColumns,
}

impl DatabaseLayout {
    fn dimensions(self, number_rows: usize) -> (usize, usize) {
        assert!(number_rows > 0);
        let mut dimension_1 = number_rows.isqrt();
        if dimension_1 * dimension_1 < number_rows {
            dimension_1 += 1;
        }
        let mut dimension_2 = number_rows.div_ceil(dimension_1);
        if matches!(self, Self::FewerColumns) {
            // Expansion performs 2^ceil(log2(dim1 + dim2)) - 1 key switches,
            // while per-column work scales with dim2: ciphertext multiplication
            // in MulPIR, or modulus switching and folding in SealPIR. Spend
            // unused expansion capacity on a taller matrix to reduce this work.
            let expansion_size = (dimension_1 + dimension_2).next_power_of_two();
            while dimension_2 > 1 {
                let candidate = dimension_2 - 1;
                let rows = number_rows.div_ceil(candidate);
                if rows > expansion_size - candidate {
                    break;
                }
                dimension_1 = rows;
                dimension_2 = candidate;
            }
        }
        (dimension_1, dimension_2)
    }
}

/// Encode a database into BFV plaintexts, returning the encoded rows and
/// layout.
#[must_use]
pub fn encode_database(
    database: &[Vec<u8>],
    par: Arc<bfv::BfvParameters>,
    level: usize,
    layout: DatabaseLayout,
) -> (Vec<bfv::Plaintext>, (usize, usize)) {
    let (rows, columns) = encoded_shape(database, &par, layout);
    let mut data = Vec::with_capacity(rows * columns);
    let shape = encode_database_with(database, par, level, layout, |pt| data.push(pt));
    (data, shape)
}

fn encoded_shape(
    database: &[Vec<u8>],
    par: &bfv::BfvParameters,
    layout: DatabaseLayout,
) -> (usize, usize) {
    assert!(!database.is_empty());
    let per_plaintext = number_elements_per_plaintext(
        par.degree(),
        par.plaintext().ilog2() as usize,
        database[0].len(),
    );
    layout.dimensions(database.len().div_ceil(per_plaintext))
}

/// Database representations supported by the PIR examples.
pub enum EncodedDatabase {
    /// Canonical NTT residues stored in machine words.
    Ntt(Vec<bfv::Plaintext>),
    /// Bit-packed canonical NTT residues.
    Packed(bfv::PackedPlaintextVec),
}

impl EncodedDatabase {
    /// Describe coefficient storage, including packing padding but excluding
    /// metadata and allocator overhead. `par` must match the encoded database.
    #[must_use]
    pub fn storage_summary(&self, par: &bfv::BfvParameters) -> String {
        let (format, bytes) = match self {
            Self::Ntt(data) => (
                "unpacked NTT",
                data.iter()
                    .map(|pt| {
                        par.degree()
                            * (par.moduli().len() - pt.level())
                            * std::mem::size_of::<u64>()
                    })
                    .sum(),
            ),
            Self::Packed(data) => ("packed NTT", data.size_bytes()),
        };
        format!(
            "Encoded database ({format}): {} ({bytes} bytes of coefficient storage)",
            indicatif::HumanBytes(bytes as u64)
        )
    }

    /// Compute one column's ciphertext/plaintext dot product.
    pub fn dot_product<'a>(
        &self,
        workspace: &mut bfv::DotProductScalarWorkspace,
        query: impl Iterator<Item = &'a bfv::Ciphertext> + Clone,
        column: usize,
        columns: usize,
    ) -> fhe::Result<bfv::Ciphertext> {
        match self {
            Self::Ntt(data) => {
                workspace.dot_product_scalar(query, data.iter().skip(column).step_by(columns))
            }
            Self::Packed(data) => workspace
                .dot_product_scalar_packed(query, data.iter().skip(column).step_by(columns)),
        }
    }
}

/// Encode and optionally pack each plaintext before proceeding to the next.
/// The packed mode never constructs a complete unpacked NTT database.
#[must_use]
pub fn prepare_database(
    database: &[Vec<u8>],
    par: Arc<bfv::BfvParameters>,
    level: usize,
    layout: DatabaseLayout,
    packed: bool,
) -> (EncodedDatabase, (usize, usize)) {
    if packed {
        let (rows, columns) = encoded_shape(database, &par, layout);
        let mut data = bfv::PackedPlaintextVec::with_capacity(&par, level, rows * columns).unwrap();
        let shape =
            encode_database_with(database, par, level, layout, |pt| data.push(&pt).unwrap());
        (EncodedDatabase::Packed(data), shape)
    } else {
        let (data, shape) = encode_database(database, par, level, layout);
        (EncodedDatabase::Ntt(data), shape)
    }
}

fn encode_database_with(
    database: &[Vec<u8>],
    par: Arc<bfv::BfvParameters>,
    level: usize,
    layout: DatabaseLayout,
    mut append: impl FnMut(bfv::Plaintext),
) -> (usize, usize) {
    assert!(!database.is_empty());

    let elements_size = database[0].len();
    let plaintext_nbits = par.plaintext().ilog2() as usize;
    let number_elements_per_plaintext =
        number_elements_per_plaintext(par.degree(), plaintext_nbits, elements_size);
    let number_rows = database.len().div_ceil(number_elements_per_plaintext);
    println!("number_rows = {number_rows}");
    println!("number_elements_per_plaintext = {number_elements_per_plaintext}");
    let (dimension_1, dimension_2) = layout.dimensions(number_rows);
    println!("dimensions = {dimension_1} {dimension_2}");
    println!("dimension = {}", dimension_1 * dimension_2);

    // The server database and its padding are public. Explicitly opt into
    // variable-time encoding so public PIR arithmetic retains its optimized
    // path without changing the constant-time default for other plaintexts.
    let variable_time = fhe_traits::VariableTime::new(fhe_traits::PublicData::assert_public());
    let public_zero = bfv::Plaintext::try_encode_vt(
        &[] as &[u64],
        bfv::Encoding::poly_at_level(level),
        &par,
        variable_time,
    )
    .unwrap();
    // Encode populated rows directly instead of allocating a zero polynomial
    // for every row only to replace it immediately.
    (0..number_rows).for_each(|i| {
        let mut serialized_plaintext = vec![0u8; number_elements_per_plaintext * elements_size];
        for j in 0..number_elements_per_plaintext {
            if let Some(pt) = database.get(j + i * number_elements_per_plaintext) {
                serialized_plaintext[j * elements_size..(j + 1) * elements_size].copy_from_slice(pt)
            }
        }
        let pt_values = transcode_from_bytes(&serialized_plaintext, plaintext_nbits);
        append(
            bfv::Plaintext::try_encode_vt(
                pt_values.as_slice(),
                bfv::Encoding::poly_at_level(level),
                &par,
                variable_time,
            )
            .unwrap(),
        );
    });
    for _ in number_rows..dimension_1 * dimension_2 {
        append(public_zero.clone());
    }
    (dimension_1, dimension_2)
}

#[cfg(test)]
mod tests {
    use super::{DatabaseLayout, encode_database, number_elements_per_plaintext};
    use fhe::bfv::{self, BfvParametersBuilder, Encoding};
    use fhe_traits::FheDecoder;
    use fhe_util::transcode_to_bytes;

    #[test]
    fn database_storage_summary_reports_format_size_and_level()
    -> Result<(), Box<dyn std::error::Error>> {
        let params = BfvParametersBuilder::new()
            .set_degree(32)
            .set_plaintext_modulus(1153)
            .set_moduli_sizes(&[40, 40])
            .build_arc()?;
        // Nine records occupy two plaintexts; level one retains one RNS row.
        let database = vec![vec![1u8; 5]; 9];
        for (level, unpacked_bytes, packed_bytes) in [(0, 1024, 648), (1, 512, 328)] {
            for (packed, mode, bytes) in [
                (false, "unpacked NTT", unpacked_bytes),
                (true, "packed NTT", packed_bytes),
            ] {
                let (encoded, shape) = super::prepare_database(
                    &database,
                    params.clone(),
                    level,
                    DatabaseLayout::FewerColumns,
                    packed,
                );
                assert_eq!(shape.0 * shape.1, 2);
                let summary = encoded.storage_summary(&params);
                assert!(summary.starts_with(&format!("Encoded database ({mode}): ")));
                assert!(summary.ends_with(&format!("({bytes} bytes of coefficient storage)")));
            }
        }
        Ok(())
    }

    #[test]
    fn pir_layout_minimizes_columns_within_the_same_expansion_rounds() {
        for number_rows in (1usize..=4096).chain([14085, 16384, 28572, 65535, 65536]) {
            let (square_rows, square_columns) = DatabaseLayout::Square.dimensions(number_rows);
            let (rows, columns) = DatabaseLayout::FewerColumns.dimensions(number_rows);
            let expansion_size = (square_rows + square_columns).next_power_of_two();
            assert!(rows * columns >= number_rows);
            assert!(columns <= square_columns);
            assert_eq!((rows + columns).next_power_of_two(), expansion_size);
            // No smaller column count can fit the same expansion budget.
            for candidate in 1..columns {
                assert!(number_rows.div_ceil(candidate) + candidate > expansion_size);
            }
        }
        assert_eq!(DatabaseLayout::FewerColumns.dimensions(14085), (174, 81));
        assert_eq!(DatabaseLayout::FewerColumns.dimensions(28572), (447, 64));
    }

    #[test]
    fn database_layouts_preserve_every_byte_and_zero_pad() -> Result<(), Box<dyn std::error::Error>>
    {
        let params = BfvParametersBuilder::new()
            .set_degree(32)
            .set_plaintext_modulus(1153)
            .set_moduli_sizes(&[40, 40])
            .build_arc()?;
        let element_size = 5;
        let database: Vec<Vec<u8>> = (0..329)
            .map(|row| {
                (0..element_size)
                    .map(|byte| (row * 31 + byte * 17) as u8)
                    .collect()
            })
            .collect();
        let bits = params.plaintext().ilog2() as usize;
        let per_plaintext = number_elements_per_plaintext(params.degree(), bits, element_size);
        for layout in [DatabaseLayout::Square, DatabaseLayout::FewerColumns] {
            let (encoded, (rows, columns)) = encode_database(&database, params.clone(), 1, layout);
            let (packed, packed_shape) =
                super::prepare_database(&database, params.clone(), 1, layout, true);
            assert_eq!(packed_shape, (rows, columns));
            let super::EncodedDatabase::Packed(packed) = packed else {
                unreachable!()
            };
            assert_eq!(packed.len(), encoded.len());
            let mut workspace = bfv::DotProductScalarWorkspace::new(&params, 1)?;
            let mut rng = rand::rng();
            let sk = bfv::SecretKey::random(&params, &mut rng);
            let ct: bfv::Ciphertext =
                fhe_traits::FheEncrypter::try_encrypt(&sk, encoded.first().unwrap(), &mut rng)?;
            assert_eq!(
                workspace.dot_product_scalar_packed(
                    std::iter::repeat_n(&ct, packed.len()),
                    packed.iter()
                )?,
                workspace
                    .dot_product_scalar(std::iter::repeat_n(&ct, encoded.len()), encoded.iter())?
            );
            assert_eq!(encoded.len(), rows * columns);
            for (row, plaintext) in encoded.iter().enumerate() {
                let coefficients = Vec::<u64>::try_decode(plaintext, Encoding::poly_at_level(1))?;
                let bytes = transcode_to_bytes(&coefficients, bits);
                for (column, element) in bytes
                    .chunks_exact(element_size)
                    .take(per_plaintext)
                    .enumerate()
                {
                    if let Some(expected) = database.get(row * per_plaintext + column) {
                        assert_eq!(element, expected);
                    } else {
                        assert!(element.iter().all(|byte| *byte == 0));
                    }
                }
            }
        }
        Ok(())
    }
}

fn main() {}
