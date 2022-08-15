#![feature(int_log)]
#![feature(int_roundings)]

use bfv::{traits::*, *};
use indicatif::HumanBytes;
use itertools::izip;
use ndarray::{Array2, Axis};
use rand::{thread_rng, RngCore};
use std::{fmt, sync::Arc, time::Duration};
use util::{transcode_backward, transcode_forward};

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

// Utility macro
macro_rules! timeit {
	($name:expr, $code:expr) => {{
		let start = std::time::Instant::now();
		let r = $code;
		println!("⏱  {}: {}", $name, DisplayDuration(start.elapsed()));
		r
	}};
}

// Utility functions for PIR
fn generate_database(database_size: usize, elements_size: usize) -> Vec<Vec<u8>> {
	let mut database = vec![vec![0u8; elements_size]; database_size];
	for i in 0..database.len() {
		database[i][..4].copy_from_slice(&(i as u32).to_le_bytes());
	}
	database
}

fn number_elements_per_plaintext(par: &Arc<BfvParameters>, elements_size: usize) -> usize {
	(par.plaintext().ilog2() as usize * par.degree()) / (elements_size * 8)
}

fn encode_database(database: &Vec<Vec<u8>>, par: &Arc<BfvParameters>) -> Array2<Plaintext> {
	let elements_size = database[0].len();
	let number_elements_per_plaintext = number_elements_per_plaintext(par, elements_size);
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
	let mut preprocessed_database =
		vec![Plaintext::zero(Encoding::Poly, &par); dimension_1 * dimension_2];
	(0..number_rows).for_each(|i| {
		let mut serialized_plaintext = vec![0u8; number_elements_per_plaintext * elements_size];
		for j in 0..number_elements_per_plaintext {
			if let Some(pt) = database.get(j + i * number_elements_per_plaintext) {
				serialized_plaintext[j * elements_size..(j + 1) * elements_size]
					.copy_from_slice(&pt)
			}
		}
		let pt_values = transcode_backward(&serialized_plaintext, par.plaintext().ilog2() as usize);
		let pt = Plaintext::try_encode(&pt_values as &[u64], Encoding::Poly, &par).unwrap();
		preprocessed_database[i] = pt.clone();
	});
	Array2::from_shape_vec((dimension_1, dimension_2), preprocessed_database).unwrap()
}

fn main() -> Result<(), String> {
	let database_size = 1 << 20;
	let elements_size = 288;

	// We use the parameters reported in Table 1 of https://eprint.iacr.org/2019/1483.pdf.
	let degree = 8192;
	let plaintext_modulus: u64 = (1 << 20) + (1 << 19) + (1 << 17) + (1 << 16) + (1 << 14) + 1;
	let moduli_sizes = [50, 55, 55];

	println!("# MulPIR with fhe.rs");

	println!("database_size = {}", database_size);
	println!("elements_size = {}", elements_size);

	// Generation of a random database
	let database = timeit!("Database generation", {
		generate_database(database_size, elements_size)
	});

	// Tower of parameters
	let (params, params_switchers) = timeit!("Parameters generation", {
		let mut params = vec![];
		let mut params_switchers = vec![];
		for i in 0..moduli_sizes.len() {
			params.push(Arc::new(
				BfvParametersBuilder::new()
					.set_degree(degree)?
					.set_plaintext_modulus(plaintext_modulus)?
					.set_ciphertext_moduli_sizes(&moduli_sizes[..moduli_sizes.len() - i])?
					.build()?,
			))
		}
		for i in 0..params.len() - 1 {
			params_switchers.push(BfvParametersSwitcher::new(&params[i], &params[i + 1])?)
		}
		(params, params_switchers)
	});

	// Database preprocessing on the server side
	let preprocessed_database = timeit!("Database preprocessing", {
		encode_database(&database, &params[1])
	});

	// Client setup
	let (sk_encrypt, mut sk, ek_expansion_serialized, ek_relin_serialized) =
		timeit!("Client setup", {
			let sk_encrypt = SecretKey::random(&params[0]);
			let dim = preprocessed_database.shape();
			let level = (dim[0] * dim[1]).next_power_of_two().ilog2().div_ceil(2) + 1;
			println!("level = {}", level);
			let ek_expansion = EvaluationKeyBuilder::new(&sk_encrypt)
				.enable_expansion(level as usize)?
				.build()?;
			let mut sk = sk_encrypt.clone();
			sk.switch_parameters(&params_switchers[0])?;
			let ek_relin = EvaluationKeyBuilder::new(&sk)
				.enable_relinearization()?
				.build()?;
			sk.switch_parameters(&params_switchers[1])?;
			let ek_expansion_serialized = ek_expansion.serialize();
			let ek_relin_serialized = ek_relin.serialize();
			(sk_encrypt, sk, ek_expansion_serialized, ek_relin_serialized)
		});
	println!(
		"📄 Evaluation key (expansion): {}",
		HumanBytes(ek_expansion_serialized.len() as u64)
	);
	println!(
		"📄 Evaluation key (relin): {}",
		HumanBytes(ek_relin_serialized.len() as u64)
	);

	// Server setup
	let (ek_expansion, ek_relin) = timeit!("Server setup", {
		(
			EvaluationKey::try_deserialize(&ek_expansion_serialized, &params[0])?,
			EvaluationKey::try_deserialize(&ek_relin_serialized, &params[1])?,
		)
	});

	// Client query
	let index = (thread_rng().next_u64() as usize) % database_size;
	let query = timeit!("Client query", {
		let dim = preprocessed_database.shape();
		let level = (dim[0] * dim[1]).next_power_of_two().ilog2().div_ceil(2) + 1;
		let query_index = index / number_elements_per_plaintext(&params[0], elements_size);
		let mut pt = vec![0u64; dim[0] + dim[1]];
		let inv = util::inverse(1 << level, plaintext_modulus).unwrap();
		pt[query_index / dim[1]] = inv;
		pt[dim[0] + (query_index % dim[1])] = inv;
		let query_pt = Plaintext::try_encode(&pt as &[u64], Encoding::Poly, &params[0])?;
		let query = sk_encrypt.encrypt(&query_pt)?;
		query.serialize()
	});
	println!("📄 Query: {}", HumanBytes(query.len() as u64));

	// Server response
	let response = timeit!("Server response", {
		let start = std::time::Instant::now();
		let query = Ciphertext::try_deserialize(&query, &params[0])?;
		let dim = preprocessed_database.shape();
		let level = (dim[0] * dim[1]).next_power_of_two().ilog2().div_ceil(2) + 1;
		let mut expanded_query = ek_expansion.expands(&query, level as usize)?;
		println!("Expand: {:?}", start.elapsed());
		expanded_query.iter_mut().for_each(|ct| {
			assert!(ct.switch_parameters(&params_switchers[0]).is_ok());
		});
		println!("Switch parameters: {:?}", start.elapsed());
		let mut out = Ciphertext::zero(&params[1]);
		izip!(
			&expanded_query[dim[0]..dim[0] + dim[1]],
			preprocessed_database.axis_iter(Axis(1))
		)
		.for_each(|(cj, column)| {
			let mut c = Ciphertext::zero(&params[1]);
			izip!(&expanded_query[..dim[0]], column.iter()).for_each(|(ci, pt)| {
				c += ci * pt;
			});
			c = mul(&c, cj, &ek_relin).unwrap();
			out += &c;
		});
		out.switch_parameters(&params_switchers[1])?;
		out.serialize()
	});
	println!("📄 Response: {}", HumanBytes(response.len() as u64));

	// Client processing
	let answer = timeit!("Client answer", {
		let response = Ciphertext::try_deserialize(&response, &params[2]).unwrap();

		let pt = sk.decrypt(&response).unwrap();
		let pt = Vec::<u64>::try_decode(&pt, Encoding::Poly).unwrap();
		let plaintext = transcode_forward(&pt, plaintext_modulus.ilog2() as usize);
		let offset = index % number_elements_per_plaintext(&params[2], elements_size);

		println!("Noise in response: {:?}", unsafe {
			sk.measure_noise(&response)
		});

		plaintext[offset * elements_size..(offset + 1) * elements_size].to_vec()
	});

	assert_eq!(&database[index], &answer);

	Ok(())
}
