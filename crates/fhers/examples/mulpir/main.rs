#![feature(int_log)]
#![feature(int_roundings)]

use bfv::{
	traits::{Decoder, Decryptor, Deserialize, Encoder, Encryptor, Serialize},
	*,
};
use itertools::izip;
use ndarray::{Array2, Axis};
use rand::{thread_rng, RngCore};
use std::sync::Arc;
use util::{transcode_backward, transcode_forward};

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

	let degree = 8192;
	let plaintext_modulus: u64 = (1 << 22) + 1;
	let moduli_sizes = [50, 55, 55];

	println!("# MulPIR with fhe.rs");

	println!("database_size = {}", database_size);
	println!("elements_size = {}", elements_size);

	// Generation of a random database
	let now = std::time::SystemTime::now();
	let database = generate_database(database_size, elements_size);
	println!("Database generation: {:?}", now.elapsed().unwrap());

	// Server parameters
	let server_params = Arc::new(
		BfvParametersBuilder::new()
			.set_degree(degree)?
			.set_plaintext_modulus(plaintext_modulus)?
			.set_ciphertext_moduli_sizes(&moduli_sizes)?
			.build()?,
	);

	// Database preprocessing on the server side
	let now = std::time::SystemTime::now();
	let preprocessed_database = encode_database(&database, &server_params);
	println!("Database preprocessing: {:?}", now.elapsed().unwrap());

	// BFV parameters
	let now = std::time::SystemTime::now();
	let client_params = Arc::new(
		BfvParametersBuilder::new()
			.set_degree(degree)?
			.set_plaintext_modulus(plaintext_modulus)?
			.set_ciphertext_moduli_sizes(&moduli_sizes)?
			.build()?,
	);
	println!("Parameters generation: {:?}", now.elapsed().unwrap());

	// Client setup
	let now = std::time::SystemTime::now();
	let mut sk = SecretKey::random(&client_params);
	let dim = preprocessed_database.shape();
	let level = (dim[0] * dim[1]).next_power_of_two().ilog2().div_ceil(2) + 1;
	println!("level = {}", level);
	let ek = EvaluationKeyBuilder::new(&sk)
		.enable_expansion(level as usize)?
		.enable_relinearization()?
		.build()?;
	let ek_serialized = ek.serialize();
	println!("Client setup: {:?}", now.elapsed().unwrap());
	println!("Evaluation key: {} B", ek_serialized.len());

	// Server setup
	let now = std::time::SystemTime::now();
	let ek = EvaluationKey::try_deserialize(&ek_serialized, &server_params)?;
	println!("Server setup: {:?}", now.elapsed().unwrap());

	// Client query
	let index = (thread_rng().next_u64() as usize) % database_size;
	println!("Index: {}", index);
	let now = std::time::SystemTime::now();
	let query_index = index / number_elements_per_plaintext(&client_params, elements_size);
	let mut pt = vec![0u64; dim[0] + dim[1]];
	let inv = util::inverse(1 << level, plaintext_modulus).unwrap();
	pt[query_index / dim[1]] = inv;
	pt[dim[0] + (query_index % dim[1])] = inv;
	let query_pt = Plaintext::try_encode(&pt as &[u64], Encoding::Poly, &client_params)?;
	let query = sk.encrypt(&query_pt)?;
	let query_serialized = query.serialize();
	println!("Client query: {:?}", now.elapsed().unwrap());
	println!("Query: {:?} B", query_serialized.len());

	// Server response
	let now = std::time::SystemTime::now();
	let query = Ciphertext::try_deserialize(&query_serialized, &server_params)?;
	let expanded_query = ek.expands(&query, level as usize)?;
	println!("Expand: {:?}", now.elapsed().unwrap());

	let mut out = Ciphertext::zero(&server_params);
	izip!(
		&expanded_query[dim[0]..dim[0] + dim[1]],
		preprocessed_database.axis_iter(Axis(1))
	)
	.for_each(|(cj, column)| {
		let mut c = Ciphertext::zero(&server_params);
		izip!(&expanded_query[..dim[0]], column.iter()).for_each(|(ci, pt)| {
			c += ci * pt;
		});
		c = mul(&c, cj, &ek).unwrap();
		out += &c;
	});
	let response = out.serialize();
	println!("Server response: {:?}", now.elapsed().unwrap());
	println!("Response: {:?} B", response.len());

	// Client processing
	let now = std::time::SystemTime::now();
	let response = Ciphertext::try_deserialize(&response, &client_params).unwrap();

	let pt = sk.decrypt(&response).unwrap();
	let pt = Vec::<u64>::try_decode(&pt, Encoding::Poly).unwrap();
	let plaintext = transcode_forward(&pt, plaintext_modulus.ilog2() as usize);
	println!("Client process: {:?}", now.elapsed().unwrap());
	println!("Noise in response: {:?}", unsafe {
		sk.measure_noise(&response)
	});

	let offset = index % number_elements_per_plaintext(&client_params, elements_size);
	// println!(
	// 	"plaintext: {:?}",
	// 	&plaintext[offset * elements_size..(offset + 1) * elements_size]
	// );
	// println!("plaintext: {:?}", &database[index]);
	assert_eq!(
		&database[index],
		&plaintext[offset * elements_size..(offset + 1) * elements_size]
	);

	Ok(())
}
