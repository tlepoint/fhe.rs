#![feature(int_log)]
#![feature(int_roundings)]
#![feature(generators, proc_macro_hygiene, stmt_expr_attributes)]

use fhers::bfv;
use fhers_traits::{
	DeserializeParametrized, FheDecoder, FheDecrypter, FheEncoder, FheEncoderVariableTime,
	FheEncrypter, Serialize,
};
use indicatif::HumanBytes;
use itertools::Itertools;
use ndarray::Axis;
use rand::{thread_rng, RngCore};
use std::{error::Error, sync::Arc};
use util::{transcode_backward, transcode_forward};
use utilities::{
	encode_database, generate_database, number_elements_per_plaintext, timeit, timeit_n,
};

fn main() -> Result<(), Box<dyn Error>> {
	let database_size = 1 << 21;
	let elements_size = 288;

	let degree = 4096;
	let plaintext_modulus: u64 = (1 << 22) + 1;
	let moduli_sizes = [33, 38, 38];

	println!("# SealPIR with fhe.rs");

	println!("database_size = {}", database_size);
	println!("elements_size = {}", elements_size);

	// Generation of a random database
	let database = timeit!("Database generation", {
		generate_database(database_size, elements_size)
	});

	// Tower of parameters
	let params = timeit!(
		"Parameters generation",
		Arc::new(
			bfv::BfvParametersBuilder::new()
				.set_degree(degree)
				.set_plaintext_modulus(plaintext_modulus)
				.set_ciphertext_moduli_sizes(&moduli_sizes)
				.build()
				.unwrap()
		)
	);

	// Database preprocessing on the server side
	let preprocessed_database = timeit!("Database preprocessing", {
		encode_database(&database, params.clone(), 1)
	});

	// Client setup
	let (sk, ek_expansion_serialized) = timeit!("Client setup", {
		let sk = bfv::SecretKey::random(&params);
		let dim = preprocessed_database.shape();
		let level = (dim[0] * dim[1]).next_power_of_two().ilog2().div_ceil(2) + 1;
		println!("expansion_level = {}", level);
		let ek_expansion = bfv::LeveledEvaluationKeyBuilder::new(&sk, 1, 0)?
			.enable_expansion(level as usize)?
			.build()?;
		let ek_expansion_serialized = ek_expansion.to_bytes();
		(sk, ek_expansion_serialized)
	});
	println!(
		"📄 Evaluation key (expansion): {}",
		HumanBytes(ek_expansion_serialized.len() as u64)
	);

	// Server setup
	let ek_expansion = timeit!(
		"Server setup",
		bfv::LeveledEvaluationKey::from_bytes(&ek_expansion_serialized, &params)?
	);

	// Client query
	let index = (thread_rng().next_u64() as usize) % database_size;
	let query = timeit!("Client query", {
		let dim = preprocessed_database.shape();
		let level = (dim[0] * dim[1]).next_power_of_two().ilog2().div_ceil(2) + 1;
		let query_index = index / number_elements_per_plaintext(params.clone(), elements_size);
		let mut pt = vec![0u64; dim[0] + dim[1]];
		let inv = util::inverse(1 << level, plaintext_modulus).unwrap();
		pt[query_index / dim[1]] = inv;
		pt[dim[0] + (query_index % dim[1])] = inv;
		let query_pt =
			bfv::Plaintext::try_encode(&pt as &[u64], bfv::Encoding::poly_at_level(1), &params)?;
		let query = sk.try_encrypt(&query_pt)?;
		query.to_bytes()
	});
	println!("📄 Query: {}", HumanBytes(query.len() as u64));

	// Server response
	let responses: Vec<Vec<u8>> = timeit_n!("Server response", 5, {
		let start = std::time::Instant::now();
		let query = bfv::Ciphertext::from_bytes(&query, &params);
		let query = query.unwrap();
		let dim = preprocessed_database.shape();
		let level = (dim[0] * dim[1]).next_power_of_two().ilog2().div_ceil(2) + 1;
		let mut expanded_query = ek_expansion.expands(&query, level as usize)?;
		expanded_query.truncate(dim[0] + dim[1]);
		println!("Expand: {:?}", start.elapsed());
		// println!("expanded_query = {:?}", expanded_query);
		// let now = std::time::Instant::now();
		let fold = preprocessed_database
			.axis_iter(Axis(1))
			.map(|column| {
				// let now = std::time::Instant::now();
				let mut c = bfv::dot_product_scalar(expanded_query[..dim[0]].iter(), column.iter())
					.unwrap();
				// println!("dot_product_scalar: {:?}", now.elapsed());
				// let now = std::time::Instant::now();
				c.mod_switch_to_last_level();
				// println!("mod_switch_to_last_level: {:?}", now.elapsed());
				let c_serialized = c.to_bytes();
				let pt_values =
					transcode_backward(&c_serialized, plaintext_modulus.ilog2() as usize);
				let r = unsafe {
					bfv::PlaintextVec::try_encode_vt(
						&pt_values as &[u64],
						bfv::Encoding::poly_at_level(1),
						&params,
					)
					.unwrap()
					.0
				};
				r
			})
			.collect_vec();
		// println!("pir_dot_product: {:?}", now.elapsed());
		// let now = std::time::Instant::now();
		let v = (0..fold[0].len())
			.map(|i| {
				let mut outi = bfv::dot_product_scalar(
					expanded_query[dim[0]..].iter(),
					fold.iter().map(|pts| pts.get(i).unwrap()).into_iter(),
				)
				.unwrap();
				outi.mod_switch_to_last_level();
				outi.to_bytes()
			})
			.collect_vec();
		// println!("last_step: {:?}", now.elapsed());
		v
	});

	println!(
		"📄 Response: {}",
		HumanBytes(responses.iter().map(|r| r.len()).sum::<usize>() as u64)
	);

	// Client processing
	let answer = timeit!("Client answer", {
		let responses = responses
			.iter()
			.map(|r| bfv::Ciphertext::from_bytes(r, &params).unwrap())
			.collect_vec();
		let decrypted_pt = responses
			.iter()
			.map(|r| sk.try_decrypt(r).unwrap())
			.collect_vec();
		let decrypted_vec = decrypted_pt
			.iter()
			.map(|pt| Vec::<u64>::try_decode(&pt, bfv::Encoding::poly_at_level(2)).unwrap())
			.collect_vec();
		let mut decrypted_ct = vec![];
		for v in &decrypted_vec {
			decrypted_ct.append(&mut transcode_forward(
				v,
				plaintext_modulus.ilog2() as usize,
			));
		}
		// println!("decrypted_ct = {:?}", &decrypted_ct[..30]);
		let mut answer = vec![];
		// There may be a few 0 bytes, let's bruteforce for now
		for i in 0..decrypted_ct.len() {
			if let Ok(ct) =
				bfv::Ciphertext::from_bytes(&decrypted_ct[..decrypted_ct.len() - i], &params)
			{
				let pt = sk.try_decrypt(&ct).unwrap();
				let pt = Vec::<u64>::try_decode(&pt, bfv::Encoding::poly_at_level(2)).unwrap();
				let plaintext = transcode_forward(&pt, plaintext_modulus.ilog2() as usize);
				let offset = index % number_elements_per_plaintext(params.clone(), elements_size);

				println!("Noise in response (ct): {:?}", unsafe {
					sk.measure_noise(&ct)
				});

				answer = plaintext[offset * elements_size..(offset + 1) * elements_size].to_vec();
				break;
			}
		}
		answer
	});

	assert_eq!(&database[index], &answer);

	Ok(())
}
