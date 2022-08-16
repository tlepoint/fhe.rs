#![feature(int_log)]
#![feature(int_roundings)]
#![feature(generators, proc_macro_hygiene, stmt_expr_attributes)]

use bfv::{traits::*, *};
use indicatif::HumanBytes;
use itertools::Itertools;
use ndarray::Axis;
use rand::{thread_rng, RngCore};
use std::sync::Arc;
use util::{transcode_backward, transcode_forward};
use utilities::{encode_database, generate_database, number_elements_per_plaintext, timeit};

fn main() -> Result<(), String> {
	let database_size = 1 << 21;
	let elements_size = 288;

	let degree = 4096;
	let plaintext_modulus: u64 = (1 << 22) + 1;
	let moduli_sizes = [36, 37, 37];

	println!("# SealPIR with fhe.rs");

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
		encode_database(&database, params[1].clone())
	});

	// Client setup
	let (sk_encrypt, mut sk_1, mut sk_2, ek_expansion_serialized) = timeit!("Client setup", {
		let sk_encrypt = SecretKey::random(&params[0]);
		let dim = preprocessed_database.shape();
		let level = (dim[0] * dim[1]).next_power_of_two().ilog2().div_ceil(2) + 1;
		println!("level = {}", level);
		let ek_expansion = EvaluationKeyBuilder::new(&sk_encrypt)
			.enable_expansion(level as usize)?
			.build()?;
		let mut sk_1 = sk_encrypt.clone();
		sk_1.switch_parameters(&params_switchers[0])?;
		let mut sk_2 = sk_1.clone();
		sk_2.switch_parameters(&params_switchers[1])?;
		let ek_expansion_serialized = ek_expansion.serialize();
		(sk_encrypt, sk_1, sk_2, ek_expansion_serialized)
	});
	println!(
		"📄 Evaluation key (expansion): {}",
		HumanBytes(ek_expansion_serialized.len() as u64)
	);

	// Server setup
	let ek_expansion = timeit!(
		"Server setup",
		EvaluationKey::try_deserialize(&ek_expansion_serialized, &params[0])?
	);

	// Client query
	let index = (thread_rng().next_u64() as usize) % database_size;
	let query = timeit!("Client query", {
		let dim = preprocessed_database.shape();
		let level = (dim[0] * dim[1]).next_power_of_two().ilog2().div_ceil(2) + 1;
		let query_index = index / number_elements_per_plaintext(params[0].clone(), elements_size);
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
	let responses: Vec<Vec<u8>> = timeit!("Server response", {
		let start = std::time::Instant::now();
		let query = Ciphertext::try_deserialize(&query, &params[0])?;
		let dim = preprocessed_database.shape();
		let level = (dim[0] * dim[1]).next_power_of_two().ilog2().div_ceil(2) + 1;
		let mut expanded_query = ek_expansion.expands(&query, level as usize)?;
		expanded_query.truncate(dim[0] + dim[1]);
		println!("Expand: {:?}", start.elapsed());
		expanded_query.iter_mut().for_each(|ct| {
			assert!(ct.switch_parameters(&params_switchers[0]).is_ok());
		});
		println!("Switch parameters: {:?}", start.elapsed());
		let fold = preprocessed_database
			.axis_iter(Axis(1))
			.map(|column| {
				let c = dot_product_scalar(expanded_query[..dim[0]].iter(), column.iter()).unwrap();
				let c_serialized = c.serialize();
				let pt_values =
					transcode_backward(&c_serialized, plaintext_modulus.ilog2() as usize);
				Vec::<Plaintext>::try_encode(&pt_values as &[u64], Encoding::Poly, &params[1])
					.unwrap()
			})
			.collect_vec();
		(0..fold[0].len())
			.map(|i| {
				let mut outi = dot_product_scalar(
					expanded_query[dim[0]..].iter(),
					fold.iter().map(|pts| pts.get(i).unwrap()).into_iter(),
				)
				.unwrap();
				assert!(outi.switch_parameters(&params_switchers[1]).is_ok());
				outi.serialize()
			})
			.collect_vec()
	});

	println!(
		"📄 Response: {}",
		HumanBytes(responses.iter().map(|r| r.len()).sum::<usize>() as u64)
	);

	// Client processing
	let answer = timeit!("Client answer", {
		let responses = responses
			.iter()
			.map(|r| Ciphertext::try_deserialize(r, params.last().unwrap()).unwrap())
			.collect_vec();
		let decrypted_pt = responses
			.iter()
			.map(|r| sk_2.decrypt(r).unwrap())
			.collect_vec();
		let decrypted_vec = decrypted_pt
			.iter()
			.map(|pt| Vec::<u64>::try_decode(&pt, Encoding::Poly).unwrap())
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
			if let Ok(ct) = Ciphertext::try_deserialize(
				&decrypted_ct[..decrypted_ct.len() - i],
				params.get(1).unwrap(),
			) {
				let pt = sk_1.decrypt(&ct).unwrap();
				let pt = Vec::<u64>::try_decode(&pt, Encoding::Poly).unwrap();
				let plaintext = transcode_forward(&pt, plaintext_modulus.ilog2() as usize);
				let offset = index
					% number_elements_per_plaintext(params.last().unwrap().clone(), elements_size);

				println!("Noise in response (ct): {:?}", unsafe {
					sk_1.measure_noise(&ct)
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
