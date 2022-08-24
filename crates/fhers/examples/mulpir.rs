#![feature(int_log)]
#![feature(int_roundings)]
#![feature(generators, proc_macro_hygiene, stmt_expr_attributes)]

use fhers::bfv::{self, Ciphertext, Plaintext};
use fhers_traits::{
	DeserializeParametrized, FheDecoder, FheDecrypter, FheEncoder, FheEncrypter, Serialize,
};
use indicatif::HumanBytes;
use rand::{thread_rng, RngCore};
use std::{error::Error, sync::Arc};
use util::transcode_to_bytes;
use utilities::{
	encode_database, generate_database, number_elements_per_plaintext, timeit, timeit_n,
};

fn main() -> Result<(), Box<dyn Error>> {
	let database_size = 1 << 21;
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
	let (preprocessed_database, (dim1, dim2)) = timeit!("Database preprocessing", {
		encode_database(&database, params.clone(), 1)
	});

	// Client setup
	let (sk, ek_expansion_serialized, ek_relin_serialized) = timeit!("Client setup", {
		let sk = bfv::SecretKey::random(&params);
		let level = (dim1 + dim2).next_power_of_two().ilog2();
		println!("level = {}", level);
		let ek_expansion = bfv::LeveledEvaluationKeyBuilder::new(&sk, 1, 0)?
			.enable_expansion(level as usize)?
			.build()?;
		let ek_relin = bfv::LeveledEvaluationKeyBuilder::new(&sk, 1, 1)?
			.enable_relinearization()?
			.build()?;
		let ek_expansion_serialized = ek_expansion.to_bytes();
		let ek_relin_serialized = ek_relin.to_bytes();
		(sk, ek_expansion_serialized, ek_relin_serialized)
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
			bfv::LeveledEvaluationKey::from_bytes(&ek_expansion_serialized, &params)?,
			bfv::LeveledEvaluationKey::from_bytes(&ek_relin_serialized, &params)?,
		)
	});

	// Client query
	let index = (thread_rng().next_u64() as usize) % database_size;
	let query = timeit!("Client query", {
		let level = (dim1 + dim2).next_power_of_two().ilog2();
		let query_index = index / number_elements_per_plaintext(params.clone(), elements_size);
		let mut pt = vec![0u64; dim1 + dim2];
		let inv = util::inverse(1 << level, plaintext_modulus).unwrap();
		pt[query_index / dim2] = inv;
		pt[dim1 + (query_index % dim2)] = inv;
		let query_pt =
			bfv::Plaintext::try_encode(&pt as &[u64], bfv::Encoding::poly_at_level(1), &params)?;
		let query = sk.try_encrypt(&query_pt)?;
		query.to_bytes()
	});
	println!("📄 Query: {}", HumanBytes(query.len() as u64));

	// Server response
	let response = timeit_n!("Server response", 5, {
		let start = std::time::Instant::now();
		let query = bfv::Ciphertext::from_bytes(&query, &params)?;
		let expanded_query = ek_expansion.expands(&query, dim1 + dim2)?;
		println!("Expand: {:?}", start.elapsed());

		let query_vec = &expanded_query[..dim1];
		let dot_product_mod_switch =
			move |i, database: &[Plaintext]| -> fhers::Result<Ciphertext> {
				let column = database.iter().skip(i).step_by(dim2);
				bfv::dot_product_scalar(query_vec.iter(), column)
			};

		let mut out = bfv::Ciphertext::zero(&params);
		for (i, ci) in expanded_query[dim1..].iter().enumerate() {
			out += &dot_product_mod_switch(i, &preprocessed_database)? * ci
		}
		ek_relin.relinearizes(&mut out)?;
		out.mod_switch_to_last_level();
		out.to_bytes()
	});
	println!("📄 Response: {}", HumanBytes(response.len() as u64));

	// Client processing
	let answer = timeit!("Client answer", {
		let response = bfv::Ciphertext::from_bytes(&response, &params).unwrap();

		let pt = sk.try_decrypt(&response).unwrap();
		let pt = Vec::<u64>::try_decode(&pt, bfv::Encoding::poly_at_level(2)).unwrap();
		let plaintext = transcode_to_bytes(&pt, plaintext_modulus.ilog2() as usize);
		let offset = index % number_elements_per_plaintext(params.clone(), elements_size);

		println!("Noise in response: {:?}", unsafe {
			sk.measure_noise(&response)
		});

		plaintext[offset * elements_size..(offset + 1) * elements_size].to_vec()
	});

	assert_eq!(&database[index], &answer);

	Ok(())
}
