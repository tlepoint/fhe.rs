#![feature(int_log)]
#![feature(int_roundings)]
#![feature(generators, proc_macro_hygiene, stmt_expr_attributes)]

use fhers::bfv;
use fhers_traits::{
	DeserializeParametrized, FheDecoder, FheDecrypter, FheEncoder, FheEncrypter,
	FheParametersSwitchable, Serialize,
};
use indicatif::HumanBytes;
use itertools::izip;
use ndarray::Axis;
use rand::{thread_rng, RngCore};
use std::{error::Error, sync::Arc};
use util::transcode_forward;
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
	let (params, params_switchers) = timeit!("Parameters generation", {
		let mut params = vec![];
		let mut params_switchers = vec![];
		for i in 0..moduli_sizes.len() {
			params.push(Arc::new(
				bfv::BfvParametersBuilder::new()
					.set_degree(degree)
					.set_plaintext_modulus(plaintext_modulus)
					.set_ciphertext_moduli_sizes(&moduli_sizes[..moduli_sizes.len() - i])
					.build()?,
			))
		}
		for i in 0..params.len() - 1 {
			params_switchers.push(bfv::BfvParametersSwitcher::new(&params[i], &params[i + 1])?)
		}
		(params, params_switchers)
	});

	// Database preprocessing on the server side
	let preprocessed_database = timeit!("Database preprocessing", {
		encode_database(&database, params[1].clone())
	});

	// Client setup
	let (sk_encrypt, mut sk, ek_expansion, ek_expansion_serialized, ek_relin, ek_relin_serialized) =
		timeit!("Client setup", {
			let sk_encrypt = bfv::SecretKey::random(&params[1]);
			let dim = preprocessed_database.shape();
			let level = (dim[0] + dim[1]).next_power_of_two().ilog2();
			println!("level = {}", level);
			let ek_expansion = bfv::EvaluationKeyBuilder::new(&sk_encrypt, &params[0])
				.enable_expansion(level as usize)?
				.build()?;
			let mut sk = sk_encrypt.clone();
			let ek_relin = bfv::EvaluationKeyBuilder::new(&sk, &params[1])
				.enable_relinearization()?
				.build()?;
			sk.switch_parameters(&params_switchers[1])?;
			let ek_expansion_serialized = ek_expansion.to_bytes();
			let ek_relin_serialized = ek_relin.to_bytes();
			(
				sk_encrypt,
				sk,
				ek_expansion,
				ek_expansion_serialized,
				ek_relin,
				ek_relin_serialized,
			)
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
	/// TODO: To fix
	// let (ek_expansion, ek_relin) = timeit!("Server setup", {
	// 	(
	// 		bfv::EvaluationKey::from_bytes(&ek_expansion_serialized, &params[0])?,
	// 		bfv::EvaluationKey::from_bytes(&ek_relin_serialized, &params[1])?,
	// 	)
	// });

	// Client query
	let index = (thread_rng().next_u64() as usize) % database_size;
	let query = timeit!("Client query", {
		let dim = preprocessed_database.shape();
		let level = (dim[0] + dim[1]).next_power_of_two().ilog2();
		let query_index = index / number_elements_per_plaintext(params[0].clone(), elements_size);
		let mut pt = vec![0u64; dim[0] + dim[1]];
		let inv = util::inverse(1 << level, plaintext_modulus).unwrap();
		pt[query_index / dim[1]] = inv;
		pt[dim[0] + (query_index % dim[1])] = inv;
		let query_pt = bfv::Plaintext::try_encode(&pt as &[u64], bfv::Encoding::Poly, &params[1])?;
		let query = sk_encrypt.try_encrypt(&query_pt)?;
		query.to_bytes()
	});
	println!("📄 Query: {}", HumanBytes(query.len() as u64));

	// Server response
	let response = timeit_n!("Server response", 5, {
		let start = std::time::Instant::now();
		let query = bfv::Ciphertext::from_bytes(&query, &params[1])?;
		let dim = preprocessed_database.shape();
		let level = (dim[0] + dim[1]).next_power_of_two().ilog2();
		let mut expanded_query = ek_expansion.expands(&query, level as usize)?;
		expanded_query.truncate(dim[0] + dim[1]);
		println!("Expand: {:?}", start.elapsed());
		// expanded_query.iter_mut().for_each(|ct| {
		// 	assert!(ct.switch_parameters(&params_switchers[0]).is_ok());
		// });
		// println!("Switch parameters: {:?}", start.elapsed());
		let mut out = bfv::Ciphertext::zero(&params[1]);
		izip!(
			&expanded_query[dim[0]..],
			preprocessed_database.axis_iter(Axis(1))
		)
		.for_each(|(cj, column)| {
			let c =
				bfv::dot_product_scalar(expanded_query[..dim[0]].iter(), column.iter()).unwrap();
			out += &c * cj;
		});
		ek_relin.relinearizes(&mut out)?;
		out.switch_parameters(&params_switchers[1])?;
		out.to_bytes()
	});
	println!("📄 Response: {}", HumanBytes(response.len() as u64));

	// Client processing
	let answer = timeit!("Client answer", {
		let response = bfv::Ciphertext::from_bytes(&response, &params[2]).unwrap();

		let pt = sk.try_decrypt(&response).unwrap();
		let pt = Vec::<u64>::try_decode(&pt, bfv::Encoding::Poly).unwrap();
		let plaintext = transcode_forward(&pt, plaintext_modulus.ilog2() as usize);
		let offset = index % number_elements_per_plaintext(params[2].clone(), elements_size);

		println!("Noise in response: {:?}", unsafe {
			sk.measure_noise(&response)
		});

		plaintext[offset * elements_size..(offset + 1) * elements_size].to_vec()
	});

	assert_eq!(&database[index], &answer);

	Ok(())
}
