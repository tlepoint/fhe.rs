#![feature(int_log)]
#![feature(int_roundings)]
#![feature(generators, proc_macro_hygiene, stmt_expr_attributes)]

use console::style;
use fhers::bfv::{self, Ciphertext, Plaintext};
use fhers_traits::{
	DeserializeParametrized, FheDecoder, FheDecrypter, FheEncoder, FheEncoderVariableTime,
	FheEncrypter, Serialize,
};
use indicatif::HumanBytes;
use itertools::Itertools;
use math::rq::{traits::TryConvertFrom, Context, Poly, Representation};
use rand::{thread_rng, RngCore};
use std::{env, error::Error, process::exit, sync::Arc};
use util::{transcode_bidirectional, transcode_to_bytes};
use utilities::{
	encode_database, generate_database, number_elements_per_plaintext, timeit, timeit_n,
};

fn print_notice_and_exit(max_element_size: usize, error: Option<String>) {
	println!(
		"{} SealPIR with fhe.rs",
		style("  overview:").magenta().bold()
	);
	println!(
		"{} sealpir [-h] [--help] [--database_size=<value>] [--element_size=<value>]",
		style("     usage:").magenta().bold()
	);
	println!(
		"{} {} must be at least 1, and {} must be between 1 and {}",
		style("constraints:").magenta().bold(),
		style("database_size").blue(),
		style("element_size").blue(),
		max_element_size
	);
	if let Some(error) = error {
		println!("{} {}", style("     error:").red().bold(), error);
	}
	exit(0);
}

fn main() -> Result<(), Box<dyn Error>> {
	let degree = 4096usize;
	let plaintext_modulus: u64 = (1 << 22) + 1;
	let moduli_sizes = [36, 36, 37];
	let max_element_size = ((plaintext_modulus.ilog2() * degree as u32) / 8) as usize;

	let args: Vec<String> = env::args().skip(1).collect();

	if args.contains(&"-h".to_string()) || args.contains(&"--help".to_string()) {
		print_notice_and_exit(max_element_size, None)
	}

	let mut database_size = 1 << 21;
	let mut elements_size = 288;
	for arg in &args {
		if arg.starts_with("--database_size") {
			let a: Vec<&str> = arg.rsplit('=').collect();
			if a.len() != 2 || a[0].parse::<usize>().is_err() {
				print_notice_and_exit(
					max_element_size,
					Some("Invalid `--database_size` command".to_string()),
				)
			} else {
				database_size = a[0].parse::<usize>().unwrap()
			}
		} else if arg.starts_with("--element_size") {
			let a: Vec<&str> = arg.rsplit('=').collect();
			if a.len() != 2 || a[0].parse::<usize>().is_err() {
				print_notice_and_exit(
					max_element_size,
					Some("Invalid `--element_size` command".to_string()),
				)
			} else {
				elements_size = a[0].parse::<usize>().unwrap()
			}
		} else {
			print_notice_and_exit(
				max_element_size,
				Some(format!("Unrecognized command: {}", arg)),
			)
		}
	}

	if elements_size > max_element_size || elements_size == 0 || database_size == 0 {
		print_notice_and_exit(
			max_element_size,
			Some("Element or database sizes out of bound".to_string()),
		)
	}

	println!("# SealPIR with fhe.rs");

	println!(
		"database of {}",
		HumanBytes((database_size * elements_size) as u64)
	);
	println!("\tdatabase_size = {}", database_size);
	println!("\telements_size = {}", elements_size);

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
				.set_moduli_sizes(&moduli_sizes)
				.build()
				.unwrap()
		)
	);

	// Database preprocessing on the server side
	let (preprocessed_database, (dim1, dim2)) = timeit!("Database preprocessing", {
		encode_database(&database, params.clone(), 1)
	});

	// Client setup
	let (sk, ek_expansion_serialized) = timeit!("Client setup", {
		let sk = bfv::SecretKey::random(&params);
		let level = (dim1 * dim2).next_power_of_two().ilog2().div_ceil(2) + 1;
		println!("expansion_level = {}", level);
		let ek_expansion = bfv::EvaluationKeyBuilder::new_leveled(&sk, 1, 0)?
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
		bfv::EvaluationKey::from_bytes(&ek_expansion_serialized, &params)?
	);

	// Client query
	let index = (thread_rng().next_u64() as usize) % database_size;
	let query = timeit!("Client query", {
		let level = (dim1 * dim2).next_power_of_two().ilog2().div_ceil(2) + 1;
		let query_index = index
			/ number_elements_per_plaintext(
				params.degree(),
				params.plaintext().ilog2() as usize,
				elements_size,
			);
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
	let responses: Vec<Vec<u8>> = timeit_n!("Server response", 5, {
		let start = std::time::Instant::now();
		let query = bfv::Ciphertext::from_bytes(&query, &params);
		let query = query.unwrap();
		let expanded_query = ek_expansion.expands(&query, dim1 + dim2)?;
		println!("Expand: {:?}", start.elapsed());

		let query_vec = &expanded_query[..dim1];
		let dot_product_mod_switch = move |i, database: &[Plaintext]| {
			let column = database.iter().skip(i).step_by(dim2);
			let mut c = bfv::dot_product_scalar(query_vec.iter(), column)?;
			c.mod_switch_to_last_level();
			Ok(c)
		};

		let dot_products = (0..dim2)
			.map(|i| dot_product_mod_switch(i, &preprocessed_database))
			.collect::<fhers::Result<Vec<Ciphertext>>>()?;

		let fold = dot_products
			.iter()
			.map(|c| {
				let mut pt_values = Vec::with_capacity(
					2 * (params.degree() * (64 - params.moduli()[0].leading_zeros() as usize))
						.div_ceil(plaintext_modulus.ilog2() as usize),
				);
				pt_values.append(&mut transcode_bidirectional(
					c.get(0).unwrap().coefficients().as_slice().unwrap(),
					64 - params.moduli()[0].leading_zeros() as usize,
					plaintext_modulus.ilog2() as usize,
				));
				pt_values.append(&mut transcode_bidirectional(
					c.get(1).unwrap().coefficients().as_slice().unwrap(),
					64 - params.moduli()[0].leading_zeros() as usize,
					plaintext_modulus.ilog2() as usize,
				));
				unsafe {
					Ok(bfv::PlaintextVec::try_encode_vt(
						&pt_values as &[u64],
						bfv::Encoding::poly_at_level(1),
						&params,
					)?
					.0)
				}
			})
			.collect::<fhers::Result<Vec<Vec<Plaintext>>>>()?;
		(0..fold[0].len())
			.map(|i| {
				let mut outi = bfv::dot_product_scalar(
					expanded_query[dim1..].iter(),
					fold.iter().map(|pts| pts.get(i).unwrap()),
				)?;
				outi.mod_switch_to_last_level();
				Ok(outi.to_bytes())
			})
			.collect::<fhers::Result<Vec<Vec<u8>>>>()?
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
			.flat_map(|pt| Vec::<u64>::try_decode(pt, bfv::Encoding::poly_at_level(2)).unwrap())
			.collect_vec();
		let expect_ncoefficients = (params.degree()
			* (64 - params.moduli()[0].leading_zeros() as usize))
			.div_ceil(plaintext_modulus.ilog2() as usize);
		assert!(decrypted_vec.len() >= 2 * expect_ncoefficients);
		let mut poly0 = transcode_bidirectional(
			&decrypted_vec[..expect_ncoefficients],
			plaintext_modulus.ilog2() as usize,
			64 - params.moduli()[0].leading_zeros() as usize,
		);
		let mut poly1 = transcode_bidirectional(
			&decrypted_vec[expect_ncoefficients..2 * expect_ncoefficients],
			plaintext_modulus.ilog2() as usize,
			64 - params.moduli()[0].leading_zeros() as usize,
		);
		assert!(poly0.len() >= params.degree());
		assert!(poly1.len() >= params.degree());
		poly0.truncate(params.degree());
		poly1.truncate(params.degree());

		let ctx = Arc::new(Context::new(&params.moduli()[..1], params.degree())?);
		let ct = Ciphertext::new(
			vec![
				Poly::try_convert_from(poly0, &ctx, true, Representation::Ntt)?,
				Poly::try_convert_from(poly1, &ctx, true, Representation::Ntt)?,
			],
			&params,
		)?;

		let pt = sk.try_decrypt(&ct).unwrap();
		let pt = Vec::<u64>::try_decode(&pt, bfv::Encoding::poly_at_level(2)).unwrap();
		let plaintext = transcode_to_bytes(&pt, plaintext_modulus.ilog2() as usize);
		let offset = index
			% number_elements_per_plaintext(
				params.degree(),
				params.plaintext().ilog2() as usize,
				elements_size,
			);

		println!("Noise in response (ct): {:?}", unsafe {
			sk.measure_noise(&ct)
		});

		plaintext[offset * elements_size..(offset + 1) * elements_size].to_vec()
	});

	assert_eq!(&database[index], &answer);

	Ok(())
}
