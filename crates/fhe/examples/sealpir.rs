// Implementation of SealPIR using the `fhe` crate.
//
// SealPIR is a Private Information Retrieval scheme that enables a client to
// retrieve a row from a database without revealing the index to the server.
// SealPIR is described in <https://eprint.iacr.org/2017/1142>.
// We use the same parameters as in Microsoft's public implementation
// <https://github.com/microsoft/SealPIR> to enable an apple-to-apple comparison.

mod util;

use console::style;
use fhe::bfv;
use fhe_math::rq::{traits::TryConvertFrom, Context, Poly, Representation};
use fhe_traits::{
    DeserializeParametrized, FheDecoder, FheDecrypter, FheEncoder, FheEncoderVariableTime,
    FheEncrypter, Serialize,
};
use fhe_util::{div_ceil, ilog2, inverse, transcode_bidirectional, transcode_to_bytes};
use indicatif::HumanBytes;
use itertools::Itertools;
use rand::{rngs::OsRng, thread_rng, RngCore};
use std::{env, error::Error, process::exit, sync::Arc};
use util::{
    encode_database, generate_database, number_elements_per_plaintext,
    timeit::{timeit, timeit_n},
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
    let plaintext_modulus = 2056193;
    let moduli_sizes = [36, 36, 37];

    // Compute what is the maximum byte-length of an element to fit within one
    // ciphertext. Each coefficient of the ciphertext polynomial can contain
    // floor(log2(plaintext_modulus)) bits.
    let max_element_size = (ilog2(plaintext_modulus) * degree) / 8;

    // This executable is a command line tool which enables to specify different
    // database and element sizes.
    let args: Vec<String> = env::args().skip(1).collect();

    // Print the help if requested.
    if args.contains(&"-h".to_string()) || args.contains(&"--help".to_string()) {
        print_notice_and_exit(max_element_size, None)
    }

    // Use the default values from <https://github.com/microsoft/SealPIR>.
    let mut database_size = 1 << 16;
    let mut elements_size = 1024;

    // Update the database size and/or element size depending on the arguments
    // provided.
    for arg in &args {
        if arg.starts_with("--database_size") {
            let a: Vec<&str> = arg.rsplit('=').collect();
            if a.len() != 2 || a[0].parse::<usize>().is_err() {
                print_notice_and_exit(
                    max_element_size,
                    Some("Invalid `--database_size` command".to_string()),
                )
            } else {
                database_size = a[0].parse::<usize>()?
            }
        } else if arg.starts_with("--element_size") {
            let a: Vec<&str> = arg.rsplit('=').collect();
            if a.len() != 2 || a[0].parse::<usize>().is_err() {
                print_notice_and_exit(
                    max_element_size,
                    Some("Invalid `--element_size` command".to_string()),
                )
            } else {
                elements_size = a[0].parse::<usize>()?
            }
        } else {
            print_notice_and_exit(
                max_element_size,
                Some(format!("Unrecognized command: {arg}")),
            )
        }
    }
    if elements_size > max_element_size || elements_size == 0 || database_size == 0 {
        print_notice_and_exit(
            max_element_size,
            Some("Element or database sizes out of bound".to_string()),
        )
    }

    // The parameters are within bound, let's go! Let's first display some
    // information about the database.
    println!("# SealPIR with fhe.rs");
    println!(
        "database of {}",
        HumanBytes((database_size * elements_size) as u64)
    );
    println!("\tdatabase_size = {database_size}");
    println!("\telements_size = {elements_size}");

    // Generation of a random database.
    let database = timeit!("Database generation", {
        generate_database(database_size, elements_size)
    });

    // Let's generate the BFV parameters structure.
    let params = timeit!(
        "Parameters generation",
        bfv::BfvParametersBuilder::new()
            .set_degree(degree)
            .set_plaintext_modulus(plaintext_modulus)
            .set_moduli_sizes(&moduli_sizes)
            .build_arc()?
    );

    // Proprocess the database on the server side: the database will be reshaped
    // so as to pack as many values as possible in every row so that it fits in one
    // ciphertext, and each element will be encoded as a polynomial in Ntt
    // representation.
    let (preprocessed_database, (dim1, dim2)) = timeit!("Database preprocessing", {
        encode_database(&database, params.clone(), 1)
    });

    // Client setup: the client generates a secret key, and an evaluation key for
    // the server will which enable to obliviously expand a ciphertext up to (dim1 +
    // dim2) values, i.e. with expansion level ceil(log2(dim1 + dim2)).
    let (sk, ek_expansion_serialized) = timeit!("Client setup", {
        let sk = bfv::SecretKey::random(&params, &mut OsRng);
        let level = ilog2((dim1 + dim2).next_power_of_two() as u64);
        println!("expansion_level = {level}");
        let ek_expansion = bfv::EvaluationKeyBuilder::new_leveled(&sk, 1, 0)?
            .enable_expansion(level)?
            .build(&mut thread_rng())?;
        let ek_expansion_serialized = ek_expansion.to_bytes();
        (sk, ek_expansion_serialized)
    });
    println!(
        "📄 Evaluation key: {}",
        HumanBytes(ek_expansion_serialized.len() as u64)
    );

    // Server setup: the server receives the evaluation key and deserializes it.
    let ek_expansion = timeit!(
        "Server setup",
        bfv::EvaluationKey::from_bytes(&ek_expansion_serialized, &params)?
    );

    // Client query: when the client wants to retrieve the `index`-th row of the
    // original database, it first computes to which row it corresponds in the
    // original database, and then encrypt a selection vector with 0 everywhere,
    // except at two indices i and (dim1 + j) such that `query_index = i * dim 2 +
    // j` where it sets the value (2^level)^(-1) modulo the plaintext space.
    // It then encodes this vector as a `polynomial` and encrypt the plaintext.
    // The ciphertext is set at level `1`, which means that one of the three moduli
    // has been dropped already; the reason is that the expansion will happen at
    // level 0 (with all three moduli) and then one of the moduli will be dropped
    // to reduce the noise.
    let index = (thread_rng().next_u64() as usize) % database_size;
    let query = timeit!("Client query", {
        let level = ilog2((dim1 + dim2).next_power_of_two() as u64);
        let query_index = index
            / number_elements_per_plaintext(
                params.degree(),
                ilog2(plaintext_modulus),
                elements_size,
            );
        let mut pt = vec![0u64; dim1 + dim2];
        let inv = inverse(1 << level, plaintext_modulus).unwrap();
        pt[query_index / dim2] = inv;
        pt[dim1 + (query_index % dim2)] = inv;
        let query_pt = bfv::Plaintext::try_encode(&pt, bfv::Encoding::poly_at_level(1), &params)?;
        let query: bfv::Ciphertext = sk.try_encrypt(&query_pt, &mut thread_rng())?;
        query.to_bytes()
    });
    println!("📄 Query: {}", HumanBytes(query.len() as u64));

    // Server response: The server receives the query, and after deserializing it,
    // performs the following steps:
    // 1- It expands the query ciphertext into `dim1 + dim2` ciphertexts.
    //    If the client created the query correctly, the server will have obtained
    //    `dim1 + dim2` ciphertexts all encrypting `0`, expect the `i`th and
    //    `dim1 + j`th ones encrypting `1`.
    // 2- It computes the inner product of the first `dim1` ciphertexts with the
    //    columns if the database viewed as a dim1 * dim2 matrix, and modulo-switch
    //    the ciphertext once.
    // 3- It parses the resulting ciphertexts as vector of plaintexts, and compute
    //    the inner product of the last `dim2` ciphertexts from step 1 with the
    //    transposed of the plaintext obtained above.
    // The operation is done `5` times to compute an average response time.
    let responses: Vec<Vec<u8>> = timeit_n!("Server response", 5, {
        let start = std::time::Instant::now();
        let query = bfv::Ciphertext::from_bytes(&query, &params)?;
        let expanded_query = ek_expansion.expands(&query, dim1 + dim2)?;
        println!("Expand: {}", DisplayDuration(start.elapsed()));

        let query_vec = &expanded_query[..dim1];
        let dot_product_mod_switch = move |i, database: &[bfv::Plaintext]| {
            let column = database.iter().skip(i).step_by(dim2);
            let mut c = bfv::dot_product_scalar(query_vec.iter(), column)?;
            c.mod_switch_to_last_level()?;
            Ok(c)
        };

        let dot_products = (0..dim2)
            .map(|i| dot_product_mod_switch(i, &preprocessed_database))
            .collect::<fhe::Result<Vec<bfv::Ciphertext>>>()?;

        let fold = dot_products
            .iter()
            .map(|c| {
                let mut pt_values = Vec::with_capacity(div_ceil(
                    2 * (params.degree() * (64 - params.moduli()[0].leading_zeros() as usize)),
                    ilog2(plaintext_modulus),
                ));
                pt_values.append(&mut transcode_bidirectional(
                    c.get(0).unwrap().coefficients().as_slice().unwrap(),
                    64 - params.moduli()[0].leading_zeros() as usize,
                    ilog2(plaintext_modulus),
                ));
                pt_values.append(&mut transcode_bidirectional(
                    c.get(1).unwrap().coefficients().as_slice().unwrap(),
                    64 - params.moduli()[0].leading_zeros() as usize,
                    ilog2(plaintext_modulus),
                ));
                unsafe {
                    Ok(bfv::PlaintextVec::try_encode_vt(
                        &pt_values,
                        bfv::Encoding::poly_at_level(1),
                        &params,
                    )?
                    .0)
                }
            })
            .collect::<fhe::Result<Vec<Vec<bfv::Plaintext>>>>()?;
        (0..fold[0].len())
            .map(|i| {
                let mut outi = bfv::dot_product_scalar(
                    expanded_query[dim1..].iter(),
                    fold.iter().map(|pts| pts.get(i).unwrap()),
                )?;
                outi.mod_switch_to_last_level()?;
                Ok(outi.to_bytes())
            })
            .collect::<fhe::Result<Vec<Vec<u8>>>>()?
    });
    println!(
        "📄 Response: {}",
        HumanBytes(responses.iter().map(|r| r.len()).sum::<usize>() as u64)
    );

    // Client processing: Upon reception of the response, the client decrypts
    // the ciphertexts and recover the "ciphertexts" which were parsed as plaintext,
    // which it decrypts too. Finally, it outputs the plaintext bytes, offset by the
    // correct value (remember the database was reshaped to maximize how many
    // elements) were embedded in a single ciphertext.
    let answer = timeit!("Client answer", {
        let responses = responses
            .iter()
            .map(|r| bfv::Ciphertext::from_bytes(r, &params).unwrap())
            .collect_vec();
        let decrypted_pt = responses
            .iter()
            .flat_map(|r| sk.try_decrypt(r))
            .collect_vec();
        let decrypted_vec = decrypted_pt
            .iter()
            .flat_map(|pt| Vec::<u64>::try_decode(pt, bfv::Encoding::poly_at_level(2)).unwrap())
            .collect_vec();
        let expect_ncoefficients = div_ceil(
            params.degree() * (64 - params.moduli()[0].leading_zeros() as usize),
            ilog2(plaintext_modulus),
        );
        assert!(decrypted_vec.len() >= 2 * expect_ncoefficients);
        let mut poly0 = transcode_bidirectional(
            &decrypted_vec[..expect_ncoefficients],
            ilog2(plaintext_modulus),
            64 - params.moduli()[0].leading_zeros() as usize,
        );
        let mut poly1 = transcode_bidirectional(
            &decrypted_vec[expect_ncoefficients..2 * expect_ncoefficients],
            ilog2(plaintext_modulus),
            64 - params.moduli()[0].leading_zeros() as usize,
        );
        assert!(poly0.len() >= params.degree());
        assert!(poly1.len() >= params.degree());
        poly0.truncate(params.degree());
        poly1.truncate(params.degree());

        let ctx = Arc::new(Context::new(&params.moduli()[..1], params.degree())?);
        let ct = bfv::Ciphertext::new(
            vec![
                Poly::try_convert_from(poly0, &ctx, true, Representation::Ntt)?,
                Poly::try_convert_from(poly1, &ctx, true, Representation::Ntt)?,
            ],
            &params,
        )?;

        let pt = sk.try_decrypt(&ct).unwrap();
        let pt = Vec::<u64>::try_decode(&pt, bfv::Encoding::poly_at_level(2))?;
        let plaintext = transcode_to_bytes(&pt, ilog2(plaintext_modulus));
        let offset = index
            % number_elements_per_plaintext(
                params.degree(),
                ilog2(plaintext_modulus),
                elements_size,
            );

        println!("Noise in response (ct): {:?}", unsafe {
            sk.measure_noise(&ct)
        });

        plaintext[offset * elements_size..(offset + 1) * elements_size].to_vec()
    });

    // Assert that the answer is indeed the `index`-th element of the initial
    // database.
    assert_eq!(&database[index], &answer);

    Ok(())
}
