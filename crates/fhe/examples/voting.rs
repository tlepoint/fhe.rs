// Implementation of multiparty voting using the `fhe` crate.

mod util;

use std::{env, error::Error, process::exit, sync::Arc};

use console::style;
use fhe::{
    bfv::{self, Ciphertext, Encoding, Plaintext, PublicKey, SecretKey},
    mbfv::{AggregateIter, CommonRandomPoly, DecryptionShare, PublicKeyShare},
};
use fhe_traits::{FheDecoder, FheEncoder, FheEncrypter};
use rand::{distributions::Uniform, prelude::Distribution, rngs::OsRng, thread_rng};
use rayon::prelude::*;
use util::timeit::timeit;

fn print_notice_and_exit(error: Option<String>) {
    println!(
        "{} Voting with fhe.rs",
        style("  overview:").magenta().bold()
    );
    println!(
        "{} voting [-h] [--help] [--num_voters=<value>] [--num_parties=<value>]",
        style("     usage:").magenta().bold()
    );
    println!(
        "{} {} and {} must be at least 1",
        style("constraints:").magenta().bold(),
        style("num_voters").blue(),
        style("num_parties").blue(),
    );
    if let Some(error) = error {
        println!("{} {}", style("     error:").red().bold(), error);
    }
    exit(0);
}

fn main() -> Result<(), Box<dyn Error>> {
    let degree = 4096;
    let plaintext_modulus: u64 = 4096;
    let moduli = vec![0xffffee001, 0xffffc4001, 0x1ffffe0001];

    // This executable is a command line tool which enables to specify
    // voter/election worker sizes.
    let args: Vec<String> = env::args().skip(1).collect();

    // Print the help if requested.
    if args.contains(&"-h".to_string()) || args.contains(&"--help".to_string()) {
        print_notice_and_exit(None)
    }

    let mut num_voters = 1000;
    let mut num_parties = 10;

    // Update the number of voters and/or number of parties depending on the
    // arguments provided.
    for arg in &args {
        if arg.starts_with("--num_voters") {
            let a: Vec<&str> = arg.rsplit('=').collect();
            if a.len() != 2 || a[0].parse::<usize>().is_err() {
                print_notice_and_exit(Some("Invalid `--num_voters` argument".to_string()))
            } else {
                num_voters = a[0].parse::<usize>()?
            }
        } else if arg.starts_with("--num_parties") {
            let a: Vec<&str> = arg.rsplit('=').collect();
            if a.len() != 2 || a[0].parse::<usize>().is_err() {
                print_notice_and_exit(Some("Invalid `--num_parties` argument".to_string()))
            } else {
                num_parties = a[0].parse::<usize>()?
            }
        } else {
            print_notice_and_exit(Some(format!("Unrecognized argument: {arg}")))
        }
    }

    if num_parties == 0 || num_voters == 0 {
        print_notice_and_exit(Some("Voter and party sizes must be nonzero".to_string()))
    }

    // The parameters are within bound, let's go! Let's first display some
    // information about the vote.
    println!("# Voting with fhe.rs");
    println!("\tnum_voters = {num_voters}");
    println!("\tnum_parties = {num_parties}");

    // Let's generate the BFV parameters structure.
    let params = timeit!(
        "Parameters generation",
        bfv::BfvParametersBuilder::new()
            .set_degree(degree)
            .set_plaintext_modulus(plaintext_modulus)
            .set_moduli(&moduli)
            .build_arc()?
    );
    let crp = CommonRandomPoly::new(&params, &mut thread_rng())?;

    // Party setup: each party generates a secret key and shares of a collective
    // public key.
    struct Party {
        sk_share: SecretKey,
        pk_share: PublicKeyShare,
    }
    let parties: Vec<Party> = timeit!("Party setup", {
        (0..num_parties)
            .into_par_iter()
            .map(|_| {
                let mut rng = OsRng;
                let sk_share = SecretKey::random(&params, &mut rng);
                let mut trng = thread_rng();
                let pk_share = PublicKeyShare::new(&sk_share, crp.clone(), &mut trng)?;
                Ok(Party { sk_share, pk_share })
            })
            .collect::<fhe::Result<Vec<_>>>()?
    });

    // Aggregation: this could be one of the parties or a separate entity. Or the
    // parties can aggregate cooperatively, in a tree-like fashion.
    let pk = timeit!("Public key aggregation", {
        let pk: PublicKey = parties.iter().map(|p| p.pk_share.clone()).aggregate()?;
        pk
    });

    // Vote casting
    let dist = Uniform::new_inclusive(0, 1);
    let votes: Vec<u64> = dist
        .sample_iter(&mut thread_rng())
        .take(num_voters)
        .collect();
    let votes_encrypted: Vec<Ciphertext> = timeit!("Vote casting", {
        votes
            .par_iter()
            .map(|v| {
                let pt = Plaintext::try_encode(&[*v], Encoding::poly(), &params)?;
                let mut rng = thread_rng();
                pk.try_encrypt(&pt, &mut rng)
            })
            .collect::<fhe::Result<Vec<_>>>()?
    });

    // Computing the tally: this can be done by anyone (party, aggregator, separate
    // computing entity).
    let params_clone = params.clone();
    let tally = timeit!("Vote tallying", {
        let sum = votes_encrypted.par_iter().cloned().reduce(
            || Ciphertext::zero(&params_clone),
            |mut acc, ct| {
                acc += &ct;
                acc
            },
        );
        Arc::new(sum)
    });

    // The result of a vote is typically public, so in this scenario the parties can
    // perform a collective decryption. If instead the result of the computation
    // should be kept private, the parties could collectively perform a
    // keyswitch to a different public key.
    let decryption_shares: Vec<DecryptionShare> = timeit!("Decryption", {
        parties
            .par_iter()
            .map(|p| {
                let mut rng = thread_rng();
                DecryptionShare::new(&p.sk_share, &tally, &mut rng)
            })
            .collect::<fhe::Result<Vec<_>>>()?
    });

    // Again, an aggregating party aggregates the decryption shares to produce the
    // decrypted plaintext.
    let tally_pt = timeit!("Decryption share aggregation", {
        let pt: Plaintext = decryption_shares.into_iter().aggregate()?;
        pt
    });
    let tally_vec = Vec::<u64>::try_decode(&tally_pt, Encoding::poly())?;
    let tally_result = tally_vec[0];

    // Show vote result
    println!("Vote result = {} / {}", tally_result, num_voters);

    let expected_tally = votes.iter().sum();
    assert_eq!(tally_result, expected_tally);

    Ok(())
}
