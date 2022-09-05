#![feature(int_log)]

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use fhe::bfv::{BfvParameters, Ciphertext, Encoding, Plaintext, RGSWCiphertext, SecretKey};
use fhe_traits::{FheEncoder, FheEncrypter};
use itertools::Itertools;
use std::time::Duration;

pub fn bfv_rgsw_benchmark(c: &mut Criterion) {
	let mut group = c.benchmark_group("bfv_rgsw");
	group.sample_size(10);
	group.warm_up_time(Duration::from_secs(1));
	group.measurement_time(Duration::from_secs(1));

	for par in &BfvParameters::default_parameters_128(20)[2..] {
		let sk = SecretKey::random(par);

		let pt1 = Plaintext::try_encode(&(1..16u64).collect_vec() as &[u64], Encoding::simd(), par)
			.unwrap();
		let pt2 = Plaintext::try_encode(&(3..39u64).collect_vec() as &[u64], Encoding::simd(), par)
			.unwrap();
		let c1: Ciphertext = sk.try_encrypt(&pt1).unwrap();
		let c2: RGSWCiphertext = sk.try_encrypt(&pt2).unwrap();
		let q = par.moduli_sizes().iter().sum::<usize>();

		group.bench_function(
			BenchmarkId::new("external", format!("n={}/log(q)={}", par.degree(), q)),
			|b| {
				b.iter(|| &c1 * &c2);
			},
		);
	}

	group.finish();
}

criterion_group!(bfv_rgsw, bfv_rgsw_benchmark);
criterion_main!(bfv_rgsw);
