use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use fhe::bfv::{dot_product_scalar, BfvParameters, Ciphertext, Encoding, Plaintext, SecretKey};
use fhe_traits::{FheEncoder, FheEncrypter};
use itertools::{izip, Itertools};
use rand::{rngs::OsRng, thread_rng};
use std::time::Duration;

pub fn bfv_benchmark(c: &mut Criterion) {
	let mut rng = thread_rng();
	let mut group = c.benchmark_group("bfv_optimized_ops");
	group.sample_size(10);
	group.warm_up_time(Duration::from_secs(1));
	group.measurement_time(Duration::from_secs(1));

	for par in &BfvParameters::default_parameters_128(20)[2..] {
		for size in [10, 128, 1000] {
			let sk = SecretKey::random(par, &mut OsRng);
			let pt1 =
				Plaintext::try_encode(&(1..16u64).collect_vec(), Encoding::poly(), par).unwrap();
			let mut c1: Ciphertext = sk.try_encrypt(&pt1, &mut rng).unwrap();

			let ct_vec = (0..size)
				.map(|i| {
					let pt =
						Plaintext::try_encode(&(i..16u64).collect_vec(), Encoding::poly(), par)
							.unwrap();
					sk.try_encrypt(&pt, &mut rng).unwrap()
				})
				.collect_vec();
			let pt_vec = (0..size)
				.map(|i| {
					Plaintext::try_encode(&(i..39u64).collect_vec(), Encoding::poly(), par).unwrap()
				})
				.collect_vec();

			group.bench_function(
				BenchmarkId::new(
					"dot_product/naive",
					format!(
						"size={}/degree={}/logq={}",
						size,
						par.degree(),
						par.moduli_sizes().iter().sum::<usize>()
					),
				),
				|b| {
					b.iter(|| izip!(&ct_vec, &pt_vec).for_each(|(cti, pti)| c1 += &(cti * pti)));
				},
			);

			group.bench_function(
				BenchmarkId::new(
					"dot_product/opt",
					format!(
						"size={}/degree={}/logq={}",
						size,
						par.degree(),
						par.moduli_sizes().iter().sum::<usize>()
					),
				),
				|b| {
					b.iter(|| dot_product_scalar(ct_vec.iter(), pt_vec.iter()));
				},
			);
		}
	}

	group.finish();
}

criterion_group!(bfv, bfv_benchmark);
criterion_main!(bfv);
