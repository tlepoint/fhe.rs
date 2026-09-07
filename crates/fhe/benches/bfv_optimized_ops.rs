// Expect indexing in benchmarks for convenience
#![expect(missing_docs, reason = "examples/benches/tests omit docs by design")]

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use fhe::bfv::{
    Ciphertext, Encoding, Parameters, Plaintext, SecretKey, evaluation::dot_product_scalar,
};

use itertools::{Itertools, izip};
use rand::rng;
use std::time::Duration;

pub fn bfv_benchmark(c: &mut Criterion) {
    let mut rng = rng();
    let mut group = c.benchmark_group("bfv_optimized_ops");
    group.sample_size(10);
    group.warm_up_time(Duration::from_secs(1));
    group.measurement_time(Duration::from_secs(1));

    for par in Parameters::profiles_128(20).unwrap() {
        let par = par.build().unwrap();
        for size in [10, 128, 1000] {
            let sk = SecretKey::generate(&par, &mut rng);
            let pt1 =
                Plaintext::encode(&par, &(1..16u64).collect_vec(), Encoding::Polynomial).unwrap();
            let mut c1: Ciphertext = sk.encrypt(&pt1, &mut rng).unwrap();

            let ct_vec = (0..size)
                .map(|i| {
                    let pt =
                        Plaintext::encode(&par, &(i..16u64).collect_vec(), Encoding::Polynomial)
                            .unwrap();
                    sk.encrypt(&pt, &mut rng).unwrap()
                })
                .collect_vec();
            let pt_vec = (0..size)
                .map(|i| {
                    Plaintext::encode(&par, &(i..39u64).collect_vec(), Encoding::Polynomial)
                        .unwrap()
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
                    b.iter(|| {
                        izip!(&ct_vec, &pt_vec).for_each(|(cti, pti)| {
                            c1.add_assign(&(cti.multiply_plaintext(pti).unwrap()))
                                .unwrap()
                        })
                    });
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
                    b.iter(|| dot_product_scalar(&ct_vec, &pt_vec));
                },
            );
        }
    }

    group.finish();
}

criterion_group!(bfv, bfv_benchmark);
criterion_main!(bfv);
