// Expect indexing in benchmarks for convenience
#![expect(missing_docs, reason = "examples/benches/tests omit docs by design")]

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use fhe::bfv::{
    Ciphertext, Encoding, Parameters, Plaintext, SecretKey, evaluation::RgswCiphertext,
};

use itertools::Itertools;
use rand::rng;
use std::time::Duration;

pub fn bfv_rgsw_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("bfv_rgsw");
    group.sample_size(10);
    group.warm_up_time(Duration::from_secs(1));
    group.measurement_time(Duration::from_secs(1));

    for par in Parameters::profiles_128(20).unwrap() {
        let par = par.build().unwrap();
        let mut rng = rng();
        let sk = SecretKey::generate(&par, &mut rng);

        let pt1 = Plaintext::encode(&par, &(1..16u64).collect_vec(), Encoding::Simd).unwrap();
        let pt2 = Plaintext::encode(&par, &(3..39u64).collect_vec(), Encoding::Simd).unwrap();
        let c1: Ciphertext = sk.encrypt(&pt1, &mut rng).unwrap();
        let c2: RgswCiphertext = sk.encrypt_rgsw(&pt2, &mut rng).unwrap();
        let q = par.moduli_sizes().iter().sum::<usize>();

        group.bench_function(
            BenchmarkId::new("external", format!("n={}/log(q)={}", par.degree(), q)),
            |b| {
                b.iter(|| c1.multiply_rgsw(&c2).unwrap());
            },
        );
    }

    group.finish();
}

criterion_group!(bfv_rgsw, bfv_rgsw_benchmark);
criterion_main!(bfv_rgsw);
