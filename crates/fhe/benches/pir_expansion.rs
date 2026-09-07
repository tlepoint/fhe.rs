//! Focused query-expansion benchmarks using the two PIR examples' parameters.

use criterion::{BenchmarkId, Criterion, SamplingMode, criterion_group, criterion_main};
use fhe::bfv::{
    Ciphertext, Encoding, ParametersBuilder, Plaintext, SecretKey, evaluation::EvaluationKey,
    evaluation::EvaluationKeyBuilder,
};

use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use std::{hint::black_box, time::Duration};

fn expansion(c: &mut Criterion) {
    let mut group = c.benchmark_group("pir_expansion");
    group.sampling_mode(SamplingMode::Flat);
    group.warm_up_time(Duration::from_millis(100));
    group.measurement_time(Duration::from_millis(600));

    for (name, degree, plaintext_modulus, moduli, size) in [
        ("mulpir", 8192, 1_785_857, [50, 55, 55], 255usize),
        ("sealpir", 4096, 2_056_193, [36, 36, 37], 511usize),
    ] {
        let par = ParametersBuilder::new()
            .degree(degree)
            .plaintext_modulus(plaintext_modulus)
            .ciphertext_modulus_bits(moduli)
            .build()
            .unwrap();
        let mut rng = ChaCha8Rng::seed_from_u64(degree as u64);
        let sk = SecretKey::generate(&par, &mut rng);
        let rounds = size.next_power_of_two().ilog2() as usize;
        let inverse = fhe_util::inverse(1 << rounds, plaintext_modulus).unwrap();
        let values: Vec<_> = (0..size)
            .map(|i| if i == 3 || i == size - 1 { inverse } else { 0 })
            .collect();
        let pt =
            Plaintext::encode_at_level(&par, values.as_slice(), Encoding::Polynomial, 1).unwrap();
        let ct: Ciphertext = sk.encrypt(&pt, &mut rng).unwrap();
        let ct = Ciphertext::from_bytes(&ct.to_bytes(), &par).unwrap();
        for key_level in [0, 1] {
            let key = EvaluationKeyBuilder::new(&sk)
                .ciphertext_level(1)
                .key_level(key_level)
                .enable_expansion(rounds)
                .build(&mut rng)
                .unwrap();
            let key = EvaluationKey::from_bytes(&key.to_bytes(), &par).unwrap();
            if key_level == 0 {
                // Match the PIR path, including deserialized fixtures. Check
                // every expanded selection bit outside the timed region.
                for (i, expanded) in key.expand(&ct, size).unwrap().iter().enumerate() {
                    let decoded = (sk.decrypt(expanded).unwrap())
                        .decode(Encoding::Polynomial)
                        .unwrap();
                    assert_eq!(
                        decoded.first().copied(),
                        Some(u64::from(i == 3 || i == size - 1))
                    );
                    assert!(decoded.iter().skip(1).all(|x| *x == 0));
                }
                group.sample_size(10);
                group.bench_function(BenchmarkId::new("full", name), |b| {
                    b.iter(|| key.expand(black_box(&ct), black_box(size)).unwrap());
                });
            }
            // One Galois step; key_level=1 is a control with no modulus drop.
            group.sample_size(30);
            group.bench_function(
                BenchmarkId::new("one_step", format!("{name}/key_level={key_level}")),
                |b| b.iter(|| key.expand(black_box(&ct), black_box(2)).unwrap()),
            );
        }
    }
    group.finish();
}

criterion_group!(benches, expansion);
criterion_main!(benches);
