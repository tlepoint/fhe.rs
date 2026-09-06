//! Focused query-expansion benchmarks using the two PIR examples' parameters.

use criterion::{BenchmarkId, Criterion, SamplingMode, criterion_group, criterion_main};
use fhe::bfv::{
    BfvParametersBuilder, Ciphertext, Encoding, EvaluationKey, EvaluationKeyBuilder, Plaintext,
    SecretKey,
};
use fhe_traits::{
    DeserializeParametrized, FheDecoder, FheDecrypter, FheEncoder, FheEncrypter, Serialize,
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
        let par = BfvParametersBuilder::new()
            .set_degree(degree)
            .set_plaintext_modulus(plaintext_modulus)
            .set_moduli_sizes(&moduli)
            .build_arc()
            .unwrap();
        let mut rng = ChaCha8Rng::seed_from_u64(degree as u64);
        let sk = SecretKey::random(&par, &mut rng);
        let rounds = size.next_power_of_two().ilog2() as usize;
        let inverse = fhe_util::inverse(1 << rounds, plaintext_modulus).unwrap();
        let values: Vec<_> = (0..size)
            .map(|i| if i == 3 || i == size - 1 { inverse } else { 0 })
            .collect();
        let pt =
            Plaintext::try_encode(values.as_slice(), Encoding::poly_at_level(1), &par).unwrap();
        let ct: Ciphertext = sk.try_encrypt(&pt, &mut rng).unwrap();
        let ct = Ciphertext::from_bytes(&ct.to_bytes(), &par).unwrap();
        for key_level in [0, 1] {
            let key = EvaluationKeyBuilder::new_leveled(&sk, 1, key_level)
                .unwrap()
                .enable_expansion(rounds)
                .unwrap()
                .build(&mut rng)
                .unwrap();
            let key = EvaluationKey::from_bytes(&key.to_bytes(), &par).unwrap();
            if key_level == 0 {
                // Match the PIR path, including deserialized fixtures. Check
                // every expanded selection bit outside the timed region.
                for (i, expanded) in key.expands(&ct, size).unwrap().iter().enumerate() {
                    let decoded = Vec::<u64>::try_decode(
                        &sk.try_decrypt(expanded).unwrap(),
                        Encoding::poly_at_level(1),
                    )
                    .unwrap();
                    assert_eq!(
                        decoded.first().copied(),
                        Some(u64::from(i == 3 || i == size - 1))
                    );
                    assert!(decoded.iter().skip(1).all(|x| *x == 0));
                }
                group.sample_size(10);
                group.bench_function(BenchmarkId::new("full", name), |b| {
                    b.iter(|| key.expands(black_box(&ct), black_box(size)).unwrap());
                });
            }
            // One Galois step; key_level=1 is a control with no modulus drop.
            group.sample_size(30);
            group.bench_function(
                BenchmarkId::new("one_step", format!("{name}/key_level={key_level}")),
                |b| b.iter(|| key.expands(black_box(&ct), black_box(2)).unwrap()),
            );
        }
    }
    group.finish();
}

criterion_group!(benches, expansion);
criterion_main!(benches);
