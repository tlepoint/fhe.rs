//! Repeated-operand multiplication, including the cost of preparing it.

use criterion::{BenchmarkId, Criterion, SamplingMode, criterion_group, criterion_main};
use fhe::bfv::{
    Ciphertext, Encoding, Parameters, Plaintext, SecretKey, evaluation::MultiplicationPlan,
};

use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use std::{hint::black_box, time::Duration};

fn multiplication(c: &mut Criterion) {
    let mut group = c.benchmark_group("bfv_multiplication");
    group.sample_size(30);
    group.sampling_mode(SamplingMode::Flat);
    group.warm_up_time(Duration::from_millis(200));
    group.measurement_time(Duration::from_secs(1));
    for par in Parameters::profiles_128(20)
        .unwrap()
        .filter(|p| matches!(p.degree(), 4096 | 8192))
    {
        let par = par.build().unwrap();
        let mut rng = ChaCha8Rng::seed_from_u64(0x0f4e + par.degree() as u64);
        let sk = SecretKey::generate(&par, &mut rng);
        let values: Vec<_> = (0..par.degree()).map(|i| i as u64 % 19).collect();
        let pt = Plaintext::encode(&par, &values, Encoding::Simd).unwrap();
        let ct: Ciphertext = sk.encrypt(&pt, &mut rng).unwrap();
        let rhs: Ciphertext = sk.encrypt(&pt, &mut rng).unwrap();
        let rhs_batch: Vec<Ciphertext> =
            (0..8).map(|_| sk.encrypt(&pt, &mut rng).unwrap()).collect();
        let ct_copy = ct.clone();
        let strategy = MultiplicationPlan::builder(&par).level(0).build().unwrap();
        let prepared = strategy.prepare_lhs(&ct).unwrap();
        let product = strategy.multiply(&ct, &rhs).unwrap();
        assert_eq!(prepared.multiply(&rhs).unwrap(), product);
        assert_eq!(
            (sk.decrypt(&product).unwrap())
                .decode(Encoding::Simd)
                .unwrap(),
            values
                .iter()
                .map(|x| x * x % par.plaintext_modulus_u64().unwrap())
                .collect::<Vec<_>>()
        );
        assert_eq!(strategy.square(&ct).unwrap(), ct.square().unwrap());
        assert_eq!(
            strategy.multiply(&ct, &ct_copy).unwrap(),
            ct.square().unwrap()
        );
        for rhs in &rhs_batch {
            assert_eq!(
                prepared.multiply(rhs).unwrap(),
                strategy.multiply(&ct, rhs).unwrap()
            );
        }
        let parameter = format!(
            "n={}/logq={}",
            par.degree(),
            par.moduli_sizes().iter().sum::<usize>()
        );
        group.bench_function(BenchmarkId::new("ordinary", &parameter), |b| {
            b.iter(|| strategy.multiply(black_box(&ct), black_box(&rhs)).unwrap());
        });
        group.bench_function(BenchmarkId::new("prepared", &parameter), |b| {
            b.iter(|| prepared.multiply(black_box(&rhs)).unwrap());
        });
        group.bench_function(BenchmarkId::new("prepare_only", &parameter), |b| {
            b.iter(|| strategy.prepare_lhs(black_box(&ct)).unwrap());
        });
        group.bench_function(BenchmarkId::new("ordinary_8", &parameter), |b| {
            b.iter(|| {
                for rhs in &rhs_batch {
                    black_box(strategy.multiply(black_box(&ct), black_box(rhs)).unwrap());
                }
            });
        });
        group.bench_function(
            BenchmarkId::new("prepared_8_including_setup", &parameter),
            |b| {
                b.iter(|| {
                    let prepared = strategy.prepare_lhs(black_box(&ct)).unwrap();
                    for rhs in &rhs_batch {
                        black_box(prepared.multiply(black_box(rhs)).unwrap());
                    }
                });
            },
        );
        group.bench_function(BenchmarkId::new("square", &parameter), |b| {
            b.iter(|| strategy.square(black_box(&ct)).unwrap());
        });
        group.bench_function(BenchmarkId::new("square_via_multiply", &parameter), |b| {
            b.iter(|| {
                strategy
                    .multiply(black_box(&ct), black_box(&ct_copy))
                    .unwrap()
            });
        });
    }
    group.finish();
}

criterion_group!(benches, multiplication);
criterion_main!(benches);
