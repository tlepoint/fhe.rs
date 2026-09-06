//! A short selection of core BFV operations on two default parameter sets.

use criterion::{BatchSize, BenchmarkId, Criterion, SamplingMode, criterion_group, criterion_main};
use fhe::bfv::{BfvParameters, Ciphertext, Encoding, Plaintext, RelinearizationKey, SecretKey};
use fhe_traits::{FheDecoder, FheDecrypter, FheEncoder, FheEncrypter};
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use std::{hint::black_box, time::Duration};

fn core_bfv(c: &mut Criterion) {
    let mut group = c.benchmark_group("bfv_core");
    group.sample_size(30);
    group.sampling_mode(SamplingMode::Flat);
    group.warm_up_time(Duration::from_millis(100));
    group.measurement_time(Duration::from_millis(600));

    for par in BfvParameters::default_parameters_128(20)
        .unwrap()
        .filter(|par| matches!(par.degree(), 4096 | 8192))
    {
        let mut rng = ChaCha8Rng::seed_from_u64(0xbf00 + par.degree() as u64);
        let sk = SecretKey::random(&par, &mut rng);
        let rk = RelinearizationKey::new(&sk, &mut rng).unwrap();
        let left: Vec<_> = (0..par.degree())
            .map(|i| (17 * i as u64 + 3) % par.plaintext())
            .collect();
        let right: Vec<_> = (0..par.degree())
            .map(|i| (31 * i as u64 + 5) % par.plaintext())
            .collect();
        let pt = Plaintext::try_encode(left.as_slice(), Encoding::simd(), &par).unwrap();
        let other_pt = Plaintext::try_encode(right.as_slice(), Encoding::simd(), &par).unwrap();
        let ct: Ciphertext = sk.try_encrypt(&pt, &mut rng).unwrap();
        let other: Ciphertext = sk.try_encrypt(&other_pt, &mut rng).unwrap();
        let product = &ct * &other;
        let mut relinearized = product.clone();
        rk.relinearizes(&mut relinearized).unwrap();
        let expected: Vec<_> = left
            .iter()
            .zip(&right)
            .map(|(a, b)| a * b % par.plaintext())
            .collect();
        // Verify fixtures before timing; the operations below reuse these
        // operands instead of accumulating noise across benchmark iterations.
        for ciphertext in [&product, &relinearized] {
            assert_eq!(
                Vec::<u64>::try_decode(&sk.try_decrypt(ciphertext).unwrap(), Encoding::simd())
                    .unwrap(),
                expected
            );
        }
        assert_eq!(
            Vec::<u64>::try_decode(&sk.try_decrypt(&ct).unwrap(), Encoding::simd()).unwrap(),
            left
        );
        assert_eq!(
            Vec::<u64>::try_decode(&sk.try_decrypt(&(&ct * &ct)).unwrap(), Encoding::simd())
                .unwrap(),
            left.iter()
                .map(|a| a * a % par.plaintext())
                .collect::<Vec<_>>()
        );

        let parameter = format!(
            "n={}/logq={}",
            par.degree(),
            par.moduli_sizes().iter().sum::<usize>()
        );
        group.bench_function(BenchmarkId::new("encrypt_sk", &parameter), |b| {
            b.iter(|| {
                let encrypted: Ciphertext = sk.try_encrypt(black_box(&pt), &mut rng).unwrap();
                encrypted
            });
        });
        group.bench_function(BenchmarkId::new("decrypt", &parameter), |b| {
            b.iter(|| sk.try_decrypt(black_box(&ct)).unwrap());
        });
        group.bench_function(BenchmarkId::new("multiply", &parameter), |b| {
            b.iter(|| black_box(&ct) * black_box(&other));
        });
        group.bench_function(BenchmarkId::new("square", &parameter), |b| {
            b.iter(|| black_box(&ct) * black_box(&ct));
        });
        group.bench_function(BenchmarkId::new("multiply_relinearize", &parameter), |b| {
            b.iter(|| {
                let mut result = black_box(&ct) * black_box(&other);
                rk.relinearizes(&mut result).unwrap();
                result
            });
        });
        group.bench_function(BenchmarkId::new("relinearize", &parameter), |b| {
            b.iter_batched_ref(
                || product.clone(),
                |input| rk.relinearizes(black_box(input)).unwrap(),
                BatchSize::PerIteration,
            );
        });
    }
    group.finish();
}

criterion_group!(benches, core_bfv);
criterion_main!(benches);
