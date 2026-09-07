//! A short selection of core BFV operations on two default parameter sets.

use criterion::{BatchSize, BenchmarkId, Criterion, SamplingMode, criterion_group, criterion_main};
use fhe::bfv::{
    BfvParameters, Ciphertext, CiphertextProductAccumulator, Encoding, Plaintext,
    RelinearizationKey, SecretKey,
};
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
        for (label, target) in [("one", 1), ("last", par.max_level())] {
            let round_trips = || {
                let mut parts = black_box(&ct).components().to_vec();
                for _ in 0..target {
                    for p in &mut parts {
                        let mut pb = p.clone().into_power_basis();
                        pb.switch_down().unwrap();
                        *p = pb.into_ntt();
                    }
                }
                Ciphertext::from_components(parts, &par).unwrap()
            };
            let retained_ntt = || {
                let mut result = black_box(&ct).clone();
                result.switch_to_level(target).unwrap();
                result
            };
            assert_eq!(retained_ntt(), round_trips());
            group.bench_function(
                BenchmarkId::new(format!("modulus_switch/{label}/round_trips"), &parameter),
                |b| b.iter(round_trips),
            );
            group.bench_function(
                BenchmarkId::new(format!("modulus_switch/{label}/retained_ntt"), &parameter),
                |b| b.iter(retained_ntt),
            );
        }
        let separate_sum = || {
            let mut sum = Ciphertext::trivial_zero(&par, 0).unwrap();
            for _ in 0..8 {
                sum += &(black_box(&ct) * black_box(&other));
            }
            sum
        };
        let fused_sum = || {
            let mut accumulator = CiphertextProductAccumulator::new(&par, 0).unwrap();
            for _ in 0..8 {
                accumulator
                    .add_product(black_box(&ct), black_box(&other))
                    .unwrap();
            }
            accumulator.finish().unwrap()
        };
        for sum in [separate_sum(), fused_sum()] {
            assert_eq!(
                Vec::<u64>::try_decode(&sk.try_decrypt(&sum).unwrap(), Encoding::simd()).unwrap(),
                expected
                    .iter()
                    .map(|x| 8 * x % par.plaintext())
                    .collect::<Vec<_>>()
            );
        }
        group.bench_function(
            BenchmarkId::new("product_sum_8/separate", &parameter),
            |b| {
                b.iter(separate_sum);
            },
        );
        group.bench_function(BenchmarkId::new("product_sum_8/fused", &parameter), |b| {
            b.iter(fused_sum);
        });
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
