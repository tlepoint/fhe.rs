// Expect indexing in benchmarks for convenience
#![expect(missing_docs, reason = "examples/benches/tests omit docs by design")]
#![expect(
    clippy::indexing_slicing,
    reason = "performance or example code relies on validated indices"
)]

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use fhe::bfv::{
    Ciphertext, Encoding, Parameters, Plaintext, PublicKey, SecretKey,
    evaluation::EvaluationKeyBuilder, evaluation::MultiplicationPlan,
    evaluation::RelinearizationKey,
};
use fhe_math::rns::{RnsContext, ScalingFactor};
use fhe_math::zq::primes::generate_prime;

use itertools::Itertools;
use num_bigint::BigUint;
use rand::rng;
use std::time::Duration;

pub fn bfv_benchmark(c: &mut Criterion) {
    let mut rng = rng();
    let mut group = c.benchmark_group("bfv");
    group.sample_size(10);
    group.warm_up_time(Duration::from_millis(600));
    group.measurement_time(Duration::from_millis(1000));

    for par in Parameters::profiles_128(20).unwrap() {
        let par = par.build().unwrap();
        let sk = SecretKey::generate(&par, &mut rng);
        let ek = if par.moduli().len() > 1 {
            Some(
                EvaluationKeyBuilder::new(&sk)
                    .enable_inner_sum()
                    .enable_column_rotation(1)
                    .enable_expansion(par.degree().ilog2() as usize)
                    .build(&mut rng)
                    .unwrap(),
            )
        } else {
            None
        };

        let rk = if par.moduli().len() > 1 {
            Some(RelinearizationKey::new(&sk, &mut rng).unwrap())
        } else {
            None
        };

        let pt1 = Plaintext::encode(&par, &(1..16u64).collect_vec(), Encoding::Simd).unwrap();
        let pt2 = Plaintext::encode(&par, &(3..39u64).collect_vec(), Encoding::Simd).unwrap();
        let mut c1: Ciphertext = sk.encrypt(&pt1, &mut rng).unwrap();
        let c2: Ciphertext = sk.encrypt(&pt2, &mut rng).unwrap();

        let q = par.moduli_sizes().iter().sum::<usize>();

        group.bench_function(
            BenchmarkId::new("keygen_sk", format!("n={}/log(q)={}", par.degree(), q)),
            |b| {
                b.iter(|| SecretKey::generate(&par, &mut rng));
            },
        );

        group.bench_function(
            BenchmarkId::new("keygen_pk", format!("n={}/log(q)={}", par.degree(), q)),
            |b| {
                b.iter(|| PublicKey::from_secret_key(&sk, &mut rng));
            },
        );

        group.bench_function(
            BenchmarkId::new("keygen_rk", format!("n={}/log(q)={}", par.degree(), q)),
            |b| {
                b.iter(|| RelinearizationKey::new(&sk, &mut rng));
            },
        );

        group.bench_function(
            BenchmarkId::new("encode_poly", format!("n={}/log(q)={}", par.degree(), q)),
            |b| {
                b.iter(|| Plaintext::encode(&par, &(1..16u64).collect_vec(), Encoding::Polynomial));
            },
        );

        group.bench_function(
            BenchmarkId::new("encode_simd", format!("n={}/log(q)={}", par.degree(), q)),
            |b| {
                b.iter(|| Plaintext::encode(&par, &(1..16u64).collect_vec(), Encoding::Simd));
            },
        );

        group.bench_function(
            BenchmarkId::new("encrypt_sk", format!("n={}/log(q)={}", par.degree(), q)),
            |b| {
                b.iter(|| {
                    let _: fhe::Result<Ciphertext> = sk.encrypt(&pt1, &mut rng);
                });
            },
        );

        group.bench_function(
            BenchmarkId::new("add_ct", format!("n={}/log(q)={}", par.degree(), q)),
            |b| {
                b.iter(|| c1 = c1.add(&c2).unwrap());
            },
        );

        group.bench_function(
            BenchmarkId::new("add_assign_ct", format!("n={}/log(q)={}", par.degree(), q)),
            |b| {
                b.iter(|| c1.add_assign(&c2).unwrap());
            },
        );

        group.bench_function(
            BenchmarkId::new("add_pt", format!("n={}/log(q)={}", par.degree(), q)),
            |b| {
                b.iter(|| c1 = c1.add_plaintext(&pt2).unwrap());
            },
        );

        group.bench_function(
            BenchmarkId::new("add_assign_pt", format!("n={}/log(q)={}", par.degree(), q)),
            |b| {
                b.iter(|| c1.add_plaintext_assign(&pt2).unwrap());
            },
        );

        group.bench_function(
            BenchmarkId::new("sub_ct", format!("n={}/log(q)={}", par.degree(), q)),
            |b| {
                b.iter(|| c1 = c1.subtract(&c2).unwrap());
            },
        );

        group.bench_function(
            BenchmarkId::new("sub_assign_ct", format!("n={}/log(q)={}", par.degree(), q)),
            |b| {
                b.iter(|| c1.subtract_assign(&c2).unwrap());
            },
        );

        group.bench_function(
            BenchmarkId::new("sub_pt", format!("n={}/log(q)={}", par.degree(), q)),
            |b| {
                b.iter(|| c1 = c1.subtract_plaintext(&pt2).unwrap());
            },
        );

        group.bench_function(
            BenchmarkId::new("sub_assign_pt", format!("n={}/log(q)={}", par.degree(), q)),
            |b| {
                b.iter(|| c1.subtract_plaintext_assign(&pt2).unwrap());
            },
        );

        group.bench_function(
            BenchmarkId::new("neg", format!("n={}/log(q)={}", par.degree(), q)),
            |b| {
                b.iter(|| c1 = -&c2);
            },
        );

        let mut c3 = c1.multiply(&c1).unwrap();
        let c3_clone = c3.clone();
        if let Some(rk) = rk.as_ref() {
            group.bench_function(
                BenchmarkId::new("relinearize", format!("n={}/log(q)={}", par.degree(), q)),
                |b| {
                    b.iter(|| {
                        assert!(rk.relinearize(&mut c3).is_ok());
                        c3 = c3_clone.clone();
                    });
                },
            );
        }

        if let Some(ek) = ek {
            group.bench_function(
                BenchmarkId::new("rotate_rows", format!("n={}/log(q)={}", par.degree(), q)),
                |b| {
                    b.iter(|| c1 = ek.rotate_rows(&c1).unwrap());
                },
            );

            group.bench_function(
                BenchmarkId::new("rotate_columns", format!("n={}/log(q)={}", par.degree(), q)),
                |b| {
                    b.iter(|| c1 = ek.rotate_columns(&c1, 1).unwrap());
                },
            );

            group.bench_function(
                BenchmarkId::new("inner_sum", format!("n={}/log(q)={}", par.degree(), q)),
                |b| {
                    b.iter(|| c1 = ek.inner_sum(&c1).unwrap());
                },
            );

            for i in 1..=par.degree().ilog2() {
                if par.degree() > 2048 && i > 4 {
                    continue; // Skip slow benchmarks
                }
                group.bench_function(
                    BenchmarkId::new(
                        format!("expand_{i}"),
                        format!("n={}/log(q)={}", par.degree(), q),
                    ),
                    |b| {
                        b.iter(|| ek.expand(&c1, 1 << i).unwrap());
                    },
                );
            }
        }

        group.bench_function(
            BenchmarkId::new("mul", format!("n={}/log(q)={}", par.degree(), q)),
            |b| {
                b.iter(|| c1.multiply(&c2).unwrap());
            },
        );

        group.bench_function(
            BenchmarkId::new("square", format!("n={}/log(q)={}", par.degree(), q)),
            |b| {
                b.iter(|| c1.multiply(&c1).unwrap());
            },
        );

        if let Some(rk) = rk.as_ref() {
            group.bench_function(
                BenchmarkId::new(
                    "mul_then_relinearize",
                    format!("n={}/log(q)={}", par.degree(), q),
                ),
                |b| {
                    b.iter(|| {
                        c3 = c1.multiply(&c2).unwrap();
                        assert!(rk.relinearize(&mut c3).is_ok());
                    });
                },
            );

            // Default multiplication method
            let multiplicator = MultiplicationPlan::with_relinearization(rk).unwrap();

            group.bench_function(
                BenchmarkId::new("mul_and_relin", format!("n={}/log(q)={}", par.degree(), q)),
                |b| {
                    b.iter(|| assert!(multiplicator.multiply(&c1, &c2).is_ok()));
                },
            );

            // Second multiplication option.
            let nmoduli = q.div_ceil(62);
            let mut extended_basis = par.moduli().to_vec();
            let mut upper_bound = u64::MAX >> 2;
            while extended_basis.len() != nmoduli + par.moduli().len() {
                upper_bound = generate_prime(62, 2 * par.degree() as u64, upper_bound).unwrap();
                if !extended_basis.contains(&upper_bound) {
                    extended_basis.push(upper_bound)
                }
            }
            let rns_q = RnsContext::new(&extended_basis[..par.moduli().len()]).unwrap();
            let rns_p = RnsContext::new(&extended_basis[par.moduli().len()..]).unwrap();
            let multiplicator = MultiplicationPlan::builder(&par)
                .extended_basis(&extended_basis)
                .scaling(fhe::bfv::evaluation::MultiplicationScaling {
                    left: ScalingFactor::one(),
                    right: ScalingFactor::new(rns_p.modulus(), rns_q.modulus()),
                    product: ScalingFactor::new(
                        &BigUint::from(par.plaintext_modulus_u64().unwrap()),
                        rns_p.modulus(),
                    ),
                })
                .relinearization(rk)
                .build()
                .unwrap();
            group.bench_function(
                BenchmarkId::new(
                    "mul_and_relin_2",
                    format!("n={}/log(q)={}", par.degree(), q),
                ),
                |b| {
                    b.iter(|| assert!(multiplicator.multiply(&c1, &c2).is_ok()));
                },
            );
        }
    }

    group.finish();
}

criterion_group!(bfv, bfv_benchmark);
criterion_main!(bfv);
