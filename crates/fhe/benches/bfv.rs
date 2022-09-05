#![feature(int_log)]
#![feature(int_roundings)]

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use fhe::{
	bfv::{
		BfvParameters, Ciphertext, Encoding, EvaluationKeyBuilder, Multiplicator, Plaintext,
		PublicKey, RelinearizationKey, SecretKey,
	},
	Result,
};
use fhe_traits::{FheEncoder, FheEncrypter};
use itertools::Itertools;
use math::rns::{RnsContext, ScalingFactor};
use math::zq::primes::generate_prime;
use num_bigint::BigUint;
use std::time::Duration;

pub fn bfv_benchmark(c: &mut Criterion) {
	let mut group = c.benchmark_group("bfv");
	group.sample_size(10);
	group.warm_up_time(Duration::from_millis(600));
	group.measurement_time(Duration::from_millis(1000));

	for par in BfvParameters::default_parameters_128(20) {
		let sk = SecretKey::random(&par);

		let ek = if par.moduli().len() > 1 {
			Some(
				EvaluationKeyBuilder::new(&sk)
					.unwrap()
					.enable_inner_sum()
					.unwrap()
					.enable_column_rotation(1)
					.unwrap()
					.enable_expansion(par.degree().ilog2() as usize)
					.unwrap()
					.build()
					.unwrap(),
			)
		} else {
			None
		};

		let rk = if par.moduli().len() > 1 {
			Some(RelinearizationKey::new(&sk).unwrap())
		} else {
			None
		};

		let pt1 =
			Plaintext::try_encode(&(1..16u64).collect_vec() as &[u64], Encoding::simd(), &par)
				.unwrap();
		let pt2 =
			Plaintext::try_encode(&(3..39u64).collect_vec() as &[u64], Encoding::simd(), &par)
				.unwrap();
		let mut c1: Ciphertext = sk.try_encrypt(&pt1).unwrap();
		let c2: Ciphertext = sk.try_encrypt(&pt2).unwrap();

		let q = par.moduli_sizes().iter().sum::<usize>();

		group.bench_function(
			BenchmarkId::new("keygen_sk", format!("n={}/log(q)={}", par.degree(), q)),
			|b| {
				b.iter(|| SecretKey::random(&par));
			},
		);

		group.bench_function(
			BenchmarkId::new("keygen_pk", format!("n={}/log(q)={}", par.degree(), q)),
			|b| {
				b.iter(|| PublicKey::new(&sk));
			},
		);

		group.bench_function(
			BenchmarkId::new("keygen_rk", format!("n={}/log(q)={}", par.degree(), q)),
			|b| {
				b.iter(|| RelinearizationKey::new(&sk));
			},
		);

		group.bench_function(
			BenchmarkId::new("encode_poly", format!("n={}/log(q)={}", par.degree(), q)),
			|b| {
				b.iter(|| {
					Plaintext::try_encode(
						&(1..16u64).collect_vec() as &[u64],
						Encoding::poly(),
						&par,
					)
				});
			},
		);

		group.bench_function(
			BenchmarkId::new("encode_simd", format!("n={}/log(q)={}", par.degree(), q)),
			|b| {
				b.iter(|| {
					Plaintext::try_encode(
						&(1..16u64).collect_vec() as &[u64],
						Encoding::simd(),
						&par,
					)
				});
			},
		);

		group.bench_function(
			BenchmarkId::new("encrypt_sk", format!("n={}/log(q)={}", par.degree(), q)),
			|b| {
				b.iter(|| {
					let _: Result<Ciphertext> = sk.try_encrypt(&pt1);
				});
			},
		);

		group.bench_function(
			BenchmarkId::new("add_ct", format!("n={}/log(q)={}", par.degree(), q)),
			|b| {
				b.iter(|| c1 = &c1 + &c2);
			},
		);

		group.bench_function(
			BenchmarkId::new("add_assign_ct", format!("n={}/log(q)={}", par.degree(), q)),
			|b| {
				b.iter(|| c1 += &c2);
			},
		);

		group.bench_function(
			BenchmarkId::new("add_pt", format!("n={}/log(q)={}", par.degree(), q)),
			|b| {
				b.iter(|| c1 = &c1 + &pt2);
			},
		);

		group.bench_function(
			BenchmarkId::new("add_assign_pt", format!("n={}/log(q)={}", par.degree(), q)),
			|b| {
				b.iter(|| c1 += &pt2);
			},
		);

		group.bench_function(
			BenchmarkId::new("sub_ct", format!("n={}/log(q)={}", par.degree(), q)),
			|b| {
				b.iter(|| c1 = &c1 - &c2);
			},
		);

		group.bench_function(
			BenchmarkId::new("sub_assign_ct", format!("n={}/log(q)={}", par.degree(), q)),
			|b| {
				b.iter(|| c1 -= &c2);
			},
		);

		group.bench_function(
			BenchmarkId::new("sub_pt", format!("n={}/log(q)={}", par.degree(), q)),
			|b| {
				b.iter(|| c1 = &c1 - &pt2);
			},
		);

		group.bench_function(
			BenchmarkId::new("sub_assign_pt", format!("n={}/log(q)={}", par.degree(), q)),
			|b| {
				b.iter(|| c1 -= &pt2);
			},
		);

		group.bench_function(
			BenchmarkId::new("neg", format!("n={}/log(q)={}", par.degree(), q)),
			|b| {
				b.iter(|| c1 = -&c2);
			},
		);

		let mut c3 = &c1 * &c1;
		let c3_clone = c3.clone();
		if let Some(rk) = rk.as_ref() {
			group.bench_function(
				BenchmarkId::new("relinearize", format!("n={}/log(q)={}", par.degree(), q)),
				|b| {
					b.iter(|| {
						assert!(rk.relinearizes(&mut c3).is_ok());
						c3 = c3_clone.clone();
					});
				},
			);
		}

		if let Some(ek) = ek {
			group.bench_function(
				BenchmarkId::new("rotate_rows", format!("n={}/log(q)={}", par.degree(), q)),
				|b| {
					b.iter(|| c1 = ek.rotates_rows(&c1).unwrap());
				},
			);

			group.bench_function(
				BenchmarkId::new("rotate_columns", format!("n={}/log(q)={}", par.degree(), q)),
				|b| {
					b.iter(|| c1 = ek.rotates_columns_by(&c1, 1).unwrap());
				},
			);

			group.bench_function(
				BenchmarkId::new("inner_sum", format!("n={}/log(q)={}", par.degree(), q)),
				|b| {
					b.iter(|| c1 = ek.computes_inner_sum(&c1).unwrap());
				},
			);

			for i in 1..par.degree().ilog2() + 1 {
				if par.degree() > 2048 && i > 4 {
					continue; // Skip slow benchmarks
				}
				group.bench_function(
					BenchmarkId::new(
						format!("expand_{}", i),
						format!("n={}/log(q)={}", par.degree(), q),
					),
					|b| {
						b.iter(|| ek.expands(&c1, 1 << i).unwrap());
					},
				);
			}
		}

		group.bench_function(
			BenchmarkId::new("mul", format!("n={}/log(q)={}", par.degree(), q)),
			|b| {
				b.iter(|| &c1 * &c2);
			},
		);

		group.bench_function(
			BenchmarkId::new("square", format!("n={}/log(q)={}", par.degree(), q)),
			|b| {
				b.iter(|| &c1 * &c1);
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
						c3 = &c1 * &c2;
						assert!(rk.relinearizes(&mut c3).is_ok());
					});
				},
			);

			// Default multiplication method
			let multiplicator = Multiplicator::default(rk).unwrap();

			group.bench_function(
				BenchmarkId::new("mul_and_relin", format!("n={}/log(q)={}", par.degree(), q)),
				|b| {
					b.iter(|| assert!(multiplicator.multiply(&c1, &c2).is_ok()));
				},
			);

			// Second multiplication option.
			let nmoduli = (q).div_ceil(62);
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
			let mut multiplicator = Multiplicator::new(
				ScalingFactor::one(),
				ScalingFactor::new(rns_p.modulus(), rns_q.modulus()),
				&extended_basis,
				ScalingFactor::new(&BigUint::from(par.plaintext()), rns_p.modulus()),
				&par,
			)
			.unwrap();
			assert!(multiplicator.enable_relinearization(rk).is_ok());
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
