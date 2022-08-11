use bfv::{
	mul, mul2,
	traits::{Encoder, Encryptor},
	BfvParameters, BfvParametersBuilder, Encoding, GaloisKey, InnerSumKey, Plaintext,
	RelinearizationKey, SecretKey,
};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use itertools::Itertools;
use math::rq::{Context, Poly, Representation};
use std::rc::Rc;

fn params() -> Vec<Rc<BfvParameters>> {
	let par_62 = BfvParametersBuilder::default()
		.polynomial_degree(16384)
		.plaintext_modulus(1153)
		.ciphertext_moduli_sizes(vec![62, 62, 62, 62, 62, 62, 62])
		.build()
		.unwrap();
	let par_50 = BfvParametersBuilder::default()
		.polynomial_degree(16384)
		.plaintext_modulus(1153)
		.ciphertext_moduli_sizes(vec![50, 50, 50, 50, 50, 50, 50])
		.build()
		.unwrap();
	vec![Rc::new(par_62), Rc::new(par_50)]
}

pub fn ops_benchmark(c: &mut Criterion) {
	let mut group = c.benchmark_group("ops");

	for par in params() {
		let sk = SecretKey::random(&par);

		let pt1 = Plaintext::try_encode(&(1..16u64).collect_vec() as &[u64], Encoding::Poly, &par)
			.unwrap();
		let pt2 = Plaintext::try_encode(&(3..39u64).collect_vec() as &[u64], Encoding::Poly, &par)
			.unwrap();
		let mut c1 = sk.encrypt(&pt1).unwrap();
		let c2 = sk.encrypt(&pt2).unwrap();

		group.bench_function(
			BenchmarkId::new(
				"add",
				format!(
					"{}/{}",
					par.degree(),
					par.moduli_sizes().iter().sum::<usize>()
				),
			),
			|b| {
				b.iter(|| c1 += &c2);
			},
		);

		group.bench_function(
			BenchmarkId::new(
				"sub",
				format!(
					"{}/{}",
					par.degree(),
					par.moduli_sizes().iter().sum::<usize>()
				),
			),
			|b| {
				b.iter(|| c1 -= &c2);
			},
		);

		group.bench_function(
			BenchmarkId::new(
				"neg",
				format!(
					"{}/{}",
					par.degree(),
					par.moduli_sizes().iter().sum::<usize>()
				),
			),
			|b| {
				b.iter(|| c1 = -&c2);
			},
		);

		let isk = InnerSumKey::new(&sk).unwrap();
		let rk = RelinearizationKey::new(&sk).unwrap();
		let gk = GaloisKey::new(&sk, par.degree() + 1).unwrap();
		let ctx = Rc::new(Context::new(par.moduli(), par.degree()).unwrap());
		let p3 = Poly::random(&ctx, Representation::PowerBasis);
		let mut p2 = Poly::random(&ctx, Representation::Ntt);
		let mut p1 = Poly::random(&ctx, Representation::Ntt);

		group.bench_function(
			BenchmarkId::new(
				"relinearize",
				format!(
					"{}/{}",
					par.degree(),
					par.moduli_sizes().iter().sum::<usize>()
				),
			),
			|b| {
				b.iter(|| rk.relinearize(&mut p1, &mut p2, &p3));
			},
		);

		group.bench_function(
			BenchmarkId::new(
				"galois",
				format!(
					"{}/{}",
					par.degree(),
					par.moduli_sizes().iter().sum::<usize>()
				),
			),
			|b| {
				b.iter(|| gk.relinearize(&mut c1));
			},
		);

		group.bench_function(
			BenchmarkId::new(
				"inner_sum",
				format!(
					"{}/{}",
					par.degree(),
					par.moduli_sizes().iter().sum::<usize>()
				),
			),
			|b| {
				b.iter(|| c1 = isk.inner_sum(&c1).unwrap());
			},
		);

		group.bench_function(
			BenchmarkId::new(
				"mul",
				format!(
					"{}/{}",
					par.degree(),
					par.moduli_sizes().iter().sum::<usize>()
				),
			),
			|b| {
				b.iter(|| mul(&c1, &c2, &rk));
			},
		);

		group.bench_function(
			BenchmarkId::new(
				"mul2",
				format!(
					"{}/{}",
					par.degree(),
					par.moduli_sizes().iter().sum::<usize>()
				),
			),
			|b| {
				b.iter(|| mul2(&c1, &c2, &rk));
			},
		);
	}

	group.finish();
}

criterion_group!(ops, ops_benchmark);
criterion_main!(ops);
