use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use crrust_math::zq::{ntt::NttOperator, Modulus};
use rand::RngCore;
use std::{rc::Rc, vec};

fn random_vector(size: usize, p: u64) -> Vec<u64> {
	let mut rng = rand::thread_rng();
	let mut v = vec![];
	for _ in 0..size {
		v.push(rng.next_u64() % p)
	}
	v
}

pub fn ntt_benchmark(c: &mut Criterion) {
	let mut group = c.benchmark_group("ntt");
	group.sample_size(50);

	let p = 4611686018326724609;

	for vector_size in [1024usize, 4096].iter() {
		let mut a = random_vector(*vector_size, p);
		let q = Modulus::new(p).unwrap();
		let op = NttOperator::new(&Rc::new(q), *vector_size).unwrap();

		group.bench_function(BenchmarkId::new("forward", vector_size), |b| {
			b.iter(|| op.forward(&mut a));
		});

		group.bench_function(BenchmarkId::new("backward", vector_size), |b| {
			b.iter(|| op.forward(&mut a));
		});
	}

	group.finish();
}

criterion_group!(ntt, ntt_benchmark);
criterion_main!(ntt);
