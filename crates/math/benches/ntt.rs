use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use math::zq::{ntt::NttOperator, Modulus};
use std::rc::Rc;

pub fn ntt_benchmark(c: &mut Criterion) {
	let mut group = c.benchmark_group("ntt");
	group.sample_size(50);

	for vector_size in [1024usize, 4096].iter() {
		for p in [4611686018326724609u64, 40961u64] {
			let p_nbits = 64 - p.leading_zeros();
			let q = Modulus::new(p).unwrap();
			let mut a = q.random_vec(*vector_size);
			let op = NttOperator::new(&Rc::new(q), *vector_size).unwrap();

			group.bench_function(
				BenchmarkId::new("forward", format!("{}/{}", vector_size, p_nbits)),
				|b| b.iter(|| op.forward(&mut a)),
			);

			group.bench_function(
				BenchmarkId::new("forward_vt", format!("{}/{}", vector_size, p_nbits)),
				|b| b.iter(|| unsafe { op.forward_vt(&mut a) }),
			);

			group.bench_function(
				BenchmarkId::new("backward", format!("{}/{}", vector_size, p_nbits)),
				|b| b.iter(|| op.backward(&mut a)),
			);

			group.bench_function(
				BenchmarkId::new("backward_vt", format!("{}/{}", vector_size, p_nbits)),
				|b| b.iter(|| unsafe { op.backward_vt(&mut a) }),
			);
		}
	}

	group.finish();
}

criterion_group!(ntt, ntt_benchmark);
criterion_main!(ntt);
