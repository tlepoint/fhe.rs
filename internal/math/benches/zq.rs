use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use math::zq::Modulus;

pub fn zq_benchmark(c: &mut Criterion) {
	let mut group = c.benchmark_group("zq");
	group.sample_size(50);

	let p = 4611686018326724609;

	for vector_size in [1024usize, 4096].iter() {
		let q = Modulus::new(p).unwrap();
		let mut a = q.random_vec(*vector_size);
		let c = q.random_vec(*vector_size);
		let c_shoup = q.shoup_vec(&c);
		let scalar = c[0];

		group.bench_function(BenchmarkId::new("add_vec", vector_size), |b| {
			b.iter(|| q.add_vec(&mut a, &c));
		});

		group.bench_function(BenchmarkId::new("sub_vec", vector_size), |b| {
			b.iter(|| q.sub_vec(&mut a, &c));
		});

		group.bench_function(BenchmarkId::new("neg_vec", vector_size), |b| {
			b.iter(|| q.neg_vec(&mut a));
		});

		group.bench_function(BenchmarkId::new("mul_vec", vector_size), |b| {
			b.iter(|| q.mul_vec(&mut a, &c));
		});

		group.bench_function(BenchmarkId::new("mul_shoup_vec", vector_size), |b| {
			b.iter(|| q.mul_shoup_vec(&mut a, &c, &c_shoup));
		});

		group.bench_function(BenchmarkId::new("scalar_mul_vec", vector_size), |b| {
			b.iter(|| q.scalar_mul_vec(&mut a, scalar));
		});

		group.bench_function(BenchmarkId::new("vt_mul_vec", vector_size), |b| unsafe {
			b.iter(|| q.mul_vec_vt(&mut a, &c));
		});
	}

	group.finish();
}

criterion_group!(zq, zq_benchmark);
criterion_main!(zq);
