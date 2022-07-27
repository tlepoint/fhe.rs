use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use crrust_math::rq::*;
use rand::RngCore;
use std::rc::Rc;

fn random_vector(size: usize, p: u64) -> Vec<u64> {
	let mut rng = rand::thread_rng();
	let mut v = vec![];
	for _ in 0..size {
		v.push(rng.next_u64() % p)
	}
	v
}

static MODULI: &[u64; 4] = &[
	4611686018326724609,
	4611686018309947393,
	4611686018282684417,
	4611686018257518593,
];

pub fn rq_benchmark(c: &mut Criterion) {
	let mut group = c.benchmark_group("rq");
	group.sample_size(50);

	for degree in &[1024usize, 4096] {
		let ctx1 = Rc::new(Context::new(&MODULI[0..1], *degree).unwrap());
		let ctx2 = Rc::new(Context::new(&MODULI[0..2], *degree).unwrap());
		let ctx3 = Rc::new(Context::new(&MODULI[0..3], *degree).unwrap());
		let ctx4 = Rc::new(Context::new(MODULI, *degree).unwrap());

		let mut p1 = <Poly as TryFrom<&[u64]>>::try_from(
			&random_vector(*degree, MODULI[0]),
			&ctx1,
			Representation::Ntt,
		)
		.ok()
		.unwrap();
		let q1 = p1.clone();

		let mut a = random_vector(*degree, MODULI[0]);
		a.append(&mut random_vector(*degree, MODULI[1]));
		let mut p2 = <Poly as TryFrom<&[u64]>>::try_from(&a, &ctx2, Representation::Ntt)
			.ok()
			.unwrap();
		let q2 = p2.clone();

		let mut a = random_vector(*degree, MODULI[0]);
		a.append(&mut random_vector(*degree, MODULI[1]));
		a.append(&mut random_vector(*degree, MODULI[2]));
		let mut p3 = <Poly as TryFrom<&[u64]>>::try_from(&a, &ctx3, Representation::Ntt)
			.ok()
			.unwrap();
		let q3 = p3.clone();

		let mut a = random_vector(*degree, MODULI[0]);
		a.append(&mut random_vector(*degree, MODULI[1]));
		a.append(&mut random_vector(*degree, MODULI[2]));
		a.append(&mut random_vector(*degree, MODULI[3]));
		let mut p4 = <Poly as TryFrom<&[u64]>>::try_from(&a, &ctx4, Representation::Ntt)
			.ok()
			.unwrap();
		let q4 = p4.clone();

		group.bench_function(BenchmarkId::new("add", format!("{}/{}", degree, 62)), |b| {
			b.iter(|| p1 += &q1);
		});

		group.bench_function(
			BenchmarkId::new("add", format!("{}/{}", degree, 2 * 62)),
			|b| {
				b.iter(|| p2 += &q2);
			},
		);

		group.bench_function(
			BenchmarkId::new("add", format!("{}/{}", degree, 3 * 62)),
			|b| {
				b.iter(|| p3 += &q3);
			},
		);

		group.bench_function(
			BenchmarkId::new("add", format!("{}/{}", degree, 4 * 62)),
			|b| {
				b.iter(|| p4 += &q4);
			},
		);
	}

	group.finish();
}

criterion_group!(rq, rq_benchmark);
criterion_main!(rq);
