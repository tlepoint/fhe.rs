use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use math::rns::{RnsContext, RnsConverter, RnsScaler};
use num_bigint::BigUint;
use rand::{thread_rng, RngCore};

pub fn rns_benchmark(c: &mut Criterion) {
	let mut group = c.benchmark_group("rns");
	group.sample_size(50);

	let q = [
		4611686018326724609u64,
		4611686018309947393,
		4611686018282684417,
	];
	let p = [
		4611686018257518593u64,
		4611686018232352769,
		4611686018171535361,
		4611686018106523649,
	];

	let mut rng = thread_rng();
	let mut x = vec![];
	for qi in &q {
		x.push(rng.next_u64() % *qi);
	}

	let rns_q = RnsContext::new(&q).unwrap();
	let rns_p = RnsContext::new(&p).unwrap();
	let converter = RnsConverter::new(&rns_q, &rns_p);
	let scaler = RnsScaler::new(
		&rns_q,
		&BigUint::from(1u64),
		&BigUint::from(46116860181065u64),
	);

	group.bench_function(
		BenchmarkId::new("converter", format!("{}->{}", q.len(), p.len())),
		|b| {
			b.iter(|| converter.convert(&x));
		},
	);

	group.bench_function(BenchmarkId::new("scaler", q.len()), |b| {
		b.iter(|| scaler.scale(&x, x.len(), true));
	});

	group.finish();
}

criterion_group!(rns, rns_benchmark);
criterion_main!(rns);
