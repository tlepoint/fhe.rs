use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use crrust_math::rq::*;
use rand::RngCore;
use std::rc::Rc;
use std::time::Duration;

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
	group.warm_up_time(Duration::from_secs(2));
	group.measurement_time(Duration::from_secs(2));

	for degree in &[1024usize, 4096] {
		for nmoduli in 1..=MODULI.len() {
			if !nmoduli.is_power_of_two() {
				continue;
			}
			let ctx = Rc::new(Context::new(&MODULI[0..nmoduli], *degree).unwrap());

			let mut a = random_vector(*degree, MODULI[0]);
			for i in 1..nmoduli {
				a.append(&mut random_vector(*degree, MODULI[i]));
			}
			let mut p = Poly::try_convert_from(&a, &ctx, Representation::Ntt)
				.ok()
				.unwrap();
			let mut q = p.clone();

			group.bench_function(
				BenchmarkId::new("add", format!("{}/{}", degree, 62 * nmoduli)),
				|b| {
					b.iter(|| p += &q);
				},
			);

			group.bench_function(
				BenchmarkId::new("sub", format!("{}/{}", degree, 62 * nmoduli)),
				|b| {
					b.iter(|| p -= &q);
				},
			);

			group.bench_function(
				BenchmarkId::new("mul", format!("{}/{}", degree, 62 * nmoduli)),
				|b| {
					b.iter(|| p *= &q);
				},
			);

			q.change_representation(Representation::NttShoup);

			group.bench_function(
				BenchmarkId::new("mul_shoup", format!("{}/{}", degree, 62 * nmoduli)),
				|b| {
					b.iter(|| p *= &q);
				},
			);

			group.bench_function(
				BenchmarkId::new(
					"change_representation/PowerBasis_to_Ntt",
					format!("{}/{}", degree, 62 * nmoduli),
				),
				|b| {
					b.iter(|| {
						unsafe {
							p.override_representation(Representation::PowerBasis);
						}
						p.change_representation(Representation::Ntt)
					});
				},
			);

			group.bench_function(
				BenchmarkId::new(
					"change_representation/Ntt_to_PowerBasis",
					format!("{}/{}", degree, 62 * nmoduli),
				),
				|b| {
					b.iter(|| {
						unsafe {
							p.override_representation(Representation::Ntt);
						}
						p.change_representation(Representation::PowerBasis)
					});
				},
			);
		}
	}

	group.finish();
}

criterion_group!(rq, rq_benchmark);
criterion_main!(rq);
