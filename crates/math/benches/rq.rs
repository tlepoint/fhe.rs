use criterion::measurement::WallTime;
use criterion::{criterion_group, criterion_main, BenchmarkGroup, BenchmarkId, Criterion};
use itertools::{izip, Itertools};
use math::rq::*;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use std::sync::Arc;
use std::time::Duration;

static MODULI: &[u64; 4] = &[
	562949954093057,
	4611686018326724609,
	4611686018309947393,
	4611686018282684417,
];

static DEGREE: &[usize] = &[1024, 2048, 4096, 8192];

fn create_group<'a>(c: &'a mut Criterion, name: String) -> BenchmarkGroup<'a, WallTime> {
	let mut group = c.benchmark_group(name);
	group.warm_up_time(Duration::from_millis(100));
	group.measurement_time(Duration::from_secs(1));
	group
}

macro_rules! bench_op {
	($c: expr, $name:expr, $op:expr, $vt:expr) => {{
		let name = if $vt {
			format!("{}_vt", $name)
		} else {
			$name.to_string()
		};
		let mut group = create_group($c, name);

		for degree in DEGREE {
			let ctx = Arc::new(Context::new(&MODULI[..1], *degree).unwrap());
			let p = Poly::random(&ctx, Representation::Ntt);
			let mut q = Poly::random(&ctx, Representation::Ntt);
			if $vt {
				unsafe { q.allow_variable_time_computations() }
			}

			group.bench_function(
				BenchmarkId::from_parameter(format!("{}/{}", degree, ctx.modulus().bits())),
				|b| {
					b.iter(|| $op(&p, &q));
				},
			);
		}
	}};
}

macro_rules! bench_op_unary {
	($c: expr, $name:expr, $op:expr, $vt:expr) => {{
		let name = if $vt {
			format!("{}_vt", $name)
		} else {
			$name.to_string()
		};
		let mut group = create_group($c, name);

		for degree in DEGREE {
			let ctx = Arc::new(Context::new(&MODULI[..1], *degree).unwrap());
			let p = Poly::random(&ctx, Representation::Ntt);
			let mut q = Poly::random(&ctx, Representation::Ntt);
			if $vt {
				unsafe { q.allow_variable_time_computations() }
			}

			group.bench_function(
				BenchmarkId::from_parameter(format!("{}/{}", degree, ctx.modulus().bits())),
				|b| {
					b.iter(|| $op(&p));
				},
			);
		}
	}};
}

macro_rules! bench_op_assign {
	($c: expr, $name:expr, $op:expr, $vt:expr) => {{
		let name = if $vt {
			format!("{}_vt", $name)
		} else {
			$name.to_string()
		};
		let mut group = create_group($c, name);

		for degree in DEGREE {
			let ctx = Arc::new(Context::new(&MODULI[..1], *degree).unwrap());
			let mut p = Poly::random(&ctx, Representation::Ntt);
			let mut q = Poly::random(&ctx, Representation::Ntt);
			if $vt {
				unsafe { q.allow_variable_time_computations() }
			}

			group.bench_function(
				BenchmarkId::from_parameter(format!("{}/{}", degree, ctx.modulus().bits())),
				|b| {
					b.iter(|| $op(&mut p, &q));
				},
			);
		}
	}};
}

pub fn rq_op_benchmark(c: &mut Criterion) {
	for vt in [false, true] {
		bench_op!(c, "rq_add", <&Poly>::add, vt);
		bench_op_assign!(c, "rq_add_assign", Poly::add_assign, vt);
		bench_op!(c, "rq_sub", <&Poly>::sub, vt);
		bench_op_assign!(c, "rq_sub_assign", Poly::sub_assign, vt);
		bench_op!(c, "rq_mul", <&Poly>::mul, vt);
		bench_op_assign!(c, "rq_mul_assign", Poly::mul_assign, vt);
		bench_op_unary!(c, "rq_neg", <&Poly>::neg, vt);
	}
}

pub fn rq_dot_product(c: &mut Criterion) {
	let mut group = create_group(c, "rq_dot_product".to_string());
	for degree in DEGREE {
		let ctx = Arc::new(Context::new(&MODULI[..1], *degree).unwrap());
		let p_vec = (0..128)
			.map(|_| Poly::random(&ctx, Representation::Ntt))
			.collect_vec();
		let mut q_vec = (0..128)
			.map(|_| Poly::random(&ctx, Representation::Ntt))
			.collect_vec();
		let mut out = Poly::zero(&ctx, Representation::Ntt);

		group.bench_function(
			BenchmarkId::from_parameter(format!("naive/{}/{}", degree, ctx.modulus().bits())),
			|b| {
				b.iter(|| izip!(p_vec.iter(), q_vec.iter()).for_each(|(pi, qi)| out += pi * qi));
			},
		);

		q_vec
			.iter_mut()
			.for_each(|qi| qi.change_representation(Representation::NttShoup));
		group.bench_function(
			BenchmarkId::from_parameter(format!("naive_shoup/{}/{}", degree, ctx.modulus().bits())),
			|b| {
				b.iter(|| izip!(p_vec.iter(), q_vec.iter()).for_each(|(pi, qi)| out += pi * qi));
			},
		);

		q_vec
			.iter_mut()
			.for_each(|qi| qi.change_representation(Representation::Ntt));
		group.bench_function(
			BenchmarkId::from_parameter(format!("opt/{}/{}", degree, ctx.modulus().bits())),
			|b| {
				b.iter(|| dot_product(p_vec.iter(), q_vec.iter()));
			},
		);

		let ctx = Arc::new(Context::new(&MODULI[1..2], *degree).unwrap());
		let p_vec = (0..128)
			.map(|_| Poly::random(&ctx, Representation::Ntt))
			.collect_vec();
		let q_vec = (0..128)
			.map(|_| Poly::random(&ctx, Representation::Ntt))
			.collect_vec();

		group.bench_function(
			BenchmarkId::from_parameter(format!("opt/{}/{}", degree, ctx.modulus().bits())),
			|b| {
				b.iter(|| dot_product(p_vec.iter(), q_vec.iter()));
			},
		);
	}
}

pub fn rq_benchmark(c: &mut Criterion) {
	let mut group = c.benchmark_group("rq");
	group.warm_up_time(Duration::from_millis(100));
	group.measurement_time(Duration::from_secs(1));

	for degree in &[1024usize, 4096] {
		for nmoduli in 1..=MODULI.len() {
			if !nmoduli.is_power_of_two() {
				continue;
			}
			let ctx = Arc::new(Context::new(&MODULI[..nmoduli], *degree).unwrap());
			let mut p = Poly::random(&ctx, Representation::Ntt);
			let mut q = Poly::random(&ctx, Representation::Ntt);
			let p_vec = (0..128)
				.map(|_| Poly::random(&ctx, Representation::Ntt))
				.collect_vec();
			let mut q_vec = (0..128)
				.map(|_| Poly::random(&ctx, Representation::Ntt))
				.collect_vec();

			group.bench_function(
				BenchmarkId::new(
					"dot_product/128/naive",
					format!("{}/{}", degree, ctx.modulus().bits()),
				),
				|b| {
					b.iter(|| izip!(&p_vec, &q_vec).for_each(|(pi, qi)| p += pi * qi));
				},
			);

			group.bench_function(
				BenchmarkId::new(
					"dot_product/128/opt",
					format!("{}/{}", degree, ctx.modulus().bits()),
				),
				|b| {
					b.iter(|| dot_product(p_vec.iter(), q_vec.iter()));
				},
			);

			q.change_representation(Representation::NttShoup);

			q_vec
				.iter_mut()
				.for_each(|qi| qi.change_representation(Representation::NttShoup));
			group.bench_function(
				BenchmarkId::new(
					"dot_product/128/shoup",
					format!("{}/{}", degree, ctx.modulus().bits()),
				),
				|b| {
					b.iter(|| izip!(&p_vec, &q_vec).for_each(|(pi, qi)| p += pi * qi));
				},
			);

			group.bench_function(
				BenchmarkId::new("mul_shoup", format!("{}/{}", degree, ctx.modulus().bits())),
				|b| {
					b.iter(|| p = &p * &q);
				},
			);

			group.bench_function(
				BenchmarkId::new(
					"mul_shoup_assign",
					format!("{}/{}", degree, ctx.modulus().bits()),
				),
				|b| {
					b.iter(|| p *= &q);
				},
			);

			group.bench_function(
				BenchmarkId::new(
					"change_representation/PowerBasis_to_Ntt",
					format!("{}/{}", degree, ctx.modulus().bits()),
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
					format!("{}/{}", degree, ctx.modulus().bits()),
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

			p.change_representation(Representation::Ntt);
			q.change_representation(Representation::Ntt);

			unsafe {
				q.allow_variable_time_computations();
				q_vec.iter_mut().for_each(|qi| {
					qi.change_representation(Representation::Ntt);
					qi.allow_variable_time_computations()
				});

				group.bench_function(
					BenchmarkId::new("add_vt", format!("{}/{}", degree, ctx.modulus().bits())),
					|b| {
						b.iter(|| p += &q);
					},
				);

				group.bench_function(
					BenchmarkId::new("sub_vt", format!("{}/{}", degree, ctx.modulus().bits())),
					|b| {
						b.iter(|| p -= &q);
					},
				);

				group.bench_function(
					BenchmarkId::new("mul_vt", format!("{}/{}", degree, ctx.modulus().bits())),
					|b| {
						b.iter(|| p *= &q);
					},
				);

				group.bench_function(
					BenchmarkId::new(
						"dot_product/128/naive_vt",
						format!("{}/{}", degree, ctx.modulus().bits()),
					),
					|b| {
						b.iter(|| izip!(&p_vec, &q_vec).for_each(|(pi, qi)| p += pi * qi));
					},
				);

				q.change_representation(Representation::NttShoup);
				q_vec
					.iter_mut()
					.for_each(|qi| qi.change_representation(Representation::NttShoup));

				group.bench_function(
					BenchmarkId::new(
						"mul_shoup_vt",
						format!("{}/{}", degree, ctx.modulus().bits()),
					),
					|b| {
						b.iter(|| p *= &q);
					},
				);

				group.bench_function(
					BenchmarkId::new(
						"dot_product/128/shoup_vt",
						format!("{}/{}", degree, ctx.modulus().bits()),
					),
					|b| {
						b.iter(|| izip!(&p_vec, &q_vec).for_each(|(pi, qi)| p += pi * qi));
					},
				);

				p.allow_variable_time_computations();

				group.bench_function(
					BenchmarkId::new(
						"change_representation/PowerBasis_to_Ntt_vt",
						format!("{}/{}", degree, ctx.modulus().bits()),
					),
					|b| {
						b.iter(|| {
							p.override_representation(Representation::PowerBasis);
							p.change_representation(Representation::Ntt)
						});
					},
				);

				group.bench_function(
					BenchmarkId::new(
						"change_representation/Ntt_to_PowerBasis_vt",
						format!("{}/{}", degree, ctx.modulus().bits()),
					),
					|b| {
						b.iter(|| {
							p.override_representation(Representation::Ntt);
							p.change_representation(Representation::PowerBasis)
						});
					},
				);
			}
		}
	}

	group.finish();
}

criterion_group!(rq, rq_op_benchmark, rq_dot_product, rq_benchmark);
criterion_main!(rq);
