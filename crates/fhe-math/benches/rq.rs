// Allow indexing in benchmarks for convenience
#![allow(clippy::indexing_slicing)]

use criterion::measurement::WallTime;
use criterion::{criterion_group, criterion_main, BenchmarkGroup, BenchmarkId, Criterion};
use fhe_math::rq::*;
use itertools::{izip, Itertools};
use rand::rng;
use std::{
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
    sync::Arc,
    time::Duration,
};

static MODULI: &[u64; 4] = &[
    562949954093057,
    4611686018326724609,
    4611686018309947393,
    4611686018282684417,
];

static DEGREE: &[usize] = &[1024, 2048, 4096, 8192];

fn create_group(c: &mut Criterion, name: String) -> BenchmarkGroup<'_, WallTime> {
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
        let mut rng = rng();

        for degree in DEGREE {
            let ctx = Arc::new(Context::new(&MODULI[..1], *degree).unwrap());
            let p = Poly::random(&ctx, Representation::Ntt, &mut rng);
            let mut q = Poly::random(&ctx, Representation::Ntt, &mut rng);
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
        let mut rng = rng();

        for degree in DEGREE {
            let ctx = Arc::new(Context::new(&MODULI[..1], *degree).unwrap());
            let p = Poly::random(&ctx, Representation::Ntt, &mut rng);
            let mut q = Poly::random(&ctx, Representation::Ntt, &mut rng);
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
        let mut rng = rng();

        for degree in DEGREE {
            let ctx = Arc::new(Context::new(&MODULI[..1], *degree).unwrap());
            let mut p = Poly::random(&ctx, Representation::Ntt, &mut rng);
            let mut q = Poly::random(&ctx, Representation::Ntt, &mut rng);
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
    let mut rng = rng();
    for degree in DEGREE {
        for i in [1, 4] {
            let ctx = Arc::new(Context::new(&MODULI[..i], *degree).unwrap());
            let p_vec = (0..256)
                .map(|_| Poly::random(&ctx, Representation::Ntt, &mut rng))
                .collect_vec();
            let mut q_vec = (0..256)
                .map(|_| Poly::random(&ctx, Representation::Ntt, &mut rng))
                .collect_vec();
            let mut out = Poly::zero(&ctx, Representation::Ntt);

            group.bench_function(
                BenchmarkId::from_parameter(format!("naive/{}/{}", degree, ctx.modulus().bits())),
                |b| {
                    b.iter(|| {
                        izip!(p_vec.iter(), q_vec.iter()).for_each(|(pi, qi)| out += &(pi * qi))
                    });
                },
            );

            q_vec
                .iter_mut()
                .for_each(|qi| qi.change_representation(Representation::NttShoup));
            group.bench_function(
                BenchmarkId::from_parameter(format!(
                    "naive_shoup/{}/{}",
                    degree,
                    ctx.modulus().bits()
                )),
                |b| {
                    b.iter(|| {
                        izip!(p_vec.iter(), q_vec.iter()).for_each(|(pi, qi)| out += &(pi * qi))
                    });
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
        }
    }
}

pub fn rq_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("rq");
    group.warm_up_time(Duration::from_millis(100));
    group.measurement_time(Duration::from_secs(1));

    let mut rng = rng();
    for degree in DEGREE {
        for nmoduli in 1..=MODULI.len() {
            if !nmoduli.is_power_of_two() {
                continue;
            }
            let ctx = Arc::new(Context::new(&MODULI[..nmoduli], *degree).unwrap());
            let mut p = Poly::random(&ctx, Representation::Ntt, &mut rng);
            let mut q = Poly::random(&ctx, Representation::Ntt, &mut rng);
            q.change_representation(Representation::NttShoup);

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
                q.change_representation(Representation::NttShoup);

                group.bench_function(
                    BenchmarkId::new(
                        "mul_shoup_vt",
                        format!("{}/{}", degree, ctx.modulus().bits()),
                    ),
                    |b| {
                        b.iter(|| p *= &q);
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
