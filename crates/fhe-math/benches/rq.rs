#![expect(missing_docs, reason = "examples/benches/tests omit docs by design")]
#![expect(
    clippy::indexing_slicing,
    reason = "performance or example code relies on validated indices"
)]
use criterion::measurement::WallTime;
use criterion::{BenchmarkGroup, BenchmarkId, Criterion, criterion_group, criterion_main};
use fhe_math::rq::{Context, Ntt, NttShoup, Poly, PowerBasis, dot_product, traits::TryConvertFrom};
use itertools::{Itertools, izip};
use rand::rng;
use std::{sync::Arc, time::Duration};

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
            let p = Poly::<Ntt>::random(&ctx, &mut rng);
            let mut q = Poly::<Ntt>::random(&ctx, &mut rng);
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
            let p = Poly::<Ntt>::random(&ctx, &mut rng);
            let mut q = Poly::<Ntt>::random(&ctx, &mut rng);
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
            let mut p = Poly::<Ntt>::random(&ctx, &mut rng);
            let mut q = Poly::<Ntt>::random(&ctx, &mut rng);
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
        bench_op!(c, "rq_add", |p, q| p + q, vt);
        bench_op_assign!(
            c,
            "rq_add_assign",
            |p: &mut Poly<Ntt>, q: &Poly<Ntt>| *p += q,
            vt
        );
        bench_op!(c, "rq_sub", |p, q| p - q, vt);
        bench_op_assign!(
            c,
            "rq_sub_assign",
            |p: &mut Poly<Ntt>, q: &Poly<Ntt>| *p -= q,
            vt
        );
        bench_op!(c, "rq_mul", |p, q| p * q, vt);
        bench_op_assign!(
            c,
            "rq_mul_assign",
            |p: &mut Poly<Ntt>, q: &Poly<Ntt>| *p *= q,
            vt
        );
        bench_op_unary!(c, "rq_neg", |p: &Poly<_>| -p, vt);
    }
}

pub fn rq_dot_product(c: &mut Criterion) {
    let mut group = create_group(c, "rq_dot_product".to_string());
    let mut rng = rng();
    for degree in DEGREE {
        for i in [1, 4] {
            let ctx = Arc::new(Context::new(&MODULI[..i], *degree).unwrap());
            let p_vec = (0..256)
                .map(|_| Poly::<Ntt>::random(&ctx, &mut rng))
                .collect_vec();
            let q_vec = (0..256)
                .map(|_| Poly::<Ntt>::random(&ctx, &mut rng))
                .collect_vec();
            let mut out = Poly::<Ntt>::zero(&ctx);

            group.bench_function(
                BenchmarkId::from_parameter(format!("naive/{}/{}", degree, ctx.modulus().bits())),
                |b| {
                    b.iter(|| {
                        izip!(p_vec.iter(), q_vec.iter()).for_each(|(pi, qi)| out += &(pi * qi))
                    });
                },
            );

            let q_vec_shoup = q_vec
                .iter()
                .cloned()
                .map(Poly::<Ntt>::into_ntt_shoup)
                .collect_vec();
            group.bench_function(
                BenchmarkId::from_parameter(format!(
                    "naive_shoup/{}/{}",
                    degree,
                    ctx.modulus().bits()
                )),
                |b| {
                    b.iter(|| {
                        izip!(p_vec.iter(), q_vec_shoup.iter())
                            .for_each(|(pi, qi)| out += &(pi * qi))
                    });
                },
            );

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
            let mut p = Poly::<Ntt>::random(&ctx, &mut rng);
            let q = Poly::<NttShoup>::random(&ctx, &mut rng);

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

            let p_pb = Poly::<PowerBasis>::random(&ctx, &mut rng);
            group.bench_function(
                BenchmarkId::new(
                    "change_representation/PowerBasis_to_Ntt",
                    format!("{}/{}", degree, ctx.modulus().bits()),
                ),
                |b| {
                    b.iter(|| {
                        let _ = p_pb.clone().into_ntt();
                    });
                },
            );

            let p_ntt = Poly::<Ntt>::random(&ctx, &mut rng);
            group.bench_function(
                BenchmarkId::new(
                    "change_representation/Ntt_to_PowerBasis",
                    format!("{}/{}", degree, ctx.modulus().bits()),
                ),
                |b| {
                    b.iter(|| {
                        let _ = p_ntt.clone().into_power_basis();
                    });
                },
            );

            unsafe {
                let mut q_vt = q.clone();
                q_vt.allow_variable_time_computations();
                let mut p_vt = p.clone();
                p_vt.allow_variable_time_computations();

                group.bench_function(
                    BenchmarkId::new(
                        "mul_shoup_vt",
                        format!("{}/{}", degree, ctx.modulus().bits()),
                    ),
                    |b| {
                        b.iter(|| p_vt *= &q_vt);
                    },
                );

                let mut p_pb_vt = Poly::<PowerBasis>::random(&ctx, &mut rng);
                p_pb_vt.allow_variable_time_computations();

                group.bench_function(
                    BenchmarkId::new(
                        "change_representation/PowerBasis_to_Ntt_vt",
                        format!("{}/{}", degree, ctx.modulus().bits()),
                    ),
                    |b| {
                        b.iter(|| {
                            let _ = p_pb_vt.clone().into_ntt();
                        });
                    },
                );

                let mut p_ntt_vt = Poly::<Ntt>::random(&ctx, &mut rng);
                p_ntt_vt.allow_variable_time_computations();
                group.bench_function(
                    BenchmarkId::new(
                        "change_representation/Ntt_to_PowerBasis_vt",
                        format!("{}/{}", degree, ctx.modulus().bits()),
                    ),
                    |b| {
                        b.iter(|| {
                            let _ = p_ntt_vt.clone().into_power_basis();
                        });
                    },
                );
            }
        }
    }

    group.finish();
}

pub fn rq_convert_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("rq_convert");
    group.warm_up_time(Duration::from_millis(100));
    group.measurement_time(Duration::from_secs(1));

    for degree in DEGREE {
        for nmoduli in 1..=MODULI.len() {
            let ctx = Arc::new(Context::new(&MODULI[..nmoduli], *degree).unwrap());
            let v: Vec<u64> = (0..*degree as u64).collect();
            let slice = v.as_slice();

            group.bench_function(
                BenchmarkId::new("try_convert_from_slice", format!("{}/{}", degree, nmoduli)),
                |b| {
                    b.iter(|| Poly::<PowerBasis>::try_convert_from(slice, &ctx, false).unwrap());
                },
            );
        }
    }
    group.finish();
}

criterion_group!(
    rq,
    rq_op_benchmark,
    rq_dot_product,
    rq_benchmark,
    rq_convert_benchmark
);
criterion_main!(rq);
