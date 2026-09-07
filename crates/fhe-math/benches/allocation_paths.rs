#![expect(missing_docs, reason = "benchmark entry point")]
use criterion::{Criterion, criterion_group, criterion_main};
use fhe_math::{
    rns::RnsContext,
    rq::{Context, NttShoup, Poly, PowerBasis},
};

use std::hint::black_box;

fn benchmarks(c: &mut Criterion) {
    let moduli = [
        4611686018326724609,
        4611686018309947393,
        4611686018282684417,
    ];
    let ctx = Context::new_arc(&moduli, 2048).unwrap();
    let pb = Poly::<PowerBasis>::random_from_seed(&ctx, [7; 32]);
    let ntt = pb.clone().into_ntt();
    let shoup = ntt.clone().into_ntt_shoup();
    let bytes = shoup.to_bytes();
    let mut group = c.benchmark_group("allocation_paths");
    group.sample_size(20);
    group.warm_up_time(std::time::Duration::from_millis(100));
    group.measurement_time(std::time::Duration::from_millis(300));
    group.bench_function("context_setup", |b| {
        b.iter(|| Context::new_arc(black_box(&moduli), 2048).unwrap())
    });
    group.bench_function("add", |b| b.iter(|| black_box(&ntt) + black_box(&ntt)));
    group.bench_function("forward", |b| b.iter(|| black_box(&pb).clone().into_ntt()));
    group.bench_function("random_shoup", |b| {
        b.iter(|| Poly::<NttShoup>::random_from_seed(black_box(&ctx), [7; 32]))
    });
    group.bench_function("shoup", |b| {
        b.iter(|| black_box(&ntt).clone().into_ntt_shoup())
    });
    group.bench_function("level", |b| {
        b.iter(|| black_box(&ctx).context_at_level(2).unwrap())
    });
    let child = ctx.context_at_level(2).unwrap();
    group.bench_function("distance", |b| {
        b.iter(|| black_box(&ctx).niterations_to(black_box(&child)).unwrap())
    });
    group.bench_function("rns_setup", |b| {
        b.iter(|| RnsContext::new(black_box(&moduli)).unwrap())
    });
    group.bench_function("serialize_pb", |b| b.iter(|| black_box(&pb).to_bytes()));
    group.bench_function("serialize_shoup", |b| {
        b.iter(|| black_box(&shoup).to_bytes())
    });
    group.bench_function("deserialize_shoup", |b| {
        b.iter(|| Poly::<NttShoup>::from_bytes(black_box(&bytes), &ctx).unwrap())
    });
    for length in [4, 16, 256] {
        let left = vec![ntt.clone(); length];
        let right = vec![ntt.clone(); length];
        group.bench_function(format!("dot_product/{length}"), |b| {
            b.iter(|| fhe_math::rq::dot_product(black_box(&left), black_box(&right)).unwrap())
        });
    }
    for length in [4, 16, 256] {
        let left = vec![ntt.clone(); length];
        let right = left.clone();
        let mut workspace = fhe_math::rq::DotProductWorkspace::new(&ctx);
        let mut out = ntt.clone();
        group.bench_function(format!("dot_product_reuse/{length}"), |b| {
            b.iter(|| {
                workspace
                    .dot_product_into(black_box(&left), black_box(&right), &mut out)
                    .unwrap();
                black_box(&out);
            })
        });
    }

    group.finish();
}
criterion_group!(benches, benchmarks);
criterion_main!(benches);
