use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use fhe_math::rns::{RnsContext, RnsScaler, ScalingFactor};
use num_bigint::BigUint;
use rand::{thread_rng, RngCore};
use std::sync::Arc;

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

    let rns_q = Arc::new(RnsContext::new(&q).unwrap());
    let rns_p = Arc::new(RnsContext::new(&p).unwrap());
    let scaler = RnsScaler::new(
        &rns_q,
        &rns_p,
        ScalingFactor::new(&BigUint::from(1u64), &BigUint::from(46116860181065u64)),
    );
    let scaler_as_converter = RnsScaler::new(&rns_q, &rns_p, ScalingFactor::one());

    let mut y = vec![0; p.len()];

    group.bench_function(
        BenchmarkId::new("scaler", format!("{}->{}", q.len(), p.len())),
        |b| {
            b.iter(|| scaler.scale((&x).into(), (&mut y).into(), 0));
        },
    );

    group.bench_function(
        BenchmarkId::new("scaler_as_converter", format!("{}->{}", q.len(), p.len())),
        |b| {
            b.iter(|| scaler_as_converter.scale((&x).into(), (&mut y).into(), 0));
        },
    );

    group.finish();
}

criterion_group!(rns, rns_benchmark);
criterion_main!(rns);
