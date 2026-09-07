#![expect(missing_docs, reason = "examples/benches/tests omit docs by design")]
use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use fhe_math::rns::{RnsContext, RnsScaler, ScalingFactor};
use num_bigint::BigUint;
use rand::{Rng as RngCore, rng};
use std::{hint::black_box, sync::Arc, time::Duration};

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

    let mut rng = rng();
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

pub fn bfv_scaling_benchmark(c: &mut Criterion) {
    use fhe_math::{
        rq::{Context, Ntt, Poly, scaler::Scaler},
        zq::primes::generate_prime,
    };

    let mut group = c.benchmark_group("bfv_scaling");
    group.sample_size(50);
    group.warm_up_time(Duration::from_millis(200));
    group.measurement_time(Duration::from_millis(600));
    let degree = 8192;
    let base = [
        generate_prime(50, 2 * degree as u64, 1 << 50).unwrap(),
        generate_prime(55, 2 * degree as u64, 1 << 55).unwrap(),
    ];
    let mut extended = base.to_vec();
    let mut upper_bound = 1 << 62;
    for _ in 0..3 {
        upper_bound = generate_prime(62, 2 * degree as u64, upper_bound).unwrap();
        extended.push(upper_bound);
    }
    let base_rns = Arc::new(RnsContext::new(&base).unwrap());
    let extended_rns = Arc::new(RnsContext::new(&extended).unwrap());
    let factor = ScalingFactor::new(&BigUint::from(2056193u64), base_rns.modulus());
    let scaler = RnsScaler::new(&extended_rns, &base_rns, factor.clone());
    let mut rng = rng();
    let input: Vec<_> = extended.iter().map(|q| rng.next_u64() % q).collect();
    let mut output = vec![0; base.len()];
    group.bench_function("rns_downscale_5_to_2", |b| {
        b.iter(|| scaler.scale(black_box(input.as_slice()).into(), (&mut output).into(), 0));
    });

    let from = Context::new_arc(&extended, degree).unwrap();
    let to = Context::new_arc(&base, degree).unwrap();
    let poly_scaler = Scaler::new(&from, &to, factor).unwrap();
    for public in [false, true] {
        let mut input = Poly::<Ntt>::random(&from, &mut rng);
        if public {
            input.allow_variable_time_computations(fhe_util::VariableTime::new(
                fhe_util::PublicData::assert_public(),
            ));
        }
        group.bench_function(format!("ntt_downscale_5_to_2/public={public}"), |b| {
            b.iter(|| black_box(&input).scale(black_box(&poly_scaler)).unwrap());
        });
    }
    group.finish();
}

criterion_group!(rns, rns_benchmark, bfv_scaling_benchmark);
criterion_main!(rns);
