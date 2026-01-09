// Allow indexing in benchmarks for convenience
#![allow(clippy::indexing_slicing)]
#![allow(missing_docs)]

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use fhe_math::{ntt::NttOperator, zq::Modulus};
use rand::rng;
use std::sync::Arc;

pub fn ntt_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("ntt");
    group.sample_size(50);
    let mut rng = rng();

    for vector_size in [1024usize, 4096].iter() {
        for p in [4611686018326724609u64, 40961u64] {
            let p_nbits = 64 - p.leading_zeros();
            let q = Modulus::new(p).unwrap();
            let mut a = q.random_vec(*vector_size, &mut rng);
            let op = NttOperator::new(&Arc::new(q), *vector_size).unwrap();

            group.bench_function(
                BenchmarkId::new("forward", format!("{vector_size}/{p_nbits}")),
                |b| b.iter(|| op.forward(&mut a)),
            );

            group.bench_function(
                BenchmarkId::new("forward_vt", format!("{vector_size}/{p_nbits}")),
                |b| b.iter(|| unsafe { op.forward_vt(a.as_mut_ptr()) }),
            );

            group.bench_function(
                BenchmarkId::new("backward", format!("{vector_size}/{p_nbits}")),
                |b| b.iter(|| op.backward(&mut a)),
            );

            group.bench_function(
                BenchmarkId::new("backward_vt", format!("{vector_size}/{p_nbits}")),
                |b| b.iter(|| unsafe { op.backward_vt(a.as_mut_ptr()) }),
            );
        }
    }

    group.finish();
}

criterion_group!(ntt, ntt_benchmark);
criterion_main!(ntt);
