use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use crrust_math::zq::Modulus;
use rand::RngCore;

fn random_vector(size: usize, p: u64) -> Vec<u64> {
    let mut rng = rand::thread_rng();
    let mut v = vec![];
    for _ in 0..size {
        v.push(rng.next_u64() % p)
    }
    v
}

pub fn zq_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("zq");
    group.sample_size(50);

    let p = 4611686018326724609;

    for vector_size in [1024usize, 4096].iter() {
        let mut a = random_vector(*vector_size, p);
        let c = random_vector(*vector_size, p);
        let q = Modulus::new(p).unwrap();

        group.bench_function(BenchmarkId::new("add_vec", vector_size), |b| {
            b.iter(|| q.add_vec(&mut a, &c));
        });

        group.bench_function(BenchmarkId::new("sub_vec", vector_size), |b| {
            b.iter(|| q.sub_vec(&mut a, &c));
        });

        group.bench_function(BenchmarkId::new("neg_vec", vector_size), |b| {
            b.iter(|| q.neg_vec(&mut a));
        });

        group.bench_function(BenchmarkId::new("mul_vec", vector_size), |b| {
            b.iter(|| q.mul_vec(&mut a, &c));
        });

        group.bench_function(BenchmarkId::new("mul_opt_vec", vector_size), |b| {
            b.iter(|| q.mul_opt_vec(&mut a, &c));
        });
    }

    group.finish();
}

criterion_group!(zq, zq_benchmark);
criterion_main!(zq);
