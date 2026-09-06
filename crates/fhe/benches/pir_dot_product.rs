//! Focused scalar dot products with the PIR examples' shapes and dense data.

use criterion::{BenchmarkId, Criterion, SamplingMode, criterion_group, criterion_main};
use fhe::bfv::{
    BfvParametersBuilder, Ciphertext, DotProductScalarWorkspace, Encoding, PackedPlaintext,
    PackedPlaintextVec, Plaintext, SecretKey,
};
use fhe_traits::{FheEncoderVariableTime, FheEncrypter, PublicData, VariableTime};
use rand::{RngExt, SeedableRng};
use rand_chacha::ChaCha8Rng;
use std::{hint::black_box, time::Duration};

fn scalar(c: &mut Criterion) {
    let mut group = c.benchmark_group("pir_dot_product");
    group.sampling_mode(SamplingMode::Flat);
    group.sample_size(10);
    group.warm_up_time(Duration::from_millis(100));
    group.measurement_time(Duration::from_millis(600));
    let public = VariableTime::new(PublicData::assert_public());
    for (name, degree, t, moduli, rows, columns) in [
        ("mulpir", 8192, 1_785_857, [50, 55, 55], 174, 81),
        ("sealpir", 4096, 2_056_193, [36, 36, 37], 447, 64),
    ] {
        let par = BfvParametersBuilder::new()
            .set_degree(degree)
            .set_plaintext_modulus(t)
            .set_moduli_sizes(&moduli)
            .build_arc()
            .unwrap();
        let encoding = Encoding::poly_at_level(1);
        let mut rng = ChaCha8Rng::seed_from_u64(degree as u64);
        let sk = SecretKey::random(&par, &mut rng);
        let selection = Plaintext::try_encode_vt(&[1][..], encoding.clone(), &par, public).unwrap();
        let query: Vec<Ciphertext> = (0..rows)
            .map(|_| sk.try_encrypt(&selection, &mut rng).unwrap())
            .collect();
        // Independent dense plaintexts: do not benefit from the examples' mostly
        // zero generated records or repeatedly reference a single allocation.
        let coefficients: Vec<Vec<u32>> = (0..rows)
            .map(|_| (0..degree).map(|_| rng.random_range(0..t) as u32).collect())
            .collect();
        let encode = |coefficients: &Vec<u32>| {
            let values: Vec<u64> = coefficients.iter().map(|&x| u64::from(x)).collect();
            Plaintext::try_encode_vt(values.as_slice(), encoding.clone(), &par, public).unwrap()
        };
        let column: Vec<_> = coefficients.iter().map(encode).collect();
        let mut workspace = DotProductScalarWorkspace::new(&par, 1).unwrap();
        let expected = workspace
            .dot_product_scalar(query.iter(), column.iter())
            .unwrap();
        group.bench_function(BenchmarkId::new("column", name), |b| {
            b.iter(|| {
                workspace
                    .dot_product_scalar(black_box(query.iter()), black_box(column.iter()))
                    .unwrap()
            });
        });
        // A compact coefficient-only database needs these transforms per query.
        let rebuilt: Vec<_> = coefficients.iter().map(encode).collect();
        assert_eq!(
            expected,
            workspace
                .dot_product_scalar(query.iter(), rebuilt.iter())
                .unwrap()
        );
        drop(rebuilt);
        group.bench_function(BenchmarkId::new("encode_and_column", name), |b| {
            b.iter(|| {
                let rebuilt: Vec<_> = black_box(&coefficients).iter().map(encode).collect();
                workspace
                    .dot_product_scalar(query.iter(), rebuilt.iter())
                    .unwrap()
            });
        });
        let packed: Vec<_> = column.iter().map(PackedPlaintext::from).collect();
        assert_eq!(
            expected,
            workspace
                .dot_product_scalar_packed(query.iter(), packed.iter())
                .unwrap()
        );
        group.bench_function(BenchmarkId::new("packed_column", name), |b| {
            b.iter(|| {
                workspace
                    .dot_product_scalar_packed(query.iter(), packed.iter())
                    .unwrap()
            });
        });
        let database: Vec<_> = (0..rows * columns)
            .map(|_| {
                let values: Vec<_> = (0..degree).map(|_| rng.random_range(0..t)).collect();
                Plaintext::try_encode_vt(values.as_slice(), encoding.clone(), &par, public).unwrap()
            })
            .collect();
        group.bench_function(BenchmarkId::new("matrix", name), |b| {
            b.iter(|| {
                for column in 0..columns {
                    black_box(
                        workspace
                            .dot_product_scalar(
                                query.iter(),
                                database.iter().skip(column).step_by(columns),
                            )
                            .unwrap(),
                    );
                }
            });
        });
        let mut packed_database =
            PackedPlaintextVec::with_capacity(&par, 1, database.len()).unwrap();
        for pt in database {
            packed_database.push(&pt).unwrap();
        }
        group.bench_function(BenchmarkId::new("packed_matrix", name), |b| {
            b.iter(|| {
                for column in 0..columns {
                    black_box(
                        workspace
                            .dot_product_scalar_packed(
                                query.iter(),
                                packed_database.iter().skip(column).step_by(columns),
                            )
                            .unwrap(),
                    );
                }
            });
        });
    }
    group.finish();
}

criterion_group!(benches, scalar);
criterion_main!(benches);
