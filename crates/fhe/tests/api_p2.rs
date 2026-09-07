//! Public P2 APIs, including one-pass inputs and immutable evaluation plans.
#![expect(clippy::indexing_slicing, reason = "fixtures have known lengths")]
use fhe::bfv::{
    Ciphertext, Encoding, Parameters, Plaintext, SecretKey,
    evaluation::{
        DotProductScalarWorkspace, MultiplicationPlan, RelinearizationKey, dot_product_scalar,
        dot_product_scalar_iter,
    },
    packing::{PackedPlaintext, PackedPlaintextBatch},
};
use fhe::{PublicData, VariableTime};
use std::cell::Cell;

fn parameters() -> fhe::Result<Parameters> {
    Parameters::builder()
        .degree(16)
        .plaintext_modulus(1153_u64)
        .ciphertext_modulus_bits([62, 62, 62])
        .build()
}

#[test]
fn chunks_are_ordinary_vectors_with_explicit_empty_and_signed_semantics() -> fhe::Result<()> {
    let par = parameters()?;
    for level in 0..=par.max_level() {
        for encoding in [Encoding::Polynomial, Encoding::Simd] {
            for length in [0_usize, 1, 15, 16, 17, 32, 35] {
                let values: Vec<_> = (0..length).map(|i| i as i64 - 10).collect();
                let chunks: Vec<Plaintext> =
                    Plaintext::encode_chunks_signed_at_level(&par, &values, encoding, level)?;
                assert_eq!(chunks.len(), length.div_ceil(par.degree()).max(1));
                let mut decoded = Vec::new();
                // Owned iteration no longer requires a custom batch wrapper.
                for pt in chunks {
                    assert_eq!(pt.level(), level);
                    decoded.extend(pt.decode_signed(encoding)?);
                }
                assert_eq!(&decoded[..length], values);
                assert!(decoded[length..].iter().all(|value| *value == 0));
            }
        }
    }
    Ok(())
}

#[test]
fn slice_and_single_pass_dot_products_preserve_results_and_timing() -> fhe::Result<()> {
    let par = parameters()?;
    let mut rng = rand::rng();
    let sk = SecretKey::generate(&par, &mut rng);
    let permission = VariableTime::new(PublicData::assert_public());
    for length in [1, 16, 17, 33] {
        // fused and long-product paths for 62-bit primes
        let ciphertexts: Vec<_> = (0..length)
            .map(|_| sk.encrypt(&Plaintext::encode(&par, &[3, 5], Encoding::Simd)?, &mut rng))
            .collect::<fhe::Result<_>>()?;
        for public in [true, false] {
            let pt = if public {
                Plaintext::encode_public(&par, &[2, 7], Encoding::Simd, permission)?
            } else {
                Plaintext::encode(&par, &[2, 7], Encoding::Simd)?
            };
            let plaintexts = vec![pt; length];
            let expected = dot_product_scalar(&ciphertexts, &plaintexts)?;
            let left_visits = Cell::new(0);
            let right_visits = Cell::new(0);
            let actual = dot_product_scalar_iter(
                ciphertexts
                    .iter()
                    .inspect(|_| left_visits.set(left_visits.get() + 1)),
                plaintexts
                    .iter()
                    .inspect(|_| right_visits.set(right_visits.get() + 1)),
            )?;
            assert_eq!(actual.to_bytes(), expected.to_bytes());
            assert_eq!((left_visits.get(), right_visits.get()), (length, length));
            let mut workspace = DotProductScalarWorkspace::new(&par, 0)?;
            let batch = PackedPlaintextBatch::try_from_iter(&par, 0, &plaintexts)?;
            let views: Vec<_> = (&batch).into_iter().collect();
            for result in [
                workspace.dot_product_scalar(&ciphertexts, &plaintexts)?,
                workspace.dot_product_scalar_refs(
                    &ciphertexts.iter().collect::<Vec<_>>(),
                    &plaintexts.iter().collect::<Vec<_>>(),
                )?,
                workspace.dot_product_scalar_packed(&ciphertexts, &views)?,
                workspace.dot_product_scalar_packed_iter(&ciphertexts, &batch)?,
            ] {
                assert_eq!(result.to_bytes(), expected.to_bytes());
                assert!(
                    result
                        .components()
                        .iter()
                        .all(|p| p.allows_variable_time_computations() == public)
                );
            }
            assert_eq!(
                sk.decrypt(&actual)?.decode(Encoding::Simd)?[..2],
                [6 * length as u64, 35 * length as u64 % 1153]
            );
            assert!(workspace.dot_product_scalar(&ciphertexts, &[]).is_err());
            assert!(
                workspace
                    .dot_product_scalar(&ciphertexts, &plaintexts[..length - 1])
                    .is_err()
            );
            assert_eq!(
                workspace.dot_product_scalar(&ciphertexts, &plaintexts)?,
                expected
            );
        }
    }
    Ok(())
}

#[test]
fn packed_batch_extension_is_atomic_and_borrowed_iteration_is_conventional() -> fhe::Result<()> {
    let par = parameters()?;
    let public = VariableTime::new(PublicData::assert_public());
    let a = Plaintext::encode_public(&par, &[3], Encoding::Polynomial, public)?;
    let b = Plaintext::encode(&par, &[5], Encoding::Polynomial)?;
    let bad = Plaintext::encode_at_level(&par, &[9], Encoding::Polynomial, 1)?;
    let mut batch = PackedPlaintextBatch::try_from_iter(&par, 0, [&a])?;
    let size = batch.size_bytes();
    assert!(batch.try_extend([&b, &bad]).is_err());
    assert_eq!(batch.len(), 1);
    assert_eq!(batch.size_bytes(), size);
    batch.try_extend([&b, &a])?;
    let mut rows = (&batch).into_iter();
    assert_eq!(rows.len(), 3);
    assert!(rows.next().is_some());
    assert!(rows.next_back().is_some());
    assert_eq!(rows.len(), 1);
    assert!(rows.next().is_some());
    assert!(rows.next().is_none());
    assert!(rows.next_back().is_none());
    assert!(batch.get(3).is_none());
    let ct = Ciphertext::trivial_zero(&par, 0)?;
    let mut workspace = DotProductScalarWorkspace::new(&par, 0)?;
    assert_eq!(
        workspace.dot_product_scalar_packed_iter([&ct; 3], &batch)?,
        workspace.dot_product_scalar_iter([&ct; 3], [&a, &b, &a])?
    );
    batch.clear();
    assert!(batch.is_empty());
    assert_eq!(batch.size_bytes(), 8);
    batch.push(&b)?;
    assert_eq!(batch.len(), 1);
    assert_eq!(PackedPlaintext::from(&b).unpack(), b);
    Ok(())
}

#[test]
fn immutable_plans_validate_named_configuration_and_share_across_threads() -> fhe::Result<()> {
    fn send_sync<T: Send + Sync>() {}
    send_sync::<Parameters>();
    send_sync::<MultiplicationPlan<'_>>();
    send_sync::<fhe::bfv::evaluation::PreparedMultiplicand<'_>>();
    send_sync::<DotProductScalarWorkspace>();
    send_sync::<PackedPlaintextBatch>();
    let par = parameters()?;
    let mut rng = rand::rng();
    let sk = SecretKey::generate(&par, &mut rng);
    let rk = RelinearizationKey::new(&sk, &mut rng)?;
    let builder = MultiplicationPlan::builder(&par)
        .relinearization(&rk)
        .modulus_switching(true);
    let plan = builder.build()?;
    let ct = sk.encrypt(&Plaintext::encode(&par, &[3, 5], Encoding::Simd)?, &mut rng)?;
    let prepared = plan.prepare_lhs(&ct)?;
    let expected = plan.square(&ct)?;
    std::thread::scope(|scope| {
        let result = scope
            .spawn(|| prepared.multiply(&ct))
            .join()
            .unwrap()
            .unwrap();
        assert_eq!(result, expected);
    });
    assert_eq!(expected.component_count(), 2);
    assert_eq!(expected.level(), 1);
    assert!(
        MultiplicationPlan::builder(&par)
            .level(par.max_level())
            .modulus_switching(true)
            .build()
            .is_err()
    );
    assert!(
        MultiplicationPlan::builder(&par)
            .level(1)
            .relinearization(&rk)
            .build()
            .is_err()
    );
    assert!(
        MultiplicationPlan::builder(&par)
            .extended_basis([])
            .build()
            .is_err()
    );
    assert_eq!(prepared.multiply(&ct)?, expected);
    Ok(())
}
