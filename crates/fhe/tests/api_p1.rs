//! Consumer-facing API contracts. Ordinary BFV operations require only fhe and
//! rand.
#![expect(clippy::indexing_slicing, reason = "tests use known nonempty fixtures")]

use fhe::bfv::{
    Ciphertext, Encoding, Parameters, Plaintext, PublicKey, SecretKey,
    evaluation::CiphertextProductAccumulator, evaluation::EvaluationKey,
    evaluation::MultiplicationPlan, evaluation::RelinearizationKey, evaluation::RgswCiphertext,
    packing::PackedPlaintext, packing::PackedPlaintextBatch,
};
use fhe::{CiphertextError, Error, ParametersError, PublicData, VariableTime};

fn parameters() -> fhe::Result<Parameters> {
    Parameters::builder()
        .degree(16)
        .plaintext_modulus(1153_u64)
        .ciphertext_modulus_bits([50, 50, 50])
        .build()
}

#[test]
fn inherent_api_round_trip_and_explicit_message_interpretation() -> fhe::Result<()> {
    let par = parameters()?;
    let mut rng = rand::rng();
    let sk = SecretKey::generate(&par, &mut rng);
    let pk = PublicKey::from_secret_key(&sk, &mut rng);
    for level in 0..=par.max_level() {
        for encoding in [Encoding::Polynomial, Encoding::Simd] {
            let a = Plaintext::encode_at_level(&par, &[20], encoding, level)?;
            let b = Plaintext::encode_signed_at_level(&par, &[-7], encoding, level)?;
            let ca = sk.encrypt(&a, &mut rng)?;
            let cb = pk.encrypt(&b, &mut rng)?;
            let product = ca.multiply(&cb)?;
            assert_eq!(product.level(), level);
            assert_eq!(product.component_count(), 3);
            assert_eq!(sk.decrypt(&ca)?, a);
            let decoded = sk.decrypt(&product)?.decode_signed(encoding)?;
            assert_eq!(decoded[0], -140);
            assert_eq!(decoded.len(), par.degree());
            assert!(decoded[1..].iter().all(|value| *value == 0));
            assert_eq!(a.parameters(), &par);
        }
    }
    Ok(())
}

#[test]
fn independent_and_imported_parameter_handles_are_compatible() -> fhe::Result<()> {
    let original = parameters()?;
    let rebuilt = parameters()?;
    let imported = Parameters::from_bytes(&original.to_bytes())?;
    let mut rng = rand::rng();
    let sk = SecretKey::generate(&original, &mut rng);
    let pk = PublicKey::from_secret_key(&sk, &mut rng);
    let bytes = sk.export_secret_bytes();
    for par in [&original.clone(), &rebuilt, &imported] {
        assert_eq!(original, *par);
        assert!(original.compatible(par));
        // This imports the same secret, rather than conflating key identity
        // with parameter compatibility.
        let imported_sk = SecretKey::from_bytes(&bytes, par)?;
        let imported_pk = PublicKey::from_bytes(&pk.to_bytes(), par)?;
        let pt = Plaintext::encode(par, &[2, 3], Encoding::Simd)?;
        let ct = sk.encrypt(&pt, &mut rng)?;
        let other = Ciphertext::from_bytes(&ct.to_bytes(), par)?;
        assert_eq!(ct, other);
        assert_eq!(
            imported_sk.decrypt(&imported_pk.encrypt(&pt, &mut rng)?)?,
            pt
        );
        assert_eq!(
            sk.decrypt(&ct.add(&other)?)?.decode(Encoding::Simd)?[..2],
            [4, 6]
        );
        assert_eq!(
            sk.decrypt(&ct.multiply(&other)?)?.decode(Encoding::Simd)?[..2],
            [4, 9]
        );
        assert_eq!(ct.multiply_plaintext(&pt)?, other.multiply_plaintext(&pt)?);
        let mut packed = PackedPlaintextBatch::with_capacity(&original, 0, 1)?;
        packed.push(&pt)?;
        assert_eq!(PackedPlaintext::from(&pt).unpack(), pt);
        let mut accumulator = CiphertextProductAccumulator::new(&original, 0)?;
        accumulator.add_product(&ct, &other)?;
        assert_eq!(accumulator.finish()?, ct.multiply(&other)?);
        let rk = RelinearizationKey::new(&imported_sk, &mut rng)?;
        let plan = MultiplicationPlan::with_relinearization(&rk)?;
        assert_eq!(
            sk.decrypt(&plan.multiply(&ct, &other)?)?
                .decode(Encoding::Simd)?[..2],
            [4, 9]
        );
    }
    Ok(())
}

#[test]
fn every_defining_parameter_setting_participates_in_compatibility() -> fhe::Result<()> {
    let par = parameters()?;
    let builder = Parameters::builder()
        .degree(par.degree())
        .plaintext_modulus(par.plaintext_modulus().clone())
        .ciphertext_moduli(par.moduli())
        .noise_variance(par.noise_variance());
    let mut reversed = par.moduli().to_vec();
    reversed.reverse();
    let changed = [
        builder.clone().degree(8).build()?,
        builder.clone().plaintext_modulus(17_u64).build()?,
        builder.clone().noise_variance(11).build()?,
        builder.clone().ciphertext_moduli(reversed).build()?,
        builder.ciphertext_modulus_bits([49, 50, 50]).build()?,
    ];
    for other in changed {
        assert_ne!(par, other);
        assert!(!par.compatible(&other));
        assert!(!other.compatible(&par));
    }
    Ok(())
}

#[test]
fn checked_arithmetic_errors_leave_receivers_unchanged() -> fhe::Result<()> {
    let par = parameters()?;
    let other = Parameters::builder()
        .degree(par.degree())
        .plaintext_modulus(1153_u64)
        .ciphertext_moduli(par.moduli())
        .noise_variance(11)
        .build()?;
    let mut rng = rand::rng();
    let sk = SecretKey::generate(&par, &mut rng);
    let pt = Plaintext::encode(&par, &[2, 3], Encoding::Simd)?;
    let ct = sk.encrypt(&pt, &mut rng)?;
    let encoded = ct.to_bytes();
    let invalid = [
        Ciphertext::trivial_zero(&other, 0)?,
        Ciphertext::trivial_zero(&par, 1)?,
    ];
    for rhs in &invalid {
        assert!(ct.add(rhs).is_err());
        assert!(ct.subtract(rhs).is_err());
        assert!(ct.multiply(rhs).is_err());
        let mut receiver = ct.clone();
        assert!(receiver.add_assign(rhs).is_err());
        assert_eq!(receiver.to_bytes(), encoded);
        assert!(receiver.subtract_assign(rhs).is_err());
        assert_eq!(receiver.to_bytes(), encoded);
        assert!(receiver.multiply_assign(rhs).is_err());
        assert_eq!(receiver.to_bytes(), encoded);
    }
    for rhs in [Plaintext::zero(&other, 0)?, Plaintext::zero(&par, 1)?] {
        let mut receiver = ct.clone();
        assert!(receiver.add_plaintext(&rhs).is_err());
        assert!(receiver.subtract_plaintext(&rhs).is_err());
        assert!(receiver.multiply_plaintext(&rhs).is_err());
        assert!(receiver.add_plaintext_assign(&rhs).is_err());
        assert_eq!(receiver.to_bytes(), encoded);
        assert!(receiver.subtract_plaintext_assign(&rhs).is_err());
        assert_eq!(receiver.to_bytes(), encoded);
        assert!(receiver.multiply_plaintext_assign(&rhs).is_err());
        assert_eq!(receiver.to_bytes(), encoded);
    }
    let squared = ct.square()?;
    let mut receiver = ct.clone();
    assert!(matches!(
        receiver.add_assign(&squared),
        Err(Error::Ciphertext(
            CiphertextError::ComponentCountMismatch { .. }
        ))
    ));
    assert!(receiver.subtract_assign(&squared).is_err());
    assert_eq!(receiver.to_bytes(), encoded);
    // Higher component counts remain supported by ordinary multiplication.
    receiver.multiply_assign(&squared)?;
    assert_eq!(receiver.component_count(), 4);
    assert_eq!(receiver, ct.multiply(&squared)?);
    assert_eq!(receiver.level(), ct.level());
    let mut receiver = ct.clone();
    receiver.square_assign()?;
    assert_eq!(receiver, squared);
    receiver = ct.clone();
    receiver.add_plaintext_assign(&pt)?;
    assert_eq!(receiver, ct.add_plaintext(&pt)?);
    receiver.subtract_plaintext_assign(&pt)?;
    assert_eq!(receiver, ct);
    receiver.multiply_plaintext_assign(&pt)?;
    assert_eq!(receiver, ct.multiply_plaintext(&pt)?);
    assert_eq!(
        Ciphertext::from_bytes(&receiver.to_bytes(), &par)?,
        receiver
    );
    Ok(())
}

#[test]
fn encoding_reduces_inputs_and_chunks_without_storing_metadata() -> fhe::Result<()> {
    let par = parameters()?;
    let t = par.plaintext_modulus_u64().unwrap();
    let permission = VariableTime::new(PublicData::assert_public());
    for encoding in [Encoding::Polynomial, Encoding::Simd] {
        let values = [t, t + 1, u64::MAX];
        let pt = Plaintext::encode(&par, &values, encoding)?;
        assert_eq!(&pt.decode(encoding)?[..3], &[0, 1, u64::MAX % t]);
        assert_eq!(
            pt,
            Plaintext::encode_public(&par, &values, encoding, permission)?
        );
        assert_eq!(
            Plaintext::encode(&par, &[], encoding)?,
            Plaintext::zero(&par, 0)?
        );
        assert!(Plaintext::encode(&par, &vec![0; par.degree() + 1], encoding).is_err());
        let values = vec![t + 3; par.degree() + 1];
        let chunks = Plaintext::encode_chunks_at_level(&par, &values, encoding, 1)?;
        assert_eq!(chunks.len(), 2);
        assert_eq!(chunks[0].decode(encoding)?, vec![3; par.degree()]);
        assert_eq!(chunks[1].decode(encoding)?[0], 3);
        assert!(chunks[1].decode(encoding)?[1..].iter().all(|x| *x == 0));
        assert!(chunks.iter().all(|pt| pt.level() == 1));
    }
    assert!(
        Plaintext::encode_at_level(&par, &[1], Encoding::Polynomial, par.max_level() + 1).is_err()
    );
    Ok(())
}

#[test]
fn consuming_builders_replace_choices_and_validate_at_build() -> fhe::Result<()> {
    let base = Parameters::builder().degree(16).plaintext_modulus(1153_u64);
    let par = base.clone().ciphertext_modulus_bits([50, 50]).build()?;
    let explicit = base
        .clone()
        .ciphertext_modulus_bits([1])
        .ciphertext_moduli(par.moduli())
        .build()?;
    let generated = base
        .ciphertext_moduli([0])
        .ciphertext_modulus_bits([50, 50])
        .build()?;
    assert_eq!(par, explicit);
    assert_eq!(par, generated);
    assert!(matches!(
        Parameters::builder().build(),
        Err(Error::ParametersError(ParametersError::MissingDegree))
    ));
    assert!(matches!(
        Parameters::builder().degree(16).build(),
        Err(Error::ParametersError(
            ParametersError::MissingPlaintextModulus
        ))
    ));
    let mut rng = rand::rng();
    let sk = SecretKey::generate(&par, &mut rng);
    assert!(
        EvaluationKey::builder(&sk)
            .ciphertext_level(par.max_level() + 1)
            .build(&mut rng)
            .is_err()
    );
    assert!(
        EvaluationKey::builder(&sk)
            .key_level(1)
            .build(&mut rng)
            .is_err()
    );
    assert!(
        EvaluationKey::builder(&sk)
            .enable_column_rotation(0)
            .build(&mut rng)
            .is_err()
    );
    assert!(
        EvaluationKey::builder(&sk)
            .enable_expansion(5)
            .build(&mut rng)
            .is_err()
    );
    let key = EvaluationKey::builder(&sk)
        .ciphertext_level(1)
        .key_level(0)
        .enable_inner_sum()
        .enable_row_rotation()
        .enable_column_rotation(1)
        .enable_expansion(4)
        .build(&mut rng)?;
    let key = EvaluationKey::from_bytes(&key.to_bytes(), &par)?;
    let pt = Plaintext::encode_at_level(&par, &[1; 16], Encoding::Simd, 1)?;
    let ct = sk.encrypt(&pt, &mut rng)?;
    assert_eq!(
        sk.decrypt(&key.inner_sum(&ct)?)?.decode(Encoding::Simd)?,
        vec![16; 16]
    );
    assert_eq!(sk.decrypt(&key.rotate_rows(&ct)?)?, pt);
    assert_eq!(sk.decrypt(&key.rotate_columns(&ct, 1)?)?, pt);
    Ok(())
}

#[test]
fn profiles_are_selected_by_degree_and_built_on_demand() -> fhe::Result<()> {
    let profiles = Parameters::profiles_128(16)?.collect::<Vec<_>>();
    assert_eq!(
        profiles.iter().map(|p| p.degree()).collect::<Vec<_>>(),
        [1024, 2048, 4096]
    );
    let selected = profiles.into_iter().find(|p| p.degree() == 2048).unwrap();
    let par = Parameters::profile_128(2048, 16)?;
    assert_eq!(par, selected.build()?);
    assert!(Parameters::profile_128(1234, 16).is_err());
    for bits in [0, 1, 64, usize::MAX] {
        assert!(Parameters::profiles_128(bits).is_err());
    }
    Ok(())
}

#[test]
fn rgsw_has_explicit_encryption_and_checked_external_product() -> fhe::Result<()> {
    let par = parameters()?;
    let mut rng = rand::rng();
    let sk = SecretKey::generate(&par, &mut rng);
    let pt = Plaintext::encode(&par, &[2, 3], Encoding::Simd)?;
    let ct = sk.encrypt(&pt, &mut rng)?;
    let rgsw = sk.encrypt_rgsw(&pt, &mut rng)?;
    let imported = RgswCiphertext::from_bytes(&rgsw.to_bytes(), &parameters()?)?;
    let product = ct.multiply_rgsw(&imported)?;
    assert_eq!(sk.decrypt(&product)?.decode(Encoding::Simd)?[..2], [4, 9]);
    assert_eq!(product.level(), ct.level());
    let mut switched = product.clone();
    switched.switch_to_level(switched.max_switchable_level())?;
    assert_eq!(
        sk.decrypt(&switched)?.decode(Encoding::Simd)?,
        sk.decrypt(&product)?.decode(Encoding::Simd)?
    );
    assert!(ct.square()?.multiply_rgsw(&rgsw).is_err());
    let mut lower = Ciphertext::trivial_zero(&par, 1)?;
    let original = lower.to_bytes();
    assert!(lower.multiply_rgsw(&rgsw).is_err());
    assert!(lower.multiply_rgsw_assign(&rgsw).is_err());
    assert_eq!(lower.to_bytes(), original);
    let mut in_place = ct.clone();
    in_place.multiply_rgsw_assign(&imported)?;
    assert_eq!(in_place, product);
    let foreign = Parameters::builder()
        .degree(16)
        .plaintext_modulus(17_u64)
        .ciphertext_moduli(par.moduli())
        .build()?;
    assert!(
        Ciphertext::trivial_zero(&foreign, 0)?
            .multiply_rgsw(&rgsw)
            .is_err()
    );
    Ok(())
}

#[test]
fn imports_p0_protobuf_fixtures_without_a_wire_format_migration() -> fhe::Result<()> {
    #[cfg(feature = "tfhe-ntt")]
    macro_rules! fixture {
        ($file:literal) => {
            include_bytes!(concat!("data/p0/tfhe/", $file))
        };
    }
    #[cfg(not(feature = "tfhe-ntt"))]
    macro_rules! fixture {
        ($file:literal) => {
            include_bytes!(concat!("data/p0/", $file))
        };
    }
    let par_bytes = fixture!("parameters.bin");
    let secret_bytes = fixture!("secret_key.bin");
    let public_bytes = fixture!("public_key.bin");
    let ct_bytes = fixture!("ciphertext.bin");
    let rk_bytes = fixture!("relinearization_key.bin");
    let rgsw_bytes = fixture!("rgsw_ciphertext.bin");
    let ek_bytes = fixture!("evaluation_key.bin");
    let par = Parameters::from_bytes(par_bytes)?;
    let sk = SecretKey::from_bytes(secret_bytes, &par)?;
    let pk = PublicKey::from_bytes(public_bytes, &par)?;
    let ct = Ciphertext::from_bytes(ct_bytes, &par)?;
    let rk = RelinearizationKey::from_bytes(rk_bytes, &par)?;
    let rgsw = RgswCiphertext::from_bytes(rgsw_bytes, &par)?;
    let ek = EvaluationKey::from_bytes(ek_bytes, &par)?;
    assert_eq!(par.to_bytes(), par_bytes);
    assert_eq!(sk.export_secret_bytes().as_slice(), secret_bytes);
    assert_eq!(pk.to_bytes(), public_bytes);
    assert_eq!(ct.to_bytes(), ct_bytes);
    assert_eq!(rk.to_bytes(), rk_bytes);
    assert_eq!(rgsw.to_bytes(), rgsw_bytes);
    let expected = Plaintext::encode(&par, &[2, 3], Encoding::Simd)?;
    assert_eq!(sk.decrypt(&ct)?, expected);
    assert_eq!(
        sk.decrypt(&pk.encrypt(&expected, &mut rand::rng())?)?,
        expected
    );
    let mut squared = ct.square()?;
    rk.relinearize(&mut squared)?;
    assert_eq!(sk.decrypt(&squared)?.decode(Encoding::Simd)?[..2], [4, 9]);
    assert_eq!(
        sk.decrypt(&ct.multiply_rgsw(&rgsw)?)?,
        sk.decrypt(&squared)?
    );
    assert_eq!(
        sk.decrypt(&ek.inner_sum(&ct)?)?.decode(Encoding::Simd)?,
        vec![5; 16]
    );
    let imported_ek = EvaluationKey::from_bytes(&ek.to_bytes(), &par)?;
    assert_eq!(imported_ek, ek);
    Ok(())
}
