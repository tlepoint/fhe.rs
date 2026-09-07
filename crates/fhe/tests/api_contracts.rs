//! Public API regressions for validated values and sensitive-data handling.
#![expect(clippy::indexing_slicing, reason = "tests use known nonempty fixtures")]

use fhe::bfv::{
    Ciphertext, Encoding, Parameters, ParametersBuilder, Plaintext, SecretKey,
    evaluation::EvaluationKeyBuilder,
};
use fhe::{PublicData, SecretDependentDiagnostics, VariableTime};
use fhe_math::rq::{Ntt, NttShoup, Poly, PowerBasis, SubstitutionExponent};
use std::sync::Arc;

fn parameters() -> fhe::Result<Parameters> {
    ParametersBuilder::new()
        .degree(16)
        .plaintext_modulus(1153_u64)
        .ciphertext_modulus_bits([62, 62])
        .build()
}

#[test]
fn plaintext_equality_is_an_equivalence_relation() -> fhe::Result<()> {
    let par = parameters()?;
    let sk = SecretKey::generate(&par, &mut rand::rng());
    let unknown = sk.decrypt(&Ciphertext::trivial_zero(&par, 0)?)?;
    let poly = Plaintext::zero(&par, 0)?;
    let simd = Plaintext::zero(&par, 0)?;
    assert_eq!(poly, unknown);
    assert_eq!(simd, unknown);
    assert_eq!(poly, simd);
    let nonzero = Plaintext::encode(&par, &[1_u64], Encoding::Polynomial)?;
    let lower = Plaintext::zero(&par, 1)?;
    let values = [
        unknown.clone(),
        unknown,
        poly.clone(),
        poly,
        simd,
        nonzero,
        lower,
    ];
    for a in &values {
        assert_eq!(a, a);
        for b in &values {
            assert_eq!(a == b, b == a);
            for c in &values {
                if a == b && b == c {
                    assert_eq!(a, c);
                }
            }
        }
    }
    Ok(())
}

#[test]
fn trivial_zeros_are_valid_at_every_level() -> fhe::Result<()> {
    let par = parameters()?;
    let mut rng = rand::rng();
    let sk = SecretKey::generate(&par, &mut rng);
    for level in 0..=par.max_level() {
        let zero = Ciphertext::trivial_zero(&par, level)?;
        assert_eq!(zero.level(), level);
        assert_eq!(zero.component_count(), 2);
        assert!(Parameters::compatible(zero.parameters(), &par));
        let restored = Ciphertext::from_bytes(&zero.to_bytes(), &par)?;
        assert_eq!(restored, zero);
        assert_eq!(restored.to_bytes(), zero.to_bytes());
        let encoding = Encoding::Polynomial;
        for value in [&zero, &zero.square()?] {
            assert_eq!(sk.decrypt(value)?.decode(encoding)?, vec![0; par.degree()],);
        }
        let pt = Plaintext::encode_at_level(&par, &[17_u64], encoding, level)?;
        let ct: Ciphertext = sk.encrypt(&pt, &mut rng)?;
        assert_eq!(ct.add(&zero).unwrap(), ct);
        assert_eq!(ct.subtract(&zero).unwrap(), ct);
    }
    assert!(Ciphertext::trivial_zero(&par, par.max_level() + 1).is_err());
    Ok(())
}

#[test]
fn component_import_rejects_invalid_structure_and_lazy_residues() -> fhe::Result<()> {
    let par = parameters()?;
    let head = par.context_at_level(0)?;
    let tail = par.context_at_level(1)?;
    assert!(Ciphertext::from_components(vec![], &par).is_err());
    assert!(Ciphertext::from_components(vec![Poly::zero(head)], &par).is_err());
    assert!(Ciphertext::from_components(vec![Poly::zero(head), Poly::zero(tail)], &par).is_err());
    let lazy = Poly::create_constant_ntt_polynomial_with_lazy_coefficients_and_variable_time(
        &[3; 16],
        head,
        VariableTime::new(PublicData::assert_public()),
    );
    assert!(!lazy.is_canonical());
    assert!(matches!(
        Ciphertext::from_components(vec![lazy; 2], &par),
        Err(fhe::Error::Ciphertext(
            fhe::error::CiphertextError::NonCanonicalPolynomial
        ))
    ));
    let lower = Ciphertext::from_components(vec![Poly::zero(tail); 3], &par)?;
    assert_eq!(lower.level(), 1);
    assert_eq!(lower.component_count(), 3);
    assert_eq!(Ciphertext::from_bytes(&lower.to_bytes(), &par)?, lower);
    Ok(())
}

#[test]
fn rebuilding_components_discards_seed_without_changing_structural_equality() -> fhe::Result<()> {
    let par = parameters()?;
    let mut rng = rand::rng();
    let sk = SecretKey::generate(&par, &mut rng);
    let pt = Plaintext::encode(&par, &[7_u64], Encoding::Polynomial)?;
    let original: Ciphertext = sk.encrypt(&pt, &mut rng)?;
    let rebuilt = Ciphertext::from_components(original.clone().into_components(), &par)?;
    assert_eq!(rebuilt, original);
    assert!(rebuilt.to_bytes().len() > original.to_bytes().len());
    for index in 0..original.component_count() {
        let mut parts = original.clone().into_components();
        parts[index] = -&parts[index];
        let modified = Ciphertext::from_components(parts, &par)?;
        let restored = Ciphertext::from_bytes(&modified.to_bytes(), &par)?;
        assert_eq!(modified, restored);
        assert_eq!(modified.to_bytes(), restored.to_bytes());
        assert_eq!(sk.decrypt(&modified)?, sk.decrypt(&restored)?);
    }
    Ok(())
}

#[test]
fn debug_redacts_secret_key_plaintext_polynomial_and_builder() -> fhe::Result<()> {
    let par = parameters()?;
    let mut rng = rand::rng();
    let a = SecretKey::generate(&par, &mut rng);
    let b = SecretKey::generate(&par, &mut rng);
    assert_eq!(format!("{a:?}"), format!("{b:?}"));
    assert_eq!(format!("{a:?}"), "SecretKey { degree: 16, .. }");
    let builder_a = EvaluationKeyBuilder::new(&a);
    let builder_b = EvaluationKeyBuilder::new(&b);
    // HashMap iteration order is unspecified; check the secret's redacted field.
    for text in [format!("{builder_a:?}"), format!("{builder_b:?}")] {
        assert!(text.contains("sk: SecretKey { degree: 16, .. }"));
        assert!(!text.contains("coeffs"));
    }
    let pa = Plaintext::encode(&par, &[123_u64], Encoding::Polynomial)?;
    let pb = Plaintext::encode(&par, &[456_u64], Encoding::Polynomial)?;
    assert_eq!(format!("{pa:?}"), format!("{pb:?}"));
    let ctx = par.context_at_level(0)?;
    let a = Poly::<PowerBasis>::from_coefficients(&[123_u64], ctx)?;
    let b = Poly::<PowerBasis>::from_coefficients(&[456_u64], ctx)?;
    assert_eq!(format!("{a:?}"), format!("{b:?}"));
    assert!(!format!("{a:?}").contains("coefficients"));
    Ok(())
}

#[test]
fn timing_permissions_and_diagnostic_acknowledgment_are_explicit() -> fhe::Result<()> {
    let par = parameters()?;
    let ctx = par.context_at_level(0)?;
    let restricted = Poly::<PowerBasis>::from_coefficients(&[3_u64], ctx)?;
    let public = Poly::<PowerBasis>::from_coefficients_with_timing(
        &[3_u64],
        ctx,
        Some(VariableTime::new(PublicData::assert_public())),
    )?;
    assert!(!restricted.allows_variable_time_computations());
    assert!(public.allows_variable_time_computations());
    assert_eq!(restricted, public);
    let sk = SecretKey::generate(&par, &mut rand::rng());
    let zero = Ciphertext::trivial_zero(&par, 0)?;
    assert_eq!(
        sk.measure_noise_vartime(&zero, SecretDependentDiagnostics::acknowledge_leakage(),)?,
        0
    );
    Ok(())
}

#[test]
fn descriptors_keep_their_validated_values() -> fhe::Result<()> {
    let par = parameters()?;
    let ctx = par.context_at_level(0)?;
    let exponent = SubstitutionExponent::new(ctx, 2 * par.degree() + 3)?;
    assert_eq!(exponent.exponent(), 3);
    assert!(SubstitutionExponent::new(ctx, 2).is_err());
    let poly = Poly::<PowerBasis>::from_coefficients(&[1_u64, 2, 3], ctx)?;
    assert_eq!(
        poly.substitute(&exponent)?.into_ntt(),
        poly.into_ntt().substitute(&exponent)?,
    );
    assert!(Arc::ptr_eq(par.context_level_at(0)?.poly_context(), ctx));
    Ok(())
}

#[test]
fn immutable_crypto_types_inherit_send_and_sync() {
    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<Parameters>();
    assert_send_sync::<Plaintext>();
    assert_send_sync::<Ciphertext>();
    assert_send_sync::<SecretKey>();
    assert_send_sync::<Poly<PowerBasis>>();
    assert_send_sync::<Poly<Ntt>>();
    assert_send_sync::<Poly<NttShoup>>();
    assert_send_sync::<EvaluationKeyBuilder<'static>>();
}
