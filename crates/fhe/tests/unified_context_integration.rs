// Integration test to verify the unified context management API works correctly
use fhe::bfv::BfvParametersBuilder;

#[test]
fn test_unified_context_api() -> Result<(), Box<dyn std::error::Error>> {
    // Create parameters using the builder
    let params = BfvParametersBuilder::new()
        .set_degree(16)
        .set_plaintext_modulus(1153)
        .set_moduli_sizes(&[62, 62])
        .build()?;

    // Test basic properties
    assert_eq!(params.degree(), 16);
    assert_eq!(params.plaintext(), 1153);
    assert_eq!(params.max_level(), 1);

    // Test unified context access
    let fhe_context = params.fhe_context();
    assert_eq!(fhe_context.polynomial_degree, 16);
    assert_eq!(fhe_context.plaintext_modulus, 1153);
    assert_eq!(fhe_context.ciphertext_contexts.len(), 2);
    assert_eq!(fhe_context.cipher_plain_contexts.len(), 2);

    // Test level-based access
    let ctx_level_0 = params.context_at_level(0)?;
    let ctx_level_1 = params.context_at_level(1)?;
    assert_eq!(ctx_level_0.moduli().len(), 2); // Full modulus chain
    assert_eq!(ctx_level_1.moduli().len(), 1); // Reduced modulus chain

    // Test cipher-plain context access
    let cp_ctx_0 = params.cipher_plain_context_at_level(0)?;
    let cp_ctx_1 = params.cipher_plain_context_at_level(1)?;

    // Verify context chaining
    assert!(cp_ctx_0.next.is_some());
    assert!(cp_ctx_1.next.is_none()); // Last level has no next

    // Verify that contexts are correctly set up
    assert_eq!(ctx_level_0.moduli().len(), 2);
    assert_eq!(ctx_level_1.moduli().len(), 1);

    Ok(())
}

#[test]
fn test_unified_context_error_handling() -> Result<(), Box<dyn std::error::Error>> {
    let params = BfvParametersBuilder::new()
        .set_degree(16)
        .set_plaintext_modulus(1153)
        .set_moduli_sizes(&[62, 62])
        .build()?;

    // Test invalid level access
    assert!(params.context_at_level(2).is_err());
    assert!(params.cipher_plain_context_at_level(2).is_err());

    Ok(())
}

#[test]
fn test_context_consistency() -> Result<(), Box<dyn std::error::Error>> {
    let params = BfvParametersBuilder::new()
        .set_degree(32)
        .set_plaintext_modulus(65537)
        .set_moduli_sizes(&[60, 60, 60])
        .build()?;

    let fhe_context = params.fhe_context();

    // Verify context chain consistency
    for level in 0..=params.max_level() {
        let ctx = params.context_at_level(level)?;
        let cp_ctx = params.cipher_plain_context_at_level(level)?;

        // Verify cipher-plain context points to correct ciphertext context
        assert_eq!(&cp_ctx.ciphertext_context, ctx);

        // Verify plaintext context is consistent
        assert_eq!(cp_ctx.plaintext_context, fhe_context.plaintext_context);

        // Verify modulus chain decreases
        let expected_moduli_count = fhe_context.ciphertext_moduli.len() - level;
        assert_eq!(ctx.moduli().len(), expected_moduli_count);
    }

    Ok(())
}
