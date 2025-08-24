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
    assert_eq!(params.context_chain().iter_chain().count(), 2);

    // Test level-based access
    let ctx_level_0 = params.context_at_level(0)?;
    let ctx_level_1 = params.context_at_level(1)?;
    assert_eq!(ctx_level_0.moduli().len(), 2); // Full modulus chain
    assert_eq!(ctx_level_1.moduli().len(), 1); // Reduced modulus chain

    // Verify context chaining
    let head = params.context_chain();
    assert!(head.next.get().is_some());
    assert!(head.next.get().unwrap().next.get().is_none());

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

    Ok(())
}

#[test]
fn test_context_consistency() -> Result<(), Box<dyn std::error::Error>> {
    let params = BfvParametersBuilder::new()
        .set_degree(32)
        .set_plaintext_modulus(65537)
        .set_moduli_sizes(&[60, 60, 60])
        .build()?;

    // Verify context chain consistency
    for level in 0..=params.max_level() {
        let ctx = params.context_at_level(level)?;
        let level_ctx = params.context_level_at(level)?;

        // Verify cipher-plain context points to correct ciphertext context
        assert_eq!(&level_ctx.poly_context, ctx);

        // Verify modulus chain decreases
        let expected_moduli_count = params.moduli().len() - level;
        assert_eq!(ctx.moduli().len(), expected_moduli_count);
    }

    Ok(())
}
