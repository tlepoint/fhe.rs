//! Checked raw-transform and decoding boundaries on either selected backend.
use fhe_math::{
    DecodeLimits, Error, PublicData, VariableTime,
    error::PolynomialSerializationError,
    ntt::NttOperator,
    rq::{Context, Ntt, NttShoup, Poly, PowerBasis},
    zq::Modulus,
};
use std::error::Error as StdError;

#[test]
fn checked_ntt_is_atomic_and_matches_the_existing_transforms() -> fhe_math::Result<()> {
    let public = VariableTime::new(PublicData::assert_public());
    for (size, modulus) in [(32, 1153), (1024, 4611686018326724609)] {
        let q = Modulus::new(modulus)?;
        let op = NttOperator::new(&q, size).unwrap();
        assert_eq!((op.size(), op.modulus()), (size, modulus));
        let input: Vec<_> = (0..size).map(|i| (i as u64 * 19) % modulus).collect();
        let mut reference = input.clone();
        op.forward(&mut reference);
        let mut checked = input.clone();
        op.try_forward(&mut checked)?;
        assert_eq!(reference, checked);
        let mut variable = input.clone();
        op.try_forward_public(&mut variable, public)?;
        assert_eq!(reference, variable);
        op.try_backward(&mut checked)?;
        op.try_backward_public(&mut variable, public)?;
        assert_eq!(input, checked);
        assert_eq!(input, variable);
        for invalid in [
            vec![0; size - 1],
            vec![0; size + 1],
            vec![modulus; size],
            vec![u64::MAX; size],
        ] {
            let mut data = invalid.clone();
            assert!(op.try_forward(&mut data).is_err());
            assert_eq!(data, invalid);
            assert!(op.try_backward(&mut data).is_err());
            assert_eq!(data, invalid);
            assert!(op.try_forward_public(&mut data, public).is_err());
            assert_eq!(data, invalid);
            assert!(op.try_backward_public(&mut data, public).is_err());
            assert_eq!(data, invalid);
        }
        // Existing safe slice methods reject wrong sizes even in release builds.
        for inverse in [false, true] {
            assert!(
                std::panic::catch_unwind(|| {
                    let mut invalid = vec![0; size - 1];
                    if inverse {
                        op.backward(&mut invalid)
                    } else {
                        op.forward(&mut invalid)
                    }
                })
                .is_err()
            );
        }
    }
    Ok(())
}

#[test]
fn polynomial_decoding_bounds_expansion_and_preserves_error_sources() -> fhe_math::Result<()> {
    let ctx = Context::new_arc(&[1153, 12289], 16)?;
    let pb = Poly::<PowerBasis>::from_coefficients(&[1, 2], &ctx)?;
    let limits = DecodeLimits {
        max_residue_bytes: 16 * 2 * 16,
        ..DecodeLimits::default()
    };
    assert_eq!(
        pb,
        Poly::<PowerBasis>::from_bytes_with_limits(&pb.to_bytes(), &ctx, &limits)?
    );
    let limited = DecodeLimits {
        max_residue_bytes: limits.max_residue_bytes - 1,
        ..limits
    };
    assert!(matches!(
        Poly::<PowerBasis>::from_bytes_with_limits(&pb.to_bytes(), &ctx, &limited),
        Err(Error::DecodeLimit(_))
    ));
    let limit = DecodeLimits {
        max_input_bytes: 0,
        ..limits
    };
    assert!(matches!(
        Poly::<Ntt>::from_bytes_with_limits(&[0x80], &ctx, &limit),
        Err(Error::DecodeLimit(_))
    ));
    assert!(matches!(
        Poly::<NttShoup>::from_bytes_with_limits(&[0x80], &ctx, &limit),
        Err(Error::DecodeLimit(_))
    ));
    for error in [
        Poly::<PowerBasis>::from_bytes(&[0x80], &ctx).unwrap_err(),
        Poly::<Ntt>::from_bytes(&[0x80], &ctx).unwrap_err(),
        Poly::<NttShoup>::from_bytes(&[0x80], &ctx).unwrap_err(),
    ] {
        assert!(matches!(
            error,
            Error::PolynomialSerialization(PolynomialSerializationError::Decode { .. })
        ));
        assert!(
            error
                .source()
                .unwrap()
                .source()
                .unwrap()
                .downcast_ref::<prost::DecodeError>()
                .is_some()
        );
    }
    Ok(())
}
