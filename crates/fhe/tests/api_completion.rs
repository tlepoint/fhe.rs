//! Remaining construction and wire-boundary contracts from the API proposal.
use fhe::bfv::{
    Ciphertext, Encoding, Parameters, Plaintext, PublicKey, SecretKey,
    evaluation::{EvaluationKey, RelinearizationKey, RgswCiphertext},
};
use fhe::{DecodeLimits, Error, error::SerializationError};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use std::error::Error as StdError;

fn parameters() -> fhe::Result<Parameters> {
    Parameters::builder()
        .degree(16)
        .plaintext_modulus(1153_u64)
        .ciphertext_modulus_bits([50, 50, 50])
        .build()
}

#[test]
fn relinearization_builder_keeps_levels_distinct_and_validates_before_rng() -> fhe::Result<()> {
    let par = parameters()?;
    let mut rng = ChaCha8Rng::seed_from_u64(47);
    let sk = SecretKey::generate(&par, &mut rng);
    let mut identical_rng = rng.clone();
    assert_eq!(
        RelinearizationKey::new(&sk, &mut rng)?,
        RelinearizationKey::builder(&sk).build(&mut identical_rng)?
    );
    for (ct_level, key_level) in [(0, 0), (1, 0), (1, 1), (2, 0), (2, 1)] {
        let key = RelinearizationKey::builder(&sk)
            .ciphertext_level(ct_level)
            .key_level(key_level)
            .build(&mut rng)?;
        let pt = Plaintext::encode_at_level(&par, &[3], Encoding::Simd, ct_level)?;
        let mut ct = sk.encrypt(&pt, &mut rng)?.square()?;
        key.relinearize(&mut ct)?;
        assert_eq!(ct.component_count(), 2);
        assert_eq!(ct.level(), ct_level);
        assert_eq!(sk.decrypt(&ct)?.decode(Encoding::Simd)?.first(), Some(&9));
    }
    for (ct_level, key_level) in [(3, 0), (0, 3), (0, 2), (0, 1)] {
        let mut copy = rng.clone();
        assert!(
            RelinearizationKey::builder(&sk)
                .ciphertext_level(ct_level)
                .key_level(key_level)
                .build(&mut rng)
                .is_err()
        );
        assert_eq!(rng.next_u64(), copy.next_u64());
    }
    Ok(())
}

#[test]
fn every_public_bfv_import_enforces_input_limits_and_keeps_decode_sources() -> fhe::Result<()> {
    let par = parameters()?;
    let limits = DecodeLimits {
        max_input_bytes: 0,
        ..DecodeLimits::default()
    };
    let invalid = &[0x80][..];
    let errors = [
        Parameters::from_bytes_with_limits(invalid, &limits).unwrap_err(),
        SecretKey::from_bytes_with_limits(invalid, &par, &limits).unwrap_err(),
        PublicKey::from_bytes_with_limits(invalid, &par, &limits).unwrap_err(),
        Ciphertext::from_bytes_with_limits(invalid, &par, &limits).unwrap_err(),
        RelinearizationKey::from_bytes_with_limits(invalid, &par, &limits).unwrap_err(),
        EvaluationKey::from_bytes_with_limits(invalid, &par, &limits).unwrap_err(),
        RgswCiphertext::from_bytes_with_limits(invalid, &par, &limits).unwrap_err(),
    ];
    for error in errors {
        assert!(matches!(error, Error::DecodeLimit(_)));
    }
    for error in [
        Parameters::from_bytes(invalid).unwrap_err(),
        SecretKey::from_bytes(invalid, &par).unwrap_err(),
        PublicKey::from_bytes(invalid, &par).unwrap_err(),
        Ciphertext::from_bytes(invalid, &par).unwrap_err(),
        RelinearizationKey::from_bytes(invalid, &par).unwrap_err(),
        EvaluationKey::from_bytes(invalid, &par).unwrap_err(),
        RgswCiphertext::from_bytes(invalid, &par).unwrap_err(),
    ] {
        assert!(matches!(
            error,
            Error::SerializationError(SerializationError::Decode { .. })
        ));
        let source = error.source().unwrap().source().unwrap();
        assert!(source.downcast_ref::<prost::DecodeError>().is_some());
    }
    // A syntactically valid field with the wrong schema wire type reaches
    // prost's decoder and retains its message/field context too.
    let error = SecretKey::from_bytes(&[13, 0, 0, 0, 0], &par).unwrap_err();
    assert!(
        error
            .source()
            .unwrap()
            .source()
            .unwrap()
            .to_string()
            .contains("SecretKey.coeffs")
    );
    Ok(())
}

#[test]
fn preflight_counts_packed_unpacked_duplicate_and_seeded_fields() -> fhe::Result<()> {
    let par = parameters()?;
    let limits = DecodeLimits {
        max_degree: 8,
        ..DecodeLimits::default()
    };
    // A huge degree is rejected before context construction; even if overwritten.
    assert!(matches!(
        Parameters::from_bytes_with_limits(&[8, 16, 8, 8], &limits),
        Err(Error::DecodeLimit(_))
    ));
    let limits = DecodeLimits {
        max_moduli: 1,
        ..DecodeLimits::default()
    };
    for bytes in [
        &[16, 17, 16, 19][..],
        &[18, 2, 17, 19][..],
        &[18, 1, 17, 16, 19][..],
    ] {
        assert!(matches!(
            Parameters::from_bytes_with_limits(bytes, &limits),
            Err(Error::DecodeLimit(_))
        ));
    }
    let limits = DecodeLimits {
        max_plaintext_bytes: 1,
        ..DecodeLimits::default()
    };
    assert!(matches!(
        Parameters::from_bytes_with_limits(&[42, 2, 1, 1], &limits),
        Err(Error::DecodeLimit(_))
    ));
    let limits = DecodeLimits {
        max_polynomials: 1,
        ..DecodeLimits::default()
    };
    assert!(matches!(
        Ciphertext::from_bytes_with_limits(&[10, 0, 10, 0], &par, &limits),
        Err(Error::DecodeLimit(_))
    ));
    let mut rng = ChaCha8Rng::seed_from_u64(19);
    let sk = SecretKey::generate(&par, &mut rng);
    let ct = sk.encrypt(&Plaintext::zero(&par, 0)?, &mut rng)?;
    // The seed also creates a polynomial, though only one is stored on the wire.
    assert!(matches!(
        Ciphertext::from_bytes_with_limits(&ct.to_bytes(), &par, &limits),
        Err(Error::DecodeLimit(_))
    ));
    let limits = DecodeLimits {
        max_polynomials: 2,
        max_residue_bytes: 2 * 16 * 3 * 16,
        ..DecodeLimits::default()
    };
    assert_eq!(
        ct,
        Ciphertext::from_bytes_with_limits(&ct.to_bytes(), &par, &limits)?
    );
    assert!(
        Ciphertext::from_bytes_with_limits(
            &ct.to_bytes(),
            &par,
            &DecodeLimits {
                max_residue_bytes: limits.max_residue_bytes - 1,
                ..limits
            }
        )
        .is_err()
    );
    let key = RelinearizationKey::builder(&sk).build(&mut rng)?;
    assert!(matches!(
        RelinearizationKey::from_bytes_with_limits(&key.to_bytes(), &par, &limits),
        Err(Error::DecodeLimit(_))
    ));
    assert_eq!(key, RelinearizationKey::from_bytes(&key.to_bytes(), &par)?);
    Ok(())
}

#[test]
fn default_wire_policy_preserves_unknown_fields_and_rejects_truncation() -> fhe::Result<()> {
    let par = parameters()?;
    let mut bytes = par.to_bytes();
    // Unknown varint, fixed32, fixed64, bytes, and a matched legacy group.
    bytes.extend_from_slice(&[
        80, 7, 93, 0, 0, 0, 0, 97, 0, 0, 0, 0, 0, 0, 0, 0, 106, 1, 1, 115, 8, 1, 116,
    ]);
    assert_eq!(Parameters::from_bytes(&bytes)?, par);
    bytes.push(0x80);
    assert!(Parameters::from_bytes(&bytes).is_err());
    assert!(Parameters::from_bytes(&[42, 0xff, 0xff, 0xff, 0xff, 0x7f]).is_err());
    Ok(())
}
