#![allow(missing_docs, clippy::indexing_slicing)]
use fhe::bfv::{BfvParametersBuilder, Ciphertext, Encoding, Plaintext, SecretKey};
use fhe_traits::{FheDecoder, FheDecrypter, FheEncoder, FheEncrypter};
use num_bigint::BigUint;
use rand::rng;
use std::error::Error;

#[test]
fn test_biguint_plaintext_encryption_decryption() -> Result<(), Box<dyn Error>> {
    let mut rng = rng();

    // Choose a large plaintext modulus: 2^127 - 1 (Mersenne prime M127)
    // 170141183460469231731687303715884105727
    let p_str = "170141183460469231731687303715884105727";
    let p = BigUint::parse_bytes(p_str.as_bytes(), 10).unwrap();

    // Create parameters
    // We need enough ciphertext moduli to support the plaintext modulus + noise.
    // p is 127 bits. Noise adds ~20-30 bits (at least).
    // So we need ~160 bits of ciphertext moduli.
    // 3 moduli of 60 bits = 180 bits.
    let params = BfvParametersBuilder::new()
        .set_degree(16)
        .set_plaintext_modulus_biguint(p.clone())
        .set_moduli_sizes(&[60, 60, 60])
        .build_arc()?;

    let sk = SecretKey::random(&params, &mut rng);

    // Create a vector of BigUint values
    let mut values = vec![BigUint::from(0u32); params.degree()];
    values[0] = BigUint::from(123456789u64);
    values[1] = p.clone() - 1u32; // -1
    values[2] = p.clone() / 2u32;

    let pt = Plaintext::try_encode(values.as_slice(), Encoding::poly(), &params)?;

    let ct: Ciphertext = sk.try_encrypt(&pt, &mut rng)?;

    let decrypted_pt = sk.try_decrypt(&ct)?;

    // Decode
    let decrypted_values: Vec<BigUint> =
        Vec::<BigUint>::try_decode(&decrypted_pt, Encoding::poly())?;

    assert_eq!(decrypted_values, values);

    Ok(())
}

#[test]
fn test_biguint_homomorphic_addition() -> Result<(), Box<dyn Error>> {
    let mut rng = rng();

    let p_str = "170141183460469231731687303715884105727";
    let p = BigUint::parse_bytes(p_str.as_bytes(), 10).unwrap();

    let params = BfvParametersBuilder::new()
        .set_degree(16)
        .set_plaintext_modulus_biguint(p.clone())
        .set_moduli_sizes(&[60, 60, 60])
        .build_arc()?;

    let sk = SecretKey::random(&params, &mut rng);

    let val1 = BigUint::from(100u32);
    let val2 = p.clone() - 50u32; // -50

    let mut vec1 = vec![BigUint::from(0u32); params.degree()];
    vec1[0] = val1.clone();

    let mut vec2 = vec![BigUint::from(0u32); params.degree()];
    vec2[0] = val2.clone();

    let pt1 = Plaintext::try_encode(vec1.as_slice(), Encoding::poly(), &params)?;
    let pt2 = Plaintext::try_encode(vec2.as_slice(), Encoding::poly(), &params)?;

    let ct1: Ciphertext = sk.try_encrypt(&pt1, &mut rng)?;
    let ct2: Ciphertext = sk.try_encrypt(&pt2, &mut rng)?;

    let ct_res = &ct1 + &ct2;

    let decrypted_pt = sk.try_decrypt(&ct_res)?;
    let decrypted_values: Vec<BigUint> =
        Vec::<BigUint>::try_decode(&decrypted_pt, Encoding::poly())?;

    // 100 + (-50) = 50
    assert_eq!(decrypted_values[0], BigUint::from(50u32));

    Ok(())
}
