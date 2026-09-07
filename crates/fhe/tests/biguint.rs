#![allow(missing_docs, clippy::indexing_slicing)]
use fhe::bfv::{
    Ciphertext, Encoding, Parameters, ParametersBuilder, Plaintext, RelinearizationKey, SecretKey,
};

use num_bigint::BigUint;
use rand::rng;
use std::error::Error;

fn parameters() -> Parameters {
    // Choose a large plaintext modulus: 2^127 - 1 (Mersenne prime M127)
    // 170141183460469231731687303715884105727
    let p_str = "170141183460469231731687303715884105727";
    let p = BigUint::parse_bytes(p_str.as_bytes(), 10).unwrap();

    // Create parameters
    ParametersBuilder::new()
        .degree(16)
        .plaintext_modulus(p.clone())
        .ciphertext_modulus_bits([60, 60, 60, 60, 60])
        .build()
        .unwrap()
}

#[test]
fn test_biguint_plaintext_encryption_decryption() -> Result<(), Box<dyn Error>> {
    let mut rng = rng();

    let params = parameters();
    let sk = SecretKey::generate(&params, &mut rng);

    // Create a vector of BigUint values
    let mut values = vec![BigUint::from(0u32); params.degree()];
    values[0] = BigUint::from(123456789u64);
    values[1] = params.plaintext_modulus() - 1u32; // -1
    values[2] = params.plaintext_modulus() / 2u32;

    let pt = Plaintext::encode_biguint(&params, values.as_slice(), Encoding::Polynomial)?;

    let ct: Ciphertext = sk.encrypt(&pt, &mut rng)?;

    let decrypted_pt = sk.decrypt(&ct)?;

    // Decode
    let decrypted_values: Vec<BigUint> = decrypted_pt.decode_biguint(Encoding::Polynomial)?;

    assert_eq!(decrypted_values, values);

    Ok(())
}

#[test]
fn test_biguint_homomorphic_addition() -> Result<(), Box<dyn Error>> {
    let mut rng = rng();

    let params = parameters();
    let sk = SecretKey::generate(&params, &mut rng);

    let val1 = BigUint::from(10u32);
    let val2 = params.plaintext_modulus() - 50u32; // -50

    let mut vec1 = vec![BigUint::from(0u32); params.degree()];
    vec1[0] = val1.clone();

    let mut vec2 = vec![BigUint::from(0u32); params.degree()];
    vec2[0] = val2.clone();

    let pt1 = Plaintext::encode_biguint(&params, vec1.as_slice(), Encoding::Polynomial)?;
    let pt2 = Plaintext::encode_biguint(&params, vec2.as_slice(), Encoding::Polynomial)?;

    let ct1: Ciphertext = sk.encrypt(&pt1, &mut rng)?;
    let ct2: Ciphertext = sk.encrypt(&pt2, &mut rng)?;

    let ct_res = ct1.add(&ct2).unwrap();

    let decrypted_pt = sk.decrypt(&ct_res)?;
    let decrypted_values: Vec<BigUint> = decrypted_pt.decode_biguint(Encoding::Polynomial)?;

    // 10 + (-50) = -40
    assert_eq!(
        decrypted_values[0],
        params.plaintext_modulus() - BigUint::from(40u32)
    );

    Ok(())
}

#[test]
fn test_biguint_multiplication_without_relin() -> Result<(), Box<dyn Error>> {
    let mut rng = rng();

    let params = parameters();
    let sk = SecretKey::generate(&params, &mut rng);

    let val1 = BigUint::from(10u32);
    let val2 = params.plaintext_modulus() - BigUint::from(20u32);

    let mut vec1 = vec![BigUint::from(0u32); params.degree()];
    vec1[0] = val1.clone();

    let mut vec2 = vec![BigUint::from(0u32); params.degree()];
    vec2[0] = val2.clone();

    let pt1 = Plaintext::encode_biguint(&params, vec1.as_slice(), Encoding::Polynomial)?;
    let pt2 = Plaintext::encode_biguint(&params, vec2.as_slice(), Encoding::Polynomial)?;

    let ct1: Ciphertext = sk.encrypt(&pt1, &mut rng)?;
    let ct2: Ciphertext = sk.encrypt(&pt2, &mut rng)?;

    let ct_res = ct1.multiply(&ct2).unwrap();

    assert_eq!(ct_res.components().len(), 3); // Degree increases

    let decrypted_pt = sk.decrypt(&ct_res)?;
    let decrypted_values: Vec<BigUint> = decrypted_pt.decode_biguint(Encoding::Polynomial)?;

    // 10 * (-20) = -200
    assert_eq!(
        decrypted_values[0],
        params.plaintext_modulus() - BigUint::from(200u32)
    );

    Ok(())
}

#[test]
fn test_biguint_multiplication_with_relin() -> Result<(), Box<dyn Error>> {
    let mut rng = rng();

    // Use default parameters with biguint
    let params = ParametersBuilder::new()
        .degree(16)
        .plaintext_modulus(1153_u64)
        .ciphertext_modulus_bits([62usize; 3])
        .build()
        .unwrap();
    let sk = SecretKey::generate(&params, &mut rng);
    let rk = RelinearizationKey::new(&sk, &mut rng)?;

    let val1 = BigUint::from(10u32);
    let val2 = params.plaintext_modulus() - BigUint::from(20u32);

    let mut vec1 = vec![BigUint::from(0u32); params.degree()];
    vec1[0] = val1.clone();

    let mut vec2 = vec![BigUint::from(0u32); params.degree()];
    vec2[0] = val2.clone();

    let pt1 = Plaintext::encode_biguint(&params, vec1.as_slice(), Encoding::Polynomial)?;
    let pt2 = Plaintext::encode_biguint(&params, vec2.as_slice(), Encoding::Polynomial)?;

    let ct1: Ciphertext = sk.encrypt(&pt1, &mut rng)?;
    let ct2: Ciphertext = sk.encrypt(&pt2, &mut rng)?;

    let mut ct_res = ct1.multiply(&ct2).unwrap();
    rk.relinearize(&mut ct_res)?;

    assert_eq!(ct_res.components().len(), 2); // Degree reduced

    let decrypted_pt = sk.decrypt(&ct_res)?;
    let decrypted_values: Vec<BigUint> = decrypted_pt.decode_biguint(Encoding::Polynomial)?;

    // 10 * (-20) = -200
    assert_eq!(
        decrypted_values[0],
        params.plaintext_modulus() - BigUint::from(200u32)
    );

    Ok(())
}

#[test]
fn test_small_modulus_with_biguint_input() -> Result<(), Box<dyn Error>> {
    let mut rng = rng();
    // Standard small modulus parameters
    let params = ParametersBuilder::new()
        .degree(16)
        .plaintext_modulus(1153_u64)
        .ciphertext_modulus_bits([62usize; 1])
        .build()
        .unwrap();
    let sk = SecretKey::generate(&params, &mut rng);

    // Let's just pick a value larger than t, but small enough to verify reduction.
    // t = 1153 (default for test_parameters(1, 16) in parameters.rs)
    let t = params.plaintext_modulus_u64().unwrap();
    let val = BigUint::from(t) + 5u32; // Should reduce to 5

    let mut values = vec![BigUint::from(0u32); params.degree()];
    values[0] = val.clone();

    let pt = Plaintext::encode_biguint(&params, values.as_slice(), Encoding::Polynomial)?;
    let ct: Ciphertext = sk.encrypt(&pt, &mut rng)?;
    let decrypted_pt = sk.decrypt(&ct)?;

    let decrypted_values: Vec<u64> = decrypted_pt.decode(Encoding::Polynomial)?;

    assert_eq!(decrypted_values[0], 5);

    Ok(())
}
