#![allow(missing_docs, clippy::indexing_slicing)]
use fhe::bfv::{
    BfvParameters, BfvParametersBuilder, Ciphertext, Encoding, Plaintext, RelinearizationKey,
    SecretKey,
};
use fhe_traits::{FheDecoder, FheDecrypter, FheEncoder as _, FheEncrypter};
use num_bigint::BigUint;
use rand::rng;
use std::{error::Error, sync::Arc};

fn parameters() -> Arc<BfvParameters> {
    // Choose a large plaintext modulus: 2^127 - 1 (Mersenne prime M127)
    // 170141183460469231731687303715884105727
    let p_str = "170141183460469231731687303715884105727";
    let p = BigUint::parse_bytes(p_str.as_bytes(), 10).unwrap();

    // Create parameters
    BfvParametersBuilder::new()
        .set_degree(16)
        .set_plaintext_modulus_biguint(p.clone())
        .set_moduli_sizes(&[60, 60, 60, 60, 60])
        .build_arc()
        .unwrap()
}

#[test]
fn test_biguint_plaintext_encryption_decryption() -> Result<(), Box<dyn Error>> {
    let mut rng = rng();

    let params = parameters();
    let sk = SecretKey::random(&params, &mut rng);

    // Create a vector of BigUint values
    let mut values = vec![BigUint::from(0u32); params.degree()];
    values[0] = BigUint::from(123456789u64);
    values[1] = params.plaintext_big() - 1u32; // -1
    values[2] = params.plaintext_big() / 2u32;

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

    let params = parameters();
    let sk = SecretKey::random(&params, &mut rng);

    let val1 = BigUint::from(10u32);
    let val2 = params.plaintext_big() - 50u32; // -50

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

    // 10 + (-50) = -40
    assert_eq!(
        decrypted_values[0],
        params.plaintext_big() - BigUint::from(40u32)
    );

    Ok(())
}

#[test]
fn test_biguint_multiplication_without_relin() -> Result<(), Box<dyn Error>> {
    let mut rng = rng();

    let params = parameters();
    let sk = SecretKey::random(&params, &mut rng);

    let val1 = BigUint::from(10u32);
    let val2 = params.plaintext_big() - BigUint::from(20u32);

    let mut vec1 = vec![BigUint::from(0u32); params.degree()];
    vec1[0] = val1.clone();

    let mut vec2 = vec![BigUint::from(0u32); params.degree()];
    vec2[0] = val2.clone();

    let pt1 = Plaintext::try_encode(vec1.as_slice(), Encoding::poly(), &params)?;
    let pt2 = Plaintext::try_encode(vec2.as_slice(), Encoding::poly(), &params)?;

    let ct1: Ciphertext = sk.try_encrypt(&pt1, &mut rng)?;
    let ct2: Ciphertext = sk.try_encrypt(&pt2, &mut rng)?;

    let ct_res = &ct1 * &ct2;

    assert_eq!(ct_res.len(), 3); // Degree increases

    let decrypted_pt = sk.try_decrypt(&ct_res)?;
    let decrypted_values: Vec<BigUint> =
        Vec::<BigUint>::try_decode(&decrypted_pt, Encoding::poly())?;

    // 10 * (-20) = -200
    assert_eq!(
        decrypted_values[0],
        params.plaintext_big() - BigUint::from(200u32)
    );

    Ok(())
}

#[test]
fn test_biguint_multiplication_with_relin() -> Result<(), Box<dyn Error>> {
    let mut rng = rng();

    // Use default parameters with biguint
    let params = BfvParametersBuilder::new()
        .set_degree(16)
        .set_plaintext_modulus(1153)
        .set_moduli_sizes(&[62usize; 3])
        .build_arc()
        .unwrap();
    let sk = SecretKey::random(&params, &mut rng);
    let rk = RelinearizationKey::new(&sk, &mut rng)?;

    let val1 = BigUint::from(10u32);
    let val2 = params.plaintext_big() - BigUint::from(20u32);

    let mut vec1 = vec![BigUint::from(0u32); params.degree()];
    vec1[0] = val1.clone();

    let mut vec2 = vec![BigUint::from(0u32); params.degree()];
    vec2[0] = val2.clone();

    let pt1 = Plaintext::try_encode(vec1.as_slice(), Encoding::poly(), &params)?;
    let pt2 = Plaintext::try_encode(vec2.as_slice(), Encoding::poly(), &params)?;

    let ct1: Ciphertext = sk.try_encrypt(&pt1, &mut rng)?;
    let ct2: Ciphertext = sk.try_encrypt(&pt2, &mut rng)?;

    let mut ct_res = &ct1 * &ct2;
    rk.relinearizes(&mut ct_res)?;

    assert_eq!(ct_res.len(), 2); // Degree reduced

    let decrypted_pt = sk.try_decrypt(&ct_res)?;
    let decrypted_values: Vec<BigUint> =
        Vec::<BigUint>::try_decode(&decrypted_pt, Encoding::poly())?;

    // 10 * (-20) = -200
    assert_eq!(
        decrypted_values[0],
        params.plaintext_big() - BigUint::from(200u32)
    );

    Ok(())
}

#[test]
fn test_small_modulus_with_biguint_input() -> Result<(), Box<dyn Error>> {
    let mut rng = rng();
    // Standard small modulus parameters
    let params = BfvParametersBuilder::new()
        .set_degree(16)
        .set_plaintext_modulus(1153)
        .set_moduli_sizes(&[62usize; 1])
        .build_arc()
        .unwrap();
    let sk = SecretKey::random(&params, &mut rng);

    // Let's just pick a value larger than t, but small enough to verify reduction.
    // t = 1153 (default for default_arc(1, 16) in parameters.rs)
    let t = params.plaintext();
    let val = BigUint::from(t) + 5u32; // Should reduce to 5

    let mut values = vec![BigUint::from(0u32); params.degree()];
    values[0] = val.clone();

    let pt = Plaintext::try_encode(values.as_slice(), Encoding::poly(), &params)?;
    let ct: Ciphertext = sk.try_encrypt(&pt, &mut rng)?;
    let decrypted_pt = sk.try_decrypt(&ct)?;

    let decrypted_values: Vec<u64> = Vec::<u64>::try_decode(&decrypted_pt, Encoding::poly())?;

    assert_eq!(decrypted_values[0], 5);

    Ok(())
}
