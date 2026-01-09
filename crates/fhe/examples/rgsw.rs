// Allow indexing in examples for simplicity
#![allow(clippy::indexing_slicing)]
#![allow(missing_docs)]

use std::error::Error;

use fhe::bfv::{BfvParameters, Ciphertext, Encoding, Plaintext, RGSWCiphertext, SecretKey};
use fhe_traits::{FheDecoder, FheDecrypter, FheEncoder, FheEncrypter, Serialize};
use rand::rng;

fn main() -> Result<(), Box<dyn Error>> {
    let mut rng = rng();
    let params = BfvParameters::default_parameters_128(20)
        .unwrap()
        .nth(2)
        .unwrap();
    let sk = SecretKey::random(&params, &mut rng);

    let v1 = vec![1u64, 2, 3, 4];
    let v2 = vec![5u64, 6, 7, 8];
    let pt1 = Plaintext::try_encode(&v1, Encoding::simd(), &params)?;
    let pt2 = Plaintext::try_encode(&v2, Encoding::simd(), &params)?;
    let ct1: Ciphertext = sk.try_encrypt(&pt1, &mut rng)?;
    let ct2: Ciphertext = sk.try_encrypt(&pt2, &mut rng)?;
    let ct2_rgsw: RGSWCiphertext = sk.try_encrypt(&pt2, &mut rng)?;

    let mut product = &ct1 * &ct2_rgsw;
    let expected = &ct1 * &ct2;

    println!("Noise in product: {}", unsafe {
        sk.measure_noise(&product)?
    });
    println!("Size of product: {} bytes", product.to_bytes().len());
    println!("Noise in expected: {}", unsafe {
        sk.measure_noise(&product)?
    });

    product.switch_to_level(product.max_switchable_level())?;
    println!("Noise in product: {}", unsafe {
        sk.measure_noise(&product)?
    });
    println!("Size of product: {} bytes", product.to_bytes().len());

    let pt_prod = sk.try_decrypt(&product)?;
    let pt_exp = sk.try_decrypt(&expected)?;
    assert_eq!(pt_prod, pt_exp);
    let decoded = Vec::<u64>::try_decode(&pt_prod, Encoding::simd())?;
    println!(
        "RGSW external product successful: {:?}",
        &decoded[..v1.len()]
    );

    Ok(())
}
