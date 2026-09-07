// Expect indexing in examples for simplicity
#![expect(missing_docs, reason = "examples/benches/tests omit docs by design")]
#![expect(
    clippy::indexing_slicing,
    reason = "performance or example code relies on validated indices"
)]

use std::error::Error;

use fhe::bfv::{
    Ciphertext, Encoding, Parameters, Plaintext, SecretKey, evaluation::RgswCiphertext,
};

use rand::rng;

fn main() -> Result<(), Box<dyn Error>> {
    let mut rng = rng();
    let params = Parameters::profile_128(4096, 20)?;
    let sk = SecretKey::generate(&params, &mut rng);

    let v1 = vec![1u64, 2, 3, 4];
    let v2 = vec![5u64, 6, 7, 8];
    let pt1 = Plaintext::encode(&params, &v1, Encoding::Simd)?;
    let pt2 = Plaintext::encode(&params, &v2, Encoding::Simd)?;
    let ct1: Ciphertext = sk.encrypt(&pt1, &mut rng)?;
    let ct2: Ciphertext = sk.encrypt(&pt2, &mut rng)?;
    let ct2_rgsw: RgswCiphertext = sk.encrypt_rgsw(&pt2, &mut rng)?;

    let mut product = ct1.multiply_rgsw(&ct2_rgsw)?;
    let expected = ct1.multiply(&ct2)?;

    println!(
        "Noise in product: {}",
        sk.measure_noise_vartime(
            &product,
            fhe::SecretDependentDiagnostics::acknowledge_leakage()
        )?
    );
    println!("Size of product: {} bytes", product.to_bytes().len());
    println!(
        "Noise in expected: {}",
        sk.measure_noise_vartime(
            &expected,
            fhe::SecretDependentDiagnostics::acknowledge_leakage()
        )?
    );

    product.switch_to_level(product.max_switchable_level())?;
    println!(
        "Noise in product: {}",
        sk.measure_noise_vartime(
            &product,
            fhe::SecretDependentDiagnostics::acknowledge_leakage()
        )?
    );
    println!("Size of product: {} bytes", product.to_bytes().len());

    let pt_prod = sk.decrypt(&product)?;
    let pt_exp = sk.decrypt(&expected)?;
    let decoded = pt_prod.decode(Encoding::Simd)?;
    // These plaintexts have different levels; compare the decoded messages.
    assert_eq!(decoded, pt_exp.decode(Encoding::Simd)?);
    println!(
        "RGSW external product successful: {:?}",
        &decoded[..v1.len()]
    );

    Ok(())
}
