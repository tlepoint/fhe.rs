// Allow indexing in examples for simplicity
#![allow(clippy::indexing_slicing)]
#![allow(missing_docs)]

use std::error::Error;

use fhe::bfv::{BfvParameters, Encoding, Plaintext, PublicKey, SecretKey};
use fhe_traits::{FheDecoder, FheDecrypter, FheEncoder, FheEncrypter};
use rand::rng;

fn main() -> Result<(), Box<dyn Error>> {
    let mut rng = rng();
    // Use default parameters
    let params = BfvParameters::default_parameters_128(16)?
        .nth(2)
        .ok_or("Could not generate parameters")?;

    // Generate keys
    let sk = SecretKey::random(&params, &mut rng);
    let pk = PublicKey::new(&sk, &mut rng);

    // ----- Without SIMD -----
    let pt_a = Plaintext::try_encode(&[3u64], Encoding::poly(), &params)?;
    let pt_b = Plaintext::try_encode(&[5u64], Encoding::poly(), &params)?;
    let ct_a = pk.try_encrypt(&pt_a, &mut rng)?;
    let ct_b = pk.try_encrypt(&pt_b, &mut rng)?;
    let ct_sum = &ct_a + &ct_b;
    let pt_sum = sk.try_decrypt(&ct_sum)?;
    let res = Vec::<u64>::try_decode(&pt_sum, Encoding::poly())?;
    println!("3 + 5 = {}", res[0]);

    // ----- With SIMD -----
    let v1 = vec![1u64, 2, 3, 4];
    let v2 = vec![5u64, 6, 7, 8];
    let pt_v1 = Plaintext::try_encode(&v1, Encoding::simd(), &params)?;
    let pt_v2 = Plaintext::try_encode(&v2, Encoding::simd(), &params)?;
    let ct_v1 = pk.try_encrypt(&pt_v1, &mut rng)?;
    let ct_v2 = pk.try_encrypt(&pt_v2, &mut rng)?;
    let ct_vsum = &ct_v1 + &ct_v2;
    let pt_vsum = sk.try_decrypt(&ct_vsum)?;
    let res_v = Vec::<u64>::try_decode(&pt_vsum, Encoding::simd())?;
    println!("{:?} + {:?} = {:?}", v1, v2, &res_v[..v1.len()]);

    Ok(())
}
