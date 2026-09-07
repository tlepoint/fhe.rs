// Expect indexing in examples for simplicity
#![expect(missing_docs, reason = "examples/benches/tests omit docs by design")]
#![expect(
    clippy::indexing_slicing,
    reason = "performance or example code relies on validated indices"
)]

use std::error::Error;

use fhe::bfv::{Encoding, Parameters, Plaintext, PublicKey, SecretKey};

use rand::rng;

fn main() -> Result<(), Box<dyn Error>> {
    let mut rng = rng();
    // Use default parameters
    let params = Parameters::profile_128(4096, 16)?;

    // Generate keys
    let sk = SecretKey::generate(&params, &mut rng);
    let pk = PublicKey::from_secret_key(&sk, &mut rng);

    // ----- Without SIMD -----
    let pt_a = Plaintext::encode(&params, &[3u64], Encoding::Polynomial)?;
    let pt_b = Plaintext::encode(&params, &[5u64], Encoding::Polynomial)?;
    let ct_a = pk.encrypt(&pt_a, &mut rng)?;
    let ct_b = pk.encrypt(&pt_b, &mut rng)?;
    let ct_sum = ct_a.add(&ct_b)?;
    let pt_sum = sk.decrypt(&ct_sum)?;
    let res = pt_sum.decode(Encoding::Polynomial)?;
    println!("3 + 5 = {}", res[0]);

    // ----- With SIMD -----
    let v1 = vec![1u64, 2, 3, 4];
    let v2 = vec![5u64, 6, 7, 8];
    let pt_v1 = Plaintext::encode(&params, &v1, Encoding::Simd)?;
    let pt_v2 = Plaintext::encode(&params, &v2, Encoding::Simd)?;
    let ct_v1 = pk.encrypt(&pt_v1, &mut rng)?;
    let ct_v2 = pk.encrypt(&pt_v2, &mut rng)?;
    let ct_vsum = ct_v1.add(&ct_v2)?;
    let pt_vsum = sk.decrypt(&ct_vsum)?;
    let res_v = pt_vsum.decode(Encoding::Simd)?;
    println!("{:?} + {:?} = {:?}", v1, v2, &res_v[..v1.len()]);

    Ok(())
}
