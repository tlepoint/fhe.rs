// Expect indexing in examples for simplicity
#![expect(missing_docs, reason = "examples/benches/tests omit docs by design")]
#![expect(
    clippy::indexing_slicing,
    reason = "performance or example code relies on validated indices"
)]

mod util;

use std::error::Error;

use fhe::bfv::{
    Ciphertext, Encoding, Parameters, Plaintext, PublicKey, SecretKey,
    evaluation::EvaluationKeyBuilder, evaluation::RelinearizationKey,
};

use rand::rng;
use util::timeit::timeit;

fn weighted_sum_plain(
    cts: &[Ciphertext],
    weights: &[u64],
    params: &Parameters,
    sk: &SecretKey,
) -> Result<u64, Box<dyn Error>> {
    let mut acc = Ciphertext::trivial_zero(params, 0)?;
    for (ct, w) in cts.iter().zip(weights.iter()) {
        let pt_w = Plaintext::encode(params, &[*w], Encoding::Polynomial)?;
        acc.add_assign(&(ct.multiply_plaintext(&pt_w)?))?;
    }
    let pt = sk.decrypt(&acc)?;
    let v = pt.decode(Encoding::Polynomial)?;
    Ok(v[0])
}

fn weighted_sum_simd(
    ct: &Ciphertext,
    weights: &Plaintext,
    ek: &fhe::bfv::evaluation::EvaluationKey,
    sk: &SecretKey,
) -> Result<u64, Box<dyn Error>> {
    let tmp = ct.multiply_plaintext(weights)?;
    let summed = ek.inner_sum(&tmp)?;
    let pt = sk.decrypt(&summed)?;
    let v = pt.decode(Encoding::Simd)?;
    Ok(v[0])
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut rng = rng();
    let params = Parameters::profile_128(4096, 20)?;
    let sk = SecretKey::generate(&params, &mut rng);
    let pk = PublicKey::from_secret_key(&sk, &mut rng);
    let ek = EvaluationKeyBuilder::new(&sk)
        .ciphertext_level(0)
        .key_level(0)
        .enable_inner_sum()
        .build(&mut rng)?;
    let rk = RelinearizationKey::new(&sk, &mut rng)?;

    // ----- Weighted sum without SIMD -----
    let values = [1u64, 2, 3];
    let weights = [4u64, 5, 6];
    timeit!("inner product (no SIMD)", {
        let cts: Vec<Ciphertext> = values
            .iter()
            .map(|v| {
                let pt = Plaintext::encode(&params, &[*v], Encoding::Polynomial)?;
                Ok(pk.encrypt(&pt, &mut rng)?)
            })
            .collect::<Result<_, Box<dyn Error>>>()?;
        let ws_plain = weighted_sum_plain(&cts, &weights, &params, &sk)?;
        println!("Weighted sum (no SIMD) = {ws_plain}");
    });

    // ----- Weighted sum with SIMD -----
    let pt_vals = Plaintext::encode(&params, &values, Encoding::Simd)?;
    let ct_vals = pk.encrypt(&pt_vals, &mut rng)?;
    let pt_ws = Plaintext::encode(&params, &weights, Encoding::Simd)?;
    timeit!("inner product (SIMD)", {
        let ws_simd = weighted_sum_simd(&ct_vals, &pt_ws, &ek, &sk)?;
        println!("Weighted sum (SIMD) = {ws_simd}");
    });

    // ----- Inner product without SIMD -----
    let v1 = [1u64, 2, 3];
    let v2 = [7u64, 8, 9];
    let ct_v1: Vec<Ciphertext> = v1
        .iter()
        .map(|v| {
            let pt = Plaintext::encode(&params, &[*v], Encoding::Polynomial)?;
            Ok(pk.encrypt(&pt, &mut rng)?)
        })
        .collect::<Result<_, Box<dyn Error>>>()?;
    let ct_v2: Vec<Ciphertext> = v2
        .iter()
        .map(|v| {
            let pt = Plaintext::encode(&params, &[*v], Encoding::Polynomial)?;
            Ok(pk.encrypt(&pt, &mut rng)?)
        })
        .collect::<Result<_, Box<dyn Error>>>()?;
    let mut acc = Ciphertext::trivial_zero(&params, 0)?;
    for (a, b) in ct_v1.iter().zip(ct_v2.iter()) {
        let mut prod = a.multiply(b)?;
        rk.relinearize(&mut prod)?;
        acc.add_assign(&prod)?;
    }
    let pt = sk.decrypt(&acc)?;
    let ip_plain = pt.decode(Encoding::Polynomial)?[0];
    println!("Inner product (no SIMD) = {ip_plain}");

    // ----- Inner product with SIMD -----
    let pt1 = Plaintext::encode(&params, &v1, Encoding::Simd)?;
    let pt2 = Plaintext::encode(&params, &v2, Encoding::Simd)?;
    let ct1 = pk.encrypt(&pt1, &mut rng)?;
    let ct2 = pk.encrypt(&pt2, &mut rng)?;
    let mut prod = ct1.multiply(&ct2)?;
    rk.relinearize(&mut prod)?;
    let summed = ek.inner_sum(&prod)?;
    let pt = sk.decrypt(&summed)?;
    let ip_simd = pt.decode(Encoding::Simd)?[0];
    println!("Inner product (SIMD) = {ip_simd}");

    // ----- Polynomial evaluation without SIMD -----
    let x = 3u64;
    let pt_x = Plaintext::encode(&params, &[x], Encoding::Polynomial)?;
    let ct_x = pk.encrypt(&pt_x, &mut rng)?;
    let mut ct_x2 = ct_x.multiply(&ct_x)?; // x^2
    rk.relinearize(&mut ct_x2)?;
    let pt_three = Plaintext::encode(&params, &[3u64], Encoding::Polynomial)?;
    let pt_two = Plaintext::encode(&params, &[2u64], Encoding::Polynomial)?;
    let pt_one = Plaintext::encode(&params, &[1u64], Encoding::Polynomial)?;
    let mut ct_res = ct_x2.multiply_plaintext(&pt_three)?;
    ct_res.add_assign(&(ct_x.multiply_plaintext(&pt_two)?))?;
    ct_res.add_plaintext_assign(&pt_one)?;
    let pt = sk.decrypt(&ct_res)?;
    let poly_plain = pt.decode(Encoding::Polynomial)?[0];
    println!("Polynomial (no SIMD) = {poly_plain}");

    // ----- Polynomial evaluation with SIMD -----
    let x_vec = [1u64, 2, 3, 4];
    let pt_xv = Plaintext::encode(&params, &x_vec, Encoding::Simd)?;
    let ct_xv = pk.encrypt(&pt_xv, &mut rng)?;
    let mut ct_xv2 = ct_xv.multiply(&ct_xv)?;
    rk.relinearize(&mut ct_xv2)?;
    let pt_three_v = Plaintext::encode(&params, &vec![3u64; x_vec.len()], Encoding::Simd)?;
    let pt_two_v = Plaintext::encode(&params, &vec![2u64; x_vec.len()], Encoding::Simd)?;
    let pt_one_v = Plaintext::encode(&params, &vec![1u64; x_vec.len()], Encoding::Simd)?;
    let mut ct_res_v = ct_xv2.multiply_plaintext(&pt_three_v)?;
    ct_res_v.add_assign(&(ct_xv.multiply_plaintext(&pt_two_v)?))?;
    ct_res_v.add_plaintext_assign(&pt_one_v)?;
    let pt = sk.decrypt(&ct_res_v)?;
    let poly_simd = pt.decode(Encoding::Simd)?;
    println!("Polynomial (SIMD) = {:?}", &poly_simd[..x_vec.len()]);

    Ok(())
}
