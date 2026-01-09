// Allow indexing in examples for simplicity
#![allow(clippy::indexing_slicing)]
#![allow(missing_docs)]

mod util;

use std::error::Error;
use std::sync::Arc;

use fhe::bfv::{
    BfvParameters, Ciphertext, Encoding, EvaluationKeyBuilder, Plaintext, PublicKey,
    RelinearizationKey, SecretKey,
};
use fhe_traits::{FheDecoder, FheDecrypter, FheEncoder, FheEncrypter};
use rand::rng;
use util::timeit::timeit;

fn weighted_sum_plain(
    cts: &[Ciphertext],
    weights: &[u64],
    params: &Arc<BfvParameters>,
    sk: &SecretKey,
) -> Result<u64, Box<dyn Error>> {
    let mut acc = Ciphertext::zero(params);
    for (ct, w) in cts.iter().zip(weights.iter()) {
        let pt_w = Plaintext::try_encode(&[*w], Encoding::poly(), params)?;
        acc += &(ct * &pt_w);
    }
    let pt = sk.try_decrypt(&acc)?;
    let v = Vec::<u64>::try_decode(&pt, Encoding::poly())?;
    Ok(v[0])
}

fn weighted_sum_simd(
    ct: &Ciphertext,
    weights: &Plaintext,
    ek: &fhe::bfv::EvaluationKey,
    sk: &SecretKey,
) -> Result<u64, Box<dyn Error>> {
    let tmp = ct * weights;
    let summed = ek.computes_inner_sum(&tmp)?;
    let pt = sk.try_decrypt(&summed)?;
    let v = Vec::<u64>::try_decode(&pt, Encoding::simd())?;
    Ok(v[0])
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut rng = rng();
    let params = BfvParameters::default_parameters_128(20)
        .unwrap()
        .nth(2) // first parameters do not support key switching
        .unwrap();
    let sk = SecretKey::random(&params, &mut rng);
    let pk = PublicKey::new(&sk, &mut rng);
    let ek = EvaluationKeyBuilder::new_leveled(&sk, 0, 0)?
        .enable_inner_sum()?
        .build(&mut rng)?;
    let rk = RelinearizationKey::new(&sk, &mut rng)?;

    // ----- Weighted sum without SIMD -----
    let values = [1u64, 2, 3];
    let weights = [4u64, 5, 6];
    timeit!("inner product (no SIMD)", {
        let cts: Vec<Ciphertext> = values
            .iter()
            .map(|v| {
                let pt = Plaintext::try_encode(&[*v], Encoding::poly(), &params)?;
                Ok(pk.try_encrypt(&pt, &mut rng)?)
            })
            .collect::<Result<_, Box<dyn Error>>>()?;
        let ws_plain = weighted_sum_plain(&cts, &weights, &params, &sk)?;
        println!("Weighted sum (no SIMD) = {ws_plain}");
    });

    // ----- Weighted sum with SIMD -----
    let pt_vals = Plaintext::try_encode(&values, Encoding::simd(), &params)?;
    let ct_vals = pk.try_encrypt(&pt_vals, &mut rng)?;
    let pt_ws = Plaintext::try_encode(&weights, Encoding::simd(), &params)?;
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
            let pt = Plaintext::try_encode(&[*v], Encoding::poly(), &params)?;
            Ok(pk.try_encrypt(&pt, &mut rng)?)
        })
        .collect::<Result<_, Box<dyn Error>>>()?;
    let ct_v2: Vec<Ciphertext> = v2
        .iter()
        .map(|v| {
            let pt = Plaintext::try_encode(&[*v], Encoding::poly(), &params)?;
            Ok(pk.try_encrypt(&pt, &mut rng)?)
        })
        .collect::<Result<_, Box<dyn Error>>>()?;
    let mut acc = Ciphertext::zero(&params);
    for (a, b) in ct_v1.iter().zip(ct_v2.iter()) {
        let mut prod = a * b;
        rk.relinearizes(&mut prod)?;
        acc += &prod;
    }
    let pt = sk.try_decrypt(&acc)?;
    let ip_plain = Vec::<u64>::try_decode(&pt, Encoding::poly())?[0];
    println!("Inner product (no SIMD) = {ip_plain}");

    // ----- Inner product with SIMD -----
    let pt1 = Plaintext::try_encode(&v1, Encoding::simd(), &params)?;
    let pt2 = Plaintext::try_encode(&v2, Encoding::simd(), &params)?;
    let ct1 = pk.try_encrypt(&pt1, &mut rng)?;
    let ct2 = pk.try_encrypt(&pt2, &mut rng)?;
    let mut prod = &ct1 * &ct2;
    rk.relinearizes(&mut prod)?;
    let summed = ek.computes_inner_sum(&prod)?;
    let pt = sk.try_decrypt(&summed)?;
    let ip_simd = Vec::<u64>::try_decode(&pt, Encoding::simd())?[0];
    println!("Inner product (SIMD) = {ip_simd}");

    // ----- Polynomial evaluation without SIMD -----
    let x = 3u64;
    let pt_x = Plaintext::try_encode(&[x], Encoding::poly(), &params)?;
    let ct_x = pk.try_encrypt(&pt_x, &mut rng)?;
    let mut ct_x2 = &ct_x * &ct_x; // x^2
    rk.relinearizes(&mut ct_x2)?;
    let pt_three = Plaintext::try_encode(&[3u64], Encoding::poly(), &params)?;
    let pt_two = Plaintext::try_encode(&[2u64], Encoding::poly(), &params)?;
    let pt_one = Plaintext::try_encode(&[1u64], Encoding::poly(), &params)?;
    let mut ct_res = &ct_x2 * &pt_three;
    ct_res += &(&ct_x * &pt_two);
    ct_res += &pt_one;
    let pt = sk.try_decrypt(&ct_res)?;
    let poly_plain = Vec::<u64>::try_decode(&pt, Encoding::poly())?[0];
    println!("Polynomial (no SIMD) = {poly_plain}");

    // ----- Polynomial evaluation with SIMD -----
    let x_vec = [1u64, 2, 3, 4];
    let pt_xv = Plaintext::try_encode(&x_vec, Encoding::simd(), &params)?;
    let ct_xv = pk.try_encrypt(&pt_xv, &mut rng)?;
    let mut ct_xv2 = &ct_xv * &ct_xv;
    rk.relinearizes(&mut ct_xv2)?;
    let pt_three_v = Plaintext::try_encode(&vec![3u64; x_vec.len()], Encoding::simd(), &params)?;
    let pt_two_v = Plaintext::try_encode(&vec![2u64; x_vec.len()], Encoding::simd(), &params)?;
    let pt_one_v = Plaintext::try_encode(&vec![1u64; x_vec.len()], Encoding::simd(), &params)?;
    let mut ct_res_v = &ct_xv2 * &pt_three_v;
    ct_res_v += &(&ct_xv * &pt_two_v);
    ct_res_v += &pt_one_v;
    let pt = sk.try_decrypt(&ct_res_v)?;
    let poly_simd = Vec::<u64>::try_decode(&pt, Encoding::simd())?;
    println!("Polynomial (SIMD) = {:?}", &poly_simd[..x_vec.len()]);

    Ok(())
}
