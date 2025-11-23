use std::cmp::min;

use fhe_math::rq::{dot_product as poly_dot_product, traits::TryConvertFrom, Poly, Representation};
use itertools::{izip, Itertools};
use ndarray::{Array, Array2};

use crate::{
    bfv::{Ciphertext, Plaintext},
    Error, Result,
};

/// Computes the Fused-Mul-Add operation `out[i] += x[i] * y[i]`
unsafe fn fma(out: &mut [u128], x: &[u64], y: &[u64]) {
    let n = out.len();
    debug_assert_eq!(x.len(), n);
    debug_assert_eq!(y.len(), n);

    macro_rules! fma_at {
        ($idx:expr) => {
            *out.get_unchecked_mut($idx) +=
                (*x.get_unchecked($idx) as u128) * (*y.get_unchecked($idx) as u128);
        };
    }

    let r = n / 16;
    for i in 0..r {
        fma_at!(16 * i);
        fma_at!(16 * i + 1);
        fma_at!(16 * i + 2);
        fma_at!(16 * i + 3);
        fma_at!(16 * i + 4);
        fma_at!(16 * i + 5);
        fma_at!(16 * i + 6);
        fma_at!(16 * i + 7);
        fma_at!(16 * i + 8);
        fma_at!(16 * i + 9);
        fma_at!(16 * i + 10);
        fma_at!(16 * i + 11);
        fma_at!(16 * i + 12);
        fma_at!(16 * i + 13);
        fma_at!(16 * i + 14);
        fma_at!(16 * i + 15);
    }

    for i in 0..n % 16 {
        fma_at!(16 * r + i);
    }
}

/// Compute the dot product between an iterator of [`Ciphertext`] and an
/// iterator of [`Plaintext`]. Returns an error if the iterator counts are 0, if
/// the parameters don't match, or if the ciphertexts have different
/// number of parts.
pub fn dot_product_scalar<'a, I, J>(ct: I, pt: J) -> Result<Ciphertext>
where
    I: Iterator<Item = &'a Ciphertext> + Clone,
    J: Iterator<Item = &'a Plaintext> + Clone,
{
    let count = min(ct.clone().count(), pt.clone().count());
    if count == 0 {
        return Err(Error::DefaultError(
            "At least one iterator is empty".to_string(),
        ));
    }
    let ct_first = ct.clone().next().unwrap();
    let ctx = ct_first[0].ctx();

    if izip!(ct.clone(), pt.clone()).any(|(cti, pti)| {
        cti.par != ct_first.par || pti.par != ct_first.par || cti.len() != ct_first.len()
    }) {
        return Err(Error::DefaultError("Mismatched parameters".to_string()));
    }
    if ct.clone().any(|cti| cti.len() != ct_first.len()) {
        return Err(Error::DefaultError(
            "Mismatched number of parts in the ciphertexts".to_string(),
        ));
    }

    let max_acc = ctx
        .moduli()
        .iter()
        .map(|qi| 1u128 << (2 * qi.leading_zeros()))
        .collect_vec();
    let min_of_max = max_acc.iter().min().unwrap();

    if count as u128 > *min_of_max {
        // Too many ciphertexts for the optimized method, instead, we call
        // `poly_dot_product`.
        let c = (0..ct_first.len())
            .map(|i| {
                poly_dot_product(
                    ct.clone().map(|cti| unsafe { cti.get_unchecked(i) }),
                    pt.clone().map(|pti| &pti.poly_ntt),
                )
                .map_err(Error::MathError)
            })
            .collect::<Result<Vec<Poly>>>()?;

        Ok(Ciphertext {
            par: ct_first.par.clone(),
            seed: None,
            c,
            level: ct_first.level,
        })
    } else {
        let mut acc = Array::zeros((ct_first.len(), ctx.moduli().len(), ct_first.par.degree()));
        for (ciphertext, plaintext) in izip!(ct, pt) {
            let pt_coefficients = plaintext.poly_ntt.coefficients();
            for (mut acci, ci) in izip!(acc.outer_iter_mut(), ciphertext.iter()) {
                let ci_coefficients = ci.coefficients();
                for (mut accij, cij, pij) in izip!(
                    acci.outer_iter_mut(),
                    ci_coefficients.outer_iter(),
                    pt_coefficients.outer_iter()
                ) {
                    unsafe {
                        fma(
                            accij.as_slice_mut().unwrap(),
                            cij.as_slice().unwrap(),
                            pij.as_slice().unwrap(),
                        )
                    }
                }
            }
        }

        // Reduce
        let mut c = Vec::with_capacity(ct_first.len());
        for acci in acc.outer_iter() {
            let mut coeffs = Array2::zeros((ctx.moduli().len(), ct_first.par.degree()));
            for (mut outij, accij, q) in izip!(
                coeffs.outer_iter_mut(),
                acci.outer_iter(),
                ctx.moduli_operators()
            ) {
                for (outij_coeff, accij_coeff) in izip!(outij.iter_mut(), accij.iter()) {
                    unsafe { *outij_coeff = q.reduce_u128_vt(*accij_coeff) }
                }
            }
            c.push(Poly::try_convert_from(
                coeffs,
                ctx,
                true,
                Representation::Ntt,
            )?)
        }

        Ok(Ciphertext {
            par: ct_first.par.clone(),
            seed: None,
            c,
            level: ct_first.level,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::dot_product_scalar;
    use crate::bfv::{BfvParameters, Ciphertext, Encoding, Plaintext, SecretKey};
    use fhe_traits::{FheEncoder, FheEncrypter};
    use itertools::{izip, Itertools};
    use rand::rng;
    use std::error::Error;

    #[test]
    fn test_dot_product_scalar() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for params in [
            BfvParameters::default_arc(1, 16),
            BfvParameters::default_arc(2, 32),
        ] {
            let sk = SecretKey::random(&params, &mut rng);
            for size in 1..128 {
                let ct = (0..size)
                    .map(|_| {
                        let v = params.plaintext.random_vec(params.degree(), &mut rng);
                        let pt = Plaintext::try_encode(&v, Encoding::simd(), &params).unwrap();
                        sk.try_encrypt(&pt, &mut rng).unwrap()
                    })
                    .collect_vec();
                let pt = (0..size)
                    .map(|_| {
                        let v = params.plaintext.random_vec(params.degree(), &mut rng);
                        Plaintext::try_encode(&v, Encoding::simd(), &params).unwrap()
                    })
                    .collect_vec();

                let r = dot_product_scalar(ct.iter(), pt.iter())?;

                let mut expected = Ciphertext::zero(&params);
                izip!(&ct, &pt).for_each(|(cti, pti)| expected += &(cti * pti));
                assert_eq!(r, expected);
            }
        }
        Ok(())
    }
}
