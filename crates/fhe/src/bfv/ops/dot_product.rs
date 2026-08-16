use fhe_math::rq::{Ntt, Poly, dot_product as poly_dot_product, traits::TryConvertFrom};
use itertools::{Itertools, izip};
use ndarray::{Array, Array2};

use crate::{
    Error, Result,
    bfv::{Ciphertext, Plaintext},
};

/// Computes the Fused-Mul-Add operation `out[i] += x[i] * y[i]`
unsafe fn fma(out: &mut [u128], x: &[u64], y: &[u64]) {
    let n = out.len();
    debug_assert_eq!(x.len(), n);
    debug_assert_eq!(y.len(), n);

    macro_rules! fma_at {
        ($idx:expr) => {
            unsafe {
                *out.get_unchecked_mut($idx) +=
                    (*x.get_unchecked($idx) as u128) * (*y.get_unchecked($idx) as u128);
            }
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
/// iterator of [`Plaintext`]. Returns an error if either iterator is empty, the
/// iterator lengths differ, the parameters or levels do not match, or the
/// ciphertexts have different numbers of parts.
pub fn dot_product_scalar<'a, I, J>(ct: I, pt: J) -> Result<Ciphertext>
where
    I: Iterator<Item = &'a Ciphertext> + Clone,
    J: Iterator<Item = &'a Plaintext> + Clone,
{
    let ct_count = ct.clone().count();
    let pt_count = pt.clone().count();
    if ct_count == 0 || pt_count == 0 {
        return Err(crate::DotProductError::EmptyInput.into());
    }
    if ct_count != pt_count {
        return Err(crate::DotProductError::OperandCountMismatch {
            ciphertexts: ct_count,
            plaintexts: pt_count,
        }
        .into());
    }
    let count = ct_count;
    let ct_first = ct
        .clone()
        .next()
        .ok_or(crate::DotProductError::EmptyInput)?;
    ct_first.validate_for(&ct_first.par)?;
    let ctx = ct_first.par.context_at_level(ct_first.level)?;

    // Variable-time reductions are permitted only when every ciphertext and
    // plaintext polynomial in the dot product has been classified as public.
    let allow_variable_time_computations =
        ct.clone()
            .zip(pt.clone())
            .take(count)
            .all(|(ciphertext, plaintext)| {
                ciphertext
                    .iter()
                    .all(Poly::allows_variable_time_computations)
                    && plaintext.poly_ntt.allows_variable_time_computations()
            });

    for (cti, pti) in ct.clone().zip(pt.clone()) {
        cti.validate_for_context(&ct_first.par, ct_first.level, ctx)?;
        pti.validate_for_context(&ct_first.par, ct_first.level, ctx)?;
        if cti.len() != ct_first.len() {
            return Err(crate::DotProductError::CiphertextPolynomialCountMismatch {
                actual: cti.len(),
                expected: ct_first.len(),
            }
            .into());
        }
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
            .collect::<Result<Vec<Poly<Ntt>>>>()?;

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
                    if allow_variable_time_computations {
                        unsafe { *outij_coeff = q.reduce_u128_vt(*accij_coeff) }
                    } else {
                        *outij_coeff = q.reduce_u128(*accij_coeff)
                    }
                }
            }
            c.push(Poly::<Ntt>::try_convert_from(
                coeffs,
                ctx,
                allow_variable_time_computations,
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
    use itertools::{Itertools, izip};
    use rand::rng;
    use std::error::Error;

    #[test]
    fn test_dot_product_scalar() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let empty_ct: Vec<Ciphertext> = Vec::new();
        let empty_pt: Vec<Plaintext> = Vec::new();
        assert!(dot_product_scalar(empty_ct.iter(), empty_pt.iter()).is_err());

        for params in [
            BfvParameters::default_arc(1, 16),
            BfvParameters::default_arc(2, 32),
        ] {
            let sk = SecretKey::random(&params, &mut rng);
            for size in 1..128 {
                let ct = (0..size)
                    .map(|_| {
                        let v = fhe_math::zq::Modulus::new(params.plaintext())
                            .unwrap()
                            .random_vec(params.degree(), &mut rng);
                        let pt = Plaintext::try_encode(&v, Encoding::simd(), &params).unwrap();
                        sk.try_encrypt(&pt, &mut rng).unwrap()
                    })
                    .collect_vec();
                let pt = (0..size)
                    .map(|_| {
                        let v = fhe_math::zq::Modulus::new(params.plaintext())
                            .unwrap()
                            .random_vec(params.degree(), &mut rng);
                        Plaintext::try_encode(&v, Encoding::simd(), &params).unwrap()
                    })
                    .collect_vec();

                let r = dot_product_scalar(ct.iter(), pt.iter())?;
                assert!(
                    r.iter()
                        .all(|poly| !poly.allows_variable_time_computations())
                );

                let mut expected = Ciphertext::zero(&params);
                izip!(&ct, &pt).for_each(|(cti, pti)| expected += &(cti * pti));
                assert_eq!(r, expected);

                let variable_time =
                    fhe_traits::VariableTime::new(fhe_traits::PublicData::assert_public());
                let mut public_pt = pt.clone();
                public_pt.iter_mut().for_each(|plaintext| {
                    plaintext
                        .poly_ntt
                        .allow_variable_time_computations(variable_time)
                });
                let public_result = dot_product_scalar(ct.iter(), public_pt.iter())?;
                assert!(
                    public_result
                        .iter()
                        .all(|poly| poly.allows_variable_time_computations())
                );
            }
        }
        Ok(())
    }

    #[test]
    fn dot_product_scalar_rejects_mismatched_inputs() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let params = BfvParameters::default_arc(1, 16);
        let sk = SecretKey::random(&params, &mut rng);
        let pt = Plaintext::try_encode(&[1u64][..], Encoding::poly(), &params)?;
        let ct = sk.try_encrypt(&pt, &mut rng)?;

        assert!(matches!(
            dot_product_scalar([&ct].into_iter(), [&pt, &pt].into_iter()),
            Err(crate::Error::DotProduct(
                crate::DotProductError::OperandCountMismatch { .. }
            ))
        ));
        assert!(matches!(
            dot_product_scalar([&Ciphertext::zero(&params)].into_iter(), [&pt].into_iter()),
            Err(crate::Error::Ciphertext(_))
        ));
        Ok(())
    }
}
