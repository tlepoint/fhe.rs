use fhe_math::rq::{DotProductWorkspace, Ntt, Poly, traits::TryConvertFrom};
use itertools::izip;
use ndarray::Array3;
use std::sync::Arc;
use zeroize::Zeroize;

use crate::{
    Error, Result,
    bfv::{BfvParameters, Ciphertext, Plaintext},
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

mod packed;

/// Compute a ciphertext/plaintext dot product with temporary scratch.
/// For repeated calls at the same level, use [`DotProductScalarWorkspace`].
pub fn dot_product_scalar<'a, 'b, I, J>(ct: I, pt: J) -> Result<Ciphertext>
where
    I: Iterator<Item = &'a Ciphertext> + Clone,
    J: Iterator<Item = &'b Plaintext> + Clone,
{
    let first = ct
        .clone()
        .next()
        .ok_or(crate::DotProductError::EmptyInput)?;
    DotProductScalarWorkspace::new(&first.par, first.level)?.dot_product_scalar(ct, pt)
}

struct ClearAccumulator<'a>(&'a mut Array3<u128>);
impl Drop for ClearAccumulator<'_> {
    fn drop(&mut self) {
        self.0.as_slice_mut().unwrap().zeroize();
    }
}

/// Reusable scratch for BFV ciphertext/plaintext dot products at a fixed level.
///
/// Both the fused fast path and the long-product fallback retain their scratch
/// across calls. Ciphertext results own their buffers independently. Scratch is
/// cleared after each fast-path call, including unwinding; resizing or dropping
/// the workspace therefore releases only cleared accumulator data.
pub struct DotProductScalarWorkspace {
    par: Arc<BfvParameters>,
    level: usize,
    min_limit: u128,
    accumulator: Array3<u128>,
    fallback: Option<DotProductWorkspace>,
}

impl DotProductScalarWorkspace {
    /// Create a workspace bound to these parameters and level. Coefficient
    /// buffers are allocated on first use and reused for matching part counts.
    pub fn new(par: &Arc<BfvParameters>, level: usize) -> Result<Self> {
        let ctx = par.context_at_level(level)?;
        let min_limit = ctx
            .moduli()
            .iter()
            .map(|q| 1u128 << (2 * q.leading_zeros()))
            .min()
            .unwrap();
        Ok(Self {
            par: par.clone(),
            level,
            min_limit,
            accumulator: Array3::zeros((0, 0, 0)),
            fallback: None,
        })
    }

    /// Compute a dot product, reusing scratch and allocating its result.
    /// Returns an error for empty/unequal inputs, parameter or level
    /// mismatches, or inconsistent ciphertext part counts. Iterators must
    /// yield the same operands when cloned. Timing permission is recomputed
    /// on every call.
    pub fn dot_product_scalar<'a, 'b, I, J>(&mut self, ct: I, pt: J) -> Result<Ciphertext>
    where
        I: Iterator<Item = &'a Ciphertext> + Clone,
        J: Iterator<Item = &'b Plaintext> + Clone,
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
        let ctx = self.par.context_at_level(self.level)?;
        ct_first.validate_for_context(&self.par, self.level, ctx)?;

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
            cti.validate_for_context(&self.par, self.level, ctx)?;
            pti.validate_for_context(&self.par, self.level, ctx)?;
            if cti.len() != ct_first.len() {
                return Err(crate::DotProductError::CiphertextPolynomialCountMismatch {
                    actual: cti.len(),
                    expected: ct_first.len(),
                }
                .into());
            }
        }

        if count as u128 > self.min_limit {
            // Too many ciphertexts for the optimized method, instead, we call
            // polynomial dot products, sharing scratch across ciphertext components.
            let workspace = self
                .fallback
                .get_or_insert_with(|| DotProductWorkspace::new(ctx));
            let c = (0..ct_first.len())
                .map(|i| {
                    workspace
                        .dot_product(
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
            let shape = (ct_first.len(), ctx.moduli().len(), self.par.degree());
            if self.accumulator.dim() != shape {
                // Previous calls always clear scratch before returning or unwinding.
                self.accumulator = Array3::zeros(shape);
            }
            let acc = ClearAccumulator(&mut self.accumulator);
            for (ciphertext, plaintext) in izip!(ct, pt) {
                let pt_coefficients = plaintext.poly_ntt.coefficients();
                for (mut acci, ci) in izip!(acc.0.outer_iter_mut(), ciphertext.iter()) {
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
            for acci in acc.0.outer_iter() {
                c.push(Poly::<Ntt>::try_convert_from(
                    acci,
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
    fn workspace_reuses_fast_scratch_and_switches_paths_and_part_counts()
    -> Result<(), Box<dyn Error>> {
        let params = crate::bfv::BfvParametersBuilder::new()
            .set_degree(16)
            .set_plaintext_modulus(1153)
            .set_moduli_sizes(&[62, 62])
            .build_arc()?;
        let mut rng = rng();
        let sk = SecretKey::random(&params, &mut rng);
        let original = Plaintext::try_encode(&[3u64, 5][..], Encoding::poly(), &params)?;
        let ct: Ciphertext = sk.try_encrypt(&original, &mut rng)?;
        let three_parts = &ct * &ct;
        let mut workspace = super::DotProductScalarWorkspace::new(&params, 0)?;
        assert!(super::DotProductScalarWorkspace::new(&params, 3).is_err());
        for ciphertext in [&ct, &three_parts, &ct] {
            let mut pointer = None;
            for length in [1, 16, 17, 33, 4, 1] {
                for public in [true, false] {
                    let mut pt = original.clone();
                    if public {
                        pt.poly_ntt.allow_variable_time_computations(
                            fhe_traits::VariableTime::new(fhe_traits::PublicData::assert_public()),
                        );
                    }
                    let actual = workspace.dot_product_scalar(
                        std::iter::repeat_n(ciphertext, length),
                        std::iter::repeat_n(&pt, length),
                    )?;
                    let mut expected = Ciphertext::zero(&params);
                    for _ in 0..length {
                        expected += &(ciphertext * &pt);
                    }
                    assert_eq!(actual, expected);
                    assert!(
                        actual
                            .iter()
                            .all(|p| p.allows_variable_time_computations() == public)
                    );
                    assert!(workspace.accumulator.iter().all(|x| *x == 0));
                    if length <= 16 {
                        if let Some(previous) = pointer {
                            assert_eq!(workspace.accumulator.as_ptr(), previous);
                        }
                        pointer = Some(workspace.accumulator.as_ptr());
                    }
                    assert!(
                        workspace
                            .dot_product_scalar(std::iter::empty(), std::iter::once(&pt))
                            .is_err()
                    );
                    assert!(
                        workspace
                            .dot_product_scalar(
                                std::iter::once(ciphertext),
                                std::iter::repeat_n(&pt, 2)
                            )
                            .is_err()
                    );
                }
            }
        }
        assert!(workspace.fallback.is_some());
        let mut lower = ct.clone();
        lower.switch_to_level(1)?;
        assert!(
            workspace
                .dot_product_scalar(std::iter::once(&lower), std::iter::once(&original))
                .is_err()
        );
        let other = BfvParameters::default_arc(1, 16);
        let foreign = Plaintext::try_encode(&[1u64][..], Encoding::poly(), &other)?;
        assert!(
            workspace
                .dot_product_scalar(std::iter::once(&ct), std::iter::once(&foreign))
                .is_err()
        );
        assert_eq!(
            workspace.dot_product_scalar(std::iter::once(&ct), std::iter::once(&original))?,
            &ct * &original
        );
        Ok(())
    }

    #[test]
    fn long_dot_product_reuses_scratch_across_ciphertext_components() -> Result<(), Box<dyn Error>>
    {
        let params = crate::bfv::BfvParametersBuilder::new()
            .set_degree(16)
            .set_plaintext_modulus(1153)
            .set_moduli_sizes(&[62, 62])
            .build_arc()?;
        let mut rng = rng();
        let sk = SecretKey::random(&params, &mut rng);
        let pt = Plaintext::try_encode(&[3u64, 7][..], Encoding::poly(), &params)?;
        let ct: Ciphertext = sk.try_encrypt(&pt, &mut rng)?;
        for length in [17, 33] {
            let ciphertexts = vec![ct.clone(); length];
            let plaintexts = vec![pt.clone(); length];
            let actual = dot_product_scalar(ciphertexts.iter(), plaintexts.iter())?;
            let mut expected = Ciphertext::zero(&params);
            for (ciphertext, plaintext) in ciphertexts.iter().zip(plaintexts.iter()) {
                expected += &(ciphertext * plaintext);
            }
            assert_eq!(actual, expected);
        }
        Ok(())
    }

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
        let ct: Ciphertext = sk.try_encrypt(&pt, &mut rng)?;

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
