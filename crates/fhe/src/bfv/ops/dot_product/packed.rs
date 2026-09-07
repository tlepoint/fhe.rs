//! Fused unpacking and accumulation of compact plaintext NTT residues.

use super::{ClearAccumulator, DotProductScalarWorkspace};
use crate::bfv::Parameters;
use crate::{
    Error, Result,
    bfv::{Ciphertext, PackedPlaintextView},
};
use fhe_math::rq::{Ntt, Poly, traits::TryConvertFrom};
use ndarray::Array3;

// Specialize the eight possible bit offsets, rather than particular moduli.
// Each group of eight residues starts at a byte boundary. Narrow residues
// fit in one unaligned word; wide residues can need a ninth byte. All shifts
// and branches depend on the public modulus, including for restricted inputs.
unsafe fn fma_packed<const REM: usize, const WIDE: bool>(
    out: &mut [u128],
    x: &[u64],
    packed: &[u8],
    bits: usize,
) {
    let mask = u64::MAX >> (64 - bits);
    let bytes = bits / 8;
    macro_rules! at {
        ($i:expr, $j:expr, $base:expr) => {
            unsafe {
                let word = u64::from_le(std::ptr::read_unaligned(
                    packed
                        .as_ptr()
                        .add($base + $j * bytes + $j * REM / 8)
                        .cast::<u64>(),
                ));
                let mut y = word >> ($j * REM % 8);
                if WIDE && bits + ($j * REM % 8) > 64 {
                    y |= u64::from(*packed.get_unchecked($base + $j * bytes + $j * REM / 8 + 8))
                        << (64 - ($j * REM % 8));
                }
                let y = y & mask;
                *out.get_unchecked_mut($i + $j) +=
                    u128::from(*x.get_unchecked($i + $j)) * u128::from(y);
            }
        };
    }
    for i in (0..out.len()).step_by(8) {
        let base = i / 8 * bits;
        at!(i, 0, base);
        at!(i, 1, base);
        at!(i, 2, base);
        at!(i, 3, base);
        at!(i, 4, base);
        at!(i, 5, base);
        at!(i, 6, base);
        at!(i, 7, base);
    }
}

// Check memory bounds even if a caller's cloned iterator changes its operands
// after validation. The workspace separately bounds unreduced sums to u128.
fn fma_packed_row(out: &mut [u128], x: &[u64], packed: &[u8], bits: usize) {
    assert_eq!(x.len(), out.len());
    assert!(out.len().is_multiple_of(8));
    assert!((1..=64).contains(&bits));
    assert!(packed.len() >= out.len() / 8 * bits + 8);
    macro_rules! run {
        ($wide:expr) => {
            unsafe {
                match bits % 8 {
                    0 => fma_packed::<0, $wide>(out, x, packed, bits),
                    1 => fma_packed::<1, $wide>(out, x, packed, bits),
                    2 => fma_packed::<2, $wide>(out, x, packed, bits),
                    3 => fma_packed::<3, $wide>(out, x, packed, bits),
                    4 => fma_packed::<4, $wide>(out, x, packed, bits),
                    5 => fma_packed::<5, $wide>(out, x, packed, bits),
                    6 => fma_packed::<6, $wide>(out, x, packed, bits),
                    7 => fma_packed::<7, $wide>(out, x, packed, bits),
                    _ => unreachable!(),
                }
            }
        };
    }
    if bits <= 57 {
        run!(false);
    } else {
        run!(true);
    }
}

impl DotProductScalarWorkspace {
    /// Compute a dot product directly from bit-packed NTT residues.
    ///
    /// No transforms or unpacked plaintext buffers are needed. Returns the same
    /// result and timing permission as [`Self::dot_product_scalar`] on unpacked
    /// inputs. Errors cover empty or unequal inputs, foreign parameters, wrong
    /// levels, and inconsistent ciphertext part counts. Cloned iterators must
    /// yield the same operands. Long products periodically reduce the shared
    /// accumulator to preserve its overflow bound.
    pub fn dot_product_scalar_packed<'a, 'b, I, J>(&mut self, ct: I, pt: J) -> Result<Ciphertext>
    where
        I: Iterator<Item = &'a Ciphertext> + Clone,
        J: Iterator + Clone,
        J::Item: Into<PackedPlaintextView<'b>>,
    {
        let pt = pt.map(Into::into);
        let count = ct.clone().count();
        let pt_count = pt.clone().count();
        if count == 0 || pt_count == 0 {
            return Err(crate::DotProductError::EmptyInput.into());
        }
        if count != pt_count {
            return Err(crate::DotProductError::OperandCountMismatch {
                ciphertexts: count,
                plaintexts: pt_count,
            }
            .into());
        }
        let first = ct.clone().next().unwrap();
        let ctx = self.par.context_at_level(self.level)?;
        first.validate_for_context(&self.par, self.level, ctx)?;
        let mut public = true;
        for (ct, pt) in ct.clone().zip(pt.clone()) {
            ct.validate_for_context(&self.par, self.level, ctx)?;
            if !Parameters::compatible(&self.par, pt.par) {
                return Err(Error::ParameterMismatch {
                    left: crate::ParameterSource::Plaintext,
                    right: crate::ParameterSource::Parameters,
                });
            }
            if pt.level != self.level {
                return Err(Error::InvalidLevel {
                    level: pt.level,
                    min_level: self.level,
                    max_level: self.level,
                });
            }
            if ct.len() != first.len() {
                return Err(crate::DotProductError::CiphertextPolynomialCountMismatch {
                    actual: ct.len(),
                    expected: first.len(),
                }
                .into());
            }
            public &= pt.public && ct.iter().all(Poly::allows_variable_time_computations);
        }
        let shape = (first.len(), ctx.moduli().len(), self.par.degree());
        if self.accumulator.dim() != shape {
            self.accumulator = Array3::zeros(shape);
        }
        let acc = ClearAccumulator(&mut self.accumulator);
        let needs_reduction = count as u128 > self.min_limit;
        // From zero, min_limit canonical products fit in u128. After a
        // reduction reserve one product's space for the previous residue.
        // The bound and schedule depend only on public moduli and input count.
        let interval = if needs_reduction {
            self.min_limit as usize - 1
        } else {
            usize::MAX
        };
        for (index, (ct, pt)) in ct.zip(pt).enumerate() {
            for (mut output, ci) in acc.0.outer_iter_mut().zip(ct.iter()) {
                let coefficients = ci.coefficients();
                let mut offset = 0;
                for ((mut out, x), q) in output
                    .outer_iter_mut()
                    .zip(coefficients.outer_iter())
                    .zip(ctx.moduli())
                {
                    let bits = (64 - q.leading_zeros()) as usize;
                    let packed = pt.coefficients.get(offset..).unwrap();
                    // Validated contexts give identical row lengths. Packing
                    // preserves canonical residues and includes eight pad bytes.
                    fma_packed_row(
                        out.as_slice_mut().unwrap(),
                        x.as_slice().unwrap(),
                        packed,
                        bits,
                    );
                    offset += bits * self.par.degree() / 8;
                }
            }
            if needs_reduction && (index + 1).is_multiple_of(interval) && index + 1 < count {
                for mut part in acc.0.outer_iter_mut() {
                    for (mut row, modulus) in part.outer_iter_mut().zip(ctx.moduli_operators()) {
                        for value in &mut row {
                            *value = if public {
                                unsafe { modulus.reduce_u128_vt(*value) }
                            } else {
                                modulus.reduce_u128(*value)
                            } as u128;
                        }
                    }
                }
            }
        }
        let c = acc
            .0
            .outer_iter()
            .map(|a| {
                Poly::<Ntt>::try_convert_from_with_timing(
                    a,
                    ctx,
                    public.then(|| crate::VariableTime::new(crate::PublicData::assert_public())),
                )
            })
            .collect::<std::result::Result<Vec<_>, _>>()?;
        Ok(Ciphertext {
            par: first.par.clone(),
            seed: None,
            c,
            level: first.level,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bfv::{Encoding, PackedPlaintext, ParametersBuilder, Plaintext};
    use crate::{PublicData, VariableTime};
    use rand::{RngExt, SeedableRng};
    use rand_chacha::ChaCha8Rng;

    #[test]
    fn fused_decoder_covers_every_width_and_final_word() {
        let mut rng = ChaCha8Rng::seed_from_u64(47);
        for bits in 1..=64 {
            let mask = u64::MAX >> (64 - bits);
            for degree in [8, 16, 128] {
                let x: Vec<u64> = (0..degree).map(|_| rng.random()).collect();
                for y in [
                    vec![0; degree],
                    vec![mask; degree],
                    (0..degree).map(|_| rng.random::<u64>() & mask).collect(),
                ] {
                    let mut bytes = fhe_util::transcode_to_bytes(&y, bits);
                    bytes.resize(bytes.len() + 8, 0);
                    let mut out = vec![7; degree];
                    fma_packed_row(&mut out, &x, &bytes, bits);
                    for ((actual, x), y) in out.iter().zip(&x).zip(&y) {
                        assert_eq!(*actual, 7 + u128::from(*x) * u128::from(*y));
                    }
                }
            }
        }
    }

    #[test]
    fn packed_dot_products_preserve_exact_results_bounds_and_permissions() -> Result<()> {
        let public = VariableTime::new(PublicData::assert_public());
        for bits in [
            13, 14, 15, 16, 17, 18, 19, 20, 36, 50, 55, 57, 58, 59, 60, 61, 62,
        ] {
            let par = ParametersBuilder::new()
                .degree(16)
                .plaintext_modulus(17_u64)
                .ciphertext_modulus_bits([bits, bits])
                .build()?;
            for level in [0, 1] {
                let ctx = par.context_at_level(level)?;
                // Worst canonical residues exercise carries and the exact
                // u128 accumulation limit, independently of encryption noise.
                let values: Vec<_> = ctx
                    .moduli()
                    .iter()
                    .flat_map(|q| vec![q - 1; par.degree()])
                    .collect();
                let worst = Poly::<Ntt>::try_convert_from_public(
                    values,
                    ctx,
                    crate::VariableTime::new(crate::PublicData::assert_public()),
                )?;
                let mut pt = Plaintext {
                    par: par.clone(),

                    poly_ntt: worst.clone(),
                };
                let mut workspace = DotProductScalarWorkspace::new(&par, level)?;
                for parts in [2, 3, 2] {
                    let mut ct = Ciphertext {
                        par: par.clone(),
                        seed: None,
                        c: vec![worst.clone(); parts],
                        level,
                    };
                    let mut pointer = None;
                    for length in [1, 15, 16, 17, 31, 33, 4] {
                        for restricted in [false, true, false] {
                            pt.poly_ntt.allow_variable_time_computations(public);
                            if restricted {
                                pt.poly_ntt.disallow_variable_time_computations();
                            }
                            let packed = PackedPlaintext::from(&pt);
                            let restored = packed.unpack();
                            assert_eq!(restored, pt);
                            assert_eq!(
                                restored.poly_ntt.allows_variable_time_computations(),
                                !restricted
                            );
                            assert_eq!(
                                packed.size_bytes(),
                                par.degree() * ctx.moduli().len() * bits / 8 + 8
                            );
                            let actual = workspace.dot_product_scalar_packed(
                                std::iter::repeat_n(&ct, length),
                                std::iter::repeat_n(&packed, length),
                            )?;
                            let mut expected = ct.multiply_plaintext(&pt).unwrap();
                            for _ in 1..length {
                                (expected)
                                    .add_assign(&(ct.multiply_plaintext(&pt).unwrap()))
                                    .unwrap();
                            }
                            assert_eq!(actual, expected);
                            assert!(
                                actual
                                    .iter()
                                    .all(|p| p.allows_variable_time_computations() != restricted)
                            );
                            assert!(workspace.accumulator.iter().all(|x| *x == 0));
                            if let Some(previous) = pointer {
                                assert_eq!(workspace.accumulator.as_ptr(), previous);
                            }
                            pointer = Some(workspace.accumulator.as_ptr());
                            assert_eq!(
                                workspace.dot_product_scalar(
                                    std::iter::repeat_n(&ct, length),
                                    std::iter::repeat_n(&pt, length)
                                )?,
                                actual
                            );
                        }
                    }
                    ct.c.last_mut()
                        .unwrap()
                        .disallow_variable_time_computations();
                    let packed = PackedPlaintext::from(&pt);
                    assert!(
                        workspace
                            .dot_product_scalar_packed(
                                std::iter::once(&ct),
                                std::iter::once(&packed)
                            )?
                            .iter()
                            .all(|p| !p.allows_variable_time_computations())
                    );
                }
            }
        }
        Ok(())
    }

    #[test]
    fn validation_leaves_workspace_reusable() -> Result<()> {
        let par = ParametersBuilder::new()
            .degree(16)
            .plaintext_modulus(17_u64)
            .ciphertext_modulus_bits([40, 40])
            .build()?;
        let pt = Plaintext::encode(&par, &[3u64][..], Encoding::Polynomial)?;
        let packed = PackedPlaintext::from(&pt);
        let ct = Ciphertext {
            par: par.clone(),
            seed: None,
            c: vec![pt.poly_ntt.clone(); 2],
            level: 0,
        };
        let mut workspace = DotProductScalarWorkspace::new(&par, 0)?;
        let compute =
            |workspace: &mut DotProductScalarWorkspace, ct: &Ciphertext, pt: &PackedPlaintext| {
                workspace.dot_product_scalar_packed(std::iter::once(ct), std::iter::once(pt))
            };
        let expected = compute(&mut workspace, &ct, &packed)?;
        assert!(
            workspace
                .dot_product_scalar_packed(std::iter::empty(), std::iter::once(&packed))
                .is_err()
        );
        assert!(
            workspace
                .dot_product_scalar_packed(std::iter::once(&ct), std::iter::repeat_n(&packed, 2))
                .is_err()
        );
        assert!(compute(&mut workspace, &Ciphertext::invalid_empty(&par), &packed).is_err());
        let lower = Plaintext::encode_at_level(&par, &[3u64][..], Encoding::Polynomial, 1)?;
        assert!(compute(&mut workspace, &ct, &PackedPlaintext::from(&lower)).is_err());
        let foreign = ParametersBuilder::new()
            .noise_variance(11)
            .degree(16)
            .plaintext_modulus(17_u64)
            .ciphertext_modulus_bits([40, 40])
            .build()?;
        let foreign = Plaintext::encode(&foreign, &[3u64][..], Encoding::Polynomial)?;
        assert!(compute(&mut workspace, &ct, &PackedPlaintext::from(&foreign)).is_err());
        let mut three = ct.clone();
        three.c.push(pt.poly_ntt.clone());
        assert!(
            workspace
                .dot_product_scalar_packed(
                    [&ct, &three].into_iter(),
                    std::iter::repeat_n(&packed, 2)
                )
                .is_err()
        );
        assert_eq!(compute(&mut workspace, &ct, &packed)?, expected);
        assert!(workspace.accumulator.iter().all(|x| *x == 0));
        let mut cleared = packed.clone();
        zeroize::Zeroize::zeroize(&mut cleared);
        assert!(cleared.coefficients.iter().all(|x| *x == 0));
        assert!(
            cleared
                .unpack()
                .poly_ntt
                .coefficients()
                .iter()
                .all(|x| *x == 0)
        );
        Ok(())
    }

    #[test]
    #[expect(
        clippy::panic,
        reason = "exercise scratch cleanup during iterator unwinding"
    )]
    fn scratch_is_cleared_when_the_input_iterator_panics() -> Result<()> {
        use std::cell::Cell;
        let par = ParametersBuilder::new()
            .degree(16)
            .plaintext_modulus(17_u64)
            .ciphertext_modulus_bits([40])
            .build()?;
        let pt = Plaintext::encode(&par, &[3u64][..], Encoding::Polynomial)?;
        let packed = PackedPlaintext::from(&pt);
        let ct = Ciphertext {
            par: par.clone(),
            seed: None,
            c: vec![pt.poly_ntt.clone(); 2],
            level: 0,
        };
        let mut workspace = DotProductScalarWorkspace::new(&par, 0)?;
        let visits = Cell::new(0);
        let ciphertexts = [&ct, &ct].into_iter().inspect(|_| {
            let visit = visits.get() + 1;
            visits.set(visit);
            // Counting: 2 visits; first operand: 1; validation: 2;
            // arithmetic: allow one accumulation, then fail on the second.
            if visit == 7 {
                panic!("iterator failed during accumulation");
            }
        });
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            workspace.dot_product_scalar_packed(ciphertexts, std::iter::repeat_n(&packed, 2))
        }));
        assert!(result.is_err());
        assert_eq!(visits.get(), 7);
        assert!(!workspace.accumulator.is_empty());
        assert!(workspace.accumulator.iter().all(|x| *x == 0));
        assert_eq!(
            workspace.dot_product_scalar_packed(std::iter::once(&ct), std::iter::once(&packed))?,
            ct.multiply_plaintext(&pt).unwrap()
        );
        Ok(())
    }
}
