//! Lift a coefficient-domain RNS component while reusing its original NTT row.

use super::{Context, Ntt, Poly};
use crate::{Error, Result, ntt::NttOperator, zq::Modulus};
use itertools::izip;
use ndarray::Axis;
use std::sync::Arc;

fn forward_row(row: &mut [u64], modulus: &Modulus, op: &NttOperator, public: bool) {
    if public {
        modulus.lazy_reduce_vec(row);
        unsafe { op.forward_vt_lazy(row.as_mut_ptr()) };
    } else {
        modulus.reduce_vec(row);
        op.forward(row);
    }
}

impl Poly<Ntt> {
    /// Lift one RNS component for multiplication by a
    /// [`Poly<super::NttShoup>`].
    ///
    /// For source modulus q at `index`, each power-basis coefficient is reduced
    /// to its representative in [0, q), then embedded into every target
    /// modulus. This is an unsigned coefficient lift, not a centered basis
    /// extension. The result may have lazy NTT coefficients: it must be
    /// multiplied by an NttShoup polynomial before other operations. Timing
    /// permission is preserved.
    ///
    /// The source NTT row is reused wherever the target contains q. Only other
    /// target rows require a forward transform. Invalid indices, mismatched
    /// degrees, and lazy source coefficients return errors.
    pub fn lift_rns_component_for_shoup(
        &self,
        index: usize,
        target: &Arc<Context>,
    ) -> Result<Poly<Ntt>> {
        let source_modulus = self
            .ctx
            .moduli
            .get(index)
            .ok_or(Error::InvalidRnsComponent {
                index,
                moduli: self.ctx.moduli.len(),
            })?;
        if self.ctx.degree != target.degree {
            return Err(Error::DegreeMismatch {
                found: target.degree,
                expected: self.ctx.degree,
            });
        }
        if self.has_lazy_coefficients {
            return Err(Error::LazyRnsComponentOperand);
        }
        let source = self.coefficients.row(index);
        let mut out = Poly::zero(target);
        out.allow_variable_time_computations = self.allow_variable_time_computations;
        out.has_lazy_coefficients = true;

        // The last output row is temporary power-basis workspace. It is
        // converted to its final value after the other rows have consumed it,
        // so no separate coefficient buffer needs allocating or clearing.
        let last_index = target.q.len() - 1;
        let (mut rows, mut last_row) = out.coefficients.view_mut().split_at(Axis(0), last_index);
        let workspace = last_row.as_slice_mut().unwrap();
        workspace.copy_from_slice(source.as_slice().unwrap());
        if last_index == 0 && target.moduli[0] == *source_modulus {
            return Ok(out);
        }
        let source_op = &self.ctx.ops[index];
        if self.allow_variable_time_computations {
            unsafe { source_op.backward_vt(workspace.as_mut_ptr()) };
        } else {
            source_op.backward(workspace);
        }
        for (mut row, modulus, op) in
            izip!(rows.outer_iter_mut(), target.q.iter(), target.ops.iter(),)
        {
            if **modulus == *source_modulus {
                // Operators are constructed deterministically for (q, degree),
                // including when the contexts were constructed independently.
                row.assign(&source);
            } else {
                let row = row.as_slice_mut().unwrap();
                row.copy_from_slice(workspace);
                forward_row(row, modulus, op, self.allow_variable_time_computations);
            }
        }
        if target.moduli[last_index] == *source_modulus {
            workspace.copy_from_slice(source.as_slice().unwrap());
        } else {
            forward_row(
                workspace,
                &target.q[last_index],
                &target.ops[last_index],
                self.allow_variable_time_computations,
            );
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rq::{PowerBasis, traits::TryConvertFrom};
    use crate::zq::primes::generate_prime;
    use ndarray::Array2;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    #[test]
    fn component_lifts_match_unsigned_coefficient_reference() -> Result<()> {
        let mut rng = ChaCha8Rng::seed_from_u64(0xd1617);
        for degree in [8, 16, 1024] {
            let moduli: Vec<_> = [18, 40, 62, 55]
                .iter()
                .map(|bits| generate_prime(*bits, 2 * degree as u64, 1 << bits).unwrap())
                .collect();
            let source = Context::new_arc(&moduli[..3], degree)?;
            let edges = Array2::from_shape_fn((3, degree), |(row, col)| {
                let q = moduli[row];
                match col % 6 {
                    0 => 0,
                    1 => 1,
                    2 => q / 2,
                    3 => q / 2 + 1,
                    4 => q - 2,
                    _ => q - 1,
                }
            });
            let inputs = [
                Poly::<PowerBasis>::try_convert_from(edges, &source, false)?,
                Poly::<PowerBasis>::random(&source, &mut rng),
            ];
            // Match at the start, middle, or workspace row; include no overlap
            // and single-modulus targets, plus independently built equal contexts.
            for target_moduli in [
                vec![moduli[0], moduli[1], moduli[2]],
                vec![moduli[2], moduli[0], moduli[1]],
                vec![moduli[0]],
                vec![moduli[3]],
                vec![moduli[2], moduli[3], moduli[0]],
                vec![moduli[1], moduli[3]],
            ] {
                let target = Context::new_arc(&target_moduli, degree)?;
                let one = Poly::<PowerBasis>::try_convert_from(&[1u64][..], &target, true)?
                    .into_ntt_shoup();
                for input in &inputs {
                    for public in [false, true] {
                        let mut transformed = input.clone().into_ntt();
                        transformed.allow_variable_time_computations = public;
                        for index in 0..source.moduli.len() {
                            let source_row = input.coefficients.row(index);
                            let expected = Poly::<PowerBasis>::try_convert_from(
                                source_row.as_slice().unwrap(),
                                &target,
                                public,
                            )?
                            .into_ntt();
                            let mut actual =
                                transformed.lift_rns_component_for_shoup(index, &target)?;
                            assert!(actual.has_lazy_coefficients);
                            assert_eq!(actual.allows_variable_time_computations(), public);
                            // Multiplication by one canonicalizes lazy rows using
                            // the same supported operation as key switching.
                            actual *= &one;
                            assert!(!actual.has_lazy_coefficients);
                            assert_eq!(actual, expected);
                            assert_eq!(actual.allows_variable_time_computations(), public);
                        }
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn component_lift_rejects_invalid_inputs() -> Result<()> {
        let ctx = Context::new_arc(&[97, 193], 16)?;
        let input = Poly::<Ntt>::zero(&ctx);
        for index in [2, usize::MAX] {
            assert_eq!(
                input.lift_rns_component_for_shoup(index, &ctx),
                Err(Error::InvalidRnsComponent { index, moduli: 2 })
            );
        }
        let other_degree = Context::new_arc(&[97, 193], 8)?;
        assert_eq!(
            input.lift_rns_component_for_shoup(0, &other_degree),
            Err(Error::DegreeMismatch {
                found: 8,
                expected: 16
            })
        );
        let lazy = Poly::create_constant_ntt_polynomial_with_lazy_coefficients_and_variable_time(
            &[3u64; 16],
            &ctx,
            fhe_traits::VariableTime::new(fhe_traits::PublicData::assert_public()),
        );
        assert_eq!(
            lazy.lift_rns_component_for_shoup(0, &ctx),
            Err(Error::LazyRnsComponentOperand)
        );
        Ok(())
    }
}
