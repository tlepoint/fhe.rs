//! Tensor products in the extended BFV multiplication basis.

use fhe_math::rq::{Ntt, Poly};
use zeroize::Zeroize;

// Clear restricted intermediates on normal return, errors, and unwinding.
// Public ciphertext buffers need no extra memory pass when discarded.
pub(super) struct Scratch<T: AsMut<[Poly<Ntt>]>>(pub(super) T);

impl<T: AsMut<[Poly<Ntt>]>> Drop for Scratch<T> {
    fn drop(&mut self) {
        for p in self.0.as_mut() {
            if !p.allows_variable_time_computations() {
                p.zeroize();
            }
        }
    }
}

// Callers supply nonempty, canonical NTT polynomials in a common context.
// Restrict every output if any input part is restricted, including parts that
// do not contribute directly to that output coefficient.
fn restrict_timing(output: &mut [Poly<Ntt>], public: bool) {
    if !public {
        for p in output {
            p.disallow_variable_time_computations();
        }
    }
}

pub(super) fn product(lhs: &[Poly<Ntt>], rhs: &[Poly<Ntt>]) -> Vec<Poly<Ntt>> {
    let public = lhs
        .iter()
        .chain(rhs)
        .all(Poly::allows_variable_time_computations);
    let mut output = if let ([a0, a1], [b0, b1]) = (lhs, rhs) {
        // Karatsuba: three pointwise products instead of four. Perform the
        // additions AFTER basis extension: centered lifting is not linear.
        let c0 = a0 * b0;
        let c2 = a1 * b1;
        let mut c1 = a0 + a1;
        let sum = Scratch([b0 + b1]);
        c1 *= &sum.0[0];
        c1 -= &c0;
        c1 -= &c2;
        vec![c0, c1, c2]
    } else {
        let mut output = vec![Poly::zero(lhs[0].ctx()); lhs.len() + rhs.len() - 1];
        if public {
            for p in &mut output {
                p.allow_variable_time_computations(crate::VariableTime::new(
                    crate::PublicData::assert_public(),
                ));
            }
        }
        for (i, a) in lhs.iter().enumerate() {
            for (j, b) in rhs.iter().enumerate() {
                let term = Scratch([a * b]);
                output[i + j] += &term.0[0];
            }
        }
        output
    };
    restrict_timing(&mut output, public);
    output
}

pub(super) fn square(input: &[Poly<Ntt>]) -> Vec<Poly<Ntt>> {
    let public = input.iter().all(Poly::allows_variable_time_computations);
    let mut output = if let [a0, a1] = input {
        let cross = Scratch([a0 * a1]);
        vec![a0 * a0, &cross.0[0] + &cross.0[0], a1 * a1]
    } else {
        let mut output = vec![Poly::zero(input[0].ctx()); 2 * input.len() - 1];
        if public {
            for p in &mut output {
                p.allow_variable_time_computations(crate::VariableTime::new(
                    crate::PublicData::assert_public(),
                ));
            }
        }
        for (i, a) in input.iter().enumerate() {
            for (j, b) in input.iter().enumerate().skip(i) {
                let term = Scratch([a * b]);
                output[i + j] += &term.0[0];
                if i != j {
                    output[i + j] += &term.0[0];
                }
            }
        }
        output
    };
    restrict_timing(&mut output, public);
    output
}

#[cfg(test)]
mod tests {
    use super::*;
    use fhe_math::rq::Context;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    // Independent full convolution, including both off-diagonal terms.
    fn schoolbook(lhs: &[Poly<Ntt>], rhs: &[Poly<Ntt>]) -> Vec<Poly<Ntt>> {
        let mut c = vec![Poly::zero(lhs[0].ctx()); lhs.len() + rhs.len() - 1];
        let public = lhs
            .iter()
            .chain(rhs)
            .all(Poly::allows_variable_time_computations);
        for p in &mut c {
            if public {
                p.allow_variable_time_computations(crate::VariableTime::new(
                    crate::PublicData::assert_public(),
                ));
            }
        }
        for (i, a) in lhs.iter().enumerate() {
            for (j, b) in rhs.iter().enumerate() {
                c[i + j] += &(a * b);
            }
        }
        c
    }

    #[test]
    fn scratch_clears_only_restricted_parts() -> crate::Result<()> {
        let ctx = Context::new_arc(&[97, 193], 16)?;
        let mut rng = ChaCha8Rng::seed_from_u64(0xc1ea2);
        let mut parts = [
            Poly::<Ntt>::random(&ctx, &mut rng),
            Poly::<Ntt>::random(&ctx, &mut rng),
        ];
        parts[1].allow_variable_time_computations(crate::VariableTime::new(
            crate::PublicData::assert_public(),
        ));
        let public = parts[1].clone();
        drop(Scratch(&mut parts));
        assert_eq!(parts[0], Poly::zero(&ctx));
        assert_eq!(parts[1], public);
        Ok(())
    }

    #[test]
    fn optimized_tensors_match_full_convolution() -> crate::Result<()> {
        let mut rng = ChaCha8Rng::seed_from_u64(0x5a0a2e);
        let ctx = Context::new_arc(&[97, 193], 16)?;
        for left_len in 2..=5 {
            for right_len in 2..=5 {
                for restricted in 0..=left_len + right_len {
                    let mut parts: Vec<_> = (0..left_len + right_len)
                        .map(|_| Poly::<Ntt>::random(&ctx, &mut rng))
                        .collect();
                    for (i, p) in parts.iter_mut().enumerate() {
                        if i != restricted {
                            p.allow_variable_time_computations(crate::VariableTime::new(
                                crate::PublicData::assert_public(),
                            ));
                        }
                    }
                    let (left, right) = parts.split_at(left_len);
                    let result = product(left, right);
                    assert!(result.iter().all(|p| p.allows_variable_time_computations()
                        == (restricted == left_len + right_len)));
                    assert!(
                        square(left)
                            .iter()
                            .all(|p| p.allows_variable_time_computations()
                                == (restricted >= left_len))
                    );
                    assert_eq!(result, schoolbook(left, right));
                    assert_eq!(square(left), schoolbook(left, left));
                }
            }
        }
        Ok(())
    }
}
