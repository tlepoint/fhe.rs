#![warn(missing_docs, unused_imports)]

//! Polynomial scaler.

use super::{Context, Poly, Representation};
use crate::{
    rns::{RnsScaler, ScalingFactor},
    Error, Result,
};
use itertools::izip;
use ndarray::{s, Array2, Axis};
use std::sync::Arc;

/// Context extender.
#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub struct Scaler {
    from: Arc<Context>,
    to: Arc<Context>,
    number_common_moduli: usize,
    scaler: RnsScaler,
}

impl Scaler {
    /// Create a scaler from a context `from` to a context `to`.
    pub fn new(from: &Arc<Context>, to: &Arc<Context>, factor: ScalingFactor) -> Result<Self> {
        if from.degree != to.degree {
            return Err(Error::Default("Incompatible degrees".to_string()));
        }

        let mut number_common_moduli = 0;
        if factor.is_one {
            for (qi, pi) in izip!(from.q.iter(), to.q.iter()) {
                if qi == pi {
                    number_common_moduli += 1
                } else {
                    break;
                }
            }
        }

        let scaler = RnsScaler::new(&from.rns, &to.rns, factor);

        Ok(Self {
            from: from.clone(),
            to: to.clone(),
            number_common_moduli,
            scaler,
        })
    }

    /// Scale a polynomial
    pub(crate) fn scale(&self, p: &Poly) -> Result<Poly> {
        if p.ctx.as_ref() != self.from.as_ref() {
            Err(Error::Default(
                "The input polynomial does not have the correct context".to_string(),
            ))
        } else {
            let mut representation = p.representation.clone();
            if representation == Representation::NttShoup {
                representation = Representation::Ntt;
            }

            let mut new_coefficients = Array2::<u64>::zeros((self.to.q.len(), self.to.degree));

            if self.number_common_moduli > 0 {
                new_coefficients
                    .slice_mut(s![..self.number_common_moduli, ..])
                    .assign(&p.coefficients.slice(s![..self.number_common_moduli, ..]));
            }

            if self.number_common_moduli < self.to.q.len() {
                if p.representation == Representation::PowerBasis {
                    izip!(
                        new_coefficients
                            .slice_mut(s![self.number_common_moduli.., ..])
                            .axis_iter_mut(Axis(1)),
                        p.coefficients.axis_iter(Axis(1))
                    )
                    .for_each(|(new_column, column)| {
                        self.scaler
                            .scale(column, new_column, self.number_common_moduli)
                    });
                } else if self.number_common_moduli < self.to.q.len() {
                    let mut p_coefficients_powerbasis = p.coefficients.clone();
                    // Backward NTT
                    if p.allow_variable_time_computations {
                        izip!(p_coefficients_powerbasis.outer_iter_mut(), p.ctx.ops.iter())
                            .for_each(|(mut v, op)| unsafe { op.backward_vt(v.as_mut_ptr()) });
                    } else {
                        izip!(p_coefficients_powerbasis.outer_iter_mut(), p.ctx.ops.iter())
                            .for_each(|(mut v, op)| op.backward(v.as_slice_mut().unwrap()));
                    }
                    // Conversion
                    izip!(
                        new_coefficients
                            .slice_mut(s![self.number_common_moduli.., ..])
                            .axis_iter_mut(Axis(1)),
                        p_coefficients_powerbasis.axis_iter(Axis(1))
                    )
                    .for_each(|(new_column, column)| {
                        self.scaler
                            .scale(column, new_column, self.number_common_moduli)
                    });
                    // Forward NTT on the second half
                    if p.allow_variable_time_computations {
                        izip!(
                            new_coefficients
                                .slice_mut(s![self.number_common_moduli.., ..])
                                .outer_iter_mut(),
                            &self.to.ops[self.number_common_moduli..]
                        )
                        .for_each(|(mut v, op)| unsafe { op.forward_vt(v.as_mut_ptr()) });
                    } else {
                        izip!(
                            new_coefficients
                                .slice_mut(s![self.number_common_moduli.., ..])
                                .outer_iter_mut(),
                            &self.to.ops[self.number_common_moduli..]
                        )
                        .for_each(|(mut v, op)| op.forward(v.as_slice_mut().unwrap()));
                    }
                }
            }

            Ok(Poly {
                ctx: self.to.clone(),
                representation,
                allow_variable_time_computations: p.allow_variable_time_computations,
                coefficients: new_coefficients,
                coefficients_shoup: None,
                has_lazy_coefficients: false,
            })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{Scaler, ScalingFactor};
    use crate::rq::{Context, Poly, Representation};
    use itertools::Itertools;
    use num_bigint::BigUint;
    use num_traits::{One, Zero};
    use rand::thread_rng;
    use std::{error::Error, sync::Arc};

    // Moduli to be used in tests.
    static Q: &[u64; 3] = &[
        4611686018282684417,
        4611686018326724609,
        4611686018309947393,
    ];

    static P: &[u64; 3] = &[
        4611686018282684417,
        4611686018309947393,
        4611686018257518593,
    ];

    #[test]
    fn scaler() -> Result<(), Box<dyn Error>> {
        let mut rng = thread_rng();
        let ntests = 100;
        let from = Arc::new(Context::new(Q, 16)?);
        let to = Arc::new(Context::new(P, 16)?);

        for numerator in &[1u64, 2, 3, 100, 1000, 4611686018326724610] {
            for denominator in &[1u64, 2, 3, 4, 100, 101, 1000, 1001, 4611686018326724610] {
                let n = BigUint::from(*numerator);
                let d = BigUint::from(*denominator);

                let scaler = Scaler::new(&from, &to, ScalingFactor::new(&n, &d))?;

                for _ in 0..ntests {
                    let mut poly = Poly::random(&from, Representation::PowerBasis, &mut rng);
                    let poly_biguint = Vec::<BigUint>::from(&poly);

                    let scaled_poly = scaler.scale(&poly)?;
                    let scaled_biguint = Vec::<BigUint>::from(&scaled_poly);

                    let expected = poly_biguint
                        .iter()
                        .map(|i| {
                            if i >= &(from.modulus() >> 1usize) {
                                if &d & BigUint::one() == BigUint::zero() {
                                    to.modulus()
                                        - (&(&(from.modulus() - i) * &n + ((&d >> 1usize) - 1u64))
                                            / &d)
                                            % to.modulus()
                                } else {
                                    to.modulus()
                                        - (&(&(from.modulus() - i) * &n + (&d >> 1)) / &d)
                                            % to.modulus()
                                }
                            } else {
                                ((i * &n + (&d >> 1)) / &d) % to.modulus()
                            }
                        })
                        .collect_vec();
                    assert_eq!(expected, scaled_biguint);

                    poly.change_representation(Representation::Ntt);
                    let mut scaled_poly = scaler.scale(&poly)?;
                    scaled_poly.change_representation(Representation::PowerBasis);
                    let scaled_biguint = Vec::<BigUint>::from(&scaled_poly);
                    assert_eq!(expected, scaled_biguint);
                }
            }
        }

        Ok(())
    }
}
