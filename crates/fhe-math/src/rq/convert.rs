//! Implementation of conversions from and to polynomials.

use super::{traits::TryConvertFrom, Context, Poly, Representation};
use crate::{
    proto::rq::{Representation as RepresentationProto, Rq},
    Error, Result,
};
use itertools::{izip, Itertools};
use ndarray::{Array2, ArrayView, Axis};
use num_bigint::BigUint;
use std::sync::Arc;
use zeroize::{Zeroize, Zeroizing};

impl From<&Poly> for Rq {
    fn from(p: &Poly) -> Self {
        assert!(!p.has_lazy_coefficients);

        let mut q = p.clone();
        if p.representation != Representation::PowerBasis {
            q.change_representation(Representation::PowerBasis);
        }

        let mut proto = Rq::default();
        match p.representation {
            Representation::PowerBasis => {
                proto.representation = RepresentationProto::Powerbasis as i32;
            }
            Representation::Ntt => {
                proto.representation = RepresentationProto::Ntt as i32;
            }
            Representation::NttShoup => {
                proto.representation = RepresentationProto::Nttshoup as i32;
            }
        }
        let mut serialization_length = 0;
        p.ctx
            .q
            .iter()
            .for_each(|qi| serialization_length += qi.serialization_length(p.ctx.degree));
        let mut serialization = Vec::with_capacity(serialization_length);

        izip!(q.coefficients.outer_iter(), p.ctx.q.iter())
            .for_each(|(v, qi)| serialization.append(&mut qi.serialize_vec(v.as_slice().unwrap())));
        proto.coefficients = serialization;
        proto.degree = p.ctx.degree as u32;
        proto.allow_variable_time = p.allow_variable_time_computations;
        proto
    }
}

impl TryConvertFrom<Vec<u64>> for Poly {
    fn try_convert_from<R>(
        mut v: Vec<u64>,
        ctx: &Arc<Context>,
        variable_time: bool,
        representation: R,
    ) -> Result<Self>
    where
        R: Into<Option<Representation>>,
    {
        let repr = representation.into();
        match repr {
            Some(Representation::Ntt) => {
                if let Ok(coefficients) = Array2::from_shape_vec((ctx.q.len(), ctx.degree), v) {
                    Ok(Self {
                        ctx: ctx.clone(),
                        representation: repr.unwrap(),
                        allow_variable_time_computations: variable_time,
                        coefficients,
                        coefficients_shoup: None,
                        has_lazy_coefficients: false,
                    })
                } else {
                    Err(Error::Default(
                        "In Ntt representation, all coefficients must be specified".to_string(),
                    ))
                }
            }
            Some(Representation::NttShoup) => {
                if let Ok(coefficients) = Array2::from_shape_vec((ctx.q.len(), ctx.degree), v) {
                    let mut p = Self {
                        ctx: ctx.clone(),
                        representation: repr.unwrap(),
                        allow_variable_time_computations: variable_time,
                        coefficients,
                        coefficients_shoup: None,
                        has_lazy_coefficients: false,
                    };
                    p.compute_coefficients_shoup();
                    Ok(p)
                } else {
                    Err(Error::Default(
                        "In NttShoup representation, all coefficients must be specified"
                            .to_string(),
                    ))
                }
            }
            Some(Representation::PowerBasis) => {
                if v.len() == ctx.q.len() * ctx.degree {
                    let coefficients =
                        Array2::from_shape_vec((ctx.q.len(), ctx.degree), v).unwrap();
                    Ok(Self {
                        ctx: ctx.clone(),
                        representation: repr.unwrap(),
                        allow_variable_time_computations: variable_time,
                        coefficients,
                        coefficients_shoup: None,
                        has_lazy_coefficients: false,
                    })
                } else if v.len() <= ctx.degree {
                    let mut out = Self::zero(ctx, repr.unwrap());
                    if variable_time {
                        unsafe {
                            izip!(out.coefficients.outer_iter_mut(), ctx.q.iter()).for_each(
                                |(mut w, qi)| {
                                    let wi = w.as_slice_mut().unwrap();
                                    wi[..v.len()].copy_from_slice(&v);
                                    qi.reduce_vec_vt(wi);
                                },
                            );
                            out.allow_variable_time_computations();
                        }
                    } else {
                        izip!(out.coefficients.outer_iter_mut(), ctx.q.iter()).for_each(
                            |(mut w, qi)| {
                                let wi = w.as_slice_mut().unwrap();
                                wi[..v.len()].copy_from_slice(&v);
                                qi.reduce_vec(wi);
                            },
                        );
                        v.zeroize();
                    }
                    Ok(out)
                } else {
                    Err(Error::Default("In PowerBasis representation, either all coefficients must be specified, or only coefficients up to the degree".to_string()))
                }
            }
            None => Err(Error::Default(
                "When converting from a vector, the representation needs to be specified"
                    .to_string(),
            )),
        }
    }
}

impl TryConvertFrom<&Rq> for Poly {
    fn try_convert_from<R>(
        value: &Rq,
        ctx: &Arc<Context>,
        variable_time: bool,
        representation: R,
    ) -> Result<Self>
    where
        R: Into<Option<Representation>>,
    {
        let repr = value
            .representation
            .try_into()
            .map_err(|_| Error::Default("Invalid representation".to_string()))?;
        let representation_from_proto = match repr {
            RepresentationProto::Powerbasis => Representation::PowerBasis,
            RepresentationProto::Ntt => Representation::Ntt,
            RepresentationProto::Nttshoup => Representation::NttShoup,
            _ => return Err(Error::Default("Unknown representation".to_string())),
        };

        let variable_time = variable_time || value.allow_variable_time;

        if let Some(r) = representation.into() as Option<Representation> {
            if r != representation_from_proto {
                return Err(Error::Default("The representation asked for does not match the representation in the serialization".to_string()));
            }
        }

        let degree = value.degree as usize;
        if degree % 8 != 0 || degree < 8 {
            return Err(Error::Default("Invalid degree".to_string()));
        }

        let mut expected_nbytes = 0;
        ctx.q
            .iter()
            .for_each(|qi| expected_nbytes += qi.serialization_length(degree));
        if value.coefficients.len() != expected_nbytes {
            return Err(Error::Default("Invalid coefficients".to_string()));
        }

        let mut power_basis_coefficients = Vec::with_capacity(ctx.q.len() * ctx.degree);
        let mut index = 0;
        for i in 0..ctx.q.len() {
            let qi = &ctx.q[i];
            let size = qi.serialization_length(degree);
            let mut v = qi.deserialize_vec(&value.coefficients[index..index + size]);
            power_basis_coefficients.append(&mut v);
            index += size;
        }

        let mut p = Poly::try_convert_from(
            power_basis_coefficients,
            ctx,
            variable_time,
            Representation::PowerBasis,
        )?;
        p.change_representation(representation_from_proto);
        Ok(p)
    }
}

impl TryConvertFrom<Array2<u64>> for Poly {
    fn try_convert_from<R>(
        a: Array2<u64>,
        ctx: &Arc<Context>,
        variable_time: bool,
        representation: R,
    ) -> Result<Self>
    where
        R: Into<Option<Representation>>,
    {
        if a.shape() != [ctx.q.len(), ctx.degree] {
            Err(Error::Default(
                "The array of coefficient does not have the correct shape".to_string(),
            ))
        } else if let Some(repr) = representation.into() {
            let mut p = Self {
                ctx: ctx.clone(),
                representation: repr,
                allow_variable_time_computations: variable_time,
                coefficients: a,
                coefficients_shoup: None,
                has_lazy_coefficients: false,
            };
            if p.representation == Representation::NttShoup {
                p.compute_coefficients_shoup()
            }
            Ok(p)
        } else {
            Err(Error::Default("When converting from a 2-dimensional array, the representation needs to be specified".to_string()))
        }
    }
}

impl<'a> TryConvertFrom<&'a [u64]> for Poly {
    fn try_convert_from<R>(
        v: &'a [u64],
        ctx: &Arc<Context>,
        variable_time: bool,
        representation: R,
    ) -> Result<Self>
    where
        R: Into<Option<Representation>>,
    {
        Poly::try_convert_from(v.to_vec(), ctx, variable_time, representation)
    }
}

impl<'a> TryConvertFrom<&'a [i64]> for Poly {
    fn try_convert_from<R>(
        v: &'a [i64],
        ctx: &Arc<Context>,
        variable_time: bool,
        representation: R,
    ) -> Result<Self>
    where
        R: Into<Option<Representation>>,
    {
        if representation.into() != Some(Representation::PowerBasis) {
            Err(Error::Default(
                "Converting signed integer require to import in PowerBasis representation"
                    .to_string(),
            ))
        } else if v.len() <= ctx.degree {
            let mut out = Self::zero(ctx, Representation::PowerBasis);
            if variable_time {
                unsafe { out.allow_variable_time_computations() }
            }
            izip!(out.coefficients.outer_iter_mut(), ctx.q.iter()).for_each(|(mut w, qi)| {
                let wi = w.as_slice_mut().unwrap();
                if variable_time {
                    unsafe { wi[..v.len()].copy_from_slice(&qi.reduce_vec_i64_vt(v)) }
                } else {
                    wi[..v.len()].copy_from_slice(Zeroizing::new(qi.reduce_vec_i64(v)).as_ref());
                }
            });
            Ok(out)
        } else {
            Err(Error::Default("In PowerBasis representation with signed integers, only `degree` coefficients can be specified".to_string()))
        }
    }
}

impl<'a> TryConvertFrom<&'a Vec<i64>> for Poly {
    fn try_convert_from<R>(
        v: &'a Vec<i64>,
        ctx: &Arc<Context>,
        variable_time: bool,
        representation: R,
    ) -> Result<Self>
    where
        R: Into<Option<Representation>>,
    {
        Poly::try_convert_from(v.as_ref() as &[i64], ctx, variable_time, representation)
    }
}

impl<'a> TryConvertFrom<&'a [BigUint]> for Poly {
    fn try_convert_from<R>(
        v: &'a [BigUint],
        ctx: &Arc<Context>,
        variable_time: bool,
        representation: R,
    ) -> Result<Self>
    where
        R: Into<Option<Representation>>,
    {
        let repr = representation.into();

        if v.len() > ctx.degree {
            Err(Error::Default(
                "The slice contains too many big integers compared to the polynomial degree"
                    .to_string(),
            ))
        } else if repr.is_some() {
            let mut coefficients = Array2::zeros((ctx.q.len(), ctx.degree));

            izip!(coefficients.axis_iter_mut(Axis(1)), v).for_each(|(mut c, vi)| {
                c.assign(&ArrayView::from(&ctx.rns.project(vi)));
            });

            let mut p = Self {
                ctx: ctx.clone(),
                representation: repr.unwrap(),
                allow_variable_time_computations: variable_time,
                coefficients,
                coefficients_shoup: None,
                has_lazy_coefficients: false,
            };

            match p.representation {
                Representation::PowerBasis => Ok(p),
                Representation::Ntt => Ok(p),
                Representation::NttShoup => {
                    p.compute_coefficients_shoup();
                    Ok(p)
                }
            }
        } else {
            Err(Error::Default(
                "When converting from a vector, the representation needs to be specified"
                    .to_string(),
            ))
        }
    }
}

impl<'a> TryConvertFrom<&'a Vec<u64>> for Poly {
    fn try_convert_from<R>(
        v: &'a Vec<u64>,
        ctx: &Arc<Context>,
        variable_time: bool,
        representation: R,
    ) -> Result<Self>
    where
        R: Into<Option<Representation>>,
    {
        Poly::try_convert_from(v.to_vec(), ctx, variable_time, representation)
    }
}

impl<'a, const N: usize> TryConvertFrom<&'a [u64; N]> for Poly {
    fn try_convert_from<R>(
        v: &'a [u64; N],
        ctx: &Arc<Context>,
        variable_time: bool,
        representation: R,
    ) -> Result<Self>
    where
        R: Into<Option<Representation>>,
    {
        Poly::try_convert_from(v.as_ref(), ctx, variable_time, representation)
    }
}

impl<'a, const N: usize> TryConvertFrom<&'a [BigUint; N]> for Poly {
    fn try_convert_from<R>(
        v: &'a [BigUint; N],
        ctx: &Arc<Context>,
        variable_time: bool,
        representation: R,
    ) -> Result<Self>
    where
        R: Into<Option<Representation>>,
    {
        Poly::try_convert_from(v.as_ref(), ctx, variable_time, representation)
    }
}

impl<'a, const N: usize> TryConvertFrom<&'a [i64; N]> for Poly {
    fn try_convert_from<R>(
        v: &'a [i64; N],
        ctx: &Arc<Context>,
        variable_time: bool,
        representation: R,
    ) -> Result<Self>
    where
        R: Into<Option<Representation>>,
    {
        Poly::try_convert_from(v.as_ref(), ctx, variable_time, representation)
    }
}

impl From<&Poly> for Vec<u64> {
    fn from(p: &Poly) -> Self {
        p.coefficients.as_slice().unwrap().to_vec()
    }
}

impl From<&Poly> for Vec<BigUint> {
    fn from(p: &Poly) -> Self {
        izip!(p.coefficients.axis_iter(Axis(1)))
            .map(|c| p.ctx.rns.lift(c))
            .collect_vec()
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        proto::rq::Rq,
        rq::{traits::TryConvertFrom, Context, Poly, Representation},
        Error as CrateError,
    };
    use num_bigint::BigUint;
    use rand::thread_rng;
    use std::{error::Error, sync::Arc};

    static MODULI: &[u64; 3] = &[1153, 4611686018326724609, 4611686018309947393];

    #[test]
    fn proto() -> Result<(), Box<dyn Error>> {
        let mut rng = thread_rng();
        for modulus in MODULI {
            let ctx = Arc::new(Context::new(&[*modulus], 16)?);
            let p = Poly::random(&ctx, Representation::PowerBasis, &mut rng);
            let proto = Rq::from(&p);
            assert_eq!(Poly::try_convert_from(&proto, &ctx, false, None)?, p);
            assert_eq!(
                Poly::try_convert_from(&proto, &ctx, false, Representation::PowerBasis)?,
                p
            );
            assert_eq!(
				Poly::try_convert_from(&proto, &ctx, false, Representation::Ntt)
					.expect_err("Should fail because of mismatched representations"),
					CrateError::Default("The representation asked for does not match the representation in the serialization".to_string())
			);
            assert_eq!(
				Poly::try_convert_from(&proto, &ctx, false, Representation::NttShoup)
					.expect_err("Should fail because of mismatched representations"),
					CrateError::Default("The representation asked for does not match the representation in the serialization".to_string())
			);
        }

        let ctx = Arc::new(Context::new(MODULI, 16)?);
        let p = Poly::random(&ctx, Representation::PowerBasis, &mut rng);
        let proto = Rq::from(&p);
        assert_eq!(Poly::try_convert_from(&proto, &ctx, false, None)?, p);
        assert_eq!(
            Poly::try_convert_from(&proto, &ctx, false, Representation::PowerBasis)?,
            p
        );
        assert_eq!(
			Poly::try_convert_from(&proto, &ctx, false, Representation::Ntt)
				.expect_err("Should fail because of mismatched representations"),
				CrateError::Default("The representation asked for does not match the representation in the serialization".to_string())
		);
        assert_eq!(
			Poly::try_convert_from(&proto, &ctx, false, Representation::NttShoup)
				.expect_err("Should fail because of mismatched representations"),
				CrateError::Default("The representation asked for does not match the representation in the serialization".to_string())
		);

        let ctx = Arc::new(Context::new(&MODULI[0..1], 16)?);
        assert_eq!(
            Poly::try_convert_from(&proto, &ctx, false, None)
                .expect_err("Should fail because of incorrect context"),
            CrateError::Default("Invalid coefficients".to_string())
        );

        Ok(())
    }

    #[test]
    fn try_convert_from_slice_zero() -> Result<(), Box<dyn Error>> {
        for modulus in MODULI {
            let ctx = Arc::new(Context::new(&[*modulus], 16)?);

            // Power Basis
            assert_eq!(
                Poly::try_convert_from(&[0u64], &ctx, false, Representation::PowerBasis)?,
                Poly::zero(&ctx, Representation::PowerBasis)
            );
            assert_eq!(
                Poly::try_convert_from(&[0i64], &ctx, false, Representation::PowerBasis)?,
                Poly::zero(&ctx, Representation::PowerBasis)
            );
            assert_eq!(
                Poly::try_convert_from(&[0u64; 16], &ctx, false, Representation::PowerBasis)?,
                Poly::zero(&ctx, Representation::PowerBasis)
            );
            assert_eq!(
                Poly::try_convert_from(&[0i64; 16], &ctx, false, Representation::PowerBasis)?,
                Poly::zero(&ctx, Representation::PowerBasis)
            );
            assert!(Poly::try_convert_from(
                &[0u64; 17], // One too many
                &ctx,
                false,
                Representation::PowerBasis,
            )
            .is_err());

            // Ntt
            assert!(Poly::try_convert_from(&[0u64], &ctx, false, Representation::Ntt).is_err());
            assert!(Poly::try_convert_from(&[0i64], &ctx, false, Representation::Ntt).is_err());
            assert_eq!(
                Poly::try_convert_from(&[0u64; 16], &ctx, false, Representation::Ntt)?,
                Poly::zero(&ctx, Representation::Ntt)
            );
            assert!(Poly::try_convert_from(&[0i64; 16], &ctx, false, Representation::Ntt).is_err());
            assert!(Poly::try_convert_from(
                &[0u64; 17], // One too many
                &ctx,
                false,
                Representation::Ntt,
            )
            .is_err());
        }

        let ctx = Arc::new(Context::new(MODULI, 16)?);
        assert_eq!(
            Poly::try_convert_from(
                Vec::<u64>::default(),
                &ctx,
                false,
                Representation::PowerBasis,
            )?,
            Poly::zero(&ctx, Representation::PowerBasis)
        );
        assert!(
            Poly::try_convert_from(Vec::<u64>::default(), &ctx, false, Representation::Ntt)
                .is_err()
        );

        assert_eq!(
            Poly::try_convert_from(&[0u64], &ctx, false, Representation::PowerBasis)?,
            Poly::zero(&ctx, Representation::PowerBasis)
        );
        assert!(Poly::try_convert_from(&[0u64], &ctx, false, Representation::Ntt).is_err());

        assert_eq!(
            Poly::try_convert_from(&[0u64; 16], &ctx, false, Representation::PowerBasis)?,
            Poly::zero(&ctx, Representation::PowerBasis)
        );
        assert!(Poly::try_convert_from(&[0u64; 16], &ctx, false, Representation::Ntt).is_err());

        assert!(
            Poly::try_convert_from(&[0u64; 17], &ctx, false, Representation::PowerBasis).is_err()
        );
        assert!(Poly::try_convert_from(&[0u64; 17], &ctx, false, Representation::Ntt).is_err());

        assert_eq!(
            Poly::try_convert_from(&[0u64; 16], &ctx, false, Representation::PowerBasis)?,
            Poly::zero(&ctx, Representation::PowerBasis)
        );
        assert_eq!(
            Poly::try_convert_from(&[0u64; 48], &ctx, false, Representation::Ntt)?,
            Poly::zero(&ctx, Representation::Ntt)
        );

        Ok(())
    }

    #[test]
    fn try_convert_from_vec_zero() -> Result<(), Box<dyn Error>> {
        for modulus in MODULI {
            let ctx = Arc::new(Context::new(&[*modulus], 16)?);
            assert_eq!(
                Poly::try_convert_from(vec![], &ctx, false, Representation::PowerBasis)?,
                Poly::zero(&ctx, Representation::PowerBasis)
            );
            assert!(Poly::try_convert_from(vec![], &ctx, false, Representation::Ntt).is_err());

            assert_eq!(
                Poly::try_convert_from(vec![0], &ctx, false, Representation::PowerBasis)?,
                Poly::zero(&ctx, Representation::PowerBasis)
            );
            assert!(Poly::try_convert_from(vec![0], &ctx, false, Representation::Ntt).is_err());

            assert_eq!(
                Poly::try_convert_from(vec![0; 16], &ctx, false, Representation::PowerBasis)?,
                Poly::zero(&ctx, Representation::PowerBasis)
            );
            assert_eq!(
                Poly::try_convert_from(vec![0; 16], &ctx, false, Representation::Ntt)?,
                Poly::zero(&ctx, Representation::Ntt)
            );

            assert!(
                Poly::try_convert_from(vec![0; 17], &ctx, false, Representation::PowerBasis)
                    .is_err()
            );
            assert!(Poly::try_convert_from(vec![0; 17], &ctx, false, Representation::Ntt).is_err());
        }

        let ctx = Arc::new(Context::new(MODULI, 16)?);
        assert_eq!(
            Poly::try_convert_from(vec![], &ctx, false, Representation::PowerBasis)?,
            Poly::zero(&ctx, Representation::PowerBasis)
        );
        assert!(Poly::try_convert_from(vec![], &ctx, false, Representation::Ntt).is_err());

        assert_eq!(
            Poly::try_convert_from(vec![0], &ctx, false, Representation::PowerBasis)?,
            Poly::zero(&ctx, Representation::PowerBasis)
        );
        assert!(Poly::try_convert_from(vec![0], &ctx, false, Representation::Ntt).is_err());

        assert_eq!(
            Poly::try_convert_from(vec![0; 16], &ctx, false, Representation::PowerBasis)?,
            Poly::zero(&ctx, Representation::PowerBasis)
        );
        assert!(Poly::try_convert_from(vec![0; 16], &ctx, false, Representation::Ntt).is_err());

        assert!(
            Poly::try_convert_from(vec![0; 17], &ctx, false, Representation::PowerBasis).is_err()
        );
        assert!(Poly::try_convert_from(vec![0; 17], &ctx, false, Representation::Ntt).is_err());

        assert_eq!(
            Poly::try_convert_from(vec![0; 48], &ctx, false, Representation::PowerBasis)?,
            Poly::zero(&ctx, Representation::PowerBasis)
        );
        assert_eq!(
            Poly::try_convert_from(vec![0; 48], &ctx, false, Representation::Ntt)?,
            Poly::zero(&ctx, Representation::Ntt)
        );

        Ok(())
    }

    #[test]
    fn biguint() -> Result<(), Box<dyn Error>> {
        let mut rng = thread_rng();
        for _ in 0..100 {
            for modulus in MODULI {
                let ctx = Arc::new(Context::new(&[*modulus], 16)?);
                let p = Poly::random(&ctx, Representation::PowerBasis, &mut rng);
                let p_coeffs = Vec::<BigUint>::from(&p);
                let q = Poly::try_convert_from(
                    p_coeffs.as_slice(),
                    &ctx,
                    false,
                    Representation::PowerBasis,
                )?;
                assert_eq!(p, q);
            }

            let ctx = Arc::new(Context::new(MODULI, 16)?);
            let p = Poly::random(&ctx, Representation::PowerBasis, &mut rng);
            let p_coeffs = Vec::<BigUint>::from(&p);
            assert_eq!(p_coeffs.len(), ctx.degree);
            let q = Poly::try_convert_from(
                p_coeffs.as_slice(),
                &ctx,
                false,
                Representation::PowerBasis,
            )?;
            assert_eq!(p, q);
        }
        Ok(())
    }
}
