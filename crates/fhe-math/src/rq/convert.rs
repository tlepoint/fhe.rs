//! Implementation of conversions from and to polynomials.

use super::{
    Context, Ntt, NttShoup, Poly, PowerBasis, Representation, RepresentationTag,
    traits::TryConvertFrom,
};
use crate::{
    Error, Result,
    proto::rq::{Representation as RepresentationProto, Rq},
};
use itertools::{Itertools, izip};
use ndarray::{Array2, ArrayView, Axis};
use num_bigint::BigUint;
use std::sync::Arc;
use zeroize::{Zeroize, Zeroizing};

impl<R: RepresentationTag> From<&Poly<R>> for Rq {
    fn from(p: &Poly<R>) -> Self {
        assert!(!p.has_lazy_coefficients);
        let q: Poly<PowerBasis> = match R::REPRESENTATION {
            Representation::PowerBasis => Poly::<PowerBasis>::from_parts(p.clone()),
            Representation::Ntt => Poly::<Ntt>::from_parts(p.clone()).into_power_basis(),
            Representation::NttShoup => Poly::<NttShoup>::from_parts(p.clone()).into_power_basis(),
        };

        let mut proto = Rq::default();
        match R::REPRESENTATION {
            Representation::PowerBasis => {
                proto.representation = RepresentationProto::Powerbasis as i32
            }
            Representation::Ntt => proto.representation = RepresentationProto::Ntt as i32,
            Representation::NttShoup => proto.representation = RepresentationProto::Nttshoup as i32,
        }
        let serialization: Vec<u8> = izip!(q.coefficients.outer_iter(), p.ctx.q.iter())
            .flat_map(|(v, qi)| qi.serialize_vec(v.as_slice().unwrap()))
            .collect();
        proto.coefficients = serialization;
        proto.degree = p.ctx.degree as u32;
        proto.allow_variable_time = p.allow_variable_time_computations;
        proto
    }
}

fn parse_proto(
    value: &Rq,
    ctx: &Arc<Context>,
    variable_time: bool,
) -> Result<(Representation, Vec<u64>, bool)> {
    let repr = value
        .representation
        .try_into()
        .map_err(|_| Error::Default("Invalid representation".to_string()))?;
    let representation_from_proto = match repr {
        RepresentationProto::Powerbasis => Representation::PowerBasis,
        RepresentationProto::Ntt => Representation::Ntt,
        RepresentationProto::Nttshoup => Representation::NttShoup,
        RepresentationProto::Unknown => {
            return Err(Error::Default("Unknown representation".to_string()));
        }
    };

    let variable_time = variable_time || value.allow_variable_time;

    let degree = value.degree as usize;
    if !degree.is_multiple_of(8) || degree < 8 {
        return Err(Error::Default("Invalid degree".to_string()));
    }

    let mut expected_nbytes = 0;
    ctx.q
        .iter()
        .for_each(|qi| expected_nbytes += qi.serialization_length(degree));
    if value.coefficients.len() != expected_nbytes {
        return Err(Error::Default("Invalid coefficients".to_string()));
    }

    let mut index = 0;
    let power_basis_coefficients: Vec<u64> = ctx
        .q
        .iter()
        .flat_map(|qi| {
            let size = qi.serialization_length(degree);
            let v = qi.deserialize_vec(&value.coefficients[index..index + size]);
            index += size;
            v
        })
        .collect();

    Ok((
        representation_from_proto,
        power_basis_coefficients,
        variable_time,
    ))
}

impl TryConvertFrom<&Rq> for Poly<PowerBasis> {
    fn try_convert_from(value: &Rq, ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        let (representation_from_proto, coefficients, variable_time) =
            parse_proto(value, ctx, variable_time)?;
        if representation_from_proto != Representation::PowerBasis {
            return Err(Error::Default(
                "The representation asked for does not match the representation in the serialization".to_string(),
            ));
        }
        Poly::<PowerBasis>::try_convert_from(coefficients, ctx, variable_time)
    }
}

impl TryConvertFrom<&Rq> for Poly<Ntt> {
    fn try_convert_from(value: &Rq, ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        let (representation_from_proto, coefficients, variable_time) =
            parse_proto(value, ctx, variable_time)?;
        if representation_from_proto != Representation::Ntt {
            return Err(Error::Default(
                "The representation asked for does not match the representation in the serialization".to_string(),
            ));
        }
        let p = Poly::<PowerBasis>::try_convert_from(coefficients, ctx, variable_time)?;
        Ok(p.into_ntt())
    }
}

impl TryConvertFrom<&Rq> for Poly<NttShoup> {
    fn try_convert_from(value: &Rq, ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        let (representation_from_proto, coefficients, variable_time) =
            parse_proto(value, ctx, variable_time)?;
        if representation_from_proto != Representation::NttShoup {
            return Err(Error::Default(
                "The representation asked for does not match the representation in the serialization".to_string(),
            ));
        }
        let p = Poly::<PowerBasis>::try_convert_from(coefficients, ctx, variable_time)?;
        Ok(p.into_ntt_shoup())
    }
}

impl TryConvertFrom<Vec<u64>> for Poly<PowerBasis> {
    fn try_convert_from(mut v: Vec<u64>, ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        if v.len() == ctx.q.len() * ctx.degree {
            let coefficients = Array2::from_shape_vec((ctx.q.len(), ctx.degree), v).unwrap();
            Ok(Self {
                ctx: ctx.clone(),
                allow_variable_time_computations: variable_time,
                coefficients,
                coefficients_shoup: None,
                has_lazy_coefficients: false,
                _repr: std::marker::PhantomData,
            })
        } else if v.len() <= ctx.degree {
            let mut out = Self::zero(ctx);
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
                izip!(out.coefficients.outer_iter_mut(), ctx.q.iter()).for_each(|(mut w, qi)| {
                    let wi = w.as_slice_mut().unwrap();
                    wi[..v.len()].copy_from_slice(&v);
                    qi.reduce_vec(wi);
                });
                v.zeroize();
            }
            Ok(out)
        } else {
            Err(Error::Default(
                "In PowerBasis representation, either all coefficients must be specified, or only coefficients up to the degree".to_string(),
            ))
        }
    }
}

impl TryConvertFrom<Vec<u64>> for Poly<Ntt> {
    fn try_convert_from(v: Vec<u64>, ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        if let Ok(coefficients) = Array2::from_shape_vec((ctx.q.len(), ctx.degree), v) {
            Ok(Self {
                ctx: ctx.clone(),
                allow_variable_time_computations: variable_time,
                coefficients,
                coefficients_shoup: None,
                has_lazy_coefficients: false,
                _repr: std::marker::PhantomData,
            })
        } else {
            Err(Error::Default(
                "In Ntt representation, all coefficients must be specified".to_string(),
            ))
        }
    }
}

impl TryConvertFrom<Vec<u64>> for Poly<NttShoup> {
    fn try_convert_from(v: Vec<u64>, ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        if let Ok(coefficients) = Array2::from_shape_vec((ctx.q.len(), ctx.degree), v) {
            let mut p = Self {
                ctx: ctx.clone(),
                allow_variable_time_computations: variable_time,
                coefficients,
                coefficients_shoup: None,
                has_lazy_coefficients: false,
                _repr: std::marker::PhantomData,
            };
            p.compute_coefficients_shoup();
            Ok(p)
        } else {
            Err(Error::Default(
                "In NttShoup representation, all coefficients must be specified".to_string(),
            ))
        }
    }
}

impl TryConvertFrom<Array2<u64>> for Poly<PowerBasis> {
    fn try_convert_from(a: Array2<u64>, ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        if a.shape() != [ctx.q.len(), ctx.degree] {
            Err(Error::Default(
                "The array of coefficient does not have the correct shape".to_string(),
            ))
        } else {
            Ok(Self {
                ctx: ctx.clone(),
                allow_variable_time_computations: variable_time,
                coefficients: a,
                coefficients_shoup: None,
                has_lazy_coefficients: false,
                _repr: std::marker::PhantomData,
            })
        }
    }
}

impl TryConvertFrom<Array2<u64>> for Poly<Ntt> {
    fn try_convert_from(a: Array2<u64>, ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        if a.shape() != [ctx.q.len(), ctx.degree] {
            Err(Error::Default(
                "The array of coefficient does not have the correct shape".to_string(),
            ))
        } else {
            Ok(Self {
                ctx: ctx.clone(),
                allow_variable_time_computations: variable_time,
                coefficients: a,
                coefficients_shoup: None,
                has_lazy_coefficients: false,
                _repr: std::marker::PhantomData,
            })
        }
    }
}

impl TryConvertFrom<Array2<u64>> for Poly<NttShoup> {
    fn try_convert_from(a: Array2<u64>, ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        if a.shape() != [ctx.q.len(), ctx.degree] {
            Err(Error::Default(
                "The array of coefficient does not have the correct shape".to_string(),
            ))
        } else {
            let mut p = Self {
                ctx: ctx.clone(),
                allow_variable_time_computations: variable_time,
                coefficients: a,
                coefficients_shoup: None,
                has_lazy_coefficients: false,
                _repr: std::marker::PhantomData,
            };
            p.compute_coefficients_shoup();
            Ok(p)
        }
    }
}

impl<'a> TryConvertFrom<&'a [u64]> for Poly<PowerBasis> {
    fn try_convert_from(v: &'a [u64], ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        Poly::<PowerBasis>::try_convert_from(v.to_vec(), ctx, variable_time)
    }
}

impl<'a> TryConvertFrom<&'a [u64]> for Poly<Ntt> {
    fn try_convert_from(v: &'a [u64], ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        Poly::<Ntt>::try_convert_from(v.to_vec(), ctx, variable_time)
    }
}

impl<'a> TryConvertFrom<&'a [u64]> for Poly<NttShoup> {
    fn try_convert_from(v: &'a [u64], ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        Poly::<NttShoup>::try_convert_from(v.to_vec(), ctx, variable_time)
    }
}

impl<'a> TryConvertFrom<&'a [i64]> for Poly<PowerBasis> {
    fn try_convert_from(v: &'a [i64], ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        if v.len() <= ctx.degree {
            let mut out = Self::zero(ctx);
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
            Err(Error::Default(
                "In PowerBasis representation with signed integers, only `degree` coefficients can be specified".to_string(),
            ))
        }
    }
}

impl<'a> TryConvertFrom<&'a Vec<i64>> for Poly<PowerBasis> {
    fn try_convert_from(v: &'a Vec<i64>, ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        Poly::try_convert_from(v.as_ref() as &[i64], ctx, variable_time)
    }
}

impl<'a> TryConvertFrom<&'a [BigUint]> for Poly<PowerBasis> {
    fn try_convert_from(v: &'a [BigUint], ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        if v.len() > ctx.degree {
            Err(Error::Default(
                "The slice contains too many big integers compared to the polynomial degree"
                    .to_string(),
            ))
        } else {
            let mut coefficients = Array2::zeros((ctx.q.len(), ctx.degree));

            izip!(coefficients.axis_iter_mut(Axis(1)), v).for_each(|(mut c, vi)| {
                c.assign(&ArrayView::from(&ctx.rns.project(vi)));
            });

            Ok(Self {
                ctx: ctx.clone(),
                allow_variable_time_computations: variable_time,
                coefficients,
                coefficients_shoup: None,
                has_lazy_coefficients: false,
                _repr: std::marker::PhantomData,
            })
        }
    }
}

impl<'a> TryConvertFrom<&'a [BigUint]> for Poly<Ntt> {
    fn try_convert_from(v: &'a [BigUint], ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        let p = Poly::<PowerBasis>::try_convert_from(v, ctx, variable_time)?;
        Ok(p.into_ntt())
    }
}

impl<'a> TryConvertFrom<&'a [BigUint]> for Poly<NttShoup> {
    fn try_convert_from(v: &'a [BigUint], ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        let p = Poly::<PowerBasis>::try_convert_from(v, ctx, variable_time)?;
        Ok(p.into_ntt_shoup())
    }
}

impl<'a> TryConvertFrom<&'a Vec<u64>> for Poly<PowerBasis> {
    fn try_convert_from(v: &'a Vec<u64>, ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        Poly::try_convert_from(v.to_vec(), ctx, variable_time)
    }
}

impl<'a> TryConvertFrom<&'a Vec<u64>> for Poly<Ntt> {
    fn try_convert_from(v: &'a Vec<u64>, ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        Poly::try_convert_from(v.to_vec(), ctx, variable_time)
    }
}

impl<'a> TryConvertFrom<&'a Vec<u64>> for Poly<NttShoup> {
    fn try_convert_from(v: &'a Vec<u64>, ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        Poly::try_convert_from(v.to_vec(), ctx, variable_time)
    }
}

impl<'a, const N: usize> TryConvertFrom<&'a [u64; N]> for Poly<PowerBasis> {
    fn try_convert_from(v: &'a [u64; N], ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        Poly::try_convert_from(v.as_ref(), ctx, variable_time)
    }
}

impl<'a, const N: usize> TryConvertFrom<&'a [u64; N]> for Poly<Ntt> {
    fn try_convert_from(v: &'a [u64; N], ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        Poly::try_convert_from(v.as_ref(), ctx, variable_time)
    }
}

impl<'a, const N: usize> TryConvertFrom<&'a [u64; N]> for Poly<NttShoup> {
    fn try_convert_from(v: &'a [u64; N], ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        Poly::try_convert_from(v.as_ref(), ctx, variable_time)
    }
}

impl<'a, const N: usize> TryConvertFrom<&'a [BigUint; N]> for Poly<PowerBasis> {
    fn try_convert_from(
        v: &'a [BigUint; N],
        ctx: &Arc<Context>,
        variable_time: bool,
    ) -> Result<Self> {
        Poly::try_convert_from(v.as_ref(), ctx, variable_time)
    }
}

impl<'a, const N: usize> TryConvertFrom<&'a [BigUint; N]> for Poly<Ntt> {
    fn try_convert_from(
        v: &'a [BigUint; N],
        ctx: &Arc<Context>,
        variable_time: bool,
    ) -> Result<Self> {
        Poly::try_convert_from(v.as_ref(), ctx, variable_time)
    }
}

impl<'a, const N: usize> TryConvertFrom<&'a [BigUint; N]> for Poly<NttShoup> {
    fn try_convert_from(
        v: &'a [BigUint; N],
        ctx: &Arc<Context>,
        variable_time: bool,
    ) -> Result<Self> {
        Poly::try_convert_from(v.as_ref(), ctx, variable_time)
    }
}

impl<'a, const N: usize> TryConvertFrom<&'a [i64; N]> for Poly<PowerBasis> {
    fn try_convert_from(v: &'a [i64; N], ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        Poly::try_convert_from(v.as_ref(), ctx, variable_time)
    }
}

impl TryFrom<&Poly<PowerBasis>> for Vec<u64> {
    type Error = Error;

    fn try_from(p: &Poly<PowerBasis>) -> Result<Self> {
        p.coefficients
            .as_slice()
            .ok_or_else(|| {
                Error::Default("Polynomial coefficients are not contiguous in memory".to_string())
            })
            .map(|slice| slice.to_vec())
    }
}

impl TryFrom<&Poly<Ntt>> for Vec<u64> {
    type Error = Error;

    fn try_from(p: &Poly<Ntt>) -> Result<Self> {
        p.coefficients
            .as_slice()
            .ok_or_else(|| {
                Error::Default("Polynomial coefficients are not contiguous in memory".to_string())
            })
            .map(|slice| slice.to_vec())
    }
}

impl TryFrom<&Poly<NttShoup>> for Vec<u64> {
    type Error = Error;

    fn try_from(p: &Poly<NttShoup>) -> Result<Self> {
        p.coefficients
            .as_slice()
            .ok_or_else(|| {
                Error::Default("Polynomial coefficients are not contiguous in memory".to_string())
            })
            .map(|slice| slice.to_vec())
    }
}

impl From<&Poly<PowerBasis>> for Vec<BigUint> {
    fn from(p: &Poly<PowerBasis>) -> Self {
        izip!(p.coefficients.axis_iter(Axis(1)))
            .map(|c| p.ctx.rns.lift(c))
            .collect_vec()
    }
}

impl From<&Poly<Ntt>> for Vec<BigUint> {
    fn from(p: &Poly<Ntt>) -> Self {
        izip!(p.coefficients.axis_iter(Axis(1)))
            .map(|c| p.ctx.rns.lift(c))
            .collect_vec()
    }
}

impl From<&Poly<NttShoup>> for Vec<BigUint> {
    fn from(p: &Poly<NttShoup>) -> Self {
        izip!(p.coefficients.axis_iter(Axis(1)))
            .map(|c| p.ctx.rns.lift(c))
            .collect_vec()
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        Error as CrateError,
        proto::rq::Rq,
        rq::{Context, Ntt, NttShoup, Poly, PowerBasis, traits::TryConvertFrom},
    };
    use num_bigint::BigUint;
    use rand::rng;
    use std::{error::Error, sync::Arc};

    static MODULI: &[u64; 3] = &[1153, 4611686018326724609, 4611686018309947393];

    #[test]
    fn proto() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        for modulus in MODULI {
            let ctx = Arc::new(Context::new(&[*modulus], 16)?);
            let p = Poly::<PowerBasis>::random(&ctx, &mut rng);
            let proto = Rq::from(&p);
            assert_eq!(
                Poly::<PowerBasis>::try_convert_from(&proto, &ctx, false)?,
                p
            );
            assert_eq!(
                Poly::<Ntt>::try_convert_from(&proto, &ctx, false).unwrap_err(),
                CrateError::Default(
                    "The representation asked for does not match the representation in the serialization".to_string()
                )
            );
            assert_eq!(
                Poly::<NttShoup>::try_convert_from(&proto, &ctx, false).unwrap_err(),
                CrateError::Default(
                    "The representation asked for does not match the representation in the serialization".to_string()
                )
            );
        }

        let ctx = Arc::new(Context::new(MODULI, 16)?);
        let p = Poly::<Ntt>::random(&ctx, &mut rng);
        let proto = Rq::from(&p);
        assert_eq!(Poly::<Ntt>::try_convert_from(&proto, &ctx, false)?, p);

        let p = Poly::<NttShoup>::random(&ctx, &mut rng);
        let proto = Rq::from(&p);
        assert_eq!(Poly::<NttShoup>::try_convert_from(&proto, &ctx, false)?, p);

        Ok(())
    }

    #[test]
    fn try_convert_from_slice_zero() -> Result<(), Box<dyn Error>> {
        for modulus in MODULI {
            let ctx = Arc::new(Context::new(&[*modulus], 16)?);

            // Power Basis
            assert_eq!(
                Poly::<PowerBasis>::try_convert_from(&[0u64], &ctx, false)?,
                Poly::<PowerBasis>::zero(&ctx)
            );
            assert_eq!(
                Poly::<PowerBasis>::try_convert_from(&[0i64], &ctx, false)?,
                Poly::<PowerBasis>::zero(&ctx)
            );
            assert_eq!(
                Poly::<PowerBasis>::try_convert_from(&[0u64; 16], &ctx, false)?,
                Poly::<PowerBasis>::zero(&ctx)
            );
            assert_eq!(
                Poly::<PowerBasis>::try_convert_from(&[0i64; 16], &ctx, false)?,
                Poly::<PowerBasis>::zero(&ctx)
            );
            assert!(Poly::<PowerBasis>::try_convert_from(&[0u64; 17], &ctx, false).is_err());

            // Ntt
            assert!(Poly::<Ntt>::try_convert_from(&[0u64], &ctx, false).is_err());
            assert!(Poly::<Ntt>::try_convert_from(&[0u64; 16], &ctx, false).is_ok());
            assert!(Poly::<Ntt>::try_convert_from(&[0u64; 17], &ctx, false).is_err());
        }

        let ctx = Arc::new(Context::new(MODULI, 16)?);
        assert_eq!(
            Poly::<PowerBasis>::try_convert_from(Vec::<u64>::default(), &ctx, false)?,
            Poly::<PowerBasis>::zero(&ctx)
        );
        assert!(Poly::<Ntt>::try_convert_from(Vec::<u64>::default(), &ctx, false).is_err());

        Ok(())
    }

    #[test]
    fn try_convert_from_vec_zero() -> Result<(), Box<dyn Error>> {
        for modulus in MODULI {
            let ctx = Arc::new(Context::new(&[*modulus], 16)?);
            assert_eq!(
                Poly::<PowerBasis>::try_convert_from(vec![], &ctx, false)?,
                Poly::<PowerBasis>::zero(&ctx)
            );
            assert!(Poly::<Ntt>::try_convert_from(vec![], &ctx, false).is_err());

            assert_eq!(
                Poly::<PowerBasis>::try_convert_from(vec![0], &ctx, false)?,
                Poly::<PowerBasis>::zero(&ctx)
            );
            assert!(Poly::<Ntt>::try_convert_from(vec![0], &ctx, false).is_err());

            assert_eq!(
                Poly::<PowerBasis>::try_convert_from(vec![0; 16], &ctx, false)?,
                Poly::<PowerBasis>::zero(&ctx)
            );
            assert_eq!(
                Poly::<Ntt>::try_convert_from(vec![0; 16], &ctx, false)?,
                Poly::<Ntt>::zero(&ctx)
            );
        }

        Ok(())
    }

    #[test]
    fn biguint() -> Result<(), Box<dyn Error>> {
        let mut rng = rng();
        let ctx = Arc::new(Context::new(MODULI, 16)?);
        let p = Poly::<PowerBasis>::random(&ctx, &mut rng);
        let values = Vec::<BigUint>::from(&p);
        let p2 = Poly::<PowerBasis>::try_convert_from(values.as_slice(), &ctx, false)?;
        assert_eq!(p, p2);
        Ok(())
    }
}
