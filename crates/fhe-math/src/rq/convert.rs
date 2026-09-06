//! Implementation of conversions from and to polynomials.

use super::{
    Context, Ntt, NttShoup, Poly, PowerBasis, Representation, RepresentationTag,
    traits::TryConvertFrom,
};
use crate::{
    Error, PolynomialSerializationError, Result,
    proto::rq::{Representation as RepresentationProto, Rq},
};
use itertools::{Itertools, izip};
use ndarray::{Array2, ArrayView, Axis};
use num_bigint::BigUint;
use std::sync::Arc;
use zeroize::{Zeroize, Zeroizing};

#[expect(
    clippy::fallible_impl_from,
    reason = "Poly constructors guarantee contiguous standard-layout coefficient rows"
)]
impl<R: RepresentationTag> From<&Poly<R>> for Rq {
    fn from(p: &Poly<R>) -> Self {
        assert!(!p.has_lazy_coefficients);
        // Only transformed coefficients need scratch. In particular, never clone
        // the Shoup cache just to discard it during conversion to power basis.
        let mut scratch = if R::REPRESENTATION == Representation::PowerBasis {
            None
        } else {
            Some(Zeroizing::new(Poly::<PowerBasis> {
                ctx: p.ctx.clone(),
                coefficients: p.coefficients.clone(),
                coefficients_shoup: None,
                has_lazy_coefficients: false,
                allow_variable_time_computations: p.allow_variable_time_computations,
                _repr: std::marker::PhantomData,
            }))
        };
        if let Some(q) = scratch.as_mut() {
            q.ntt_backward();
        }
        let coefficients = scratch
            .as_ref()
            .map_or(&p.coefficients, |q| &q.coefficients);

        let mut proto = Rq::default();
        match R::REPRESENTATION {
            Representation::PowerBasis => {
                proto.representation = RepresentationProto::Powerbasis as i32
            }
            Representation::Ntt => proto.representation = RepresentationProto::Ntt as i32,
            Representation::NttShoup => proto.representation = RepresentationProto::Nttshoup as i32,
        }
        let capacity = p
            .ctx
            .q
            .iter()
            .map(|qi| qi.serialization_length(p.ctx.degree))
            .sum();
        let mut serialization = Vec::with_capacity(capacity);
        for (row, qi) in coefficients.outer_iter().zip(p.ctx.q.iter()) {
            qi.serialize_into(row.as_slice().unwrap(), &mut serialization);
        }
        proto.coefficients = serialization;
        proto.degree = p.ctx.degree as u32;
        // Timing policy is local execution state, not serialized data. In
        // particular, untrusted bytes must not authorize variable-time work.
        proto.allow_variable_time = false;
        proto
    }
}

fn parse_proto(
    value: &Rq,
    ctx: &Arc<Context>,
    variable_time: bool,
) -> Result<(Representation, Poly<PowerBasis>)> {
    let repr = value.representation.try_into().map_err(|_| {
        PolynomialSerializationError::InvalidRepresentation {
            value: value.representation,
        }
    })?;
    let representation_from_proto = match repr {
        RepresentationProto::Powerbasis => Representation::PowerBasis,
        RepresentationProto::Ntt => Representation::Ntt,
        RepresentationProto::Nttshoup => Representation::NttShoup,
        RepresentationProto::Unknown => {
            return Err(PolynomialSerializationError::UnknownRepresentation.into());
        }
    };

    let degree = value.degree as usize;
    if !degree.is_multiple_of(8) || degree < 8 {
        return Err(PolynomialSerializationError::InvalidDegree { degree }.into());
    }

    if degree != ctx.degree {
        return Err(Error::DegreeMismatch {
            found: degree,
            expected: ctx.degree,
        });
    }

    let mut expected_nbytes = 0;
    ctx.q
        .iter()
        .for_each(|qi| expected_nbytes += qi.serialization_length(degree));
    if value.coefficients.len() != expected_nbytes {
        return Err(PolynomialSerializationError::InvalidCoefficientCount {
            actual: value.coefficients.len(),
            expected: expected_nbytes,
        }
        .into());
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

    for (row, modulus) in power_basis_coefficients
        .chunks_exact(degree)
        .zip(ctx.q.iter())
    {
        if row.iter().any(|coefficient| *coefficient >= **modulus) {
            return Err(PolynomialSerializationError::NonCanonicalCoefficient {
                modulus: **modulus,
            }
            .into());
        }
    }

    // Only this parser constructs directly from wire coefficients: the exact
    // shape and canonical residues have been validated above. Public raw
    // constructors must still normalize arbitrary input.
    let coefficients =
        Array2::from_shape_vec((ctx.q.len(), degree), power_basis_coefficients).unwrap();
    let polynomial = from_canonical_coefficients(coefficients, ctx, variable_time);
    Ok((representation_from_proto, polynomial))
}

impl TryConvertFrom<&Rq> for Poly<PowerBasis> {
    fn try_convert_from(value: &Rq, ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        let (representation_from_proto, p) = parse_proto(value, ctx, variable_time)?;
        if representation_from_proto != Representation::PowerBasis {
            return Err(PolynomialSerializationError::RepresentationMismatch {
                found: representation_from_proto,
                expected: Representation::PowerBasis,
            }
            .into());
        }
        Ok(p)
    }
}

impl TryConvertFrom<&Rq> for Poly<Ntt> {
    fn try_convert_from(value: &Rq, ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        let (representation_from_proto, p) = parse_proto(value, ctx, variable_time)?;
        if representation_from_proto != Representation::Ntt {
            return Err(PolynomialSerializationError::RepresentationMismatch {
                found: representation_from_proto,
                expected: Representation::Ntt,
            }
            .into());
        }
        Ok(p.into_ntt())
    }
}

impl TryConvertFrom<&Rq> for Poly<NttShoup> {
    fn try_convert_from(value: &Rq, ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        let (representation_from_proto, p) = parse_proto(value, ctx, variable_time)?;
        if representation_from_proto != Representation::NttShoup {
            return Err(PolynomialSerializationError::RepresentationMismatch {
                found: representation_from_proto,
                expected: Representation::NttShoup,
            }
            .into());
        }
        Ok(p.into_ntt_shoup())
    }
}

// Preserve logical element order and reduce all RNS rows in constant time.
fn canonical_coefficients(a: Array2<u64>, ctx: &Context) -> Array2<u64> {
    let mut a = if a.is_standard_layout() {
        a
    } else {
        a.as_standard_layout().into_owned()
    };
    for (mut row, modulus) in a.outer_iter_mut().zip(ctx.q.iter()) {
        modulus.reduce_vec(row.as_slice_mut().unwrap());
    }
    a
}

// Common initialization for canonical, standard-layout coefficients. Callers
// must validate dimensions and either normalize or reject noncanonical input.
fn from_canonical_coefficients<R: RepresentationTag>(
    coefficients: Array2<u64>,
    ctx: &Arc<Context>,
    variable_time: bool,
) -> Poly<R> {
    debug_assert_eq!(coefficients.shape(), [ctx.q.len(), ctx.degree]);
    debug_assert!(coefficients.is_standard_layout());
    let mut p = Poly {
        ctx: ctx.clone(),
        coefficients,
        coefficients_shoup: None,
        has_lazy_coefficients: false,
        allow_variable_time_computations: variable_time,
        _repr: std::marker::PhantomData,
    };
    if R::REPRESENTATION == Representation::NttShoup {
        p.compute_coefficients_shoup();
    }
    p
}

fn from_rns_array<R: RepresentationTag>(
    a: Array2<u64>,
    ctx: &Arc<Context>,
    variable_time: bool,
) -> Result<Poly<R>> {
    if a.shape() != [ctx.q.len(), ctx.degree] {
        return Err(Error::InvalidCoefficientShape {
            actual_rows: a.nrows(),
            actual_columns: a.ncols(),
            expected_rows: ctx.q.len(),
            expected_columns: ctx.degree,
        });
    }
    Ok(from_canonical_coefficients(
        canonical_coefficients(a, ctx),
        ctx,
        variable_time,
    ))
}

fn from_rns_vec<R: RepresentationTag>(
    v: Vec<u64>,
    ctx: &Arc<Context>,
    variable_time: bool,
) -> Result<Poly<R>> {
    let actual = v.len();
    let coefficients = Array2::from_shape_vec((ctx.q.len(), ctx.degree), v).map_err(|_| {
        Error::InvalidCoefficientCount {
            representation: R::REPRESENTATION,
            actual,
            degree: ctx.degree,
            moduli: ctx.q.len(),
        }
    })?;
    from_rns_array(coefficients, ctx, variable_time)
}

impl TryConvertFrom<Vec<u64>> for Poly<PowerBasis> {
    fn try_convert_from(mut v: Vec<u64>, ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        if v.len() == ctx.q.len() * ctx.degree {
            from_rns_vec(v, ctx, variable_time)
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
                    out.allow_variable_time_computations(fhe_traits::VariableTime::new(
                        fhe_traits::PublicData::assert_public(),
                    ));
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
            Err(Error::InvalidCoefficientCount {
                representation: Representation::PowerBasis,
                actual: v.len(),
                degree: ctx.degree,
                moduli: ctx.q.len(),
            })
        }
    }
}

impl TryConvertFrom<Vec<u64>> for Poly<Ntt> {
    fn try_convert_from(v: Vec<u64>, ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        from_rns_vec(v, ctx, variable_time)
    }
}

impl TryConvertFrom<Vec<u64>> for Poly<NttShoup> {
    fn try_convert_from(v: Vec<u64>, ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        from_rns_vec(v, ctx, variable_time)
    }
}

impl TryConvertFrom<Array2<u64>> for Poly<PowerBasis> {
    fn try_convert_from(a: Array2<u64>, ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        from_rns_array(a, ctx, variable_time)
    }
}

impl TryConvertFrom<Array2<u64>> for Poly<Ntt> {
    fn try_convert_from(a: Array2<u64>, ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        from_rns_array(a, ctx, variable_time)
    }
}

impl TryConvertFrom<Array2<u64>> for Poly<NttShoup> {
    fn try_convert_from(a: Array2<u64>, ctx: &Arc<Context>, variable_time: bool) -> Result<Self> {
        from_rns_array(a, ctx, variable_time)
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
                out.allow_variable_time_computations(fhe_traits::VariableTime::new(
                    fhe_traits::PublicData::assert_public(),
                ));
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
            Err(Error::TooManyCoefficients {
                actual: v.len(),
                maximum: ctx.degree,
            })
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
            Err(Error::TooManyCoefficients {
                actual: v.len(),
                maximum: ctx.degree,
            })
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
            .ok_or(Error::NonContiguousCoefficients)
            .map(|slice| slice.to_vec())
    }
}

impl TryFrom<&Poly<Ntt>> for Vec<u64> {
    type Error = Error;

    fn try_from(p: &Poly<Ntt>) -> Result<Self> {
        p.coefficients
            .as_slice()
            .ok_or(Error::NonContiguousCoefficients)
            .map(|slice| slice.to_vec())
    }
}

impl TryFrom<&Poly<NttShoup>> for Vec<u64> {
    type Error = Error;

    fn try_from(p: &Poly<NttShoup>) -> Result<Self> {
        p.coefficients
            .as_slice()
            .ok_or(Error::NonContiguousCoefficients)
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
        Error as CrateError, PolynomialSerializationError,
        proto::rq::Rq,
        rq::{Context, Ntt, NttShoup, Poly, PowerBasis, traits::TryConvertFrom},
    };
    use num_bigint::BigUint;
    use rand::rng;
    use std::{error::Error, sync::Arc};

    static MODULI: &[u64; 3] = &[1153, 4611686018326724609, 4611686018309947393];

    #[test]
    fn shared_constructors_reject_invalid_shapes_and_lengths() -> Result<(), Box<dyn Error>> {
        let ctx = Context::new_arc(&[1153, 2017], 16)?;
        macro_rules! check {
            ($repr:ty) => {
                for shape in [(1, 32), (2, 15), (3, 16)] {
                    let error =
                        Poly::<$repr>::try_convert_from(ndarray::Array2::zeros(shape), &ctx, false)
                            .unwrap_err();
                    assert!(matches!(
                        error,
                        crate::Error::InvalidCoefficientShape { .. }
                    ));
                }
                for length in [17, 31, 33] {
                    let error =
                        Poly::<$repr>::try_convert_from(vec![0; length], &ctx, false).unwrap_err();
                    assert!(matches!(
                        error,
                        crate::Error::InvalidCoefficientCount { .. }
                    ));
                }
            };
        }
        check!(PowerBasis);
        check!(Ntt);
        check!(NttShoup);
        // Power-basis short inputs still denote a polynomial, not flattened RNS rows.
        let short = Poly::<PowerBasis>::try_convert_from(vec![2018], &ctx, false)?;
        assert_eq!(
            Vec::<BigUint>::from(&short).first(),
            Some(&BigUint::from(2018u64))
        );
        Ok(())
    }

    #[test]
    fn owned_array_layouts_preserve_values_and_ntt_safety() -> Result<(), Box<dyn Error>> {
        use ndarray::{Array2, Axis};
        for moduli in [&MODULI[..1], &MODULI[..2]] {
            let ctx = Context::new_arc(moduli, 16)?;
            let base =
                Array2::from_shape_fn((moduli.len(), 16), |(row, col)| (row * 16 + col) as u64);
            let mut reverse_columns = base.clone();
            reverse_columns.invert_axis(Axis(1));
            let mut reverse_rows = base.clone();
            reverse_rows.invert_axis(Axis(0));
            let transposed =
                Array2::from_shape_fn((16, moduli.len()), |(col, row)| (row * 16 + col) as u64)
                    .reversed_axes();
            for array in [reverse_columns, reverse_rows, transposed] {
                let standard = array.as_standard_layout().into_owned();
                for public in [false, true] {
                    let pb = Poly::<PowerBasis>::try_convert_from(array.clone(), &ctx, public)?;
                    let expected =
                        Poly::<PowerBasis>::try_convert_from(standard.clone(), &ctx, false)?;
                    assert!(pb.coefficients().is_standard_layout());
                    assert_eq!(pb.clone().into_ntt().into_power_basis(), expected);
                    assert_eq!(pb.into_ntt_shoup().into_power_basis(), expected);
                    let ntt = Poly::<Ntt>::try_convert_from(array.clone(), &ctx, public)?;
                    let expected = Poly::<Ntt>::try_convert_from(standard.clone(), &ctx, false)?;
                    assert!(ntt.coefficients().is_standard_layout());
                    assert_eq!(ntt.into_power_basis(), expected.into_power_basis());
                    let shoup = Poly::<NttShoup>::try_convert_from(array.clone(), &ctx, public)?;
                    let expected =
                        Poly::<NttShoup>::try_convert_from(standard.clone(), &ctx, false)?;
                    assert!(shoup.coefficients().is_standard_layout());
                    assert_eq!(shoup.into_power_basis(), expected.into_power_basis());
                }
            }
        }
        Ok(())
    }

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
                CrateError::PolynomialSerialization(
                    PolynomialSerializationError::RepresentationMismatch {
                        found: crate::rq::Representation::PowerBasis,
                        expected: crate::rq::Representation::Ntt,
                    }
                )
            );
            assert_eq!(
                Poly::<NttShoup>::try_convert_from(&proto, &ctx, false).unwrap_err(),
                CrateError::PolynomialSerialization(
                    PolynomialSerializationError::RepresentationMismatch {
                        found: crate::rq::Representation::PowerBasis,
                        expected: crate::rq::Representation::NttShoup,
                    }
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

    #[test]
    fn wire_requires_matching_degree_and_canonical_coefficients() -> Result<(), Box<dyn Error>> {
        use fhe_traits::{DeserializeWithContext, Serialize};
        use prost::Message;
        let small = Context::new_arc(&[1153], 8)?;
        let large = Context::new_arc(&[1153], 16)?;
        let bytes = Poly::<PowerBasis>::zero(&small).to_bytes();
        assert!(matches!(
            Poly::<PowerBasis>::from_bytes(&bytes, &large),
            Err(CrateError::DegreeMismatch {
                found: 8,
                expected: 16
            })
        ));
        let bytes = Poly::<PowerBasis>::zero(&large).to_bytes();
        assert!(Poly::<PowerBasis>::from_bytes(&bytes, &small).is_err());
        for representation in [1, 2, 3] {
            let proto = Rq {
                representation,
                degree: 16,
                coefficients: vec![255; 22],
                allow_variable_time: false,
            };
            let bytes = proto.encode_to_vec();
            let error = match representation {
                1 => Poly::<PowerBasis>::from_bytes(&bytes, &large).unwrap_err(),
                2 => Poly::<Ntt>::from_bytes(&bytes, &large).unwrap_err(),
                _ => Poly::<NttShoup>::from_bytes(&bytes, &large).unwrap_err(),
            };
            assert!(matches!(
                error,
                CrateError::PolynomialSerialization(
                    PolynomialSerializationError::NonCanonicalCoefficient { modulus: 1153 }
                )
            ));
        }
        Ok(())
    }

    #[test]
    fn raw_rns_constructors_reduce_every_row() -> Result<(), Box<dyn Error>> {
        use ndarray::Array2;
        for moduli in [&MODULI[..1], &MODULI[..2]] {
            let ctx = Context::new_arc(moduli, 16)?;
            let values = vec![u64::MAX; moduli.len() * 16];
            let expected =
                Array2::from_shape_fn((moduli.len(), 16), |(row, _)| u64::MAX % moduli[row]);
            let array = Array2::from_shape_vec((moduli.len(), 16), values.clone())?;
            for public in [false, true] {
                let pb = Poly::<PowerBasis>::try_convert_from(values.clone(), &ctx, public)?;
                assert_eq!(pb.coefficients(), expected);
                assert_eq!(
                    (&pb + &Poly::<PowerBasis>::zero(&ctx)).coefficients(),
                    expected
                );
                assert_eq!(
                    Poly::<PowerBasis>::try_convert_from(array.clone(), &ctx, public)?
                        .coefficients(),
                    expected
                );
                assert_eq!(
                    Poly::<Ntt>::try_convert_from(values.clone(), &ctx, public)?.coefficients(),
                    expected
                );
                assert_eq!(
                    Poly::<Ntt>::try_convert_from(array.clone(), &ctx, public)?.coefficients(),
                    expected
                );
                let shoup = Poly::<NttShoup>::try_convert_from(values.clone(), &ctx, public)?;
                assert_eq!(shoup.coefficients(), expected);
                assert_eq!(
                    Poly::<NttShoup>::try_convert_from(array.clone(), &ctx, public)?,
                    shoup
                );
                assert_eq!(pb.clone().into_ntt().into_power_basis(), pb);
            }
        }
        Ok(())
    }
}
