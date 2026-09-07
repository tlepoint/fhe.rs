//! Implementation of conversions from and to polynomials.

use super::{
    Context, Ntt, NttShoup, Poly, PowerBasis, Representation, RepresentationTag, wire::FromProto,
};
use crate::{
    Error, PolynomialSerializationError, Result,
    proto::rq::{Representation as RepresentationProto, Rq},
};
use itertools::{Itertools, izip};
use ndarray::{Array2, ArrayView, Axis};
use num_bigint::BigUint;
use std::sync::Arc;
use zeroize::Zeroizing;

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
            let v = qi
                .deserialize_vec(&value.coefficients[index..index + size])
                .unwrap();
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

impl FromProto<&Rq> for Poly<PowerBasis> {
    fn from_proto_with_timing(
        value: &Rq,
        ctx: &Arc<Context>,
        permission: Option<fhe_util::VariableTime>,
    ) -> Result<Self> {
        let variable_time = permission.is_some();
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

impl FromProto<&Rq> for Poly<Ntt> {
    fn from_proto_with_timing(
        value: &Rq,
        ctx: &Arc<Context>,
        permission: Option<fhe_util::VariableTime>,
    ) -> Result<Self> {
        let variable_time = permission.is_some();
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

impl FromProto<&Rq> for Poly<NttShoup> {
    fn from_proto_with_timing(
        value: &Rq,
        ctx: &Arc<Context>,
        permission: Option<fhe_util::VariableTime>,
    ) -> Result<Self> {
        let variable_time = permission.is_some();
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

impl<R: RepresentationTag> Poly<R> {
    /// Import raw RNS residues in representation `R`; no transform is
    /// performed. Shape must be `(number of moduli, degree)`. Rows are
    /// reduced modulo their respective primes and copied into standard
    /// layout if necessary.
    pub fn from_rns_residues(values: Array2<u64>, ctx: &Arc<Context>) -> Result<Self> {
        Self::from_rns_residues_with_timing(values, ctx, None)
    }

    /// Import raw RNS residues with optional permission for variable-time work.
    /// Supplying a token asserts the residues are public.
    pub fn from_rns_residues_with_timing(
        values: Array2<u64>,
        ctx: &Arc<Context>,
        permission: Option<fhe_util::VariableTime>,
    ) -> Result<Self> {
        from_rns_array(values, ctx, permission.is_some())
    }

    /// Import a contiguous, modulus-major sequence of raw RNS residues in `R`.
    /// Length must be exactly `number of moduli * degree`; no transform occurs.
    pub fn from_rns_slice(values: &[u64], ctx: &Arc<Context>) -> Result<Self> {
        Self::from_rns_slice_with_timing(values, ctx, None)
    }

    /// Import a raw residue slice with optional public-data timing permission.
    pub fn from_rns_slice_with_timing(
        values: &[u64],
        ctx: &Arc<Context>,
        permission: Option<fhe_util::VariableTime>,
    ) -> Result<Self> {
        from_rns_vec(values.to_vec(), ctx, permission.is_some())
    }
}

impl Poly<Ntt> {
    /// Reduce wide NTT accumulators into canonical residues. Accepts strided
    /// views with shape `(number of moduli, degree)`. No NTT is performed.
    pub fn from_wide_ntt_residues(
        a: ndarray::ArrayView2<'_, u128>,
        ctx: &Arc<Context>,
    ) -> Result<Self> {
        Self::from_wide_ntt_residues_with_timing(a, ctx, None)
    }

    /// Reduce wide NTT residues with optional public-data timing permission.
    pub fn from_wide_ntt_residues_with_timing(
        a: ndarray::ArrayView2<'_, u128>,
        ctx: &Arc<Context>,
        permission: Option<fhe_util::VariableTime>,
    ) -> Result<Self> {
        let variable_time = permission.is_some();
        if a.shape() != [ctx.q.len(), ctx.degree] {
            return Err(Error::InvalidCoefficientShape {
                actual_rows: a.nrows(),
                actual_columns: a.ncols(),
                expected_rows: ctx.q.len(),
                expected_columns: ctx.degree,
            });
        }
        let mut coefficients = Array2::zeros((ctx.q.len(), ctx.degree));
        for ((mut output, input), modulus) in coefficients
            .outer_iter_mut()
            .zip(a.outer_iter())
            .zip(ctx.q.iter())
        {
            for (out, value) in output.iter_mut().zip(input.iter()) {
                *out = if variable_time {
                    unsafe { modulus.reduce_u128_vt(*value) }
                } else {
                    modulus.reduce_u128(*value)
                };
            }
        }
        Ok(from_canonical_coefficients(
            coefficients,
            ctx,
            variable_time,
        ))
    }
}

impl Poly<PowerBasis> {
    /// Encode integer coefficients in power basis, reducing each coefficient
    /// modulo every RNS prime and padding to the degree. Rejects longer input;
    /// its length never changes the interpretation to raw residues.
    pub fn from_coefficients(v: &[u64], ctx: &Arc<Context>) -> Result<Self> {
        Self::from_coefficients_with_timing(v, ctx, None)
    }

    /// Encode coefficients with optional public-data timing permission.
    pub fn from_coefficients_with_timing(
        v: &[u64],
        ctx: &Arc<Context>,
        permission: Option<fhe_util::VariableTime>,
    ) -> Result<Self> {
        if v.len() > ctx.degree {
            return Err(Error::TooManyCoefficients {
                actual: v.len(),
                maximum: ctx.degree,
            });
        }
        let mut out = Self::zero(ctx);
        out.allow_variable_time_computations = permission.is_some();
        for (mut row, modulus) in out.coefficients.outer_iter_mut().zip(ctx.q.iter()) {
            let row = row.as_slice_mut().unwrap();
            row[..v.len()].copy_from_slice(v);
            if permission.is_some() {
                unsafe { modulus.reduce_vec_vt(row) };
            } else {
                modulus.reduce_vec(row);
            }
        }
        Ok(out)
    }

    /// Encode signed integer coefficients, reducing and zero-padding to degree.
    pub fn from_signed_coefficients(v: &[i64], ctx: &Arc<Context>) -> Result<Self> {
        Self::from_signed_coefficients_with_timing(v, ctx, None)
    }

    /// Encode signed coefficients with optional public-data timing permission.
    pub fn from_signed_coefficients_with_timing(
        v: &[i64],
        ctx: &Arc<Context>,
        permission: Option<fhe_util::VariableTime>,
    ) -> Result<Self> {
        let variable_time = permission.is_some();
        if v.len() <= ctx.degree {
            let mut out = Self::zero(ctx);
            if variable_time {
                out.allow_variable_time_computations(fhe_util::VariableTime::new(
                    fhe_util::PublicData::assert_public(),
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

    /// Encode arbitrary unsigned integer coefficients, reducing and padding to
    /// degree. Big-integer operations do not provide constant-time guarantees.
    pub fn from_biguint_coefficients(v: &[BigUint], ctx: &Arc<Context>) -> Result<Self> {
        Self::from_biguint_coefficients_with_timing(v, ctx, None)
    }

    /// Encode big-integer coefficients with optional public-data timing
    /// permission.
    pub fn from_biguint_coefficients_with_timing(
        v: &[BigUint],
        ctx: &Arc<Context>,
        permission: Option<fhe_util::VariableTime>,
    ) -> Result<Self> {
        let variable_time = permission.is_some();
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
        rq::{Context, Ntt, NttShoup, Poly, PowerBasis, wire::FromProto},
    };
    use num_bigint::BigUint;
    use rand::rng;
    use std::{error::Error, sync::Arc};

    static MODULI: &[u64; 3] = &[1153, 4611686018326724609, 4611686018309947393];

    #[test]
    fn wide_ntt_conversion_reduces_once_and_accepts_strided_views() -> Result<(), Box<dyn Error>> {
        let ctx = Context::new_arc(&[1153, 4611686018326724609], 16)?;
        let a = ndarray::Array2::from_shape_fn((2, 16), |(r, c)| match c % 4 {
            0 => 0,
            1 => u128::MAX,
            2 => ctx.moduli()[r] as u128,
            _ => ctx.moduli()[r] as u128 + 1,
        });
        for view in [a.view(), a.slice(ndarray::s![.., ..;-1])] {
            for public in [false, true] {
                let p = Poly::<Ntt>::from_wide_ntt_residues_with_timing(
                    view,
                    &ctx,
                    (public).then(|| {
                        fhe_util::VariableTime::new(fhe_util::PublicData::assert_public())
                    }),
                )?;
                assert_eq!(p.allows_variable_time_computations(), public);
                for ((output, input), modulus) in p
                    .coefficients
                    .outer_iter()
                    .zip(view.outer_iter())
                    .zip(ctx.moduli())
                {
                    for (actual, value) in output.iter().zip(input) {
                        assert_eq!(*actual, (*value % *modulus as u128) as u64);
                    }
                }
                assert_eq!(p.clone().into_power_basis().into_ntt(), p);
            }
        }
        let wrong = ndarray::Array2::<u128>::zeros((1, 32));
        assert!(matches!(
            Poly::<Ntt>::from_wide_ntt_residues(wrong.view(), &ctx),
            Err(crate::Error::InvalidCoefficientShape { .. })
        ));
        Ok(())
    }

    #[test]
    fn shared_constructors_reject_invalid_shapes_and_lengths() -> Result<(), Box<dyn Error>> {
        let ctx = Context::new_arc(&[1153, 2017], 16)?;
        macro_rules! check {
            ($repr:ty) => {
                for shape in [(1, 32), (2, 15), (3, 16)] {
                    let error =
                        Poly::<$repr>::from_rns_residues(ndarray::Array2::zeros(shape), &ctx)
                            .unwrap_err();
                    assert!(matches!(
                        error,
                        crate::Error::InvalidCoefficientShape { .. }
                    ));
                }
                for length in [17, 31, 33] {
                    let error = Poly::<$repr>::from_rns_slice(&vec![0; length], &ctx).unwrap_err();
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
        let short = Poly::<PowerBasis>::from_coefficients(&[2018], &ctx)?;
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
                    let pb = Poly::<PowerBasis>::from_rns_residues_with_timing(
                        array.clone(),
                        &ctx,
                        (public).then(|| {
                            fhe_util::VariableTime::new(fhe_util::PublicData::assert_public())
                        }),
                    )?;
                    let expected = Poly::<PowerBasis>::from_rns_residues(standard.clone(), &ctx)?;
                    assert!(pb.coefficients().is_standard_layout());
                    assert_eq!(pb.clone().into_ntt().into_power_basis(), expected);
                    assert_eq!(pb.into_ntt_shoup().into_power_basis(), expected);
                    let ntt = Poly::<Ntt>::from_rns_residues_with_timing(
                        array.clone(),
                        &ctx,
                        (public).then(|| {
                            fhe_util::VariableTime::new(fhe_util::PublicData::assert_public())
                        }),
                    )?;
                    let expected = Poly::<Ntt>::from_rns_residues(standard.clone(), &ctx)?;
                    assert!(ntt.coefficients().is_standard_layout());
                    assert_eq!(ntt.into_power_basis(), expected.into_power_basis());
                    let shoup = Poly::<NttShoup>::from_rns_residues_with_timing(
                        array.clone(),
                        &ctx,
                        (public).then(|| {
                            fhe_util::VariableTime::new(fhe_util::PublicData::assert_public())
                        }),
                    )?;
                    let expected = Poly::<NttShoup>::from_rns_residues(standard.clone(), &ctx)?;
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
            assert_eq!(Poly::<PowerBasis>::from_proto(&proto, &ctx)?, p);
            assert_eq!(
                Poly::<Ntt>::from_proto(&proto, &ctx).unwrap_err(),
                CrateError::PolynomialSerialization(
                    PolynomialSerializationError::RepresentationMismatch {
                        found: crate::rq::Representation::PowerBasis,
                        expected: crate::rq::Representation::Ntt,
                    }
                )
            );
            assert_eq!(
                Poly::<NttShoup>::from_proto(&proto, &ctx).unwrap_err(),
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
        assert_eq!(Poly::<Ntt>::from_proto(&proto, &ctx)?, p);

        let p = Poly::<NttShoup>::random(&ctx, &mut rng);
        let proto = Rq::from(&p);
        assert_eq!(Poly::<NttShoup>::from_proto(&proto, &ctx)?, p);

        Ok(())
    }

    #[test]
    fn try_convert_from_slice_zero() -> Result<(), Box<dyn Error>> {
        for modulus in MODULI {
            let ctx = Arc::new(Context::new(&[*modulus], 16)?);

            // Power Basis
            assert_eq!(
                Poly::<PowerBasis>::from_coefficients(&[0u64], &ctx)?,
                Poly::<PowerBasis>::zero(&ctx)
            );
            assert_eq!(
                Poly::<PowerBasis>::from_signed_coefficients(&[0i64], &ctx)?,
                Poly::<PowerBasis>::zero(&ctx)
            );
            assert_eq!(
                Poly::<PowerBasis>::from_coefficients(&[0u64; 16], &ctx)?,
                Poly::<PowerBasis>::zero(&ctx)
            );
            assert_eq!(
                Poly::<PowerBasis>::from_signed_coefficients(&[0i64; 16], &ctx)?,
                Poly::<PowerBasis>::zero(&ctx)
            );
            assert!(Poly::<PowerBasis>::from_coefficients(&[0u64; 17], &ctx).is_err());

            // Ntt
            assert!(Poly::<Ntt>::from_rns_slice(&[0u64], &ctx).is_err());
            assert!(Poly::<Ntt>::from_rns_slice(&[0u64; 16], &ctx).is_ok());
            assert!(Poly::<Ntt>::from_rns_slice(&[0u64; 17], &ctx).is_err());
        }

        let ctx = Arc::new(Context::new(MODULI, 16)?);
        assert_eq!(
            Poly::<PowerBasis>::from_coefficients(&Vec::<u64>::default(), &ctx)?,
            Poly::<PowerBasis>::zero(&ctx)
        );
        assert!(Poly::<Ntt>::from_rns_slice(&Vec::<u64>::default(), &ctx).is_err());

        Ok(())
    }

    #[test]
    fn try_convert_from_vec_zero() -> Result<(), Box<dyn Error>> {
        for modulus in MODULI {
            let ctx = Arc::new(Context::new(&[*modulus], 16)?);
            assert_eq!(
                Poly::<PowerBasis>::from_coefficients(&[], &ctx)?,
                Poly::<PowerBasis>::zero(&ctx)
            );
            assert!(Poly::<Ntt>::from_rns_slice(&[], &ctx).is_err());

            assert_eq!(
                Poly::<PowerBasis>::from_coefficients(&[0], &ctx)?,
                Poly::<PowerBasis>::zero(&ctx)
            );
            assert!(Poly::<Ntt>::from_rns_slice(&[0], &ctx).is_err());

            assert_eq!(
                Poly::<PowerBasis>::from_coefficients(&[0; 16], &ctx)?,
                Poly::<PowerBasis>::zero(&ctx)
            );
            assert_eq!(
                Poly::<Ntt>::from_rns_slice(&[0; 16], &ctx)?,
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
        let p2 = Poly::<PowerBasis>::from_biguint_coefficients(values.as_slice(), &ctx)?;
        assert_eq!(p, p2);
        Ok(())
    }

    #[test]
    fn wire_requires_matching_degree_and_canonical_coefficients() -> Result<(), Box<dyn Error>> {
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
                let pb = Poly::<PowerBasis>::from_rns_slice_with_timing(
                    &values,
                    &ctx,
                    (public).then(|| {
                        fhe_util::VariableTime::new(fhe_util::PublicData::assert_public())
                    }),
                )?;
                assert_eq!(pb.coefficients(), expected);
                assert_eq!(
                    (&pb + &Poly::<PowerBasis>::zero(&ctx)).coefficients(),
                    expected
                );
                assert_eq!(
                    Poly::<PowerBasis>::from_rns_residues_with_timing(
                        array.clone(),
                        &ctx,
                        (public).then(|| fhe_util::VariableTime::new(
                            fhe_util::PublicData::assert_public()
                        ))
                    )?
                    .coefficients(),
                    expected
                );
                assert_eq!(
                    Poly::<Ntt>::from_rns_slice_with_timing(
                        &values,
                        &ctx,
                        (public).then(|| fhe_util::VariableTime::new(
                            fhe_util::PublicData::assert_public()
                        ))
                    )?
                    .coefficients(),
                    expected
                );
                assert_eq!(
                    Poly::<Ntt>::from_rns_residues_with_timing(
                        array.clone(),
                        &ctx,
                        (public).then(|| fhe_util::VariableTime::new(
                            fhe_util::PublicData::assert_public()
                        ))
                    )?
                    .coefficients(),
                    expected
                );
                let shoup = Poly::<NttShoup>::from_rns_slice_with_timing(
                    &values,
                    &ctx,
                    (public).then(|| {
                        fhe_util::VariableTime::new(fhe_util::PublicData::assert_public())
                    }),
                )?;
                assert_eq!(shoup.coefficients(), expected);
                assert_eq!(
                    Poly::<NttShoup>::from_rns_residues_with_timing(
                        array.clone(),
                        &ctx,
                        (public).then(|| fhe_util::VariableTime::new(
                            fhe_util::PublicData::assert_public()
                        ))
                    )?,
                    shoup
                );
                assert_eq!(pb.clone().into_ntt().into_power_basis(), pb);
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod boundary_tests {
    use super::*;

    #[test]
    fn coefficient_count_never_changes_interpretation_to_raw_residues() -> Result<()> {
        let ctx = Context::new_arc(&[1153, 2017], 16)?;
        let residues = vec![1; 32];
        assert!(matches!(
            Poly::<PowerBasis>::from_coefficients(&residues, &ctx),
            Err(Error::TooManyCoefficients {
                actual: 32,
                maximum: 16
            })
        ));
        let raw = Poly::<Ntt>::from_rns_slice(&residues, &ctx)?;
        assert_eq!(
            raw.coefficients().iter().copied().collect::<Vec<_>>(),
            residues
        );
        let coefficient = Poly::<PowerBasis>::from_coefficients(&[1], &ctx)?.into_ntt();
        assert_eq!(raw, coefficient);
        assert!(Poly::<Ntt>::from_rns_slice(&[1], &ctx).is_err());
        assert!(Poly::<PowerBasis>::from_signed_coefficients(&[1; 17], &ctx).is_err());
        assert!(
            Poly::<PowerBasis>::from_biguint_coefficients(&vec![BigUint::from(1_u64); 17], &ctx)
                .is_err()
        );
        Ok(())
    }
}
