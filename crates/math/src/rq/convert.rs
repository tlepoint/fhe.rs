//! Implementation of conversions from and to polynomials.

use super::{traits::TryConvertFrom, traits::Unsigned, Context, Poly, Representation};
use fhers_protos::protos::rq as proto_rq;
use itertools::{izip, Itertools};
use ndarray::{Array2, ArrayView, Axis};
use num_bigint::BigUint;
use protobuf::EnumOrUnknown;
use std::sync::Arc;
use zeroize::{Zeroize, Zeroizing};

impl Unsigned for u8 {}
impl Unsigned for u16 {}
impl Unsigned for u32 {}
impl Unsigned for u64 {}

impl<T: Unsigned> TryConvertFrom<T> for Poly {
	type Error = String;

	fn try_convert_from<R>(
		value: T,
		ctx: &Arc<Context>,
		representation: R,
	) -> Result<Self, Self::Error>
	where
		R: Into<Option<Representation>>,
	{
		let repr = representation.into();
		match repr {
			Some(Representation::PowerBasis) => Poly::try_convert_from(&[value], ctx, repr),
			_ => Err(
				"Converting from constant values is only possible in PowerBasis representation"
					.to_string(),
			),
		}
	}
}

impl From<&Poly> for proto_rq::Rq {
	fn from(p: &Poly) -> Self {
		let mut proto = proto_rq::Rq::new();
		match p.representation {
			Representation::PowerBasis => {
				proto.representation = EnumOrUnknown::new(proto_rq::rq::Representation::POWERBASIS);
			}
			Representation::Ntt => {
				proto.representation = EnumOrUnknown::new(proto_rq::rq::Representation::NTT);
			}
			Representation::NttShoup => {
				proto.representation = EnumOrUnknown::new(proto_rq::rq::Representation::NTTSHOUP);
			}
		}
		let mut serialization_length = 0;
		izip!(&p.ctx.q)
			.for_each(|qi| serialization_length += qi.serialization_length(p.ctx.degree));
		let mut serialization = Vec::with_capacity(serialization_length);

		izip!(p.coefficients.outer_iter(), &p.ctx.q)
			.for_each(|(v, qi)| serialization.append(&mut qi.serialize_vec(v.as_slice().unwrap())));
		proto.coefficients = serialization;
		proto.degree = p.ctx.degree as u32;
		proto
	}
}

impl TryConvertFrom<Vec<u64>> for Poly {
	type Error = String;

	fn try_convert_from<R>(
		mut v: Vec<u64>,
		ctx: &Arc<Context>,
		representation: R,
	) -> Result<Self, Self::Error>
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
						allow_variable_time_computations: false,
						coefficients,
						coefficients_shoup: None,
					})
				} else {
					Err("In Ntt representation, all coefficients must be specified".to_string())
				}
			}
			Some(Representation::NttShoup) => {
				if let Ok(coefficients) = Array2::from_shape_vec((ctx.q.len(), ctx.degree), v) {
					let mut p = Self {
						ctx: ctx.clone(),
						representation: repr.unwrap(),
						allow_variable_time_computations: false,
						coefficients,
						coefficients_shoup: None,
					};
					p.compute_coefficients_shoup();
					Ok(p)
				} else {
					Err(
						"In NttShoup representation, all coefficients must be specified"
							.to_string(),
					)
				}
			}
			Some(Representation::PowerBasis) => {
				if v.len() == ctx.q.len() * ctx.degree {
					let coefficients =
						Array2::from_shape_vec((ctx.q.len(), ctx.degree), v).unwrap();
					Ok(Self {
						ctx: ctx.clone(),
						representation: repr.unwrap(),
						allow_variable_time_computations: false,
						coefficients,
						coefficients_shoup: None,
					})
				} else if v.len() <= ctx.degree {
					let mut out = Self::zero(ctx, repr.unwrap());
					izip!(out.coefficients.outer_iter_mut(), &ctx.q).for_each(|(mut w, qi)| {
						let wi = w.as_slice_mut().unwrap();
						wi[..v.len()].copy_from_slice(&v);
						qi.reduce_vec(wi);
					});
					v.zeroize();
					Ok(out)
				} else {
					Err("In PowerBasis representation, either all coefficients must be specified, or only coefficients up to the degree".to_string())
				}
			}
			None => Err(
				"When converting from a vector, the representation needs to be specified"
					.to_string(),
			),
		}
	}
}

impl TryConvertFrom<&proto_rq::Rq> for Poly {
	type Error = String;

	fn try_convert_from<R>(
		value: &proto_rq::Rq,
		ctx: &Arc<Context>,
		representation: R,
	) -> Result<Self, Self::Error>
	where
		R: Into<Option<Representation>>,
	{
		let repr = value.representation.enum_value_or_default();
		let representation_from_proto = match repr {
			proto_rq::rq::Representation::POWERBASIS => Representation::PowerBasis,
			proto_rq::rq::Representation::NTT => Representation::Ntt,
			proto_rq::rq::Representation::NTTSHOUP => Representation::NttShoup,
			_ => return Err("Unknown representation".to_string()),
		};

		if (representation.into() as Option<Representation>)
			.is_some_and(|r| *r != representation_from_proto)
		{
			return Err("The representation asked for does not match the representation in the serialization".to_string());
		}

		let degree = value.degree as usize;
		if degree % 8 != 0 || degree < 8 {
			return Err("Invalid degree".to_string());
		}

		let mut expected_nbytes = 0;
		izip!(&ctx.q).for_each(|qi| expected_nbytes += qi.serialization_length(degree));
		if value.coefficients.len() != expected_nbytes {
			return Err("Invalid coefficients".to_string());
		}

		let mut coefficients = Vec::with_capacity(ctx.q.len() * ctx.degree);
		let mut index = 0;
		for i in 0..ctx.q.len() {
			let qi = &ctx.q[i];
			let size = qi.serialization_length(degree);
			let mut v = qi.deserialize_vec(&value.coefficients[index..index + size]);
			coefficients.append(&mut v);
			index += size;
		}

		Poly::try_convert_from(coefficients, ctx, representation_from_proto)
	}
}

impl TryConvertFrom<Array2<u64>> for Poly {
	type Error = String;

	fn try_convert_from<R>(
		a: Array2<u64>,
		ctx: &Arc<Context>,
		representation: R,
	) -> Result<Self, Self::Error>
	where
		R: Into<Option<Representation>>,
	{
		if a.shape() != [ctx.q.len(), ctx.degree] {
			Err("The array of coefficient does not have the correct shape".to_string())
		} else if let Some(repr) = representation.into() {
			let mut p = Self {
				ctx: ctx.clone(),
				representation: repr,
				allow_variable_time_computations: false,
				coefficients: a,
				coefficients_shoup: None,
			};
			if p.representation == Representation::NttShoup {
				p.compute_coefficients_shoup()
			}
			Ok(p)
		} else {
			Err("When converting from a 2-dimensional array, the representation needs to be specified".to_string())
		}
	}
}

impl<T: Unsigned> TryConvertFrom<&[T]> for Poly {
	type Error = String;

	fn try_convert_from<R>(
		v: &[T],
		ctx: &Arc<Context>,
		representation: R,
	) -> Result<Self, Self::Error>
	where
		R: Into<Option<Representation>>,
	{
		let v_clone: Vec<u64> = v.iter().map(|vi| (*vi).into()).collect_vec();
		Poly::try_convert_from(v_clone, ctx, representation)
	}
}

impl TryConvertFrom<&[i64]> for Poly {
	type Error = String;

	fn try_convert_from<R>(
		v: &[i64],
		ctx: &Arc<Context>,
		representation: R,
	) -> Result<Self, Self::Error>
	where
		R: Into<Option<Representation>>,
	{
		if representation.into() != Some(Representation::PowerBasis) {
			Err(
				"Converting signed integer require to import in PowerBasis representation"
					.to_string(),
			)
		} else if v.len() <= ctx.degree {
			let mut out = Self::zero(ctx, Representation::PowerBasis);
			izip!(out.coefficients.outer_iter_mut(), &ctx.q).for_each(|(mut w, qi)| {
				let wi = w.as_slice_mut().unwrap();
				wi[..v.len()].copy_from_slice(Zeroizing::new(qi.reduce_vec_i64(v)).as_ref());
			});
			Ok(out)
		} else {
			Err("In PowerBasis representation with signed integers, only `degree` coefficients can be specified".to_string())
		}
	}
}

impl TryConvertFrom<&[BigUint]> for Poly {
	type Error = String;

	fn try_convert_from<R>(
		v: &[BigUint],
		ctx: &Arc<Context>,
		representation: R,
	) -> Result<Self, Self::Error>
	where
		R: Into<Option<Representation>>,
	{
		let repr = representation.into();

		if v.len() > ctx.degree {
			Err(
				"The slice contains too many big integers compared to the polynomial degree"
					.to_string(),
			)
		} else if repr.is_some() {
			let mut coefficients = Array2::zeros((ctx.q.len(), ctx.degree));

			izip!(coefficients.axis_iter_mut(Axis(1)), v).for_each(|(mut c, vi)| {
				c.assign(&ArrayView::from(&ctx.rns.project(vi)));
			});

			let mut p = Self {
				ctx: ctx.clone(),
				representation: repr.unwrap(),
				allow_variable_time_computations: false,
				coefficients,
				coefficients_shoup: None,
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
			Err(
				"When converting from a vector, the representation needs to be specified"
					.to_string(),
			)
		}
	}
}

impl<T: Unsigned> TryConvertFrom<&Vec<T>> for Poly {
	type Error = String;

	fn try_convert_from<R>(
		v: &Vec<T>,
		ctx: &Arc<Context>,
		representation: R,
	) -> Result<Self, Self::Error>
	where
		R: Into<Option<Representation>>,
	{
		let v_clone: Vec<u64> = v.iter().map(|vi| (*vi).into()).collect_vec();
		Poly::try_convert_from(v_clone, ctx, representation)
	}
}

impl<T: Unsigned, const N: usize> TryConvertFrom<&[T; N]> for Poly {
	type Error = String;

	fn try_convert_from<R>(
		v: &[T; N],
		ctx: &Arc<Context>,
		representation: R,
	) -> Result<Self, Self::Error>
	where
		R: Into<Option<Representation>>,
	{
		Poly::try_convert_from(v.as_ref(), ctx, representation)
	}
}

impl<const N: usize> TryConvertFrom<&[BigUint; N]> for Poly {
	type Error = String;

	fn try_convert_from<R>(
		v: &[BigUint; N],
		ctx: &Arc<Context>,
		representation: R,
	) -> Result<Self, Self::Error>
	where
		R: Into<Option<Representation>>,
	{
		Poly::try_convert_from(v.as_ref(), ctx, representation)
	}
}

impl<const N: usize> TryConvertFrom<&[i64; N]> for Poly {
	type Error = String;

	fn try_convert_from<R>(
		v: &[i64; N],
		ctx: &Arc<Context>,
		representation: R,
	) -> Result<Self, Self::Error>
	where
		R: Into<Option<Representation>>,
	{
		Poly::try_convert_from(v.as_ref(), ctx, representation)
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
			.map(|c| p.ctx.rns.lift(&c))
			.collect_vec()
	}
}

#[cfg(test)]
mod tests {
	use crate::rq::{traits::TryConvertFrom, Context, Poly, Representation};
	use fhers_protos::protos::rq as proto_rq;
	use num_bigint::BigUint;
	use std::sync::Arc;

	static MODULI: &[u64; 3] = &[1153, 4611686018326724609, 4611686018309947393];

	#[test]
	fn test_proto() -> Result<(), String> {
		for modulus in MODULI {
			let ctx = Arc::new(Context::new(&[*modulus], 8)?);
			let p = Poly::random(&ctx, Representation::PowerBasis);
			let proto = proto_rq::Rq::from(&p);
			assert!(Poly::try_convert_from(&proto, &ctx, None).is_ok_and(|q| q == &p));
			assert!(Poly::try_convert_from(&proto, &ctx, None).is_ok_and(|q| q == &p));
			assert_eq!(
				Poly::try_convert_from(&proto, &ctx, Representation::Ntt)
					.expect_err("Should fail because of mismatched representations"),
				"The representation asked for does not match the representation in the serialization"
			);
			assert_eq!(
				Poly::try_convert_from(&proto, &ctx, Representation::NttShoup)
					.expect_err("Should fail because of mismatched representations"),
				"The representation asked for does not match the representation in the serialization"
			);
		}

		let ctx = Arc::new(Context::new(MODULI, 8)?);
		let p = Poly::random(&ctx, Representation::PowerBasis);
		let proto = proto_rq::Rq::from(&p);
		assert!(Poly::try_convert_from(&proto, &ctx, None).is_ok_and(|q| q == &p));
		assert!(Poly::try_convert_from(&proto, &ctx, None).is_ok_and(|q| q == &p));
		assert_eq!(
			Poly::try_convert_from(&proto, &ctx, Representation::Ntt)
				.expect_err("Should fail because of mismatched representations"),
			"The representation asked for does not match the representation in the serialization"
		);
		assert_eq!(
			Poly::try_convert_from(&proto, &ctx, Representation::NttShoup)
				.expect_err("Should fail because of mismatched representations"),
			"The representation asked for does not match the representation in the serialization"
		);

		let ctx = Arc::new(Context::new(&MODULI[0..1], 8)?);
		assert_eq!(
			Poly::try_convert_from(&proto, &ctx, None)
				.expect_err("Should fail because of incorrect context"),
			"Invalid coefficients"
		);

		Ok(())
	}

	#[test]
	fn test_try_convert_from_slice_zero() -> Result<(), String> {
		for modulus in MODULI {
			let ctx = Arc::new(Context::new(&[*modulus], 8)?);

			// Power Basis
			let p = Poly::try_convert_from(&[0u64], &ctx, Representation::PowerBasis);
			assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
			let p = Poly::try_convert_from(&[0i64], &ctx, Representation::PowerBasis);
			assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
			let p = Poly::try_convert_from(&[0u64; 8], &ctx, Representation::PowerBasis);
			assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
			let p = Poly::try_convert_from(&[0i64; 8], &ctx, Representation::PowerBasis);
			assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
			let p = Poly::try_convert_from(
				&[0u64; 9], // One too many
				&ctx,
				Representation::PowerBasis,
			);
			assert!(p.is_err());

			// Ntt
			let p = Poly::try_convert_from(&[0u64], &ctx, Representation::Ntt);
			assert!(p.is_err());
			let p = Poly::try_convert_from(&[0i64], &ctx, Representation::Ntt);
			assert!(p.is_err());
			let p = Poly::try_convert_from(&[0u64; 8], &ctx, Representation::Ntt);
			assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::Ntt)));
			let p = Poly::try_convert_from(&[0i64; 8], &ctx, Representation::Ntt);
			assert!(p.is_err());
			let p = Poly::try_convert_from(
				&[0u64; 9], // One too many
				&ctx,
				Representation::Ntt,
			);
			assert!(p.is_err());
		}

		let ctx = Arc::new(Context::new(MODULI, 8)?);
		let mut p = Poly::try_convert_from(Vec::<u64>::default(), &ctx, Representation::PowerBasis);
		assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
		p = Poly::try_convert_from(Vec::<u64>::default(), &ctx, Representation::Ntt);
		assert!(p.is_err());

		p = Poly::try_convert_from(&[0u64], &ctx, Representation::PowerBasis);
		assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
		p = Poly::try_convert_from(&[0u64], &ctx, Representation::Ntt);
		assert!(p.is_err());

		p = Poly::try_convert_from(&[0u64; 8], &ctx, Representation::PowerBasis);
		assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
		p = Poly::try_convert_from(&[0u64; 8], &ctx, Representation::Ntt);
		assert!(p.is_err());

		p = Poly::try_convert_from(&[0u64; 9], &ctx, Representation::PowerBasis);
		assert!(p.is_err());
		p = Poly::try_convert_from(&[0u64; 9], &ctx, Representation::Ntt);
		assert!(p.is_err());

		p = Poly::try_convert_from(&[0u64; 24], &ctx, Representation::PowerBasis);
		assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
		p = Poly::try_convert_from(&[0u64; 24], &ctx, Representation::Ntt);
		assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::Ntt)));

		Ok(())
	}

	#[test]
	fn test_try_convert_from_vec_zero() -> Result<(), String> {
		for modulus in MODULI {
			let ctx = Arc::new(Context::new(&[*modulus], 8)?);
			let mut p = Poly::try_convert_from(vec![], &ctx, Representation::PowerBasis);
			assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
			p = Poly::try_convert_from(vec![], &ctx, Representation::Ntt);
			assert!(p.is_err());

			p = Poly::try_convert_from(vec![0], &ctx, Representation::PowerBasis);
			assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
			p = Poly::try_convert_from(vec![0], &ctx, Representation::Ntt);
			assert!(p.is_err());

			p = Poly::try_convert_from(vec![0; 8], &ctx, Representation::PowerBasis);
			assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
			p = Poly::try_convert_from(vec![0; 8], &ctx, Representation::Ntt);
			assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::Ntt)));

			p = Poly::try_convert_from(vec![0; 9], &ctx, Representation::PowerBasis);
			assert!(p.is_err());
			p = Poly::try_convert_from(vec![0; 9], &ctx, Representation::Ntt);
			assert!(p.is_err());
		}

		let ctx = Arc::new(Context::new(MODULI, 8)?);
		let mut p = Poly::try_convert_from(vec![], &ctx, Representation::PowerBasis);
		assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
		p = Poly::try_convert_from(vec![], &ctx, Representation::Ntt);
		assert!(p.is_err());

		p = Poly::try_convert_from(vec![0], &ctx, Representation::PowerBasis);
		assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
		p = Poly::try_convert_from(vec![0], &ctx, Representation::Ntt);
		assert!(p.is_err());

		p = Poly::try_convert_from(vec![0; 8], &ctx, Representation::PowerBasis);
		assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
		p = Poly::try_convert_from(vec![0; 8], &ctx, Representation::Ntt);
		assert!(p.is_err());

		p = Poly::try_convert_from(vec![0; 9], &ctx, Representation::PowerBasis);
		assert!(p.is_err());
		p = Poly::try_convert_from(vec![0; 9], &ctx, Representation::Ntt);
		assert!(p.is_err());

		p = Poly::try_convert_from(vec![0; 24], &ctx, Representation::PowerBasis);
		assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
		p = Poly::try_convert_from(vec![0; 24], &ctx, Representation::Ntt);
		assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::Ntt)));

		Ok(())
	}

	#[test]
	fn test_try_convert_from_u64_zero() -> Result<(), String> {
		for modulus in MODULI {
			let ctx = Arc::new(Context::new(&[*modulus], 8)?);
			let mut p = <Poly as TryConvertFrom<u64>>::try_convert_from(
				0,
				&ctx,
				Representation::PowerBasis,
			);
			assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
			p = <Poly as TryConvertFrom<u64>>::try_convert_from(0, &ctx, Representation::Ntt);
			assert!(p.is_err());
		}

		let ctx = Arc::new(Context::new(MODULI, 8)?);
		let mut p =
			<Poly as TryConvertFrom<u64>>::try_convert_from(0, &ctx, Representation::PowerBasis);
		assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
		p = <Poly as TryConvertFrom<u64>>::try_convert_from(0, &ctx, Representation::Ntt);
		assert!(p.is_err());

		Ok(())
	}

	#[test]
	fn test_biguint() -> Result<(), String> {
		for _ in 0..100 {
			for modulus in MODULI {
				let ctx = Arc::new(Context::new(&[*modulus], 8)?);
				let p = Poly::random(&ctx, Representation::PowerBasis);
				let p_coeffs = Vec::<BigUint>::from(&p);
				let q =
					Poly::try_convert_from(p_coeffs.as_slice(), &ctx, Representation::PowerBasis)?;
				assert_eq!(p, q);
			}

			let ctx = Arc::new(Context::new(MODULI, 8)?);
			let p = Poly::random(&ctx, Representation::PowerBasis);
			let p_coeffs = Vec::<BigUint>::from(&p);
			assert_eq!(p_coeffs.len(), ctx.degree);
			let q = Poly::try_convert_from(p_coeffs.as_slice(), &ctx, Representation::PowerBasis)?;
			assert_eq!(p, q);
		}
		Ok(())
	}
}
