#![warn(missing_docs, unused_imports)]

//! Polynomials in R_q\[x\] = (ZZ_q1 x ... x ZZ_qn)\[x\] where the qi's are prime moduli in zq.

pub mod extender;
pub mod scaler;
pub mod traits;

use crate::{
	rns::RnsContext,
	zq::{ntt::NttOperator, Modulus},
};
use fhers_protos::protos::rq as proto_rq;
use itertools::{izip, Itertools};
use ndarray::{Array2, ArrayView, ArrayView2, Axis};
use num_bigint::BigUint;
use protobuf::EnumOrUnknown;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use std::{
	ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
	rc::Rc,
};
use traits::{TryConvertFrom, Unsigned};
use util::sample_vec_cbd;
use zeroize::{Zeroize, Zeroizing};

/// Struct that holds the context associated with elements in rq.
#[derive(Default, Debug, Clone, PartialEq)]
pub struct Context {
	q: Vec<Modulus>,
	rns: RnsContext,
	ops: Vec<NttOperator>,
	degree: usize,
}

impl Context {
	/// Creates a context from a list of moduli and a polynomial degree.
	///
	/// Returns an error if the moduli are not primes less than 62 bits which supports the NTT of size `degree`.
	pub fn new(moduli: &[u64], degree: usize) -> Result<Self, String> {
		if !degree.is_power_of_two() || degree < 8 {
			Err("The degree is not a power of two larger or equal to 8".to_string())
		} else {
			let mut q = Vec::with_capacity(moduli.len());
			let rns = RnsContext::new(moduli)?;
			let mut ops = Vec::with_capacity(moduli.len());
			for modulus in moduli {
				if let Some(qi) = Modulus::new(*modulus) {
					if let Some(op) = NttOperator::new(&qi, degree) {
						q.push(qi);
						ops.push(op);
					} else {
						return Err("Impossible to construct a Ntt operator".to_string());
					}
				} else {
					return Err("The modulus is invalid".to_string());
				}
			}

			Ok(Self {
				q,
				rns,
				ops,
				degree,
			})
		}
	}

	/// Returns the modulus as a BigUint
	/// TODO: To test
	pub fn modulus(&self) -> &BigUint {
		self.rns.modulus()
	}
}

/// Possible representations of the underlying polynomial.
#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub enum Representation {
	/// This is the list of coefficients ci, such that the polynomial is c0 + c1 * x + ... + c_(degree - 1) * x^(degree - 1)
	#[default]
	PowerBasis,
	/// This is the NTT representation of the PowerBasis representation.
	Ntt,
	/// This is a "Shoup" representation of the Ntt representation used for faster multiplication.
	NttShoup,
}

/// Struct that holds a polynomial for a specific context.
#[derive(Default, Debug, Clone, PartialEq)]
pub struct Poly {
	ctx: Rc<Context>,
	representation: Representation,
	allow_variable_time_computations: bool,
	coefficients: Array2<u64>,
	coefficients_shoup: Option<Array2<u64>>,
}

impl Poly {
	/// Creates a polynomial holding the constant 0.
	pub fn zero(ctx: &Rc<Context>, representation: Representation) -> Self {
		Self {
			ctx: ctx.clone(),
			representation: representation.clone(),
			allow_variable_time_computations: false,
			coefficients: Array2::zeros((ctx.q.len(), ctx.degree)),
			coefficients_shoup: if representation == Representation::NttShoup {
				Some(Array2::zeros((ctx.q.len(), ctx.degree)))
			} else {
				None
			},
		}
	}

	/// # Safety
	///
	/// Enable variable time computations when this polynomial is involved.
	/// By default, this is marked as unsafe, but is usually safe when only public data is processed.
	/// TODO: To test
	pub unsafe fn allow_variable_time_computations(&mut self) {
		self.allow_variable_time_computations = true
	}

	/// Disable variable time computations when this polynomial is involved.
	pub fn disallow_variable_time_computations(&mut self) {
		self.allow_variable_time_computations = false
	}

	/// Change the representation of the underlying polynomial.
	///
	/// Panics if the change of representation is illegal.
	pub fn change_representation(&mut self, to: Representation) {
		// If we are already in the correct representation, returns immediately.
		if self.representation == to {
			return;
		}

		// TODO: Should we use `match` instead?
		if self.representation == Representation::PowerBasis && to == Representation::Ntt {
			if self.allow_variable_time_computations {
				izip!(self.coefficients.outer_iter_mut(), &self.ctx.ops)
					.for_each(|(mut v, op)| unsafe { op.forward_vt(v.as_slice_mut().unwrap()) });
			} else {
				izip!(self.coefficients.outer_iter_mut(), &self.ctx.ops)
					.for_each(|(mut v, op)| op.forward(v.as_slice_mut().unwrap()));
			}
		} else if self.representation == Representation::Ntt && to == Representation::PowerBasis {
			if self.allow_variable_time_computations {
				izip!(self.coefficients.outer_iter_mut(), &self.ctx.ops)
					.for_each(|(mut v, op)| unsafe { op.backward_vt(v.as_slice_mut().unwrap()) });
			} else {
				izip!(self.coefficients.outer_iter_mut(), &self.ctx.ops)
					.for_each(|(mut v, op)| op.backward(v.as_slice_mut().unwrap()));
			}
		} else if self.representation == Representation::PowerBasis
			&& to == Representation::NttShoup
		{
			if self.allow_variable_time_computations {
				izip!(self.coefficients.outer_iter_mut(), &self.ctx.ops)
					.for_each(|(mut v, op)| unsafe { op.forward_vt(v.as_slice_mut().unwrap()) });
			} else {
				izip!(self.coefficients.outer_iter_mut(), &self.ctx.ops)
					.for_each(|(mut v, op)| op.forward(v.as_slice_mut().unwrap()));
			}
			self.compute_coefficients_shoup();
		} else if self.representation == Representation::Ntt && to == Representation::NttShoup {
			self.compute_coefficients_shoup();
		} else if self.representation == Representation::NttShoup && to == Representation::Ntt {
			// We are not sure whether this polynomial was sensitive or not, so for security, we zeroize the Shoup coefficients.
			self.coefficients_shoup
				.as_mut()
				.unwrap()
				.as_slice_mut()
				.unwrap()
				.zeroize();
			self.coefficients_shoup = None;
		} else if self.representation == Representation::NttShoup
			&& to == Representation::PowerBasis
		{
			if self.allow_variable_time_computations {
				izip!(self.coefficients.outer_iter_mut(), &self.ctx.ops)
					.for_each(|(mut v, op)| unsafe { op.backward_vt(v.as_slice_mut().unwrap()) });
			} else {
				izip!(self.coefficients.outer_iter_mut(), &self.ctx.ops)
					.for_each(|(mut v, op)| op.backward(v.as_slice_mut().unwrap()));
			}
			// We are not sure whether this polynomial was sensitive or not, so for security, we zeroize the Shoup coefficients.
			self.coefficients_shoup
				.as_mut()
				.unwrap()
				.as_slice_mut()
				.unwrap()
				.zeroize();
			self.coefficients_shoup = None;
		} else {
			panic!(
				"Invalid change of representation from {:?} to {:?}",
				self.representation, to
			)
		}
		self.representation = to;
	}

	/// Compute the Shoup representation of the coefficients.
	fn compute_coefficients_shoup(&mut self) {
		let mut coefficients_shoup = Array2::zeros((self.ctx.q.len(), self.ctx.degree));
		izip!(
			coefficients_shoup.outer_iter_mut(),
			self.coefficients.outer_iter(),
			&self.ctx.q
		)
		.for_each(|(mut v_shoup, v, qi)| {
			v_shoup
				.as_slice_mut()
				.unwrap()
				.copy_from_slice(&qi.shoup_vec(v.as_slice().unwrap()))
		});
		self.coefficients_shoup = Some(coefficients_shoup)
	}

	/// # Safety
	///
	/// Override the internal representation to a given representation.
	/// If the `to` representation is NttShoup, the coefficients are still computed correctly to avoid being in an unstable state.
	pub unsafe fn override_representation(&mut self, to: Representation) {
		if to == Representation::NttShoup {
			self.compute_coefficients_shoup()
		}
		self.representation = to;
	}

	/// Generate a random polynomial.
	pub fn random(ctx: &Rc<Context>, representation: Representation) -> Self {
		let mut p = Poly::zero(ctx, representation);
		izip!(p.coefficients.outer_iter_mut(), &ctx.q).for_each(|(mut v, qi)| {
			v.as_slice_mut()
				.unwrap()
				.copy_from_slice(&qi.random_vec(ctx.degree))
		});
		if p.representation == Representation::NttShoup {
			p.compute_coefficients_shoup()
		}
		p
	}

	/// Generate a random polynomial deterministically from a seed.
	pub fn random_from_seed(
		ctx: &Rc<Context>,
		representation: Representation,
		seed: <ChaCha8Rng as SeedableRng>::Seed,
	) -> Self {
		let mut rng = ChaCha8Rng::from_seed(seed);
		let mut p = Poly::zero(ctx, representation);
		izip!(p.coefficients.outer_iter_mut(), &ctx.q).for_each(|(mut v, qi)| {
			let mut seed_for_vec = <ChaCha8Rng as SeedableRng>::Seed::default();
			rng.fill(&mut seed_for_vec);
			v.as_slice_mut()
				.unwrap()
				.copy_from_slice(&qi.random_vec_from_seed(ctx.degree, seed_for_vec))
		});
		if p.representation == Representation::NttShoup {
			p.compute_coefficients_shoup()
		}
		p
	}

	/// Generate a small polynomial.
	///
	/// Returns an error if the variance does not belong to [1, ..., 16].
	/// TODO: To test
	pub fn small(
		ctx: &Rc<Context>,
		representation: Representation,
		variance: usize,
	) -> Result<Self, String> {
		if !(1..=16).contains(&variance) {
			Err("The variance should be an integer between 1 and 16".to_string())
		} else {
			let mut coeffs = sample_vec_cbd(ctx.degree, variance)?;
			let p = Poly::try_convert_from(coeffs.as_ref() as &[i64], ctx, representation)?;
			coeffs.zeroize();
			Ok(p)
		}
	}

	/// Access the polynomial coefficients in RNS representation.
	/// TODO: To test?
	pub fn coefficients(&self) -> ArrayView2<u64> {
		self.coefficients.view()
	}
}

impl Zeroize for Poly {
	fn zeroize(&mut self) {
		self.coefficients.as_slice_mut().unwrap().zeroize();
		if let Some(s) = self.coefficients_shoup.as_mut() {
			s.as_slice_mut().unwrap().zeroize();
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

impl TryConvertFrom<&proto_rq::Rq> for Poly {
	type Error = String;

	fn try_convert_from<R>(
		value: &proto_rq::Rq,
		ctx: &Rc<Context>,
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
			let v = qi.deserialize_vec(&value.coefficients[index..index + size]);
			if v == None {
				return Err("Could not deserialize the polynomial coefficients".to_string());
			}
			coefficients.append(&mut v.unwrap());
			index += size;
		}

		Poly::try_convert_from(coefficients, ctx, representation_from_proto)
	}
}

impl Unsigned for u8 {}
impl Unsigned for u16 {}
impl Unsigned for u32 {}
impl Unsigned for u64 {}

impl<T: Unsigned> TryConvertFrom<T> for Poly {
	type Error = String;

	fn try_convert_from<R>(
		value: T,
		ctx: &Rc<Context>,
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

impl TryConvertFrom<Vec<u64>> for Poly {
	type Error = String;

	fn try_convert_from<R>(
		mut v: Vec<u64>,
		ctx: &Rc<Context>,
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

impl TryConvertFrom<Array2<u64>> for Poly {
	type Error = String;

	fn try_convert_from<R>(
		a: Array2<u64>,
		ctx: &Rc<Context>,
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
		ctx: &Rc<Context>,
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
		ctx: &Rc<Context>,
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
		ctx: &Rc<Context>,
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
				// c.clone_from(&ctx.rns.project(vi));
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
		ctx: &Rc<Context>,
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
		ctx: &Rc<Context>,
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
		ctx: &Rc<Context>,
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
		ctx: &Rc<Context>,
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

impl AddAssign<&Poly> for Poly {
	fn add_assign(&mut self, p: &Poly) {
		assert_ne!(
			self.representation,
			Representation::NttShoup,
			"Cannot add to a polynomial in NttShoup representation"
		);
		assert_eq!(
			self.representation, p.representation,
			"Incompatible representations"
		);
		debug_assert_eq!(self.ctx, p.ctx, "Incompatible contexts");
		if self.allow_variable_time_computations || p.allow_variable_time_computations {
			izip!(
				self.coefficients.outer_iter_mut(),
				p.coefficients.outer_iter(),
				&self.ctx.q
			)
			.for_each(|(mut v1, v2, qi)| unsafe {
				qi.add_vec_vt(v1.as_slice_mut().unwrap(), v2.as_slice().unwrap())
			});
		} else {
			izip!(
				self.coefficients.outer_iter_mut(),
				p.coefficients.outer_iter(),
				&self.ctx.q
			)
			.for_each(|(mut v1, v2, qi)| {
				qi.add_vec(v1.as_slice_mut().unwrap(), v2.as_slice().unwrap())
			});
		}
	}
}

impl Add<&Poly> for &Poly {
	type Output = Poly;
	fn add(self, p: &Poly) -> Poly {
		let mut q = self.clone();
		q += p;
		q
	}
}

impl SubAssign<&Poly> for Poly {
	fn sub_assign(&mut self, p: &Poly) {
		assert_ne!(
			self.representation,
			Representation::NttShoup,
			"Cannot subtract from a polynomial in NttShoup representation"
		);
		assert_eq!(
			self.representation, p.representation,
			"Incompatible representations"
		);
		debug_assert_eq!(self.ctx, p.ctx, "Incompatible contexts");
		if self.allow_variable_time_computations || p.allow_variable_time_computations {
			izip!(
				self.coefficients.outer_iter_mut(),
				p.coefficients.outer_iter(),
				&self.ctx.q
			)
			.for_each(|(mut v1, v2, qi)| unsafe {
				qi.sub_vec_vt(v1.as_slice_mut().unwrap(), v2.as_slice().unwrap())
			});
		} else {
			izip!(
				self.coefficients.outer_iter_mut(),
				p.coefficients.outer_iter(),
				&self.ctx.q
			)
			.for_each(|(mut v1, v2, qi)| {
				qi.sub_vec(v1.as_slice_mut().unwrap(), v2.as_slice().unwrap())
			});
		}
	}
}

impl Sub<&Poly> for &Poly {
	type Output = Poly;
	fn sub(self, p: &Poly) -> Poly {
		let mut q = self.clone();
		q -= p;
		q
	}
}

impl MulAssign<&Poly> for Poly {
	fn mul_assign(&mut self, p: &Poly) {
		assert_ne!(
			self.representation,
			Representation::NttShoup,
			"Cannot multiply to a polynomial in NttShoup representation"
		);
		assert_eq!(
			self.representation,
			Representation::Ntt,
			"Multiplication requires an Ntt representation."
		);
		debug_assert_eq!(self.ctx, p.ctx, "Incompatible contexts");

		match p.representation {
			Representation::Ntt => {
				if self.allow_variable_time_computations || p.allow_variable_time_computations {
					unsafe {
						izip!(
							self.coefficients.outer_iter_mut(),
							p.coefficients.outer_iter(),
							&self.ctx.q
						)
						.for_each(|(mut v1, v2, qi)| {
							qi.mul_vec_vt(v1.as_slice_mut().unwrap(), v2.as_slice().unwrap())
						});
					}
				} else {
					izip!(
						self.coefficients.outer_iter_mut(),
						p.coefficients.outer_iter(),
						&self.ctx.q
					)
					.for_each(|(mut v1, v2, qi)| {
						qi.mul_vec(v1.as_slice_mut().unwrap(), v2.as_slice().unwrap())
					});
				}
			}
			Representation::NttShoup => {
				if self.allow_variable_time_computations || p.allow_variable_time_computations {
					izip!(
						self.coefficients.outer_iter_mut(),
						p.coefficients.outer_iter(),
						p.coefficients_shoup.as_ref().unwrap().outer_iter(),
						&self.ctx.q
					)
					.for_each(|(mut v1, v2, v2_shoup, qi)| unsafe {
						qi.mul_shoup_vec_vt(
							v1.as_slice_mut().unwrap(),
							v2.as_slice().unwrap(),
							v2_shoup.as_slice().unwrap(),
						)
					});
				} else {
					izip!(
						self.coefficients.outer_iter_mut(),
						p.coefficients.outer_iter(),
						p.coefficients_shoup.as_ref().unwrap().outer_iter(),
						&self.ctx.q
					)
					.for_each(|(mut v1, v2, v2_shoup, qi)| {
						qi.mul_shoup_vec(
							v1.as_slice_mut().unwrap(),
							v2.as_slice().unwrap(),
							v2_shoup.as_slice().unwrap(),
						)
					});
				}
			}
			_ => {
				panic!("Multiplication requires a multipliand in Ntt or NttShoup representation.")
			}
		}
	}
}

impl MulAssign<&BigUint> for Poly {
	fn mul_assign(&mut self, p: &BigUint) {
		let v: Vec<BigUint> = vec![p.clone()];
		let mut q = Poly::try_convert_from(
			v.as_ref() as &[BigUint],
			&self.ctx,
			self.representation.clone(),
		)
		.unwrap();
		q.change_representation(Representation::Ntt);
		if self.allow_variable_time_computations {
			unsafe {
				izip!(
					self.coefficients.outer_iter_mut(),
					q.coefficients.outer_iter(),
					&self.ctx.q
				)
				.for_each(|(mut v1, v2, qi)| {
					qi.mul_vec_vt(v1.as_slice_mut().unwrap(), v2.as_slice().unwrap())
				});
			}
		} else {
			izip!(
				self.coefficients.outer_iter_mut(),
				q.coefficients.outer_iter(),
				&self.ctx.q
			)
			.for_each(|(mut v1, v2, qi)| {
				qi.mul_vec(v1.as_slice_mut().unwrap(), v2.as_slice().unwrap())
			});
		}
	}
}

impl Mul<&Poly> for &Poly {
	type Output = Poly;
	fn mul(self, p: &Poly) -> Poly {
		match self.representation {
			Representation::NttShoup => {
				// TODO: To test, and do the same thing for add, sub, and neg
				let mut q = p.clone();
				if q.representation == Representation::NttShoup {
					q.coefficients_shoup
						.as_mut()
						.unwrap()
						.as_slice_mut()
						.unwrap()
						.zeroize();
					unsafe { q.override_representation(Representation::Ntt) }
				}
				q *= self;
				q
			}
			_ => {
				let mut q = self.clone();
				q *= p;
				q
			}
		}
	}
}

impl Mul<&BigUint> for &Poly {
	type Output = Poly;
	fn mul(self, p: &BigUint) -> Poly {
		let mut q = self.clone();
		q *= p;
		q
	}
}

impl Mul<&Poly> for &BigUint {
	type Output = Poly;
	fn mul(self, p: &Poly) -> Poly {
		p * self
	}
}

impl Neg for &Poly {
	type Output = Poly;

	fn neg(self) -> Poly {
		let mut out = self.clone();
		if self.allow_variable_time_computations {
			izip!(out.coefficients.outer_iter_mut(), &out.ctx.q)
				.for_each(|(mut v1, qi)| unsafe { qi.neg_vec_vt(v1.as_slice_mut().unwrap()) });
		} else {
			izip!(out.coefficients.outer_iter_mut(), &out.ctx.q)
				.for_each(|(mut v1, qi)| qi.neg_vec(v1.as_slice_mut().unwrap()));
		}
		out
	}
}

#[cfg(test)]
mod tests {
	use super::{Context, Poly, Representation, TryConvertFrom};
	use crate::zq::{ntt::supports_ntt, Modulus};
	use fhers_protos::protos::rq as proto_rq;
	use num_bigint::BigUint;
	use num_traits::Zero;
	use rand::{thread_rng, Rng, SeedableRng};
	use rand_chacha::ChaCha8Rng;
	use std::rc::Rc;

	// Moduli to be used in tests.
	static MODULI: &[u64; 3] = &[1153, 4611686018326724609, 4611686018309947393];

	#[test]
	fn test_context_constructor() {
		for modulus in MODULI {
			assert!(Context::new(&[*modulus], 8).is_ok());
			if supports_ntt(*modulus, 128) {
				assert!(Context::new(&[*modulus], 128).is_ok());
			} else {
				assert!(Context::new(&[*modulus], 128).is_err());
			}
		}

		assert!(Context::new(MODULI, 8).is_ok());
		assert!(Context::new(MODULI, 128).is_err());
	}

	#[test]
	fn test_poly_zero() -> Result<(), String> {
		let reference = &[
			BigUint::zero(),
			BigUint::zero(),
			BigUint::zero(),
			BigUint::zero(),
			BigUint::zero(),
			BigUint::zero(),
			BigUint::zero(),
			BigUint::zero(),
		];

		for modulus in MODULI {
			let ctx = Rc::new(Context::new(&[*modulus], 8)?);
			let p = Poly::zero(&ctx, Representation::PowerBasis);
			let q = Poly::zero(&ctx, Representation::Ntt);
			assert_ne!(p, q);
			assert_eq!(Vec::<u64>::from(&p), &[0; 8]);
			assert_eq!(Vec::<u64>::from(&q), &[0; 8]);
		}

		let ctx = Rc::new(Context::new(MODULI, 8)?);
		let p = Poly::zero(&ctx, Representation::PowerBasis);
		let q = Poly::zero(&ctx, Representation::Ntt);
		assert_ne!(p, q);
		assert_eq!(Vec::<u64>::from(&p), [0; 24]);
		assert_eq!(Vec::<u64>::from(&q), [0; 24]);
		assert_eq!(Vec::<BigUint>::from(&p), reference);
		assert_eq!(Vec::<BigUint>::from(&q), reference);

		Ok(())
	}

	#[test]
	fn test_try_convert_from_slice_zero() -> Result<(), String> {
		for modulus in MODULI {
			let ctx = Rc::new(Context::new(&[*modulus], 8)?);

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

		let ctx = Rc::new(Context::new(MODULI, 8)?);
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
			let ctx = Rc::new(Context::new(&[*modulus], 8)?);
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

		let ctx = Rc::new(Context::new(MODULI, 8)?);
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
			let ctx = Rc::new(Context::new(&[*modulus], 8)?);
			let mut p = <Poly as TryConvertFrom<u64>>::try_convert_from(
				0,
				&ctx,
				Representation::PowerBasis,
			);
			assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
			p = <Poly as TryConvertFrom<u64>>::try_convert_from(0, &ctx, Representation::Ntt);
			assert!(p.is_err());
		}

		let ctx = Rc::new(Context::new(MODULI, 8)?);
		let mut p =
			<Poly as TryConvertFrom<u64>>::try_convert_from(0, &ctx, Representation::PowerBasis);
		assert!(p.is_ok_and(|pp| pp == &Poly::zero(&ctx, Representation::PowerBasis)));
		p = <Poly as TryConvertFrom<u64>>::try_convert_from(0, &ctx, Representation::Ntt);
		assert!(p.is_err());

		Ok(())
	}

	#[test]
	fn test_add() -> Result<(), String> {
		for _ in 0..100 {
			for modulus in MODULI {
				let ctx = Rc::new(Context::new(&[*modulus], 8)?);
				let m = Modulus::new(*modulus).unwrap();

				let p = Poly::random(&ctx, Representation::PowerBasis);
				let q = Poly::random(&ctx, Representation::PowerBasis);
				let r = &p + &q;
				assert_eq!(r.representation, Representation::PowerBasis);
				let mut a = Vec::<u64>::from(&p);
				m.add_vec(&mut a, &Vec::<u64>::from(&q));
				assert_eq!(Vec::<u64>::from(&r), a);

				let p = Poly::random(&ctx, Representation::Ntt);
				let q = Poly::random(&ctx, Representation::Ntt);
				let r = &p + &q;
				assert_eq!(r.representation, Representation::Ntt);
				let mut a = Vec::<u64>::from(&p);
				m.add_vec(&mut a, &Vec::<u64>::from(&q));
				assert_eq!(Vec::<u64>::from(&r), a);
			}

			let ctx = Rc::new(Context::new(MODULI, 8)?);
			let p = Poly::random(&ctx, Representation::PowerBasis);
			let q = Poly::random(&ctx, Representation::PowerBasis);
			let mut a = Vec::<u64>::from(&p);
			let b = Vec::<u64>::from(&q);
			for i in 0..MODULI.len() {
				let m = Modulus::new(MODULI[i]).unwrap();
				m.add_vec(&mut a[i * 8..(i + 1) * 8], &b[i * 8..(i + 1) * 8])
			}
			let r = &p + &q;
			assert_eq!(r.representation, Representation::PowerBasis);
			assert_eq!(Vec::<u64>::from(&r), a);
		}
		Ok(())
	}

	#[test]
	fn test_sub() -> Result<(), String> {
		for _ in 0..100 {
			for modulus in MODULI {
				let ctx = Rc::new(Context::new(&[*modulus], 8)?);
				let m = Modulus::new(*modulus).unwrap();

				let p = Poly::random(&ctx, Representation::PowerBasis);
				let q = Poly::random(&ctx, Representation::PowerBasis);
				let r = &p - &q;
				assert_eq!(r.representation, Representation::PowerBasis);
				let mut a = Vec::<u64>::from(&p);
				m.sub_vec(&mut a, &Vec::<u64>::from(&q));
				assert_eq!(Vec::<u64>::from(&r), a);

				let p = Poly::random(&ctx, Representation::Ntt);
				let q = Poly::random(&ctx, Representation::Ntt);
				let r = &p - &q;
				assert_eq!(r.representation, Representation::Ntt);
				let mut a = Vec::<u64>::from(&p);
				m.sub_vec(&mut a, &Vec::<u64>::from(&q));
				assert_eq!(Vec::<u64>::from(&r), a);
			}

			let ctx = Rc::new(Context::new(MODULI, 8)?);
			let p = Poly::random(&ctx, Representation::PowerBasis);
			let q = Poly::random(&ctx, Representation::PowerBasis);
			let mut a = Vec::<u64>::from(&p);
			let b = Vec::<u64>::from(&q);
			for i in 0..MODULI.len() {
				let m = Modulus::new(MODULI[i]).unwrap();
				m.sub_vec(&mut a[i * 8..(i + 1) * 8], &b[i * 8..(i + 1) * 8])
			}
			let r = &p - &q;
			assert_eq!(r.representation, Representation::PowerBasis);
			assert_eq!(Vec::<u64>::from(&r), a);
		}
		Ok(())
	}

	#[test]
	fn test_mul() -> Result<(), String> {
		for _ in 0..100 {
			for modulus in MODULI {
				let ctx = Rc::new(Context::new(&[*modulus], 8)?);
				let m = Modulus::new(*modulus).unwrap();

				let p = Poly::random(&ctx, Representation::Ntt);
				let q = Poly::random(&ctx, Representation::Ntt);
				let r = &p * &q;
				assert_eq!(r.representation, Representation::Ntt);
				let mut a = Vec::<u64>::from(&p);
				m.mul_vec(&mut a, &Vec::<u64>::from(&q));
				assert_eq!(Vec::<u64>::from(&r), a);
			}

			let ctx = Rc::new(Context::new(MODULI, 8)?);
			let p = Poly::random(&ctx, Representation::Ntt);
			let q = Poly::random(&ctx, Representation::Ntt);
			let mut a = Vec::<u64>::from(&p);
			let b = Vec::<u64>::from(&q);
			for i in 0..MODULI.len() {
				let m = Modulus::new(MODULI[i]).unwrap();
				m.mul_vec(&mut a[i * 8..(i + 1) * 8], &b[i * 8..(i + 1) * 8])
			}
			let r = &p * &q;
			assert_eq!(r.representation, Representation::Ntt);
			assert_eq!(Vec::<u64>::from(&r), a);
		}
		Ok(())
	}

	#[test]
	fn test_mul_shoup() -> Result<(), String> {
		for _ in 0..100 {
			for modulus in MODULI {
				let ctx = Rc::new(Context::new(&[*modulus], 8)?);
				let m = Modulus::new(*modulus).unwrap();

				let p = Poly::random(&ctx, Representation::Ntt);
				let q = Poly::random(&ctx, Representation::NttShoup);
				let r = &p * &q;
				assert_eq!(r.representation, Representation::Ntt);
				let mut a = Vec::<u64>::from(&p);
				m.mul_vec(&mut a, &Vec::<u64>::from(&q));
				assert_eq!(Vec::<u64>::from(&r), a);
			}

			let ctx = Rc::new(Context::new(MODULI, 8)?);
			let p = Poly::random(&ctx, Representation::Ntt);
			let q = Poly::random(&ctx, Representation::NttShoup);
			let mut a = Vec::<u64>::from(&p);
			let b = Vec::<u64>::from(&q);
			for i in 0..MODULI.len() {
				let m = Modulus::new(MODULI[i]).unwrap();
				m.mul_vec(&mut a[i * 8..(i + 1) * 8], &b[i * 8..(i + 1) * 8])
			}
			let r = &p * &q;
			assert_eq!(r.representation, Representation::Ntt);
			assert_eq!(Vec::<u64>::from(&r), a);
		}
		Ok(())
	}

	#[test]
	fn test_neg() -> Result<(), String> {
		for _ in 0..100 {
			for modulus in MODULI {
				let ctx = Rc::new(Context::new(&[*modulus], 8)?);
				let m = Modulus::new(*modulus).unwrap();

				let p = Poly::random(&ctx, Representation::PowerBasis);
				let r = -&p;
				assert_eq!(r.representation, Representation::PowerBasis);
				let mut a = Vec::<u64>::from(&p);
				m.neg_vec(&mut a);
				assert_eq!(Vec::<u64>::from(&r), a);

				let p = Poly::random(&ctx, Representation::Ntt);
				let r = -&p;
				assert_eq!(r.representation, Representation::Ntt);
				let mut a = Vec::<u64>::from(&p);
				m.neg_vec(&mut a);
				assert_eq!(Vec::<u64>::from(&r), a);
			}

			let ctx = Rc::new(Context::new(MODULI, 8)?);
			let p = Poly::random(&ctx, Representation::PowerBasis);
			let mut a = Vec::<u64>::from(&p);
			for i in 0..MODULI.len() {
				let m = Modulus::new(MODULI[i]).unwrap();
				m.neg_vec(&mut a[i * 8..(i + 1) * 8])
			}
			let r = -&p;
			assert_eq!(r.representation, Representation::PowerBasis);
			assert_eq!(Vec::<u64>::from(&r), a);
		}
		Ok(())
	}

	#[test]
	fn test_random() -> Result<(), String> {
		for _ in 0..100 {
			let mut seed = <ChaCha8Rng as SeedableRng>::Seed::default();
			thread_rng().fill(&mut seed);

			for modulus in MODULI {
				let ctx = Rc::new(Context::new(&[*modulus], 8)?);
				let p = Poly::random_from_seed(&ctx, Representation::Ntt, seed.clone());
				let q = Poly::random_from_seed(&ctx, Representation::Ntt, seed.clone());
				assert_eq!(p, q);
			}

			let ctx = Rc::new(Context::new(MODULI, 8)?);
			let p = Poly::random_from_seed(&ctx, Representation::Ntt, seed.clone());
			let q = Poly::random_from_seed(&ctx, Representation::Ntt, seed.clone());
			assert_eq!(p, q);

			thread_rng().fill(&mut seed);
			let p = Poly::random_from_seed(&ctx, Representation::Ntt, seed);
			assert_ne!(p, q);

			let r = Poly::random(&ctx, Representation::Ntt);
			assert_ne!(p, r);
			assert_ne!(q, r);
		}
		Ok(())
	}

	#[test]
	fn test_proto() -> Result<(), String> {
		for modulus in MODULI {
			let ctx = Rc::new(Context::new(&[*modulus], 8)?);
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

		let ctx = Rc::new(Context::new(MODULI, 8)?);
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

		let ctx = Rc::new(Context::new(&MODULI[0..1], 8)?);
		assert_eq!(
			Poly::try_convert_from(&proto, &ctx, None)
				.expect_err("Should fail because of incorrect context"),
			"Invalid coefficients"
		);

		Ok(())
	}

	#[test]
	fn test_biguint() -> Result<(), String> {
		for _ in 0..100 {
			for modulus in MODULI {
				let ctx = Rc::new(Context::new(&[*modulus], 8)?);
				let p = Poly::random(&ctx, Representation::PowerBasis);
				let p_coeffs = Vec::<BigUint>::from(&p);
				let q =
					Poly::try_convert_from(p_coeffs.as_slice(), &ctx, Representation::PowerBasis)?;
				assert_eq!(p, q);
			}

			let ctx = Rc::new(Context::new(MODULI, 8)?);
			let p = Poly::random(&ctx, Representation::PowerBasis);
			let p_coeffs = Vec::<BigUint>::from(&p);
			assert_eq!(p_coeffs.len(), ctx.degree);
			let q = Poly::try_convert_from(p_coeffs.as_slice(), &ctx, Representation::PowerBasis)?;
			assert_eq!(p, q);
		}
		Ok(())
	}
}
