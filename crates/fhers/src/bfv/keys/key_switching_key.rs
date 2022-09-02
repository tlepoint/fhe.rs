//! Key-switching keys for the BFV encryption scheme

use crate::bfv::{
	proto::bfv::KeySwitchingKey as KeySwitchingKeyProto,
	traits::TryConvertFrom as BfvTryConvertFrom, BfvParameters, SecretKey,
};
use crate::{Error, Result};
use fhers_traits::{DeserializeWithContext, Serialize};
use itertools::{izip, Itertools};
use math::rq::traits::TryConvertFrom;
use math::rq::Context;
use math::{
	rns::RnsContext,
	rq::{Poly, Representation},
};
use rand::{thread_rng, Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use std::sync::Arc;
use zeroize::Zeroize;

/// Key switching key for the BFV encryption scheme.
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct KeySwitchingKey {
	/// The parameters of the underlying BFV encryption scheme.
	pub(crate) par: Arc<BfvParameters>,

	/// The seed that generated the polynomials c1.
	pub(crate) seed: Option<<ChaCha8Rng as SeedableRng>::Seed>,

	/// The key switching elements c0.
	pub(crate) c0: Box<[Poly]>,

	/// The key switching elements c1.
	pub(crate) c1: Box<[Poly]>,

	/// The context of the polynomials that will be key switched.
	pub(crate) ciphertext_level: usize,
	pub(crate) ctx_ciphertext: Arc<Context>,

	/// The context of the key switching key.
	pub(crate) ksk_level: usize,
	pub(crate) ctx_ksk: Arc<Context>,
}

impl KeySwitchingKey {
	/// Generate a [`KeySwitchingKey`] to this [`SecretKey`] from a polynomial
	/// `from`.
	pub fn new(
		sk: &SecretKey,
		from: &Poly,
		ciphertext_level: usize,
		ksk_level: usize,
	) -> Result<Self> {
		let ctx_ksk = sk.par.ctx_at_level(ksk_level)?;
		let ctx_ciphertext = sk.par.ctx_at_level(ciphertext_level)?;

		if ctx_ksk.moduli().len() == 1 {
			return Err(Error::DefaultError(
				"These parameters do not support key switching".to_string(),
			));
		}

		if from.ctx() != ctx_ksk {
			return Err(Error::DefaultError(
				"Incorrect context for polynomial from".to_string(),
			));
		}

		let mut seed = <ChaCha8Rng as SeedableRng>::Seed::default();
		thread_rng().fill(&mut seed);
		let c1 = Self::generate_c1(ctx_ksk, seed, ctx_ciphertext.moduli().len());
		let c0 = Self::generate_c0(sk, from, &c1)?;

		Ok(Self {
			par: sk.par.clone(),
			seed: Some(seed),
			c0: c0.into_boxed_slice(),
			c1: c1.into_boxed_slice(),
			ciphertext_level,
			ctx_ciphertext: ctx_ciphertext.clone(),
			ksk_level,
			ctx_ksk: ctx_ksk.clone(),
		})
	}

	/// Generate the c1's from the seed
	fn generate_c1(
		ctx: &Arc<Context>,
		seed: <ChaCha8Rng as SeedableRng>::Seed,
		size: usize,
	) -> Vec<Poly> {
		let mut c1 = Vec::with_capacity(size);
		let mut rng = ChaCha8Rng::from_seed(seed);
		(0..size).for_each(|_| {
			let mut seed_i = <ChaCha8Rng as SeedableRng>::Seed::default();
			rng.fill(&mut seed_i);
			let mut a = Poly::random_from_seed(ctx, Representation::NttShoup, seed_i);
			unsafe { a.allow_variable_time_computations() }
			c1.push(a);
		});
		c1
	}

	/// Generate the c0's from the c1's and the secret key
	fn generate_c0(sk: &SecretKey, from: &Poly, c1: &[Poly]) -> Result<Vec<Poly>> {
		if c1.is_empty() {
			return Err(Error::DefaultError("Empty number of c1's".to_string()));
		}
		if from.representation() != &Representation::PowerBasis {
			return Err(Error::DefaultError(
				"Unexpected representation for from".to_string(),
			));
		}

		let size = c1.len();

		let mut s = Poly::try_convert_from(
			&sk.s_coefficients as &[i64],
			c1[0].ctx(),
			false,
			Representation::PowerBasis,
		)?;
		s.change_representation(Representation::Ntt);

		let rns = RnsContext::new(&sk.par.moduli[..size])?;
		let c0 = c1
			.iter()
			.enumerate()
			.map(|(i, c1i)| {
				let mut a = c1i.clone();
				a.disallow_variable_time_computations();
				let mut a_s = &a * &s;
				a_s.change_representation(Representation::PowerBasis);

				let mut b =
					Poly::small(a.ctx(), Representation::PowerBasis, sk.par.variance).unwrap();
				b -= &a_s;

				let gi = rns.get_garner(i).unwrap();
				let mut g_i_from = gi * from;
				b += &g_i_from;

				a.zeroize();
				a_s.zeroize();
				g_i_from.zeroize();

				// It is now safe to enable variable time computations.
				unsafe { b.allow_variable_time_computations() }
				b.change_representation(Representation::NttShoup);
				b
			})
			.collect_vec();

		s.zeroize();
		Ok(c0)
	}

	/// Key switch a polynomial.
	pub fn key_switch(&self, p: &Poly) -> Result<(Poly, Poly)> {
		if p.ctx().as_ref() != self.ctx_ciphertext.as_ref() {
			return Err(Error::DefaultError(
				"The input polynomial does not have the correct context.".to_string(),
			));
		}
		if p.representation() != &Representation::PowerBasis {
			return Err(Error::DefaultError("Incorrect representation".to_string()));
		}

		let mut c0 = Poly::zero(&self.ctx_ksk, Representation::Ntt);
		let mut c1 = Poly::zero(&self.ctx_ksk, Representation::Ntt);
		for (c2_i_coefficients, c0_i, c1_i) in izip!(
			p.coefficients().outer_iter(),
			self.c0.iter(),
			self.c1.iter()
		) {
			let mut c2_i = unsafe {
				Poly::create_constant_ntt_polynomial_with_lazy_coefficients_and_variable_time(
					c2_i_coefficients.as_slice().unwrap(),
					&self.ctx_ksk,
				)
			};
			c0 += &(&c2_i * c0_i);
			c2_i *= c1_i;
			c1 += &c2_i;
		}
		Ok((c0, c1))
	}
}

impl From<&KeySwitchingKey> for KeySwitchingKeyProto {
	fn from(value: &KeySwitchingKey) -> Self {
		let mut ksk = KeySwitchingKeyProto::new();
		if let Some(seed) = value.seed.as_ref() {
			ksk.seed = seed.to_vec();
		} else {
			ksk.c1.reserve_exact(value.c1.len());
			for c1 in value.c1.iter() {
				ksk.c1.push(c1.to_bytes())
			}
		}
		ksk.c0.reserve_exact(value.c0.len());
		for c0 in value.c0.iter() {
			ksk.c0.push(c0.to_bytes())
		}
		ksk.ciphertext_level = value.ciphertext_level as u32;
		ksk.ksk_level = value.ksk_level as u32;
		ksk
	}
}

impl BfvTryConvertFrom<&KeySwitchingKeyProto> for KeySwitchingKey {
	fn try_convert_from(value: &KeySwitchingKeyProto, par: &Arc<BfvParameters>) -> Result<Self> {
		let ciphertext_level = value.ciphertext_level as usize;
		let ksk_level = value.ksk_level as usize;
		let ctx_ksk = par.ctx_at_level(ksk_level)?;
		let ctx_ciphertext = par.ctx_at_level(ciphertext_level)?;

		if value.c0.len() != ctx_ciphertext.moduli().len() {
			return Err(Error::DefaultError(
				"Incorrect number of values in c0".to_string(),
			));
		}

		let seed = if value.seed.is_empty() {
			if value.c1.len() != ctx_ciphertext.moduli().len() {
				return Err(Error::DefaultError(
					"Incorrect number of values in c1".to_string(),
				));
			}
			None
		} else {
			let unwrapped = <ChaCha8Rng as SeedableRng>::Seed::try_from(value.seed.clone());
			if unwrapped.is_err() {
				return Err(Error::DefaultError("Invalid seed".to_string()));
			}
			Some(unwrapped.unwrap())
		};

		let c1 = if let Some(seed) = seed {
			Self::generate_c1(ctx_ksk, seed, value.c0.len())
		} else {
			value
				.c1
				.iter()
				.map(|c1i| Poly::from_bytes(c1i, ctx_ksk).map_err(Error::MathError))
				.collect::<Result<Vec<Poly>>>()?
		};

		let c0 = value
			.c0
			.iter()
			.map(|c0i| Poly::from_bytes(c0i, ctx_ksk).map_err(Error::MathError))
			.collect::<Result<Vec<Poly>>>()?;

		Ok(Self {
			par: par.clone(),
			seed,
			c0: c0.into_boxed_slice(),
			c1: c1.into_boxed_slice(),
			ciphertext_level,
			ctx_ciphertext: ctx_ciphertext.clone(),
			ksk_level,
			ctx_ksk: ctx_ksk.clone(),
		})
	}
}

#[cfg(test)]
mod tests {
	use crate::bfv::{
		keys::key_switching_key::KeySwitchingKey,
		proto::bfv::KeySwitchingKey as KeySwitchingKeyProto, traits::TryConvertFrom, BfvParameters,
		SecretKey,
	};
	use math::{
		rns::RnsContext,
		rq::{traits::TryConvertFrom as TryConvertFromPoly, Poly, Representation},
	};
	use num_bigint::BigUint;
	use std::{error::Error, sync::Arc};

	#[test]
	fn constructor() -> Result<(), Box<dyn Error>> {
		for params in [
			Arc::new(BfvParameters::default(6, 8)),
			Arc::new(BfvParameters::default(3, 8)),
		] {
			let sk = SecretKey::random(&params);
			let ctx = params.ctx_at_level(0)?;
			let p = Poly::small(ctx, Representation::PowerBasis, 10)?;
			let ksk = KeySwitchingKey::new(&sk, &p, 0, 0);
			assert!(ksk.is_ok());
		}
		Ok(())
	}

	#[test]
	fn key_switch() -> Result<(), Box<dyn Error>> {
		for params in [Arc::new(BfvParameters::default(6, 8))] {
			for _ in 0..100 {
				let sk = SecretKey::random(&params);
				let ctx = params.ctx_at_level(0)?;
				let mut p = Poly::small(ctx, Representation::PowerBasis, 10)?;
				let ksk = KeySwitchingKey::new(&sk, &p, 0, 0)?;
				let mut s = Poly::try_convert_from(
					&sk.s_coefficients as &[i64],
					ctx,
					false,
					Representation::PowerBasis,
				)
				.map_err(crate::Error::MathError)?;
				s.change_representation(Representation::Ntt);

				let mut input = Poly::random(ctx, Representation::PowerBasis);
				let (c0, c1) = ksk.key_switch(&input)?;

				let mut c2 = &c0 + &(&c1 * &s);
				c2.change_representation(Representation::PowerBasis);

				input.change_representation(Representation::Ntt);
				p.change_representation(Representation::Ntt);
				let mut c3 = &input * &p;
				c3.change_representation(Representation::PowerBasis);

				let rns = RnsContext::new(&params.moduli)?;
				Vec::<BigUint>::from(&(&c2 - &c3)).iter().for_each(|b| {
					assert!(std::cmp::min(b.bits(), (rns.modulus() - b).bits()) <= 70)
				});
			}
		}
		Ok(())
	}

	#[test]
	fn proto_conversion() -> Result<(), Box<dyn Error>> {
		for params in [
			Arc::new(BfvParameters::default(6, 8)),
			Arc::new(BfvParameters::default(3, 8)),
		] {
			let sk = SecretKey::random(&params);
			let ctx = params.ctx_at_level(0)?;
			let p = Poly::small(ctx, Representation::PowerBasis, 10)?;
			let ksk = KeySwitchingKey::new(&sk, &p, 0, 0)?;
			let ksk_proto = KeySwitchingKeyProto::from(&ksk);
			assert_eq!(ksk, KeySwitchingKey::try_convert_from(&ksk_proto, &params)?);
		}
		Ok(())
	}
}
