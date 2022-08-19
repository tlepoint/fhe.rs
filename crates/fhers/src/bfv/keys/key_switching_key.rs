//! Key-switching keys for the BFV encryption scheme

use crate::bfv::{Ciphertext, BfvParametersSwitcher};
use crate::bfv::{
	proto::bfv::KeySwitchingKey as KeySwitchingKeyProto,
	traits::TryConvertFrom as BfvTryConvertFrom, BfvParameters, SecretKey,
};
use crate::{Error, Result};
use fhers_traits::{DeserializeWithContext, Serialize};
use itertools::izip;
use math::{
	rns::RnsContext,
	rq::{Poly, Representation},
};
use rand::{thread_rng, Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use std::sync::Arc;
use zeroize::Zeroize;

/// Key switching key for the BFV encryption scheme.
#[derive(Debug, PartialEq, Eq)]
pub struct KeySwitchingKey {
	/// The parameters of the underlying BFV encryption scheme.
	pub(crate) par: Arc<BfvParameters>,

	/// The seed that generated the polynomials c1.
	pub(crate) seed: <ChaCha8Rng as SeedableRng>::Seed,

	/// The key switching elements c0.
	pub(crate) c0: Vec<Poly>,

	/// The key switching elements c1.
	pub(crate) c1: Vec<Poly>,
}

impl KeySwitchingKey {
	/// Generate a [`KeySwitchingKey`] to this [`SecretKey`] from a polynomial
	/// `from`.
	pub fn new(sk: &SecretKey, from: &Poly) -> Result<Self> {
		if sk.par.ciphertext_moduli.len() == 1 {
			return Err(Error::DefaultError(
				"These parameters do not support key switching".to_string(),
			));
		}

		let mut seed = <ChaCha8Rng as SeedableRng>::Seed::default();
		thread_rng().fill(&mut seed);
		let c1 = Self::generate_c1(&sk.par, seed, sk.par.ciphertext_moduli.len());
		let c0 = Self::generate_c0(sk, from, &c1)?;

		Ok(Self {
			par: sk.par.clone(),
			seed,
			c0,
			c1,
		})
	}

	/// Generate the c1's from the seed
	fn generate_c1(
		par: &Arc<BfvParameters>,
		seed: <ChaCha8Rng as SeedableRng>::Seed,
		size: usize,
	) -> Vec<Poly> {
		let mut c1 = Vec::with_capacity(size);
		let mut rng = ChaCha8Rng::from_seed(seed);
		(0..size).for_each(|_| {
			let mut seed_i = <ChaCha8Rng as SeedableRng>::Seed::default();
			rng.fill(&mut seed_i);
			let mut a = Poly::random_from_seed(&par.ctx, Representation::NttShoup, seed_i);
			unsafe { a.allow_variable_time_computations() }
			c1.push(a);
		});
		c1
	}

	/// Generate the c0's from the c1's and the secret key
	fn generate_c0(sk: &SecretKey, from: &Poly, c1: &[Poly]) -> Result<Vec<Poly>> {
		if c1.len() != sk.par.ciphertext_moduli.len() {
			return Err(Error::DefaultError("Invalid number of c1".to_string()));
		}
		let rns = RnsContext::new(&sk.par.ciphertext_moduli)?;
		let mut c0 = Vec::with_capacity(sk.par.ciphertext_moduli.len());
		for (i, c1i) in c1.iter().enumerate().take(sk.par.ciphertext_moduli.len()) {
			let mut a = c1i.clone();
			a.disallow_variable_time_computations();
			let mut a_s = &a * &sk.s[0];
			a_s.change_representation(Representation::PowerBasis);

			let mut b = Poly::small(&sk.par.ctx, Representation::PowerBasis, sk.par.variance)?;
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
			c0.push(b);
		}
		Ok(c0)
	}

	/// Key switch a polynomial.
	pub fn key_switch(&self, p: &Poly, acc_0: &mut Poly, acc_1: &mut Poly) -> Result<()> {
		// TODO: Check representation of input polynomials
		for (c2_i_coefficients, c0_i, c1_i) in
			izip!(p.coefficients().outer_iter(), &self.c0, &self.c1)
		{
			let mut c2_i = unsafe {
				Poly::create_constant_ntt_polynomial_with_lazy_coefficients_and_variable_time(
					c2_i_coefficients.as_slice().unwrap(),
					&self.par.ctx,
				)
			};
			*acc_0 += &(&c2_i * c0_i);
			c2_i *= c1_i;
			*acc_1 += &c2_i;
		}
		Ok(())
	}

	pub fn key_switch_and_params_switch(&self, p: &Poly, acc_0: &mut Poly, acc_1: &mut Poly, switcher: &BfvParametersSwitcher) -> Result<()> {
		if self.par != switcher.from {
			return Err(Error::DefaultError("Invalid switcher".to_string()))
		}
		// let now = std::time::Instant::now();
		let mut c0 = Poly::zero(&self.par.ctx, Representation::Ntt);
		let mut c1 = Poly::zero(&self.par.ctx, Representation::Ntt);
		self.key_switch(p, &mut c0, &mut c1)?;
		// println!("key_switch_no_scale {:?}", now.elapsed());
		// let now = std::time::Instant::now();
		c0.change_representation(Representation::PowerBasis);
		c1.change_representation(Representation::PowerBasis);
		c0.mod_switch_down(acc_0.ctx())?;
		c1.mod_switch_down(acc_1.ctx())?;
		c0.change_representation(Representation::Ntt);
		c1.change_representation(Representation::Ntt);
		// println!("mod_switch_down {:?}", now.elapsed());
		*acc_0 += &c0;
		*acc_1 += &c1;
		// *acc_0 += &switcher.scaler.scale(&c0).map_err(Error::MathError)?;
		// *acc_1 += &switcher.scaler.scale(&c1).map_err(Error::MathError)?;

		Ok(())
	}
}

impl From<&KeySwitchingKey> for KeySwitchingKeyProto {
	fn from(value: &KeySwitchingKey) -> Self {
		let mut ksk = KeySwitchingKeyProto::new();
		ksk.seed = value.seed.to_vec();
		for c0 in &value.c0 {
			ksk.c0.push(c0.to_bytes())
		}
		ksk
	}
}

impl BfvTryConvertFrom<&KeySwitchingKeyProto> for KeySwitchingKey {
	fn try_convert_from(value: &KeySwitchingKeyProto, par: &Arc<BfvParameters>) -> Result<Self> {
		let seed = <ChaCha8Rng as SeedableRng>::Seed::try_from(value.seed.clone());
		if seed.is_err() {
			return Err(Error::DefaultError("Invalid seed".to_string()));
		}
		if value.c0.len() != par.ciphertext_moduli.len() {
			return Err(Error::DefaultError("Invalid number of c0".to_string()));
		}

		let seed = seed.unwrap();
		let c1 = Self::generate_c1(par, seed, par.ciphertext_moduli.len());
		let mut c0 = Vec::with_capacity(par.ciphertext_moduli.len());
		for c0i in &value.c0 {
			c0.push(Poly::from_bytes(c0i, &par.ctx)?)
		}

		Ok(Self {
			par: par.clone(),
			seed,
			c0,
			c1,
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
		rq::{Poly, Representation},
	};
	use num_bigint::BigUint;
	use std::{error::Error, sync::Arc};

	#[test]
	fn test_constructor() -> Result<(), Box<dyn Error>> {
		for params in [Arc::new(BfvParameters::default(2))] {
			let sk = SecretKey::random(&params);
			let p = Poly::small(&params.ctx, Representation::PowerBasis, 10)?;
			let ksk = KeySwitchingKey::new(&sk, &p);
			assert!(ksk.is_ok());
		}
		Ok(())
	}

	#[test]
	fn test_key_switch() -> Result<(), Box<dyn Error>> {
		for params in [Arc::new(BfvParameters::default(2))] {
			for _ in 0..100 {
				let sk = SecretKey::random(&params);
				let mut s = Poly::small(&params.ctx, Representation::PowerBasis, 10)?;
				let ksk = KeySwitchingKey::new(&sk, &s)?;

				let mut input = Poly::random(&params.ctx, Representation::PowerBasis);
				let mut c0 = Poly::zero(&params.ctx, Representation::Ntt);
				let mut c1 = Poly::zero(&params.ctx, Representation::Ntt);
				ksk.key_switch(&input, &mut c0, &mut c1)?;

				let mut c2 = &c0 + &(&c1 * &sk.s[0]);
				c2.change_representation(Representation::PowerBasis);

				input.change_representation(Representation::Ntt);
				s.change_representation(Representation::Ntt);
				let mut c3 = &input * &s;
				c3.change_representation(Representation::PowerBasis);

				let rns = RnsContext::new(&params.ciphertext_moduli)?;
				Vec::<BigUint>::from(&(&c2 - &c3)).iter().for_each(|b| {
					assert!(std::cmp::min(b.bits(), (rns.modulus() - b).bits()) <= 70)
				});
			}
		}
		Ok(())
	}

	#[test]
	fn test_proto_conversion() -> Result<(), Box<dyn Error>> {
		for params in [Arc::new(BfvParameters::default(2))] {
			let sk = SecretKey::random(&params);
			let p = Poly::small(&params.ctx, Representation::PowerBasis, 10)?;
			let ksk = KeySwitchingKey::new(&sk, &p)?;
			let ksk_proto = KeySwitchingKeyProto::from(&ksk);
			assert_eq!(ksk, KeySwitchingKey::try_convert_from(&ksk_proto, &params)?);
		}
		Ok(())
	}
}
