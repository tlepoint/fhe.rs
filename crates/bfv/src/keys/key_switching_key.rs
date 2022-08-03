//! Key-switching keys for the BFV encryption scheme

use crate::{BfvParameters, SecretKey};
use itertools::izip;
use math::{
	rns::RnsContext,
	rq::{traits::TryConvertFrom, Poly, Representation},
};
use ndarray::ArrayView2;
use rand::{thread_rng, Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use std::rc::Rc;
use zeroize::Zeroize;

/// Key switching key for the BFV encryption scheme.
#[derive(Debug, PartialEq)]
pub struct KeySwitchingKey {
	/// The parameters of the underlying BFV encryption scheme.
	pub(crate) par: Rc<BfvParameters>,

	/// The seed that generated the polynomials c1.
	pub(crate) seed: Option<<ChaCha8Rng as SeedableRng>::Seed>,

	/// The key switching elements c0.
	pub(crate) c0: Vec<Poly>,

	/// The key switching elements c1.
	pub(crate) c1: Vec<Poly>,
}

impl KeySwitchingKey {
	/// Generate a [`KeySwitchingKey`] to this [`SecretKey`] from a polynomial `from`.
	pub fn new(sk: &SecretKey, from: &Poly) -> Result<Self, String> {
		let mut c0 = Vec::with_capacity(sk.par.ciphertext_moduli.len());
		let mut c1 = Vec::with_capacity(sk.par.ciphertext_moduli.len());

		let rns = RnsContext::new(&sk.par.ciphertext_moduli)?;

		let mut seed = <ChaCha8Rng as SeedableRng>::Seed::default();
		thread_rng().fill(&mut seed);
		let mut rng = ChaCha8Rng::from_seed(seed);

		for i in 0..sk.par.ciphertext_moduli.len() {
			let mut seed_i = <ChaCha8Rng as SeedableRng>::Seed::default();
			rng.fill(&mut seed_i);

			let mut b = Poly::small(&sk.par.ctx, Representation::PowerBasis, sk.par.variance)?;
			let mut a = Poly::random_from_seed(&sk.par.ctx, Representation::Ntt, seed_i);

			let mut a_s = &a * &sk.s;
			a_s.change_representation(Representation::PowerBasis);
			b -= &a_s;

			let gi = rns.get_garner(i).unwrap();
			let mut g_i_from = gi * from;
			b += &g_i_from;

			a_s.zeroize();
			g_i_from.zeroize();

			// It is now safe to enable variable time computations.
			unsafe { a.allow_variable_time_computations() }
			unsafe { b.allow_variable_time_computations() }

			a.change_representation(Representation::NttShoup);
			b.change_representation(Representation::NttShoup);

			c0.push(b);
			c1.push(a);
		}

		Ok(Self {
			par: sk.par.clone(),
			seed: Some(seed),
			c0,
			c1,
		})
	}

	/// Key switch a polynomial.
	pub fn key_switch(
		&self,
		p_coefficients: &ArrayView2<u64>,
		acc_0: &mut Poly,
		acc_1: &mut Poly,
	) -> Result<(), String> {
		// TODO: Check representation of input polynomials
		for (c2_i_coefficients, c0_i, c1_i) in
			izip!(p_coefficients.outer_iter(), &self.c0, &self.c1)
		{
			let mut c2_i = Poly::try_convert_from(
				c2_i_coefficients.as_slice().unwrap(),
				&self.par.ctx,
				Representation::PowerBasis,
			)?;
			unsafe { c2_i.allow_variable_time_computations() }
			c2_i.change_representation(Representation::Ntt);
			*acc_0 += &(&c2_i * c0_i);
			c2_i *= c1_i;
			*acc_1 += &c2_i;
		}
		Ok(())
	}
}

#[cfg(test)]
mod tests {
	use crate::{keys::key_switching_key::KeySwitchingKey, BfvParameters, SecretKey};
	use math::{
		rns::RnsContext,
		rq::{Poly, Representation},
	};
	use num_bigint::BigUint;
	use std::rc::Rc;

	#[test]
	fn test_constructor() -> Result<(), String> {
		for params in [
			Rc::new(BfvParameters::default_one_modulus()),
			Rc::new(BfvParameters::default_two_moduli()),
		] {
			let sk = SecretKey::random(&params);
			let p = Poly::small(&params.ctx, Representation::PowerBasis, 10)?;
			let ksk = KeySwitchingKey::new(&sk, &p);
			assert!(ksk.is_ok());
		}
		Ok(())
	}

	#[test]
	fn test_key_switch() -> Result<(), String> {
		for params in [Rc::new(BfvParameters::default_two_moduli())] {
			for _ in 0..100 {
				let sk = SecretKey::random(&params);
				let mut s = Poly::small(&params.ctx, Representation::PowerBasis, 10)?;
				let ksk = KeySwitchingKey::new(&sk, &s)?;

				let mut input = Poly::random(&params.ctx, Representation::PowerBasis);
				let mut c0 = Poly::zero(&params.ctx, Representation::Ntt);
				let mut c1 = Poly::zero(&params.ctx, Representation::Ntt);
				ksk.key_switch(&input.coefficients(), &mut c0, &mut c1)?;

				let mut c2 = &c0 + &(&c1 * &sk.s);
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
}
