use std::sync::Arc;

use math::{
	rns::ScalingFactor,
	rq::{scaler::Scaler, Context, Representation},
	zq::primes::generate_prime,
};
use num_bigint::BigUint;

use crate::{
	bfv::{keys::RelinearizationKey, BfvParameters, Ciphertext},
	Error, Result,
};

/// Multiplicator that implements a strategy for multiplying. In particular, the
/// following information can be specified:
/// - Whether `lhs` must be scaled;
/// - Whether `rhs` must be scaled;
/// - The basis at which the multiplication will occur;
/// - The scaling factor after multiplication;
/// - Whether relinearization should be used.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Multiplicator {
	par: Arc<BfvParameters>,
	pub(crate) extender_lhs: Scaler,
	pub(crate) extender_rhs: Scaler,
	pub(crate) down_scaler: Scaler,
	pub(crate) base_ctx: Arc<Context>,
	pub(crate) mul_ctx: Arc<Context>,
	rk: Option<RelinearizationKey>,
	mod_switch: bool,
	level: usize,
}

impl Multiplicator {
	/// Construct a multiplicator using custom scaling factors and extended
	/// basis.
	pub fn new(
		lhs_scaling_factor: ScalingFactor,
		rhs_scaling_factor: ScalingFactor,
		extended_basis: &[u64],
		post_mul_scaling_factor: ScalingFactor,
		par: &Arc<BfvParameters>,
	) -> Result<Self> {
		Self::new_leveled(
			lhs_scaling_factor,
			rhs_scaling_factor,
			extended_basis,
			post_mul_scaling_factor,
			0,
			par,
		)
	}

	#[cfg(feature = "leveled_bfv")]
	#[doc(cfg(feature = "leveled_bfv"))]
	/// Construct a multiplicator using custom scaling factors and extended
	/// basis at a given level.
	pub fn new_leveled(
		lhs_scaling_factor: ScalingFactor,
		rhs_scaling_factor: ScalingFactor,
		extended_basis: &[u64],
		post_mul_scaling_factor: ScalingFactor,
		level: usize,
		par: &Arc<BfvParameters>,
	) -> Result<Self> {
		let base_ctx = par.ctx_at_level(level)?;
		let mul_ctx = Arc::new(Context::new(extended_basis, par.degree())?);
		let extender_lhs = Scaler::new(base_ctx, &mul_ctx, lhs_scaling_factor)?;
		let extender_rhs = Scaler::new(base_ctx, &mul_ctx, rhs_scaling_factor)?;
		let down_scaler = Scaler::new(&mul_ctx, base_ctx, post_mul_scaling_factor)?;
		Ok(Self {
			par: par.clone(),
			extender_lhs,
			extender_rhs,
			down_scaler,
			base_ctx: base_ctx.clone(),
			mul_ctx,
			rk: None,
			mod_switch: false,
			level,
		})
	}

	/// Default multiplication strategy using relinearization.
	pub fn default(rk: &RelinearizationKey) -> Result<Self> {
		let ctx = rk.ksk.par.ctx_at_level(rk.ksk.ciphertext_level)?;

		let modulus_size = rk.ksk.par.moduli_sizes()[..ctx.moduli().len()]
			.iter()
			.sum::<usize>();
		let n_moduli = (modulus_size + 60).div_ceil(62);

		let mut extended_basis = Vec::with_capacity(ctx.moduli().len() + n_moduli);
		extended_basis.append(&mut ctx.moduli().to_vec());
		let mut upper_bound = 1 << 62;
		while extended_basis.len() != ctx.moduli().len() + n_moduli {
			upper_bound = generate_prime(62, 2 * rk.ksk.par.degree() as u64, upper_bound).unwrap();
			if !extended_basis.contains(&upper_bound) && !ctx.moduli().contains(&upper_bound) {
				extended_basis.push(upper_bound)
			}
		}

		let mut multiplicator = Multiplicator::new_leveled(
			ScalingFactor::one(),
			ScalingFactor::one(),
			&extended_basis,
			ScalingFactor::new(
				&BigUint::from(rk.ksk.par.plaintext.modulus()),
				ctx.modulus(),
			),
			rk.ksk.ciphertext_level,
			&rk.ksk.par,
		)?;

		multiplicator.enable_relinearization(rk)?;
		Ok(multiplicator)
	}

	/// Enable relinearization after multiplication.
	pub fn enable_relinearization(&mut self, rk: &RelinearizationKey) -> Result<()> {
		let rk_ctx = self.par.ctx_at_level(rk.ksk.ciphertext_level)?;
		if rk_ctx != &self.base_ctx {
			return Err(Error::DefaultError(
				"Invalid relinearization key context".to_string(),
			));
		}
		self.rk = Some(rk.clone());
		Ok(())
	}

	#[cfg(feature = "leveled_bfv")]
	#[doc(cfg(feature = "leveled_bfv"))]
	/// Enable modulus switching after multiplication (and relinearization, if
	/// applicable).
	pub fn enable_mod_switching(&mut self) -> Result<()> {
		if self.par.ctx_at_level(self.par.max_level())? == &self.base_ctx {
			Err(Error::DefaultError(
				"Cannot modulo switch as this is already the last level".to_string(),
			))
		} else {
			self.mod_switch = true;
			Ok(())
		}
	}

	/// Multiply two ciphertexts using the defined multiplication strategy.
	pub fn multiply(&self, lhs: &Ciphertext, rhs: &Ciphertext) -> Result<Ciphertext> {
		if lhs.par != self.par || rhs.par != self.par {
			return Err(Error::DefaultError(
				"Ciphertexts do not have the same parameters".to_string(),
			));
		}
		if lhs.level != self.level || rhs.level != self.level {
			return Err(Error::DefaultError(
				"Ciphertexts are not at expected level".to_string(),
			));
		}
		if lhs.c.len() != 2 || rhs.c.len() != 2 {
			return Err(Error::DefaultError(
				"Multiplication can only be performed on ciphertexts of size 2".to_string(),
			));
		}

		// Extend
		// let mut now = std::time::SystemTime::now();
		let c00 = lhs.c[0].scale(&self.extender_lhs)?;
		let c01 = lhs.c[1].scale(&self.extender_lhs)?;
		let c10 = rhs.c[0].scale(&self.extender_rhs)?;
		let c11 = rhs.c[1].scale(&self.extender_rhs)?;
		// println!("Extend: {:?}", now.elapsed().unwrap());

		// Multiply
		// now = std::time::SystemTime::now();
		let mut c0 = &c00 * &c10;
		let mut c1 = &c00 * &c11;
		c1 += &(&c01 * &c10);
		let mut c2 = &c01 * &c11;
		c0.change_representation(Representation::PowerBasis);
		c1.change_representation(Representation::PowerBasis);
		c2.change_representation(Representation::PowerBasis);
		// println!("Multiply: {:?}", now.elapsed().unwrap());

		// Scale
		// now = std::time::SystemTime::now();
		let c0 = c0.scale(&self.down_scaler)?;
		let c1 = c1.scale(&self.down_scaler)?;
		let c2 = c2.scale(&self.down_scaler)?;
		// println!("Scale: {:?}", now.elapsed().unwrap());

		let mut c = vec![c0, c1, c2];

		// Relinearize
		if let Some(rk) = self.rk.as_ref() {
			// now = std::time::SystemTime::now();
			let (mut c0r, mut c1r) = rk.relinearizes_poly(&c[2])?;
			if c0r.ctx() != c[0].ctx() {
				c0r.change_representation(Representation::PowerBasis);
				c1r.change_representation(Representation::PowerBasis);
				c0r.mod_switch_down_to(c[0].ctx())?;
				c1r.mod_switch_down_to(c[1].ctx())?;
			} else {
				c[0].change_representation(Representation::Ntt);
				c[1].change_representation(Representation::Ntt);
			}
			c[0] += &c0r;
			c[1] += &c1r;
			c.truncate(2);
			// println!("Relinearize: {:?}", now.elapsed().unwrap());
		}

		// We construct a ciphertext, but it may not have the right representation for
		// the polynomials yet.
		let mut c = Ciphertext {
			par: self.par.clone(),
			seed: None,
			c,
			level: self.level,
		};

		if self.mod_switch {
			// now = std::time::SystemTime::now();
			c.mod_switch_to_next_level();
		// println!("Modulo switch: {:?}", now.elapsed().unwrap());
		} else {
			// now = std::time::SystemTime::now();
			c.c.iter_mut()
				.for_each(|p| p.change_representation(Representation::Ntt));
			// println!("Convert to NTT: {:?}", now.elapsed().unwrap());
		}

		Ok(c)
	}
}

#[cfg(test)]
mod tests {
	use crate::bfv::{
		BfvParameters, Ciphertext, Encoding, Plaintext, RelinearizationKey, SecretKey,
	};
	use fhe_traits::{FheDecoder, FheDecrypter, FheEncoder, FheEncrypter};
	use math::{
		rns::{RnsContext, ScalingFactor},
		zq::primes::generate_prime,
	};
	use num_bigint::BigUint;
	use std::{error::Error, sync::Arc};

	use super::Multiplicator;

	#[test]
	fn mul() -> Result<(), Box<dyn Error>> {
		let par = Arc::new(BfvParameters::default(3, 8));
		for _ in 0..30 {
			// We will encode `values` in an Simd format, and check that the product is
			// computed correctly.
			let values = par.plaintext.random_vec(par.degree());
			let mut expected = values.clone();
			par.plaintext.mul_vec(&mut expected, &values);

			let sk = SecretKey::random(&par);
			let rk = RelinearizationKey::new(&sk)?;
			let pt = Plaintext::try_encode(&values as &[u64], Encoding::simd(), &par)?;
			let ct1 = sk.try_encrypt(&pt)?;
			let ct2 = sk.try_encrypt(&pt)?;

			let mut multiplicator = Multiplicator::default(&rk)?;
			let ct3 = multiplicator.multiply(&ct1, &ct2)?;
			println!("Noise: {}", unsafe { sk.measure_noise(&ct3)? });
			let pt = sk.try_decrypt(&ct3)?;
			assert_eq!(Vec::<u64>::try_decode(&pt, Encoding::simd())?, expected);

			multiplicator.enable_mod_switching()?;
			let ct3 = multiplicator.multiply(&ct1, &ct2)?;
			assert_eq!(ct3.level, 1);
			println!("Noise: {}", unsafe { sk.measure_noise(&ct3)? });
			let pt = sk.try_decrypt(&ct3)?;
			assert_eq!(Vec::<u64>::try_decode(&pt, Encoding::simd())?, expected);
		}
		Ok(())
	}

	#[test]
	fn mul_at_level() -> Result<(), Box<dyn Error>> {
		let par = Arc::new(BfvParameters::default(3, 8));
		for _ in 0..15 {
			for level in 0..2 {
				let values = par.plaintext.random_vec(par.degree());
				let mut expected = values.clone();
				par.plaintext.mul_vec(&mut expected, &values);

				let sk = SecretKey::random(&par);
				let rk = RelinearizationKey::new_leveled(&sk, level, level)?;
				let pt =
					Plaintext::try_encode(&values as &[u64], Encoding::simd_at_level(level), &par)?;
				let ct1: Ciphertext = sk.try_encrypt(&pt)?;
				let ct2: Ciphertext = sk.try_encrypt(&pt)?;
				assert_eq!(ct1.level, level);
				assert_eq!(ct2.level, level);

				let mut multiplicator = Multiplicator::default(&rk).unwrap();
				let ct3 = multiplicator.multiply(&ct1, &ct2).unwrap();
				println!("Noise: {}", unsafe { sk.measure_noise(&ct3)? });
				let pt = sk.try_decrypt(&ct3)?;
				assert_eq!(Vec::<u64>::try_decode(&pt, Encoding::simd())?, expected);

				multiplicator.enable_mod_switching()?;
				let ct3 = multiplicator.multiply(&ct1, &ct2)?;
				assert_eq!(ct3.level, level + 1);
				println!("Noise: {}", unsafe { sk.measure_noise(&ct3)? });
				let pt = sk.try_decrypt(&ct3)?;
				assert_eq!(Vec::<u64>::try_decode(&pt, Encoding::simd())?, expected);
			}
		}
		Ok(())
	}

	#[test]
	fn mul_no_relin() -> Result<(), Box<dyn Error>> {
		let par = Arc::new(BfvParameters::default(6, 8));
		for _ in 0..30 {
			// We will encode `values` in an Simd format, and check that the product is
			// computed correctly.
			let values = par.plaintext.random_vec(par.degree());
			let mut expected = values.clone();
			par.plaintext.mul_vec(&mut expected, &values);

			let sk = SecretKey::random(&par);
			let rk = RelinearizationKey::new(&sk)?;
			let pt = Plaintext::try_encode(&values as &[u64], Encoding::simd(), &par)?;
			let ct1 = sk.try_encrypt(&pt)?;
			let ct2 = sk.try_encrypt(&pt)?;

			let mut multiplicator = Multiplicator::default(&rk)?;
			// Remove the relinearization key.
			multiplicator.rk = None;
			let ct3 = multiplicator.multiply(&ct1, &ct2)?;
			println!("Noise: {}", unsafe { sk.measure_noise(&ct3)? });
			let pt = sk.try_decrypt(&ct3)?;
			assert_eq!(Vec::<u64>::try_decode(&pt, Encoding::simd())?, expected);

			multiplicator.enable_mod_switching()?;
			let ct3 = multiplicator.multiply(&ct1, &ct2)?;
			assert_eq!(ct3.level, 1);
			println!("Noise: {}", unsafe { sk.measure_noise(&ct3)? });
			let pt = sk.try_decrypt(&ct3)?;
			assert_eq!(Vec::<u64>::try_decode(&pt, Encoding::simd())?, expected);
		}
		Ok(())
	}

	#[test]
	fn different_mul_strategy() -> Result<(), Box<dyn Error>> {
		// Implement the second multiplication strategy from <https://eprint.iacr.org/2021/204>

		let par = Arc::new(BfvParameters::default(3, 8));
		let mut extended_basis = par.moduli().to_vec();
		extended_basis
			.push(generate_prime(62, 2 * par.degree() as u64, extended_basis[2]).unwrap());
		extended_basis
			.push(generate_prime(62, 2 * par.degree() as u64, extended_basis[3]).unwrap());
		extended_basis
			.push(generate_prime(62, 2 * par.degree() as u64, extended_basis[4]).unwrap());
		let rns = RnsContext::new(&extended_basis[3..])?;

		for _ in 0..30 {
			// We will encode `values` in an Simd format, and check that the product is
			// computed correctly.
			let values = par.plaintext.random_vec(par.degree());
			let mut expected = values.clone();
			par.plaintext.mul_vec(&mut expected, &values);

			let sk = SecretKey::random(&par);
			let pt = Plaintext::try_encode(&values as &[u64], Encoding::simd(), &par)?;
			let ct1 = sk.try_encrypt(&pt)?;
			let ct2 = sk.try_encrypt(&pt)?;

			let mut multiplicator = Multiplicator::new(
				ScalingFactor::one(),
				ScalingFactor::new(rns.modulus(), par.ctx[0].modulus()),
				&extended_basis,
				ScalingFactor::new(&BigUint::from(par.plaintext()), rns.modulus()),
				&par,
			)?;

			let ct3 = multiplicator.multiply(&ct1, &ct2)?;
			println!("Noise: {}", unsafe { sk.measure_noise(&ct3)? });
			let pt = sk.try_decrypt(&ct3)?;
			assert_eq!(Vec::<u64>::try_decode(&pt, Encoding::simd())?, expected);

			multiplicator.enable_mod_switching()?;
			let ct3 = multiplicator.multiply(&ct1, &ct2)?;
			assert_eq!(ct3.level, 1);
			println!("Noise: {}", unsafe { sk.measure_noise(&ct3)? });
			let pt = sk.try_decrypt(&ct3)?;
			assert_eq!(Vec::<u64>::try_decode(&pt, Encoding::simd())?, expected);
		}

		Ok(())
	}
}
