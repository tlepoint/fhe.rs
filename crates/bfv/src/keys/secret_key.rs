//! Secret keys for the BFV encryption scheme

use crate::{
	ciphertext::Ciphertext,
	parameters::BfvParameters,
	plaintext::Plaintext,
	traits::{Decryptor, Encryptor},
	Error, Result,
};
use itertools::Itertools;
use math::{
	rq::{traits::TryConvertFrom, Poly, Representation},
	zq::Modulus,
};
use num_bigint::BigUint;
use rand::{thread_rng, Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use std::sync::Arc;
use util::sample_vec_cbd;
use zeroize::{Zeroize, ZeroizeOnDrop};

/// Secret key for the BFV encryption scheme.
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct SecretKey {
	pub(crate) par: Arc<BfvParameters>,
	pub(crate) s: Vec<Poly>,
	pub(crate) s_coefficients: Vec<i64>,
}

impl Zeroize for SecretKey {
	fn zeroize(&mut self) {
		self.s.zeroize();
		self.s_coefficients.zeroize();
	}
}

impl ZeroizeOnDrop for SecretKey {}

impl SecretKey {
	/// Generate a random [`SecretKey`].
	pub fn random(par: &Arc<BfvParameters>) -> Self {
		let s_coefficients = sample_vec_cbd(par.degree(), par.variance).unwrap();
		Self::new(s_coefficients, par)
	}

	/// Generate a [`SecretKey`] from its coefficients.
	pub(crate) fn new(s_coefficients: Vec<i64>, par: &Arc<BfvParameters>) -> Self {
		let mut s = Poly::try_convert_from(
			&s_coefficients as &[i64],
			&par.ctx,
			false,
			Representation::PowerBasis,
		)
		.unwrap();
		s.change_representation(Representation::NttShoup);
		let mut s2 = &s * &s;
		s2.change_representation(Representation::NttShoup);
		Self {
			par: par.clone(),
			s: vec![s, s2],
			s_coefficients,
		}
	}

	/// Measure the noise in a [`Ciphertext`].
	///
	/// # Safety
	///
	/// This operations may run in a variable time depending on the value of the
	/// noise.
	pub unsafe fn measure_noise(&mut self, ct: &Ciphertext) -> Result<usize> {
		let plaintext = self.decrypt(ct)?;
		let mut m = plaintext.encode()?;

		// Let's disable variable time computations
		let mut c = ct.c[0].clone();
		c.disallow_variable_time_computations();

		for i in 1..ct.c.len() {
			if self.s.len() < i {
				self.s
					.push(self.s.last().unwrap() * self.s.first().unwrap());
				debug_assert_eq!(self.s.len(), i)
			}
			let mut cis = ct.c[i].clone();
			cis.disallow_variable_time_computations();
			cis *= &self.s[i - 1];
			c += &cis;
			cis.zeroize();
		}
		c -= &m;
		c.change_representation(Representation::PowerBasis);

		let ciphertext_modulus = self.par.ctx.modulus();
		let mut noise = 0usize;
		for coeff in Vec::<BigUint>::from(&c) {
			noise = std::cmp::max(
				noise,
				std::cmp::min(coeff.bits(), (ciphertext_modulus - &coeff).bits()) as usize,
			)
		}

		c.zeroize();
		m.zeroize();

		Ok(noise)
	}
}

impl Encryptor for SecretKey {
	fn encrypt(&self, pt: &Plaintext) -> Result<Ciphertext> {
		assert_eq!(self.par, pt.par);

		let mut seed = <ChaCha8Rng as SeedableRng>::Seed::default();
		thread_rng().fill(&mut seed);

		let mut a = Poly::random_from_seed(&self.par.ctx, Representation::Ntt, seed);
		let mut a_s = &a * &self.s[0];

		let mut b = Poly::small(&self.par.ctx, Representation::Ntt, self.par.variance).unwrap();
		b -= &a_s;

		let mut m = pt.encode()?;
		b += &m;

		// Zeroize the temporary variables holding sensitive information.
		a_s.zeroize();
		m.zeroize();

		// It is now safe to enable variable time computations.
		unsafe {
			a.allow_variable_time_computations();
			b.allow_variable_time_computations()
		}

		Ok(Ciphertext {
			par: self.par.clone(),
			seed: Some(seed),
			c: vec![b, a],
		})
	}
}

impl Decryptor for SecretKey {
	fn decrypt(&mut self, ct: &Ciphertext) -> Result<Plaintext> {
		if self.par != ct.par {
			Err(Error::DefaultError(
				"Incompatible BFV parameters".to_string(),
			))
		} else {
			let mut c = ct.c[0].clone();
			c.disallow_variable_time_computations();

			for i in 1..ct.c.len() {
				if self.s.len() < i {
					self.s
						.push(self.s.last().unwrap() * self.s.first().unwrap());
					debug_assert_eq!(self.s.len(), i)
				}
				let mut cis = ct.c[i].clone();
				cis.disallow_variable_time_computations();
				cis *= &self.s[i - 1];
				c += &cis;
				cis.zeroize();
			}
			c.change_representation(Representation::PowerBasis);

			let mut d = self.par.scaler.scale(&c, false)?;

			// TODO: Can we handle plaintext moduli that are BigUint?
			let mut v = Vec::<u64>::from(&d)
				.iter_mut()
				.map(|vi| *vi + self.par.plaintext.modulus())
				.collect_vec();
			let mut w = v[..self.par.degree()].to_vec();
			let q = Modulus::new(self.par.ciphertext_moduli[0]).unwrap();
			q.reduce_vec(&mut w);
			self.par.plaintext.reduce_vec(&mut w);

			let mut poly = Poly::try_convert_from(
				&w as &[u64],
				&self.par.ctx,
				false,
				Representation::PowerBasis,
			)?;
			poly.change_representation(Representation::Ntt);

			let pt = Plaintext {
				par: self.par.clone(),
				value: w,
				encoding: None,
				poly_ntt: poly,
			};

			// Zeroize the temporary variables potentially holding sensitive information.
			c.zeroize();
			d.zeroize();
			v.zeroize();

			Ok(pt)
		}
	}
}

#[cfg(test)]
mod tests {
	use super::SecretKey;
	use crate::{
		parameters::BfvParameters,
		traits::{Decryptor, Encoder, Encryptor},
		Encoding, Plaintext,
	};
	use math::rq::Representation;
	use std::{error::Error, sync::Arc};

	#[test]
	fn test_keygen() {
		let params = Arc::new(BfvParameters::default(1));
		let sk = SecretKey::random(&params);
		assert_eq!(sk.par, params);

		let mut s = sk.s[0].clone();
		s.change_representation(Representation::PowerBasis);
		let coefficients = Vec::<u64>::from(&s);
		coefficients.iter().for_each(|ci| {
			// Check that this is a small polynomial
			assert!(
				*ci <= 2 * sk.par.variance as u64
					|| *ci >= (sk.par.ciphertext_moduli[0] - 2 * sk.par.variance as u64)
			)
		})
	}

	#[test]
	fn test_encrypt_decrypt() -> Result<(), Box<dyn Error>> {
		for params in [
			Arc::new(BfvParameters::default(1)),
			Arc::new(BfvParameters::default(2)),
		] {
			for _ in 0..1 {
				let mut sk = SecretKey::random(&params);

				let pt = Plaintext::try_encode(
					&params.plaintext.random_vec(params.degree()) as &[u64],
					Encoding::Poly,
					&params,
				)?;
				let ct = sk.encrypt(&pt)?;
				let pt2 = sk.decrypt(&ct);

				println!("Noise: {}", unsafe { sk.measure_noise(&ct)? });
				assert!(pt2.is_ok_and(|pt2| pt2 == &pt));
			}
		}

		Ok(())
	}
}
