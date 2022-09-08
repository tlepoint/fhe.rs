//! Secret keys for the BFV encryption scheme

use crate::bfv::{BfvParameters, Ciphertext, Plaintext};
use crate::{Error, Result};
use fhe_math::{
	rq::{traits::TryConvertFrom, Poly, Representation},
	zq::Modulus,
};
use fhe_traits::{FheDecrypter, FheEncrypter, FheParametrized};
use fhe_util::sample_vec_cbd;
use itertools::Itertools;
use num_bigint::BigUint;
use rand::{thread_rng, Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use std::sync::Arc;
use zeroize::{Zeroize, ZeroizeOnDrop, Zeroizing};

/// Secret key for the BFV encryption scheme.
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct SecretKey {
	pub(crate) par: Arc<BfvParameters>,
	pub(crate) coeffs: Box<[i64]>,
}

impl Zeroize for SecretKey {
	fn zeroize(&mut self) {
		self.coeffs.zeroize();
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
		Self {
			par: par.clone(),
			coeffs: s_coefficients.into_boxed_slice(),
		}
	}

	/// Measure the noise in a [`Ciphertext`].
	///
	/// # Safety
	///
	/// This operations may run in a variable time depending on the value of the
	/// noise.
	pub unsafe fn measure_noise(&self, ct: &Ciphertext) -> Result<usize> {
		let plaintext = Zeroizing::new(self.try_decrypt(ct)?);
		let m = Zeroizing::new(plaintext.to_poly()?);

		// Let's create a secret key with the ciphertext context
		let mut s = Zeroizing::new(Poly::try_convert_from(
			&self.coeffs as &[i64],
			ct.c[0].ctx(),
			false,
			Representation::PowerBasis,
		)?);
		s.change_representation(Representation::Ntt);
		let mut si = s.clone();

		// Let's disable variable time computations
		let mut c = Zeroizing::new(ct.c[0].clone());
		c.disallow_variable_time_computations();

		for i in 1..ct.c.len() {
			let mut cis = Zeroizing::new(ct.c[i].clone());
			cis.disallow_variable_time_computations();
			*cis.as_mut() *= si.as_ref();
			*c.as_mut() += &cis;
			*si.as_mut() *= s.as_ref();
		}
		*c.as_mut() -= &m;
		c.change_representation(Representation::PowerBasis);

		let ciphertext_modulus = ct.c[0].ctx().modulus();
		let mut noise = 0usize;
		for coeff in Vec::<BigUint>::from(c.as_ref()) {
			noise = std::cmp::max(
				noise,
				std::cmp::min(coeff.bits(), (ciphertext_modulus - &coeff).bits()) as usize,
			)
		}

		Ok(noise)
	}

	pub(crate) fn encrypt_poly(&self, p: &Poly) -> Result<Ciphertext> {
		if p.representation() != &Representation::Ntt {
			return Err(Error::MathError(fhe_math::Error::IncorrectRepresentation(
				p.representation().clone(),
				Representation::Ntt,
			)));
		}

		let level = self.par.level_of_ctx(p.ctx())?;

		let mut seed = <ChaCha8Rng as SeedableRng>::Seed::default();
		thread_rng().fill(&mut seed);

		// Let's create a secret key with the ciphertext context
		let mut s = Zeroizing::new(Poly::try_convert_from(
			&self.coeffs as &[i64],
			p.ctx(),
			false,
			Representation::PowerBasis,
		)?);
		s.change_representation(Representation::Ntt);

		let mut a = Poly::random_from_seed(p.ctx(), Representation::Ntt, seed);
		let a_s = Zeroizing::new(&a * s.as_ref());

		let mut b = Poly::small(p.ctx(), Representation::Ntt, self.par.variance)
			.map_err(Error::MathError)?;
		b -= &a_s;
		b += p;

		// It is now safe to enable variable time computations.
		unsafe {
			a.allow_variable_time_computations();
			b.allow_variable_time_computations()
		}

		Ok(Ciphertext {
			par: self.par.clone(),
			seed: Some(seed),
			c: vec![b, a],
			level,
		})
	}
}

impl FheParametrized for SecretKey {
	type Parameters = BfvParameters;
}

impl FheEncrypter<Plaintext, Ciphertext> for SecretKey {
	type Error = Error;

	fn try_encrypt(&self, pt: &Plaintext) -> Result<Ciphertext> {
		assert_eq!(self.par, pt.par);
		let m = Zeroizing::new(pt.to_poly()?);
		self.encrypt_poly(m.as_ref())
	}
}

impl FheDecrypter<Plaintext, Ciphertext> for SecretKey {
	type Error = Error;

	fn try_decrypt(&self, ct: &Ciphertext) -> Result<Plaintext> {
		if self.par != ct.par {
			Err(Error::DefaultError(
				"Incompatible BFV parameters".to_string(),
			))
		} else {
			// Let's create a secret key with the ciphertext context
			let mut s = Zeroizing::new(Poly::try_convert_from(
				&self.coeffs as &[i64],
				ct.c[0].ctx(),
				false,
				Representation::PowerBasis,
			)?);
			s.change_representation(Representation::Ntt);
			let mut si = s.clone();

			let mut c = Zeroizing::new(ct.c[0].clone());
			c.disallow_variable_time_computations();

			for i in 1..ct.c.len() {
				let mut cis = Zeroizing::new(ct.c[i].clone());
				cis.disallow_variable_time_computations();
				*cis.as_mut() *= si.as_ref();
				*c.as_mut() += &cis;
				*si.as_mut() *= s.as_ref();
			}
			c.change_representation(Representation::PowerBasis);

			let d = Zeroizing::new(c.scale(&self.par.scalers[ct.level])?);

			// TODO: Can we handle plaintext moduli that are BigUint?
			let v = Zeroizing::new(
				Vec::<u64>::from(d.as_ref())
					.iter_mut()
					.map(|vi| *vi + self.par.plaintext.modulus())
					.collect_vec(),
			);
			let mut w = v[..self.par.degree()].to_vec();
			let q = Modulus::new(self.par.moduli[0]).map_err(Error::MathError)?;
			q.reduce_vec(&mut w);
			self.par.plaintext.reduce_vec(&mut w);

			let mut poly = Poly::try_convert_from(
				&w as &[u64],
				ct.c[0].ctx(),
				false,
				Representation::PowerBasis,
			)?;
			poly.change_representation(Representation::Ntt);

			let pt = Plaintext {
				par: self.par.clone(),
				value: w.into_boxed_slice(),
				encoding: None,
				poly_ntt: poly,
				level: ct.level,
			};

			Ok(pt)
		}
	}
}

#[cfg(test)]
mod tests {
	use super::SecretKey;
	use crate::bfv::{parameters::BfvParameters, Encoding, Plaintext};
	use fhe_traits::{FheDecrypter, FheEncoder, FheEncrypter};
	use std::{error::Error, sync::Arc};

	#[test]
	fn keygen() {
		let params = Arc::new(BfvParameters::default(1, 8));
		let sk = SecretKey::random(&params);
		assert_eq!(sk.par, params);

		sk.coeffs.iter().for_each(|ci| {
			// Check that this is a small polynomial
			assert!((*ci).abs() <= 2 * sk.par.variance as i64)
		})
	}

	#[test]
	fn encrypt_decrypt() -> Result<(), Box<dyn Error>> {
		for params in [
			Arc::new(BfvParameters::default(1, 8)),
			Arc::new(BfvParameters::default(6, 8)),
		] {
			for level in 0..params.max_level() {
				for _ in 0..20 {
					let sk = SecretKey::random(&params);

					let pt = Plaintext::try_encode(
						&params.plaintext.random_vec(params.degree()) as &[u64],
						Encoding::poly_at_level(level),
						&params,
					)?;
					let ct = sk.try_encrypt(&pt)?;
					let pt2 = sk.try_decrypt(&ct)?;

					println!("Noise: {}", unsafe { sk.measure_noise(&ct)? });
					assert_eq!(pt2, pt);
				}
			}
		}

		Ok(())
	}
}
