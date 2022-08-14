//! Secret keys for the BFV encryption scheme

use crate::{
	ciphertext::Ciphertext,
	parameters::BfvParameters,
	plaintext::Plaintext,
	traits::{Decryptor, Encryptor},
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
	pub(crate) s_minimized: Vec<Poly>,
}

impl Zeroize for SecretKey {
	fn zeroize(&mut self) {
		self.s.iter_mut().for_each(|si| si.zeroize());
	}
}

impl ZeroizeOnDrop for SecretKey {}

impl SecretKey {
	/// Generate a random [`SecretKey`].
	pub fn random(par: &Arc<BfvParameters>) -> Self {
		let mut s_coefficients = sample_vec_cbd(par.degree(), par.variance).unwrap();
		let mut s = Poly::try_convert_from(
			&s_coefficients as &[i64],
			&par.ctx,
			Representation::PowerBasis,
		)
		.unwrap();
		let mut s_minimized = Poly::try_convert_from(
			&s_coefficients as &[i64],
			&par.plaintext_ctx,
			Representation::PowerBasis,
		)
		.unwrap();
		s_coefficients.zeroize();
		s.change_representation(Representation::Ntt);
		s_minimized.change_representation(Representation::Ntt);
		let mut s2 = &s * &s;
		let mut s2_minimized = &s_minimized * &s_minimized;
		// TODO: Can I multiply in NttShoup representation directly?
		s.change_representation(Representation::NttShoup);
		s2.change_representation(Representation::NttShoup);
		s_minimized.change_representation(Representation::NttShoup);
		s2_minimized.change_representation(Representation::NttShoup);
		Self {
			par: par.clone(),
			s: vec![s, s2],
			s_minimized: vec![s_minimized, s2_minimized],
		}
	}

	/// Measure the noise in a [`Ciphertext`].
	///
	/// # Safety
	///
	/// This operations may run in a variable time depending on the value of the noise.
	pub unsafe fn measure_noise(&mut self, ct: &Ciphertext) -> Result<usize, String> {
		let plaintext = self.decrypt(ct)?;
		let mut m = plaintext.encode(ct.minimized)?;

		// Let's disable variable time computations
		let mut c = ct.c[0].clone();
		c.disallow_variable_time_computations();

		let s = if ct.minimized {
			&mut self.s_minimized
		} else {
			&mut self.s
		};
		for i in 1..ct.c.len() {
			if s.len() < i {
				s.push(s.last().unwrap() * s.first().unwrap());
				debug_assert_eq!(s.len(), i)
			}
			let mut cis = ct.c[i].clone();
			cis.disallow_variable_time_computations();
			cis *= &s[i - 1];
			c += &cis;
			cis.zeroize();
		}
		c -= &m;
		c.change_representation(Representation::PowerBasis);

		let ciphertext_modulus = if ct.minimized {
			self.par.plaintext_ctx.modulus()
		} else {
			self.par.ctx.modulus()
		};
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
	type Error = String;

	fn encrypt(&self, pt: &Plaintext) -> Result<Ciphertext, Self::Error> {
		let mut seed = <ChaCha8Rng as SeedableRng>::Seed::default();
		thread_rng().fill(&mut seed);

		let mut a = Poly::random_from_seed(&self.par.ctx, Representation::Ntt, seed);
		let mut a_s = &a * &self.s[0];

		let mut b = Poly::small(&self.par.ctx, Representation::Ntt, self.par.variance).unwrap();
		b -= &a_s;

		let mut m = pt.encode(false)?;
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
			minimized: false,
		})
	}
}

impl Decryptor for SecretKey {
	type Error = String;

	fn decrypt(&mut self, ct: &Ciphertext) -> Result<Plaintext, Self::Error> {
		if self.par != ct.par {
			Err("Incompatible BFV parameters".to_string())
		} else {
			let mut c = ct.c[0].clone();
			c.disallow_variable_time_computations();

			let s = if ct.minimized {
				&mut self.s_minimized
			} else {
				&mut self.s
			};
			for i in 1..ct.c.len() {
				if s.len() < i {
					s.push(s.last().unwrap() * s.first().unwrap());
					debug_assert_eq!(s.len(), i)
				}
				let mut cis = ct.c[i].clone();
				cis.disallow_variable_time_computations();
				cis *= &s[i - 1];
				c += &cis;
				cis.zeroize();
			}
			c.change_representation(Representation::PowerBasis);

			let mut d = if ct.minimized {
				self.par.scaler_minimized.scale(&c, false)?
			} else {
				self.par.scaler.scale(&c, false)?
			};

			// TODO: Can we handle plaintext moduli that are BigUint?
			let mut v = Vec::<u64>::from(&d)
				.iter_mut()
				.map(|vi| *vi + self.par.plaintext.modulus())
				.collect_vec();
			let mut w = v[..self.par.degree()].to_vec();
			let q = Modulus::new(self.par.ciphertext_moduli[0]).unwrap();
			q.reduce_vec(&mut w);
			self.par.plaintext.reduce_vec(&mut w);

			let mut poly =
				Poly::try_convert_from(&w as &[u64], &self.par.ctx, Representation::PowerBasis)?;
			poly.change_representation(Representation::Ntt);

			let pt = Plaintext {
				par: self.par.clone(),
				value: unsafe { self.par.plaintext.center_vec_vt(&w[..self.par.degree()]) },
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
	use std::sync::Arc;

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
	fn test_encrypt_decrypt() -> Result<(), String> {
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
				let mut ct = sk.encrypt(&pt)?;
				let pt2 = sk.decrypt(&ct);

				println!("Noise: {}", unsafe { sk.measure_noise(&ct)? });
				assert!(pt2.is_ok_and(|pt2| pt2 == &pt));

				ct.minimizes();
				let pt2 = sk.decrypt(&ct);

				println!("Noise: {}", unsafe { sk.measure_noise(&ct)? });
				assert!(pt2.is_ok_and(|pt2| pt2 == &pt));
			}
		}

		Ok(())
	}
}
