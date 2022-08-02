//! Keys for the BFV encryption scheme

use crate::{
	ciphertext::Ciphertext,
	parameters::BfvParameters,
	plaintext::Plaintext,
	traits::{Decryptor, Encryptor},
};
use math::rq::{traits::TryConvertFrom, Poly, Representation};
use rand::{thread_rng, Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use std::rc::Rc;
use zeroize::{Zeroize, ZeroizeOnDrop};

/// Secret key for the BFV encryption scheme.
#[derive(Debug, PartialEq)]
pub struct SecretKey {
	par: Rc<BfvParameters>,
	s: Poly,
}

impl Zeroize for SecretKey {
	fn zeroize(&mut self) {
		self.s.zeroize();
	}
}

impl ZeroizeOnDrop for SecretKey {}

impl SecretKey {
	/// Generate a random [`SecretKey`].
	pub fn random(par: &Rc<BfvParameters>) -> Self {
		let mut s = Poly::small(par.ctx(), Representation::PowerBasis, par.variance()).unwrap();
		s.change_representation(Representation::NttShoup);
		Self {
			par: par.clone(),
			s,
		}
	}
}

impl Encryptor for SecretKey {
	type Error = String;

	fn encrypt(&self, pt: &Plaintext) -> Result<Ciphertext, Self::Error> {
		let mut seed = <ChaCha8Rng as SeedableRng>::Seed::default();
		thread_rng().fill(&mut seed);

		let mut a = Poly::random_from_seed(self.par.ctx(), Representation::Ntt, seed);

		let mut b = Poly::small(
			self.par.ctx(),
			Representation::PowerBasis,
			self.par.variance(),
		)
		.unwrap();
		b.change_representation(Representation::Ntt);
		let mut a_s = &a * &self.s;
		b -= &a_s;

		let mut m = Poly::try_convert_from(pt, self.par.ctx(), Representation::PowerBasis)?;
		m.change_representation(Representation::Ntt);
		m *= self.par.delta();
		b += &m;

		a_s.zeroize();
		m.zeroize();

		// It is now safe to enable variable time computations.
		unsafe { a.allow_variable_time_computations() }
		unsafe { b.allow_variable_time_computations() }

		Ok(Ciphertext {
			par: self.par.clone(),
			seed: Some(seed),
			c0: b,
			c1: a,
		})
	}
}

impl Decryptor for SecretKey {
	type Error = String;

	fn decrypt(&self, ct: &Ciphertext) -> Result<Plaintext, Self::Error> {
		if self.par != ct.par {
			Err("Incompatible BFV parameters".to_string())
		} else {
			// Let's disable variable time computations
			let mut c0 = ct.c0.clone();
			let mut c1 = ct.c1.clone();
			c0.disallow_variable_time_computations();
			c1.disallow_variable_time_computations();

			let mut c1_s = &c1 * &self.s;
			let mut c = &c0 + &c1_s;
			c.change_representation(Representation::PowerBasis);
			let mut d = self.par.scaler().scale(&c, false)?;
			let mut v = Vec::<u64>::from(&d);
			self.par.plaintext().reduce_vec(&mut v);
			let pt = Plaintext {
				par: self.par.clone(),
				value: v[..self.par.degree()].to_vec(),
			};

			c1_s.zeroize();
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
	use std::rc::Rc;

	#[test]
	fn test_keygen() {
		let params = Rc::new(BfvParameters::default_one_modulus());
		let sk = SecretKey::random(&params);
		assert_eq!(sk.par, params);

		let mut s = sk.s.clone();
		s.change_representation(Representation::PowerBasis);
		let coefficients = Vec::<u64>::from(&s);
		coefficients.iter().for_each(|ci| {
			// Check that this is a small polynomial
			assert!(
				*ci <= 2 * sk.par.variance() as u64
					|| *ci >= (sk.par.moduli()[0] - 2 * sk.par.variance() as u64)
			)
		})
	}

	#[test]
	fn test_encrypt_decrypt() {
		let params = Rc::new(BfvParameters::default_one_modulus());
		let sk = SecretKey::random(&params);

		let pt = Plaintext::try_encode(&[1, 2, 3, 4, 5, 6, 7, 8], Encoding::Poly, &params).unwrap();
		let ct = sk.encrypt(&pt).unwrap();
		let pt2 = sk.decrypt(&ct);

		assert!(pt2.is_ok_and(|pt2| pt2 == &pt));

		let params = Rc::new(BfvParameters::default_two_moduli());
		let sk = SecretKey::random(&params);

		let pt = Plaintext::try_encode(&[1, 2, 3, 4, 5, 6, 7, 8], Encoding::Poly, &params).unwrap();
		let ct = sk.encrypt(&pt).unwrap();
		let pt2 = sk.decrypt(&ct);

		assert!(pt2.is_ok_and(|pt2| pt2 == &pt));
	}
}
