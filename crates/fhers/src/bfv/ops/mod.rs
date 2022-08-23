//! Operations over ciphertexts

#[cfg(feature = "optimized_ops")]
mod dot_product;

#[cfg(feature = "optimized_ops")]
mod mul;

#[cfg(feature = "optimized_ops")]
pub use dot_product::dot_product_scalar;

#[cfg(feature = "optimized_ops")]
pub use mul::{mul_relin, mul_relin_2};

use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use itertools::{izip, Itertools};
use math::rq::{Poly, Representation};

use super::{Ciphertext, Plaintext};

impl Add<&Ciphertext> for &Ciphertext {
	type Output = Ciphertext;

	fn add(self, rhs: &Ciphertext) -> Ciphertext {
		let mut self_clone = self.clone();
		self_clone += rhs;
		self_clone
	}
}

impl AddAssign<&Ciphertext> for Ciphertext {
	fn add_assign(&mut self, rhs: &Ciphertext) {
		assert_eq!(self.par, rhs.par);

		if self.c.is_empty() {
			*self = rhs.clone()
		} else {
			assert_eq!(self.level, rhs.level);
			assert_eq!(self.c.len(), rhs.c.len());
			izip!(&mut self.c, &rhs.c).for_each(|(c1i, c2i)| *c1i += c2i);
			self.seed = None
		}
	}
}

impl AddAssign<Ciphertext> for Ciphertext {
	fn add_assign(&mut self, rhs: Ciphertext) {
		assert_eq!(self.par, rhs.par);

		if self.c.is_empty() {
			*self = rhs
		} else {
			assert_eq!(self.level, rhs.level);
			assert_eq!(self.c.len(), rhs.c.len());
			izip!(&mut self.c, rhs.c).for_each(|(c1i, c2i)| *c1i += c2i);
			self.seed = None
		}
	}
}

impl Sub<&Ciphertext> for &Ciphertext {
	type Output = Ciphertext;

	fn sub(self, rhs: &Ciphertext) -> Ciphertext {
		let mut self_clone = self.clone();
		self_clone -= rhs;
		self_clone
	}
}

impl SubAssign<&Ciphertext> for Ciphertext {
	fn sub_assign(&mut self, rhs: &Ciphertext) {
		assert_eq!(self.par, rhs.par);

		if self.c.is_empty() {
			*self = -rhs
		} else {
			assert_eq!(self.level, rhs.level);
			assert_eq!(self.c.len(), rhs.c.len());
			izip!(&mut self.c, &rhs.c).for_each(|(c1i, c2i)| *c1i -= c2i);
			self.seed = None
		}
	}
}

impl Neg for &Ciphertext {
	type Output = Ciphertext;

	fn neg(self) -> Ciphertext {
		let c = self.c.iter().map(|c1i| -c1i).collect_vec();
		Ciphertext {
			par: self.par.clone(),
			seed: None,
			c,
			level: self.level,
		}
	}
}

impl MulAssign<&Plaintext> for Ciphertext {
	fn mul_assign(&mut self, rhs: &Plaintext) {
		assert_eq!(self.par, rhs.par);
		if !self.c.is_empty() {
			assert_eq!(self.level, rhs.level);
			self.c.iter_mut().for_each(|ci| *ci *= &rhs.poly_ntt);
		}
		self.seed = None
	}
}

impl Mul<&Plaintext> for &Ciphertext {
	type Output = Ciphertext;

	fn mul(self, rhs: &Plaintext) -> Self::Output {
		let mut self_clone = self.clone();
		self_clone *= rhs;
		self_clone
	}
}

impl Mul<&Ciphertext> for &Ciphertext {
	type Output = Ciphertext;

	fn mul(self, rhs: &Ciphertext) -> Self::Output {
		assert_eq!(self.par, rhs.par);

		if self.c.is_empty() {
			return self.clone();
		}
		assert_eq!(self.level, rhs.level);

		let mp = &self.par.mul_1_params[self.level];

		// Scale all ciphertexts
		// let mut now = std::time::SystemTime::now();
		let self_c = self
			.c
			.iter()
			.map(|ci| ci.scale(&mp.extender_self).unwrap())
			.collect_vec();
		let other_c = rhs
			.c
			.iter()
			.map(|ci| ci.scale(&mp.extender_self).unwrap())
			.collect_vec();
		// println!("Extend: {:?}", now.elapsed().unwrap());

		// Multiply
		// now = std::time::SystemTime::now();
		let mut c = vec![Poly::zero(&mp.to, Representation::Ntt); self_c.len() + other_c.len() - 1];
		for i in 0..self_c.len() {
			for j in 0..other_c.len() {
				c[i + j] += &(&self_c[i] * &other_c[j])
			}
		}
		// println!("Multiply: {:?}", now.elapsed().unwrap());

		// Scale
		// now = std::time::SystemTime::now();
		let c = c
			.iter_mut()
			.map(|ci| {
				ci.change_representation(Representation::PowerBasis);
				let mut ci = ci.scale(&mp.down_scaler).unwrap();
				ci.change_representation(Representation::Ntt);
				ci
			})
			.collect_vec();
		// println!("Scale: {:?}", now.elapsed().unwrap());

		Ciphertext {
			par: self.par.clone(),
			seed: None,
			c,
			level: rhs.level,
		}
	}
}
