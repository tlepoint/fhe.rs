//! Key-switching keys for the BFV encryption scheme

use crate::BfvParameters;
use itertools::izip;
use math::rq::{Poly, Representation};
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use std::rc::Rc;

/// Secret key for the BFV encryption scheme.
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
	/// Key switch a polynomial.
	pub fn key_switch(&self, poly: &Poly) -> (Poly, Poly) {
		let mut c0 = Poly::zero(self.par.ctx(), Representation::Ntt);
		let mut c1 = Poly::zero(self.par.ctx(), Representation::Ntt);
		izip!(&self.c0, &self.c1).for_each(|(c0_i, c1_i)| {
			c0 += &(poly * c0_i);
			c1 += &(poly * c1_i);
		});
		(c0, c1)
	}
}

#[cfg(test)]
mod tests {
	use crate::{
		traits::{Decoder, Decryptor},
		BfvParameters, Ciphertext, Encoding, SecretKey,
	};
	use math::rq::{Poly, Representation};
	use std::rc::Rc;

	#[test]
	fn test_key_switching() {
		for params in [Rc::new(BfvParameters::default_two_moduli())] {
			let sk = SecretKey::random(&params);

			let mut p = Poly::small(params.ctx(), Representation::PowerBasis, 10).unwrap();
			let ksk = sk.key_switching_new(&p).unwrap();

			// Let's construct a ciphertext by key-switching p itself.
			// We should obtain a ciphertext with small noise. We print the noise and verifies that
			// the ciphertext decrypts to 0.
			p.change_representation(Representation::Ntt);
			let (c0, c1) = ksk.key_switch(&p);
			let ct = Ciphertext {
				par: params.clone(),
				seed: None,
				c0,
				c1,
			};

			println!("Noise: {}", unsafe { sk.measure_noise(&ct).unwrap() });
			let pt = sk.decrypt(&ct).unwrap();
			assert_eq!(
				Vec::<u64>::try_decode(&pt, Encoding::Poly).unwrap(),
				&[0; 8]
			)
		}
	}
}
