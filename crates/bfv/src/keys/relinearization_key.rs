//! Relinearization keys for the BFV encryption scheme

use super::key_switching_key::KeySwitchingKey;
use crate::SecretKey;
use math::rq::{Poly, Representation};
use ndarray::ArrayView2;
use zeroize::Zeroize;

/// Relinearization key for the BFV encryption scheme.
/// A relinearization key is a special type of key switching key,
/// which switch from `s^2` to `s` where `s` is the secret key.
#[derive(Debug, PartialEq)]
pub struct RelinearizationKey {
	ksk: KeySwitchingKey,
}

impl RelinearizationKey {
	/// Generate a [`RelinearizationKey`] from a [`SecretKey`].
	pub fn new(sk: &SecretKey) -> Result<Self, String> {
		let mut s_squared = &sk.s * &sk.s;
		s_squared.change_representation(Representation::PowerBasis);
		let ksk = KeySwitchingKey::new(sk, &s_squared)?;
		s_squared.zeroize();
		Ok(Self { ksk })
	}

	/// Relinearize an "extended" ciphertext (c0, c1, c2) into a [`Ciphertext`]
	pub fn relinearize(
		&self,
		c0: &mut Poly,
		c1: &mut Poly,
		c2_coefficients: &ArrayView2<u64>,
	) -> Result<(), String> {
		self.ksk.key_switch(c2_coefficients, c0, c1)
	}
}

#[cfg(test)]
mod tests {
	use math::rq::{Poly, Representation};
	use std::rc::Rc;

	use super::RelinearizationKey;
	use crate::{
		traits::{Decoder, Decryptor},
		BfvParameters, Ciphertext, Encoding, SecretKey,
	};

	#[test]
	fn test_relinearization() -> Result<(), String> {
		for params in [Rc::new(BfvParameters::default_two_moduli())] {
			for _ in 0..100 {
				let sk = SecretKey::random(&params);
				let rk = RelinearizationKey::new(&sk)?;

				// Let's generate manually an "extended" ciphertext (c0 = e - c1 * s - c2 * s^2, c1, c2) encrypting 0.
				let mut c2 = Poly::random(&params.ctx, Representation::Ntt);
				let mut c1 = Poly::random(&params.ctx, Representation::Ntt);
				let mut c0 = Poly::small(&params.ctx, Representation::PowerBasis, 16)?;
				c0.change_representation(Representation::Ntt);
				c0 -= &(&c1 * &sk.s);
				c0 -= &(&(&c2 * &sk.s) * &sk.s);

				// Relinearize the extended ciphertext!
				c2.change_representation(Representation::PowerBasis);
				rk.relinearize(&mut c0, &mut c1, &c2.coefficients())?;

				let ct = Ciphertext {
					par: params.clone(),
					seed: None,
					c0,
					c1,
				};

				// Print the noise and decrypt
				println!("Noise: {}", unsafe {
					sk.measure_noise(&ct, Encoding::Poly)?
				});
				let pt = sk.decrypt(&ct)?;
				assert!(Vec::<u64>::try_decode(&pt, Encoding::Poly).is_ok_and(|v| v == &[0u64; 8]))
			}
		}
		Ok(())
	}
}
