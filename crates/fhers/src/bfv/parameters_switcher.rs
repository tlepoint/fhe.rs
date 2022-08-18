//! Switching parameters

use crate::bfv::{BfvParameters, Ciphertext, Plaintext, SecretKey};
use crate::{Error, Result};
use itertools::Itertools;
use math::{
	rns::ScalingFactor,
	rq::{scaler::Scaler, traits::TryConvertFrom, Poly, Representation},
};
use std::sync::Arc;
use zeroize::Zeroize;

/// Switcher that enable to switch the underlying [`BfvParameters`]
pub struct BfvParametersSwitcher {
	pub(crate) from: Arc<BfvParameters>,
	pub(crate) to: Arc<BfvParameters>,
	scaler: Scaler,
}

impl BfvParametersSwitcher {
	/// Create a switcher between two sets of parameters
	pub fn new(from: &Arc<BfvParameters>, to: &Arc<BfvParameters>) -> Result<Self> {
		Ok(Self {
			from: from.clone(),
			to: to.clone(),
			scaler: Scaler::new(
				&from.ctx,
				&to.ctx,
				ScalingFactor::new(to.ctx.modulus(), from.ctx.modulus()),
			)?,
		})
	}
}

/// Indicate that the [`BfvParameters`] of Self can be switched.
pub trait ParametersSwitchable {
	/// Switch the underlying [`BfvParameters`]
	fn switch_parameters(&mut self, switcher: &BfvParametersSwitcher) -> Result<()>;
}

impl ParametersSwitchable for Ciphertext {
	fn switch_parameters(&mut self, switcher: &BfvParametersSwitcher) -> Result<()> {
		if self.par != switcher.from {
			Err(Error::DefaultError("Mismatched parameters".to_string()))
		} else {
			self.seed = None;
			self.c = self
				.c
				.iter()
				.map(|ci| switcher.scaler.scale(ci, false).unwrap())
				.collect_vec();
			self.par = switcher.to.clone();
			Ok(())
		}
	}
}

impl ParametersSwitchable for SecretKey {
	fn switch_parameters(&mut self, switcher: &BfvParametersSwitcher) -> Result<()> {
		if self.par != switcher.from {
			Err(Error::DefaultError("Mismatched parameters".to_string()))
		} else {
			let s_coefficients = self.s_coefficients.clone();
			self.zeroize();
			*self = SecretKey::new(s_coefficients, &switcher.to);
			Ok(())
		}
	}
}

impl ParametersSwitchable for Plaintext {
	fn switch_parameters(&mut self, switcher: &BfvParametersSwitcher) -> Result<()> {
		if self.par != switcher.from {
			Err(Error::DefaultError("Mismatched parameters".to_string()))
		} else {
			self.poly_ntt = Poly::try_convert_from(
				&self.value as &[u64],
				&switcher.to.ctx,
				false,
				Representation::PowerBasis,
			)?;
			self.poly_ntt.change_representation(Representation::Ntt);
			self.par = switcher.to.clone();
			Ok(())
		}
	}
}

#[cfg(test)]
mod tests {
	use super::{BfvParametersSwitcher, ParametersSwitchable};
	use crate::bfv::{BfvParameters, Encoding, Plaintext, SecretKey};
	use fhers_traits::{FheDecoder, FheDecrypter, FheEncoder, FheEncrypter};
	use std::{error::Error, sync::Arc};

	#[test]
	fn test_parameters_switcher() -> Result<(), Box<dyn Error>> {
		let from = Arc::new(BfvParameters::default(5));
		let to = Arc::new(BfvParameters::default(2));

		let switcher = BfvParametersSwitcher::new(&from, &to)?;

		let v = from.plaintext.random_vec(from.degree());
		let mut sk = SecretKey::random(&from);
		let mut pt = Plaintext::try_encode(&v as &[u64], Encoding::Poly, &from)?;
		let mut ct = sk.try_encrypt(&pt)?;

		sk.switch_parameters(&switcher)?;
		ct.switch_parameters(&switcher)?;
		pt.switch_parameters(&switcher)?;

		assert_eq!(sk.par, to);
		assert_eq!(ct.par, to);
		assert_eq!(pt.par, to);

		println!("Noise: {}", unsafe { sk.measure_noise(&ct)? });
		let pt_decrypted = sk.try_decrypt(&ct)?;
		assert_eq!(pt, pt_decrypted);
		assert_eq!(Vec::<u64>::try_decode(&pt, Encoding::Poly)?, v);

		Ok(())
	}
}
