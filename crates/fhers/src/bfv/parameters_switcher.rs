//! Switching parameters

use crate::bfv::{BfvParameters, Ciphertext, Plaintext, SecretKey};
use crate::{Error, Result};
use fhers_traits::{FheParametersSwitchable, FheParametrized};
use itertools::Itertools;
use math::rq::Context;
use math::{
	rns::ScalingFactor,
	rq::{scaler::Scaler, traits::TryConvertFrom, Poly, Representation},
};
use std::sync::Arc;
use zeroize::Zeroize;

/// Switcher that enable to switch the underlying [`BfvParameters`]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BfvParametersSwitcher {
	pub(crate) from: Arc<BfvParameters>,
	pub(crate) to: Arc<BfvParameters>,
	pub(crate) scaler: Scaler,
	mod_switch_down_ctxs: Vec<Arc<Context>>,
}

impl FheParametrized for BfvParametersSwitcher {
	type Parameters = BfvParameters;
}

impl BfvParametersSwitcher {
	/// Create a switcher between two sets of parameters
	pub fn new(from: &Arc<BfvParameters>, to: &Arc<BfvParameters>) -> Result<Self> {
		let iter_mod_switch_down = if to.ctx.moduli().len() <= from.ctx.moduli().len()
			&& to.ctx.moduli() == &from.ctx.moduli()[..to.ctx.moduli().len()]
		{
			// We can do mod switch down instead of scaling
			from.ctx.moduli().len() - to.ctx.moduli().len()
		} else {
			0
		};
		let mut mod_switch_down_ctxs = vec![];
		for i in 0..iter_mod_switch_down {
			let ctx = Context::new(
				&from.ctx.moduli()[..from.ctx.moduli().len() - i - 1],
				from.degree(),
			)?;
			mod_switch_down_ctxs.push(Arc::new(ctx))
		}

		Ok(Self {
			from: from.clone(),
			to: to.clone(),
			scaler: Scaler::new(
				&from.ctx,
				&to.ctx,
				ScalingFactor::new(to.ctx.modulus(), from.ctx.modulus()),
			)?,
			mod_switch_down_ctxs,
		})
	}
}

impl FheParametersSwitchable<BfvParametersSwitcher> for Ciphertext {
	fn switch_parameters(&mut self, switcher: &BfvParametersSwitcher) -> Result<()> {
		if self.par != switcher.from {
			Err(Error::DefaultError("Mismatched parameters".to_string()))
		} else {
			self.seed = None;
			if !switcher.mod_switch_down_ctxs.is_empty() {
				self.c.iter_mut().for_each(|ci| {
					ci.change_representation(Representation::PowerBasis);
					for new_context in &switcher.mod_switch_down_ctxs {
						assert!(ci.mod_switch_down(new_context).is_ok())
					}
					ci.change_representation(Representation::Ntt)
				})
			} else {
				self.c = self
					.c
					.iter()
					.map(|ci| switcher.scaler.scale(ci).unwrap())
					.collect_vec();
			}
			self.par = switcher.to.clone();
			Ok(())
		}
	}

	type Error = Error;
}

impl FheParametersSwitchable<BfvParametersSwitcher> for SecretKey {
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

	type Error = Error;
}

impl FheParametersSwitchable<BfvParametersSwitcher> for Plaintext {
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
	type Error = Error;
}

#[cfg(test)]
mod tests {
	use super::BfvParametersSwitcher;
	use crate::bfv::{BfvParameters, Encoding, Plaintext, SecretKey};
	use fhers_traits::{
		FheDecoder, FheDecrypter, FheEncoder, FheEncrypter, FheParametersSwitchable,
	};
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
