//! Ciphertext type in the BFV encryption scheme.

use crate::bfv::{context_chain::ContextLevel, parameters::BfvParameters, traits::TryConvertFrom};
use crate::proto::bfv::Ciphertext as CiphertextProto;
use crate::{Error, Result};
use fhe_math::rq::{Poly, Representation};
use fhe_traits::{
    DeserializeParametrized, DeserializeWithContext, FheCiphertext, FheParametrized, Serialize,
};
use prost::Message;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use std::ops::{Deref, DerefMut};
use std::sync::Arc;

/// A ciphertext encrypting a plaintext.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Ciphertext {
    /// The parameters of the underlying BFV encryption scheme.
    pub(crate) par: Arc<BfvParameters>,

    /// The seed that generated the polynomial c1 in a fresh ciphertext.
    pub(crate) seed: Option<<ChaCha8Rng as SeedableRng>::Seed>,

    /// The ciphertext elements.
    pub(crate) c: Vec<Poly>,

    /// The ciphertext level
    pub(crate) level: usize,
}

impl Deref for Ciphertext {
    type Target = [Poly];

    fn deref(&self) -> &Self::Target {
        &self.c
    }
}

impl DerefMut for Ciphertext {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.c
    }
}

impl Ciphertext {
    /// Modulo switch the ciphertext to the last level.
    pub fn mod_switch_to_last_level(&mut self) -> Result<()> {
        self.level = self.par.max_level();
        let last_ctx = self.par.context_at_level(self.level)?;
        self.seed = None;
        for ci in self.c.iter_mut() {
            if ci.ctx() != last_ctx {
                ci.change_representation(Representation::PowerBasis);
                ci.mod_switch_down_to(last_ctx)?;
                ci.change_representation(Representation::Ntt);
            }
        }
        Ok(())
    }

    /// Truncate the underlying vector of polynomials.
    pub(crate) fn truncate(&mut self, len: usize) {
        self.c.truncate(len)
    }

    /// Modulo switch the ciphertext to the next level.
    pub fn mod_switch_to_next_level(&mut self) -> Result<()> {
        if self.level < self.par.max_level() {
            self.seed = None;
            for ci in self.c.iter_mut() {
                ci.change_representation(Representation::PowerBasis);
                ci.mod_switch_down_next()?;
                ci.change_representation(Representation::Ntt);
            }
            self.level += 1
        }
        Ok(())
    }

    /// Get the context level for this ciphertext
    pub fn context_level(&self) -> Arc<ContextLevel> {
        // safe unwrap: parameters always contain context chain
        self.par.context_level_at(self.level).unwrap()
    }

    /// Switch to the next level in the chain
    pub fn switch_down(&mut self) -> Result<()> {
        self.mod_switch_to_next_level()
    }

    /// Switch to a specific level (only moving down)
    pub fn switch_to_level(&mut self, target_level: usize) -> Result<()> {
        if target_level < self.level {
            return Err(Error::DefaultError(format!(
                "Cannot switch to a higher level: current {}, target {}",
                self.level, target_level
            )));
        }
        if target_level > self.par.max_level() {
            return Err(Error::DefaultError(format!(
                "Invalid level: {target_level}"
            )));
        }
        while self.level < target_level {
            self.mod_switch_to_next_level()?;
        }
        Ok(())
    }

    /// Get the deepest level this ciphertext can reach
    pub fn max_switchable_level(&self) -> usize {
        self.par.max_level()
    }

    /// Clone the ciphertext and switch to a specific level
    pub fn clone_at_level(&self, level: usize) -> Result<Self> {
        let mut cloned = self.clone();
        cloned.switch_to_level(level)?;
        Ok(cloned)
    }

    /// Create a ciphertext from a vector of polynomials.
    /// A ciphertext must contain at least two polynomials, and all polynomials
    /// must be in Ntt representation and with the same context.
    pub fn new(c: Vec<Poly>, par: &Arc<BfvParameters>) -> Result<Self> {
        if c.len() < 2 {
            return Err(Error::TooFewValues(c.len(), 2));
        }

        let ctx = c[0].ctx();
        let level = par.level_of_context(ctx)?;

        // Check that all polynomials have the expected representation and context.
        for ci in c.iter() {
            if ci.representation() != &Representation::Ntt {
                return Err(Error::MathError(fhe_math::Error::IncorrectRepresentation(
                    ci.representation().clone(),
                    Representation::Ntt,
                )));
            }
            if ci.ctx() != ctx {
                return Err(Error::MathError(fhe_math::Error::InvalidContext));
            }
        }

        Ok(Self {
            par: par.clone(),
            seed: None,
            c,
            level,
        })
    }
}

impl FheCiphertext for Ciphertext {}

impl FheParametrized for Ciphertext {
    type Parameters = BfvParameters;
}

impl Serialize for Ciphertext {
    fn to_bytes(&self) -> Vec<u8> {
        CiphertextProto::from(self).encode_to_vec()
    }
}

impl DeserializeParametrized for Ciphertext {
    fn from_bytes(bytes: &[u8], par: &Arc<BfvParameters>) -> Result<Self> {
        if let Ok(ctp) = Message::decode(bytes) {
            Ciphertext::try_convert_from(&ctp, par)
        } else {
            Err(Error::SerializationError)
        }
    }

    type Error = Error;
}

impl Ciphertext {
    /// Generate the zero ciphertext.
    pub fn zero(par: &Arc<BfvParameters>) -> Self {
        Self {
            par: par.clone(),
            seed: None,
            c: Default::default(),
            level: 0,
        }
    }
}

/// Conversions from and to protobuf.
impl From<&Ciphertext> for CiphertextProto {
    fn from(ct: &Ciphertext) -> Self {
        let mut proto = CiphertextProto::default();
        for i in 0..ct.len() - 1 {
            proto.c.push(ct[i].to_bytes())
        }
        if let Some(seed) = ct.seed {
            proto.seed = seed.to_vec()
        } else {
            proto.c.push(ct[ct.len() - 1].to_bytes())
        }
        proto.level = ct.level as u32;
        proto
    }
}

impl TryConvertFrom<&CiphertextProto> for Ciphertext {
    fn try_convert_from(value: &CiphertextProto, par: &Arc<BfvParameters>) -> Result<Self> {
        if value.c.is_empty() || (value.c.len() == 1 && value.seed.is_empty()) {
            return Err(Error::DefaultError("Not enough polynomials".to_string()));
        }

        if value.level as usize > par.max_level() {
            return Err(Error::DefaultError("Invalid level".to_string()));
        }

        let ctx = par.context_at_level(value.level as usize)?;

        let mut c = Vec::with_capacity(value.c.len() + 1);
        for cip in &value.c {
            c.push(Poly::from_bytes(cip, ctx)?)
        }

        let mut seed = None;
        if !value.seed.is_empty() {
            let try_seed = <ChaCha8Rng as SeedableRng>::Seed::try_from(value.seed.clone())
                .map_err(|_| {
                    Error::MathError(fhe_math::Error::InvalidSeedSize(
                        value.seed.len(),
                        <ChaCha8Rng as SeedableRng>::Seed::default().len(),
                    ))
                })?;
            seed = Some(try_seed);
            let mut c1 = Poly::random_from_seed(ctx, Representation::Ntt, try_seed);
            unsafe { c1.allow_variable_time_computations() }
            c.push(c1)
        }

        Ok(Ciphertext {
            par: par.clone(),
            seed,
            c,
            level: value.level as usize,
        })
    }
}

#[cfg(test)]
mod tests {
    use crate::bfv::{
        traits::TryConvertFrom, BfvParameters, Ciphertext, Encoding, Plaintext, SecretKey,
    };
    use crate::proto::bfv::Ciphertext as CiphertextProto;
    use fhe_traits::FheDecrypter;
    use fhe_traits::{DeserializeParametrized, FheEncoder, FheEncrypter, Serialize};
    use rand::thread_rng;
    use std::error::Error;

    #[test]
    fn proto_conversion() -> Result<(), Box<dyn Error>> {
        let mut rng = thread_rng();
        for params in [
            BfvParameters::default_arc(1, 16),
            BfvParameters::default_arc(6, 16),
        ] {
            let sk = SecretKey::random(&params, &mut rng);
            let v = params.plaintext.random_vec(params.degree(), &mut rng);
            let pt = Plaintext::try_encode(&v, Encoding::simd(), &params)?;
            let ct = sk.try_encrypt(&pt, &mut rng)?;
            let ct_proto = CiphertextProto::from(&ct);
            assert_eq!(ct, Ciphertext::try_convert_from(&ct_proto, &params)?);

            let ct = &ct * &ct;
            let ct_proto = CiphertextProto::from(&ct);
            assert_eq!(ct, Ciphertext::try_convert_from(&ct_proto, &params)?)
        }
        Ok(())
    }

    #[test]
    fn serialize() -> Result<(), Box<dyn Error>> {
        let mut rng = thread_rng();
        for params in [
            BfvParameters::default_arc(1, 16),
            BfvParameters::default_arc(6, 16),
        ] {
            let sk = SecretKey::random(&params, &mut rng);
            let v = params.plaintext.random_vec(params.degree(), &mut rng);
            let pt = Plaintext::try_encode(&v, Encoding::simd(), &params)?;
            let ct: Ciphertext = sk.try_encrypt(&pt, &mut rng)?;
            let ct_bytes = ct.to_bytes();
            assert_eq!(ct, Ciphertext::from_bytes(&ct_bytes, &params)?);
        }
        Ok(())
    }

    #[test]
    fn new() -> Result<(), Box<dyn Error>> {
        let mut rng = thread_rng();
        for params in [
            BfvParameters::default_arc(1, 16),
            BfvParameters::default_arc(6, 16),
        ] {
            let sk = SecretKey::random(&params, &mut rng);
            let v = params.plaintext.random_vec(params.degree(), &mut rng);
            let pt = Plaintext::try_encode(&v, Encoding::simd(), &params)?;
            let ct: Ciphertext = sk.try_encrypt(&pt, &mut rng)?;
            let mut ct3 = &ct * &ct;

            let c0 = &ct3[0];
            let c1 = &ct3[1];
            let c2 = &ct3[2];

            assert_eq!(
                ct3,
                Ciphertext::new(vec![c0.clone(), c1.clone(), c2.clone()], &params)?
            );
            assert_eq!(ct3.level, 0);

            ct3.mod_switch_to_last_level()?;

            let c0 = ct3.first().unwrap();
            let c1 = ct3.get(1).unwrap();
            let c2 = ct3.get(2).unwrap();
            assert_eq!(
                ct3,
                Ciphertext::new(vec![c0.clone(), c1.clone(), c2.clone()], &params)?
            );
            assert_eq!(ct3.level, params.max_level());
        }

        Ok(())
    }

    #[test]
    fn mod_switch_to_last_level() -> Result<(), Box<dyn Error>> {
        let mut rng = thread_rng();
        for params in [
            BfvParameters::default_arc(1, 16),
            BfvParameters::default_arc(6, 16),
        ] {
            let sk = SecretKey::random(&params, &mut rng);
            let v = params.plaintext.random_vec(params.degree(), &mut rng);
            let pt = Plaintext::try_encode(&v, Encoding::simd(), &params)?;
            let mut ct: Ciphertext = sk.try_encrypt(&pt, &mut rng)?;

            assert_eq!(ct.level, 0);
            ct.mod_switch_to_last_level()?;
            assert_eq!(ct.level, params.max_level());

            let decrypted = sk.try_decrypt(&ct)?;
            assert_eq!(decrypted.value, pt.value);
        }

        Ok(())
    }
}
