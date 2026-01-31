//! Ciphertext type in the BFV encryption scheme.

use crate::bfv::{parameters::BfvParameters, traits::TryConvertFrom};
use crate::proto::bfv::Ciphertext as CiphertextProto;
use crate::{Error, Result, SerializationError};
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
    /// Create a ciphertext from a vector of polynomials.
    /// A ciphertext must contain at least two polynomials, and all polynomials
    /// must be in Ntt representation and with the same context.
    #[expect(clippy::expect_used, reason = "bounds are validated before use")]
    pub fn new(c: Vec<Poly>, par: &Arc<BfvParameters>) -> Result<Self> {
        if c.len() < 2 {
            return Err(Error::TooFewValues {
                actual: c.len(),
                minimum: 2,
            });
        }

        let ctx = c
            .first()
            .expect("c has at least 2 elements due to length check above")
            .ctx();
        let level = par.level_of_context(ctx)?;

        // Check that all polynomials have the expected representation and context.
        for ci in c.iter() {
            if ci.representation() != &Representation::Ntt {
                return Err(Error::MathError(fhe_math::Error::IncorrectRepresentation(
                    *ci.representation(),
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

    /// Truncate the underlying vector of polynomials.
    pub(crate) fn truncate(&mut self, len: usize) {
        self.c.truncate(len)
    }

    /// Switch to the next level in the chain
    pub fn switch_down(&mut self) -> Result<()> {
        if self.level < self.max_switchable_level() {
            self.seed = None;
            for ci in self.c.iter_mut() {
                ci.change_representation(Representation::PowerBasis);
                ci.switch_down()?;
                ci.change_representation(Representation::Ntt);
            }
            self.level += 1
        }
        Ok(())
    }

    /// Switch to a specific level (only moving down)
    pub fn switch_to_level(&mut self, target_level: usize) -> Result<()> {
        if target_level < self.level {
            return Err(Error::InvalidLevel {
                level: target_level,
                min_level: self.level,
                max_level: self.max_switchable_level(),
            });
        }
        if target_level > self.max_switchable_level() {
            return Err(Error::InvalidLevel {
                level: target_level,
                min_level: self.level,
                max_level: self.max_switchable_level(),
            });
        }
        while self.level < target_level {
            self.switch_down()?;
        }
        Ok(())
    }

    /// Get the deepest level this ciphertext can reach
    #[must_use]
    pub fn max_switchable_level(&self) -> usize {
        self.par.max_level()
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
        let ctp = Message::decode(bytes).map_err(|_| {
            Error::SerializationError(SerializationError::ProtobufError {
                message: "Ciphertext decode".into(),
            })
        })?;
        Ciphertext::try_convert_from(&ctp, par)
    }

    type Error = Error;
}

impl Ciphertext {
    /// Generate the zero ciphertext.
    #[must_use]
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

        // Split the ciphertext polynomials into all-but-last and last
        match ct.c.split_last() {
            None => {
                // Empty ciphertext - this should not happen as new() requires
                // at least 2 polys but we handle it gracefully
            }
            Some((last, rest)) => {
                // Serialize all but the last polynomial
                for poly in rest {
                    proto.c.push(poly.to_bytes());
                }

                // Handle the last polynomial based on whether we have a seed
                if let Some(seed) = ct.seed {
                    proto.seed = seed.to_vec();
                } else {
                    proto.c.push(last.to_bytes());
                }
            }
        }

        proto.level = ct.level as u32;
        proto
    }
}

impl TryConvertFrom<&CiphertextProto> for Ciphertext {
    fn try_convert_from(value: &CiphertextProto, par: &Arc<BfvParameters>) -> Result<Self> {
        if value.c.is_empty() || (value.c.len() == 1 && value.seed.is_empty()) {
            return Err(Error::InvalidCiphertext {
                reason: "Not enough polynomials".into(),
            });
        }

        if value.level as usize > par.max_level() {
            return Err(Error::InvalidLevel {
                level: value.level as usize,
                min_level: 0,
                max_level: par.max_level(),
            });
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
    use crate::Error as FheError;
    use crate::bfv::{
        BfvParameters, Ciphertext, Encoding, Plaintext, SecretKey, traits::TryConvertFrom,
    };
    use crate::proto::bfv::Ciphertext as CiphertextProto;
    use fhe_traits::FheDecrypter;
    use fhe_traits::{DeserializeParametrized, FheEncoder, FheEncrypter, Serialize};
    use rand::rng;
    use std::error::Error as StdError;

    #[test]
    fn proto_conversion() -> Result<(), Box<dyn StdError>> {
        let mut rng = rng();
        for params in [
            BfvParameters::default_arc(1, 16),
            BfvParameters::default_arc(6, 16),
        ] {
            let sk = SecretKey::random(&params, &mut rng);
            let v = fhe_math::zq::Modulus::new(params.plaintext()).unwrap().random_vec(params.degree(), &mut rng);
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
    fn serialize() -> Result<(), Box<dyn StdError>> {
        let mut rng = rng();
        for params in [
            BfvParameters::default_arc(1, 16),
            BfvParameters::default_arc(6, 16),
        ] {
            let sk = SecretKey::random(&params, &mut rng);
            let v = fhe_math::zq::Modulus::new(params.plaintext()).unwrap().random_vec(params.degree(), &mut rng);
            let pt = Plaintext::try_encode(&v, Encoding::simd(), &params)?;
            let ct: Ciphertext = sk.try_encrypt(&pt, &mut rng)?;
            let ct_bytes = ct.to_bytes();
            assert_eq!(ct, Ciphertext::from_bytes(&ct_bytes, &params)?);
        }
        Ok(())
    }

    #[test]
    fn new() -> Result<(), Box<dyn StdError>> {
        let mut rng = rng();
        for params in [
            BfvParameters::default_arc(1, 16),
            BfvParameters::default_arc(6, 16),
        ] {
            let sk = SecretKey::random(&params, &mut rng);
            let v = fhe_math::zq::Modulus::new(params.plaintext()).unwrap().random_vec(params.degree(), &mut rng);
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

            ct3.switch_to_level(ct3.max_switchable_level())?;

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
    fn switch_to_last_level() -> Result<(), Box<dyn StdError>> {
        let mut rng = rng();
        for params in [
            BfvParameters::default_arc(1, 16),
            BfvParameters::default_arc(6, 16),
        ] {
            let sk = SecretKey::random(&params, &mut rng);
            let v = fhe_math::zq::Modulus::new(params.plaintext()).unwrap().random_vec(params.degree(), &mut rng);
            let pt = Plaintext::try_encode(&v, Encoding::simd(), &params)?;
            let mut ct: Ciphertext = sk.try_encrypt(&pt, &mut rng)?;

            assert_eq!(ct.level, 0);
            ct.switch_to_level(ct.max_switchable_level())?;
            assert_eq!(ct.level, params.max_level());

            let decrypted = sk.try_decrypt(&ct)?;
            assert_eq!(decrypted.value, pt.value);
        }

        Ok(())
    }

    #[test]
    #[expect(clippy::panic, reason = "panic indicates violated internal invariant")]
    fn switch_to_level_invalid() -> Result<(), Box<dyn StdError>> {
        let mut rng = rng();
        let params = BfvParameters::default_arc(2, 16);
        let sk = SecretKey::random(&params, &mut rng);
        let v = fhe_math::zq::Modulus::new(params.plaintext()).unwrap().random_vec(params.degree(), &mut rng);
        let pt = Plaintext::try_encode(&v, Encoding::simd(), &params)?;
        let mut ct: Ciphertext = sk.try_encrypt(&pt, &mut rng)?;

        // Move to level 1
        ct.switch_down()?;
        assert_eq!(ct.level, 1);

        // Target level smaller than current
        match ct.switch_to_level(0) {
            Err(FheError::InvalidLevel {
                level,
                min_level,
                max_level,
            }) => {
                assert_eq!(level, 0);
                assert_eq!(min_level, 1);
                assert_eq!(max_level, params.max_level());
            }
            _ => panic!("expected InvalidLevel error"),
        }

        // Target level larger than max
        let too_high = params.max_level() + 1;
        match ct.switch_to_level(too_high) {
            Err(FheError::InvalidLevel {
                level,
                min_level,
                max_level,
            }) => {
                assert_eq!(level, too_high);
                assert_eq!(min_level, 1);
                assert_eq!(max_level, params.max_level());
            }
            _ => panic!("expected InvalidLevel error"),
        }

        Ok(())
    }
}
