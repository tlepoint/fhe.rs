//! Ciphertext type in the BFV encryption scheme.

use crate::bfv::{parameters::Parameters, wire::FromProto};
use crate::proto::bfv::Ciphertext as CiphertextProto;
use crate::{Error, Result, error::SerializationError};
use fhe_math::rq::{Context, Ntt, Poly};

use prost::Message;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use std::sync::Arc;

/// A ciphertext encrypting a plaintext.
///
/// Polynomial components are immutable through the public API.
/// ```compile_fail
/// use fhe::bfv::Ciphertext;
/// use fhe_math::rq::{Poly, Ntt};
/// fn replace(ct: &mut Ciphertext, polynomial: Poly<Ntt>) {
///     ct.components()[0] = polynomial;
/// }
/// ```
/// The ciphertext is not a slice or a mutable smart pointer to its components.
/// ```compile_fail
/// use fhe::bfv::Ciphertext;
/// use fhe_math::rq::{Poly, Ntt};
/// fn replace(ct: &mut Ciphertext, polynomial: Poly<Ntt>) {
///     ct[0] = polynomial;
/// }
/// ```
///
/// Equality compares parameters, level and polynomial components. It ignores
/// the compression seed and does not test equality of encrypted messages.
#[derive(Debug, Clone, Eq)]
pub struct Ciphertext {
    /// The parameters of the underlying BFV encryption scheme.
    pub(crate) par: Parameters,

    /// The seed that generated the polynomial c1 in a fresh ciphertext.
    pub(crate) seed: Option<<ChaCha8Rng as SeedableRng>::Seed>,

    /// The ciphertext elements.
    pub(crate) c: Vec<Poly<Ntt>>,

    /// The ciphertext level
    pub(crate) level: usize,
}

impl PartialEq for Ciphertext {
    fn eq(&self, other: &Self) -> bool {
        self.par == other.par && self.level == other.level && self.c == other.c
    }
}

impl Ciphertext {
    /// Create a ciphertext from a vector of polynomials.
    /// A ciphertext must contain at least two polynomials, and all polynomials
    /// must have canonical Ntt residues and the same context. The context must
    /// belong to the supplied parameters; the level is inferred from it.
    #[expect(clippy::expect_used, reason = "bounds are validated before use")]
    pub fn from_components(c: Vec<Poly<Ntt>>, par: &Parameters) -> Result<Self> {
        if c.len() < 2 {
            return Err(crate::error::CiphertextError::TooFewPolynomials {
                actual: c.len(),
                minimum: 2,
            }
            .into());
        }

        let ctx = c
            .first()
            .expect("c has at least 2 elements due to length check above")
            .ctx();
        let level = par.level_of_context(ctx)?;

        // Check that all polynomials have the expected context.
        for ci in c.iter() {
            if !ci.is_canonical() {
                return Err(crate::error::CiphertextError::NonCanonicalPolynomial.into());
            }
            if ci.ctx() != ctx {
                return Err(
                    crate::error::CiphertextError::PolynomialContextMismatch { level }.into(),
                );
            }
        }

        Ok(Self {
            par: par.clone(),
            seed: None,
            c,
            level,
        })
    }

    /// Validate the structure and context of a ciphertext used as an input.
    #[inline]
    pub(crate) fn validate_for(&self, par: &Parameters) -> Result<()> {
        if !Parameters::compatible(&self.par, par) {
            return Err(Error::ParameterMismatch {
                left: crate::error::ParameterSource::Ciphertext,
                right: crate::error::ParameterSource::Parameters,
            });
        }
        let expected_ctx = par.context_at_level(self.level)?;
        self.validate_context(self.level, expected_ctx)
    }

    /// Validate against a context that the caller has already resolved.
    #[inline]
    pub(crate) fn validate_for_context(
        &self,
        par: &Parameters,
        expected_level: usize,
        expected_ctx: &Arc<Context>,
    ) -> Result<()> {
        if !Parameters::compatible(&self.par, par) {
            return Err(Error::ParameterMismatch {
                left: crate::error::ParameterSource::Ciphertext,
                right: crate::error::ParameterSource::Parameters,
            });
        }
        self.validate_context(expected_level, expected_ctx)
    }

    #[inline]
    fn validate_context(&self, expected_level: usize, expected_ctx: &Arc<Context>) -> Result<()> {
        if self.c.iter().any(|poly| !poly.is_canonical()) {
            return Err(crate::error::CiphertextError::NonCanonicalPolynomial.into());
        }
        if self.c.len() < 2 {
            return Err(crate::error::CiphertextError::TooFewPolynomials {
                actual: self.c.len(),
                minimum: 2,
            }
            .into());
        }
        if self.level != expected_level {
            return Err(Error::InvalidLevel {
                level: self.level,
                min_level: expected_level,
                max_level: expected_level,
            });
        }
        if self
            .c
            .iter()
            .any(|poly| !Arc::ptr_eq(poly.ctx(), expected_ctx) && poly.ctx() != expected_ctx)
        {
            return Err(crate::error::CiphertextError::PolynomialContextMismatch {
                level: expected_level,
            }
            .into());
        }
        Ok(())
    }

    /// Truncate the underlying vector of polynomials.
    pub(crate) fn truncate(&mut self, len: usize) {
        self.seed = None;
        self.c.truncate(len)
    }

    /// Switch to the next level in the chain.
    ///
    /// Returns an error if the ciphertext is already at the last level.
    pub fn switch_down(&mut self) -> Result<()> {
        if self.level >= self.max_switchable_level() {
            return Err(fhe_math::Error::NoMoreContext.into());
        }

        self.switch_to_level(self.level + 1)
    }

    /// Switch to a specific level (only moving down).
    ///
    /// One drop preserves the surviving NTT rows. Several drops use one
    /// conversion to power basis and back, retaining the same sequence of
    /// rounding operations. Invalid inputs return before changing any part.
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
        self.validate_for(&self.par)?;
        if self.level != target_level {
            let target = self.par.context_at_level(target_level)?;
            self.seed = None;
            for ci in &mut self.c {
                ci.switch_down_to(target)?;
            }
            self.level = target_level;
        }
        Ok(())
    }

    /// Get the deepest level this ciphertext can reach
    #[must_use]
    pub fn max_switchable_level(&self) -> usize {
        self.par.max_level()
    }
}

impl Ciphertext {
    /// Serialize in the existing protobuf wire format.
    #[must_use]
    pub fn to_bytes(&self) -> Vec<u8> {
        CiphertextProto::from(self).encode_to_vec()
    }
}

impl Ciphertext {
    /// Import validated protobuf bytes, binding contextual values to the
    /// supplied parameters.
    pub fn from_bytes(bytes: &[u8], par: &Parameters) -> Result<Self> {
        Self::from_bytes_with_limits(bytes, par, &crate::DecodeLimits::default())
    }

    /// Import with explicit resource bounds checked before allocation.
    pub fn from_bytes_with_limits(
        bytes: &[u8],
        par: &Parameters,
        limits: &crate::DecodeLimits,
    ) -> Result<Self> {
        crate::bfv::wire::preflight(
            bytes,
            crate::error::SerializedObject::Ciphertext,
            Some(par),
            limits,
        )?;
        let ctp = Message::decode(bytes).map_err(|source| {
            Error::SerializationError(SerializationError::Decode {
                object: crate::error::SerializedObject::Ciphertext,
                source,
            })
        })?;
        Ciphertext::from_proto(&ctp, par, limits)
    }
}

impl Ciphertext {
    /// Construct a public, deterministic zero at a validated level.
    ///
    /// This has two polynomial components and can be evaluated, decrypted and
    /// serialized like any other ciphertext. It does not hide its value;
    /// encrypt a zero plaintext with a key and RNG to obtain a randomized zero.
    pub fn trivial_zero(par: &Parameters, level: usize) -> Result<Self> {
        let ctx = par.context_at_level(level)?;
        let mut zero = Poly::<Ntt>::zero(ctx);
        zero.allow_variable_time_computations(crate::VariableTime::new(
            crate::PublicData::assert_public(),
        ));
        Self::from_components(vec![zero; 2], par)
    }

    /// Borrow the ciphertext's parameters.
    #[must_use]
    pub fn parameters(&self) -> &Parameters {
        &self.par
    }

    /// Return the modulus-chain level (zero is the full chain).
    #[must_use]
    pub const fn level(&self) -> usize {
        self.level
    }

    /// Return the number of polynomial components.
    #[must_use]
    pub fn component_count(&self) -> usize {
        self.c.len()
    }

    /// Borrow the polynomial components for advanced arithmetic or inspection.
    /// Mutation requires rebuilding through [`Self::from_components`].
    #[must_use]
    pub fn components(&self) -> &[Poly<Ntt>] {
        &self.c
    }

    /// Consume this ciphertext, discarding its compression seed.
    /// Reconstruct modified components with [`Self::from_components`].
    #[must_use]
    pub fn into_components(self) -> Vec<Poly<Ntt>> {
        self.c
    }

    pub(crate) fn len(&self) -> usize {
        self.c.len()
    }
    pub(crate) fn iter(&self) -> std::slice::Iter<'_, Poly<Ntt>> {
        self.c.iter()
    }
    pub(crate) fn iter_mut(&mut self) -> std::slice::IterMut<'_, Poly<Ntt>> {
        self.seed = None;
        self.c.iter_mut()
    }

    #[cfg(test)]
    pub(crate) fn invalid_empty(par: &Parameters) -> Self {
        Self {
            par: par.clone(),
            seed: None,
            c: vec![],
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
                // Only malformed crate-internal test fixtures can be empty;
                // public construction requires at least two components.
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

impl FromProto<&CiphertextProto> for Ciphertext {
    fn from_proto(
        value: &CiphertextProto,
        par: &Parameters,
        limits: &crate::DecodeLimits,
    ) -> Result<Self> {
        if value.c.is_empty() || (value.c.len() == 1 && value.seed.is_empty()) {
            return Err(Error::SerializationError(
                SerializationError::InvalidCiphertextPolynomialCount {
                    actual: value.c.len(),
                    seed_present: !value.seed.is_empty(),
                },
            ));
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
            c.push(Poly::<Ntt>::from_bytes_with_limits(cip, ctx, limits)?)
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
            let mut c1 = Poly::<Ntt>::random_from_seed(ctx, try_seed);
            c1.allow_variable_time_computations(crate::VariableTime::new(
                crate::PublicData::assert_public(),
            ));
            c.push(c1)
        }

        // Ciphertexts are public once received. Grant timing permission only
        // at this trusted type boundary; polynomial wire data cannot grant it.
        let variable_time = crate::VariableTime::new(crate::PublicData::assert_public());
        c.iter_mut()
            .for_each(|ci| ci.allow_variable_time_computations(variable_time));

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
    use crate::bfv::{Ciphertext, Encoding, Parameters, Plaintext, SecretKey, wire::FromProto};
    use crate::proto::bfv::Ciphertext as CiphertextProto;

    use rand::rng;
    use std::error::Error as StdError;

    #[test]
    fn proto_conversion() -> Result<(), Box<dyn StdError>> {
        let mut rng = rng();
        for params in [
            Parameters::test_parameters(1, 16),
            Parameters::test_parameters(6, 16),
        ] {
            let sk = SecretKey::generate(&params, &mut rng);
            let v = fhe_math::zq::Modulus::new(params.plaintext_modulus_u64().unwrap())
                .unwrap()
                .random_vec(params.degree(), &mut rng);
            let pt = Plaintext::encode(&params, &v, Encoding::Simd)?;
            let ct = sk.encrypt(&pt, &mut rng)?;
            let ct_proto = CiphertextProto::from(&ct);
            assert_eq!(
                ct,
                Ciphertext::from_proto(&ct_proto, &params, &crate::DecodeLimits::default())?
            );

            let ct = ct.multiply(&ct).unwrap();
            let ct_proto = CiphertextProto::from(&ct);
            assert_eq!(
                ct,
                Ciphertext::from_proto(&ct_proto, &params, &crate::DecodeLimits::default())?
            )
        }
        Ok(())
    }

    #[test]
    fn serialize() -> Result<(), Box<dyn StdError>> {
        let mut rng = rng();
        for params in [
            Parameters::test_parameters(1, 16),
            Parameters::test_parameters(6, 16),
        ] {
            let sk = SecretKey::generate(&params, &mut rng);
            let v = fhe_math::zq::Modulus::new(params.plaintext_modulus_u64().unwrap())
                .unwrap()
                .random_vec(params.degree(), &mut rng);
            let pt = Plaintext::encode(&params, &v, Encoding::Simd)?;
            let ct: Ciphertext = sk.encrypt(&pt, &mut rng)?;
            let ct_bytes = ct.to_bytes();
            assert_eq!(ct, Ciphertext::from_bytes(&ct_bytes, &params)?);
        }
        Ok(())
    }

    #[test]
    fn new() -> Result<(), Box<dyn StdError>> {
        let mut rng = rng();
        for params in [
            Parameters::test_parameters(1, 16),
            Parameters::test_parameters(6, 16),
        ] {
            let sk = SecretKey::generate(&params, &mut rng);
            let v = fhe_math::zq::Modulus::new(params.plaintext_modulus_u64().unwrap())
                .unwrap()
                .random_vec(params.degree(), &mut rng);
            let pt = Plaintext::encode(&params, &v, Encoding::Simd)?;
            let ct: Ciphertext = sk.encrypt(&pt, &mut rng)?;
            let mut ct3 = ct.multiply(&ct).unwrap();

            let c0 = &ct3.c[0];
            let c1 = &ct3.c[1];
            let c2 = &ct3.c[2];

            assert_eq!(
                ct3,
                Ciphertext::from_components(vec![c0.clone(), c1.clone(), c2.clone()], &params)?
            );
            assert_eq!(ct3.level, 0);

            ct3.switch_to_level(ct3.max_switchable_level())?;

            let c0 = ct3.c.first().unwrap();
            let c1 = ct3.c.get(1).unwrap();
            let c2 = ct3.c.get(2).unwrap();
            assert_eq!(
                ct3,
                Ciphertext::from_components(vec![c0.clone(), c1.clone(), c2.clone()], &params)?
            );
            assert_eq!(ct3.level, params.max_level());
        }

        Ok(())
    }

    #[test]
    fn ntt_switching_matches_the_old_round_trip_at_every_level() -> Result<(), Box<dyn StdError>> {
        use fhe_math::rq::{Ntt, Poly};
        use rand::SeedableRng;
        use rand_chacha::ChaCha8Rng;
        let par = Parameters::test_parameters(4, 16);
        let mut rng = ChaCha8Rng::seed_from_u64(0x57017c4);
        for level in 0..=par.max_level() {
            for parts in 2..=4 {
                for restricted in [false, true] {
                    let polynomials = (0..parts)
                        .map(|i| {
                            let mut p =
                                Poly::<Ntt>::random(par.context_at_level(level).unwrap(), &mut rng);
                            if !restricted || i != 1 {
                                p.allow_variable_time_computations(crate::VariableTime::new(
                                    crate::PublicData::assert_public(),
                                ));
                            }
                            p
                        })
                        .collect();
                    let original = Ciphertext::from_components(polynomials, &par)?;
                    for target in level..=par.max_level() {
                        let mut expected = original.clone();
                        while expected.level < target {
                            for p in &mut expected.c {
                                let mut pb = p.clone().into_power_basis();
                                pb.switch_down()?;
                                *p = pb.into_ntt();
                            }
                            expected.seed = None;
                            expected.level += 1;
                        }
                        let mut actual = original.clone();
                        actual.switch_to_level(target)?;
                        assert_eq!(actual, expected);
                        assert_eq!(actual.to_bytes(), expected.to_bytes());
                        for (i, p) in actual.iter().enumerate() {
                            assert_eq!(
                                p.allows_variable_time_computations(),
                                !restricted || i != 1
                            );
                        }
                        if target == level + 1 {
                            let mut single = original.clone();
                            single.switch_down()?;
                            assert_eq!(single, expected);
                        }
                    }
                }
            }
        }
        // A malformed later part must not leave the earlier part switched.
        let ctx = par.context_at_level(0)?;
        let mut malformed = Ciphertext::from_components(vec![Poly::zero(ctx); 2], &par)?;
        malformed.c[1] = Poly::zero(par.context_at_level(1)?);
        let saved = malformed.clone();
        assert!(malformed.switch_down().is_err());
        assert_eq!(malformed, saved);
        assert!(malformed.switch_to_level(2).is_err());
        assert_eq!(malformed, saved);
        let mut zero = Ciphertext::trivial_zero(&par, 0)?;
        zero.switch_to_level(par.max_level())?;
        assert_eq!(zero.component_count(), 2);
        assert_eq!(zero.level, par.max_level());
        Ok(())
    }

    #[test]
    fn switch_to_last_level() -> Result<(), Box<dyn StdError>> {
        let mut rng = rng();
        for params in [
            Parameters::test_parameters(1, 16),
            Parameters::test_parameters(6, 16),
        ] {
            let sk = SecretKey::generate(&params, &mut rng);
            let v = fhe_math::zq::Modulus::new(params.plaintext_modulus_u64().unwrap())
                .unwrap()
                .random_vec(params.degree(), &mut rng);
            let pt = Plaintext::encode(&params, &v, Encoding::Simd)?;
            let mut ct: Ciphertext = sk.encrypt(&pt, &mut rng)?;

            assert_eq!(ct.level, 0);
            ct.switch_to_level(ct.max_switchable_level())?;
            assert_eq!(ct.level, params.max_level());

            let decrypted = sk.decrypt(&ct)?;
            assert_eq!(
                decrypted.decode(Encoding::Simd)?,
                pt.decode(Encoding::Simd)?
            );
        }

        Ok(())
    }

    #[test]
    fn switch_down_from_last_level_returns_error() -> Result<(), Box<dyn StdError>> {
        let mut rng = rng();
        let params = Parameters::test_parameters(2, 16);
        let sk = SecretKey::generate(&params, &mut rng);
        let pt = Plaintext::encode(&params, &[1u64][..], Encoding::Polynomial)?;
        let mut ct: Ciphertext = sk.encrypt(&pt, &mut rng)?;
        ct.switch_to_level(params.max_level())?;

        assert!(matches!(
            ct.switch_down(),
            Err(FheError::MathError(fhe_math::Error::NoMoreContext))
        ));
        Ok(())
    }

    #[test]
    #[expect(clippy::panic, reason = "panic indicates violated internal invariant")]
    fn switch_to_level_invalid() -> Result<(), Box<dyn StdError>> {
        let mut rng = rng();
        let params = Parameters::test_parameters(2, 16);
        let sk = SecretKey::generate(&params, &mut rng);
        let v = fhe_math::zq::Modulus::new(params.plaintext_modulus_u64().unwrap())
            .unwrap()
            .random_vec(params.degree(), &mut rng);
        let pt = Plaintext::encode(&params, &v, Encoding::Simd)?;
        let mut ct: Ciphertext = sk.encrypt(&pt, &mut rng)?;

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

    #[test]
    fn reconstructed_parts_discard_seed() -> Result<(), Box<dyn StdError>> {
        let params = Parameters::test_parameters(2, 16);
        let mut rng = rng();
        let sk = SecretKey::generate(&params, &mut rng);
        let pt = Plaintext::encode(&params, &[3u64], Encoding::Polynomial)?;
        let original: Ciphertext = sk.encrypt(&pt, &mut rng)?;
        assert!(original.seed.is_some());
        for index in [0, 1] {
            let mut parts = original.clone().into_components();
            parts[index] = -&parts[index];
            let ct = Ciphertext::from_components(parts, &params)?;
            assert!(ct.seed.is_none());
            let restored = Ciphertext::from_bytes(&ct.to_bytes(), &params)?;
            assert_eq!(ct, restored);
            assert_eq!(sk.decrypt(&ct)?, sk.decrypt(&restored)?);
        }
        let mut ct = original;
        ct.iter_mut().for_each(|poly| *poly = -&*poly);
        assert!(ct.seed.is_none());
        assert_eq!(ct, Ciphertext::from_bytes(&ct.to_bytes(), &params)?);
        Ok(())
    }
}
