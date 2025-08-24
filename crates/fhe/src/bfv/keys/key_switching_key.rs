//! Key-switching keys for the BFV encryption scheme

use crate::bfv::{traits::TryConvertFrom as BfvTryConvertFrom, BfvParameters, SecretKey};
use crate::proto::bfv::KeySwitchingKey as KeySwitchingKeyProto;
use crate::{Error, Result};
use fhe_math::rq::traits::TryConvertFrom;
use fhe_math::rq::Context;
use fhe_math::{
    rns::RnsContext,
    rq::{Poly, Representation},
};
use fhe_traits::{DeserializeWithContext, Serialize};
use itertools::{izip, Itertools};
use num_bigint::BigUint;
use rand::{CryptoRng, Rng, RngCore, SeedableRng};
use rand_chacha::ChaCha8Rng;
use std::sync::Arc;
use zeroize::{Zeroize, Zeroizing};

/// Key switching key for the BFV encryption scheme.
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct KeySwitchingKey {
    /// The parameters of the underlying BFV encryption scheme.
    pub(crate) par: Arc<BfvParameters>,

    /// The (optional) seed that generated the polynomials c1.
    pub(crate) seed: Option<<ChaCha8Rng as SeedableRng>::Seed>,

    /// The key switching elements c0.
    pub(crate) c0: Box<[Poly]>,

    /// The key switching elements c1.
    pub(crate) c1: Box<[Poly]>,

    /// The level and context of the polynomials that will be key switched.
    pub(crate) ciphertext_level: usize,
    pub(crate) ctx_ciphertext: Arc<Context>,

    /// The level and context of the key switching key.
    pub(crate) ksk_level: usize,
    pub(crate) ctx_ksk: Arc<Context>,

    // For level with only one modulus, we will use basis
    pub(crate) log_base: usize,
}

impl KeySwitchingKey {
    /// Generate a [`KeySwitchingKey`] to this [`SecretKey`] from a polynomial
    /// `from`.
    pub fn new<R: RngCore + CryptoRng>(
        sk: &SecretKey,
        from: &Poly,
        ciphertext_level: usize,
        ksk_level: usize,
        rng: &mut R,
    ) -> Result<Self> {
        let ctx_ksk = sk.par.context_at_level(ksk_level)?;
        let ctx_ciphertext = sk.par.context_at_level(ciphertext_level)?;

        if from.ctx() != ctx_ksk {
            return Err(Error::DefaultError(
                "Incorrect context for polynomial from".to_string(),
            ));
        }

        let mut seed = <ChaCha8Rng as SeedableRng>::Seed::default();
        rng.fill(&mut seed);

        if ctx_ksk.moduli().len() == 1 {
            let modulus = ctx_ksk.moduli().first().unwrap();
            let log_modulus = modulus.next_power_of_two().ilog2() as usize;
            let log_base = log_modulus / 2;

            let c1 = Self::generate_c1(ctx_ksk, seed, log_modulus.div_ceil(log_base));
            let c0 = Self::generate_c0_decomposition(sk, from, &c1, rng, log_base)?;

            Ok(Self {
                par: sk.par.clone(),
                seed: Some(seed),
                c0: c0.into_boxed_slice(),
                c1: c1.into_boxed_slice(),
                ciphertext_level,
                ctx_ciphertext: ctx_ciphertext.clone(),
                ksk_level,
                ctx_ksk: ctx_ksk.clone(),
                log_base,
            })
        } else {
            let c1 = Self::generate_c1(ctx_ksk, seed, ctx_ciphertext.moduli().len());
            let c0 = Self::generate_c0(sk, from, &c1, rng)?;

            Ok(Self {
                par: sk.par.clone(),
                seed: Some(seed),
                c0: c0.into_boxed_slice(),
                c1: c1.into_boxed_slice(),
                ciphertext_level,
                ctx_ciphertext: ctx_ciphertext.clone(),
                ksk_level,
                ctx_ksk: ctx_ksk.clone(),
                log_base: 0,
            })
        }
    }

    /// Generate the c1's from the seed
    fn generate_c1(
        ctx: &Arc<Context>,
        seed: <ChaCha8Rng as SeedableRng>::Seed,
        size: usize,
    ) -> Vec<Poly> {
        let mut c1 = Vec::with_capacity(size);
        let mut rng = ChaCha8Rng::from_seed(seed);
        (0..size).for_each(|_| {
            let mut seed_i = <ChaCha8Rng as SeedableRng>::Seed::default();
            rng.fill(&mut seed_i);
            let mut a = Poly::random_from_seed(ctx, Representation::NttShoup, seed_i);
            unsafe { a.allow_variable_time_computations() }
            c1.push(a);
        });
        c1
    }

    /// Generate the c0's from the c1's and the secret key
    fn generate_c0<R: RngCore + CryptoRng>(
        sk: &SecretKey,
        from: &Poly,
        c1: &[Poly],
        rng: &mut R,
    ) -> Result<Vec<Poly>> {
        if c1.is_empty() {
            return Err(Error::DefaultError("Empty number of c1's".to_string()));
        }
        if from.representation() != &Representation::PowerBasis {
            return Err(Error::DefaultError(
                "Unexpected representation for from".to_string(),
            ));
        }

        let size = c1.len();

        let mut s = Zeroizing::new(Poly::try_convert_from(
            sk.coeffs.as_ref(),
            c1[0].ctx(),
            false,
            Representation::PowerBasis,
        )?);
        s.change_representation(Representation::Ntt);

        let rns = RnsContext::new(&sk.par.moduli[..size])?;
        let c0 = c1
            .iter()
            .enumerate()
            .map(|(i, c1i)| {
                let mut a_s = Zeroizing::new(c1i.clone());
                a_s.disallow_variable_time_computations();
                a_s.change_representation(Representation::Ntt);
                *a_s.as_mut() *= s.as_ref();
                a_s.change_representation(Representation::PowerBasis);

                let mut b =
                    Poly::small(a_s.ctx(), Representation::PowerBasis, sk.par.variance, rng)?;
                b -= &a_s;

                let gi = rns.get_garner(i).unwrap();
                let g_i_from = Zeroizing::new(gi * from);
                b += &g_i_from;

                // It is now safe to enable variable time computations.
                unsafe { b.allow_variable_time_computations() }
                b.change_representation(Representation::NttShoup);
                Ok(b)
            })
            .collect::<Result<Vec<Poly>>>()?;

        Ok(c0)
    }

    /// Generate the c0's from the c1's and the secret key
    fn generate_c0_decomposition<R: RngCore + CryptoRng>(
        sk: &SecretKey,
        from: &Poly,
        c1: &[Poly],
        rng: &mut R,
        log_base: usize,
    ) -> Result<Vec<Poly>> {
        if c1.is_empty() {
            return Err(Error::DefaultError("Empty number of c1's".to_string()));
        }

        if from.representation() != &Representation::PowerBasis {
            return Err(Error::DefaultError(
                "Unexpected representation for from".to_string(),
            ));
        }

        let mut s = Zeroizing::new(Poly::try_convert_from(
            sk.coeffs.as_ref(),
            c1[0].ctx(),
            false,
            Representation::PowerBasis,
        )?);
        s.change_representation(Representation::Ntt);

        let c0 = c1
            .iter()
            .enumerate()
            .map(|(i, c1i)| {
                let mut a_s = Zeroizing::new(c1i.clone());
                a_s.disallow_variable_time_computations();
                a_s.change_representation(Representation::Ntt);
                *a_s.as_mut() *= s.as_ref();
                a_s.change_representation(Representation::PowerBasis);

                let mut b =
                    Poly::small(a_s.ctx(), Representation::PowerBasis, sk.par.variance, rng)?;
                b -= &a_s;

                let power = BigUint::from(1u64 << (i * log_base));
                b += &(from * &power);

                // It is now safe to enable variable time computations.
                unsafe { b.allow_variable_time_computations() }
                b.change_representation(Representation::NttShoup);
                Ok(b)
            })
            .collect::<Result<Vec<Poly>>>()?;

        Ok(c0)
    }

    /// Key switch a polynomial.
    pub fn key_switch(&self, p: &Poly) -> Result<(Poly, Poly)> {
        if self.log_base != 0 {
            return self.key_switch_decomposition(p);
        }

        if p.ctx().as_ref() != self.ctx_ciphertext.as_ref() {
            return Err(Error::DefaultError(
                "The input polynomial does not have the correct context.".to_string(),
            ));
        }
        if p.representation() != &Representation::PowerBasis {
            return Err(Error::DefaultError("Incorrect representation".to_string()));
        }

        let mut c0 = Poly::zero(&self.ctx_ksk, Representation::Ntt);
        let mut c1 = Poly::zero(&self.ctx_ksk, Representation::Ntt);
        let p_coefficients = p.coefficients();
        for (c2_i_coefficients, c0_i, c1_i) in
            izip!(p_coefficients.outer_iter(), self.c0.iter(), self.c1.iter())
        {
            let mut c2_i = unsafe {
                Poly::create_constant_ntt_polynomial_with_lazy_coefficients_and_variable_time(
                    c2_i_coefficients.as_slice().unwrap(),
                    &self.ctx_ksk,
                )
            };
            c0 += &(&c2_i * c0_i);
            c2_i *= c1_i;
            c1 += &c2_i;
        }
        Ok((c0, c1))
    }

    /// Key switch a polynomial, writing the result in-place.
    pub fn key_switch_assign(&self, p: &Poly, c0: &mut Poly, c1: &mut Poly) -> Result<()> {
        if self.log_base != 0 {
            let (k0, k1) = self.key_switch_decomposition(p)?;
            *c0 = k0;
            *c1 = k1;
            return Ok(());
        }

        if p.ctx().as_ref() != self.ctx_ciphertext.as_ref() {
            return Err(Error::DefaultError(
                "The input polynomial does not have the correct context.".to_string(),
            ));
        }
        if p.representation() != &Representation::PowerBasis {
            return Err(Error::DefaultError("Incorrect representation".to_string()));
        }

        if c0.ctx().as_ref() != self.ctx_ksk.as_ref() || c0.representation() != &Representation::Ntt
        {
            *c0 = Poly::zero(&self.ctx_ksk, Representation::Ntt);
        } else {
            c0.zeroize();
        }

        if c1.ctx().as_ref() != self.ctx_ksk.as_ref() || c1.representation() != &Representation::Ntt
        {
            *c1 = Poly::zero(&self.ctx_ksk, Representation::Ntt);
        } else {
            c1.zeroize();
        }

        let p_coefficients = p.coefficients();
        for (c2_i_coefficients, c0_i, c1_i) in
            izip!(p_coefficients.outer_iter(), self.c0.iter(), self.c1.iter())
        {
            let mut c2_i = unsafe {
                Poly::create_constant_ntt_polynomial_with_lazy_coefficients_and_variable_time(
                    c2_i_coefficients.as_slice().unwrap(),
                    &self.ctx_ksk,
                )
            };
            *c0 += &(&c2_i * c0_i);
            c2_i *= c1_i;
            *c1 += &c2_i;
        }
        Ok(())
    }

    /// Key switch a polynomial.
    fn key_switch_decomposition(&self, p: &Poly) -> Result<(Poly, Poly)> {
        if p.ctx().as_ref() != self.ctx_ciphertext.as_ref() {
            return Err(Error::DefaultError(
                "The input polynomial does not have the correct context.".to_string(),
            ));
        }
        if p.representation() != &Representation::PowerBasis {
            return Err(Error::DefaultError("Incorrect representation".to_string()));
        }

        let log_modulus = p
            .ctx()
            .moduli()
            .first()
            .unwrap()
            .next_power_of_two()
            .ilog2() as usize;

        let mut coefficients = p.coefficients().to_slice().unwrap().to_vec();
        let mut c2i = vec![];
        let mask = (1u64 << self.log_base) - 1;
        (0..log_modulus.div_ceil(self.log_base)).for_each(|_| {
            c2i.push(coefficients.iter().map(|c| c & mask).collect_vec());
            coefficients.iter_mut().for_each(|c| *c >>= self.log_base);
        });

        let mut c0 = Poly::zero(&self.ctx_ksk, Representation::Ntt);
        let mut c1 = Poly::zero(&self.ctx_ksk, Representation::Ntt);
        for (c2_i_coefficients, c0_i, c1_i) in izip!(c2i.iter(), self.c0.iter(), self.c1.iter()) {
            let mut c2_i = unsafe {
                Poly::create_constant_ntt_polynomial_with_lazy_coefficients_and_variable_time(
                    c2_i_coefficients.as_slice(),
                    &self.ctx_ksk,
                )
            };
            c0 += &(&c2_i * c0_i);
            c2_i *= c1_i;
            c1 += &c2_i;
        }
        Ok((c0, c1))
    }
}

impl From<&KeySwitchingKey> for KeySwitchingKeyProto {
    fn from(value: &KeySwitchingKey) -> Self {
        let mut ksk = KeySwitchingKeyProto::default();
        if let Some(seed) = value.seed.as_ref() {
            ksk.seed = seed.to_vec();
        } else {
            ksk.c1.reserve_exact(value.c1.len());
            for c1 in value.c1.iter() {
                ksk.c1.push(c1.to_bytes())
            }
        }
        ksk.c0.reserve_exact(value.c0.len());
        for c0 in value.c0.iter() {
            ksk.c0.push(c0.to_bytes())
        }
        ksk.ciphertext_level = value.ciphertext_level as u32;
        ksk.ksk_level = value.ksk_level as u32;
        ksk.log_base = value.log_base as u32;
        ksk
    }
}

impl BfvTryConvertFrom<&KeySwitchingKeyProto> for KeySwitchingKey {
    fn try_convert_from(value: &KeySwitchingKeyProto, par: &Arc<BfvParameters>) -> Result<Self> {
        let ciphertext_level = value.ciphertext_level as usize;
        let ksk_level = value.ksk_level as usize;
        let ctx_ksk = par.context_at_level(ksk_level)?;
        let ctx_ciphertext = par.context_at_level(ciphertext_level)?;

        let c0_size: usize;
        let log_base = value.log_base as usize;
        if log_base != 0 {
            if ksk_level != par.max_level() || ciphertext_level != par.max_level() {
                return Err(Error::DefaultError(
                    "A decomposition size is specified but the levels are not maximal".to_string(),
                ));
            } else {
                let log_modulus: usize =
                    par.moduli().first().unwrap().next_power_of_two().ilog2() as usize;
                c0_size = log_modulus.div_ceil(log_base);
            }
        } else {
            c0_size = ctx_ciphertext.moduli().len();
        }

        if value.c0.len() != c0_size {
            return Err(Error::DefaultError(
                "Incorrect number of values in c0".to_string(),
            ));
        }

        let seed = if value.seed.is_empty() {
            if value.c1.len() != c0_size {
                return Err(Error::DefaultError(
                    "Incorrect number of values in c1".to_string(),
                ));
            }
            None
        } else {
            let unwrapped = <ChaCha8Rng as SeedableRng>::Seed::try_from(value.seed.clone());
            if unwrapped.is_err() {
                return Err(Error::DefaultError("Invalid seed".to_string()));
            }
            Some(unwrapped.unwrap())
        };

        let c1 = if let Some(seed) = seed {
            Self::generate_c1(ctx_ksk, seed, value.c0.len())
        } else {
            value
                .c1
                .iter()
                .map(|c1i| Poly::from_bytes(c1i, ctx_ksk).map_err(Error::MathError))
                .collect::<Result<Vec<Poly>>>()?
        };

        let c0 = value
            .c0
            .iter()
            .map(|c0i| Poly::from_bytes(c0i, ctx_ksk).map_err(Error::MathError))
            .collect::<Result<Vec<Poly>>>()?;

        Ok(Self {
            par: par.clone(),
            seed,
            c0: c0.into_boxed_slice(),
            c1: c1.into_boxed_slice(),
            ciphertext_level,
            ctx_ciphertext: ctx_ciphertext.clone(),
            ksk_level,
            ctx_ksk: ctx_ksk.clone(),
            log_base: value.log_base as usize,
        })
    }
}

#[cfg(test)]
mod tests {
    use crate::bfv::{
        keys::key_switching_key::KeySwitchingKey, traits::TryConvertFrom, BfvParameters, SecretKey,
    };
    use crate::proto::bfv::KeySwitchingKey as KeySwitchingKeyProto;
    use fhe_math::{
        rns::RnsContext,
        rq::{traits::TryConvertFrom as TryConvertFromPoly, Poly, Representation},
    };
    use num_bigint::BigUint;
    use rand::thread_rng;
    use std::error::Error;

    #[test]
    fn constructor() -> Result<(), Box<dyn Error>> {
        let mut rng = thread_rng();
        for params in [
            BfvParameters::default_arc(6, 16),
            BfvParameters::default_arc(3, 16),
        ] {
            let sk = SecretKey::random(&params, &mut rng);
            let ctx = params.context_at_level(0)?;
            let p = Poly::small(ctx, Representation::PowerBasis, 10, &mut rng)?;
            let ksk = KeySwitchingKey::new(&sk, &p, 0, 0, &mut rng);
            assert!(ksk.is_ok());
        }
        Ok(())
    }

    #[test]
    fn constructor_last_level() -> Result<(), Box<dyn Error>> {
        let mut rng = thread_rng();
        for params in [
            BfvParameters::default_arc(6, 16),
            BfvParameters::default_arc(3, 16),
        ] {
            let level = params.moduli().len() - 1;
            let sk = SecretKey::random(&params, &mut rng);
            let ctx = params.context_at_level(level)?;
            let p = Poly::small(ctx, Representation::PowerBasis, 10, &mut rng)?;
            let ksk = KeySwitchingKey::new(&sk, &p, level, level, &mut rng);
            assert!(ksk.is_ok());
        }
        Ok(())
    }

    #[test]
    fn key_switch() -> Result<(), Box<dyn Error>> {
        let mut rng = thread_rng();
        for params in [BfvParameters::default_arc(6, 16)] {
            for _ in 0..100 {
                let sk = SecretKey::random(&params, &mut rng);
                let ctx = params.context_at_level(0)?;
                let mut p = Poly::small(ctx, Representation::PowerBasis, 10, &mut rng)?;
                let ksk = KeySwitchingKey::new(&sk, &p, 0, 0, &mut rng)?;
                let mut s = Poly::try_convert_from(
                    sk.coeffs.as_ref(),
                    ctx,
                    false,
                    Representation::PowerBasis,
                )
                .map_err(crate::Error::MathError)?;
                s.change_representation(Representation::Ntt);

                let mut input = Poly::random(ctx, Representation::PowerBasis, &mut rng);
                let (c0, c1) = ksk.key_switch(&input)?;

                let mut c2 = &c0 + &(&c1 * &s);
                c2.change_representation(Representation::PowerBasis);

                input.change_representation(Representation::Ntt);
                p.change_representation(Representation::Ntt);
                let mut c3 = &input * &p;
                c3.change_representation(Representation::PowerBasis);

                let rns = RnsContext::new(&params.moduli)?;
                Vec::<BigUint>::from(&(&c2 - &c3)).iter().for_each(|b| {
                    assert!(std::cmp::min(b.bits(), (rns.modulus() - b).bits()) <= 70)
                });
            }
        }
        Ok(())
    }

    #[test]
    fn key_switch_assign_matches() -> Result<(), Box<dyn Error>> {
        let mut rng = thread_rng();
        {
            let params = BfvParameters::default_arc(6, 16);
            let sk = SecretKey::random(&params, &mut rng);
            let ctx = params.context_at_level(0)?;
            let p = Poly::small(ctx, Representation::PowerBasis, 10, &mut rng)?;
            let ksk = KeySwitchingKey::new(&sk, &p, 0, 0, &mut rng)?;
            let input = Poly::random(ctx, Representation::PowerBasis, &mut rng);

            let (c0, c1) = ksk.key_switch(&input)?;

            let mut a0 = Poly::zero(&ksk.ctx_ksk, Representation::Ntt);
            let mut a1 = Poly::zero(&ksk.ctx_ksk, Representation::Ntt);
            ksk.key_switch_assign(&input, &mut a0, &mut a1)?;

            assert_eq!(c0, a0);
            assert_eq!(c1, a1);
        }
        Ok(())
    }

    #[test]
    fn key_switch_decomposition() -> Result<(), Box<dyn Error>> {
        let mut rng = thread_rng();
        for params in [BfvParameters::default_arc(6, 16)] {
            for _ in 0..100 {
                let sk = SecretKey::random(&params, &mut rng);
                let ctx = params.context_at_level(5)?;
                let mut p = Poly::small(ctx, Representation::PowerBasis, 10, &mut rng)?;
                let ksk = KeySwitchingKey::new(&sk, &p, 5, 5, &mut rng)?;
                let mut s = Poly::try_convert_from(
                    sk.coeffs.as_ref(),
                    ctx,
                    false,
                    Representation::PowerBasis,
                )
                .map_err(crate::Error::MathError)?;
                s.change_representation(Representation::Ntt);

                let mut input = Poly::random(ctx, Representation::PowerBasis, &mut rng);
                let (c0, c1) = ksk.key_switch(&input)?;

                let mut c2 = &c0 + &(&c1 * &s);
                c2.change_representation(Representation::PowerBasis);

                input.change_representation(Representation::Ntt);
                p.change_representation(Representation::Ntt);
                let mut c3 = &input * &p;
                c3.change_representation(Representation::PowerBasis);

                let rns = RnsContext::new(ctx.moduli())?;
                Vec::<BigUint>::from(&(&c2 - &c3)).iter().for_each(|b| {
                    assert!(
                        std::cmp::min(b.bits(), (rns.modulus() - b).bits())
                            <= (rns.modulus().bits() / 2) + 10
                    )
                });
            }
        }
        Ok(())
    }

    #[test]
    fn proto_conversion() -> Result<(), Box<dyn Error>> {
        let mut rng = thread_rng();
        for params in [
            BfvParameters::default_arc(6, 16),
            BfvParameters::default_arc(3, 16),
        ] {
            let sk = SecretKey::random(&params, &mut rng);
            let ctx = params.context_at_level(0)?;
            let p = Poly::small(ctx, Representation::PowerBasis, 10, &mut rng)?;
            let ksk = KeySwitchingKey::new(&sk, &p, 0, 0, &mut rng)?;
            let ksk_proto = KeySwitchingKeyProto::from(&ksk);
            assert_eq!(ksk, KeySwitchingKey::try_convert_from(&ksk_proto, &params)?);
        }
        Ok(())
    }
}
