#![warn(missing_docs, unused_imports)]
// Allow indexing/slicing in this performance-critical polynomial arithmetic module.
// This code heavily uses RNS representation and requires fast array access.
#![allow(clippy::indexing_slicing)]

//! Polynomials in R_q\[x\] = (ZZ_q1 x ... x ZZ_qn)\[x\] where the qi's are
//! prime moduli in zq.

mod context;
mod convert;
mod ops;
mod serialize;

pub mod scaler;
pub mod switcher;
pub mod traits;
use self::{scaler::Scaler, switcher::Switcher, traits::TryConvertFrom};
use crate::{zq::Modulus, Error, Result};
pub use context::Context;
use fhe_util::sample_vec_cbd;
use itertools::{izip, Itertools};
use ndarray::{s, Array2, ArrayView2, Axis};
pub use ops::dot_product;
use rand::{CryptoRng, RngCore, SeedableRng};
use rand_chacha::ChaCha8Rng;
use sha2::{Digest, Sha256};
use std::sync::Arc;
use zeroize::{Zeroize, Zeroizing};

/// Possible representations of the underlying polynomial.
#[derive(Default, Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Representation {
    /// This is the list of coefficients ci, such that the polynomial is c0 + c1
    /// * x + ... + c_(degree - 1) * x^(degree - 1)
    #[default]
    PowerBasis,
    /// This is the NTT representation of the PowerBasis representation.
    Ntt,
    /// This is a "Shoup" representation of the Ntt representation used for
    /// faster multiplication.
    NttShoup,
}

/// An exponent for a substitution.
#[derive(Debug, PartialEq, Eq)]
pub struct SubstitutionExponent {
    /// The value of the exponent.
    pub exponent: usize,

    ctx: Arc<Context>,
    power_bitrev: Vec<usize>,
}

impl SubstitutionExponent {
    /// Creates a substitution element from an exponent.
    /// Returns an error if the exponent is even modulo 2 * degree.
    pub fn new(ctx: &Arc<Context>, exponent: usize) -> Result<Self> {
        let exponent = exponent % (2 * ctx.degree);
        if exponent & 1 == 0 {
            return Err(Error::Default(
                "The exponent should be odd modulo 2 * degree".to_string(),
            ));
        }
        let mut power = (exponent - 1) / 2;
        let mask = ctx.degree - 1;
        let power_bitrev = (0..ctx.degree)
            .map(|_| {
                let r = (power & mask).reverse_bits() >> (ctx.degree.leading_zeros() + 1);
                power += exponent;
                r
            })
            .collect_vec();
        Ok(Self {
            ctx: ctx.clone(),
            exponent,
            power_bitrev,
        })
    }
}

/// Struct that holds a polynomial for a specific context.
#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub struct Poly {
    ctx: Arc<Context>,
    representation: Representation,
    has_lazy_coefficients: bool,
    allow_variable_time_computations: bool,
    coefficients: Array2<u64>,
    coefficients_shoup: Option<Array2<u64>>,
}

// Implements zeroization of polynomials
impl Zeroize for Poly {
    fn zeroize(&mut self) {
        if let Some(coeffs) = self.coefficients.as_slice_mut() {
            coeffs.zeroize()
        }
        self.zeroize_shoup()
    }
}

impl AsRef<Poly> for Poly {
    fn as_ref(&self) -> &Poly {
        self
    }
}

impl AsMut<Poly> for Poly {
    fn as_mut(&mut self) -> &mut Poly {
        self
    }
}

impl Poly {
    /// Creates a polynomial holding the constant 0.
    #[must_use]
    pub fn zero(ctx: &Arc<Context>, representation: Representation) -> Self {
        Self {
            ctx: ctx.clone(),
            representation,
            allow_variable_time_computations: false,
            has_lazy_coefficients: false,
            coefficients: Array2::zeros((ctx.q.len(), ctx.degree)),
            coefficients_shoup: if representation == Representation::NttShoup {
                Some(Array2::zeros((ctx.q.len(), ctx.degree)))
            } else {
                None
            },
        }
    }

    /// Enable variable time computations when this polynomial is involved.
    ///
    /// # Safety
    ///
    /// By default, this is marked as unsafe, but is usually safe when only
    /// public data is processed.
    pub unsafe fn allow_variable_time_computations(&mut self) {
        self.allow_variable_time_computations = true
    }

    /// Disable variable time computations when this polynomial is involved.
    pub fn disallow_variable_time_computations(&mut self) {
        self.allow_variable_time_computations = false
    }

    /// Current representation of the polynomial.
    #[must_use]
    pub const fn representation(&self) -> &Representation {
        &self.representation
    }

    /// Zeroize the shoup coefficients
    fn zeroize_shoup(&mut self) {
        if let Some(coeffs_shoup) = self
            .coefficients_shoup
            .as_mut()
            .and_then(|f| f.as_slice_mut())
        {
            coeffs_shoup.zeroize()
        }
    }

    /// Change the representation of the underlying polynomial.
    pub fn change_representation(&mut self, to: Representation) {
        if self.representation == to {
            return;
        }

        match (&self.representation, &to) {
            (Representation::PowerBasis, Representation::Ntt) => self.ntt_forward(),
            (Representation::PowerBasis, Representation::NttShoup) => {
                self.ntt_forward();
                self.compute_coefficients_shoup()
            }
            (Representation::Ntt, Representation::PowerBasis) => self.ntt_backward(),
            (Representation::Ntt, Representation::NttShoup) => self.compute_coefficients_shoup(),
            (Representation::NttShoup, Representation::PowerBasis) => {
                self.zeroize_shoup();
                self.coefficients_shoup = None;
                self.ntt_backward()
            }
            (Representation::NttShoup, Representation::Ntt) => {
                self.zeroize_shoup();
                self.coefficients_shoup = None;
            }
            _ => unreachable!(),
        }

        self.representation = to;
    }

    /// Compute the Shoup representation of the coefficients.
    fn compute_coefficients_shoup(&mut self) {
        let mut coefficients_shoup = Array2::zeros((self.ctx.q.len(), self.ctx.degree));
        izip!(
            coefficients_shoup.outer_iter_mut(),
            self.coefficients.outer_iter(),
            self.ctx.q.iter()
        )
        .for_each(|(mut v_shoup, v, qi)| {
            v_shoup
                .as_slice_mut()
                .unwrap()
                .copy_from_slice(&qi.shoup_vec(v.as_slice().unwrap()))
        });
        self.coefficients_shoup = Some(coefficients_shoup)
    }

    /// Override the internal representation to a given representation.
    ///
    /// # Safety
    ///
    /// Prefer the `change_representation` function to safely modify the
    /// polynomial representation. If the `to` representation is NttShoup, the
    /// coefficients are still computed correctly to avoid being in an unstable
    /// state. If we override a polynomial with Shoup coefficients, we zeroize
    /// them.
    pub unsafe fn override_representation(&mut self, to: Representation) {
        if self.coefficients_shoup.is_some() {
            self.zeroize_shoup();
            self.coefficients_shoup = None
        }
        if to == Representation::NttShoup {
            self.compute_coefficients_shoup()
        }
        self.representation = to;
    }

    /// Generate a random polynomial.
    pub fn random<R: RngCore + CryptoRng>(
        ctx: &Arc<Context>,
        representation: Representation,
        rng: &mut R,
    ) -> Self {
        let mut p = Poly::zero(ctx, representation);
        izip!(p.coefficients.outer_iter_mut(), ctx.q.iter()).for_each(|(mut v, qi)| {
            v.as_slice_mut()
                .unwrap()
                .copy_from_slice(&qi.random_vec(ctx.degree, rng))
        });
        if p.representation == Representation::NttShoup {
            p.compute_coefficients_shoup()
        }
        p
    }

    /// Generate a random polynomial deterministically from a seed.
    #[must_use]
    pub fn random_from_seed(
        ctx: &Arc<Context>,
        representation: Representation,
        seed: <ChaCha8Rng as SeedableRng>::Seed,
    ) -> Self {
        // Let's hash the seed into a ChaCha8Rng seed.
        let mut hasher = Sha256::new();
        hasher.update(seed);
        let mut prng =
            ChaCha8Rng::from_seed(<ChaCha8Rng as SeedableRng>::Seed::from(hasher.finalize()));
        let mut p = Poly::zero(ctx, representation);
        izip!(p.coefficients.outer_iter_mut(), ctx.q.iter()).for_each(|(mut v, qi)| {
            v.as_slice_mut()
                .unwrap()
                .copy_from_slice(&qi.random_vec(ctx.degree, &mut prng))
        });
        if p.representation == Representation::NttShoup {
            p.compute_coefficients_shoup()
        }
        p
    }

    /// Generate a small polynomial and convert into the specified
    /// representation.
    ///
    /// Returns an error if the variance does not belong to [1, ..., 16].
    pub fn small<T: RngCore + CryptoRng>(
        ctx: &Arc<Context>,
        representation: Representation,
        variance: usize,
        rng: &mut T,
    ) -> Result<Self> {
        if !(1..=16).contains(&variance) {
            return Err(Error::Default(
                "The variance should be an integer between 1 and 16".to_string(),
            ));
        }

        let coeffs = Zeroizing::new(
            sample_vec_cbd(ctx.degree, variance, rng).map_err(|e| Error::Default(e.to_string()))?,
        );
        let mut p = Poly::try_convert_from(
            coeffs.as_ref() as &[i64],
            ctx,
            false,
            Representation::PowerBasis,
        )?;
        if representation != Representation::PowerBasis {
            p.change_representation(representation);
        }
        Ok(p)
    }

    /// Access the polynomial coefficients in RNS representation.
    #[must_use]
    pub fn coefficients(&self) -> ArrayView2<'_, u64> {
        self.coefficients.view()
    }

    /// Computes the forward Ntt on the coefficients
    fn ntt_forward(&mut self) {
        if self.allow_variable_time_computations {
            izip!(self.coefficients.outer_iter_mut(), self.ctx.ops.iter())
                .for_each(|(mut v, op)| unsafe { op.forward_vt(v.as_mut_ptr()) });
        } else {
            izip!(self.coefficients.outer_iter_mut(), self.ctx.ops.iter())
                .for_each(|(mut v, op)| op.forward(v.as_slice_mut().unwrap()));
        }
    }

    /// Computes the backward Ntt on the coefficients
    fn ntt_backward(&mut self) {
        if self.allow_variable_time_computations {
            izip!(self.coefficients.outer_iter_mut(), self.ctx.ops.iter())
                .for_each(|(mut v, op)| unsafe { op.backward_vt(v.as_mut_ptr()) });
        } else {
            izip!(self.coefficients.outer_iter_mut(), self.ctx.ops.iter())
                .for_each(|(mut v, op)| op.backward(v.as_slice_mut().unwrap()));
        }
    }

    /// Substitute x by x^i in a polynomial.
    /// In PowerBasis representation, i can be any integer that is not a
    /// multiple of 2 * degree. In Ntt and NttShoup representation, i can be any
    /// odd integer that is not a multiple of 2 * degree.
    pub fn substitute(&self, i: &SubstitutionExponent) -> Result<Poly> {
        let mut q = Poly::zero(&self.ctx, self.representation);
        if self.allow_variable_time_computations {
            unsafe { q.allow_variable_time_computations() }
        }
        match self.representation {
            Representation::Ntt | Representation::NttShoup => {
                izip!(
                    q.coefficients.outer_iter_mut(),
                    self.coefficients.outer_iter()
                )
                .for_each(|(mut q_row, p_row)| {
                    for (j, k) in izip!(self.ctx.bitrev.iter(), i.power_bitrev.iter()) {
                        q_row[*j] = p_row[*k]
                    }
                });
                if self.representation == Representation::NttShoup {
                    izip!(
                        q.coefficients_shoup.as_mut().unwrap().outer_iter_mut(),
                        self.coefficients_shoup.as_ref().unwrap().outer_iter()
                    )
                    .for_each(|(mut q_row, p_row)| {
                        for (j, k) in izip!(self.ctx.bitrev.iter(), i.power_bitrev.iter()) {
                            q_row[*j] = p_row[*k]
                        }
                    });
                }
            }
            Representation::PowerBasis => {
                let mut power = 0usize;
                let mask = self.ctx.degree - 1;
                for j in 0..self.ctx.degree {
                    izip!(
                        self.ctx.q.iter(),
                        q.coefficients.slice_mut(s![.., power & mask]),
                        self.coefficients.slice(s![.., j])
                    )
                    .for_each(|(qi, qij, pij)| {
                        if power & self.ctx.degree != 0 {
                            *qij = qi.sub(*qij, *pij)
                        } else {
                            *qij = qi.add(*qij, *pij)
                        }
                    });
                    power += i.exponent
                }
            }
        }

        Ok(q)
    }

    /// Create a polynomial which can only be multiplied by a polynomial in
    /// NttShoup representation. All other operations may panic.
    ///
    /// # Safety
    /// This operation also creates a polynomial that allows variable time
    /// operations.
    #[must_use]
    pub unsafe fn create_constant_ntt_polynomial_with_lazy_coefficients_and_variable_time(
        power_basis_coefficients: &[u64],
        ctx: &Arc<Context>,
    ) -> Self {
        let mut coefficients = Array2::zeros((ctx.q.len(), ctx.degree));
        izip!(coefficients.outer_iter_mut(), ctx.q.iter(), ctx.ops.iter()).for_each(
            |(mut p, qi, op)| {
                p.as_slice_mut()
                    .unwrap()
                    .clone_from_slice(power_basis_coefficients);
                qi.lazy_reduce_vec(p.as_slice_mut().unwrap());
                op.forward_vt_lazy(p.as_mut_ptr());
            },
        );
        Self {
            ctx: ctx.clone(),
            representation: Representation::Ntt,
            allow_variable_time_computations: true,
            coefficients,
            coefficients_shoup: None,
            has_lazy_coefficients: true,
        }
    }

    /// Modulus switch down the polynomial by dividing and rounding each
    /// coefficient by the last modulus in the chain, then drops the last
    /// modulus, as described in Algorithm 2 of <https://eprint.iacr.org/2018/931.pdf>.
    ///
    /// Returns an error if there is no next context or if the representation
    /// is not PowerBasis.
    pub fn switch_down(&mut self) -> Result<()> {
        if self.ctx.next_context.is_none() {
            return Err(Error::NoMoreContext);
        }

        if self.representation != Representation::PowerBasis {
            return Err(Error::IncorrectRepresentation(
                self.representation,
                Representation::PowerBasis,
            ));
        }

        // Unwrap the next_context.
        let next_context = self.ctx.next_context.as_ref().unwrap();

        let q_len = self.ctx.q.len();
        let q_last = self.ctx.q.last().unwrap();
        let q_last_div_2 = (**q_last) / 2;

        // Add (q_last - 1) / 2 to change from flooring to rounding
        let (mut q_new_polys, mut q_last_poly) =
            self.coefficients.view_mut().split_at(Axis(0), q_len - 1);

        let add: fn(&Modulus, u64, u64) -> u64 = if self.allow_variable_time_computations {
            |qi, a, b| unsafe { qi.add_vt(a, b) }
        } else {
            |qi, a, b| qi.add(a, b)
        };
        let reduce: unsafe fn(&Modulus, u64) -> u64 = if self.allow_variable_time_computations {
            |qi, a| unsafe { qi.reduce_vt(a) }
        } else {
            |qi, a| qi.reduce(a)
        };

        q_last_poly
            .iter_mut()
            .for_each(|coeff| *coeff = add(q_last, *coeff, q_last_div_2));
        izip!(
            q_new_polys.outer_iter_mut(),
            self.ctx.q.iter(),
            self.ctx.inv_last_qi_mod_qj.iter(),
            self.ctx.inv_last_qi_mod_qj_shoup.iter(),
        )
        .for_each(|(coeffs, qi, inv, inv_shoup)| {
            let q_last_div_2_mod_qi = **qi - unsafe { reduce(qi, q_last_div_2) }; // Up to qi.modulus()
            for (coeff, q_last_coeff) in izip!(coeffs, q_last_poly.iter()) {
                // (x mod q_last - q_L/2) mod q_i
                let tmp = qi.lazy_reduce(*q_last_coeff) + q_last_div_2_mod_qi; // Up to 3 * qi.modulus()

                // ((x mod q_i) - (x mod q_last) + (q_L/2 mod q_i)) mod q_i
                // = (x - x mod q_last + q_L/2) mod q_i
                *coeff += 3 * (**qi) - tmp; // Up to 4 * qi.modulus()

                // q_last^{-1} * (x - x mod q_last) mod q_i
                *coeff = qi.mul_shoup(*coeff, *inv, *inv_shoup);
            }
        });

        // Remove the last row, and update the context.
        if !self.allow_variable_time_computations {
            q_last_poly.as_slice_mut().unwrap().zeroize();
        }
        self.coefficients.remove_index(Axis(0), q_len - 1);
        self.ctx = next_context.clone();

        Ok(())
    }

    /// Modulo switch down to a smaller context.
    ///
    /// Returns an error if there is the provided context is not a child of the
    /// current context, or if the polynomial is not in PowerBasis
    /// representation.
    pub fn switch_down_to(&mut self, context: &Arc<Context>) -> Result<()> {
        let niterations = self.ctx.niterations_to(context)?;
        for _ in 0..niterations {
            self.switch_down()?;
        }
        assert_eq!(&self.ctx, context);
        Ok(())
    }

    /// Modulo switch to another context. The target context needs not to be
    /// related to the current context.
    pub fn switch(&self, switcher: &Switcher) -> Result<Poly> {
        switcher.switch(self)
    }

    /// Scale a polynomial using a scaler.
    pub fn scale(&self, scaler: &Scaler) -> Result<Poly> {
        scaler.scale(self)
    }

    /// Returns the context of the underlying polynomial
    #[must_use]
    pub fn ctx(&self) -> &Arc<Context> {
        &self.ctx
    }

    /// Multiplies a polynomial in PowerBasis representation by x^(-power).
    pub fn multiply_inverse_power_of_x(&mut self, power: usize) -> Result<()> {
        if self.representation != Representation::PowerBasis {
            return Err(Error::IncorrectRepresentation(
                self.representation,
                Representation::PowerBasis,
            ));
        }

        let shift = ((self.ctx.degree << 1) - power) % (self.ctx.degree << 1);
        let mask = self.ctx.degree - 1;
        let mut new_coefficients = Array2::zeros((self.ctx.q.len(), self.ctx.degree));
        izip!(
            new_coefficients.outer_iter_mut(),
            self.coefficients.outer_iter(),
            self.ctx.q.iter()
        )
        .for_each(|(mut new_coeffs, orig_coeffs, qi)| {
            for k in 0..self.ctx.degree {
                let index = shift + k;
                if index & self.ctx.degree == 0 {
                    new_coeffs[index & mask] = orig_coeffs[k];
                } else {
                    new_coeffs[index & mask] = qi.neg(orig_coeffs[k]);
                }
            }
        });
        self.coefficients = new_coefficients;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::{switcher::Switcher, Context, Poly, Representation};
    use crate::{rq::SubstitutionExponent, zq::Modulus};
    use fhe_util::variance;
    use itertools::Itertools;
    use num_bigint::BigUint;
    use num_traits::{One, Zero};
    use rand::{Rng, SeedableRng};
    use rand_chacha::ChaCha8Rng;
    use std::{error::Error, sync::Arc};

    // Moduli to be used in tests.
    const MODULI: &[u64; 5] = &[
        1153,
        4611686018326724609,
        4611686018309947393,
        4611686018232352769,
        4611686018171535361,
    ];

    #[test]
    fn poly_zero() -> Result<(), Box<dyn Error>> {
        let reference = &[
            BigUint::zero(),
            BigUint::zero(),
            BigUint::zero(),
            BigUint::zero(),
            BigUint::zero(),
            BigUint::zero(),
            BigUint::zero(),
            BigUint::zero(),
            BigUint::zero(),
            BigUint::zero(),
            BigUint::zero(),
            BigUint::zero(),
            BigUint::zero(),
            BigUint::zero(),
            BigUint::zero(),
            BigUint::zero(),
        ];

        for modulus in MODULI {
            let ctx = Arc::new(Context::new(&[*modulus], 16)?);
            let p = Poly::zero(&ctx, Representation::PowerBasis);
            let q = Poly::zero(&ctx, Representation::Ntt);
            assert_ne!(p, q);
            assert_eq!(Vec::<u64>::try_from(&p).unwrap(), &[0; 16]);
            assert_eq!(Vec::<u64>::try_from(&q).unwrap(), &[0; 16]);
        }

        let ctx = Arc::new(Context::new(MODULI, 16)?);
        let p = Poly::zero(&ctx, Representation::PowerBasis);
        let q = Poly::zero(&ctx, Representation::Ntt);
        assert_ne!(p, q);
        assert_eq!(Vec::<u64>::try_from(&p).unwrap(), [0; 16 * MODULI.len()]);
        assert_eq!(Vec::<u64>::try_from(&q).unwrap(), [0; 16 * MODULI.len()]);
        assert_eq!(Vec::<BigUint>::from(&p), reference);
        assert_eq!(Vec::<BigUint>::from(&q), reference);

        Ok(())
    }

    #[test]
    fn ctx() -> Result<(), Box<dyn Error>> {
        for modulus in MODULI {
            let ctx = Arc::new(Context::new(&[*modulus], 16)?);
            let p = Poly::zero(&ctx, Representation::PowerBasis);
            assert_eq!(p.ctx(), &ctx);
        }

        let ctx = Arc::new(Context::new(MODULI, 16)?);
        let p = Poly::zero(&ctx, Representation::PowerBasis);
        assert_eq!(p.ctx(), &ctx);

        Ok(())
    }

    #[test]
    fn random() -> Result<(), Box<dyn Error>> {
        let mut rng = rand::rng();
        for _ in 0..100 {
            let mut seed = <ChaCha8Rng as SeedableRng>::Seed::default();
            rand::rng().fill(&mut seed);

            for modulus in MODULI {
                let ctx = Arc::new(Context::new(&[*modulus], 16)?);
                let p = Poly::random_from_seed(&ctx, Representation::Ntt, seed);
                let q = Poly::random_from_seed(&ctx, Representation::Ntt, seed);
                assert_eq!(p, q);
            }

            let ctx = Arc::new(Context::new(MODULI, 16)?);
            let p = Poly::random_from_seed(&ctx, Representation::Ntt, seed);
            let q = Poly::random_from_seed(&ctx, Representation::Ntt, seed);
            assert_eq!(p, q);

            rand::rng().fill(&mut seed);
            let p = Poly::random_from_seed(&ctx, Representation::Ntt, seed);
            assert_ne!(p, q);

            let r = Poly::random(&ctx, Representation::Ntt, &mut rng);
            assert_ne!(p, r);
            assert_ne!(q, r);
        }
        Ok(())
    }

    #[test]
    fn coefficients() -> Result<(), Box<dyn Error>> {
        let mut rng = rand::rng();
        for _ in 0..50 {
            for modulus in MODULI {
                let ctx = Arc::new(Context::new(&[*modulus], 16)?);
                let p = Poly::random(&ctx, Representation::Ntt, &mut rng);
                let p_coefficients = Vec::<u64>::try_from(&p).unwrap();
                assert_eq!(p_coefficients, p.coefficients().as_slice().unwrap())
            }

            let ctx = Arc::new(Context::new(MODULI, 16)?);
            let p = Poly::random(&ctx, Representation::Ntt, &mut rng);
            let p_coefficients = Vec::<u64>::try_from(&p).unwrap();
            assert_eq!(p_coefficients, p.coefficients().as_slice().unwrap())
        }
        Ok(())
    }

    #[test]
    fn modulus() -> Result<(), Box<dyn Error>> {
        for modulus in MODULI {
            let modulus_biguint = BigUint::from(*modulus);
            let ctx = Arc::new(Context::new(&[*modulus], 16)?);
            assert_eq!(ctx.modulus(), &modulus_biguint)
        }

        let mut modulus_biguint = BigUint::one();
        MODULI.iter().for_each(|m| modulus_biguint *= *m);
        let ctx = Arc::new(Context::new(MODULI, 16)?);
        assert_eq!(ctx.modulus(), &modulus_biguint);

        Ok(())
    }

    #[test]
    fn allow_variable_time_computations() -> Result<(), Box<dyn Error>> {
        let mut rng = rand::rng();
        for modulus in MODULI {
            let ctx = Arc::new(Context::new(&[*modulus], 16)?);
            let mut p = Poly::random(&ctx, Representation::default(), &mut rng);
            assert!(!p.allow_variable_time_computations);

            unsafe { p.allow_variable_time_computations() }
            assert!(p.allow_variable_time_computations);

            let q = p.clone();
            assert!(q.allow_variable_time_computations);

            p.disallow_variable_time_computations();
            assert!(!p.allow_variable_time_computations);
        }

        let ctx = Arc::new(Context::new(MODULI, 16)?);
        let mut p = Poly::random(&ctx, Representation::default(), &mut rng);
        assert!(!p.allow_variable_time_computations);

        unsafe { p.allow_variable_time_computations() }
        assert!(p.allow_variable_time_computations);

        let q = p.clone();
        assert!(q.allow_variable_time_computations);

        // Allowing variable time propagates.
        let mut p = Poly::random(&ctx, Representation::Ntt, &mut rng);
        unsafe { p.allow_variable_time_computations() }
        let mut q = Poly::random(&ctx, Representation::Ntt, &mut rng);

        assert!(!q.allow_variable_time_computations);
        q *= &p;
        assert!(q.allow_variable_time_computations);

        q.disallow_variable_time_computations();
        q += &p;
        assert!(q.allow_variable_time_computations);

        q.disallow_variable_time_computations();
        q -= &p;
        assert!(q.allow_variable_time_computations);

        q = -&p;
        assert!(q.allow_variable_time_computations);

        Ok(())
    }

    #[test]
    fn change_representation() -> Result<(), Box<dyn Error>> {
        let mut rng = rand::rng();
        let ctx = Arc::new(Context::new(MODULI, 16)?);

        let mut p = Poly::random(&ctx, Representation::default(), &mut rng);
        assert_eq!(p.representation, Representation::default());
        assert_eq!(p.representation(), &Representation::default());

        p.change_representation(Representation::PowerBasis);
        assert_eq!(p.representation, Representation::PowerBasis);
        assert_eq!(p.representation(), &Representation::PowerBasis);
        assert!(p.coefficients_shoup.is_none());
        let q = p.clone();

        p.change_representation(Representation::Ntt);
        assert_eq!(p.representation, Representation::Ntt);
        assert_eq!(p.representation(), &Representation::Ntt);
        assert_ne!(p.coefficients, q.coefficients);
        assert!(p.coefficients_shoup.is_none());
        let q_ntt = p.clone();

        p.change_representation(Representation::NttShoup);
        assert_eq!(p.representation, Representation::NttShoup);
        assert_eq!(p.representation(), &Representation::NttShoup);
        assert_ne!(p.coefficients, q.coefficients);
        assert!(p.coefficients_shoup.is_some());
        let q_ntt_shoup = p.clone();

        p.change_representation(Representation::PowerBasis);
        assert_eq!(p, q);

        p.change_representation(Representation::NttShoup);
        assert_eq!(p, q_ntt_shoup);

        p.change_representation(Representation::Ntt);
        assert_eq!(p, q_ntt);

        p.change_representation(Representation::PowerBasis);
        assert_eq!(p, q);

        Ok(())
    }

    #[test]
    fn override_representation() -> Result<(), Box<dyn Error>> {
        let mut rng = rand::rng();
        let ctx = Arc::new(Context::new(MODULI, 16)?);

        let mut p = Poly::random(&ctx, Representation::PowerBasis, &mut rng);
        assert_eq!(p.representation(), &p.representation);
        let q = p.clone();

        unsafe { p.override_representation(Representation::Ntt) }
        assert_eq!(p.representation, Representation::Ntt);
        assert_eq!(p.representation(), &p.representation);
        assert_eq!(p.coefficients, q.coefficients);
        assert!(p.coefficients_shoup.is_none());

        unsafe { p.override_representation(Representation::NttShoup) }
        assert_eq!(p.representation, Representation::NttShoup);
        assert_eq!(p.representation(), &p.representation);
        assert_eq!(p.coefficients, q.coefficients);
        assert!(p.coefficients_shoup.is_some());

        unsafe { p.override_representation(Representation::PowerBasis) }
        assert_eq!(p, q);

        unsafe { p.override_representation(Representation::NttShoup) }
        assert!(p.coefficients_shoup.is_some());

        unsafe { p.override_representation(Representation::Ntt) }
        assert!(p.coefficients_shoup.is_none());

        Ok(())
    }

    #[test]
    fn small() -> Result<(), Box<dyn Error>> {
        let mut rng = rand::rng();
        for modulus in MODULI {
            let ctx = Arc::new(Context::new(&[*modulus], 16)?);
            let q = Modulus::new(*modulus).unwrap();

            let e = Poly::small(&ctx, Representation::PowerBasis, 0, &mut rng);
            assert!(e.is_err());
            assert_eq!(
                e.unwrap_err().to_string(),
                "The variance should be an integer between 1 and 16"
            );
            let e = Poly::small(&ctx, Representation::PowerBasis, 17, &mut rng);
            assert!(e.is_err());
            assert_eq!(
                e.unwrap_err().to_string(),
                "The variance should be an integer between 1 and 16"
            );

            for i in 1..=16 {
                let p = Poly::small(&ctx, Representation::PowerBasis, i, &mut rng)?;
                let coefficients = p.coefficients().to_slice().unwrap();
                let v = unsafe { q.center_vec_vt(coefficients) };

                assert!(v.iter().map(|vi| vi.abs()).max().unwrap() <= 2 * i as i64);
            }
        }

        // Generate a very large polynomial to check the variance (here equal to 8).
        let ctx = Arc::new(Context::new(&[4611686018326724609], 1 << 18)?);
        let q = Modulus::new(4611686018326724609).unwrap();
        let mut rng = rand::rng();
        let p = Poly::small(&ctx, Representation::PowerBasis, 16, &mut rng)?;
        let coefficients = p.coefficients().to_slice().unwrap();
        let v = unsafe { q.center_vec_vt(coefficients) };
        assert!(v.iter().map(|vi| vi.abs()).max().unwrap() <= 32);
        assert_eq!(variance(&v).round(), 16.0);

        Ok(())
    }

    #[test]
    fn substitute() -> Result<(), Box<dyn Error>> {
        let mut rng = rand::rng();
        for modulus in MODULI {
            let ctx = Arc::new(Context::new(&[*modulus], 16)?);
            let p = Poly::random(&ctx, Representation::PowerBasis, &mut rng);
            let mut p_ntt = p.clone();
            p_ntt.change_representation(Representation::Ntt);
            let mut p_ntt_shoup = p.clone();
            p_ntt_shoup.change_representation(Representation::NttShoup);
            let p_coeffs = Vec::<u64>::try_from(&p).unwrap();

            // Substitution by a multiple of 2 * degree, or even numbers, should fail
            assert!(SubstitutionExponent::new(&ctx, 0).is_err());
            assert!(SubstitutionExponent::new(&ctx, 2).is_err());
            assert!(SubstitutionExponent::new(&ctx, 16).is_err());

            // Substitution by 1 leaves the polynomials unchanged
            assert_eq!(p, p.substitute(&SubstitutionExponent::new(&ctx, 1)?)?);
            assert_eq!(
                p_ntt,
                p_ntt.substitute(&SubstitutionExponent::new(&ctx, 1)?)?
            );
            assert_eq!(
                p_ntt_shoup,
                p_ntt_shoup.substitute(&SubstitutionExponent::new(&ctx, 1)?)?
            );

            // Substitution by 3
            let mut q = p.substitute(&SubstitutionExponent::new(&ctx, 3)?)?;
            let mut v = vec![0u64; 16];
            for i in 0..16 {
                v[(3 * i) % 16] = if ((3 * i) / 16) & 1 == 1 && p_coeffs[i] > 0 {
                    *modulus - p_coeffs[i]
                } else {
                    p_coeffs[i]
                };
            }
            assert_eq!(&Vec::<u64>::try_from(&q).unwrap(), &v);

            let q_ntt = p_ntt.substitute(&SubstitutionExponent::new(&ctx, 3)?)?;
            q.change_representation(Representation::Ntt);
            assert_eq!(q, q_ntt);

            let q_ntt_shoup = p_ntt_shoup.substitute(&SubstitutionExponent::new(&ctx, 3)?)?;
            q.change_representation(Representation::NttShoup);
            assert_eq!(q, q_ntt_shoup);

            // 11 = 3^(-1) % 16
            assert_eq!(
                p,
                p.substitute(&SubstitutionExponent::new(&ctx, 3)?)?
                    .substitute(&SubstitutionExponent::new(&ctx, 11)?)?
            );
            assert_eq!(
                p_ntt,
                p_ntt
                    .substitute(&SubstitutionExponent::new(&ctx, 3)?)?
                    .substitute(&SubstitutionExponent::new(&ctx, 11)?)?
            );
            assert_eq!(
                p_ntt_shoup,
                p_ntt_shoup
                    .substitute(&SubstitutionExponent::new(&ctx, 3)?)?
                    .substitute(&SubstitutionExponent::new(&ctx, 11)?)?
            );
        }

        let ctx = Arc::new(Context::new(MODULI, 16)?);
        let p = Poly::random(&ctx, Representation::PowerBasis, &mut rng);
        let mut p_ntt = p.clone();
        p_ntt.change_representation(Representation::Ntt);
        let mut p_ntt_shoup = p.clone();
        p_ntt_shoup.change_representation(Representation::NttShoup);

        assert_eq!(
            p,
            p.substitute(&SubstitutionExponent::new(&ctx, 3)?)?
                .substitute(&SubstitutionExponent::new(&ctx, 11)?)?
        );
        assert_eq!(
            p_ntt,
            p_ntt
                .substitute(&SubstitutionExponent::new(&ctx, 3)?)?
                .substitute(&SubstitutionExponent::new(&ctx, 11)?)?
        );
        assert_eq!(
            p_ntt_shoup,
            p_ntt_shoup
                .substitute(&SubstitutionExponent::new(&ctx, 3)?)?
                .substitute(&SubstitutionExponent::new(&ctx, 11)?)?
        );

        Ok(())
    }

    #[test]
    fn switch_down() -> Result<(), Box<dyn Error>> {
        let mut rng = rand::rng();
        let ntests = 100;
        let ctx = Arc::new(Context::new(MODULI, 16)?);

        for _ in 0..ntests {
            // If the polynomial has incorrect representation, an error is returned
            let e = Poly::random(&ctx, Representation::Ntt, &mut rng).switch_down();
            assert!(e.is_err());
            assert_eq!(
                e.unwrap_err(),
                crate::Error::IncorrectRepresentation(
                    Representation::Ntt,
                    Representation::PowerBasis
                )
            );

            // Otherwise, no error happens and the coefficients evolve as expected.
            let mut p = Poly::random(&ctx, Representation::PowerBasis, &mut rng);
            let mut reference = Vec::<BigUint>::from(&p);
            let mut current_ctx = ctx.clone();
            assert_eq!(p.ctx, current_ctx);
            while current_ctx.next_context.is_some() {
                let denominator = current_ctx.modulus().clone();
                current_ctx = current_ctx.next_context.as_ref().unwrap().clone();
                let numerator = current_ctx.modulus().clone();
                assert!(p.switch_down().is_ok());
                assert_eq!(p.ctx, current_ctx);
                let p_biguint = Vec::<BigUint>::from(&p);
                assert_eq!(
                    p_biguint,
                    reference
                        .iter()
                        .map(
                            |b| (((b * &numerator) + (&denominator >> 1)) / &denominator)
                                % current_ctx.modulus()
                        )
                        .collect_vec()
                );
                reference.clone_from(&p_biguint);
            }
        }
        Ok(())
    }

    #[test]
    fn switch_down_to() -> Result<(), Box<dyn Error>> {
        let mut rng = rand::rng();
        let ntests = 100;
        let ctx1 = Arc::new(Context::new(MODULI, 16)?);
        let ctx2 = Arc::new(Context::new(&MODULI[..2], 16)?);

        for _ in 0..ntests {
            let mut p = Poly::random(&ctx1, Representation::PowerBasis, &mut rng);
            let reference = Vec::<BigUint>::from(&p);

            p.switch_down_to(&ctx2)?;

            assert_eq!(p.ctx, ctx2);
            assert_eq!(
                Vec::<BigUint>::from(&p),
                reference
                    .iter()
                    .map(|b| ((b * ctx2.modulus()) + (ctx1.modulus() >> 1)) / ctx1.modulus())
                    .collect_vec()
            );
        }

        Ok(())
    }

    #[test]
    fn switch() -> Result<(), Box<dyn Error>> {
        let mut rng = rand::rng();
        let ntests = 100;
        let ctx1 = Arc::new(Context::new(&MODULI[..2], 16)?);
        let ctx2 = Arc::new(Context::new(&MODULI[3..], 16)?);
        let switcher = Switcher::new(&ctx1, &ctx2)?;
        for _ in 0..ntests {
            let p = Poly::random(&ctx1, Representation::PowerBasis, &mut rng);
            let reference = Vec::<BigUint>::from(&p);

            let q = p.switch(&switcher)?;

            assert_eq!(q.ctx, ctx2);
            assert_eq!(
                Vec::<BigUint>::from(&q),
                reference
                    .iter()
                    .map(|b| ((b * ctx2.modulus()) + (ctx1.modulus() >> 1)) / ctx1.modulus())
                    .collect_vec()
            );
        }
        Ok(())
    }

    #[test]
    fn mul_x_power() -> Result<(), Box<dyn Error>> {
        let mut rng = rand::rng();
        let ctx = Arc::new(Context::new(MODULI, 16)?);
        let e = Poly::random(&ctx, Representation::Ntt, &mut rng).multiply_inverse_power_of_x(1);
        assert!(e.is_err());
        assert_eq!(
            e.unwrap_err(),
            crate::Error::IncorrectRepresentation(Representation::Ntt, Representation::PowerBasis)
        );

        let mut p = Poly::random(&ctx, Representation::PowerBasis, &mut rng);
        let q = p.clone();

        p.multiply_inverse_power_of_x(0)?;
        assert_eq!(p, q);

        p.multiply_inverse_power_of_x(1)?;
        assert_ne!(p, q);

        p.multiply_inverse_power_of_x(2 * ctx.degree - 1)?;
        assert_eq!(p, q);

        p.multiply_inverse_power_of_x(ctx.degree)?;
        assert_eq!(
            Vec::<BigUint>::from(&p)
                .iter()
                .map(|c| ctx.modulus() - c)
                .collect_vec(),
            Vec::<BigUint>::from(&q)
        );

        Ok(())
    }
}
