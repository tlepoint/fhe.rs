#![warn(missing_docs, unused_imports)]
// Allow indexing in this performance-critical RNS implementation
#![allow(clippy::indexing_slicing)]

//! Residue-Number System operations.

use crate::{zq::Modulus, Error, Result};
use itertools::{izip, Itertools};
use ndarray::ArrayView1;
use num_bigint::BigUint;
use num_bigint_dig::{BigInt as BigIntDig, BigUint as BigUintDig, ExtendedGcd, ModInverse};
use num_traits::{cast::ToPrimitive, One, Zero};
use std::{cmp::Ordering, fmt::Debug};

mod scaler;

pub use scaler::{RnsScaler, ScalingFactor};

/// Context for a Residue Number System.
#[derive(Default, Clone, PartialEq, Eq)]
pub struct RnsContext {
    moduli_u64: Vec<u64>,
    moduli: Vec<Modulus>,
    q_tilde: Vec<u64>,
    q_tilde_shoup: Vec<u64>,
    q_star: Vec<BigUint>,
    garner: Vec<BigUint>,
    product: BigUint,
}

impl Debug for RnsContext {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("RnsContext")
            .field("moduli_u64", &self.moduli_u64)
            // .field("moduli", &self.moduli)
            // .field("q_tilde", &self.q_tilde)
            // .field("q_tilde_shoup", &self.q_tilde_shoup)
            // .field("q_star", &self.q_star)
            // .field("garner", &self.garner)
            .field("product", &self.product)
            .finish()
    }
}

impl RnsContext {
    /// Create a RNS context from a list of moduli.
    ///
    /// Returns an error if the list is empty, or if the moduli are no coprime.
    pub fn new(moduli_u64: &[u64]) -> Result<Self> {
        if moduli_u64.is_empty() {
            Err(Error::Default("The list of moduli is empty".to_string()))
        } else {
            let mut product = BigUint::one();
            let mut product_dig = BigUintDig::one();

            for i in 0..moduli_u64.len() {
                // Return an error if the moduli are not coprime.
                for j in 0..moduli_u64.len() {
                    if i != j {
                        let (d, _, _) = BigUintDig::from(moduli_u64[i])
                            .extended_gcd(&BigUintDig::from(moduli_u64[j]));
                        if d.cmp(&BigIntDig::from(1)) != Ordering::Equal {
                            return Err(Error::Default("The moduli are not coprime".to_string()));
                        }
                    }
                }

                product *= &BigUint::from(moduli_u64[i]);
                product_dig *= &BigUintDig::from(moduli_u64[i]);
            }

            #[allow(clippy::type_complexity)]
            let (moduli, q_tilde, q_tilde_shoup, q_star, garner): (
                Vec<Modulus>,
                Vec<u64>,
                Vec<u64>,
                Vec<BigUint>,
                Vec<BigUint>,
            ) = moduli_u64
                .iter()
                .map(|modulus| {
                    let m = Modulus::new(*modulus)?;
                    let q_star_i = &product / modulus;
                    let q_tilde_i = (&product_dig / modulus)
                        .mod_inverse(&BigUintDig::from(*modulus))
                        .unwrap()
                        .to_u64()
                        .unwrap();
                    let garner_i = &q_star_i * q_tilde_i;
                    let q_tilde_shoup_i = m.shoup(q_tilde_i);
                    Ok((m, q_tilde_i, q_tilde_shoup_i, q_star_i, garner_i))
                })
                .collect::<Result<Vec<_>>>()?
                .into_iter()
                .multiunzip();

            Ok(Self {
                moduli_u64: moduli_u64.to_owned(),
                moduli,
                q_tilde,
                q_tilde_shoup,
                q_star,
                garner,
                product,
            })
        }
    }

    /// Returns the product of the moduli used when creating the RNS context.
    #[must_use]
    pub const fn modulus(&self) -> &BigUint {
        &self.product
    }

    /// Project a BigUint into its rests.
    #[must_use]
    pub fn project(&self, a: &BigUint) -> Vec<u64> {
        self.moduli_u64
            .iter()
            .map(|modulus| (a % modulus).to_u64().unwrap())
            .collect()
    }

    /// Lift rests into a BigUint.
    ///
    /// Aborts if the number of rests is different than the number of moduli in
    /// debug mode.
    #[must_use]
    pub fn lift(&self, rests: ArrayView1<u64>) -> BigUint {
        let mut result = BigUint::zero();
        izip!(rests.iter(), self.garner.iter())
            .for_each(|(r_i, garner_i)| result += garner_i * *r_i);
        result % &self.product
    }

    /// Getter for the i-th garner coefficient.
    #[must_use]
    pub fn get_garner(&self, i: usize) -> Option<&BigUint> {
        self.garner.get(i)
    }
}

#[cfg(test)]
mod tests {

    use std::error::Error;

    use super::RnsContext;
    use ndarray::ArrayView1;
    use num_bigint::BigUint;
    use rand::RngCore;

    #[test]
    fn constructor() {
        assert!(RnsContext::new(&[2]).is_ok());
        assert!(RnsContext::new(&[2, 3]).is_ok());
        assert!(RnsContext::new(&[4, 15, 1153]).is_ok());

        let e = RnsContext::new(&[]);
        assert!(e.is_err());
        assert_eq!(e.unwrap_err().to_string(), "The list of moduli is empty");
        let e = RnsContext::new(&[2, 4]);
        assert!(e.is_err());
        assert_eq!(e.unwrap_err().to_string(), "The moduli are not coprime");
        let e = RnsContext::new(&[2, 3, 5, 30]);
        assert!(e.is_err());
        assert_eq!(e.unwrap_err().to_string(), "The moduli are not coprime");
    }

    #[test]
    fn garner() -> Result<(), Box<dyn Error>> {
        let rns = RnsContext::new(&[4, 15, 1153])?;

        for i in 0..3 {
            let gi = rns.get_garner(i);
            assert!(gi.is_some());
            assert_eq!(gi.unwrap(), &rns.garner[i]);
        }
        assert!(rns.get_garner(3).is_none());

        Ok(())
    }

    #[test]
    fn modulus() -> Result<(), Box<dyn Error>> {
        let mut rns = RnsContext::new(&[2])?;
        debug_assert_eq!(rns.modulus(), &BigUint::from(2u64));

        rns = RnsContext::new(&[2, 5])?;
        debug_assert_eq!(rns.modulus(), &BigUint::from(2u64 * 5));

        rns = RnsContext::new(&[4, 15, 1153])?;
        debug_assert_eq!(rns.modulus(), &BigUint::from(4u64 * 15 * 1153));

        Ok(())
    }

    #[test]
    fn project_lift() -> Result<(), Box<dyn Error>> {
        let ntests = 100;
        let rns = RnsContext::new(&[4, 15, 1153])?;
        let product = 4u64 * 15 * 1153;

        let mut rests = rns.project(&BigUint::from(0u64));
        assert_eq!(&rests, &[0u64, 0, 0]);
        assert_eq!(rns.lift(ArrayView1::from(&rests)), BigUint::from(0u64));

        rests = rns.project(&BigUint::from(4u64));
        assert_eq!(&rests, &[0u64, 4, 4]);
        assert_eq!(rns.lift(ArrayView1::from(&rests)), BigUint::from(4u64));

        rests = rns.project(&BigUint::from(15u64));
        assert_eq!(&rests, &[3u64, 0, 15]);
        assert_eq!(rns.lift(ArrayView1::from(&rests)), BigUint::from(15u64));

        rests = rns.project(&BigUint::from(1153u64));
        assert_eq!(&rests, &[1u64, 13, 0]);
        assert_eq!(rns.lift(ArrayView1::from(&rests)), BigUint::from(1153u64));

        rests = rns.project(&BigUint::from(product - 1));
        assert_eq!(&rests, &[3u64, 14, 1152]);
        assert_eq!(
            rns.lift(ArrayView1::from(&rests)),
            BigUint::from(product - 1)
        );

        let mut rng = rand::rng();

        for _ in 0..ntests {
            let b = BigUint::from(rng.next_u64() % product);
            rests = rns.project(&b);
            assert_eq!(rns.lift(ArrayView1::from(&rests)), b);
        }

        Ok(())
    }
}
