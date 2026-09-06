use itertools::Itertools;
use num_bigint::BigUint;
use std::{fmt::Debug, sync::Arc};

use crate::{Error, Result, ntt::NttOperator, rns::RnsContext, zq::Modulus};

/// Struct that holds the context associated with elements in rq.
#[derive(Default, Clone, PartialEq, Eq)]
pub struct Context {
    pub(crate) moduli: Box<[u64]>,
    pub(crate) q: Box<[Modulus]>,
    pub(crate) rns: Arc<RnsContext>,
    pub(crate) ops: Box<[Arc<NttOperator>]>,
    pub(crate) degree: usize,
    pub(crate) bitrev: Arc<[usize]>,
    pub(crate) inv_last_qi_mod_qj: Box<[u64]>,
    pub(crate) inv_last_qi_mod_qj_shoup: Box<[u64]>,
    pub(crate) next_context: Option<Arc<Context>>,
}

impl Debug for Context {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Context")
            .field("moduli", &self.moduli)
            // .field("q", &self.q)
            // .field("rns", &self.rns)
            // .field("ops", &self.ops)
            // .field("degree", &self.degree)
            // .field("bitrev", &self.bitrev)
            // .field("inv_last_qi_mod_qj", &self.inv_last_qi_mod_qj)
            // .field("inv_last_qi_mod_qj_shoup", &self.inv_last_qi_mod_qj_shoup)
            .field("next_context", &self.next_context)
            .finish()
    }
}

impl Context {
    /// Creates a context from a list of moduli and a polynomial degree.
    ///
    /// Returns an error if the moduli are not primes less than 62 bits which
    /// supports the NTT of size `degree`.
    pub fn new(moduli: &[u64], degree: usize) -> Result<Self> {
        if !degree.is_power_of_two() || degree < 8 {
            Err(Error::InvalidPolynomialDegree { degree, minimum: 8 })
        } else {
            let rns = Arc::new(RnsContext::new(moduli)?);
            let (q, ops): (Vec<Modulus>, Vec<Arc<NttOperator>>) = moduli
                .iter()
                .map(|modulus| {
                    let qi = Modulus::new(*modulus)?;
                    NttOperator::new(&qi, degree)
                        .ok_or(Error::NttOperatorUnavailable {
                            modulus: *modulus,
                            degree,
                        })
                        .map(|op| (qi, Arc::new(op)))
                })
                .collect::<Result<Vec<(Modulus, Arc<NttOperator>)>>>()?
                .into_iter()
                .unzip();
            let bitrev = (0..degree)
                .map(|j| j.reverse_bits() >> (degree.leading_zeros() + 1))
                .collect_vec();

            Self::from_shared_operators(moduli, &q, &ops, degree, bitrev.into(), rns)
        }
    }

    // NTT operators and the bit-reversal permutation are immutable and shared
    // across the whole chain. CRT products and switching constants vary by level.
    fn from_shared_operators(
        moduli: &[u64],
        q: &[Modulus],
        ops: &[Arc<NttOperator>],
        degree: usize,
        bitrev: Arc<[usize]>,
        rns: Arc<RnsContext>,
    ) -> Result<Self> {
        let mut inv_last_qi_mod_qj = vec![];
        let mut inv_last_qi_mod_qj_shoup = vec![];
        let q_last = moduli.last().unwrap();
        for qi in &q[..q.len() - 1] {
            let inv = qi.inv(qi.reduce(*q_last)).unwrap();
            inv_last_qi_mod_qj.push(inv);
            inv_last_qi_mod_qj_shoup.push(qi.shoup(inv));
        }
        let next_context = if moduli.len() >= 2 {
            let prefix = moduli.len() - 1;
            let child_rns = Arc::new(RnsContext::new(&moduli[..prefix])?);
            Some(Arc::new(Self::from_shared_operators(
                &moduli[..prefix],
                &q[..prefix],
                &ops[..prefix],
                degree,
                bitrev.clone(),
                child_rns,
            )?))
        } else {
            None
        };
        Ok(Self {
            moduli: moduli.into(),
            q: q.into(),
            ops: ops.into(),
            degree,
            bitrev,
            rns,
            inv_last_qi_mod_qj: inv_last_qi_mod_qj.into_boxed_slice(),
            inv_last_qi_mod_qj_shoup: inv_last_qi_mod_qj_shoup.into_boxed_slice(),
            next_context,
        })
    }

    /// Creates a context in an `Arc`.
    pub fn new_arc(moduli: &[u64], degree: usize) -> Result<Arc<Self>> {
        Self::new(moduli, degree).map(Arc::new)
    }

    /// Returns the modulus as a BigUint.
    #[must_use]
    pub fn modulus(&self) -> &BigUint {
        self.rns.modulus()
    }

    /// Returns a reference to the moduli in this context.
    #[must_use]
    pub fn moduli(&self) -> &[u64] {
        &self.moduli
    }

    /// Returns a reference to the moduli as Modulus in this context.
    #[must_use]
    pub fn moduli_operators(&self) -> &[Modulus] {
        &self.q
    }

    /// Returns the number of iterations to switch to a children context.
    /// Returns an error if the context provided is not a child context.
    pub fn niterations_to(&self, context: &Arc<Context>) -> Result<usize> {
        if std::ptr::eq(context.as_ref(), self) || context.as_ref() == self {
            return Ok(0);
        }

        let mut current = self;
        let mut niterations = 0;
        while let Some(next) = &current.next_context {
            niterations += 1;
            if next == context {
                return Ok(niterations);
            }
            current = next;
        }
        Err(Error::ContextNotReachable)
    }

    /// Returns a shared context after `i` iterations, including this Arc at
    /// level zero.
    pub fn context_at_level(self: &Arc<Self>, i: usize) -> Result<Arc<Self>> {
        if i >= self.moduli.len() {
            Err(Error::InvalidContextLevel {
                level: i,
                max_level: self.moduli.len().saturating_sub(1),
            })
        } else {
            let mut current = self;
            for _ in 0..i {
                current = current.next_context.as_ref().unwrap();
            }
            Ok(current.clone())
        }
    }
}

#[cfg(test)]
mod tests {
    use std::{error::Error, sync::Arc};

    use crate::ntt::supports_ntt;
    use crate::rq::Context;

    const MODULI: &[u64; 5] = &[
        1153,
        4611686018326724609,
        4611686018309947393,
        4611686018232352769,
        4611686018171535361,
    ];

    #[test]
    fn levels_share_tables_and_keep_independent_crt_data() -> Result<(), Box<dyn Error>> {
        let root = Context::new_arc(MODULI, 16)?;
        for level in 0..MODULI.len() {
            let child = root.context_at_level(level)?;
            let separate = Context::new_arc(&MODULI[..MODULI.len() - level], 16)?;
            assert_eq!(child, separate);
            assert!(Arc::ptr_eq(&root.bitrev, &child.bitrev));
            for (original, shared) in root.ops.iter().zip(child.ops.iter()) {
                assert!(Arc::ptr_eq(original, shared));
            }
            if level != 0 {
                assert!(!Arc::ptr_eq(&root.rns, &child.rns));
            }
            let poly = crate::rq::Poly::<crate::rq::PowerBasis>::random_from_seed(&child, [3; 32]);
            assert_eq!(poly.clone().into_ntt().into_power_basis(), poly);
        }
        let child = root.context_at_level(2)?;
        let reference = Context::new_arc(&MODULI[..MODULI.len() - 2], 16)?;
        drop(root);
        assert_eq!(child, reference);
        Ok(())
    }

    #[test]
    fn level_lookup_reuses_children_and_accepts_equal_contexts() -> Result<(), Box<dyn Error>> {
        let ctx = Context::new_arc(MODULI, 16)?;
        assert!(Arc::ptr_eq(&ctx.context_at_level(0)?, &ctx));
        let mut child = ctx.clone();
        for level in 1..MODULI.len() {
            child = child.next_context.as_ref().unwrap().clone();
            assert!(Arc::ptr_eq(&ctx.context_at_level(level)?, &child));
            let separate = Context::new_arc(&MODULI[..MODULI.len() - level], 16)?;
            assert_eq!(ctx.niterations_to(&separate)?, level);
        }
        assert!(ctx.context_at_level(MODULI.len()).is_err());
        assert!(ctx.niterations_to(&Context::new_arc(MODULI, 8)?).is_err());
        Ok(())
    }

    #[test]
    fn context_constructor() {
        for modulus in MODULI {
            // modulus is = 1 modulo 2 * 8
            assert!(Context::new(&[*modulus], 16).is_ok());

            if supports_ntt(*modulus, 128) {
                assert!(Context::new(&[*modulus], 128).is_ok());
            } else {
                assert!(Context::new(&[*modulus], 128).is_err());
            }
        }

        // All moduli in MODULI are = 1 modulo 2 * 8
        assert!(Context::new(MODULI, 16).is_ok());

        // This should fail since 1153 != 1 moduli 2 * 128
        assert!(Context::new(MODULI, 128).is_err());
    }

    #[test]
    fn next_context() -> Result<(), Box<dyn Error>> {
        // A context should have a children pointing to a context with one less modulus.
        let context = Arc::new(Context::new(MODULI, 16)?);
        assert_eq!(
            context.next_context,
            Some(Arc::new(Context::new(&MODULI[..MODULI.len() - 1], 16)?))
        );

        // We can go down the chain of the MODULI.len() - 1 context's.
        let mut number_of_children = 0;
        let mut current = context;
        while current.next_context.is_some() {
            number_of_children += 1;
            current = current.next_context.as_ref().unwrap().clone();
        }
        assert_eq!(number_of_children, MODULI.len() - 1);

        Ok(())
    }

    #[test]
    fn niterations_to() -> Result<(), Box<dyn Error>> {
        // A context should have a children pointing to a context with one less modulus.
        let context = Arc::new(Context::new(MODULI, 16)?);

        assert_eq!(context.niterations_to(&context).ok(), Some(0));

        assert_eq!(
            context
                .niterations_to(&Arc::new(Context::new(&MODULI[1..], 16)?))
                .err(),
            Some(crate::Error::ContextNotReachable)
        );

        for i in 1..MODULI.len() {
            assert_eq!(
                context
                    .niterations_to(&Arc::new(Context::new(&MODULI[..MODULI.len() - i], 16)?))
                    .ok(),
                Some(i)
            );
        }

        Ok(())
    }
}
