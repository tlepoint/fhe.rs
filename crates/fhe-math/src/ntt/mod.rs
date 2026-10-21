//! Number-Theoretic Transform in ZZ_q.

mod predicates;
pub(crate) use predicates::supports_ntt;

mod native;

#[cfg(feature = "tfhe-ntt")]
mod tfhe;

#[cfg(not(feature = "tfhe-ntt"))]
pub use native::NttOperator;
#[cfg(feature = "tfhe-ntt")]
pub use tfhe::NttOperator;

impl NttOperator {
    fn validate_input(&self, input: &[u64]) -> crate::Result<()> {
        if input.len() != self.size() {
            return Err(crate::Error::NttLengthMismatch {
                actual: input.len(),
                expected: self.size(),
            });
        }
        // Visit every coefficient: do not reveal the first noncanonical index.
        let invalid = input
            .iter()
            .fold(0_u64, |bad, &x| bad | u64::from(x >= self.modulus()));
        if invalid != 0 {
            return Err(crate::Error::NonCanonicalNttInput {
                modulus: self.modulus(),
            });
        }
        Ok(())
    }

    /// Forward transform of canonical residues, with validation before
    /// mutation.
    pub fn try_forward(&self, input: &mut [u64]) -> crate::Result<()> {
        self.validate_input(input)?;
        self.forward(input);
        Ok(())
    }

    /// Inverse transform of canonical residues, with validation before
    /// mutation.
    pub fn try_backward(&self, input: &mut [u64]) -> crate::Result<()> {
        self.validate_input(input)?;
        self.backward(input);
        Ok(())
    }

    /// Checked variable-time forward transform of explicitly public residues.
    pub fn try_forward_public(
        &self,
        input: &mut [u64],
        _permission: crate::VariableTime,
    ) -> crate::Result<()> {
        self.validate_input(input)?;
        // Length and canonical range were checked; the mutable slice is live
        // and exclusive for the full duration of the raw kernel.
        unsafe { self.forward_vt(input.as_mut_ptr()) };
        Ok(())
    }

    /// Checked variable-time inverse transform of explicitly public residues.
    pub fn try_backward_public(
        &self,
        input: &mut [u64],
        _permission: crate::VariableTime,
    ) -> crate::Result<()> {
        self.validate_input(input)?;
        unsafe { self.backward_vt(input.as_mut_ptr()) };
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use rand::rng;

    use super::{NttOperator, supports_ntt};
    use crate::zq::Modulus;

    #[test]
    fn constructor() {
        for size in [32, 1024] {
            for p in [1153, 4611686018326724609] {
                let q = Modulus::new(p).unwrap();
                let supports_ntt = supports_ntt(p, size);

                let op = NttOperator::new(&q, size);

                if supports_ntt {
                    assert!(op.is_some());
                } else {
                    assert!(op.is_none());
                }
            }
        }
    }

    #[test]
    fn bijection() {
        let ntests = 100;
        let mut rng = rng();

        for size in [32, 1024] {
            for p in [1153, 4611686018326724609] {
                let q = Modulus::new(p).unwrap();

                if supports_ntt(p, size) {
                    let op = NttOperator::new(&q, size).unwrap();

                    for _ in 0..ntests {
                        let mut a = q.random_vec(size, &mut rng);
                        let a_clone = a.clone();
                        let mut b = a.clone();

                        op.forward(&mut a);
                        assert_ne!(a, a_clone);

                        unsafe { op.forward_vt(b.as_mut_ptr()) }
                        assert_eq!(a, b);

                        op.backward(&mut a);
                        assert_eq!(a, a_clone);

                        unsafe { op.backward_vt(b.as_mut_ptr()) }
                        assert_eq!(a, b);
                    }
                }
            }
        }
    }

    #[test]
    fn forward_lazy() {
        let ntests = 100;
        let mut rng = rng();

        for size in [32, 1024] {
            for p in [1153, 4611686018326724609] {
                let q = Modulus::new(p).unwrap();

                if supports_ntt(p, size) {
                    let op = NttOperator::new(&q, size).unwrap();

                    for _ in 0..ntests {
                        let mut a = q.random_vec(size, &mut rng);
                        let mut a_lazy = a.clone();

                        op.forward(&mut a);

                        unsafe {
                            op.forward_vt_lazy(a_lazy.as_mut_ptr());
                            q.reduce_vec(&mut a_lazy);
                        }

                        assert_eq!(a, a_lazy);
                    }
                }
            }
        }
    }
}
