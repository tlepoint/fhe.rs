use tfhe_ntt::prime64::Plan;

use crate::zq::Modulus;

use super::native;

/// Number-Theoretic Transform operator.
#[derive(Debug, Clone)]
pub struct NttOperator {
    modulus: u64,
    size: usize,
    backend: Backend,
}

#[derive(Debug, Clone)]
enum Backend {
    Tfhe(Plan),
    Native(native::NttOperator),
}

impl PartialEq for NttOperator {
    fn eq(&self, other: &Self) -> bool {
        // Construction deterministically selects the backend for these parameters.
        self.modulus == other.modulus && self.size == other.size
    }
}

impl Eq for NttOperator {}

impl NttOperator {
    /// Create an NTT operator given a modulus for a specific size.
    ///
    /// Aborts if the size is not a power of 2 that is >= 8 in debug mode.
    /// Returns None if the modulus does not support the NTT for this specific
    /// size.
    #[must_use]
    pub fn new(p: &Modulus, size: usize) -> Option<Self> {
        if !super::supports_ntt(p.p, size) {
            return None;
        }
        let backend = match Plan::try_new(size, p.p) {
            Some(plan) => Backend::Tfhe(plan),
            None => Backend::Native(native::NttOperator::new(p, size)?),
        };
        Some(Self {
            modulus: p.p,
            size,
            backend,
        })
    }

    /// Compute the forward NTT in place.
    /// Aborts if a is not of the size handled by the operator.
    pub fn forward(&self, a: &mut [u64]) {
        match &self.backend {
            Backend::Tfhe(tfhe_operator) => {
                tfhe_operator.fwd(a);
            }
            Backend::Native(op) => op.forward(a),
        }
    }

    /// Compute the backward NTT in place.
    /// Aborts if a is not of the size handled by the operator.
    pub fn backward(&self, a: &mut [u64]) {
        match &self.backend {
            Backend::Tfhe(tfhe_operator) => {
                tfhe_operator.inv(a);
                tfhe_operator.normalize(a);
            }
            Backend::Native(op) => op.backward(a),
        }
    }

    /// Compute the forward NTT in place in variable time in a lazily fashion.
    /// This means that the output coefficients may be up to 4 times the
    /// modulus.
    ///
    /// # Safety
    /// This function assumes that a_ptr points to at least `size` elements.
    /// This function is not constant time and its timing may reveal information
    /// about the value being reduced.
    pub(crate) unsafe fn forward_vt_lazy(&self, a_ptr: *mut u64) {
        match &self.backend {
            Backend::Tfhe(tfhe_operator) => {
                let a = unsafe { std::slice::from_raw_parts_mut(a_ptr, tfhe_operator.ntt_size()) };
                tfhe_operator.fwd(a);
            }
            Backend::Native(op) => unsafe { op.forward_vt_lazy(a_ptr) },
        }
    }

    /// Compute the forward NTT in place in variable time.
    ///
    /// # Safety
    /// This function assumes that a_ptr points to at least `size` elements.
    /// This function is not constant time and its timing may reveal information
    /// about the value being reduced.
    pub unsafe fn forward_vt(&self, a_ptr: *mut u64) {
        match &self.backend {
            Backend::Tfhe(tfhe_operator) => {
                let a = unsafe { std::slice::from_raw_parts_mut(a_ptr, tfhe_operator.ntt_size()) };
                tfhe_operator.fwd(a);
            }
            Backend::Native(op) => unsafe { op.forward_vt(a_ptr) },
        }
    }

    /// Compute the backward NTT in place in variable time.
    ///
    /// # Safety
    /// This function assumes that a_ptr points to at least `size` elements.
    /// This function is not constant time and its timing may reveal information
    /// about the value being reduced.
    pub unsafe fn backward_vt(&self, a_ptr: *mut u64) {
        match &self.backend {
            Backend::Tfhe(tfhe_operator) => {
                let a = unsafe { std::slice::from_raw_parts_mut(a_ptr, tfhe_operator.ntt_size()) };
                tfhe_operator.inv(a);
                tfhe_operator.normalize(a);
            }
            Backend::Native(op) => unsafe { op.backward_vt(a_ptr) },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{Backend, NttOperator, native};
    use crate::zq::Modulus;

    #[test]
    fn selected_backend_preserves_transforms_and_ring_products() {
        for (modulus, size) in [
            (1153, 8),
            (1153, 16),
            (1153, 32),
            (4611686018326724609, 1024),
        ] {
            let q = Modulus::new(modulus).unwrap();
            let selected = NttOperator::new(&q, size).unwrap();
            // The pinned TFHE backend supports sizes >= 16; size 8 uses native.
            assert_eq!(matches!(&selected.backend, Backend::Native(_)), size == 8);
            let reference = native::NttOperator::new(&q, size).unwrap();
            assert_eq!(selected, NttOperator::new(&q, size).unwrap());
            assert_eq!(selected, selected.clone());
            for input in [
                vec![modulus - 1; size],
                (0..size).map(|i| i as u64 % modulus).collect(),
            ] {
                let mut expected = input.clone();
                // Reproduce the previous wrapper's backend selection. Different
                // roots mean native and TFHE transform vectors need not match.
                if let Some(plan) = tfhe_ntt::prime64::Plan::try_new(size, modulus) {
                    plan.fwd(&mut expected);
                } else {
                    reference.forward(&mut expected);
                }
                let mut actual = input.clone();
                selected.forward(&mut actual);
                assert_eq!(actual, expected);
                assert!(actual.iter().all(|x| *x < modulus));
                selected.backward(&mut actual);
                assert_eq!(actual, input);
                let mut vt = input.clone();
                unsafe {
                    selected.forward_vt(vt.as_mut_ptr());
                }
                assert_eq!(vt, expected);
                unsafe {
                    selected.backward_vt(vt.as_mut_ptr());
                }
                assert_eq!(vt, input);
                let mut lazy = input.clone();
                unsafe {
                    selected.forward_vt_lazy(lazy.as_mut_ptr());
                }
                assert!(lazy.iter().all(|x| *x < 4 * modulus));
                q.reduce_vec(&mut lazy);
                assert_eq!(lazy, expected);

                // Compare the ring operation after inverse transformation,
                // where the two backends must agree despite different roots.
                let mut native_product = input.clone();
                reference.forward(&mut native_product);
                for x in &mut native_product {
                    *x = q.mul(*x, *x);
                }
                reference.backward(&mut native_product);
                let mut selected_product = expected;
                for x in &mut selected_product {
                    *x = q.mul(*x, *x);
                }
                selected.backward(&mut selected_product);
                assert_eq!(selected_product, native_product);
            }
        }
        let q = Modulus::new(1153).unwrap();
        assert!(NttOperator::new(&q, 128).is_none());
        assert_ne!(NttOperator::new(&q, 8), NttOperator::new(&q, 16));
        assert_ne!(
            NttOperator::new(&q, 16),
            NttOperator::new(&Modulus::new(2017).unwrap(), 16)
        );
    }
}
