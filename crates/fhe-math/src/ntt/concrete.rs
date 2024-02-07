use concrete_ntt::prime64::Plan;

use crate::zq::Modulus;

use super::native;

/// Number-Theoretic Transform operator.
#[derive(Debug, Clone)]
pub struct NttOperator {
    concrete_operator: Option<Plan>,
    native_operator: native::NttOperator,
}

impl PartialEq for NttOperator {
    fn eq(&self, other: &Self) -> bool {
        self.native_operator == other.native_operator
    }
}

impl Eq for NttOperator {}

impl NttOperator {
    /// Create an NTT operator given a modulus for a specific size.
    ///
    /// Aborts if the size is not a power of 2 that is >= 8 in debug mode.
    /// Returns None if the modulus does not support the NTT for this specific
    /// size.
    pub fn new(p: &Modulus, size: usize) -> Option<Self> {
        let native_operator = native::NttOperator::new(p, size)?;
        let concrete_operator = Plan::try_new(size, p.p);
        Some(Self {
            concrete_operator,
            native_operator,
        })
    }

    /// Compute the forward NTT in place.
    /// Aborts if a is not of the size handled by the operator.
    pub fn forward(&self, a: &mut [u64]) {
        if let Some(ref concrete_operator) = self.concrete_operator {
            concrete_operator.fwd(a);
        } else {
            self.native_operator.forward(a);
        }
    }

    /// Compute the backward NTT in place.
    /// Aborts if a is not of the size handled by the operator.
    pub fn backward(&self, a: &mut [u64]) {
        if let Some(ref concrete_operator) = self.concrete_operator {
            concrete_operator.inv(a);
            concrete_operator.normalize(a);
        } else {
            self.native_operator.backward(a);
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
        if let Some(ref concrete_operator) = self.concrete_operator {
            let a = std::slice::from_raw_parts_mut(a_ptr, concrete_operator.ntt_size());
            concrete_operator.fwd(a);
        } else {
            self.native_operator.forward_vt_lazy(a_ptr);
        }
    }

    /// Compute the forward NTT in place in variable time.
    ///
    /// # Safety
    /// This function assumes that a_ptr points to at least `size` elements.
    /// This function is not constant time and its timing may reveal information
    /// about the value being reduced.
    pub unsafe fn forward_vt(&self, a_ptr: *mut u64) {
        if let Some(ref concrete_operator) = self.concrete_operator {
            let a = std::slice::from_raw_parts_mut(a_ptr, concrete_operator.ntt_size());
            concrete_operator.fwd(a);
        } else {
            self.native_operator.forward_vt(a_ptr);
        }
    }

    /// Compute the backward NTT in place in variable time.
    ///
    /// # Safety
    /// This function assumes that a_ptr points to at least `size` elements.
    /// This function is not constant time and its timing may reveal information
    /// about the value being reduced.
    pub unsafe fn backward_vt(&self, a_ptr: *mut u64) {
        if let Some(ref concrete_operator) = self.concrete_operator {
            let a = std::slice::from_raw_parts_mut(a_ptr, concrete_operator.ntt_size());
            concrete_operator.inv(a);
            concrete_operator.normalize(a);
        } else {
            self.native_operator.backward_vt(a_ptr);
        }
    }
}
