//! The production NTT butterflies, verified under explicit operand and twiddle
//! bounds. Construction of the operator and the transform loops are separate
//! proof obligations.

#![cfg_attr(verus_keep_ghost, allow(unused_imports))]

use super::NttOperator;
use crate::zq::Modulus;
#[cfg(verus_keep_ghost)]
use crate::zq::verified::{shoup_quotient, valid_shoup};
#[cfg(verus_keep_ghost)]
use vstd::arithmetic::div_mod::{
    lemma_add_mod_noop, lemma_mod_mod, lemma_mod_multiples_vanish, lemma_mul_mod_noop_left,
    lemma_sub_mod_noop,
};
use vstd::prelude::*;

verus! {

/// Replace the reduced intermediates in a forward butterfly by their inputs.
proof fn forward_congruence(a: int, b: int, w: int, p: int, reduced: int, product: int)
    requires
        p > 0,
        reduced == a % (2 * p),
        product % p == (b * w) % p,
    ensures
        (reduced + product) % p == (a + b * w) % p,
        (reduced + 2 * p - product) % p == (a - b * w) % p,
{
    lemma_mod_mod(a, p, 2);
    lemma_add_mod_noop(reduced, product, p);
    lemma_add_mod_noop(a, b * w, p);
    lemma_sub_mod_noop(reduced, product, p);
    lemma_sub_mod_noop(a, b * w, p);
    lemma_mod_multiples_vanish(2, reduced - product, p);
}

/// The inverse butterfly's shifted difference represents a - b modulo p.
proof fn inverse_congruence(a: int, b: int, z: int, p: int)
    requires p > 0,
    ensures
        ((a + b) % (2 * p)) % p == (a + b) % p,
        ((2 * p + a - b) * z) % p == ((a - b) * z) % p,
{
    lemma_mod_mod(a + b, p, 2);
    lemma_mod_multiples_vanish(2, a - b, p);
    lemma_mul_mod_noop_left(2 * p + a - b, z, p);
    lemma_mul_mod_noop_left(a - b, z, p);
}

impl NttOperator {
    /// NTT butterfly.
    pub(super) fn butterfly(
        &self, x: &mut u64, y: &mut u64, w: u64, w_shoup: u64,
    )
        requires
            self.p_twice == 2 * self.p.value(),
            *old(x) < 4 * self.p.value(),
            *old(y) < 4 * self.p.value(),
            valid_shoup(w, self.p.value(), w_shoup),
        ensures
            *final(x) < 4 * self.p.value(),
            *final(y) < 4 * self.p.value(),
            *final(x) as int % self.p.value() as int
                == (*old(x) as int + *old(y) as int * w as int) % self.p.value() as int,
            *final(y) as int % self.p.value() as int
                == (*old(x) as int - *old(y) as int * w as int) % self.p.value() as int,
    {
        let p = self.p.as_u64();
        debug_assert!(*x < 4 * p);
        debug_assert!(*y < 4 * p);
        debug_assert!(w < p);
        proof { shoup_quotient(w, p, w_shoup); }
        debug_assert!(self.p.shoup(w) == w_shoup);

        let ghost a = *x as int;
        let ghost b = *y as int;
        *x = Modulus::reduce1(*x, self.p_twice);
        let t = self.p.lazy_mul_shoup(*y, w, w_shoup);
        proof { forward_congruence(a, b, w as int, p as int, *x as int, t as int); }
        *y = *x + self.p_twice - t;
        *x += t;

        debug_assert!(*x < 4 * p);
        debug_assert!(*y < 4 * p);
    }

    /// NTT butterfly in variable time.
    ///
    /// # Safety
    /// Timing may reveal information about the operands.
    pub(super) unsafe fn butterfly_vt(
        &self, x: &mut u64, y: &mut u64, w: u64, w_shoup: u64,
    )
        requires
            self.p_twice == 2 * self.p.value(),
            *old(x) < 4 * self.p.value(),
            *old(y) < 4 * self.p.value(),
            valid_shoup(w, self.p.value(), w_shoup),
        ensures
            *final(x) < 4 * self.p.value(),
            *final(y) < 4 * self.p.value(),
            *final(x) as int % self.p.value() as int
                == (*old(x) as int + *old(y) as int * w as int) % self.p.value() as int,
            *final(y) as int % self.p.value() as int
                == (*old(x) as int - *old(y) as int * w as int) % self.p.value() as int,
    {
        let p = self.p.as_u64();
        debug_assert!(*x < 4 * p);
        debug_assert!(*y < 4 * p);
        debug_assert!(w < p);
        proof { shoup_quotient(w, p, w_shoup); }
        debug_assert!(self.p.shoup(w) == w_shoup);

        let ghost a = *x as int;
        let ghost b = *y as int;
        *x = unsafe { Modulus::reduce1_vt(*x, self.p_twice) };
        let t = self.p.lazy_mul_shoup(*y, w, w_shoup);
        proof { forward_congruence(a, b, w as int, p as int, *x as int, t as int); }
        *y = *x + self.p_twice - t;
        *x += t;

        debug_assert!(*x < 4 * p);
        debug_assert!(*y < 4 * p);
    }

    /// Inverse NTT butterfly.
    pub(super) fn inv_butterfly(
        &self, x: &mut u64, y: &mut u64, z: u64, z_shoup: u64,
    )
        requires
            self.p_twice == 2 * self.p.value(),
            *old(x) < 2 * self.p.value(),
            *old(y) < 2 * self.p.value(),
            valid_shoup(z, self.p.value(), z_shoup),
        ensures
            *final(x) < 2 * self.p.value(),
            *final(y) < 2 * self.p.value(),
            *final(x) as int % self.p.value() as int
                == (*old(x) as int + *old(y) as int) % self.p.value() as int,
            *final(y) as int % self.p.value() as int
                == ((*old(x) as int - *old(y) as int) * z as int) % self.p.value() as int,
    {
        let p = self.p.as_u64();
        debug_assert!(*x < 2 * p);
        debug_assert!(*y < 2 * p);
        debug_assert!(z < p);
        proof { shoup_quotient(z, p, z_shoup); }
        debug_assert!(self.p.shoup(z) == z_shoup);

        proof { inverse_congruence(*x as int, *y as int, z as int, p as int); }
        let t = *x;
        *x = Modulus::reduce1(*y + t, self.p_twice);
        *y = self.p.lazy_mul_shoup(self.p_twice + t - *y, z, z_shoup);

        debug_assert!(*x < 2 * p);
        debug_assert!(*y < 2 * p);
    }

    /// Inverse NTT butterfly in variable time.
    ///
    /// # Safety
    /// Timing may reveal information about the operands.
    pub(super) unsafe fn inv_butterfly_vt(
        &self, x: &mut u64, y: &mut u64, z: u64, z_shoup: u64,
    )
        requires
            self.p_twice == 2 * self.p.value(),
            *old(x) < 2 * self.p.value(),
            *old(y) < 2 * self.p.value(),
            valid_shoup(z, self.p.value(), z_shoup),
        ensures
            *final(x) < 2 * self.p.value(),
            *final(y) < 2 * self.p.value(),
            *final(x) as int % self.p.value() as int
                == (*old(x) as int + *old(y) as int) % self.p.value() as int,
            *final(y) as int % self.p.value() as int
                == ((*old(x) as int - *old(y) as int) * z as int) % self.p.value() as int,
    {
        let p = self.p.as_u64();
        debug_assert!(*x < 2 * p);
        debug_assert!(*y < 2 * p);
        debug_assert!(z < p);
        proof { shoup_quotient(z, p, z_shoup); }
        debug_assert!(self.p.shoup(z) == z_shoup);

        proof { inverse_congruence(*x as int, *y as int, z as int, p as int); }
        let t = *x;
        *x = unsafe { Modulus::reduce1_vt(*y + t, self.p_twice) };
        *y = self.p.lazy_mul_shoup(self.p_twice + t - *y, z, z_shoup);

        debug_assert!(*x < 2 * p);
        debug_assert!(*y < 2 * p);
    }

}

}
