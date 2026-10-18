//! Arithmetic kernels checked by Verus and compiled into the ordinary Rust
//! library.
//!
//! Verify this file directly; see `verification/README.md` at the repository
//! root.

use vstd::prelude::*;

#[path = "verified/add.rs"]
mod addition;
pub(crate) use addition::{add_mod, add_mod_vt, add_vec, add_vec_vt};

#[path = "verified/sub.rs"]
mod subtraction;
#[cfg(test)]
pub(crate) use subtraction::sub_mod_vt;
pub(crate) use subtraction::{sub_mod, sub_vec, sub_vec_vt};

#[path = "verified/neg.rs"]
mod negation;
#[cfg(test)]
pub(crate) use negation::neg_mod_vt;
pub(crate) use negation::{neg_mod, neg_vec, neg_vec_vt};

#[path = "verified/mul.rs"]
mod multiplication;
pub(crate) use multiplication::{
    lazy_mul_shoup, mul_mod, mul_opt, mul_shoup, mul_shoup_vec, mul_shoup_vec_vt, mul_vec,
    mul_vec_vt, scalar_mul_vec, scalar_mul_vec_vt, shoup, shoup_into, shoup_vec,
};
#[cfg(test)]
pub(crate) use multiplication::{mul_mod_vt, mul_opt_vt, mul_shoup_vt};

#[path = "verified/center.rs"]
mod centering;
pub(crate) use centering::{center, center_vec};

#[path = "verified/pow.rs"]
mod exponentiation;
pub(crate) use exponentiation::pow_mod;

#[path = "verified/random.rs"]
mod sampling;
pub(crate) use sampling::{ResidueSampler, random_into, random_vec};

#[path = "verified/reduce.rs"]
mod reduction;
pub(crate) use reduction::{
    reduce, reduce_opt, reduce_opt_u128, reduce_opt_u128_vt, reduce_opt_vt, reduce_u128,
    reduce_u128_vt, reduce_vec, reduce_vec_i64, reduce_vec_i64_vt, reduce_vec_new,
    reduce_vec_new_vt, reduce_vec_vt, reduce_vt, reduce1, reduce1_vt, reduce3, reduce3_vt,
};
#[cfg(test)]
pub(crate) use reduction::{reduce_i64, reduce_i64_vt};

#[path = "verified/lazy_reduce.rs"]
mod lazy_reduce;
pub use lazy_reduce::supports_opt;
pub(crate) use lazy_reduce::{
    lazy_reduce, lazy_reduce_opt, lazy_reduce_opt_u128, lazy_reduce_u128, lazy_reduce_vec,
};

verus! {

/// Select `on_true` when `cond` is true and `on_false` otherwise.
#[inline]
pub(crate) const fn const_time_cond_select(
    on_true: u64,
    on_false: u64,
    cond: bool,
) -> (result: u64)
    ensures
        result == if cond { on_true } else { on_false },
{
    let mask = -(cond as i64) as u64;
    let diff = on_true ^ on_false;
    let result = (diff & mask) ^ on_false;

    proof {
        assert(result == if cond { on_true } else { on_false }) by (bit_vector)
            requires
                mask == -(cond as i64) as u64,
                diff == on_true ^ on_false,
                result == (diff & mask) ^ on_false,
        ;
    }

    result
}

}
