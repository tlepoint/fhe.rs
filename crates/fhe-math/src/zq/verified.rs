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

#[path = "verified/reduce.rs"]
mod reduction;
pub(crate) use reduction::{
    reduce, reduce_opt, reduce_opt_u128, reduce_opt_u128_vt, reduce_opt_vt, reduce_u128,
    reduce_u128_vt, reduce_vec, reduce_vec_i64, reduce_vec_i64_vt, reduce_vec_new,
    reduce_vec_new_vt, reduce_vec_vt, reduce_vt, reduce3, reduce3_vt,
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

/// Reduce a value below `2*p` to its canonical residue.
///
/// The wider modulus bound includes the NTT's calls with twice a 62-bit modulus.
#[inline]
pub(crate) const fn reduce1(x: u64, p: u64) -> (result: u64)
    requires
        0 < p < 0x8000_0000_0000_0000u64,
        x < 2 * p,
    ensures
        result < p,
        result as int == (x as int) % (p as int),
{
    let result = const_time_cond_select(x, x.wrapping_sub(p), x < p);

    proof {
        // The selection contract supplies the two cases; one subtraction suffices.
        assert((x as int) % (p as int) == if x < p {
            x as int
        } else {
            x as int - p as int
        }) by (nonlinear_arith)
            requires
                0 < p,
                x < 2 * p,
        ;
    }

    result
}

/// Reduce a value below `2*p` using the target's variable-time implementation.
///
/// # Safety
/// Timing may reveal information about `x`.
#[cfg(any(target_os = "macos", target_feature = "avx2"))]
#[inline]
pub(crate) const unsafe fn reduce1_vt(x: u64, p: u64) -> (result: u64)
    requires
        0 < p < 0x8000_0000_0000_0000u64,
        x < 2 * p,
    ensures
        result < p,
        result as int == (x as int) % (p as int),
{
    proof {
        assert((x as int) % (p as int) == if x >= p {
            x as int - p as int
        } else {
            x as int
        }) by (nonlinear_arith)
            requires
                0 < p,
                x < 2 * p,
        ;
    }

    if x >= p { x - p } else { x }
}

/// Use the verified selector-based reduction on targets without the branch path.
///
/// # Safety
/// Callers must permit variable timing, as on the other targets.
#[cfg(all(not(target_os = "macos"), not(target_feature = "avx2")))]
#[inline]
pub(crate) const unsafe fn reduce1_vt(x: u64, p: u64) -> (result: u64)
    requires
        0 < p < 0x8000_0000_0000_0000u64,
        x < 2 * p,
    ensures
        result < p,
        result as int == (x as int) % (p as int),
{
    reduce1(x, p)
}

}
