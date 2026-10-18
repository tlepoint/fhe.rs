//! Arithmetic kernels checked by Verus and compiled into the ordinary Rust
//! library.
//!
//! Verify this file directly; see `verification/README.md` at the repository
//! root.

use vstd::prelude::*;

#[path = "verified/lazy_reduce.rs"]
mod lazy_reduce;
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

/// Add canonical residues for a modulus accepted by `Modulus::new`.
///
/// The contract proves the result and arithmetic safety under these input bounds.
/// It does not prove constant-time machine code or validate callers' inputs.
#[inline]
pub(crate) const fn add_mod(a: u64, b: u64, p: u64) -> (result: u64)
    requires
        2 <= p < 0x4000_0000_0000_0000u64,
        a < p,
        b < p,
    ensures
        result < p,
        result as int == (a as int + b as int) % (p as int),
{
    reduce1(a + b, p)
}

/// Add canonical residues using the verified variable-time reducer.
///
/// # Safety
/// Timing may reveal information about the operands.
#[inline]
pub(crate) const unsafe fn add_mod_vt(a: u64, b: u64, p: u64) -> (result: u64)
    requires
        2 <= p < 0x4000_0000_0000_0000u64,
        a < p,
        b < p,
    ensures
        result < p,
        result as int == (a as int + b as int) % (p as int),
{
    unsafe { reduce1_vt(a + b, p) }
}

}
