//! Scalar kernels checked by Verus and compiled into the ordinary Rust library.
//!
//! Verify this file directly; see `verification/README.md` at the repository
//! root.

use vstd::prelude::*;

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

}
