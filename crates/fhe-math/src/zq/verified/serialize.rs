//! Modulus-width serialization, composed with fhe-util's checked packing.

#[cfg(verus_keep_ghost)]
use super::MODULUS_LIMIT;
#[cfg(verus_keep_ghost)]
use super::lazy_reduce::leading_zero_bounds;
use super::packing;
#[cfg(verus_keep_ghost)]
use vstd::arithmetic::{div_mod::*, power2::*};
#[cfg(verus_keep_ghost)]
use vstd::bits::{lemma_u64_pow2_no_overflow, lemma_u64_shl_is_mul};
use vstd::prelude::*;
#[cfg(verus_keep_ghost)]
use vstd::std_specs::bits::u64_leading_zeros;

verus! {

pub open spec fn valid_modulus(p: u64) -> bool {
    2 <= p < MODULUS_LIMIT
}

pub open spec fn width(p: u64) -> int {
    64 - u64_leading_zeros((p - 1) as u64)
}

/// Exactly ceil(log2(p)), including powers of two and the modulus two.
proof fn width_bounds(p: u64)
    requires valid_modulus(p),
    ensures
        1 <= width(p) <= 62,
        pow2((width(p) - 1) as nat) < p <= pow2(width(p) as nat),
{
    leading_zero_bounds((p - 1) as u64);
    if p == 2 {
        leading_zero_bounds(0);
        reveal(u64_leading_zeros);
    }
}

/// Number of bits needed to encode every canonical residue.
pub(crate) const fn serialization_width(p: u64) -> (r: usize)
    requires valid_modulus(p),
    ensures r == width(p), 1 <= r <= 62,
        pow2((r - 1) as nat) < p <= pow2(r as nat),
{
    proof { width_bounds(p); }
    64 - (p - 1).leading_zeros() as usize
}

/// Canonical residues always fit; fitting is not a canonicality check.
pub proof fn canonical_words_fit(a: Seq<u64>, p: u64)
    requires valid_modulus(p), forall|i: int| 0 <= i < a.len() ==> #[trigger] a[i] < p,
    ensures packing::valid_words(a, width(p)),
{
    width_bounds(p);
    let w = width(p) as u64;
    lemma_u64_pow2_no_overflow(w as nat);
    lemma_u64_shl_is_mul(1, w);
    assert forall|i: int| 0 <= i < a.len() implies packing::fits(#[trigger] a[i], width(p)) by {
        let x = a[i];
        assert(x & !(u64::MAX >> ((64 - w) as u64)) == 0) by (bit_vector)
            requires 1 <= w <= 62, x < (1u64 << w);
    }
    assert forall|i: int| 0 <= i < a.len() implies packing::fits(#[trigger] packing::Word::value(a[i]), width(p)) by {
        assert(packing::fits(a[i], width(p)));
    }
}

/// Checked byte count for the public API's groups of eight words.
/// Dividing first avoids an intermediate bit-count overflow.
pub(crate) const fn serialization_length(p: u64, size: usize) -> (r: Option<usize>)
    requires valid_modulus(p),
    ensures
        r.is_some() <==> size % 8 == 0 && size * width(p) / 8 <= usize::MAX,
        r matches Some(bytes) ==> bytes == size * width(p) / 8,
{
    if !size.is_multiple_of(8) {
        return None;
    }
    let w = serialization_width(p);
    proof {
        lemma_fundamental_div_mod(size as int, 8);
        assert(w * (size as int / 8) == size * (w as int) / 8) by (nonlinear_arith)
            requires size as int % 8 == 0;
    }
    w.checked_mul(size / 8)
}

/// Append at a byte boundary. The preconditions prove the existing unwrap safe.
pub(crate) fn serialize_into(a: &[u64], p: u64, out: &mut Vec<u8>)
    requires
        valid_modulus(p),
        packing::valid_words(a@, width(p)),
        a.len() * width(p) <= usize::MAX,
        old(out).len() + (a.len() * width(p) + 7) / 8 <= usize::MAX,
    ensures
        final(out)@.take(old(out).len() as int) == old(out)@,
        packing::packed(a@, width(p), final(out)@.skip(old(out).len() as int), 8),
{
    let w = serialization_width(p);
    packing::transcode_to_bytes_into(a, w, out).unwrap();
}

/// Pack words at the modulus' width, preserving the packing error semantics.
pub(crate) fn serialize_vec(a: &[u64], p: u64)
    -> (r: Result<Vec<u8>, packing::TranscodeError>)
    requires valid_modulus(p),
    ensures match r {
        Ok(b) => packing::pack_error(a@, width(p), 0).is_none()
            && packing::packed(a@, width(p), b@, 8),
        Err(e) => packing::pack_error(a@, width(p), 0) == Some(e),
    },
{
    packing::transcode_to_bytes(a, serialization_width(p))
}

/// Decode all byte bits as data, without reducing or validating against p.
pub(crate) fn deserialize_vec(b: &[u8], p: u64)
    -> (r: Result<Vec<u64>, packing::TranscodeError>)
    requires valid_modulus(p),
    ensures match r {
        Ok(a) => packing::unpack_error(b@, width(p)).is_none()
            && packing::packed(b@, 8, a@, width(p)),
        Err(e) => packing::unpack_error(b@, width(p)) == Some(e),
    },
{
    packing::transcode_from_bytes(b, serialization_width(p))
}

/// Serializing canonical words cannot fail the width or value checks.
pub proof fn canonical_serialization_succeeds(a: Seq<u64>, p: u64)
    requires valid_modulus(p), a.len() * width(p) <= usize::MAX,
        forall|i: int| 0 <= i < a.len() ==> #[trigger] a[i] < p,
    ensures packing::pack_error(a, width(p), 0).is_none(),
{
    canonical_words_fit(a, p);
    width_bounds(p);
    assert((a.len() * width(p) + 7) / 8 <= usize::MAX) by (nonlinear_arith)
        requires 0 <= a.len() * width(p) <= usize::MAX;
}

proof fn aligned_bits(count: int, w: int)
    requires count >= 0, count % 8 == 0,
    ensures (count * w) % 8 == 0,
        (count * w + 7) / 8 == count * w / 8,
{
    lemma_fundamental_div_mod(count, 8);
    let blocks = count / 8;
    assert(count * w == 8 * (blocks * w)) by (nonlinear_arith)
        requires count == 8 * blocks;
    lemma_fundamental_div_mod_converse(count * w, 8, blocks * w, 0);
    lemma_fundamental_div_mod_converse(count * w + 7, 8, blocks * w, 7);
}

/// The length API agrees with the encoded length for complete groups of eight.
pub proof fn serialization_length_matches_packing(a: Seq<u64>, p: u64, b: Seq<u8>)
    requires valid_modulus(p), a.len() % 8 == 0, packing::packed(a, width(p), b, 8),
    ensures b.len() == a.len() * width(p) / 8,
{
    aligned_bits(a.len() as int, width(p));
}

/// Stream decoding recovers the original prefix. Byte-aligned data has no extras.
pub proof fn serialization_round_trip(a: Seq<u64>, p: u64, b: Seq<u8>, result: Seq<u64>)
    requires valid_modulus(p), packing::valid_words(a, width(p)),
        packing::packed(a, width(p), b, 8), packing::packed(b, 8, result, width(p)),
    ensures
        result.len() >= a.len(),
        result.take(a.len() as int) == a,
        (a.len() * width(p)) % 8 == 0 ==> result == a,
        a.len() % 8 == 0 ==> result == a,
{
    width_bounds(p);
    let w = width(p);
    assert(result.len() >= a.len()) by (nonlinear_arith)
        requires 1 <= w, b.len() == (a.len() * w + 7) / 8,
            result.len() == (b.len() * 8 + w - 1) / w;
    let prefix = result.take(a.len() as int);
    assert forall|j: int| 0 <= j < a.len() implies
        packing::word_at(b, 8, j * w, #[trigger] prefix[j], w) by {
        assert(packing::word_at(b, 8, j * w, packing::Word::value(result[j]), w));
    }
    packing::proofs::packing_round_trip(a, w, b, prefix);
    if a.len() % 8 == 0 {
        aligned_bits(a.len() as int, w);
    }
    if (a.len() * w) % 8 == 0 {
        assert(result.len() == a.len()) by (nonlinear_arith)
            requires 1 <= w, (a.len() * w) % 8 == 0,
                b.len() == (a.len() * w + 7) / 8,
                result.len() == (b.len() * 8 + w - 1) / w;
        assert(result =~= a);
    }
}

}
