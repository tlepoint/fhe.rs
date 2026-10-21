//! Shared bit-packing implementation, verified directly with Verus.

// Verus generates enum accessors in the ghost Rust view without doc comments.
#![cfg_attr(verus_keep_ghost, allow(missing_docs))]

#[cfg(verus_keep_ghost)]
use vstd::arithmetic::div_mod::*;
use vstd::prelude::*;
#[cfg(verus_keep_ghost)]
#[path = "packing/proofs.rs"]
pub mod proofs;
#[cfg(verus_keep_ghost)]
use proofs::*;

verus! {

/// Invalid input to a bit-transcoding operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum TranscodeError {
    /// Word widths must lie in 1..=64.
    InvalidWidth(usize),
    /// At least one word would be truncated at the requested width.
    ValueTooWide,
    /// A computed input or output length cannot be represented by `usize`.
    LengthOverflow,
    /// The byte length does not match the requested word count and width.
    LengthMismatch,
    /// Unused bits in the last byte must be zero.
    NonZeroPadding,
}

/// Bit k of a machine word, with zero outside the word.
pub open spec fn bit(x: u128, k: int) -> bool {
    0 <= k < 128 && (x >> (k as u32)) & 1 == 1
}

/// Storage words supported by the packing loop.
pub trait Word: Sized {
    /// The unsigned value represented by this storage word.
    spec fn value(self) -> u64;
    /// The storage word's bit width.
    spec fn width() -> int;
    /// Read the storage word without changing its value.
    fn read(&self) -> (r: u64)
        ensures r == self.value(), Self::width() == 8 || Self::width() == 64;
    /// Store a value that fits the storage word's width.
    fn write(x: u64) -> (r: Self)
        requires fits(x, Self::width()),
        ensures r.value() == x;
}

impl Word for u64 {
    open spec fn value(self) -> u64 { self }
    open spec fn width() -> int { 64 }
    #[inline]
    fn read(&self) -> u64 { *self }
    #[inline]
    fn write(x: u64) -> u64 { x }
}

impl Word for u8 {
    open spec fn value(self) -> u64 { self as u64 }
    open spec fn width() -> int { 8 }
    #[inline]
    fn read(&self) -> u64 { *self as u64 }
    #[inline]
    fn write(x: u64) -> u8 {
        proof {
            assert(x <= 255) by (bit_vector)
                requires x & !(u64::MAX >> 56) == 0;
        }
        x as u8
    }
}

/// Whether a value fits a supported bit width without truncation.
pub open spec fn fits(x: u64, width: int) -> bool {
    1 <= width <= 64 && x & !(u64::MAX >> ((64 - width) as u32)) == 0
}

/// A supported width and a sequence of words that all fit it.
pub open spec fn valid_words<T: Word>(a: Seq<T>, width: int) -> bool {
    1 <= width <= 64 && forall|i: int| 0 <= i < a.len() ==> fits(#[trigger] a[i].value(), width)
}

/// The LSB-first input bit stream, extended with zeros beyond its end.
pub open spec fn stream<T: Word>(a: Seq<T>, width: int, k: int) -> bool {
    0 <= k < a.len() * width && bit(a[k / width].value() as u128, k % width)
}

/// One output word contains the specified slice of the zero-extended input stream.
pub open spec fn word_at<T: Word>(a: Seq<T>, iw: int, offset: int, x: u64, ow: int) -> bool {
    forall|k: int| 0 <= k < 64 ==> #[trigger] bit(x as u128, k)
        == (k < ow && stream(a, iw, offset + k))
}

/// Exact output length and bit contents for a stream repartitioned into words.
pub open spec fn packed<T: Word, U: Word>(a: Seq<T>, iw: int, b: Seq<U>, ow: int) -> bool {
    b.len() == (a.len() * iw + ow - 1) / (ow as int)
    && forall|j: int| 0 <= j < b.len() ==> word_at(a, iw, j * ow, #[trigger] b[j].value(), ow)
}

spec fn buffer<T: Word>(a: Seq<T>, iw: int, offset: int, pool: u128, n: int) -> bool {
    forall|k: int| 0 <= k < 128 ==> #[trigger] bit(pool, k)
        == (k < n && stream(a, iw, offset + k))
}

/// One buffer loop for byte packing, byte unpacking, and word transcoding.
#[expect(clippy::indexing_slicing, reason = "indices and bounds are proved by Verus")]
fn transcode<T: Word, U: Word>(a: &[T], iw: usize, ow: usize, out: &mut Vec<U>)
    requires valid_words(a@, iw as int), 1 <= iw <= 64, 1 <= ow <= 64,
        ow <= U::width() <= 64,
        a.len() * iw <= usize::MAX,
        old(out).len() + (a.len() * iw + ow - 1) / (ow as int) <= usize::MAX,
    ensures
        final(out)@.take(old(out)@.len() as int) == old(out)@,
        packed(a@, iw as int, final(out)@.skip(old(out)@.len() as int), ow as int),
{
    let count = ceil_words(a.len() * iw, ow);
    let ghost start = out.len();
    let ghost prefix = out@;
    out.reserve(count);
    let mask = (u64::MAX >> (64 - ow)) as u128;
    let mut i = 0usize;
    let mut pool = 0u128;
    let mut n = 0usize;
    proof { assert(forall|k: u32| #![auto] k < 128 ==> (0u128 >> k) & 1 == 0) by (bit_vector); }
    while i < a.len()
        invariant
            valid_words(a@, iw as int), 1 <= iw <= 64, 1 <= ow <= 64,
            ow <= U::width() <= 64, a.len() * iw <= usize::MAX,
            count as int == (a.len() * iw + ow - 1) / (ow as int),
            start == prefix.len(), start + count <= usize::MAX,
            0 <= i <= a.len(), 0 <= n < ow,
            start <= out.len(), (out.len() - start) * ow + n == i * iw,
            out@.take(start as int) == prefix,
            forall|j: int| start <= j < out.len() ==>
                word_at(a@, iw as int, (j - start) * ow, #[trigger] out@[j].value(), ow as int),
            buffer(a@, iw as int, (out.len() - start) * ow, pool, n as int),
            mask == (u64::MAX >> (64 - ow)) as u128,
        decreases a.len() - i,
    {
        let value = a[i].read();
        proof { append_buffer(a@, iw as u32, i as int, (out.len() - start) * ow, pool, n as u32); }
        pool |= (value as u128) << n;
        n += iw;
        i += 1;
        proof { assert((i - 1) * iw + iw == i * iw) by (nonlinear_arith); }
        while n >= ow
            invariant
                valid_words(a@, iw as int), 1 <= iw <= 64, 1 <= ow <= 64,
                ow <= U::width() <= 64, a.len() * iw <= usize::MAX,
                count as int == (a.len() * iw + ow - 1) / (ow as int),
                start == prefix.len(), start + count <= usize::MAX,
                0 <= i <= a.len(), 0 <= n < 128,
                start <= out.len(), (out.len() - start) * ow + n == i * iw,
                out@.take(start as int) == prefix,
                forall|j: int| start <= j < out.len() ==>
                    word_at(a@, iw as int, (j - start) * ow, #[trigger] out@[j].value(), ow as int),
                buffer(a@, iw as int, (out.len() - start) * ow, pool, n as int),
                mask == (u64::MAX >> (64 - ow)) as u128,
            decreases n,
        {
            proof { take_buffer(a@, iw as int, (out.len() - start) * ow, pool, n as u32, ow as u32); }
            let ghost len = out.len();
            proof {
                assert(pool & mask <= u64::MAX as u128) by (bit_vector)
                    requires mask <= u64::MAX as u128;
                fits_widen((pool & mask) as u64, ow as u32, U::width() as u32); }
            out.push(U::write((pool & mask) as u64));
            pool >>= ow;
            n -= ow;
            proof {
                assert((len - start) * ow + ow == (out.len() - start) * ow) by (nonlinear_arith)
                    requires out.len() == len + 1;
            }
        }
    }
    if n > 0 {
        proof {
            // Once all input is consumed, the missing bits of the last word are zero.
            assert(buffer(a@, iw as int, (out.len() - start) * ow, pool, ow as int));
            take_buffer(a@, iw as int, (out.len() - start) * ow, pool, ow as u32, ow as u32);
            buffer_bound(a@, iw as int, (out.len() - start) * ow, pool, n as u32);
            assert(pool == pool & mask) by (bit_vector)
                requires pool >> (n as u32) == 0, n < ow, 1 <= ow <= 64,
                    mask == (u64::MAX >> ((64 - ow) as u32)) as u128;
        }
        proof {
            assert(pool <= u64::MAX as u128) by (bit_vector)
                requires pool == pool & mask, mask <= u64::MAX as u128;
            fits_widen(pool as u64, ow as u32, U::width() as u32); }
        out.push(U::write(pool as u64));
    }
    proof {
        assert(out.len() - start == count) by (nonlinear_arith)
            requires 1 <= ow, 0 <= n < ow,
                (out.len() - start - if n > 0 { 1int } else { 0int }) * ow + n == a.len() * iw,
                count == (a.len() * iw + ow - 1) / (ow as int);
    }
}

/// Overflow-free ceiling division for the checked bit lengths.
fn ceil_words(bits: usize, width: usize) -> (r: usize)
    requires 1 <= width <= 64,
    ensures r as int == (bits as int + width as int - 1) / width as int,
        r <= bits,
{
    let q = bits / width;
    let r = bits % width;
    proof {
        lemma_fundamental_div_mod(bits as int, width as int);
        assert(q as int + (if r > 0 { 1int } else { 0int }) <= bits) by (nonlinear_arith)
            requires width >= 1, q >= 0, r >= 0, bits == q * width + r;
        assert((bits as int + width as int - 1) / width as int ==
            q as int + if r > 0 { 1int } else { 0int }) by (nonlinear_arith)
            requires 1 <= width, 0 <= r < width, bits == q * width + r;
    }
    q + if r > 0 { 1 } else { 0 }
}

/// Input validation error, in the public API's precedence order.
pub open spec fn words_error(a: Seq<u64>, width: int) -> Option<TranscodeError> {
    if !(1 <= width <= 64) { Some(TranscodeError::InvalidWidth(width as usize)) }
    else if a.len() * width > usize::MAX { Some(TranscodeError::LengthOverflow) }
    else if !valid_words(a, width) { Some(TranscodeError::ValueTooWide) }
    else { None }
}

/// Packing error, including overflow when appending to an existing prefix.
pub open spec fn pack_error(a: Seq<u64>, width: int, start: int) -> Option<TranscodeError> {
    if words_error(a, width).is_some() { words_error(a, width) }
    else if start + (a.len() * width + 7) / 8 > usize::MAX { Some(TranscodeError::LengthOverflow) }
    else { None }
}

fn validate_width(width: usize) -> (r: Result<(), TranscodeError>)
    ensures match r {
        Ok(()) => 1 <= width <= 64,
        Err(e) => !(1 <= width <= 64) && e == TranscodeError::InvalidWidth(width),
    },
{
    if !(1..=64).contains(&width) {
        return Err(TranscodeError::InvalidWidth(width));
    }
    Ok(())
}

#[expect(clippy::indexing_slicing, reason = "indices and bounds are proved by Verus")]
fn validate_words(values: &[u64], width: usize) -> (r: Result<usize, TranscodeError>)
    ensures match r {
        Ok(bits) => words_error(values@, width as int).is_none() && bits == values.len() * width,
        Err(e) => words_error(values@, width as int) == Some(e),
    },
{
    validate_width(width)?;
    let bits = values.len().checked_mul(width).ok_or(TranscodeError::LengthOverflow)?;
    let mut used = 0u64;
    let mut i = 0usize;
    proof { assert(0u64 & !(u64::MAX >> ((64 - width) as u32)) == 0) by (bit_vector); }
    // Scan all words, preserving the previous aggregate validation behavior.
    while i < values.len()
        invariant 1 <= width <= 64, bits == values.len() * width, bits <= usize::MAX,
            0 <= i <= values.len(),
            fits(used, width as int) == (forall|j: int| 0 <= j < i ==> fits(#[trigger] values@[j], width as int)),
        decreases values.len() - i,
    {
        let word = values[i];
        proof {
            assert(((used | word) & !(u64::MAX >> ((64 - width) as u32)) == 0)
                == (used & !(u64::MAX >> ((64 - width) as u32)) == 0
                    && word & !(u64::MAX >> ((64 - width) as u32)) == 0)) by (bit_vector);
        }
        used |= word;
        i += 1;
    }
    proof {
        if valid_words(values@, width as int) {
            assert forall|j: int| 0 <= j < i implies fits(#[trigger] values@[j], width as int) by {
                assert(fits(values@[j].value(), width as int));
            }
        } else if fits(used, width as int) {
            assert forall|j: int| 0 <= j < values.len() implies fits(#[trigger] values@[j].value(), width as int) by {
                assert(fits(values@[j], width as int));
            }
        }
        assert(valid_words(values@, width as int) == fits(used, width as int));
    }
    if used & !(u64::MAX >> (64 - width)) != 0 {
        return Err(TranscodeError::ValueTooWide);
    }
    Ok(bits)
}

/// Pack words least-significant-bit first, padding the final byte with zeros.
/// Rejects invalid widths, overflowing lengths, and words that would truncate.
pub fn transcode_to_bytes(a: &[u64], nbits: usize) -> (r: Result<Vec<u8>, TranscodeError>)
    ensures match r {
        Ok(b) => pack_error(a@, nbits as int, 0).is_none() && packed(a@, nbits as int, b@, 8),
        Err(e) => pack_error(a@, nbits as int, 0) == Some(e),
    },
{
    let mut out = Vec::new();
    transcode_to_bytes_into(a, nbits, &mut out)?;
    proof { assert(out@.skip(0) =~= out@); }
    Ok(out)
}

/// Appends packed `nbits`-bit numbers to `out`, padding the final byte with
/// zeros. Each call starts at a byte boundary and preserves the existing
/// prefix. Invalid input leaves `out` unchanged, including its allocation.
pub fn transcode_to_bytes_into(a: &[u64], nbits: usize, out: &mut Vec<u8>) -> (r: Result<(), TranscodeError>)
    ensures match r {
        Ok(()) => pack_error(a@, nbits as int, old(out).len() as int).is_none()
            && final(out)@.take(old(out).len() as int) == old(out)@
            && packed(a@, nbits as int, final(out)@.skip(old(out).len() as int), 8),
        Err(e) => pack_error(a@, nbits as int, old(out).len() as int) == Some(e)
            && final(out)@ == old(out)@,
    },
{
    let bits = validate_words(a, nbits)?;
    let nbytes = ceil_words(bits, 8);
    out.len().checked_add(nbytes).ok_or(TranscodeError::LengthOverflow)?;
    transcode(a, nbits, 8, out);
    Ok(())
}

/// Width or input bit-length error for stream decoding.
pub open spec fn unpack_error(b: Seq<u8>, width: int) -> Option<TranscodeError> {
    if !(1 <= width <= 64) { Some(TranscodeError::InvalidWidth(width as usize)) }
    else if b.len() * 8 > usize::MAX { Some(TranscodeError::LengthOverflow) }
    else { None }
}

proof fn byte_words(b: Seq<u8>)
    ensures valid_words(b, 8),
{
    assert forall|i: int| 0 <= i < b.len() implies fits(#[trigger] b[i].value(), 8) by {
        let x = b[i];
        assert((x as u64) & !(u64::MAX >> 56) == 0) by (bit_vector);
    }
}

/// Interpret every input bit as a stream of `nbits`-bit words,
/// least-significant bit first. Zero-extend a partial final word. Byte padding
/// is interpreted as data; use [`transcode_from_bytes_exact`] when the original
/// word count is known.
pub fn transcode_from_bytes(b: &[u8], nbits: usize) -> (r: Result<Vec<u64>, TranscodeError>)
    ensures match r {
        Ok(a) => unpack_error(b@, nbits as int).is_none() && packed(b@, 8, a@, nbits as int),
        Err(e) => unpack_error(b@, nbits as int) == Some(e),
    },
{
    validate_width(nbits)?;
    let bits = b.len().checked_mul(8).ok_or(TranscodeError::LengthOverflow)?;
    let count = ceil_words(bits, nbits);
    let mut out = Vec::with_capacity(count);
    proof { byte_words(b@); }
    transcode(b, 8, nbits, &mut out);
    proof { assert(out@.skip(0) =~= out@); }
    Ok(out)
}

/// Whether the unused high bits of the final byte are zero.
pub open spec fn padding_zero(b: Seq<u8>, bits: int) -> bool {
    bits % 8 == 0 || (b.len() > 0 && b[b.len() - 1] >> ((bits % 8) as u32) == 0)
}

/// Exact decoding error, including length mismatch and nonzero padding.
pub open spec fn exact_error(b: Seq<u8>, width: int, count: int) -> Option<TranscodeError> {
    if !(1 <= width <= 64) { Some(TranscodeError::InvalidWidth(width as usize)) }
    else if count * width > usize::MAX { Some(TranscodeError::LengthOverflow) }
    else if b.len() != (count * width + 7) / 8 { Some(TranscodeError::LengthMismatch) }
    else if !padding_zero(b, count * width) { Some(TranscodeError::NonZeroPadding) }
    else { unpack_error(b, width) }
}

/// Exactly count words decoded from the byte stream at the requested width.
pub open spec fn decoded(b: Seq<u8>, width: int, count: int, a: Seq<u64>) -> bool {
    a.len() == count && forall|j: int| 0 <= j < count ==>
        word_at(b, 8, j * width, #[trigger] a[j], width)
}

/// Decode exactly `count` packed words, rejecting extra/missing bytes and
/// nonzero padding. This is the inverse of [`transcode_to_bytes`] when the
/// caller supplies the original word count, including for empty inputs.
#[expect(clippy::indexing_slicing, reason = "the final byte's existence is proved by Verus")]
pub fn transcode_from_bytes_exact(b: &[u8], nbits: usize, count: usize)
    -> (r: Result<Vec<u64>, TranscodeError>)
    ensures match r {
        Ok(a) => exact_error(b@, nbits as int, count as int).is_none()
            && decoded(b@, nbits as int, count as int, a@),
        Err(e) => exact_error(b@, nbits as int, count as int) == Some(e),
    },
{
    validate_width(nbits)?;
    let bits = count.checked_mul(nbits).ok_or(TranscodeError::LengthOverflow)?;
    let nbytes = ceil_words(bits, 8);
    if b.len() != nbytes {
        return Err(TranscodeError::LengthMismatch);
    }
    if bits % 8 != 0 {
        proof {
            assert(b.len() > 0) by (nonlinear_arith)
                requires bits > 0, b.len() == (bits as int + 7) / 8;
        }
        if b[b.len() - 1] >> (bits % 8) != 0 {
            return Err(TranscodeError::NonZeroPadding);
        }
    }
    let mut words = transcode_from_bytes(b, nbits)?;
    proof {
        assert(count <= words.len()) by (nonlinear_arith)
            requires 1 <= nbits, bits == count * nbits,
                b.len() == (bits as int + 7) / 8,
                words.len() == (b.len() * 8 + nbits - 1) / (nbits as int);
    }
    let ghost all = words@;
    words.truncate(count);
    proof {
        assert forall|j: int| 0 <= j < count implies word_at(b@, 8, j * nbits, #[trigger] words@[j], nbits as int) by {
            assert(word_at(b@, 8, j * nbits, all[j].value(), nbits as int));
        }
        assert(exact_error(b@, nbits as int, count as int).is_none());
    }
    Ok(words)
}

/// Word-transcoding error, checking input validity before the output width.
pub open spec fn transcode_error(a: Seq<u64>, iw: int, ow: int) -> Option<TranscodeError> {
    if words_error(a, iw).is_some() { words_error(a, iw) }
    else if !(1 <= ow <= 64) { Some(TranscodeError::InvalidWidth(ow as usize)) }
    else { None }
}

/// Transcode `input_nbits`-bit words into `output_nbits`-bit words, LSB first.
/// Reject truncation and invalid widths. Zero-extend a partial final word;
/// the output does not retain the original bit length.
pub fn transcode_bidirectional(a: &[u64], input_nbits: usize, output_nbits: usize)
    -> (r: Result<Vec<u64>, TranscodeError>)
    ensures match r {
        Ok(b) => transcode_error(a@, input_nbits as int, output_nbits as int).is_none()
            && packed(a@, input_nbits as int, b@, output_nbits as int),
        Err(e) => transcode_error(a@, input_nbits as int, output_nbits as int) == Some(e),
    },
{
    let bits = validate_words(a, input_nbits)?;
    validate_width(output_nbits)?;
    let count = ceil_words(bits, output_nbits);
    let mut out = Vec::with_capacity(count);
    transcode(a, input_nbits, output_nbits, &mut out);
    proof { assert(out@.skip(0) =~= out@); }
    Ok(out)
}

} // verus!

#[cfg(test)]
mod tests {
    #[test]
    fn ceiling_lengths_cover_machine_boundaries() {
        for width in 1usize..=64 {
            for bits in [0, 1, 7, 8, 63, 64, 65, usize::MAX - 1, usize::MAX] {
                assert_eq!(super::ceil_words(bits, width), bits.div_ceil(width));
            }
        }
    }
}
