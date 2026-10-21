//! Bit-buffer lemmas and the packing/decoding round-trip theorem.

use super::*;

verus! {

pub(super) proof fn bit_append(x: u128, y: u64, n: u32, w: u32, k: u32)
    requires n < 64, 1 <= w <= 64, k < 128,
        x >> n == 0,
        y & !(u64::MAX >> ((64 - w) as u32)) == 0,
    ensures
        bit(x | ((y as u128) << n), k as int) ==
            if k < n { bit(x, k as int) }
            else if k < n + w { bit(y as u128, (k - n) as int) }
            else { false },
{
    assert(((x | ((y as u128) << n)) >> k) & 1 ==
        if k < n { (x >> k) & 1 }
        else if k < n + w { ((y as u128) >> ((k - n) as u32)) & 1 }
        else { 0 }) by (bit_vector)
        requires n < 64, 1 <= w <= 64, k < 128, x >> n == 0,
            y & !(u64::MAX >> ((64 - w) as u32)) == 0;
}

pub(super) proof fn bit_shift(x: u128, w: u32, k: u32)
    requires 1 <= w <= 64, k < 128,
    ensures bit(x >> w, k as int) == (k + w < 128 && bit(x, (k + w) as int)),
{
    assert(((x >> w) >> k) & 1 == if k + w < 128 { (x >> ((k + w) as u32)) & 1 } else { 0 })
        by (bit_vector);
}

pub(super) proof fn bit_mask(x: u128, w: u32, k: u32)
    requires 1 <= w <= 64, k < 64,
    ensures bit((x & ((u64::MAX >> ((64 - w) as u32)) as u128)) as u64 as u128, k as int)
        == (k < w && bit(x, k as int)),
{
    assert((((x & ((u64::MAX >> ((64 - w) as u32)) as u128)) as u64 as u128) >> k) & 1
        == if k < w { (x >> k) & 1 } else { 0 }) by (bit_vector)
        requires 1 <= w <= 64, k < 64;
}

pub(super) proof fn buffer_bound<T: Word>(a: Seq<T>, iw: int, offset: int, pool: u128, n: u32)
    requires n < 128, buffer(a, iw, offset, pool, n as int),
    ensures pool >> n == 0,
{
    if pool >> n != 0 {
        // Equality of all bits implies equality of the machine words.
        assert(forall|k: u32| #![auto] k < 128 ==> (pool >> k) & 1 ==
            ((pool & (u128::MAX >> ((128 - n) as u32))) >> k) & 1) by {
            assert forall|k: u32| #![auto] k < 128 implies (pool >> k) & 1 ==
                ((pool & (u128::MAX >> ((128 - n) as u32))) >> k) & 1 by {
                if k >= n { assert(!bit(pool, k as int)); }
                assert((pool >> k) & 1 == ((pool & (u128::MAX >> ((128 - n) as u32))) >> k) & 1)
                    by (bit_vector)
                    requires k < 128, n < 128, k >= n ==> (pool >> k) & 1 != 1;
            }
        }
        let masked = pool & (u128::MAX >> ((128 - n) as u32));
        assert(pool == masked) by (bit_vector)
            requires forall|k: u32| #![auto] k < 128 ==> (pool >> k) & 1 == (masked >> k) & 1;
        assert(pool >> n == 0) by (bit_vector)
            requires n < 128, pool == pool & (u128::MAX >> ((128 - n) as u32));
    }
}

pub(super) proof fn append_buffer<T: Word>(a: Seq<T>, iw: u32, i: int, offset: int, pool: u128, n: u32)
    requires valid_words(a, iw as int), 1 <= iw <= 64,
        0 <= i < a.len(), n < 64, offset >= 0, offset + n == i * iw,
        buffer(a, iw as int, offset, pool, n as int),
    ensures buffer(a, iw as int, offset, pool | ((a[i].value() as u128) << n), (n + iw) as int),
{
    buffer_bound(a, iw as int, offset, pool, n);
    let y = a[i].value();
    assert forall|k: int| 0 <= k < 128 implies
        #[trigger] bit(pool | ((y as u128) << n), k)
            == (k < n + iw && stream(a, iw as int, offset + k)) by {
        bit_append(pool, y, n, iw, k as u32);
        if n <= k < n + iw {
            lemma_fundamental_div_mod_converse(offset + k, iw as int, i, k - n);
            assert(offset + k < a.len() * iw) by (nonlinear_arith)
                requires offset + n == i * iw, k < n + iw, i < a.len(), iw > 0;
        }
    }
}

pub(super) proof fn take_buffer<T: Word>(a: Seq<T>, iw: int, offset: int, pool: u128, n: u32, ow: u32)
    requires 1 <= ow <= 64, ow <= n < 128, buffer(a, iw, offset, pool, n as int),
    ensures buffer(a, iw, offset + ow, pool >> ow, (n - ow) as int),
        word_at(a, iw, offset, (pool & ((u64::MAX >> ((64 - ow) as u32)) as u128)) as u64, ow as int),
        fits((pool & ((u64::MAX >> ((64 - ow) as u32)) as u128)) as u64, ow as int),
{
    assert forall|k: int| 0 <= k < 128 implies #[trigger] bit(pool >> ow, k)
        == (k < n - ow && stream(a, iw, offset + ow + k)) by {
        bit_shift(pool, ow, k as u32);
    }
    assert forall|k: int| 0 <= k < 64 implies
        #[trigger] bit((pool & ((u64::MAX >> ((64 - ow) as u32)) as u128)) as u64 as u128, k)
        == (k < ow && stream(a, iw, offset + k)) by {
        bit_mask(pool, ow, k as u32);
    }
    assert(((pool & ((u64::MAX >> ((64 - ow) as u32)) as u128)) as u64)
        & !(u64::MAX >> ((64 - ow) as u32)) == 0) by (bit_vector);
}

pub(super) proof fn fits_widen(x: u64, small: u32, large: u32)
    requires 1 <= small <= large <= 64, fits(x, small as int),
    ensures fits(x, large as int),
{
    assert(x & !(u64::MAX >> ((64 - large) as u32)) == 0) by (bit_vector)
        requires 1 <= small <= large <= 64, x & !(u64::MAX >> ((64 - small) as u32)) == 0;
}

/// Repartitioning a stream preserves every bit, including the zero extension.
pub proof fn packed_stream<T: Word, U: Word>(a: Seq<T>, iw: int, b: Seq<U>, ow: int, k: int)
    requires 1 <= iw <= 64, 1 <= ow <= 64, packed(a, iw, b, ow),
    ensures stream(b, ow, k) == stream(a, iw, k),
{
    if 0 <= k < b.len() * ow {
        lemma_fundamental_div_mod(k, ow);
        lemma_mod_bound(k, ow);
        assert(0 <= k / ow < b.len()) by (nonlinear_arith)
            requires 0 <= k < b.len() * ow, ow > 0;
        assert(word_at(a, iw, (k / ow) * ow, b[k / ow].value(), ow));
        assert(bit(b[k / ow].value() as u128, k % ow) == stream(a, iw, k));
    } else {
        assert(b.len() * ow >= a.len() * iw) by (nonlinear_arith)
            requires ow > 0, b.len() == (a.len() * iw + ow - 1) / ow;
    }
}

pub(super) proof fn fits_bit(x: u64, width: u32, k: u32)
    requires 1 <= width <= 64, width <= k < 64, fits(x, width as int),
    ensures !bit(x as u128, k as int),
{
    assert(((x as u128) >> k) & 1 == 0) by (bit_vector)
        requires 1 <= width <= 64, width <= k < 64,
            x & !(u64::MAX >> ((64 - width) as u32)) == 0;
}

pub(super) proof fn equal_bits(x: u64, y: u64)
    requires forall|k: int| 0 <= k < 64 ==> #[trigger] bit(x as u128, k) == bit(y as u128, k),
    ensures x == y,
{
    let xx = x as u128;
    let yy = y as u128;
    assert(forall|k: u32| #![auto] k < 64 ==> (xx >> k) & 1 == (yy >> k) & 1) by {
        assert forall|k: u32| #![auto] k < 64 implies (xx >> k) & 1 == (yy >> k) & 1 by {
            assert(bit(x as u128, k as int) == bit(y as u128, k as int));
            assert((xx >> k) & 1 == (yy >> k) & 1) by (bit_vector)
                requires (((xx >> k) & 1) == 1) == (((yy >> k) & 1) == 1);
        }
    }
    assert(x == y) by (bit_vector)
        requires xx == x as u128, yy == y as u128, forall|k: u32| #![auto] k < 64 ==> (xx >> k) & 1 == (yy >> k) & 1;
}

/// Exact decoding of packed words recovers the original sequence.
pub proof fn packing_round_trip(a: Seq<u64>, width: int, b: Seq<u8>, result: Seq<u64>)
    requires valid_words(a, width), packed(a, width, b, 8), decoded(b, width, a.len() as int, result),
    ensures result == a,
{
    assert forall|j: int| 0 <= j < a.len() implies #[trigger] result[j] == a[j] by {
        assert(word_at(b, 8, j * width, result[j], width));
        assert(fits(a[j].value(), width));
        assert forall|k: int| 0 <= k < 64 implies #[trigger] bit(result[j] as u128, k) == bit(a[j] as u128, k) by {
            if k < width {
                packed_stream(a, width, b, 8, j * width + k);
                lemma_fundamental_div_mod_converse(j * width + k, width, j, k);
                assert(j * width + k < a.len() * width) by (nonlinear_arith)
                    requires j < a.len(), k < width, width > 0;
            } else {
                fits_bit(a[j], width as u32, k as u32);
            }
        }
        equal_bits(result[j], a[j]);
    }
    assert(result =~= a);
}

/// Packing supplies the exact decoder's byte length and zero padding.
pub proof fn packing_is_decodable(a: Seq<u64>, width: int, b: Seq<u8>)
    requires valid_words(a, width), packed(a, width, b, 8), b.len() * 8 <= usize::MAX,
    ensures exact_error(b, width, a.len() as int).is_none(),
{
    let bits = a.len() * width;
    assert(bits <= b.len() * 8) by (nonlinear_arith)
        requires b.len() == (bits + 7) / 8;
    if bits % 8 != 0 {
        lemma_fundamental_div_mod(bits, 8);
        assert(b.len() > 0 && (b.len() - 1) * 8 == bits - bits % 8) by (nonlinear_arith)
            requires bits >= 0, bits % 8 != 0, b.len() == (bits + 7) / 8;
        let last = b[b.len() - 1];
        let r = (bits % 8) as u32;
        let last128 = last as u128;
        assert forall|k: u32| #![auto] r <= k < 8 implies (last128 >> k) & 1 != 1 by {
            assert(word_at(a, width, (b.len() - 1) * 8, last.value(), 8));
            assert(!bit(last as u128, k as int));
        }
        assert(last >> r == 0) by (bit_vector)
            requires last128 == last as u128, 0 < r < 8, forall|k: u32| #![auto] r <= k < 8 ==> (last128 >> k) & 1 != 1;
    }
}

}
