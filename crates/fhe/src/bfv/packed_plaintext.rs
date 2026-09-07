//! Compact in-memory storage of plaintext NTT coefficients.

use super::{Parameters, Plaintext};
use fhe_math::rq::{Ntt, Poly};
use zeroize::Zeroize;

/// Compact in-memory plaintext storage with bit-packed NTT residues.
///
/// Each residue uses the bit width of its ciphertext modulus. For example, two
/// 36-bit moduli occupy 72 bits per coefficient instead of 128. Packing is
/// exact: [`Self::unpack`] restores the same plaintext, level, and
/// local timing permission without any transforms.
/// [`super::evaluation::DotProductScalarWorkspace::dot_product_scalar_packed_iter`]
/// consumes this storage directly. This type is not a wire format.
///
/// Packing still requires the initial plaintext NTT. It trades additional
/// decoding work during dot products for lower retained database memory.
#[derive(Clone)]
pub struct PackedPlaintext {
    pub(crate) par: Parameters,

    pub(crate) level: usize,
    pub(crate) public: bool,
    pub(crate) coefficients: Vec<u8>,
}

/// A borrowed row of packed NTT residues for scalar dot products.
#[derive(Clone, Copy)]
pub struct PackedPlaintextView<'a> {
    pub(crate) par: &'a Parameters,
    pub(crate) level: usize,
    pub(crate) public: bool,
    pub(crate) coefficients: &'a [u8],
}

impl<'a> From<&'a PackedPlaintext> for PackedPlaintextView<'a> {
    fn from(pt: &'a PackedPlaintext) -> Self {
        Self {
            par: &pt.par,
            level: pt.level,
            public: pt.public,
            coefficients: &pt.coefficients,
        }
    }
}

/// Contiguous packed NTT storage for plaintexts at the same parameters and
/// level.
///
/// A single allocation avoids per-plaintext allocator size-class padding. The
/// collection retains each row's timing permission; its borrowed views can be
/// consumed directly by scalar dot products. Stored bytes are cleared on drop.
pub struct PackedPlaintextBatch {
    par: Parameters,
    level: usize,
    row_bytes: usize,
    coefficients: Vec<u8>,
    public: Vec<bool>,
}

impl PackedPlaintextBatch {
    /// Reserve coefficient storage for `capacity` plaintexts at this level.
    pub fn with_capacity(par: &Parameters, level: usize, capacity: usize) -> crate::Result<Self> {
        let ctx = par.context_at_level(level)?;
        let row_bytes = ctx
            .moduli_operators()
            .iter()
            .map(|q| q.serialization_length(par.degree()))
            .sum::<usize>();
        let mut coefficients =
            Vec::with_capacity(capacity.saturating_mul(row_bytes).saturating_add(8));
        coefficients.resize(8, 0);
        Ok(Self {
            par: par.clone(),
            level,
            row_bytes,
            coefficients,
            public: Vec::with_capacity(capacity),
        })
    }

    /// Pack borrowed plaintexts with these shared parameters and level. Empty
    /// input creates an empty batch. Context validation is fallible.
    pub fn try_from_iter<'a>(
        par: &Parameters,
        level: usize,
        plaintexts: impl IntoIterator<Item = &'a Plaintext>,
    ) -> crate::Result<Self> {
        let mut batch = Self::with_capacity(par, level, 0)?;
        batch.try_extend(plaintexts)?;
        Ok(batch)
    }

    /// Pack and append a plaintext. Parameter/level errors leave storage
    /// unchanged.
    pub fn push(&mut self, pt: &Plaintext) -> crate::Result<()> {
        let ctx = self.par.context_at_level(self.level)?;
        pt.validate_for_context(&self.par, self.level, ctx)?;
        self.reserve_rows(1);
        self.append_validated(pt);
        Ok(())
    }

    /// Snapshot and validate all inputs before appending. A returned error
    /// leaves the batch unchanged, including its allocation. Each iterator
    /// is used once.
    pub fn try_extend<'a>(
        &mut self,
        plaintexts: impl IntoIterator<Item = &'a Plaintext>,
    ) -> crate::Result<()> {
        let plaintexts: Vec<_> = plaintexts.into_iter().collect();
        let ctx = self.par.context_at_level(self.level)?;
        for pt in &plaintexts {
            pt.validate_for_context(&self.par, self.level, ctx)?;
        }
        self.reserve_rows(plaintexts.len());
        for pt in plaintexts {
            self.append_validated(pt);
        }
        Ok(())
    }

    fn reserve_rows(&mut self, additional: usize) {
        let needed = self
            .coefficients
            .len()
            .saturating_add(additional.saturating_mul(self.row_bytes));
        if needed > self.coefficients.capacity() {
            let mut replacement =
                Vec::with_capacity(needed.max(self.coefficients.capacity().saturating_mul(2)));
            replacement.extend_from_slice(&self.coefficients);
            self.coefficients.as_mut_slice().zeroize();
            self.coefficients = replacement;
        }
        self.public.reserve(additional);
    }

    fn append_validated(&mut self, pt: &Plaintext) {
        self.coefficients.truncate(self.coefficients.len() - 8);
        append_coefficients(pt, &mut self.coefficients);
        self.coefficients.resize(self.coefficients.len() + 8, 0);
        self.public
            .push(pt.poly_ntt.allows_variable_time_computations());
    }

    /// Clear all rows and their bytes, retaining capacity for reuse.
    pub fn clear(&mut self) {
        self.zeroize();
        self.coefficients.truncate(8);
        self.public.clear();
    }

    /// Number of stored plaintexts.
    #[must_use]
    pub fn len(&self) -> usize {
        self.public.len()
    }

    /// Whether the collection contains no plaintexts.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.public.is_empty()
    }

    /// Packed coefficient bytes, including the single eight-byte padding
    /// suffix.
    #[must_use]
    pub fn size_bytes(&self) -> usize {
        self.coefficients.len()
    }

    /// Borrow one row, without allocating or unpacking coefficients.
    #[must_use]
    pub fn get(&self, index: usize) -> Option<PackedPlaintextView<'_>> {
        let &public = self.public.get(index)?;
        let start = index * self.row_bytes;
        Some(PackedPlaintextView {
            par: &self.par,
            level: self.level,
            public,
            coefficients: self.coefficients.get(start..start + self.row_bytes + 8)?,
        })
    }

    /// Iterate over borrowed rows without allocating or unpacking coefficients.
    #[must_use]
    pub fn iter(&self) -> PackedPlaintextIter<'_> {
        PackedPlaintextIter {
            batch: self,
            indices: 0..self.len(),
        }
    }
}

/// Borrowed iteration over the validated rows of a [`PackedPlaintextBatch`].
#[derive(Clone)]
pub struct PackedPlaintextIter<'a> {
    batch: &'a PackedPlaintextBatch,
    indices: std::ops::Range<usize>,
}

impl<'a> Iterator for PackedPlaintextIter<'a> {
    type Item = PackedPlaintextView<'a>;
    fn next(&mut self) -> Option<Self::Item> {
        self.batch.get(self.indices.next()?)
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.indices.size_hint()
    }
}
impl DoubleEndedIterator for PackedPlaintextIter<'_> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.batch.get(self.indices.next_back()?)
    }
}
impl ExactSizeIterator for PackedPlaintextIter<'_> {}
impl std::iter::FusedIterator for PackedPlaintextIter<'_> {}
impl<'a> IntoIterator for &'a PackedPlaintextBatch {
    type Item = PackedPlaintextView<'a>;
    type IntoIter = PackedPlaintextIter<'a>;
    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

impl Zeroize for PackedPlaintextBatch {
    fn zeroize(&mut self) {
        self.coefficients.as_mut_slice().zeroize();
    }
}

impl Drop for PackedPlaintextBatch {
    fn drop(&mut self) {
        self.zeroize();
    }
}

// Emit packed words eight bytes at a time. The row length is a multiple
// of eight, so the final fragment always contains complete bytes.
fn append_row(row: &[u64], bits: usize, out: &mut Vec<u8>) {
    debug_assert!(row.len().is_multiple_of(8));
    debug_assert!((1..=64).contains(&bits));
    out.reserve(row.len() * bits / 8);
    let mut buffer = 0u128;
    let mut buffered_bits = 0;
    for &value in row {
        buffer |= u128::from(value) << buffered_bits;
        buffered_bits += bits;
        if buffered_bits >= 64 {
            out.extend_from_slice(&(buffer as u64).to_le_bytes());
            buffer >>= 64;
            buffered_bits -= 64;
        }
    }
    if buffered_bits != 0 {
        out.extend_from_slice(
            (buffer as u64)
                .to_le_bytes()
                .get(..buffered_bits / 8)
                .unwrap(),
        );
    }
    buffer.zeroize();
}

fn append_coefficients(pt: &Plaintext, coefficients: &mut Vec<u8>) {
    for (row, q) in pt
        .poly_ntt
        .coefficients()
        .outer_iter()
        .zip(pt.poly_ntt.ctx().moduli_operators())
    {
        append_row(
            row.as_slice().unwrap(),
            (64 - q.leading_zeros()) as usize,
            coefficients,
        );
    }
}

impl Zeroize for PackedPlaintext {
    fn zeroize(&mut self) {
        self.coefficients.as_mut_slice().zeroize();
    }
}

impl Drop for PackedPlaintext {
    fn drop(&mut self) {
        self.zeroize();
    }
}

impl From<&Plaintext> for PackedPlaintext {
    fn from(pt: &Plaintext) -> Self {
        let mut coefficients = Vec::with_capacity(
            pt.poly_ntt
                .ctx()
                .moduli_operators()
                .iter()
                .map(|q| q.serialization_length(pt.par.degree()))
                .sum::<usize>()
                + 8,
        );
        append_coefficients(pt, &mut coefficients);
        // Permit bounded unaligned word reads at the end of the final row.
        coefficients.resize(coefficients.len() + 8, 0);
        Self {
            par: pt.par.clone(),

            level: pt.level(),
            public: pt.poly_ntt.allows_variable_time_computations(),
            coefficients,
        }
    }
}

impl PackedPlaintext {
    /// Restore the canonical plaintext without an NTT.
    #[must_use]
    pub fn unpack(&self) -> Plaintext {
        let ctx = self.par.context_at_level(self.level).unwrap();
        let degree = self.par.degree();
        let mut coefficients = vec![0; degree * ctx.moduli().len()];
        let mut offset = 0;
        for (row, q) in coefficients
            .chunks_exact_mut(degree)
            .zip(ctx.moduli_operators())
        {
            let bits = (64 - q.leading_zeros()) as usize;
            let mask = u64::MAX >> (64 - bits);
            for (i, value) in row.iter_mut().enumerate() {
                let bit = i * bits;
                let byte = offset + bit / 8;
                let word = u64::from_le_bytes(
                    self.coefficients
                        .get(byte..byte + 8)
                        .unwrap()
                        .try_into()
                        .unwrap(),
                );
                *value = (word >> (bit % 8)) & mask;
                if bits + bit % 8 > 64 {
                    *value |=
                        u64::from(*self.coefficients.get(byte + 8).unwrap()) << (64 - bit % 8);
                    *value &= mask;
                }
            }
            offset += bits * degree / 8;
        }
        Plaintext {
            par: self.par.clone(),

            poly_ntt: Poly::<Ntt>::from_rns_residues_with_timing(
                ndarray::Array2::from_shape_vec(
                    (ctx.moduli().len(), self.par.degree()),
                    coefficients,
                )
                .unwrap(),
                ctx,
                (self.public).then(|| crate::VariableTime::new(crate::PublicData::assert_public())),
            )
            .unwrap(),
        }
    }

    /// Packed coefficient bytes, including eight padding bytes; excludes
    /// metadata.
    #[must_use]
    pub fn size_bytes(&self) -> usize {
        self.coefficients.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bfv::{
        Ciphertext, Encoding, ParametersBuilder, evaluation::DotProductScalarWorkspace,
    };
    use crate::{PublicData, VariableTime};

    #[test]
    fn word_packing_matches_byte_packing_for_every_width() {
        use rand::{RngExt, SeedableRng};
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(8);
        for bits in 1..=64 {
            let mask = u64::MAX >> (64 - bits);
            for degree in [8, 16, 128] {
                for row in [
                    vec![0; degree],
                    vec![mask; degree],
                    (0..degree).map(|_| rng.random::<u64>() & mask).collect(),
                ] {
                    let mut actual = vec![1, 2, 3];
                    append_row(&row, bits, &mut actual);
                    let mut expected = vec![1, 2, 3];
                    fhe_util::transcode_to_bytes_into(&row, bits, &mut expected).unwrap();
                    assert_eq!(actual, expected);
                }
            }
        }
    }

    #[test]
    fn contiguous_storage_grows_validates_and_preserves_each_rows_permission() -> crate::Result<()>
    {
        let par = ParametersBuilder::new()
            .degree(16)
            .plaintext_modulus(17_u64)
            .ciphertext_modulus_bits([36, 62])
            .build()?;
        let mut data = PackedPlaintextBatch::with_capacity(&par, 0, 2)?;
        assert!(data.is_empty());
        assert_eq!(data.iter().len(), 0);
        let mut originals = Vec::new();
        for i in 0..5 {
            let mut pt = Plaintext::encode(&par, &[i as u64, 16, 0, 1][..], Encoding::Polynomial)?;
            if i % 2 == 0 {
                pt.poly_ntt
                    .allow_variable_time_computations(VariableTime::new(
                        PublicData::assert_public(),
                    ));
            }
            data.push(&pt)?;
            originals.push(pt);
        }
        assert_eq!(data.len(), 5);
        assert_eq!(data.size_bytes(), 5 * par.degree() * (36 + 62) / 8 + 8);
        let pt = originals.first().unwrap();
        let ct = Ciphertext {
            par: par.clone(),
            seed: None,
            c: vec![pt.poly_ntt.clone(); 2],
            level: 0,
        };
        let mut workspace = DotProductScalarWorkspace::new(&par, 0)?;
        for (view, original) in data.iter().rev().zip(originals.iter().rev()) {
            assert_eq!(
                view.public,
                original.poly_ntt.allows_variable_time_computations()
            );
            assert_eq!(
                workspace
                    .dot_product_scalar_packed_iter(std::iter::once(&ct), std::iter::once(view))?,
                ct.multiply_plaintext(original).unwrap()
            );
        }
        assert_eq!(
            workspace.dot_product_scalar_packed_iter(
                std::iter::repeat_n(&ct, 3),
                data.iter().step_by(2)
            )?,
            workspace.dot_product_scalar_iter(
                std::iter::repeat_n(&ct, 3),
                originals.iter().step_by(2)
            )?
        );
        let previous = data.coefficients.clone();
        let lower = Plaintext::encode_at_level(&par, &[3u64][..], Encoding::Polynomial, 1)?;
        assert!(data.push(&lower).is_err());
        let foreign = ParametersBuilder::new()
            .noise_variance(11)
            .degree(16)
            .plaintext_modulus(17_u64)
            .ciphertext_modulus_bits([36, 62])
            .build()?;
        assert!(
            data.push(&Plaintext::encode(
                &foreign,
                &[3u64][..],
                Encoding::Polynomial
            )?)
            .is_err()
        );
        assert_eq!(data.coefficients, previous);
        assert_eq!(data.len(), 5);
        data.zeroize();
        assert_eq!(data.len(), 5);
        assert!(data.coefficients.iter().all(|x| *x == 0));
        let zero =
            workspace.dot_product_scalar_packed_iter(std::iter::repeat_n(&ct, 5), data.iter())?;
        assert!(
            zero.iter()
                .all(|p| p.coefficients().iter().all(|x| *x == 0))
        );
        assert!(PackedPlaintextBatch::with_capacity(&par, 9, 0).is_err());
        Ok(())
    }
}
