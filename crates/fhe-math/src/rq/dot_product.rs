//! Polynomial dot products with bounded accumulation and reusable storage.

use super::{Context, Ntt, Poly};
use crate::{Error, Result};
use std::sync::Arc;
use zeroize::Zeroize;

struct ClearAccumulator<'a>(&'a mut [u128]);

impl Drop for ClearAccumulator<'_> {
    fn drop(&mut self) {
        self.0.zeroize();
    }
}

struct Prepared {
    ctx: Arc<Context>,
    count: usize,
    public: bool,
}

// Visit each input once to count and validate it. Keep error precedence:
// empty inputs, unequal lengths, incompatible contexts, then lazy operands.
fn prepare<'a, 'b>(
    p: impl Iterator<Item = &'a Poly<Ntt>>,
    q: impl Iterator<Item = &'b Poly<Ntt>>,
) -> Result<Prepared> {
    let mut ctx = None;
    let mut left = 0;
    let mut right = 0;
    let mut compatible = true;
    let mut public = true;
    let mut lazy = false;
    for poly in p {
        let expected = ctx.get_or_insert_with(|| poly.ctx.clone());
        compatible &= *expected == poly.ctx;
        public &= poly.allow_variable_time_computations;
        lazy |= poly.has_lazy_coefficients;
        left += 1;
    }
    for poly in q {
        compatible &= ctx.as_ref().is_some_and(|expected| expected == &poly.ctx);
        public &= poly.allow_variable_time_computations;
        lazy |= poly.has_lazy_coefficients;
        right += 1;
    }
    if left == 0 || right == 0 {
        return Err(Error::EmptyDotProduct);
    }
    if left != right {
        return Err(Error::DotProductLengthMismatch { left, right });
    }
    if !compatible {
        return Err(Error::PolynomialContextMismatch);
    }
    if lazy {
        return Err(Error::LazyDotProductOperand);
    }
    Ok(Prepared {
        ctx: ctx.unwrap(),
        count: left,
        public,
    })
}

/// Reusable storage for dot products in a fixed polynomial context.
///
/// The accumulator and reduction schedule are allocated once. Use
/// [`Self::dot_product_into`] to reuse the result allocation as well.
/// Accumulator contents are zeroized after every computation, including panic
/// unwinding, so dropping the workspace releases only cleared storage. The
/// workspace contains no borrowed inputs.
pub struct DotProductWorkspace {
    ctx: Arc<Context>,
    accumulator: Vec<u128>,
    counts: Vec<u128>,
    limits: Vec<u128>,
    min_limit: u128,
}

impl DotProductWorkspace {
    /// Allocate workspace for `ctx`. Equivalent separately allocated contexts
    /// are accepted by the computation methods.
    #[must_use]
    pub fn new(ctx: &Arc<Context>) -> Self {
        let limits: Vec<_> = ctx
            .q
            .iter()
            .map(|q| 1u128 << (2 * q.leading_zeros()))
            .collect();
        let min_limit = limits.iter().copied().min().unwrap_or(1);
        Self {
            ctx: ctx.clone(),
            accumulator: vec![0; ctx.q.len() * ctx.degree],
            counts: vec![1; ctx.q.len()],
            limits,
            min_limit,
        }
    }

    /// Compute a dot product, allocating only its result polynomial.
    /// Returns an error for empty or unequal inputs, foreign contexts, or lazy
    /// coefficients. Iterators must yield the same operands when cloned.
    pub fn dot_product<'a, 'b, I, J>(&mut self, p: I, q: J) -> Result<Poly<Ntt>>
    where
        I: Iterator<Item = &'a Poly<Ntt>> + Clone,
        J: Iterator<Item = &'b Poly<Ntt>> + Clone,
    {
        let prepared = prepare(p.clone(), q.clone())?;
        if prepared.ctx != self.ctx {
            return Err(Error::PolynomialContextMismatch);
        }
        let mut out = Poly::zero(&self.ctx);
        self.compute(p, q, &mut out, prepared);
        Ok(out)
    }

    /// Overwrite `out` with a dot product without allocating coefficient
    /// buffers. Inputs and output must match this workspace's context.
    /// Returns the same input errors as [`Self::dot_product`]; validation
    /// errors leave `out` unchanged. Timing permission is recomputed from
    /// all inputs on every call. Iterators must yield the same operands
    /// when cloned.
    pub fn dot_product_into<'a, 'b, I, J>(&mut self, p: I, q: J, out: &mut Poly<Ntt>) -> Result<()>
    where
        I: Iterator<Item = &'a Poly<Ntt>> + Clone,
        J: Iterator<Item = &'b Poly<Ntt>> + Clone,
    {
        let prepared = prepare(p.clone(), q.clone())?;
        if prepared.ctx != self.ctx || out.ctx != self.ctx {
            return Err(Error::PolynomialContextMismatch);
        }
        self.compute(p, q, out, prepared);
        Ok(())
    }

    fn compute<'a, 'b>(
        &mut self,
        p: impl Iterator<Item = &'a Poly<Ntt>>,
        q: impl Iterator<Item = &'b Poly<Ntt>>,
        out: &mut Poly<Ntt>,
        prepared: Prepared,
    ) {
        // The buffer starts at zero and this guard clears it on every exit,
        // including an unwind from an input iterator. No secret scratch survives.
        let acc = ClearAccumulator(self.accumulator.as_mut_slice());
        self.counts.fill(1);
        let needs_reduction = prepared.count as u128 > self.min_limit;
        if needs_reduction {
            for (pi, qi) in p.zip(q) {
                fma(
                    acc.0,
                    pi.coefficients.as_slice().unwrap(),
                    qi.coefficients.as_slice().unwrap(),
                );
                for (((row, modulus), count), limit) in acc
                    .0
                    .chunks_exact_mut(self.ctx.degree)
                    .zip(self.ctx.q.iter())
                    .zip(self.counts.iter_mut())
                    .zip(self.limits.iter())
                {
                    *count += 1;
                    if *count == *limit {
                        for value in row {
                            *value = if prepared.public {
                                (unsafe { modulus.reduce_u128_vt(*value) }) as u128
                            } else {
                                modulus.reduce_u128(*value) as u128
                            };
                        }
                        *count = 1;
                    }
                }
            }
        } else {
            for (pi, qi) in p.zip(q) {
                fma(
                    acc.0,
                    pi.coefficients.as_slice().unwrap(),
                    qi.coefficients.as_slice().unwrap(),
                );
            }
        }
        // With q < 2^(64-leading_zeros(q)), at most 2^(2*leading_zeros(q))
        // products fit in u128 when starting from zero. Periodic reduction uses
        // one fewer product to reserve room for the previous reduced residue.
        for ((mut output, row), modulus) in out
            .coefficients
            .outer_iter_mut()
            .zip(acc.0.chunks_exact(self.ctx.degree))
            .zip(self.ctx.q.iter())
        {
            for (value, accumulated) in output.iter_mut().zip(row) {
                *value = if prepared.public {
                    unsafe { modulus.reduce_u128_vt(*accumulated) }
                } else {
                    modulus.reduce_u128(*accumulated)
                };
            }
        }
        out.allow_variable_time_computations = prepared.public;
        out.has_lazy_coefficients = false;
    }
}

/// Compute a dot product of NTT polynomials, allocating workspace and output.
/// Returns an error for empty or unequal inputs, incompatible contexts, or lazy
/// coefficients. Iterators must yield the same operands when cloned. For
/// repeated calls, use [`DotProductWorkspace`] to reuse storage.
pub fn dot_product<'a, 'b, I, J>(p: I, q: J) -> Result<Poly<Ntt>>
where
    I: Iterator<Item = &'a Poly<Ntt>> + Clone,
    J: Iterator<Item = &'b Poly<Ntt>> + Clone,
{
    let prepared = prepare(p.clone(), q.clone())?;
    let mut workspace = DotProductWorkspace::new(&prepared.ctx);
    let mut out = Poly::zero(&prepared.ctx);
    workspace.compute(p, q, &mut out, prepared);
    Ok(out)
}
/// Computes the Fused-Mul-Add operation `out[i] += x[i] * y[i]`
///
/// Uses safe slice chunk APIs (Rust 1.88+) to process elements in chunks of 16
/// for better performance through loop unrolling, while maintaining safety.
fn fma(out: &mut [u128], x: &[u64], y: &[u64]) {
    let n = out.len();
    assert_eq!(x.len(), n);
    assert_eq!(y.len(), n);

    // Process complete chunks of 16 elements using safe chunk APIs
    let (out_chunks, out_remainder) = out.as_chunks_mut::<16>();
    let (x_chunks, x_remainder) = x.as_chunks::<16>();
    let (y_chunks, y_remainder) = y.as_chunks::<16>();

    for ((out_chunk, x_chunk), y_chunk) in out_chunks.iter_mut().zip(x_chunks).zip(y_chunks) {
        for i in 0..16 {
            out_chunk[i] += (x_chunk[i] as u128) * (y_chunk[i] as u128);
        }
    }

    // Process any remaining elements
    for ((out_elem, x_elem), y_elem) in out_remainder.iter_mut().zip(x_remainder).zip(y_remainder) {
        *out_elem += (*x_elem as u128) * (*y_elem as u128);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rq::traits::TryConvertFrom;

    fn context() -> Arc<Context> {
        Context::new_arc(&[1153, 4611686018326724609], 16).unwrap()
    }

    #[test]
    fn workspace_reuse_preserves_bounds_and_resets_policy() {
        let ctx = context();
        let values: Vec<u64> = ctx
            .q
            .iter()
            .flat_map(|q| vec![**q - 1; ctx.degree])
            .collect();
        let worst = Poly::<Ntt>::try_convert_from(values, &ctx, false).unwrap();
        let mut workspace = DotProductWorkspace::new(&ctx);
        let mut out = Poly::zero(&ctx);
        let output_pointer = out.coefficients.as_ptr();
        let scratch_pointer = workspace.accumulator.as_ptr();
        // The large modulus permits 16 products from zero, or 15 after reduction.
        for length in [1, 15, 16, 17, 30, 31, 32, 33, 256, 4, 1] {
            for public in [true, false] {
                let mut left = vec![worst.clone(); length];
                let mut right = left.clone();
                for poly in left.iter_mut().chain(right.iter_mut()) {
                    poly.allow_variable_time_computations = true;
                }
                right.last_mut().unwrap().allow_variable_time_computations = public;
                out.has_lazy_coefficients = true;
                workspace
                    .dot_product_into(left.iter(), right.iter(), &mut out)
                    .unwrap();
                assert_eq!(out.allows_variable_time_computations(), public);
                assert!(!out.has_lazy_coefficients);
                for (row, q) in out.coefficients.outer_iter().zip(ctx.q.iter()) {
                    assert!(row.iter().all(|x| *x == q.reduce(length as u64)));
                }
                assert_eq!(out, dot_product(left.iter(), right.iter()).unwrap());
                assert_eq!(
                    out,
                    workspace.dot_product(left.iter(), right.iter()).unwrap()
                );
                assert_eq!(out.coefficients.as_ptr(), output_pointer);
                assert_eq!(workspace.accumulator.as_ptr(), scratch_pointer);
                assert!(workspace.accumulator.iter().all(|x| *x == 0));
            }
        }
    }

    #[test]
    fn validation_errors_leave_output_unchanged_and_workspace_reusable() {
        let ctx = context();
        let p = Poly::<Ntt>::random_from_seed(&ctx, [1; 32]);
        let mut out = p.clone();
        let mut workspace = DotProductWorkspace::new(&ctx);
        assert_eq!(
            workspace.dot_product_into(std::iter::empty(), std::iter::once(&p), &mut out),
            Err(Error::EmptyDotProduct)
        );
        assert_eq!(
            workspace.dot_product_into(std::iter::once(&p), std::iter::empty(), &mut out),
            Err(Error::EmptyDotProduct)
        );
        assert_eq!(
            workspace.dot_product_into([&p].into_iter(), [&p, &p].into_iter(), &mut out),
            Err(Error::DotProductLengthMismatch { left: 1, right: 2 })
        );
        let foreign_ctx = Context::new_arc(&[1153], 16).unwrap();
        let foreign = Poly::zero(&foreign_ctx);
        assert_eq!(
            workspace.dot_product_into([&p].into_iter(), [&foreign].into_iter(), &mut out),
            Err(Error::PolynomialContextMismatch)
        );
        assert_eq!(
            workspace.dot_product_into([&foreign].into_iter(), [&foreign].into_iter(), &mut out),
            Err(Error::PolynomialContextMismatch)
        );
        let mut lazy = p.clone();
        lazy.has_lazy_coefficients = true;
        assert_eq!(
            workspace.dot_product_into([&p].into_iter(), [&lazy].into_iter(), &mut out),
            Err(Error::LazyDotProductOperand)
        );
        assert_eq!(
            dot_product([&lazy].into_iter(), [&p].into_iter()),
            Err(Error::LazyDotProductOperand)
        );
        assert_eq!(out, p);
        assert_eq!(
            out.allows_variable_time_computations(),
            p.allows_variable_time_computations()
        );
        let mut foreign_output = foreign.clone();
        assert_eq!(
            workspace.dot_product_into([&p].into_iter(), [&p].into_iter(), &mut foreign_output),
            Err(Error::PolynomialContextMismatch)
        );
        assert_eq!(foreign_output, foreign);
        let separate = context();
        let equivalent = Poly::<Ntt>::random_from_seed(&separate, [1; 32]);
        workspace
            .dot_product_into([&equivalent].into_iter(), [&p].into_iter(), &mut out)
            .unwrap();
        assert_eq!(out, &p * &p);
        assert!(workspace.accumulator.iter().all(|x| *x == 0));
    }

    #[test]
    fn inputs_are_visited_once_for_validation_and_once_for_arithmetic() {
        use std::cell::Cell;
        let ctx = context();
        let operands = vec![Poly::<Ntt>::random_from_seed(&ctx, [3; 32]); 4];
        let left_visits = Cell::new(0);
        let right_visits = Cell::new(0);
        dot_product(
            operands
                .iter()
                .filter(|_| true)
                .inspect(|_| left_visits.set(left_visits.get() + 1)),
            operands
                .iter()
                .filter(|_| true)
                .inspect(|_| right_visits.set(right_visits.get() + 1)),
        )
        .unwrap();
        assert_eq!(left_visits.get(), 8);
        assert_eq!(right_visits.get(), 8);
    }

    #[test]
    #[expect(clippy::panic, reason = "exercise scratch cleanup during unwinding")]
    fn accumulator_is_cleared_when_an_iterator_panics() {
        let ctx = context();
        let p = Poly::<Ntt>::random_from_seed(&ctx, [9; 32]);
        let mut workspace = DotProductWorkspace::new(&ctx);
        let mut out = p.clone();
        let prepared = prepare(std::iter::once(&p), std::iter::once(&p)).unwrap();
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            workspace.compute(
                std::iter::once(&p).chain(std::iter::once_with(|| panic!("iterator failed"))),
                [&p, &p].into_iter(),
                &mut out,
                prepared,
            );
        }));
        assert!(result.is_err());
        assert!(workspace.accumulator.iter().all(|x| *x == 0));
        assert_eq!(out, p);
        workspace
            .dot_product_into([&p].into_iter(), [&p].into_iter(), &mut out)
            .unwrap();
        assert_eq!(out, &p * &p);
    }
}
