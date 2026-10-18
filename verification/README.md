# Verifying `fhe-math` with Verus

The verified kernels in
[`crates/fhe-math/src/zq/verified.rs`](../crates/fhe-math/src/zq/verified.rs)
form compositional proofs: `const_time_cond_select`, then `reduce1`, then
`add_mod`; and `reduce1_vt`, then `add_mod_vt`. The variable-time reducer uses
either a branch or the verified `reduce1`, according to the target. Each caller
uses its helper's proved contract. `Modulus::add`, `Modulus::add_vt`, both
reduction wrappers, and the other conditional-selection callers use these
same functions in ordinary builds. Verus checks the executable
implementations and their contracts; the proof is not a separate copy of the
algorithms. Normal Cargo builds erase the proof code through `vstd`.

The included [`verified/lazy_reduce.rs`](../crates/fhe-math/src/zq/verified/lazy_reduce.rs)
also proves every `lazy_reduce*` kernel, including the vector loop. The
`Modulus` methods delegate to these verified implementations.

## Run the proof

Use the pinned [Verus 0.2026.10.04.426d8b0 release](https://github.com/verus-lang/verus/releases/tag/release/0.2026.10.04.426d8b0)
and its required Rust toolchain:

```sh
rustup toolchain install 1.98.1 --profile minimal
```

Download and unpack the release archive for your platform, following the
[Verus installation instructions](https://github.com/verus-lang/verus/blob/426d8b01e7ffb910a36c15e1f870ff985bb61eef/INSTALL.md).
Add its directory to `PATH`, or invoke that directory's `verus` executable
directly. From the repository root, run:

```sh
verus --version
verus --crate-type lib --edition=2024 --no-cheating crates/fhe-math/src/zq/verified.rs
```

The expected result is `28 verified, 0 errors` on macOS or an AVX2 target, and
`27 verified, 0 errors` elsewhere, including arithmetic lemmas, local assertions,
and the vector loop. On x86 Linux, check both target-dependent implementations with:

```sh
verus --crate-type lib --edition=2024 --no-cheating -C target-feature=-avx2 crates/fhe-math/src/zq/verified.rs
verus --crate-type lib --edition=2024 --no-cheating -C target-feature=+avx2 crates/fhe-math/src/zq/verified.rs
```

`--no-cheating` rejects assumptions, admits, and unchecked
external bodies in our proof source. Verus, its bundled standard-library
specifications, Z3, Rust, and the compiler remain trusted infrastructure.

This checks the arithmetic module and its submodule directly, without asking
Verus to analyze the crate's unrelated dependencies or algorithms. It needs neither Cargo dependency
compilation nor `protoc`. It is not a `cargo verus verify -p fhe-math` integration
or a verification of the whole crate. The pinned `vstd` dependency is used by
ordinary Rust builds; the standalone verifier uses its bundled `vstd` instead.
The current Cargo pin is `0.0.0-2026-10-11-0230`, while the verifier release is
from October 4. Rerun both verification and Cargo checks when updating either
pin; the proof and macro-erased production builds must both work.
The verifier's compiler is separate from the library's Rust 1.91.1 MSRV.

## Contracts and proof order

1. **`const_time_cond_select(on_true, on_false, cond)`** has no preconditions.
   Its result equals `if cond { on_true } else { on_false }` for all `u64`
   values and both Boolean values. The original implementation is retained:
   `-(cond as i64) as u64` builds the mask, and XOR/AND select the word. A
   bitvector proof establishes the selection contract, including the signed
   mask's conversion to `u64`.
2. **`reduce1(x, p)`** assumes `0 < p < 2^63` and `x < 2*p`. It calls the
   verified selector with `x`, `x.wrapping_sub(p)`, and `x < p`. The selector's
   contract and a remainder lemma establish `result < p` and
   `result = x mod p`. This bound deliberately includes `p = 1` and the NTT's
   calls with twice a 62-bit modulus. The wrapping subtraction is intentional;
   when `x < p`, the selector discards its result.
3. **`add_mod(a, b, p)`** assumes `2 <= p < 2^62`, `a < p`, and `b < p`.
   Its body is just `reduce1(a + b, p)`. Verus establishes that the sum cannot
   overflow and satisfies the reduction preconditions, then obtains a
   canonical result equal to `(a + b) mod p` from the reduction contract.
4. **`reduce1_vt(x, p)`** has the same contract as `reduce1`. On macOS or
   targets with AVX2, its body is `if x >= p { x - p } else { x }`, with a
   remainder lemma proving correctness and Verus checking subtraction safety.
   On other targets it calls `reduce1` and uses that verified contract. The
   target selection is unchanged from the original implementation.
5. **`add_mod_vt(a, b, p)`** has the same contract as `add_mod`. Its body is
   `unsafe { reduce1_vt(a + b, p) }`. Verus proves the sum cannot overflow and
   satisfies the reduction preconditions, then uses the reducer's contract to
   establish the canonical modular sum. The existing `unsafe` timing boundary
   and `const fn` signatures are retained for both variable-time functions.

### Lazy reduction

All lazy reducers establish `0 <= result < 2*p` and
`result mod p = input mod p`. They intentionally permit a noncanonical result
in `[p, 2*p)`. Verus also checks intermediate arithmetic, shift counts, casts,
and slice accesses. The scalar kernels retain `const fn` support and the
original multiplication/shift algorithms.

Let `B = 2^64` and `mu = hi*B + lo`. The `valid_barrett` contract requires
`2 <= p < 2^62` and `mu*p <= B^2 < (mu+1)*p`, equivalently
`mu = floor(2^128/p)`. These are the reciprocal limbs constructed by
`Modulus::new`; the bounds are assumptions on the kernel parameters, not
additional runtime calculations.

| Kernel | Additional assumptions | Proof |
| --- | --- | --- |
| `lazy_reduce` | Any `u64` input; `valid_barrett` | The two limb products compute `floor(a*mu/B^2)`. |
| `lazy_reduce_u128` | Any `u128` input; `valid_barrett` | The four limb products compute the same quotient without overflowing their intermediate sums. |
| `lazy_reduce_opt` | Any `u64` input; `valid_opt` | The shift quotient gives `0 <= a-q*p < 2*p`. No reciprocal is needed. |
| `lazy_reduce_opt_u128` | `a < p^2`; both predicates | The optimized quotient omits a low product; the modulus condition bounds that approximation error. |
| `lazy_reduce_vec` | `valid_barrett`; also `valid_opt` for its optimized branch | Loop invariants preserve unreduced elements and establish the scalar contract for every processed element; length is unchanged. |

`valid_opt(p, s)` specifies the actual optimization condition: `s` is the
leading-zero count, `t = 2^s`, and
`(t^3 + 1)*B < t^3*(t+1)*p`. This is the integer inequality in
`primes::supports_opt`; primality is not required. The proof derives
`hi = t` and `p^2 + B*lo < B^2` from this condition and the reciprocal bounds.
Those inequalities justify dropping the low product in the optimized
128-bit quotient. A separate shift lemma proves that `a << s` does not discard
bits for the allowed inputs.

The shared Barrett remainder lemma proves that the approximate quotient never
exceeds `floor(a/p)` and is at most one below it. It then obtains congruence by
subtracting a multiple of `p`. The limb lemma `hi + lo <= B` supplies the bound
needed for the full-width intermediate sum. The vector loop invokes the
scalar contracts, with a const Boolean selecting the implementation so the
existing dispatch remains outside the loop.

All arithmetic in contracts is interpreted over mathematical integers. The
proofs cover the full stated input domains, including even and composite
moduli, and preserve `const fn` use.

The reduction and addition proofs assume valid inputs. `Modulus::new` checks
the modulus; the addition and reduction wrappers retain their debug
assertions. The optimized lazy wrappers now also assert `supports_opt` in
debug builds, and the 128-bit variant retains its `a < p^2` assertion.
Verus does not verify `Modulus::new`, the BigUint implementation of
`primes::supports_opt`, the wrappers, or their callers, and does not insert
runtime checks. Correct reciprocal limbs and optimization metadata remain
the boundary between the constructor and these kernel contracts. Although
other operations reuse the verified kernels, their contracts have not been
proved. These proofs do not establish constant-time machine code, NTT/RNS
correctness, or FHE security.

## Tests and CI

Ordinary Rust tests exercise both selection outcomes on boundary and
alternating-bit words, exhaust reduction for moduli 1 through 64 and addition
for moduli 2 through 64, and check reduction/addition boundaries at the maximum
supported widths, for both the ordinary and variable-time variants. They also
check the selection, reduction, and addition kernels remain usable in `const`
expressions. Existing randomized arithmetic tests exercise their production
callers, including rejection of
invalid operands by `add_vt` in debug builds:

```sh
cargo test -p fhe-math zq::tests
```

Lazy-reduction tests add exhaustive small inputs, full-width inputs including
`u128::MAX`, powers of two, even/composite moduli, and both sides of the strict
optimization threshold at each modulus width. They exercise all four scalar
kernels in const expressions, both vector implementations (including empty
slices), and invalid optimized inputs in debug builds. Randomized tests use
Rust's `%` as an independent congruence oracle for all four scalar kernels.

The [Verus workflow](../.github/workflows/verus.yml) installs the pinned
checksum-checked Linux release and runs the proof on pull requests and pushes,
with AVX2 disabled and enabled to cover both `reduce1_vt` implementations.
The ordinary test, formatting, Clippy, and MSRV jobs continue to check Rust
builds independently. When changing the kernel, rerun both the proof and these
tests. To sanity-check the contracts, removing the mask's negation must fail
the selector proof; changing the reduction comparison `x < p` to `x <= p`
must fail the reduction proof; dropping `b` from `reduce1(a + b, p)` must fail
the addition proof. For the variable-time path, changing `x >= p` to `x > p`
must fail the branch reducer's proof at `x = p`; dropping `b` from
`reduce1_vt(a + b, p)` must fail the variable-time addition proof.
For lazy reduction, removing the high-product contribution, corrupting either
optimized quotient shift, or leaving vector entries unreduced must fail
verification as well.
