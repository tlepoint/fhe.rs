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
also proves every `lazy_reduce*` kernel, including the vector loop, and the
`primes::supports_opt` eligibility check. The `Modulus` methods delegate to
these verified implementations.

[`verified/reduce.rs`](../crates/fhe-math/src/zq/verified/reduce.rs) completes
the `reduce*` family: canonical unsigned and signed reductions, optimized
reductions, all six vector variants, and the native NTT's `reduce3` helpers.
Both ordinary and variable-time paths use these shared implementations.

[`verified/add.rs`](../crates/fhe-math/src/zq/verified/add.rs),
[`verified/sub.rs`](../crates/fhe-math/src/zq/verified/sub.rs), and
[`verified/neg.rs`](../crates/fhe-math/src/zq/verified/neg.rs) each contain the
corresponding scalar kernels and in-place vector loops, including ordinary
and variable-time variants. All six vector operations are verified.

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

The expected result is `77 verified, 0 errors` on macOS or an AVX2 target, and
`76 verified, 0 errors` elsewhere, including arithmetic lemmas, local assertions,
and vector loops. On x86 Linux, check both target-dependent implementations with:

```sh
verus --crate-type lib --edition=2024 --no-cheating -C target-feature=-avx2 crates/fhe-math/src/zq/verified.rs
verus --crate-type lib --edition=2024 --no-cheating -C target-feature=+avx2 crates/fhe-math/src/zq/verified.rs
```

`--no-cheating` rejects assumptions, admits, and unchecked
external bodies in our proof source. Verus, its bundled standard-library
specifications, Z3, Rust, and the compiler remain trusted infrastructure.

This checks the arithmetic module and its submodules directly, without asking
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

`valid_opt(p, s)` includes the modulus and bit-length bounds together with
`optimization_condition(p, s)`: `s` is the leading-zero count, `t = 2^s`, and
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

### Optimization eligibility

`primes::supports_opt(p)` has no preconditions. For every `u64`, Verus proves
that its result is exactly `optimization_condition(p, p.leading_zeros())`.
For `2 <= p < 2^62`, it also proves that the result is equivalent to
`valid_opt(p, p.leading_zeros())`, connecting the public check to the optimized
reducers' precondition. A recursive lemma proves the bit-length bounds from
the standard-library specification of `leading_zeros`.

The executable check uses bounded `u128` arithmetic instead of `BigUint`.
For `s > 0`, dividing the original strict inequality by `t^3` gives
`(t+1)*p > B + floor(B/t^3)`. The correction is `2^(64-3*s)` when `3*s < 64`
and zero otherwise. Verus proves that this correction is the floor, that the
comparison is equivalent, and that all shifts and arithmetic are safe. The
public API path is retained by re-exporting the verified function directly.

Zero now returns `false`, fixing the previous shift by 64. Inputs with the
top bit set still return `false`; eligible 63-bit inputs remain supported even
though `Modulus` itself only accepts moduli below `2^62`. The condition does
not require primality.

### Canonical reduction

Every canonical reducer proves `0 <= result < p` and
`result = input mod p`, with modulo interpreted over mathematical integers.
For signed inputs this is the nonnegative Euclidean remainder, matching
Rust's `rem_euclid`, including at `i64::MIN`.

| Kernels (including each `_vt` variant) | Assumptions | Composition |
| --- | --- | --- |
| `reduce`, `reduce_u128` | `valid_barrett`; any input of the indicated width | The corresponding lazy reducer produces a value below `2*p`; `reduce1` or `reduce1_vt` makes it canonical. |
| `reduce_opt` | `valid_opt`; any `u64` input | Optimized lazy reduction followed by one subtraction. |
| `reduce_opt_u128` | `valid_barrett`, `valid_opt`, and `a < p^2` | Optimized 128-bit lazy reduction followed by one subtraction. |
| `reduce_i64` | `valid_barrett`; any `i64` input | Add `p*2^64`, then apply `reduce_u128`. |
| `reduce3` | `0 < p < 2^62`, `p_twice = 2*p`, and `a < 4*p` | Reduce first modulo `2*p`, then modulo `p`. |
| `reduce_vec`, `reduce_vec_new`, `reduce_vec_i64` | `valid_barrett`; unrestricted entries | Apply the appropriate canonical scalar contract to every element. |

The signed helper proves that the original signed shift and addition
`((p as i128) << 64) + (a as i128)` fit in `i128`, remain nonnegative, and
preserve the residue on conversion to `u128`. It never negates the input.
The NTT helpers use the identity `(a mod (2*p)) mod p = a mod p`; the full
NTT and construction of its `p_twice` field remain outside the proof.

All six canonical vector loops are checked, including bounds, termination,
length, element order, and canonical residues. In-place loops preserve the
unprocessed suffix; allocating loops reserve the input length and establish
the contract for each pushed value. The read-only inputs remain borrowed
immutably. The `Modulus` wrappers retain architecture dispatch around the
verified loops. Scalar kernels remain `const fn`; variable-time APIs retain
their `unsafe` timing boundary. The private signed scalar and optimized lazy
`Modulus` wrappers are now test conveniences: production vector/canonical
kernels compose those scalar helpers directly inside the verified module.

### Addition, subtraction, and negation

These kernels require `2 <= p < 2^62` and canonical inputs below `p`; they
do not require primality, Barrett reciprocals, or optimization metadata.
`add_mod` and `add_mod_vt` prove `result < p` and
`result = (a+b) mod p`. `sub_mod` and `sub_mod_vt` prove `result < p` and
`result = (a-b) mod p`. `neg_mod` and `neg_mod_vt` prove `result < p` and
`result = (-a) mod p`. These specifications use Euclidean modulo over
mathematical integers, including when the difference is negative.

The executable scalar formulas remain `reduce1(a+b, p)`,
`reduce1(a+p-b, p)`, and `reduce1(p-a, p)`, using `reduce1_vt` for the
variable-time versions. Verus proves that their additions and subtractions
fit, that the intermediate
values are below `2*p`, and that adding `p` preserves the residue. The final
reduction contract establishes canonicality, including `a-b = 0` and `-0`.

`add_vec`, `add_vec_vt`, `sub_vec`, and `sub_vec_vt` additionally require equal
lengths. Their loop invariants prove the corresponding scalar contract for
each updated entry while preserving the unprocessed suffix and total length.
The same properties hold for `neg_vec` and `neg_vec_vt`. The variable-time
addition and subtraction kernels retain complete blocks of 16 followed by a
tail; explicit index loops replace the original chunk iterators so every loop
and slice access is checked.
The `Modulus` wrappers retain architecture dispatch, scalar `const fn`
signatures, and the variable-time APIs' `unsafe` timing boundary.

Vector wrappers explicitly check canonical inputs in debug builds, along with
the existing addition and subtraction length checks. This also makes
variable-time addition and subtraction reject noncanonical operands as
documented: checking only the intermediate reduction bound could previously
miss some invalid operands.
There is no new public scalar subtraction API; `sub_mod_vt` is an internal
helper for the existing vector operation. The private `Modulus::neg_vt`
wrapper is retained for tests, while production vectors call the verified
scalar kernel directly. Polynomial addition, subtraction, and negation delegate
to these modular vector operations; their context and representation handling
remain outside the formal proof.

### Proof boundary

All arithmetic in contracts is interpreted over mathematical integers. The
proofs cover the full stated input domains, including even and composite
moduli, and preserve `const fn` use.

The arithmetic proofs assume valid inputs. `Modulus::new` checks
the modulus; the scalar and vector wrappers retain their documented debug
assertions. The optimized lazy wrappers now also assert `supports_opt` in
debug builds, and the 128-bit variant retains its `a < p^2` assertion.
Verus does not verify `Modulus::new`, the wrappers, or their callers, and does
not insert runtime checks. Correct reciprocal limbs and storage/use of the
proved optimization metadata remain the boundary between the constructor and
these kernel contracts. Although other operations reuse the verified kernels,
their contracts have not been proved. These proofs do not establish
constant-time machine code, NTT/RNS
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

Vector addition tests use an independent `u128` remainder oracle, covering
sums below, equal to, and above the modulus, including `2*p - 2`. They check
empty inputs, exact blocks, tails around multiples of 16, and sentinel-guarded
subslices, with moduli through `2^62 - 1`. Debug tests reject noncanonical
operands in either input and mismatched lengths. The randomized vector test
also uses the independent remainder oracle for both timing variants.

Lazy-reduction tests add exhaustive small inputs, full-width inputs including
`u128::MAX`, powers of two, even/composite moduli, and both sides of the strict
optimization threshold at each modulus width. They exercise all four scalar
kernels in const expressions, both vector implementations (including empty
slices), and invalid optimized inputs in debug builds. Randomized tests use
Rust's `%` as an independent congruence oracle for all four scalar kernels.

Canonical reduction tests cover every scalar and vector variant, including
`u128::MAX`, `i64::MIN`, negative multiples of the modulus, and both sides of
the optimization threshold. They exhaust small unsigned and signed inputs,
check all twelve new scalar kernels in const expressions, and compare vector
results to `%` or `i128::rem_euclid` independently. Empty, singleton, and
lengths around loop boundaries are included; sentinels check in-place slice
boundaries. Optimized wrappers reject unsupported moduli and `a >= p^2` in
debug builds. Native NTT tests cover `reduce3` and `reduce3_vt` at all four
modulus boundaries and reject `a = 4*p` in debug builds.

Subtraction and negation tests exhaust canonical operands for moduli 2 through
64, cover the largest 62-bit modulus and even/composite moduli, and check all
four scalar kernels in const expressions. Independent `i128::rem_euclid`
oracles cover negative differences and zero negation. Vector tests cover
empty inputs, exact blocks, tails around multiples of 16, and sentinel-guarded
subslices. Debug tests reject noncanonical operands in either subtraction
input, invalid negation inputs, and mismatched lengths. Existing randomized
scalar and vector tests use the same independent signed-remainder oracle:

```sh
cargo test -p fhe-math subtraction_negation
```

Eligibility tests evaluate the original NFL inequality with `BigUint` as an
independent oracle. They exhaust all 16-bit inputs, cover every power-of-two
boundary and optimization threshold, check the exact equality at `3*2^61`,
and compare random inputs across the full `u64` domain:

```sh
cargo test -p fhe-math zq::primes::tests::supports_opt
```

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
Dropping the right-hand operand from a vector addition update must fail its
loop invariant, including both the block and tail paths of `add_vec_vt`.
For lazy reduction, removing the high-product contribution, corrupting either
optimized quotient shift, or leaving vector entries unreduced must fail
verification as well.
Changing the strict comparison in `supports_opt` to `>=` must fail at
`p = 3*2^61`, where the original NFL inequality is an equality.
Canonical reduction checks also reject omitting the final subtraction,
subtracting the signed input instead of adding it in `lift_i64`, or pushing
zero instead of the computed residue in an allocating vector loop.
Subtraction and negation checks reject reversing the subtraction operands,
omitting the reduction of `p - a` at `a = 0`, and corrupting either the
16-element block update or the tail update in `sub_vec_vt`.
