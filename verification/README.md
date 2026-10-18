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
It also contains `reduce1` and both target-specific `reduce1_vt` implementations,
re-exported by `verified.rs` for the other arithmetic kernels and callers.
Both ordinary and variable-time paths use these shared implementations.

[`verified/add.rs`](../crates/fhe-math/src/zq/verified/add.rs),
[`verified/sub.rs`](../crates/fhe-math/src/zq/verified/sub.rs), and
[`verified/neg.rs`](../crates/fhe-math/src/zq/verified/neg.rs) each contain the
corresponding scalar kernels and in-place vector loops, including ordinary
and variable-time variants. All six vector operations are verified.

[`verified/mul.rs`](../crates/fhe-math/src/zq/verified/mul.rs) proves ordinary,
optimized, and Shoup multiplication, including their variable-time variants,
all six multiplication vector loops, and both Shoup precomputation loops.

[`verified/center.rs`](../crates/fhe-math/src/zq/verified/center.rs) proves scalar
and vector centering, including the signed range, congruence, and midpoint
behavior for odd and even moduli.

[`verified/pow.rs`](../crates/fhe-math/src/zq/verified/pow.rs) proves modular
exponentiation by composing multiplication contracts while tracking the
processed binary prefix of the exponent.

[`verified/random.rs`](../crates/fhe-math/src/zq/verified/random.rs) proves the
`random_vec` and `random_into` loops against an abstract sampler contract,
including canonical bounds, length, and exact correspondence to the sequence
of samples drawn. The adapter to `rand::Uniform` is a trusted boundary.

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

The expected result is `122 verified, 0 errors` on macOS or an AVX2 target, and
`121 verified, 0 errors` elsewhere, including arithmetic lemmas, local assertions,
and vector loops. On x86 Linux, check both target-dependent implementations with:

```sh
verus --crate-type lib --edition=2024 --no-cheating -C target-feature=-avx2 crates/fhe-math/src/zq/verified.rs
verus --crate-type lib --edition=2024 --no-cheating -C target-feature=+avx2 crates/fhe-math/src/zq/verified.rs
```

`--no-cheating` rejects assumptions, admits, and unchecked
external bodies in our proof source. Verus, its bundled standard-library
specifications, Z3, Rust, and the compiler remain trusted infrastructure.
For sampling, this command checks generic loops against the `ResidueSampler`
trait contract; it does not check the external `UniformSampler` implementation.

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

### Multiplication and Shoup precomputation

The ordinary and optimized scalar multiplication kernels require canonical
operands `a, b < p` and `valid_barrett(p, hi, lo)`. The optimized variants also
require `valid_opt(p, s)`. Verus proves that the `u128` product is exact and
below `p^2`, then composes the corresponding canonical reduction contract to
obtain `result < p` and `result = (a*b) mod p`. The variable-time versions
use the reducers with variable-time final subtraction.

Shoup multiplication has a wider first-operand domain: `a` may be any `u64`,
including unreduced NTT intermediates. Its second operand is canonical. With
`B = 2^64`, `valid_shoup(b, p, s)` requires `2 <= p < 2^62`, `b < p`, and
`s*p <= b*B < (s+1)*p`. Thus `s = floor(b*B/p)`. The shared `shoup` function
proves this predicate and the exact quotient, including that the left shift
and division fit in `u128` and the quotient fits in `u64`.

`lazy_mul_shoup` retains the original two products, high-word shift, and
subtraction. The proof establishes that `q = floor(a*s/B)` never exceeds
`floor(a*b/p)` and is at most one below it. Consequently, `a*b - q*p` is
nonnegative, below `2*p`, and congruent to `a*b`. This also proves that the
subtraction cannot underflow and the final cast loses no bits. `mul_shoup`
and `mul_shoup_vt` compose that result with `reduce1` and `reduce1_vt` to
obtain canonical residues. Shoup kernels need neither Barrett metadata nor
the optimization condition.

The vector kernels prove each entry's product, preserve length and element
order, and retain the unprocessed suffix. `mul_vec` and `mul_vec_vt` require
equal lengths and canonical inputs; a const parameter selects the ordinary
or optimized scalar kernel, preserving dispatch outside the loop.
`mul_shoup_vec` and its variable-time variant accept arbitrary entries in
the left slice, with equal-length canonical right operands and valid Shoup
words. `scalar_mul_vec` and its variable-time variant similarly accept any
`u64` entries and precompute their canonical scalar's Shoup word once.
`shoup_vec` and `shoup_into` prove exact precomputation and the `valid_shoup`
predicate for every output entry; `shoup_into` also requires equal lengths.

All scalar kernels remain `const fn`. The public methods retain architecture
dispatch and their debug input checks. The private variable-time scalar
wrappers and optimized variable-time reduction wrapper are test conveniences;
production vector kernels compose the verified functions directly. The
vector Shoup and scalar-multiplication documentation now states their existing
full-width left-operand support. Exponentiation composes these contracts as
described below. Other multiplication callers, including polynomial operations
and NTT butterflies, remain outside the proof.

### Centering

`center(a, p)` requires `2 <= p < 2^62` and `a < p`. It proves the exact
integer result: `a` when `2*a < p`, and `a-p` otherwise. The result satisfies
`-p <= 2*result < p` and `result mod p = a`, with Euclidean modulo. This is
the mathematical interval `[-p/2, p/2)`, or equivalently the integer interval
`[-floor(p/2), ceil(p/2))`. For even moduli, the midpoint maps to `-p/2`;
for odd moduli, `floor(p/2)` stays positive.

The executable scalar body retains `(p+1) >> 1`, signed wrapping subtraction,
and the verified conditional selector. Verus proves the threshold equals
`ceil(p/2)`, that the signed subtraction fits, and that casting a signed word
to `u64` and back preserves it. The unsigned representation of a negative
difference is intentional; the final cast recovers that signed value.

`center_vec` requires canonical input entries and proves the scalar contract
for every output, together with output length and element order. It borrows
the input immutably and reserves output space for the input length.
An explicit index loop replaces the iterator so its bounds, termination, and
processed prefix are verified. The public wrapper retains architecture
dispatch and debug canonical-input checks; scalar centering remains `const fn`.
Neither operation needs Barrett metadata, primality, or optimization eligibility.

### Exponentiation

`pow_mod(a, n, p, hi, lo)` requires `valid_barrett(p, hi, lo)`, `a < p`, and
`n < p`, retaining the public `Modulus::pow` input bounds. It proves
`result < p` and `result = a^n mod p`, where exponentiation is the standard
mathematical integer power from `vstd::arithmetic::power`. Exponent zero
returns one, including `0^0`; no primality assumption is needed.

The executable algorithm retains the two small-exponent cases and the
left-to-right square-and-multiply loop. The leading-zero lemma identifies
the highest set bit and proves the initial exponent prefix is one. At the
start of each iteration, the accumulator equals
`a^(n >> (i+1)) mod p`. Squaring doubles this processed exponent, and the
conditional multiplication adds the next bit. A bitvector lemma connects
that update to `n >> i`; a power-addition lemma and modular multiplication
identities connect it to the previously verified `mul_mod` contract.

Verus also checks that the signed loop index stays between `-1` and `60`,
every executed shift is valid, the accumulator stays canonical, and the loop
terminates. When the index reaches `-1`, the processed prefix is the whole
exponent. The public wrapper retains its debug assertions and variable-time
behavior. `inv` calls this verified exponentiation, but its inverse guarantee
and primality check remain outside the proof.

### Random vectors

`random_vec` and `random_into` share their production loops with Verus and
require a sampler for a modulus `2 <= p < 2^62`. `random_vec` returns exactly
`size` entries; `random_into` preserves its slice's length and overwrites every
entry, without requiring the original entries to be canonical. Both prove
that every output is below `p` and that the output sequence is exactly the
sequence drawn during the call, in order. Empty outputs require no draws.

The `ResidueSampler` trait describes a fixed modulus and a ghost history of
returned samples. Each `sample` call must preserve the modulus, return a value
below it, and append that value to the history. The vector loop tracks
`history = initial_history + output`; the in-place loop tracks the same
identity for the initialized prefix and preserves the unprocessed suffix.
Thus the proofs rule out omitted, extra, reordered, or substituted samples,
in addition to checking indexing, arithmetic, and loop progress. The history
erases completely in Rust builds; it allocates no runtime storage.

**Trusted sampling boundary:** `UniformSampler` in `zq/mod.rs` implements the
trait outside the standalone verified module. It holds the constructor's
`Uniform::new(0, p)` and the caller's RNG, and calls `rng.sample(distribution)`
once per trait call. Its compliance with the range and history contract is
assumed by the generic proofs, not verified by `--no-cheating`. The wrappers
and distribution construction are also outside the proof. No `assume`,
`admit`, or unchecked external body is added to the verified source.

The proofs establish functional properties conditional on this boundary;
they do not establish statistical uniformity, independence, cryptographic
unpredictability, or termination of `rand`'s rejection sampler. Loop termination
assumes each sample call returns, and allocation failure is outside the
contracts, as with the other allocating kernels. The `CryptoRng` bound and
`rand`'s existing sampling algorithm are retained. Seeded compatibility tests
check identical values and RNG consumption against the previous iterator
implementation, including rejection and empty outputs.

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
The sampling kernels additionally depend on the unverified sampler adapter
and its contract described above.

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

Multiplication tests exhaust small moduli, cover full-width Shoup operands
through `u64::MAX`, and check powers of two, even/composite moduli, and both
sides of each optimization threshold. Independent `u128` division and
remainder oracles validate Shoup words, lazy congruence, and canonical
products. All eight scalar kernels are exercised in const expressions.
Vector tests cover every variant, empty inputs, loop boundaries, and guarded
subslices; debug tests reject invalid operands, mismatched lengths, unsupported
optimized moduli, and incorrect Shoup words. Randomized Shoup tests retain
unrestricted first operands and use the same independent oracles:

```sh
cargo test -p fhe-math multiplication
```

Centering tests exhaust all canonical operands for small odd and even moduli
and exercise both sides of the midpoint at every supported bit width, through
`2^62 - 1`. An independent `i128::rem_euclid` oracle checks the result and its
congruence. Tests also cover const evaluation, empty and singleton vectors,
longer vectors, guarded input subslices, and debug rejection of noncanonical
entries. Randomized scalar and vector tests use the same independent oracle:

```sh
cargo test -p fhe-math center
```

Exponentiation tests exhaust canonical bases and permitted exponents for
moduli 2 through 64 using independent `u128` repeated multiplication. Boundary
tests and randomized tests compare against `BigUint::modpow`, including
full-width exponents, zero and one, all leading-bit positions, values adjacent
to powers of two, alternating bits, and the largest permitted exponent.
They cover odd, even, prime, and composite moduli and check debug rejection
of both invalid bases and invalid exponents:

```sh
cargo test -p fhe-math pow
```

Sampling tests compare both APIs with the original `Uniform` iterator using
seeded ChaCha generators. They cover zero and singleton lengths, block
boundaries, guarded subslices, small moduli, power-of-two boundaries, and the
largest modulus. Exact generator positions check RNG consumption, including
cases that exercise rejection. Mixed sequences of allocating and in-place
calls must match one uninterrupted call. Property tests vary moduli, seeds,
lengths, and the slice's initial contents; no test relies on two random
vectors happening to differ:

```sh
cargo test -p fhe-math zq::tests::random
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
Multiplication checks must reject dropping an operand, returning an incorrect
Shoup word, corrupting the lazy Shoup quotient, or skipping a vector product.
