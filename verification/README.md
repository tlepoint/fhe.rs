# Verifying arithmetic and packing with Verus

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

The constructor and kernel contracts share `MODULUS_LIMIT = 1u64 << 62`,
the exclusive upper bound for supported moduli. It leaves two high bits
available for lazy arithmetic; changing it requires revisiting those proofs.

[`verified/modulus.rs`](../crates/fhe-math/src/zq/verified/modulus.rs) proves the
range check and arithmetic initialization used by `Modulus::new`, establishing
the metadata predicates required by the reduction kernels.

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

[`ntt/native/butterflies.rs`](../crates/fhe-math/src/ntt/native/butterflies.rs)
proves all four native NTT butterfly methods, composing one-subtraction
reduction and lazy Shoup multiplication to establish modular results and
the bounds on lazy outputs.

[`verified/center.rs`](../crates/fhe-math/src/zq/verified/center.rs) proves scalar
and vector centering, including the signed range, congruence, and midpoint
behavior for odd and even moduli.

[`verified/pow.rs`](../crates/fhe-math/src/zq/verified/pow.rs) proves modular
exponentiation by composing multiplication contracts while tracking the
processed binary prefix of the exponent.

[`verified/inv.rs`](../crates/fhe-math/src/zq/verified/inv.rs) proves inversion
using that exponentiation. It includes a proof of Fermat's little theorem
from mathematical primality; the external primality detector remains outside
the verified module.

[`verified/random.rs`](../crates/fhe-math/src/zq/verified/random.rs) proves the
`random_vec` and `random_into` loops against an abstract sampler contract,
including canonical bounds, length, and exact correspondence to the sequence
of samples drawn. The adapter to `rand::Uniform` is a trusted boundary.

[`fhe-util/src/packing.rs`](../crates/fhe-util/src/packing.rs) verifies all five
public bit-transcoding functions, their validation, and their shared buffer
loop. Its [proof helpers](../crates/fhe-util/src/packing/proofs.rs) also prove
exact packing/decoding round trips. These are the implementations re-exported
by `fhe-util`, including the checked public boundaries.

[`verified/serialize.rs`](../crates/fhe-math/src/zq/verified/serialize.rs) composes
those packing contracts across the crate boundary to prove modulus widths,
serialized lengths, and encoding/decoding round trips.

[`zq/primes.rs`](../crates/fhe-math/src/zq/primes.rs) verifies the public
`generate_prime` search, including input validation, candidate alignment,
termination, maximality, and the exact absence of eligible primes on `None`.
It reuses the existing runtime primality contract.

## Run the proof

Use the pinned [Verus 0.2026.10.04.426d8b0 release](https://github.com/verus-lang/verus/releases/tag/release/0.2026.10.04.426d8b0)
and its required Rust toolchain:

```sh
rustup toolchain install 1.98.1 --profile minimal
```

Download and unpack the release archive for your platform, following the
[Verus installation instructions](https://github.com/verus-lang/verus/blob/426d8b01e7ffb910a36c15e1f870ff985bb61eef/INSTALL.md).
Add its directory to `PATH` so both `verus` and `cargo-verus` are available.
Install `protoc` as for ordinary `fhe-math` builds. From the repository root:

```sh
verus --version
cargo verus verify --fwd-verus-args-to roots -p fhe-util -p fhe-math --locked -- --no-cheating
```

Both crates opt in through `[package.metadata.verus] verify = true`.
Cargo Verus checks `fhe-util`, exports its contracts as `.vir` metadata, and
imports them when checking `fhe-math`. Serialization uses the ordinary
`fhe_util::packing` dependency path in both proof and production builds.
There is no source include outside either crate's package. Source releases
include the annotations; downstream Verus builds regenerate their metadata
with the selected toolchain. Ordinary Cargo builds erase the ghost code.
See the [Cargo Verus guide](https://verus-lang.github.io/verus/guide/cargo_verus.html).

Before publishing, release the updated `fhe-util` first and update the exact
dependency version in `fhe-math`. Published consumers need the release containing
the public packing module and its specifications. To check the source-package
layout locally, extract both `cargo package` archives into a temporary workspace
and patch its `fhe-util` dependency to the extracted package, then run the same
verification command there. This does not require files from the original checkout.

The expected project results are `65 verified, 0 errors` for `fhe-util` and
`253 verified, 0 errors` for `fhe-math` on macOS or an AVX2 target (`252`
elsewhere). Counts include proof lemmas, local assertions, and loops, not
just API functions. The first run also checks the pinned `vstd` dependency.
On x86 Linux, check both reducer implementations:

```sh
cargo verus verify --fwd-verus-args-to roots -p fhe-util -p fhe-math --locked -- --no-cheating -C target-feature=-avx2
cargo verus verify --fwd-verus-args-to roots -p fhe-util -p fhe-math --locked -- --no-cheating -C target-feature=+avx2
```

Select **both** project crates as roots when forwarding `--no-cheating` to
roots. It rejects assumptions, admits, and unchecked external bodies in the
verified project code. Do not forward it to `vstd`: the standard-library
specifications intentionally contain axioms and external specifications.
Verus, those specifications, Z3, Rust, and the compiler remain trusted.
The workspace enables `vstd`'s `std` feature so Cargo builds include its
`Vec` and allocation specifications, as the standalone release does.

Cargo integration verifies the code inside `verus!`, not the whole library.
The `Modulus` implementation and its public wrappers are included. The external
`Runtime` implementation in `zq/backend.rs` supplies explicit contracts for
CPU dispatch, uniform sampling, cloning runtime state, and primality detection.
`--no-cheating` checks the generic implementation; it does not verify that
external implementation's compliance with its trait contracts.

The current Cargo pin is `vstd = 0.0.0-2026-10-11-0230`, while the verifier
release is from October 4. Both work together in the Cargo verification run.
Rerun verification and Cargo checks when changing either pin. The verifier's
compiler is separate from the library's Rust 1.91.1 MSRV.

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

### Modulus construction

`new_modulus(p)` has no preconditions. It returns `Some(parameters)` exactly
when `2 <= p < 2^62`, with the original `p` and a proved
`ModulusParameters::valid()` predicate. That predicate establishes
`valid_barrett(p, hi, lo)`, the exact leading-zero count and its bounds, and
`supports_opt == valid_opt(p, leading_zeros)`. Thus both the ordinary and
optimized reduction metadata are established by the shared constructor code.

The reciprocal helper computes `q = u128::MAX / p` and
`r = u128::MAX % p`, then returns `q + (r == p - 1)`. Incrementing the dividend
from `2^128 - 1` to `2^128` increments its quotient exactly when the remainder
was `p - 1`; this handles powers of two as well as other supported moduli.
The proof checks division and subtraction safety, shows the addition cannot
overflow, and establishes the exact result `floor(2^128 / p)` and both
Barrett inequalities. The existing limb-splitting lemma then proves the
stored high and low words reconstruct this reciprocal. The leading-zero
lemma and verified `supports_opt` contract establish the remaining fields.

The public `Modulus::new` delegates to this initializer, maps `None` to the
existing `InvalidModulus(p)` error, and copies the returned arithmetic fields
into `Modulus`. It no longer uses `BigUint` to construct the reciprocal.
The wrapper's error mapping and every field copy are now verified too.

### Persistent Modulus invariant and public wrappers

[`zq/modulus.rs`](../crates/fhe-math/src/zq/modulus.rs) attaches a Verus
`type_invariant` to the actual stored object. It combines the Barrett reciprocal
predicate, exact leading-zero count, optimization eligibility, and equality
between the modulus and its sampling distribution's bound. All fields are
private. The NTT reads the modulus through a verified const accessor, and
`Deref` only exposes an immutable reference.

`Modulus::new` returns `Ok` exactly for supported moduli, preserving the input
value and optimization predicate. It returns exactly `InvalidModulus(p)`
otherwise. Construction establishes the invariant; `Clone` preserves the
complete logical state. `PartialEq` compares modulus values and `Deref`
returns that value with its range guarantee.

All production scalar and vector wrappers, including random generation and
serialization, now have checked contracts. Each arithmetic wrapper opens the
object invariant and proves its kernel's metadata preconditions. Callers still
must supply valid operands, Shoup words, compatible vector lengths, and the
stated bounds for optimized operations. Existing debug checks are proved safe
under those preconditions. The serialization-length wrapper requires aligned,
representable lengths; its documented panic cases are covered by runtime tests.

The concrete public type remains `zq::Modulus`, now an alias for an internal
`Modulus<Runtime>` implementation. Its generic `Backend` parameter allows Verus
to check the real wrapper bodies against explicit dependency contracts. There
is no dynamic dispatch or alternative arithmetic implementation. The production
backend retains `pulp::Arch` dispatch, `rand::Uniform`, and `fhe_util::is_prime`.
The trusted contract requires dispatch to satisfy its callback's postcondition,
samples to lie below the distribution bound, cloning to preserve backend state,
and the primality result to match mathematical primality. These are external
dependency assumptions, not new proofs of those dependencies.

Mutable dispatch callbacks take and return the slice borrow so Verus can connect
the callback's result to the caller's final slice. Arithmetic and vector loops
remain in checked code. [`zq/clients.rs`](../crates/fhe-math/src/zq/clients.rs)
checks composition through the concrete public alias: construction, cloning,
equality, dereferencing, scalar and vector operations, and optimization selection.
These clients supply no Barrett metadata assumptions.

### Prime generation

`generate_prime(num_bits, modulo, upper_bound)` has no caller preconditions.
It returns the greatest prime `p` satisfying
`2^(num_bits - 1) <= p < upper_bound` and `p % modulo == 1`. Supported searches
have `10 <= num_bits <= 62`, `modulo >= 2`, and
`upper_bound <= 2^num_bits`. `None` means no eligible prime exists for a
supported search; unsupported parameters also return `None`. Zero and other
empty upper bounds return `None` without subtraction or division errors.
Oversized upper bounds now return `None` in both debug and release builds,
replacing the former debug-only assertion.

The public `valid_prime_search`, `eligible_prime`, and `prime_search_result`
specifications describe the complete result. Successful results have the exact
bit width, satisfy the strict upper bound and congruence, and are maximal.
The result contract rules out a false `None` when any eligible prime exists.

The initial candidate is `start - ((start - 1) % modulo)`, where
`start = upper_bound - 1`. The proof checks every subtraction and shows this
is the largest candidate below the upper bound in residue class one. This
replaces the previous one-by-one alignment loop, avoiding a potentially large
scan when `modulo` is large. Subsequent candidates still decrease by `modulo`,
preserving the primes selected by valid existing callers.

A proved residue-spacing lemma shows that distinct candidates in the same
class differ by at least `modulo`. The search invariant states that every
eligible prime is at most the current candidate. Each failed primality check
and subtraction preserves that invariant; returning a prime proves maximality,
and leaving the bit-width interval proves absence. The candidate strictly
decreases, proving loop termination provided each primality check returns.

Only primality detection uses the existing `Backend::is_prime` assumption.
Both the public wrapper and generic search body are verified; no new trusted
implementation or number-theoretic axiom is introduced. The concrete public
API also has a proof client that feeds its result to `Modulus::new` and
preserves the complete prime-search contract.

### Lazy reduction

All lazy reducers establish `0 <= result < 2*p` and
`result mod p = input mod p`. They intentionally permit a noncanonical result
in `[p, 2*p)`. Verus also checks intermediate arithmetic, shift counts, casts,
and slice accesses. The scalar kernels retain `const fn` support and the
original multiplication/shift algorithms.

Let `B = 2^64` and `mu = hi*B + lo`. The `valid_barrett` contract requires
`2 <= p < 2^62` and `mu*p <= B^2 < (mu+1)*p`, equivalently
`mu = floor(2^128/p)`. The shared `new_modulus` initializer proves these
bounds for the reciprocal limbs used by `Modulus::new`. They remain
preconditions on standalone reduction kernels, without additional runtime
checks in those kernels.

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
described below. The native NTT butterflies also compose them; other
multiplication callers, including polynomial operations, remain outside the proof.

### Native NTT butterflies

The four production methods in `ntt/native/butterflies.rs` have these contracts,
where `x` and `y` in the formulas denote the inputs:

| Methods | Input bounds | Output residues modulo `p` | Output bounds |
| --- | --- | --- | --- |
| `butterfly`, `butterfly_vt` | `x, y < 4p` | `(x + wy, x - wy)` | Both `< 4p` |
| `inv_butterfly`, `inv_butterfly_vt` | `x, y < 2p` | `(x + y, z(x - y))` | Both `< 2p` |

All four require `p_twice = 2p` and a canonical twiddle with its exact Shoup
word (`valid_shoup`). The verified `Modulus` accessor supplies `2 <= p < 2^62`.
Neither primality nor roots of unity are needed for these local identities.
The proofs also establish that executable additions and subtractions cannot
overflow or underflow, and that the existing debug assertions hold.

The forward proof uses `(x mod 2p) mod p = x mod p`, the lazy Shoup product's
congruence, and the fact that adding `2p` preserves the difference modulo `p`.
The inverse proof uses the same reduction identity for the sum and removes
the `2p` shift before multiplying the difference by the inverse twiddle.
Ordinary and variable-time methods share these algebra lemmas; Verus checks
each executable method separately. The shared Shoup quotient lemma also
justifies the runtime debug check against `Modulus::shoup`.

`NttOperator`'s actual private fields are visible to Verus, and normal builds
call these same methods. Its constructor, derived cloning, twiddle tables,
and transform loops are not yet verified. In particular, callers must still
establish the explicit `p_twice`, operand, and Shoup preconditions; there is
no assumed NTT type invariant. These local proofs do not establish the
correctness or constant-time machine code of a complete transform.

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
behavior. `inv` composes this exponentiation with the number-theoretic proof
described below.

### Inversion

`inv_mod(a, p, hi, lo, p_is_prime)` requires `valid_barrett`, a primality flag
equal to the mathematical predicate `prime(p)`, and `a < p` when the flag is
true. It returns `None` exactly when the flag is false or `a = 0`. Every
`Some(r)` satisfies `0 < r < p`, `(a*r) mod p = 1`, and
`r = a^(p-2) mod p`. This includes `p = 2` and the zero exponent. Composite
moduli still return `None` even for invertible operands. As before, invalid
operands are rejected in debug builds for prime moduli, while composite
moduli return `None` without inspecting the operand's canonicality.

The predicate `prime(p)` means `p > 1` and no integer `d` with `1 < d < p`
divides `p`. It does not assume Fermat's theorem or the existence of inverses.
The proof builds those facts in stages:

1. A ghost Euclidean recursion constructs a positive common divisor and
   Bezout coefficients. For `0 < a < p` and prime `p`, that divisor must be
   one, giving a modular inverse witness and cancellation of nonzero factors.
2. Recursive sequence lemmas prove that removing a factor, permuting distinct
   factors, and scaling all factors have the expected effects on their
   product modulo `p`. A product of nonzero residues remains nonzero.
3. Multiplication by `a` maps `1, ..., p-1` to distinct nonzero residues.
   Equal finite cardinalities make this a permutation. Equality of the
   products gives `a^(p-1) * product = product` modulo `p`; cancellation proves
   Fermat's little theorem.
4. The existing `pow_mod` contract and the power-addition identity give
   `a * a^(p-2) = 1` modulo `p`, establishing the executable kernel's result.

All these lemmas are checked with `--no-cheating` and erase in ordinary Rust
builds. The runtime algorithm remains the existing primality gate followed by
Fermat exponentiation; the Euclidean recursion and sequence products run only
in the proof.

**Primality boundary:** the verified wrapper obtains its flag from the
backend's primality contract and passes it to the kernel. The arithmetic
metadata comes from the checked object invariant. The inverse guarantee remains
conditional on `fhe_util::is_prime` and its `num-bigint-dig` dependency correctly
recognizing mathematical primes. Their implementation is not proved here.
Fermat's theorem itself is proved without an unchecked axiom.

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

**Trusted sampling boundary:** `Runtime` in `zq/backend.rs` stores the
constructor's `Uniform::new(0, p)` and calls `rng.sample(distribution)` once per
backend draw. Its range and distribution-construction contracts are assumed.
The `Sampler` bridge is verified: it appends each actual returned value to its
ghost history. The public wrappers now prove that this adapter uses the object's
modulus and that their outputs have the promised length and canonical bounds.

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

The arithmetic proofs assume valid operands. The checked constructor and
persistent invariant now establish the arithmetic metadata used by every
production wrapper. The wrappers retain their documented debug assertions;
Verus does not insert runtime validation or verify all callers elsewhere in
the library. It also does not establish constant-time machine code, NTT/RNS
correctness, or FHE security. Test-only scalar convenience wrappers are outside
the normal library proof, while their shared scalar kernels are verified.

Runtime dependency contracts remain the explicit boundary: CPU dispatch must
execute the callback according to its contract, uniform sampling must honor
its range and configured modulus, and primality detection must match its
mathematical predicate. Statistical uniformity and RNG security are not proved.

### Checked bit packing and unpacking

`fhe-util` re-exports `transcode_to_bytes`, `transcode_to_bytes_into`,
`transcode_from_bytes`, `transcode_from_bytes_exact`, and
`transcode_bidirectional` from the verified packing module. Their contracts
have no caller preconditions: they specify successful output and the exact
error for invalid input, including validation order. The `TranscodeError`
type and its variants are unchanged.

The specification defines an LSB-first stream by bit `k` of input word
`k / input_width`, at position `k % input_width`, and extends that stream
with zeros past its end. Each output word contains exactly its corresponding
slice of this stream, with zero upper bits. The output length is
`ceil(input_length * input_width / output_width)`. This specifies individual
bits independently of the buffer algorithm, including widths 1 and 64,
empty inputs, partial bytes, and partial output words.

One generic loop replaces the three copies of the `u128` buffer algorithm.
Its `u8` and `u64` storage adapters are verified in the same module. Before
loading a word, fewer than `output_width` bits remain in the buffer; loading
at most 64 more bits leaves at most 127 bits, so the buffer cannot truncate
data. The invariant relates emitted words and every buffered bit to the
input stream. Bounds, shifts, arithmetic, and both loops' termination are
checked. A verified ceiling-division helper computes quotient plus a nonzero
remainder flag without adding to the potentially maximal input bit length.

Validation still scans every input word with an aggregate OR, rejecting any
word that exceeds its declared width. Invalid widths and unrepresentable bit
lengths retain their earlier precedence. Append checks the total output
length before reserving or writing, proves that errors preserve the output
sequence, and proves successful calls preserve the old prefix and begin at
a byte boundary. The executable code does not touch allocation on an error;
pointer/capacity preservation is covered by Rust tests, since the standard
Verus vector model exposes contents rather than allocation identity.

Stream decoding consumes every byte bit as data and zero-extends the final
word. Exact decoding additionally checks the byte count and rejects nonzero
padding before truncating to the declared word count. Its contract retains
the existing overflow error if the rounded byte bit length cannot fit in
`usize`. `packing_is_decodable` proves the size and padding checks accept
packed data when that rounded length is representable; `packing_round_trip`
proves that successful exact decoding then returns the original sequence.
Word transcoding similarly preserves the zero-extended stream but does not
retain the original bit count.

There is no trusted packing adapter or unchecked project proof body. As with
the arithmetic kernels, Verus, Z3, the compiler, and standard-library specs
remain trusted. Allocation failure, allocation-capacity limits, and constant-time
machine-code behavior are outside these contracts. Packing does not validate
residues against a modulus. The modulus-specific composition is described below.

### Modulus serialization

[`verified/serialize.rs`](../crates/fhe-math/src/zq/verified/serialize.rs)
composes `fhe_util::packing` contracts across the crate boundary. Its shared
functions implement `Modulus::serialize_vec`, `serialize_into`, and
`deserialize_vec`, plus the checked computation used by `serialization_length`.
They assume a valid modulus `2 <= p < 2^62`; they need no Barrett metadata,
optimization flag, or primality assumption.

The width is proved to be `ceil(log2(p))`, between 1 and 62, with
`2^(width-1) < p <= 2^width`. Every canonical residue therefore fits, including
at powers of two where the width is determined by `p-1`. Packing and stream
decoding inherit the exact bit layout, lengths, padding, and errors from
`fhe-util`. Fitting noncanonical words remain accepted, and decoding does not
reduce or validate residues. The append kernel proves its existing unwrap
safe for fitting words and representable input/output lengths, and preserves
the existing byte prefix.

The round-trip theorem proves decoding recovers the original prefix for
arbitrary word counts. If the encoded bit count is a multiple of eight,
there are no extra words and the full sequences are equal. A word count
that is a multiple of eight always satisfies this alignment condition.
This preserves stream decoding's treatment of byte padding as data.

The checked length computation returns `Some(bytes)` exactly when the word
count is a multiple of eight and `count * width / 8` fits in `usize`. It now
computes `width.checked_mul(count / 8)`: the previous multiply-before-divide
expression could overflow even when the final byte count was representable,
and could wrap in release builds. The public `const fn` retains its alignment
assertion and explicitly panics on an unrepresentable byte count. The proof
also establishes agreement between this byte count and the packing output.

The serialization wrappers now obtain their modulus bounds from the persistent
object invariant, and their field access and delegation are verified. The
length wrapper's panic mapping is tested in debug and release builds. These
proofs do not yet cover polynomial/protobuf serialization, RNS row concatenation, representation
conversion, or decode resource limits.

## Tests and CI

Prime-search tests compare every upper bound for 10- and 11-bit searches
against an independent trial-division oracle across small moduli and selected
large steps. They also cover invalid parameters, empty ranges, strict bounds,
large-step alignment and exhaustion, and the existing NFL prime sequence.
Invalid-input and boundary cases run in release builds as well.
Disposable mutation checks confirmed that Verus rejects an unjustified `None`,
a wrong returned prime, and an incorrect candidate step.


The public-wrapper regression clients are checked by the same Cargo Verus
command. Runtime tests also compare cloned metadata, equality and dereferencing,
dispatched reduction, and identical seeded RNG streams across boundary moduli.
Disposable mutation checks confirmed rejection of an incorrect constructor
reciprocal limb, a flipped clone optimization flag, a duplicated scalar operand,
and substitution of subtraction for addition inside a vector callback.


Modulus serialization tests cover every supported width, power-of-two
boundaries, canonical and fitting noncanonical words, arbitrary byte streams,
partial bytes, append prefixes and error atomicity. Independent bit-at-a-time
oracles check encoding and decoding. Length tests reach the largest
representable byte counts without allocating buffers, exercise const evaluation,
and require overflow rejection in both debug and release builds:

```sh
cargo test -p fhe-math zq::tests::serial
cargo test --release -p fhe-math zq::tests::serial
```

Packing tests compare all 4,096 width pairs against an independent bit-at-a-time
oracle, including empty inputs, word boundaries, full-width words, and partial
outputs. Property tests vary words, lengths, and widths. Additional tests cover
exact decoding, every padding-bit position, append boundaries, error precedence,
unchanged output allocations on validation errors, and ceiling lengths through
`usize::MAX`:

```sh
cargo test -p fhe-util
```

Constructor tests exhaust all valid 16-bit moduli and cover every supported
power of two, adjacent values, optimization thresholds, and the maximum
62-bit modulus. They compare the stored reciprocal against the original
`BigUint` division and the optimization flag against the original NFL
inequality. Property tests check exact acceptance/error behavior over the
full `u64` domain and compare metadata for random valid moduli:

```sh
cargo test -p fhe-math zq::tests::constructor
```

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

Direct butterfly tests exhaust all lazy inputs and canonical twiddles for
moduli 2 through 13. Boundary and seeded randomized cases cover even/composite
moduli and the largest supported 62-bit values. An independent signed `i128`
remainder oracle checks both output formulas, lazy bounds, and equality of
ordinary and variable-time outputs. Debug tests reject out-of-range operands,
noncanonical twiddles, and incorrect Shoup words for all four methods. Existing
native NTT round-trip and reference ring-product tests exercise their callers.

```sh
cargo test -p fhe-math butterflies
```

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

Inversion tests exhaust all canonical operands for moduli 2 through 256,
using trial division and exhaustive inverse search as independent oracles.
Boundary and property tests compare with `fhe_util::inverse`, which uses
`num-bigint-dig`'s extended Euclidean inversion rather than Fermat
exponentiation. They include `p = 2`, zero, one, `p-1`, power-of-two operand
boundaries, and large 62-bit primes. Composite tests include Carmichael
numbers, invertible operands such as one, and noncanonical inputs; separate
debug tests check rejection of invalid operands for primes:

```sh
cargo test -p fhe-math zq::tests::inv
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
through Cargo with both project crates selected as roots and `--no-cheating`.
It disables and enables AVX2 to cover both `reduce1_vt` implementations.
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
