# Proposal: simplify the public API and strengthen Rust contracts

Reviewed on 2026-09-07 against commit `f2ac193` and the current working tree.
Scope: the four workspace crates, their public interfaces, BFV examples,
selected implementations and tests, and experimental multiparty BFV.
This document proposes a breaking release; it does not implement the redesign.
Examples under “Proposed” describe the target design; see the migration notes
below for the implemented APIs.

Implementation update: the P0 contract work is now implemented. See the
[P0 migration notes](api-p0-migration.md) for the actual changes and current
signatures. P1 is also implemented; see the [P1 migration notes](api-p1-migration.md).
The analysis below records the original findings; P2 remains a proposal.

## Recommendation

Keep the arithmetic implementation and simplify the interface around it.
The current tree is already substantially more disciplined than a typical
first Rust project: it has typed polynomial representations, structured errors,
validated parameter construction, explicit variable-time permissions, and
useful allocation-conscious APIs. A wholesale rewrite would discard that work.

The largest improvement would be to make ordinary BFV use require only `fhe`,
concrete types, and inherent methods. At the same time, fix public contracts
that currently permit invalid states or surprising behavior. Those fixes matter
more than changing method spellings.

| Priority | Change | Reason | Relative effort |
| --- | --- | --- | --- |
| P0 | Fix plaintext equality; remove invalid default/sentinel states | Existing trait and value contracts are inconsistent | Small–medium |
| P0 | Close mutation paths around validated objects | Safe callers can invalidate cached or contextual invariants | Medium |
| P0 | Redact secret-bearing `Debug` output; complete timing API cleanup | Ordinary language features should not reveal coefficients or obscure timing permission | Small–medium |
| P1 | Make common operations inherent; retire speculative FHE traits | Removes imports, output-type ambiguity, and conversion scaffolding | Medium |
| P1 | Hide shared parameter ownership and define compatibility | Removes `Arc` from routine use and resolves identity surprises | Medium–large |
| P1 | Provide a consistent checked arithmetic interface | Normal runtime mismatches should be recoverable | Medium |
| P1 | Separate encoding, level, and decoded integer representation | Removes optional-encoding ambiguity and repeated configuration | Medium |
| P1 | Simplify builders and serialization entry points | Fewer ways to construct the same object; clearer validation boundaries | Medium |
| P2 | Organize advanced operations and simplify collection APIs | Retains performance without making it the entry-level interface | Medium |
| P2 | Tighten the math and multiparty extension surfaces | Make deliberate extension points and protocol states explicit | Medium |

“P0” denotes contract fixes to do first, not a security severity rating. The
effort estimates describe breadth of changes, not measured implementation time.

## 1. Fix value contracts before reorganizing the API

### Plaintext equality is not transitive

In [`Plaintext::eq`](../crates/fhe/src/bfv/plaintext.rs), encodings are compared
only when both are `Some`. Decryption constructs a plaintext with `encoding:
None`. For the same zero polynomial:

```text
A: encoding = Some(Poly)
B: encoding = None
C: encoding = Some(Simd)

A == B, B == C, but A != C
```

This was reproduced through public APIs. It violates the transitivity required
by [`PartialEq`](https://doc.rust-lang.org/std/cmp/trait.PartialEq.html), even
without the additional `Eq` implementation. Removing only `Eq` would not fix it.

**Proposed:** remove encoding metadata from `Plaintext` as described in section
5. Define equality as equality of the represented polynomial, level, and
compatible parameters. Until that change lands, compare the entire
`Option<Encoding>` exactly; use a separately named comparison if tests need to
ignore metadata. Do not use unknown metadata as a wildcard in `==`.

Also specify that ciphertext equality is structural and does not mean “encrypts
the same message.” Audit whether local serialization seed/cache state should
participate in that structural equality. Parameter equality should compare the
defining values, not precomputed tables or memory addresses.

### The zero ciphertext is an accumulator sentinel

[`Ciphertext::new`](../crates/fhe/src/bfv/ciphertext.rs) requires at least two
polynomials, but `Ciphertext::zero` creates an empty vector. Operators have
special branches for this value. Decryption, checked multiplication, and
deserialization reject it. In particular:

```rust,ignore
let zero = Ciphertext::zero(&params);
assert!(Ciphertext::from_bytes(&zero.to_bytes(), &params).is_err());
```

The serialization round-trip failure and decryption rejection were reproduced.
This is a representation of “no accumulated value,” exposed as a ciphertext.

**Proposed:** remove the empty ciphertext state. Use `Option<Ciphertext>` inside
an additive accumulator or initialize a reduction from its first input. Keep
the existing [`CiphertextProductAccumulator`](../crates/fhe/src/bfv/ops/product_accumulator.rs)
for its distinct fused product-sum semantics. Define empty reductions explicitly:
return `None`, return `EmptyInput`, or require parameters and a level to construct
a valid trivial zero. Do not invent a parameterless `Default`/`Sum` identity.

If a public zero constructor remains, call it `trivial_zero(parameters, level)`
and construct a valid pair of zero polynomials. Document that this is a public,
deterministic value; encrypting zero with an RNG is a separate operation. Every
ordinary ciphertext should then serialize and deserialize successfully.

### Some `Default` implementations bypass validation

[`rq::Context`](../crates/fhe-math/src/rq/context.rs) and
[`RnsContext`](../crates/fhe-math/src/rns/mod.rs) derive `Default`, although their
constructors reject zero degree or empty moduli. Likewise,
[`ScalingFactor`](../crates/fhe-math/src/rns/scaler.rs) derives a default with
zero denominator, which its constructor rejects. `Scaler`, `RnsScaler`, and
internal multiplication parameters also use derived defaults.

**Proposed:** remove defaults that produce uninitialized mathematical objects.
Use `Option<T>` for internal construction slots. If `ScalingFactor` needs a
default, implement a valid identity explicitly. A builder may implement
`Default` while producing a validation error from `build()`; that is a different
contract from default-constructing an invalid finished object.

## 2. Make validated objects stay valid

[`Ciphertext`](../crates/fhe/src/bfv/ciphertext.rs) dereferences to
`[Poly<Ntt>]` and implements `DerefMut`. A caller can replace one component with
a polynomial from another degree or modulus chain. This was reproduced; a
subsequent checked operation rejects the mutated ciphertext. The implementation
already invalidates the compressed seed on mutable access, which is good, but
that does not preserve the rest of the ciphertext invariants.

There are smaller versions of the same problem:

- [`SubstitutionExponent::exponent`](../crates/fhe-math/src/rq/mod.rs) is public
  and mutable although construction precomputes a permutation from it.
- [`ContextLevel::poly_context`](../crates/fhe/src/bfv/context/level.rs) is public.
  A caller can modify a cloned `ContextLevel` independently of its other cached
  data. The shared parameter table itself is only exposed immutably.
- The public `CipherPlainContext` export has no public constructor or useful
  public methods; it exposes an implementation concept without a user operation.

**Proposed:**

- Remove `Ciphertext`'s `Deref` and `DerefMut`. Provide `level()`,
  `component_count()`, and `parameters()` as inherent accessors. `Ciphertext`
  currently has no public `level()` accessor despite tracking the value.
- Provide an explicitly advanced, immutable `components()` view only where
  downstream algorithms need it. Keep mutation crate-private.
- Replace the generic `new(Vec<Poly<Ntt>>, ...)` with a checked
  `from_components(...)` import path. If editable components are necessary,
  consume the ciphertext into owned parts and reconstruct through that path;
  do not return unrestricted mutable references or validate only in `Drop`.
- Make cached descriptors' fields private and provide getters. Make
  `CipherPlainContext` crate-private and either hide `ContextLevel` or expose
  an immutable level view with a clear purpose.

`Deref` changes method lookup and implicitly exposes the target interface;
the standard library explicitly calls out that API commitment in its
[`Deref` documentation](https://doc.rust-lang.org/std/ops/trait.Deref.html).
A ciphertext should have a ciphertext interface. A transparent collection
wrapper can still reasonably offer slice access; it needs a separate assessment.

## 3. Make the common path concrete and inherent

[`bfv_basic.rs`](../crates/fhe/examples/bfv_basic.rs) imports `FheEncoder`,
`FheDecoder`, `FheEncrypter`, and `FheDecrypter` just to encode, encrypt, decrypt,
and decode. Other examples use a wildcard import from `fhe_traits`.
The root README requires both crates, while the `fhe` README's dependency
snippet lists only `fhe` even though its example imports `fhe_traits`.

[`fhe-traits`](../crates/fhe-traits/src/lib.rs) contains several abstractions
with little current payoff:

- `FheParameters` and `FhePlaintextEncoding` are empty markers with one
  implementation each in this workspace.
- `FheParametrized` supplies an associated type but no parameter accessor.
- `FheParametersSwitchable` has no implementation in the workspace.
- `FheCiphertext` requires serialization and parameterized deserialization,
  coupling arithmetic objects to persistence.
- `FheDecoder` puts decoding on `Vec`, so discoverability starts at the output
  container rather than the plaintext.
- Secret-key encryption has both BFV and RGSW implementations; the trait's
  result parameter makes a routine call depend on an output annotation.

**Proposed common interface:**

```rust,ignore
use fhe::bfv::{Encoding, Parameters, Plaintext, PublicKey, SecretKey};

fn main() -> fhe::Result<()> {
    // Explicit arithmetic configuration, not a new security-selection policy.
    let params = Parameters::builder()
        .degree(2048)
        .plaintext_modulus(1024_u64)
        .ciphertext_moduli([0x3fffffff000001])
        .build()?;
    let mut rng = rand::rng();
    let sk = SecretKey::generate(&params, &mut rng);
    let pk = PublicKey::from_secret_key(&sk, &mut rng);

    let a = Plaintext::encode(&params, &[20_u64], Encoding::Polynomial)?;
    let b = Plaintext::encode_signed(&params, &[-7_i64], Encoding::Polynomial)?;
    let ca = sk.encrypt(&a, &mut rng)?;
    let cb = pk.encrypt(&b, &mut rng)?;
    let product = ca.multiply(&cb)?; // Unrelinearized; level is unchanged.
    let decoded = sk.decrypt(&product)?.decode_signed(Encoding::Polynomial)?;
    assert_eq!(decoded[0], -140);
    Ok(())
}
```

Make these real inherent methods, not merely a prelude that hides the existing
imports. Use `encrypt_rgsw` or a constructor on `RgswCiphertext` for the alternate
output. Retain explicit caller-owned RNGs and their cryptographic bounds.

Use slices for numeric inputs: `&[u64]`, `&[i64]`, and `&[BigUint]`. Ordinary
argument coercions already accept borrowed vectors and arrays. This eliminates
the forwarding implementations for `&Vec<T>` and `&[T; N]` in
[`plaintext.rs`](../crates/fhe/src/bfv/plaintext.rs). Prefer a small set of named
numeric methods over a new public “all plaintext integers” trait until real
generic callers justify one. Decode on `Plaintext`, with return types stated by
`decode`, `decode_signed`, and `decode_biguint`.

Retire the unused/generic FHE hierarchy in the breaking release unless an
identified downstream abstraction needs it. A second future scheme can first
have its own concrete API. If a generic encryptor eventually becomes useful,
derive a small capability trait from actual consumers and keep serialization
independent. The choice to use inherent methods follows the discoverability
guidance in the [Rust API Guidelines](https://rust-lang.github.io/api-guidelines/predictability.html#functions-with-a-clear-receiver-are-methods-c-method).

**Crate dependency constraint:** `fhe-math` already depends on `fhe-traits` for
timing tokens. Moving those tokens directly into `fhe` would introduce a cycle.
Move them into a small shared support module in `fhe-util`, or keep a deliberately
small shared crate, and re-export them from `fhe` and `fhe-math`. Rename/remove
`fhe-traits` only after that dependency has been dealt with. Keep `fhe-math` as
the independently usable low-level crate; do not collapse all packages simply
to make the README shorter.

## 4. Hide `Arc` and choose a parameter compatibility contract

Nearly every constructor accepts `&Arc<BfvParameters>`. The builder offers both
`build` and `build_arc`, but ordinary keys and plaintexts require the latter.
This makes sharing strategy part of the user API without offering a meaningful
choice in the common path.

There is also a semantic inconsistency. BFV validation and operators require
`Arc::ptr_eq`, whereas some polynomial/context operations accept value-equal
contexts. Two independently built `BfvParameters` compare equal but a key from
one rejects a plaintext from the other. This was reproduced.

**Proposed:** use a cheap, immutable owning handle:

```rust,ignore
#[derive(Clone)]
pub struct Parameters(Arc<ParametersInner>); // The field and inner type are private.

impl ParametersBuilder {
    pub fn build(self) -> Result<Parameters>;
}
```

Public constructors borrow `&Parameters`; retained objects clone the handle.
Document that `Parameters::clone` shares precomputed data. Do not replace all
borrows with `&T` mechanically: construction really does retain shared ownership,
and the wrapper is what makes the signature honest and convenient.

Recommend that equal defining parameter values be compatible, including when
built or deserialized separately. Use shared identity as a fast path, followed
by comparison of a canonical parameter specification. Include the exact ordered
moduli, degree, plaintext modulus, noise configuration, and any other
compatibility-relevant settings. Exclude NTT caches and scratch buffers.
Do this once per operation/import boundary, not once per coefficient. Ensure
the low-level checks accept the same relation, and benchmark the fallback.

A fingerprint can accelerate rejection or identify wire parameters, but must
not silently replace the actual compatibility contract. Parameter compatibility
also does not establish that ciphertexts use the same secret key.

The lower-effort alternative is explicit session identity: separately created
contexts remain incompatible and all imports bind into one supplied context.
That can be coherent, but it must be named, documented, and applied consistently.
Do not preserve today's mixture accidentally. Resolve this decision before the
serialization migration and before proliferating additional context types.

Do not add a monolithic `Bfv` object owning secrets, evaluation keys, RNG, and
scratch. Parameters, keys, and per-call mutable workspaces have different
ownership and sharing requirements.

## 5. Separate encoding from level and numeric output

[`Encoding`](../crates/fhe/src/bfv/encoding.rs) wraps a private enum and a level,
with `poly`, `simd`, `poly_at_level`, and `simd_at_level` constructors. Plaintexts
also carry their actual level in their polynomial context. Encoding metadata
is optional and disappears on decryption. Decoding accepts
`Into<Option<Encoding>>`, introducing both inference and override behavior.

**Proposed:**

```rust,ignore
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Encoding {
    Polynomial,
    Simd,
}

impl Plaintext {
    pub fn encode(params: &Parameters, values: &[u64], encoding: Encoding)
        -> Result<Self>; // Level 0.
    pub fn encode_at_level(
        params: &Parameters, values: &[u64], encoding: Encoding, level: usize,
    ) -> Result<Self>;
    pub fn decode(&self, encoding: Encoding) -> Result<Vec<u64>>;
    pub fn level(&self) -> usize;
}
```

Keep one authoritative level in the value. Start with checked `usize` level
arguments; introduce a small validated level descriptor only if it removes
repeated validation or demonstrable confusion. A bare newtype alone would not
prove compatibility with a particular parameter set.

Remove stored encoding from `Plaintext` and require an explicit interpretation
when decoding. A decrypted polynomial cannot reliably reveal how the sender
encoded it. This design deliberately lets applications interpret the same
polynomial as coefficients or slots; their protocol supplies the intended
meaning. If preserving that intent is essential, use an application wrapper or
an explicit encoded-message wrapper, rather than optional wildcard metadata.

Document these semantics alongside the methods:

- Polynomial encoding gives negacyclic polynomial arithmetic; SIMD gives slot
  arithmetic and has parameter restrictions. For the current NTT path the
  modulus must be prime and congruent to 1 modulo `2 * degree`; the shorter
  `Encoding::simd` rustdoc currently omits the factor of two. The actual check is
  in [`supports_ntt`](../crates/fhe-math/src/ntt/mod.rs).
- Level zero uses the full chain; increasing the level drops moduli. Neither
  decoding nor a different encoding argument changes a plaintext's level.
- Define reduction of unsigned, signed, and `BigUint` inputs, zero padding,
  maximum length, centered signed output, and overflow errors explicitly.
  Preserve numerical behavior during the API refactor and test boundaries.
- Decoding returns `degree` values; the original input length is not recoverable
  from the ciphertext. Any truncation is an explicit application decision.
- Big-integer support does not imply that the existing SIMD implementation
  supports arbitrary large plaintext moduli. Preserve capability checks.

Avoid making every ciphertext generic over encoding or level. Those parameters
would spread through keys, serialization, and heterogeneous collections without
proving the caller's message interpretation after a network round trip.

## 6. Make checked arithmetic the ordinary interface

[`ops/mod.rs`](../crates/fhe/src/bfv/ops/mod.rs) implements `+`, `-`, `*`, and
assignment operators using assertions for parameter, level, or component-count
mismatches. Multiplication has acquired `try_mul` and `square`, and
[`Multiplicator`](../crates/fhe/src/bfv/ops/mul.rs) has checked methods, but
addition/subtraction/plaintext operations do not have an equivalent complete
interface. Examples still lead with panicking operators.

**Proposed:** provide one checked kernel per operation, with inherent methods
such as `add`, `add_assign`, `subtract`, `subtract_assign`, `multiply`, and
`multiply_plaintext`, returning `Result`. Use longer names or `_plaintext`
suffixes to disambiguate operand kinds; do not introduce another generic
operator framework. Keep allocating and in-place methods distinct. Validate
all recoverable preconditions before changing the receiver so a returned error
leaves it unchanged. Existing switching/accumulator contracts are useful models.

Recommend removing binary arithmetic and assignment operator implementations
from the main ciphertext type in the breaking release. Runtime levels and
parameters make failure ordinary, and the required checked calls should be
visible in examples. Keeping documented panicking operators as thin wrappers
is a viable compromise if downstream usage strongly favors notation; it should
be an explicit choice after implementing the checked methods. Do not make `+`
return `Result<Ciphertext>` merely to preserve the glyph.

Use `try_` when distinguishing a fallible alternative from a retained infallible
method, or where standard conversion traits require it. A lone fallible method
can simply be `encrypt`, `decode`, or `multiply`; the return type states failure.
Do not turn invariant-protected key generation into a spurious `Result` just to
make every function look alike. Verify the internal invariants before retaining
infallible constructors.

Specify multiplication precisely: unrelinearized multiplication of `m` and `n`
components yields `m + n - 1` components at the same level; relinearization and
modulus switching are explicit additional operations. The current general
`Ciphertext::try_mul` accepts more than two components, whereas the configured
`Multiplicator` validates exactly two. Preserve that distinction or deliberately
generalize the plan; a rename must not silently narrow the existing capability.

Retain `square()` and its dedicated algorithm. Retain
`CiphertextProductAccumulator`'s bound on accumulated products and its single
final rounding. It is not interchangeable with repeatedly calling `multiply`
and adding separately rounded results.

## 7. Simplify construction and names

The parameter builder is already validating a substantial set of constraints;
retain that validation. The improvements below concern representing choices
and giving callers a predictable construction path.

| Current API | Proposed API / action |
| --- | --- |
| `BfvParameters`, `BfvParametersBuilder` inside `bfv` | `Parameters`, `ParametersBuilder`; qualify with `bfv::` where needed |
| `BfvParametersBuilder::new()` | Primary entry point `Parameters::builder()` |
| `set_degree`, `set_variance` | `degree`, `noise_variance` on a consuming builder |
| `set_plaintext_modulus` / `_biguint` | One `plaintext_modulus(impl Into<BigUint>)`; examples use typed integer literals |
| `set_moduli` / `set_moduli_sizes` | `ciphertext_moduli` / `ciphertext_modulus_bits`; store one internal `ModuliSpec` enum |
| `build()` / `build_arc()` | One `build() -> Result<Parameters>` returning the shared handle |
| `plaintext()` (panics for large moduli) / `plaintext_big()` | `plaintext_modulus() -> &BigUint`, `plaintext_modulus_u64() -> Option<u64>` |
| `default_parameters_128(bits)?.nth(index)` | Select a documented named profile or an explicit degree; build only the requested profile |
| `SecretKey::random` | `SecretKey::generate` (a naming preference, not a correctness fix) |
| `PublicKey::new` | `PublicKey::from_secret_key` |
| `relinearizes`, `rotates_rows`, `rotates_columns_by`, `computes_inner_sum`, `expands` | `relinearize`, `rotate_rows`, `rotate_columns`, `inner_sum`, `expand` |
| `RGSWCiphertext` | `RgswCiphertext` |
| `Multiplicator::default(rk)` | `MultiplicationPlan::with_relinearization(rk)`; avoid a non-`Default` method named `default` |

Consuming builders are a project choice, not the only idiomatic Rust pattern.
Their benefit here is a clear transition into a finished value, fewer copies
of configured vectors, and consistent usage across parameter and key builders.
Represent unspecified required fields as `Option`, rather than using zero as
an implicit missing-value marker. For mutually exclusive modulus settings,
document that the most recent setter replaces the prior selection, or expose a
single explicit enum setter. Eliminate the invalid “both specified” state.

[`EvaluationKeyBuilder::new`](../crates/fhe/src/bfv/keys/evaluation_key.rs)
returns `Result` but unconditionally returns `Ok`, and it clones a full secret
key into the builder. Prefer a builder borrowing `&SecretKey`, consumed by
`build(rng)`. Validate configuration in setters only where needed for subsequent
steps; otherwise validate in `build`. Preserve distinct ciphertext and
evaluation-key levels, using named configuration rather than adjacent integer
arguments in `new_leveled` constructors.

The existing default-parameter iterator sorts candidates, so its order is not
random. The problem is that it eagerly builds all surviving parameter sets and
`.nth(2)` couples selection to filtering. Separate lightweight profile listing
from expensive construction. Preserve and document the existing security
assumptions; ordinary builder validation is not a security estimator.

## 8. Keep serialization explicit and context-aware

[`fhe-traits`](../crates/fhe-traits/src/lib.rs) has `Serialize`, `Deserialize`,
`DeserializeParametrized`, and `DeserializeWithContext`, with inconsistent
`from_bytes` / `try_deserialize` naming. There are two unrelated
`TryConvertFrom` traits in
[`bfv::traits` (now private wire helpers)](../crates/fhe/src/bfv/wire.rs) and
[`rq::traits`](../crates/fhe-math/src/rq/traits.rs). The `fhe` crate also publicly
exports generated protobuf types through [`proto`](../crates/fhe/src/proto/mod.rs).

**Proposed:** expose inherent `to_bytes` and `from_bytes` methods for supported
objects. Parameterized values take `&Parameters` or the math context explicitly.
Use private wire DTO conversion functions and make `fhe::proto` private, or put
it in a deliberately versioned interoperability module if consumers need raw
protobuf access. Naming a custom trait `Serialize` does not make it Serde's trait.

Use standard `From`/`TryFrom` for actual unambiguous conversions whose input
contains everything required. Additional-context conversions can be inherent
constructors; a tuple or a private conversion-input type can support `TryFrom`
internally where generic composition is useful. The extra arguments do not
require a public replacement for `TryFrom`.
[`TryFrom`](https://doc.rust-lang.org/std/convert/trait.TryFrom.html) is the
standard fallible-conversion contract; it need not be forced onto every decoder.

Do not blindly derive Serde on cryptographic structs: their fields include
contexts, caches, seeds, and local execution policy. A Serde adapter, if needed,
should use a validated wire representation and an explicit context-binding
mechanism. API compatibility and stored-byte compatibility are separate release
decisions. Prefer keeping the existing protobuf wire encoding while the Rust
API changes, then version any necessary format changes separately.

Define supported object/version identification, parameter binding, allocation
limits, invalid-field errors, and trailing-data policy at the wire boundary.
Retain the validation already present in the deserializers. Keep the current
rule that untrusted polynomial bytes cannot grant variable-time permission;
the ciphertext import boundary classifies its public data separately. Packed
plaintext storage is explicitly an in-memory optimization, not a wire format.

Secret-key export deserves an explicit method such as `export_secret_bytes`
returning `Zeroizing<Vec<u8>>`, with temporary buffers also reviewed for cleanup.
This preserves the ability to persist keys while making secret ownership
visible. It does not imply encryption or authentication of exported bytes.

Keep the existing structured `fhe::Error` and `Result` rather than replacing
them with strings or an application-oriented catch-all error. The enums already
use `thiserror` and `#[non_exhaustive]`. Consider moving detailed classification
enums under a public `error` module instead of exporting all of them at the
root. Preserve useful source errors where possible; protobuf decode errors are
currently reduced to an object classification. Separately replace
`fhe-util::sample_vec_cbd`'s `&'static str` error with a small typed error if that
function remains public. Do not remove useful diagnostics merely to shrink an
enum.

## 9. Make sensitive ownership and timing permission consistent

[`SecretKey`](../crates/fhe/src/bfv/keys/secret_key.rs) derives `Debug`, which
prints its `coeffs` field; this was confirmed without displaying the coefficients
in the review output. `Plaintext` also derives `Debug` through its polynomial.
Existing drop-time zeroization does not prevent values being copied into logs.

**Proposed:** implement redacted `Debug` for secrets, plaintexts, and builders
that contain secrets. Show public metadata such as degree and level, not
coefficients. Keep explicit inspection/export methods separate. Audit `Clone`
on secret-bearing objects and eliminate unnecessary copies such as the
evaluation-key builder's owned secret. Do not replace every key with `Arc` by
default: shared secret ownership delays destruction. Document the chosen
lifetime and retain zeroization on errors, normal drops, and unwinding.

The timing API is partly modernized:

- `PublicData` and `VariableTime` make public-data classification explicit.
- `FheEncoderVariableTime` and `Poly::try_convert_from_public` accept tokens.
- But public `rq::traits::TryConvertFrom` still accepts `variable_time: bool`,
  allowing callers to select the same behavior with an unexplained `true`.
- `SecretKey::measure_noise` is `unsafe` with a timing-only safety explanation.
- Low-level NTT `_vt` methods use raw pointers and have genuine memory-access
  preconditions; they need separate treatment.

Complete the migration: a default conversion without a timing argument, and a
named public-data path requiring the token. Remove the public boolean route.
Keep cryptographic classification distinct from Rust memory safety. After
confirming there is no additional memory-safety precondition, expose noise
measurement as an explicitly diagnostic, variable-time operation with an
acknowledgment specific to revealing secret-dependent diagnostics. Do not ask
callers to pretend that the secret-dependent noise is public input merely to
obtain the existing `PublicData` token.

Where useful, put checked slice-based wrappers around raw-pointer NTT kernels.
Validate lengths and coefficient-range requirements before the unsafe call;
retain `unsafe` on the raw kernels. Removing `unsafe` from a timing-only contract
does not justify removing real pointer safety requirements elsewhere. None of
these API changes establishes constant-time behavior of arbitrary `BigUint`
operations; document the actual supported timing guarantees.

There are manual `unsafe impl Send` declarations for `BfvParameters` and
`Plaintext`. Try removing them and let auto traits follow their fields. Check
both NTT backends; if a backend prevents automatic derivation, investigate the
actual constraint and keep any necessary unsafe proof at that boundary.
Add compile-time assertions for the intended `Send`/`Sync` promises, rather than
asserting thread safety without an explanation. This follows Rust's
[automatic `Send`/`Sync` rules](https://doc.rust-lang.org/nomicon/send-and-sync.html).

## 10. Preserve advanced paths with a smaller surrounding API

The public [`bfv` exports](../crates/fhe/src/bfv/mod.rs) mix common objects,
internal context concepts, packing, reusable scratch, and configurable
multiplication. These advanced types serve real use cases and should remain
available, with clearer grouping.

- Keep `bfv::{Parameters, Encoding, Plaintext, Ciphertext, SecretKey, PublicKey}`
  as the introductory surface. Put advanced packing and reusable evaluation
  machinery in documented submodules. Name them by function rather than
  creating an unrestricted “internals” escape hatch.
- Rename `Multiplicator` to `MultiplicationPlan`. Keep it immutable after build,
  with explicit configuration for scaling, relinearization, and modulus
  switching. A builder or config struct is clearer than the current custom
  constructor's sequence of scaling factors and basis arguments.
- Avoid deep-copying relinearization keys into plans by default. Borrow an
  immutable key, or deliberately share its immutable backing storage when an
  owned plan is needed. Keep mutable scratch separate from the plan.
- Preserve `PreparedMultiplicand`'s useful contract: it owns a snapshot and
  borrows the strategy, preventing strategy changes while the cache is used.
  Preserve the asymmetric-scaling behavior and cleanup of cached coefficients.
- Keep reusable dot-product workspaces with `&mut self`; they make exclusive
  scratch ownership explicit. Do not hide them behind global or thread-local
  mutable state. Allocating convenience methods can construct temporary scratch.

[`dot_product_scalar`](../crates/fhe/src/bfv/ops/dot_product.rs) requires
`Iterator + Clone` and repeatedly traverses cloned iterators, with a documented
requirement that every clone yield the same inputs. That is an unusually strong
semantic contract for a public iterator bound.

Prefer slices for the common checked dot product. For general borrowed inputs,
accept an `IntoIterator` and collect references once, then validate and evaluate
the same snapshot. Keep an advanced streaming accumulator if avoiding that
small allocation is important. Do not merely replace `Clone` with
`ExactSizeIterator`: a length promise does not make an iterator replayable.
Preserve count mismatch errors, validation-before-mutation, timing propagation,
and packed views. The choice of validation and snapshot strategy also matters
for unsafe arithmetic kernels.

[`PlaintextVec`](../crates/fhe/src/bfv/plaintext_vec.rs) primarily exists to
implement the encoding trait, and the single-plaintext encoder currently builds
one then clones its first element. Replace it with an inherent `encode_chunks`
returning `Vec<Plaintext>`, or a `PlaintextBatch` only if it owns useful shared
metadata/length invariants. Encode a single plaintext directly using the common
chunk kernel. Define whether an empty batch input yields zero rows or one zero
plaintext; today it yields one.

Keep [`PackedPlaintextVec`](../crates/fhe/src/bfv/packed_plaintext.rs), possibly
renamed `PackedPlaintextBatch`, because its contiguous storage and validated
shared context are substantive. Add conventional borrowed iteration where
useful. Do not add unrestricted `FromIterator` or `Extend` if inserts can fail
context validation; use `try_from_iter`/`try_extend` instead. Preserve buffer
reuse and zeroization when simplifying wrappers.

Coordinate names with the existing [bounded Rayon pool proposal](bounded-rayon-pool.md).
That is a separate proposal, not an implemented API. The redesign here should
continue to permit caller-owned pools and independent mutable workspaces; it
does not require introducing an executor or implicit parallel execution.

## 11. Tighten math extension points and retain multiparty typestate

[`Poly<PowerBasis>`, `Poly<Ntt>`, and `Poly<NttShoup>`](../crates/fhe-math/src/rq/mod.rs)
already prevent many representation mistakes statically. Keep the consuming
`into_ntt` / `into_power_basis` transitions and explicit borrowed conversions.
This is a good use of Rust's type system.

`RepresentationTag` and `ScaleRepresentation` are publicly implementable even
though algorithms know a fixed set of representations and legal transitions.
Seal these traits unless downstream representation implementations are an
intentional supported feature. Runtime `Representation` can remain for
inspection and wire validation; its existence does not mean the typed design
should be undone. Avoid adding unrelated typestate parameters for every runtime
setting.

Other math improvements should be narrow: replace the duplicated conversion
trait with descriptive constructors; distinguish coefficients from raw NTT/RNS
residues in their names; keep validated range/layout information at the boundary
of unsafe kernels. Slice conveniences can reduce unnecessary `ndarray` exposure
where contiguous input is required, but retain explicit multidimensional/strided
views where they are actually useful. Audit the public bit-transcoding helpers'
invalid-width, truncation, and padding behavior; prefer checked boundary
functions and named internal prevalidated kernels to hidden debug-only checks.

In [`mbfv`](../crates/fhe/src/mbfv/mod.rs), the sealed `Round` markers and
`RelinKeyShare<R>` already distinguish protocol rounds. Retain that design.
[`RelinKeyShare`](../crates/fhe/src/mbfv/relin_key_gen.rs) nevertheless stores
`last_round: Option<Arc<RelinKeyShare<R1Aggregated>>>` for all rounds and can emit
`MissingRelinearizationRoundOneShare`. Store required data in the corresponding
state, for example in distinct share structs or associated round-specific
payload types, so a constructed round-two share always has its dependency.

Consider consuming generator state across rounds if the protocol requires
one-use ephemeral state; decide this from protocol semantics, not aesthetics.
Retain distinct ciphertext/evaluation-key levels and validation of share
compatibility. Round marker types alone do not prove shares belong to one
protocol session.

`Aggregate` is more defensible than the empty FHE markers because it expresses
a real fallible operation. It can remain a scoped trait or become inherent
`aggregate(shares)` methods. The blanket `AggregateIter` convenience trait is
optional, and the `Aggregate<Result<S>>` adapter currently collects to a vector
before aggregation. Prefer explicit streaming error propagation if it preserves
the protocol's validation behavior. Do not substitute `FromIterator` on the
aggregate object when aggregation can fail. Keep this work behind the existing
experimental feature and preserve its documented protocol limitations.

## Implementation order and acceptance criteria

1. **Contract fixes.** Add regression tests for equality, zero handling,
   redacted debugging, context/descriptor mutation, and invalid defaults.
   Close mutation before relying on stronger constructor invariants elsewhere.
   Preserve existing arithmetic, timing, and serialization tests.
2. **Concrete facade and parameters.** Introduce inherent methods, the shared
   `Parameters` handle, and the chosen compatibility relation. Port the basic
   example and README first. Ensure an external consumer can use the full basic
   flow with `fhe` and `rand`, without `fhe-traits` or `fhe-math` imports.
3. **Encoding and checked evaluation.** Remove optional encoding metadata and
   the empty ciphertext sentinel. Complete checked operations and their
   in-place error guarantees. Port signed/large-plaintext tests, RGSW, voting,
   and PIR examples with explicit semantics.
4. **Construction and serialization.** Consolidate builders, expose inherent
   byte methods, hide wire DTOs, and remove duplicate conversion traits.
   Preserve stored-byte compatibility by default; if a format must change,
   specify its version and migration fixtures before changing the writer.
5. **Advanced and package cleanup.** Move advanced types into focused modules,
   simplify batching/iteration, seal representation traits, and retire obsolete
   traits after moving shared timing tokens. Tighten multiparty state separately.
6. **Release the breaking API.** Remove transitional aliases in this release
   rather than maintaining two permanent facades. Update all four crate READMEs,
   rustdoc, examples, benchmarks, and migration notes together.

The following checks should gate implementation:

| Area | Required evidence |
| --- | --- |
| Equality | Reflexivity, symmetry, transitivity, metadata policy, and level/parameter distinctions |
| Validity | Compile-fail coverage for sealed traits and removed mutation access; checked component imports reject mismatches |
| Arithmetic | All operand kinds, levels and part counts; allocating/in-place agreement; receiver unchanged on returned errors |
| Parameters | Clones and separately built/deserialized equal specs follow the chosen policy; unequal specs are rejected |
| Encoding | Empty/short/full/oversized inputs; signed boundaries; `BigUint` overflow; unsupported SIMD; output padding and all supported levels |
| Persistence | Round trips for every public constructible value; malformed input and limits; old format fixtures if compatibility is promised; seed invalidation |
| Sensitive data | Redacted debug; timing permission does not come from untrusted bytes; cleanup of retained and temporary buffers |
| Performance | Existing multiplication, prepared-product, packed-dot-product, allocation-path, and PIR benchmarks with setup costs and peak retained memory distinguished |

Use small external-consumer tests and rustdoc examples to test discoverability
and coercion behavior. Test important rejected programs with compile-fail
examples instead of relying on prose. Verify intended auto traits for both NTT
backends. Compare arithmetic outputs as appropriate to each operation: exact
bytes for deterministic representation-preserving changes, decoded results for
randomized encryption, and existing exact-reference tests for fused rounding.

Before committing implementation changes, follow `AGENTS.md`:

```sh
cargo test
cargo +nightly fmt --all
cargo clippy --all-targets -- -D warnings
```

Also test the optional APIs and the documented MSRV explicitly:

```sh
cargo test --workspace --all-features
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo +1.91.1 check --workspace --all-targets
RUSTDOCFLAGS='-D warnings' cargo doc --workspace --all-features --no-deps
```

## Review evidence and limits

The recommendations above are based on the current implementations, not on an
assumption that the original author's experience describes today's code.
Source links point to the relevant files; symbols are named so the references
remain useful when line numbers change.

For this documentation review:

- The default workspace tests were run with `cargo test --offline` and passed.
- `cargo test --offline --workspace --all-features` also passed, covering the
  accelerated NTT backend and experimental multiparty feature. These runs used
  the installed default nightly toolchain; they do not establish MSRV support.
- Five temporary external-consumer probes passed, confirming the current
  non-transitive plaintext comparison, the zero sentinel's failed round trip
  and decryption, rejection of equal independent parameter sets, mutable
  component invalidation, and secret coefficients appearing in `Debug`.
  These probes assert the existing behavior; they are not fixes or regression
  tests for the proposed behavior.
- No performance comparisons, constant-time verification, or proof of backend
  auto traits were performed. Recommendations about those properties include
  explicit implementation gates rather than claims of measured improvement.
- No library source was changed. Proposed signatures are intentionally marked
  `ignore`; convert them to compiling tests as the API is implemented.

The design choices—such as removing operators, accepting equal independently
built parameters, and using consuming builders—are recommendations for this
crate. They are not requirements imposed by the Rust language. The specific
trait-contract issue with equality and the observed invalid-state behaviors
are separate, reproducible findings.
