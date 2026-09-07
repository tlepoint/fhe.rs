# P0 API contract changes

The later [P1 migration](api-p1-migration.md) supersedes the P0 signatures below.

This implements the P0 work in the [API proposal](api-simplification-proposal.md).
The changes intentionally break Rust source compatibility. The parameter
ownership model, encoding API, ordinary operators, and protobuf wire format
remain as before; their broader redesign belongs to P1/P2.

## Ciphertexts

`Ciphertext::zero` is removed. Use
`Ciphertext::trivial_zero(&parameters, level)?` to construct a valid two-component
zero at an explicit level. It can be decrypted, evaluated, switched and
serialized. Its value is public and deterministic; encrypt a zero plaintext
with a key and RNG when a randomized encryption is required.

Empty accumulators should be `Option<Ciphertext>` or start from their first
input. A two-component trivial zero cannot stand in for an arbitrary level or
component count. Addition still requires matching levels and component counts.
The fused product accumulator retains its separate, single-rounding semantics.

`Ciphertext::new(parts, &parameters)` becomes
`Ciphertext::from_components(parts, &parameters)`. This checks the component
count, matching contexts, membership in the parameters' chain, and canonical
residues. Lazy NTT polynomials are rejected with
`CiphertextError::NonCanonicalPolynomial`.

`Ciphertext` no longer implements `Deref` or `DerefMut`:

| Old usage | Replacement |
| --- | --- |
| `ct.len()` | `ct.component_count()` |
| `ct.iter()` | `ct.components().iter()` |
| `ct[i]` for inspection | `ct.components()[i]` |
| Mutable indexing / slice mutation | Consume components, modify them, then reconstruct with validation |
| Inferring the level from a component context | `ct.level()` |
| No direct parameter accessor | `ct.parameters()` |

For example:

```rust,ignore
let parameters = ct.parameters().clone();
let mut parts = ct.into_components();
parts[0] = -&parts[0];
let ct = Ciphertext::from_components(parts, &parameters)?;
```

Reconstruction discards the old compression seed. Unmodified component access
is immutable, so it does not invalidate the seed. Internal arithmetic continues
to invalidate seeds when coefficients change, including relinearization.

Ciphertext equality compares parameters, level and components, excluding the
compression seed. Two equivalent seeded and expanded representations compare
equal even when their byte encodings differ. Equality does not determine whether
different randomized ciphertexts encrypt the same message.

Expansion builds its output incrementally instead of allocating empty
ciphertext placeholders. Its output order, scaling and timing permissions are
preserved, including requests whose size is not a power of two.

## Plaintext equality and debugging

Plaintext equality includes the entire `Option<Encoding>`. Unknown encoding
metadata after decryption is distinct from known encoding metadata. For
encryption/decryption checks, compare decoded values using the intended encoding:

```rust,ignore
let decrypted = sk.try_decrypt(&ct)?;
assert_eq!(
    Vec::<u64>::try_decode(&decrypted, encoding.clone())?,
    Vec::<u64>::try_decode(&plaintext, encoding)?,
);
```

This removes the wildcard comparison that violated transitivity. Encoding
metadata itself is retained until the P1 encoding redesign.

`Debug` on `SecretKey`, `Plaintext`, and `Poly<R>` now shows public metadata
without coefficients. This also redacts nested secret keys in evaluation-key
builders. Explicit decoding, coefficient inspection and serialization retain
their existing meanings.

`EvaluationKeyBuilder<'a>` borrows `&'a SecretKey` instead of cloning secret
coefficients. Keep the secret key alive while building. The builder no longer
implements `Zeroize` / `ZeroizeOnDrop`; it owns only public configuration and
a borrow. The secret key retains its existing zeroization on drop.

## Contexts and mathematical values

`ContextLevel::poly_context` is now accessed through `poly_context()` and cannot
be replaced through the public API. `SubstitutionExponent::exponent` becomes
`exponent()`, keeping its value consistent with the cached permutation.
`CipherPlainContext` is crate-private.

The invalid `Default` implementations are removed from `rq::Context`,
`RnsContext`, `Scaler`, `Switcher`, `RnsScaler`, and `ScalingFactor`. Construct
them through their explicit constructors; `ScalingFactor::one()` remains the
identity. The unused internal multiplication-parameter default is also removed.

`BfvParameters` and `Plaintext` use automatically derived thread-safety traits
instead of manual `unsafe impl Send`. Compile-time tests check the intended
`Send + Sync` properties with both NTT backends.

## Timing permission

The math conversion trait no longer accepts a boolean:

```rust,ignore
// Default: variable-time processing disabled.
let private = Poly::<PowerBasis>::try_convert_from(values, &context)?;

// The caller explicitly classifies these values as public.
let public = Poly::<PowerBasis>::try_convert_from_public(
    values,
    &context,
    VariableTime::new(PublicData::assert_public()),
)?;
```

Implementations of `rq::traits::TryConvertFrom<T>` now implement
`try_convert_from_with_timing(value, context, Option<VariableTime>)`. The default
`try_convert_from` method forwards `None`. The optional-token method also permits
forwarding existing permission without introducing a public boolean shortcut.
Wire data still cannot grant variable-time permission to a polynomial.

`SecretKey::measure_noise` is replaced by the safe, explicitly diagnostic
`measure_noise_vartime`:

```rust,ignore
let noise = sk.measure_noise_vartime(
    &ct,
    fhe_traits::SecretDependentDiagnostics::acknowledge_leakage(),
)?;
```

The acknowledgment accepts secret-dependent output and timing leakage; it does
not assert that the secret-dependent noise is public input. The method first
validates through decryption. Raw-pointer NTT kernels retain their memory-safety
requirements. These interface changes do not establish new constant-time
guarantees for the arithmetic or arbitrary big-integer operations.

## Verification

Regression coverage includes equality laws, trivial zeros at all levels,
validated component imports, seeded round trips and relinearization, partial
expansion bytes and permissions, debug redaction, timing-token behavior, and
automatic thread-safety traits. Compile-fail examples cover removed mutation,
invalid defaults, the old boolean conversion, missing diagnostic acknowledgment,
and a builder outliving its borrowed key.

Validation completed: `cargo test` passed 263 tests; the workspace all-feature
suite passed 274 tests, including experimental multiparty BFV and the optional
NTT backend. `cargo +nightly fmt --all`, default and all-feature Clippy with
warnings denied, and workspace all-feature rustdoc with warnings denied also
passed. These checks used the installed nightly toolchain. This refactor does
not add a performance or constant-time audit.
