# P1 API migration

The [P2 migration notes](api-p2-migration.md) describe subsequent advanced API,
batch, math-constructor, and experimental multiparty changes.

Implemented after the P0 contract commit `0eb11ce`. This is a breaking Rust API
change. Existing protobuf bytes remain readable with the same NTT backend; fixtures produced by that
commit exercise parameter, key, ciphertext, and evaluation-key imports.

## Ordinary use

Encoding, encryption, decryption, arithmetic, and supported serialization are
inherent methods. Applications do not need an FHE trait import or `fhe-math`
for ordinary BFV operations.

```rust
use fhe::bfv::{Encoding, Parameters, Plaintext, PublicKey, SecretKey};

fn main() -> fhe::Result<()> {
    let parameters = Parameters::builder()
        .degree(2048)
        .plaintext_modulus(1024_u64)
        .ciphertext_moduli([0x3fffffff000001])
        .build()?;
    let mut rng = rand::rng();
    let sk = SecretKey::generate(&parameters, &mut rng);
    let pk = PublicKey::from_secret_key(&sk, &mut rng);
    let a = Plaintext::encode(&parameters, &[20], Encoding::Polynomial)?;
    let b = Plaintext::encode_signed(&parameters, &[-7], Encoding::Polynomial)?;
    let ca = sk.encrypt(&a, &mut rng)?;
    let cb = pk.encrypt(&b, &mut rng)?;
    let product = ca.multiply(&cb)?;
    assert_eq!(sk.decrypt(&product)?.decode_signed(Encoding::Polynomial)?[0], -140);
    Ok(())
}
```

The speculative `Fhe*` hierarchy and custom serialization traits have been
removed, along with the `fhe-traits` workspace package. Remove that dependency
from applications. Timing and diagnostic tokens now live in `fhe-util` and are
re-exported by `fhe` and `fhe-math`.

## Shared parameters and builders

`BfvParameters` and `BfvParametersBuilder` become `Parameters` and
`ParametersBuilder` in `fhe::bfv`. A `Parameters` value privately owns shared
precomputation. Cloning it is cheap; constructors take `&Parameters` and retain
a cloned handle. Do not wrap parameters in an additional `Arc`.

Parameter compatibility uses identity as a fast path, then compares degree,
ordered ciphertext moduli, plaintext modulus, and noise variance. Independently
built and deserialized equivalent parameters are compatible. Derived caches are
excluded. This applies to arithmetic, key operations, imports, prepared
multiplication, accumulators, and packed plaintexts. Compatible parameters do
not imply that two ciphertexts use the same secret key.

Parameter builders consume and return themselves; `build(self)` produces the
shared handle. Required fields are explicitly optional until configured. The
most recent modulus setter replaces the previous choice.

| Before | Now |
| --- | --- |
| `set_degree(n)` | `degree(n)` |
| `set_variance(v)` | `noise_variance(v)` |
| `set_plaintext_modulus(t)` / `_biguint(t)` | `plaintext_modulus(t)` accepting `Into<BigUint>` |
| `set_moduli(q)` | `ciphertext_moduli(q)` |
| `set_moduli_sizes(bits)` | `ciphertext_modulus_bits(bits)` |
| `build_arc()` | `build()` |
| `plaintext()` | `plaintext_modulus_u64() -> Option<u64>` |
| `plaintext_big()` | `plaintext_modulus() -> &BigUint` |

Use typed literals, such as `1024_u64`, with the generic plaintext setter.
Builder validation retains the existing arithmetic restrictions; it does not
estimate security for arbitrary configurations.

Use `Parameters::profile_128(degree, plaintext_bits)` to build one existing
preselected profile. `Parameters::profiles_128(bits)` lists lightweight
`ParameterProfile` descriptors that can be filtered before calling `build()`.
Availability depends on the plaintext prime; selecting by the former `.nth()`
position is unnecessary. The existing profile values and security assumptions
are retained.

Evaluation-key configuration borrows the secret and is also consuming:

```rust,ignore
let key = EvaluationKey::builder(&sk)
    .ciphertext_level(1)
    .key_level(0)
    .enable_row_rotation()
    .enable_column_rotation(1)
    .enable_inner_sum()
    .enable_expansion(4)
    .build(&mut rng)?;
```

`EvaluationKeyBuilder::new(&sk)` is infallible. The setters return `Self`;
`build` validates levels, rotation steps, and expansion size before generating
keys. There is no `new_leveled` constructor for this builder. Reassign a builder
when configuring it conditionally, or clone its public configuration to build
multiple keys; cloning the builder does not copy secret coefficients.

## Encoding, decoding, and numeric types

`Encoding` is the copyable enum `Encoding::Polynomial | Encoding::Simd`. It
contains no level. `Plaintext` stores its polynomial and parameters, without
optional encoding metadata. Its level comes from the polynomial context.
Equality therefore agrees before and after decryption.

- `Plaintext::encode(&parameters, &[u64], encoding)` uses level zero.
- `encode_at_level(&parameters, values, encoding, level)` selects a level.
- `encode_signed` / `encode_signed_at_level` accept `&[i64]`.
- `encode_biguint` / `encode_biguint_at_level` accept `&[BigUint]`.
- `encode_public` / `encode_public_at_level` accept unsigned input and a
  `VariableTime` token, re-exported from `fhe`.
- `Plaintext::zero(&parameters, level)` needs no encoding choice.
- `plaintext.decode(encoding)`, `decode_signed(encoding)`, and
  `decode_biguint(encoding)` return explicitly chosen integer representations.

Inputs are reduced modulo the plaintext modulus before conversion and padded
with zeros. A single plaintext rejects more than `degree` inputs. This also
fixes unsigned/large inputs above the modulus being transformed before their
plaintext reduction. Signed decoding uses centered representatives, with an
even modulus's midpoint represented negatively. Values that do not fit the
requested output integer type return a structured error.

Decoding always requires an interpretation and returns exactly `degree` values;
the input length is not recoverable. Polynomial and SIMD interpretations may be
requested on the same plaintext, and decoding does not change its level. SIMD
still requires a machine-word prime plaintext modulus congruent to 1 modulo
`2 * degree`; BigUint input does not extend that backend's capabilities.

`PlaintextVec` retains its existing chunked collection role with inherent
`encode`, `encode_biguint`, and `encode_public`, plus their `_at_level` forms.
Empty input produces one zero plaintext and the final chunk is padded.
Packed storage retains the polynomial, level, and local timing permission.
Collection and advanced-module redesigns remain P2 work.

## Checked arithmetic and clearer names

Ciphertext binary arithmetic and assignment operators are removed. Unary
negation remains infallible. Use methods returning `fhe::Result`:

| Operand | Allocating | In place |
| --- | --- | --- |
| Ciphertext | `add`, `subtract`, `multiply` | `add_assign`, `subtract_assign`, `multiply_assign` |
| Plaintext | `add_plaintext`, `subtract_plaintext`, `multiply_plaintext` | corresponding `_assign` methods |
| RGSW ciphertext | `multiply_rgsw` | `multiply_rgsw_assign` |
| Self | `square` | `square_assign` |

Errors leave an in-place receiver unchanged, including its compressed seed.
Addition and subtraction require equal component counts. Unrelinearized
multiplication accepts `m` and `n` components and produces `m + n - 1` at the
same level. Relinearization and modulus switching remain explicit.

`Multiplicator` becomes `MultiplicationPlan`;
`Multiplicator::default(key)` becomes
`MultiplicationPlan::with_relinearization(key)`. Configured multiplication still
requires two-component operands. Prepared operands, dedicated squaring,
workspaces, and bounded product accumulation retain their optimized paths.
Accumulating products with one final rounding is still distinct from adding
separately rounded products.

`SecretKey::random` becomes `generate`; `PublicKey::new` becomes
`from_secret_key`. `RGSWCiphertext` becomes `RgswCiphertext`, and alternate
secret-key encryption is explicitly `encrypt_rgsw`. Evaluation verbs become
`relinearize`, `rotate_rows`, `rotate_columns`, `inner_sum`, and `expand`.

## Serialization boundaries

Supported objects have inherent `to_bytes` and `from_bytes`. Contextual BFV
imports take `&Parameters`; math polynomial imports still take an explicit math
context. Parameter imports are self-contained. Generated `fhe::proto` types are
private, and BFV DTO conversions use a private `FromProto` helper. The public
math numeric conversion extension trait remains for the P2 math-surface review.

The existing protobuf format has no new envelope or authentication. Native and
`tfhe-ntt` builds are not generally interoperable with each other's existing
payloads. This was reproduced on the P0 commit as well as P1; fixtures from each
backend verify preservation separately. Keep the backend consistent when
persisting or exchanging these bytes. Adding a portable, versioned format is a
separate change. The caller
selects the expected object type and supplies its parameter binding. Unknown
protobuf fields remain accepted; malformed trailing bytes are rejected, while
trailing valid protobuf fields follow protobuf merge semantics. Import validates
stored dimensions, representations, levels, component counts, and parameter
constraints. Untrusted polynomial bytes cannot grant variable-time permission.
There is no newly introduced global byte budget: callers receiving network data
must bound the input before decoding. Packed plaintexts remain an in-memory
optimization rather than a wire format.

Secret keys use `export_secret_bytes() -> Zeroizing<Vec<u8>>` instead of
`to_bytes()`. The returned buffer is unencrypted secret material and is cleared
on drop. Temporary protobuf coefficient vectors are also cleared on drop,
including failed and partially decoded imports. Caller-owned input bytes remain
the caller's responsibility. `SecretKey::from_bytes(bytes, &parameters)` retains
validated import.

## Compatibility measurement

A local release benchmark at degree 4096 (20-bit plaintext profile) measured
checked allocating addition at about 12.3 microseconds with shared parameters
and 57.0 microseconds with independently constructed equal parameters. The
fallback also reaches existing polynomial-context value comparisons. These are
single-machine measurements, not performance guarantees. Prefer cloning a
shared `Parameters` handle within an application; independently built handles
remain a supported interchange boundary. Reproduce with
`cargo bench -p fhe --bench bfv_core -- parameter_compatibility`.
