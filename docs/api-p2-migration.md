# P2 migration: advanced Rust APIs

P1 was committed as `ca81f68`. This document describes the subsequent P2 changes
from the [API proposal](api-simplification-proposal.md). They intentionally break
Rust source compatibility. Read the [P1 migration](api-p1-migration.md) first for
the common BFV flow, parameter compatibility, encoding, and checked arithmetic.
The proposal's [coverage audit](api-simplification-proposal.md#implementation-coverage-2026-09-07)
lists remaining work and optional suggestions across all three priorities;
this migration does not claim that every proposal item is complete.

## Imports

The introductory types remain in `fhe::bfv`: `Parameters`, `Encoding`,
`Plaintext`, `Ciphertext`, `SecretKey`, and `PublicKey`. Parameter builders and
profiles also remain there. Advanced types have dedicated public modules:

| Module | Exports |
| --- | --- |
| `bfv::evaluation` | `EvaluationKey`, `EvaluationKeyBuilder`, `RelinearizationKey`, `RgswCiphertext`, `MultiplicationPlan`, `MultiplicationPlanBuilder`, `MultiplicationScaling`, `PreparedMultiplicand`, `CiphertextProductAccumulator`, `DotProductScalarWorkspace`, `dot_product_scalar`, `dot_product_scalar_iter` |
| `bfv::packing` | `PackedPlaintext`, `PackedPlaintextView`, `PackedPlaintextBatch`, `PackedPlaintextIter` |
| `bfv::context` | Immutable `ContextLevel` descriptors |

Old root exports are removed. The context module exposes descriptors, not mutable
caches or an unrestricted internal API. No executor or automatic parallelism is
introduced; workspaces remain independent mutable values suitable for use in
caller-owned pools. The [bounded Rayon pool proposal](bounded-rayon-pool.md)
remains separate work.

## Immutable multiplication plans

Configure a plan before building it:

```rust
use fhe::bfv::{Parameters, evaluation::{MultiplicationPlan, RelinearizationKey}};

fn plan<'key>(parameters: &Parameters, key: &'key RelinearizationKey)
    -> fhe::Result<MultiplicationPlan<'key>>
{
    MultiplicationPlan::builder(parameters)
        .level(0)
        .relinearization(key)
        .modulus_switching(true)
        .build()
}
```

The default uses the parameters' precomputed BFV multiplication basis and
scaling, without relinearization or modulus switching. The custom `new` and
`new_leveled` constructors are replaced by `.extended_basis(moduli)` and
`.scaling(MultiplicationScaling { left, right, product })`. Each setting can be
overridden independently. Custom factors and basis sizes are still an expert
facility; the builder checks structural validity, not the suitability of a
rounding strategy for a particular circuit.

`MultiplicationPlan::with_relinearization(key)` remains a convenience for the
default strategy at the key's ciphertext level. Replace
`without_relinearization(parameters, level)` with
`builder(parameters).level(level).build()`.

There are no plan mutators. Rebuild from a cloned builder when comparing
configurations. Validation rejects incompatible keys, mismatched ciphertext
levels, invalid bases, and switching past the last level. Key and ciphertext
levels remain distinct; the key level need not equal the input level.

Plans borrow their relinearization keys instead of cloning key polynomials.
A stored plan consequently carries a lifetime and its key must outlive it.
`PreparedMultiplicand` continues to own a coefficient snapshot and borrow the
plan. Prepared coefficients are cleared on drop. Asymmetric scaling, square
specialization, two-component validation, and separate mutable workspaces are
preserved. Compile-time tests cover plan immutability and the intended
`Send + Sync` properties on both NTT backends.

## Dot products: slices and one-pass iterators

The BFV free function and workspace methods use slices:

```rust
use fhe::bfv::{Ciphertext, Plaintext, evaluation::DotProductScalarWorkspace};

fn evaluate(workspace: &mut DotProductScalarWorkspace,
            query: &[Ciphertext], rows: &[Plaintext]) -> fhe::Result<Ciphertext> {
    workspace.dot_product_scalar(query, rows)
}
```

Workspace `_refs` methods accept slices of references without allocating an
operand list. General `_iter` methods accept `IntoIterator`, collect references
once, and validate and evaluate that same snapshot. They no longer require
`Clone` or a promise that repeated traversals yield identical values. Packing
follows the same convention: `dot_product_scalar_packed` accepts ciphertext
slices and packed-view slices; `_packed_refs` and `_packed_iter` handle reference
slices and general borrowed inputs respectively. The latter accepts `&batch` or
borrowed `PackedPlaintext` values directly.

The math `rq::dot_product` and `DotProductWorkspace` follow the same conventions,
including `_refs`, `_iter`, and `dot_product_into` methods. Slice entrypoints do
not allocate operand lists; `_iter` entrypoints retain two lists of references
(or packed views) for the duration of the call. Ciphertexts and plaintext
coefficients are not cloned. In-place math output and workspace storage are
unchanged on validation errors. Public iterators that panic do so while inputs
are being collected, before arithmetic starts.

Empty inputs, unequal counts, incompatible contexts or levels, and inconsistent
ciphertext component counts still return errors. Timing permission is recomputed
from every operand. Fused accumulation, periodic reduction, the long-product
fallback, packed decoding, and scratch cleanup retain their existing bounds.

## Ordinary and packed batches

`PlaintextVec` is removed. Use `Plaintext::encode_chunks`, returning
`Vec<Plaintext>`, or the `encode_chunks_signed`, `encode_chunks_biguint`, and
`encode_chunks_public` variants. Each has an `_at_level` form.

Empty input produces **one zero plaintext**, preserving the previous behavior.
The final chunk is zero-padded to the polynomial degree. Single-plaintext
encoding calls the common chunk kernel directly. Each returned plaintext clears
its polynomial on drop, including when an enclosing ordinary vector is dropped.

`PackedPlaintextVec` becomes `PackedPlaintextBatch`. It retains its shared
parameters, level, contiguous coefficient storage, and per-row timing policy.
Use `try_from_iter(parameters, level, plaintexts)`, `try_extend(plaintexts)`,
`push`, `get`, `iter`, or borrowed `for row in &batch` iteration. An empty packed
batch has zero rows. `try_extend` snapshots and validates all rows before
appending; errors preserve both contents and allocation. `clear` clears bytes
while retaining capacity. Growth clears the old coefficient allocation before
releasing it. There is no infallible `FromIterator` or `Extend` implementation.

## Math constructors and checked transcoding

`RepresentationTag` and `ScaleRepresentation` are sealed. The supported markers
remain `PowerBasis`, `Ntt`, and `NttShoup`; consuming representation transitions
are unchanged. The public `rq::traits::TryConvertFrom` trait is removed. Private
wire decoding uses a contextual helper that is not a downstream extension point.

| Input meaning | Constructor |
| --- | --- |
| Unsigned integer coefficients | `Poly::<PowerBasis>::from_coefficients(&values, context)` |
| Signed integer coefficients | `Poly::<PowerBasis>::from_signed_coefficients(&values, context)` |
| Arbitrary unsigned integer coefficients | `Poly::<PowerBasis>::from_biguint_coefficients(&values, context)` |
| Raw residue matrix in representation `R` | `Poly::<R>::from_rns_residues(array, context)` |
| Contiguous modulus-major raw residues in `R` | `Poly::<R>::from_rns_slice(&values, context)` |
| Wide NTT accumulators, including strided views | `Poly::<Ntt>::from_wide_ntt_residues(view, context)` |

Each has a `_with_timing` variant accepting `Option<VariableTime>`; `Some(token)`
is explicit public-data classification. Ordinary constructors use `None`.
Big-integer arithmetic still carries no constant-time guarantee.

Integer coefficient input always means coefficients, is reduced modulo each
prime, and is padded to the degree. More than `degree` coefficients is an error;
length no longer silently changes the interpretation to raw RNS rows. Raw
residue constructors require exact shape/length, reduce each modulus row, and
normalize owned arrays into standard layout. They do not perform an NTT. The
unused conversion from `&Plaintext` through the removed math trait is also gone;
BFV arithmetic and explicit message decoding remain inherent methods.

`fhe-util` bit-transcoding functions now return `Result<_, TranscodeError>`.
Invalid widths, size overflow, and input words that would truncate are errors.
Packing appends only after validation. The final byte is zero-padded.
`transcode_from_bytes_exact(bytes, width, count)` rejects extra/missing bytes and
nonzero padding. The stream decoder `transcode_from_bytes(bytes, width)` uses
every input bit and zero-extends the final partial word; it cannot infer the
original count from padded bytes. `transcode_bidirectional` likewise preserves
all input bits and zero-extends its final word. Masking an input must now be an
explicit caller decision. Modulus `serialize_vec` and `deserialize_vec` propagate
transcoding errors and document that they do not validate residues against the
modulus itself. Polynomial import still checks canonical residues separately.

The public centered-binomial sampler now returns `InvalidVariance { actual }`
instead of a string error. Its algorithm and RNG consumption are unchanged.

## Experimental multiparty state

`RelinKeyShare<R>` retains sealed round markers. Its dependency is now an
associated round-specific type: `()` for first-round shares and their aggregate,
and a required `Arc<RelinKeyShare<R1Aggregated>>` for second-round shares.
`MissingRelinearizationRoundOneShare` is removed because the state cannot be
constructed without its dependency.

The generator transitions are consuming:

```rust,ignore
let generator = RelinKeyGenerator::new(&secret_share, &crp, &mut rng)?;
let (round_one_share, round_two_state) = generator.round_1(&mut rng)?;
// Exchange and aggregate all first-round shares here.
let round_two_share = round_two_state.round_2(&aggregate, &mut rng)?;
```

The design treats one generator as one protocol execution: its sampled `u` is
retained between the two rounds and cleared when the second-round state is
consumed or dropped. Retransmit by cloning an already generated share, and
create a fresh generator for a new execution. This choice follows the ephemeral
key flow also described by the authors' [Lattigo protocol documentation](https://pkg.go.dev/github.com/tuneinsight/lattigo/v6/multiparty#RelinearizationKeyGenProtocol).
It is an API ownership decision, not a new cryptographic security proof.

Relinearization generation and aggregation reject incompatible parameters,
polynomial contexts, or share shapes. Second-round aggregation also rejects
shares with different first-round aggregates; structurally equal independent
aggregates are accepted. These checks and round markers do not authenticate
participants or establish application session identity. Generation remains at
level zero; this change does not introduce new multiparty level support.

`Aggregate` and `AggregateIter` remain scoped to `mbfv`. The
`Aggregate<Result<S>>` adapter streams shares instead of collecting them. It
stops on the first input or aggregation error it encounters; an observed input
error takes precedence over a partial aggregate result. It need not consume the
rest of a failing input, unlike the previous eager adapter.

Everything remains behind `experimental-mbfv`. Existing warnings about incomplete
common-reference-string generation, noise flooding, and lack of an independent
audit still apply.

## Persistence and verification

The protobuf writers and accepted canonical representation are unchanged. The
P0 binary fixtures continue to cover import, re-export, and evaluation on their
matching backends. The pre-existing limitation on importing stored polynomial
representations across NTT backends remains; see the P1 notes.

Regression coverage includes slice/iterator arithmetic agreement, single-pass
visits, timing propagation, empty and signed chunking at every level, packed
batch error atomicity, numeric/raw constructor distinctions, bit widths 1–64,
padding and truncation rejection, immutable plans and sealed traits, and
multiparty round compatibility. Existing tensor, square, asymmetric/prepared
multiplication, serialization, and scratch-cleanup tests are retained.

The following checks passed:

```sh
cargo test
cargo test --workspace --all-features
cargo +nightly fmt --all
cargo clippy --all-targets -- -D warnings
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo +1.91.1 check --workspace --all-targets --all-features
RUSTDOCFLAGS='-D warnings' cargo doc --workspace --all-features --no-deps
```

The final iterator-unwind and fused error-adapter regressions were additionally
checked with the native packed-dot-product tests and
`cargo test -p fhe --all-features --lib`.

Criterion smoke runs (`cargo bench -p <crate> --bench <name> -- --test`) passed:
14 multiplication/preparation cases, 17 math allocation-path cases, and six PIR
expansion cases. Benchmarks use slice entrypoints where contiguous inputs are
available; iterator and packed-view paths remain exercised for strided databases.
These smoke runs verify execution and the benchmarks' arithmetic assertions;
they are not statistically meaningful before/after performance measurements.
Full dense PIR database throughput and peak process memory were not measured.

Both `mulpir` and `sealpir` completed retrieval assertions in release mode with
128 records of 256 bytes, using packed and unpacked storage. Retained coefficient
storage was 262,144 bytes unpacked in both examples, versus 215,048 bytes for
packed MulPIR and 147,464 bytes for packed SealPIR (including padding). These
figures exclude parameter/key caches, allocator capacity, and transient scratch.
The examples separately report parameter generation, preprocessing, query,
response, and decoding time; no timing comparison is claimed here.
