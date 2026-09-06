# fhe-math: optimizations and simplifications

Reviewed the source and existing benchmarks at commit `20506cd`, including the recent correctness fixes. The initial review identified implementation opportunities without measured speedups. Items 1–8 and the forward portion of item 10 have since been implemented; their original rationale remains below, with measurements in the implementation results sections. References below are repository-relative source paths and symbol names.

Let **N** denote polynomial degree and **L** the number of RNS moduli. Start with allocations and context construction; treat arithmetic kernel changes as benchmark-driven experiments.

## Suggested order

| Order | Opportunity | Main benefit | Effort / risk |
| --- | --- | --- | --- |
| 1 | Fill Shoup and random buffers directly | Fewer allocations and copies | Small / low |
| 2 | Borrow contexts during traversal | Avoid deep table clones | Small / low |
| 3 | Use word-sized coprimality checks | Cheaper context setup | Small / low |
| 4 | Borrow coefficients during serialization | Less copying and temporary storage | Medium / low |
| 5 | Consolidate polynomial arithmetic and constructors | Smaller invariant-maintenance surface | Medium / moderate |
| 6 | Share NTT tables across context levels | Lower setup time and retained memory | Medium / moderate |
| 7 | Keep only the selected NTT backend | Lower setup time and retained memory | Medium / moderate |
| 8 | Streamline dot-product validation and scratch use | Less setup overhead and allocation | Medium / moderate |
| 9 | Add a tiled polynomial-scaling kernel | Better memory locality | Large / needs measurement |
| 10 | Fuse native NTT final passes | Fewer full-buffer passes | Medium / arithmetic risk |
| 11 | Traverse substitutions by contiguous row | Better locality and simpler indexing | Medium / moderate |
| 12 | Improve benchmark coverage and comparability | Reliable prioritization | Small–medium / low |

## 1. Fill Shoup and random buffers directly

**Evidence:** `crates/fhe-math/src/rq/mod.rs`, `compute_coefficients_shoup`, `random`, and `random_from_seed`; `crates/fhe-math/src/zq/mod.rs`, `shoup_vec` and `random_vec`.

Shoup computation allocates the destination matrix, then allocates a temporary vector for every row and copies it into that matrix. Random sampling follows the same temporary-vector-and-copy pattern. For `Poly<NttShoup>::random`, `Poly::zero` also allocates a zero Shoup matrix that is subsequently replaced by `compute_coefficients_shoup`.

Add internal `shoup_into(input, output)` and `random_into(output, rng)` primitives. Fill matrix rows directly and reuse an existing correctly sized Shoup matrix. Keep allocating public helpers as thin wrappers if useful to callers. `random_from_seed` can delegate to `random` after hashing and initializing its PRNG.

**Expected benefit:** remove L temporary allocations and O(LN) copying per operation, plus the redundant initial Shoup allocation. Actual runtime improvement depends on how much division or sampling dominates.

**Validate:** compare Shoup values against the existing scalar reference; preserve exact seeded output and RNG consumption order with fixed-seed regression vectors. Benchmark NTT-to-Shoup conversion and random generation separately, including allocation counts. Preserve secret-buffer cleanup when changing ownership.

## 2. Borrow contexts while traversing the level chain

**Evidence:** `crates/fhe-math/src/rq/context.rs`, `niterations_to` and `context_at_level`.

Both start with `Arc::new(self.clone())`. `Context::clone` copies its boxed operator and bit-reversal arrays; native NTT operators themselves contain boxed tables. `niterations_to` needs no owned root at all, and `context_at_level(i > 0)` immediately discards that copied root after moving to its child.

Traverse using borrowed `&Context` references. For a nonzero level, clone only the destination `Arc`. Returning level zero without copying would require an Arc-aware entry point, such as a method with `self: &Arc<Self>`; retain compatibility or introduce an additional method rather than silently changing the public API.

Also avoid explicit `as_ref()` comparisons of large derived contexts on frequent paths, such as `rq::Scaler::scale`, when an Arc equality comparison can first recognize shared identity. Preserve structural equality for independently constructed equivalent contexts; pointer equality alone would change behavior.

**Validate:** level-zero, intermediate, and unreachable contexts; separately allocated equivalent contexts. Measure allocations and latency at multiple chain lengths.

## 3. Use word-sized coprimality checks during RNS setup

**Evidence:** `crates/fhe-math/src/rns/mod.rs`, `RnsContext::new`.

The constructor runs arbitrary-precision extended GCD for every ordered pair of distinct u64 moduli, so each pair is checked twice. Only the GCD is needed for this validation.

Use a u64 Euclidean GCD and examine each unordered pair once. This keeps O(L²) validation but removes big-integer allocation, unnecessary Bézout coefficient computation, and duplicate work. Reuse the validated result when building a chain if construction is later reorganized.

Do not replace the modular inverse used for CRT precomputation with a prime-only shortcut: `RnsContext` supports coprime composite moduli too. Removing the second big-integer dependency is a separate investigation because it currently supplies inversion as well as GCD.

**Validate:** composite coprime moduli, repeated moduli, shared factors, and invalid moduli. Preserve useful error reporting; benchmark construction rather than steady-state arithmetic.

## 4. Reduce serialization and deserialization copies

**Evidence:** `crates/fhe-math/src/rq/convert.rs`, `From<&Poly<R>> for Rq`, `parse_proto`, and the `TryConvertFrom<&Rq>` implementations.

Serialization clones even a power-basis polynomial. For Shoup input it clones the auxiliary Shoup matrix although serialization only needs coefficients transformed into power basis. Packing also allocates a byte vector per modulus before collecting the combined payload.

Borrow power-basis coefficients directly. For transformed input, copy only the coefficient matrix needed for the inverse NTT. Precompute payload capacity and add a packing routine that appends into the final payload buffer. The existing wire format stores power-basis coefficients with a representation tag; preserve that format.

Deserialization already checks canonical residues, then routes them through public constructors that reduce the values again. An internal constructor accepting validated, standard-layout coefficients could avoid this second pass. Keep public raw constructors normalizing arbitrary residues, and keep wire parsing rejecting noncanonical residues. A private validated wrapper can make the boundary explicit without exposing an unchecked public API.

**Validate:** byte-for-byte compatibility for all three representations, malformed-input regressions, degree checks, and local timing-policy behavior. Benchmark serialization/deserialization separately from NTT cost and track peak temporary memory.

## 5. Consolidate representation-independent arithmetic and constructors

**Evidence:** `crates/fhe-math/src/rq/ops.rs`, duplicated PowerBasis/Ntt addition, subtraction, scalar multiplication, and negation; `crates/fhe-math/src/rq/convert.rs`, raw Vec/Array2 constructors.

The arithmetic implementations repeat row iteration, timing-policy propagation, and modulus dispatch. Constructors repeat dimensions, normalization, field initialization, and Shoup setup.

Extract small private row kernels and a common validated-storage initializer. Keep public trait implementations explicit and thin, or use a narrowly scoped macro where that is clearer. Prefer static dispatch so an abstraction does not introduce per-coefficient indirect calls.

Do not blanket-implement mutable arithmetic for every representation tag: changing NttShoup coefficients requires updating its cached Shoup values. Preserve the distinct lazy-coefficient preconditions and the rule that one secret operand forces constant-time execution.

**Expected benefit:** easier review and fewer places to miss an invariant; runtime should remain neutral unless common helpers also remove allocations.

**Validate:** existing arithmetic property tests plus mixed timing policies, lazy-input restrictions, layout normalization, and Shoup-cache correctness. Check representative release benchmarks for abstraction overhead.

## 6. Share immutable NTT tables across context levels

**Evidence:** `crates/fhe-math/src/rq/context.rs`, recursive `Context::new`; `crates/fhe-math/src/ntt/native.rs`, NTT table fields and constructor.

Each context recursively constructs the shorter modulus prefix. An L-modulus chain therefore constructs L(L+1)/2 modulus-specific NTT operators, although there are only L distinct `(modulus, degree)` pairs. Every level also owns an identical degree-sized bit-reversal table. Native table storage across the chain consequently grows as O(L²N).

Build immutable per-modulus operators once and share them across levels using Arc-backed storage or a shared table collection with a prefix length. Share bit-reversal data by degree. Keep level-specific CRT products and last-modulus inverses separate. This can reduce NTT table storage to O(LN), without implying that all context metadata becomes linear.

**Tradeoff:** additional indirection and a broader internal layout change. Prefer per-chain sharing before adding a global cache with eviction and synchronization concerns.

**Validate:** compare all levels against independently constructed contexts, round-trip transforms, and modulus switching. Measure cold construction, retained memory, and steady-state NTT throughput.

## 7. Store only the selected NTT backend

**Evidence:** `crates/fhe-math/src/ntt/tfhe.rs`, `NttOperator::new` and backend dispatch.

With the feature enabled, the wrapper constructs a native operator and also attempts to construct a TFHE plan. When the plan is present, the wrapper routes the exposed transforms to it, leaving the native tables retained principally for fallback structure and equality.

Consider an enum containing either a supported TFHE plan or the native operator. Retain a small mathematical identity `(modulus, size)` for equality. Construct native tables only when plan construction fails. Coordinate this with context-table sharing rather than duplicating ownership redesigns.

**Validate:** both supported-plan and fallback parameters; native/backend transform agreement, normalization, lazy-output range contracts, equality, and construction failure behavior. Establish timing suitability from the pinned backend implementation before changing which backend handles secret data. Measure setup and memory, not just transform speed.

## 8. Streamline dot-product preparation and reuse scratch

**Evidence:** `crates/fhe-math/src/rq/ops.rs`, `dot_product` and `fma`.

The function clones and traverses its iterators repeatedly for lengths, the first element, context validation, timing policy, and accumulation. Every invocation allocates a u128 accumulator, per-modulus counters and limits, and a result matrix.

First combine metadata checks into one preparation pass while retaining the general iterator API and exact errors. Consider a slice-oriented fast path for common callers. For repeated key-switch workloads, offer an internal workspace or `dot_product_into` variant that reuses buffers; reset all state on every invocation. Per-modulus accumulation limits depend on public context data and could be precomputed.

The raw-pointer reduction loop can also be tested against a safe row/chunk implementation. Keep the simpler version if generated code and benchmarks are equivalent.

**Tradeoffs:** collecting arbitrary iterators adds an allocation, so it is not automatically an improvement. A reusable workspace needs explicit secret-memory cleanup. Do not perform variable-time reductions before determining that every input permits them, or loosen overflow bounds to reduce the number of reductions.

**Validate:** empty/mismatched iterators, foreign contexts, mixed timing policies, and lengths immediately below/at/above accumulation thresholds. Compare against modular multiply-and-add; benchmark short and long products with and without scratch reuse.

## 9. Tile polynomial RNS scaling for locality

**Evidence:** `crates/fhe-math/src/rq/scaler.rs`, `Scaler::scale`; `crates/fhe-math/src/rns/scaler.rs`, `RnsScaler::scale`.

Polynomial storage is modulus-major, but scaling visits one coefficient column at a time. Its residues are separated by N words, and each coefficient invokes scalar scaling with repeated dimension checks and passes over precomputed data.

Prototype a batch kernel processing a small tile of adjacent coefficients, keeping rounding accumulators per coefficient while iterating through modulus rows. Hoist public shape checks to the batch boundary. Compare this with copying a tile into coefficient-major scratch; a full-matrix transpose may cost more than it saves.

Retain the existing common-prefix shortcut for identity scaling and the rule that only new destination rows require forward NTT. Reusable scratch for the source inverse transform is another candidate for repeated calls.

**Validate:** exact BigUint references, signed rounding ties, identity/general factors, partially shared bases, and PowerBasis/Ntt inputs. Benchmark complete polynomial scaling: the existing scalar RNS benchmark cannot establish gains from layout changes. Preserve branch-free secret-dependent rounding.

## 10. Fuse native NTT final passes where beneficial

**Evidence:** `crates/fhe-math/src/ntt/native.rs`, `forward_vt`, `forward_vt_lazy`, `backward`, and `backward_vt`.

Variable-time forward NTT currently runs the lazy transform and then scans the entire output to reduce it. Inverse NTT performs a separate final normalization scan. The constant-time forward path already integrates reduction into its last butterfly stage.

Try a canonical-output variant that reduces during the final variable-time forward stage while retaining the existing lazy variant. For inverse NTT, test applying normalization as each final-stage output is produced. Removing a memory pass does not guarantee a speedup: instruction scheduling and vectorization may become worse. Arithmetic fusion into twiddle constants would need a separate range proof.

**Validate:** round trips and reference polynomial multiplication across degrees and modulus widths, especially near the modulus limit. Check every intermediate range and inspect generated code for constant-time paths. Benchmark each direction independently with the native backend before comparing feature-enabled builds.

## 11. Make substitution access more contiguous

**Evidence:** `crates/fhe-math/src/rq/mod.rs`, `SubstitutionExponent` and `Poly::substitute`.

Power-basis substitution walks coefficient columns across all moduli, using ndarray slices for each column. NTT substitution accesses both source and destination through bit-reversal indices.

For power basis, move the modulus-row loop outside and operate on contiguous row slices. Precompute destination/sign mappings in the exponent object if it is reused enough to amortize their storage. For NTT, precompute a source index for each sequential destination index; then destination writes can be linear. Reuse the same mapping for Shoup data.

**Tradeoff:** additional O(N) mapping storage and setup. Power-basis even exponents can map multiple inputs to the same output, so preserve accumulation rather than assuming every substitution is a permutation.

**Validate:** odd and even exponents, repeated destinations, signs, exponent periodicity, and all representations. Retain foreign-context rejection. Benchmark one-shot versus repeated use of an exponent.

## 12. Make benchmarks support these decisions

**Evidence:** `crates/fhe-math/benches/rq.rs`, `rns.rs`, and `ntt.rs`.

The current suite covers useful arithmetic cases, but the main binary-operation macros use one modulus, dot products use length 256, and the RNS benchmark measures a single residue vector. In the dot-product comparison, naive variants repeatedly add into an existing output while the optimized version returns a fresh output. Several transform benchmarks intentionally include cloning in the timed operation.

Add context construction/traversal, Shoup conversion, sampling, serialization, substitution, complete polynomial scaling, and short dot-product benchmarks. Separate end-to-end allocation-inclusive measurements from kernels with preallocated scratch. Reset naive dot-product outputs so both variants compute the same operation; use batched setup where the aim is to exclude cloning/reset costs. Keep a separate allocation-inclusive comparison.

Use `std::hint::black_box` consistently for inputs and outputs and inspect suspiciously small results. Cover degrees 1024–8192 plus larger supported degrees, one and multiple moduli, several modulus widths, and both timing policies. Run native and `tfhe-ntt` configurations on the same machine. Record CPU, toolchain, feature set, allocation counts, and uncertainty; do not report a universal percentage from one parameter set.

## Implementation acceptance criteria

Keep the recent boundary checks, canonical-residue rules, layout normalization, and timing-policy propagation intact. Branching on public dimensions or modulus parameters is different from branching on secret coefficients. Allocation and ownership changes must preserve zeroization behavior; new scratch containing secrets needs an explicit cleanup policy.

For each implemented change, add targeted regression/property tests and benchmark the affected operation before and after on the same machine. Run `cargo test --workspace`, tests with `--all-features`, `cargo +nightly fmt --all`, and `cargo clippy --workspace --all-targets --all-features -- -D warnings`. Arithmetic/range changes additionally need release-mode tests. The remaining proposals are unimplemented and have no measured performance gains.


## Implementation results: items 1–4

Implemented direct Shoup/random row filling and Shoup-buffer reuse, borrowed context traversal, u64 coprimality validation, and serialization buffer reuse. Wire parsing constructs polynomials directly after validating shape and canonical residues. Transformed serialization scratch is zeroized when dropped. The existing `context_at_level(0)` borrowed-receiver API still clones the root; nonzero levels reuse existing child Arcs.

Added `crates/fhe-math/benches/allocation_paths.rs` and regression tests for RNG stream compatibility, Shoup storage reuse, context identity/equivalence, composite and invalid moduli, byte packing with existing prefixes, and legacy serialization payloads. Default/all-feature workspace tests, release math tests, nightly formatting, and all-target/all-feature Clippy with warnings denied passed.

The following are Criterion timing estimates from the same macOS arm64 host, using `rustc 1.99.0-nightly (d453bdd8f 2026-08-14)`, the native NTT backend, degree 2048, and three approximately 62-bit moduli. Runs used 20 samples, 100 ms warmup, and 300 ms measurement per case. These short runs establish local evidence, not general performance guarantees. Allocation counts were not instrumented; removed allocations were identified from the code paths.

| Operation | Before | After |
| --- | --- | --- |
| Seeded random Shoup polynomial | 115.31 µs | 111.10 µs |
| Clone NTT polynomial and convert to Shoup | 70.82 µs | 69.93 µs |
| Lookup level 2 | 3.814 µs | 4.002 ns |
| Distance to level 2 | 3.791 µs | 3.581 ns |
| Construct three-modulus RNS context | 4.092 µs | 2.898 µs |
| Serialize power basis | 127.59 µs | 33.54 µs |
| Serialize Shoup | 155.46 µs | 60.76 µs |
| Deserialize Shoup | 182.68 µs | 165.46 µs |

Criterion found clear improvements in random generation, child traversal, serialization, and Shoup deserialization in this run. Shoup conversion and RNS setup comparisons were classified within the noise threshold; avoid drawing firm speedup conclusions for those two from these samples. Serialization and conversion timings include allocation; Shoup conversion also includes cloning the input. Context timings refer to an existing chain and shared child, not construction or comparison against a separately allocated equivalent context.

To reproduce a comparison, install this benchmark on the baseline revision before modifying production code, then run both revisions with the same absolute result directory:

```sh
CRITERION_HOME=/tmp/fhe-allocation-bench cargo bench -p fhe-math --bench allocation_paths -- --save-baseline before
# Apply items 1–4, then:
CRITERION_HOME=/tmp/fhe-allocation-bench cargo bench -p fhe-math --bench allocation_paths -- --baseline before
```


## Implementation results: items 5–7

Shared arithmetic is generated only for PowerBasis and Ntt; NttShoup mutation is deliberately excluded. Raw constructors share shape validation, normalization, and storage initialization, with Shoup caches computed only for the matching representation. Wire parsing still rejects noncanonical data before using the common initializer.

Context levels share per-modulus NTT operators through Arcs and share one bit-reversal array. Level-specific RNS products and switching constants remain separate. The TFHE wrapper stores exactly one backend, constructing native tables only when a TFHE plan is unavailable. Size 8 exercises that fallback. Backend selection, forward-output ordering, inverse normalization, and lazy-output ranges remain unchanged. Native and TFHE roots can differ, so tests compare decoded ring products across backends rather than requiring equal transform vectors.

**0.2.0 API change:** `Context::context_at_level` now takes `self: &Arc<Self>`. Level zero returns an Arc clone of the existing root, eliminating the copy retained in items 1–4. All repository callers use Arc contexts. This supersedes the compatibility limitation recorded above; backward compatibility with 0.1.1 is not required.

Tests verify table identity across levels, independent CRT metadata, contexts that outlive their root, equivalent independently constructed contexts, arithmetic values and timing policies, lazy-input rejection, constructor boundaries, both backend paths, and ring-product equivalence. Default/all-feature workspace tests, release math tests with all features, formatting, and Clippy passed.

Extended the allocation benchmark with context setup, addition, and forward transformation. Compared against `4114b96` with the benchmark additions, on the same macOS arm64 host/toolchain and degree-2048, three-modulus parameters used above:

| Operation | Native before → after | TFHE feature before → after |
| --- | --- | --- |
| Construct context chain | 1.532 ms → 0.883 ms | 1.816 ms → 0.460 ms |
| Clone and forward transform | 30.29 µs → 30.68 µs | 29.92 µs → 29.58 µs |
| Allocate and add polynomials | 1.698 µs → 1.718 µs | 1.699 µs → 1.730 µs |

Setup improved clearly in these short runs. Transform changes and native addition were within noise/no-change classifications. TFHE addition initially showed about a 2% slowdown; a repeat measured 1.716 µs and Criterion classified the change within its noise threshold. These samples do not establish arithmetic speedups. Table sharing is also asserted directly by tests; retained-memory bytes were not profiled. NTT table storage across an L-modulus chain becomes O(LN), while level metadata and Arc lists still have their own costs.


## Implementation results: item 8

Added `rq::DotProductWorkspace`, with fixed-context scratch, precomputed per-modulus accumulation limits, and reusable counters. `workspace.dot_product(...)` allocates a result while reusing scratch; `workspace.dot_product_into(..., &mut out)` reuses both. The original `rq::dot_product` remains the allocation-inclusive convenience function. The BFV long-dot-product fallback now shares one workspace across ciphertext components.

Validation counts and inspects each input once, followed by one arithmetic pass, without collecting iterator contents. It preserves empty/length/context error precedence and checks all timing permissions before any reduction. Cloned iterators must produce the same operands. Equivalent separately allocated contexts remain accepted. Output and workspace contexts are checked before mutation; validation errors leave the output unchanged.

The accumulator is initialized to zero and guarded throughout computation. The guard zeroizes it after success or panic unwinding; workspace drop therefore releases already-cleared storage. Buffers never grow. This replaces the previous unzeroized temporary accumulator. Cleanup is not repeated through `Zeroizing<Vec<_>>`, whose element-by-element and full-capacity wiping added avoidable overhead in an intermediate implementation.

Periodic reduction uses safe slice iteration and retains the original per-modulus bounds and short-product fast path. Lazy NTT inputs are now rejected with `LazyDotProductOperand`: their larger coefficient range invalidates the canonical-input accumulation bound. Convert them to canonical NTT coefficients before using dot products.

Usage with an existing `Arc<Context>` and slices of NTT polynomials:

```rust
let mut workspace = fhe_math::rq::DotProductWorkspace::new(&ctx);
let mut output = fhe_math::rq::Poly::<fhe_math::rq::Ntt>::zero(&ctx);
workspace.dot_product_into(left.iter(), right.iter(), &mut output)?;
// Reuse workspace and output for the next set of inputs.
```

Tests cover 15/16/17 and repeated reduction boundaries with maximum canonical residues, mixed modulus widths, changing lengths and timing policies, equivalent/foreign contexts, lazy inputs, unchanged output on errors, unchanged buffer addresses, exact iterator visit counts, and cleanup/recovery after an iterator panic. A BFV regression forces the fallback with 62-bit moduli and lengths 17 and 33.

Native-backend Criterion samples on the same macOS arm64 host/toolchain, degree 2048 and three 62-bit moduli, comparing against `2090337` with the benchmark additions:

| Terms | Previous fresh call | New fresh call | New workspace + output reuse |
| --- | --- | --- | --- |
| 4 | 22.55 µs | 23.85 µs | 21.27 µs |
| 16 | 55.57 µs | 58.91 µs | 53.93 µs |
| 256 | 890.60 µs | 894.94 µs | 895.29 µs |

These short runs show modest short-product gains from reuse, but a roughly 6–8% regression for fresh short calls with the new cleanup guarantee. Long-product changes were within noise; no long-product speedup is claimed. Buffer identity is checked by tests, but total allocator activity was not instrumented. Benchmark filters are `allocation_paths/dot_product/` and `allocation_paths/dot_product_reuse/`; the latter has no pre-change API equivalent. Retain an absolute `CRITERION_HOME` when comparing revisions.

Validation passed: default and all-feature workspace tests, release math tests with all features, nightly formatting, and all-target/all-feature Clippy with warnings denied.


### Item 8 follow-up: BFV fast path and PIR examples

The initial BFV integration only covered the long-product fallback. For the million-entry, 288-byte workloads, SealPIR uses 170-term columns with a 2^56 fallback threshold, and MulPIR uses 119-term columns with a 262,144 threshold. Neither reaches that fallback.

Added `bfv::DotProductScalarWorkspace::new(&params, level)` and its `dot_product_scalar` method. It preserves the fused ciphertext/plaintext fast loop and retains its accumulator between calls; the long path lazily retains a polynomial workspace too. Scratch is cleared after each computation and on unwinding. Parameter/level validation, part-count checks, and timing-permission checks remain in place. Results have independent coefficient buffers; a changed ciphertext part count resizes the fast accumulator on the next fast-path call.

The fast path now converts an `ArrayView2<u128>` directly into canonical `Poly<Ntt>` coefficients using a checked conversion in `fhe-math`. This removes the second reduction previously performed by the u64 constructor. The new conversion validates shape, accepts strided input views, and reduces arbitrary u128 values with the selected timing policy.

SealPIR and MulPIR each create one BFV workspace outside the timed response loop and reuse it across database columns and response repetitions. SealPIR also reuses it for its final folded dot products. The free `bfv::dot_product_scalar` function is a convenience wrapper over the same implementation with temporary workspace.

Both exact commands completed successfully and verified the requested database entry:

```sh
cargo run --example sealpir --release -- --database-size 1000000 --element-size 288
cargo run --example mulpir --release -- --database-size 1000000 --element-size 288
```

| Server response (five-response average per invocation) | Before BFV fast-path integration | After |
| --- | --- | --- |
| SealPIR | 576.2 ms | 572.2 ms |
| MulPIR | 954.0 ms | 939.9 ms |

These are single sequential before/after invocations on the same host, not a statistical performance guarantee. Observed changes are modest (approximately 0.7% and 1.5%). The workspace is now used by the actual hot paths, but bulk polynomial arithmetic and database reads remain, and query expansion alone takes about 264–273 ms. This integration does not imply a large end-to-end speedup.

New tests cover repeated fast/fallback transitions, 2/3-part ciphertexts, buffer identity, public/secret policy changes, invalid levels and parameters, and full-width/strided u128 conversion. Default and all-feature workspace tests, release math tests, nightly formatting, and all-target/all-feature Clippy passed.


## Implementation results: item 10

Native `forward_vt` now reduces each output in its final butterfly stage. A
compile-time canonical/lazy choice shares the earlier stages, and a separate
adjacent-pair final loop lets the compiler optimize that stage independently.
`forward_vt_lazy` retains its [0, 4p) output contract. No buffers, twiddle constants,
public APIs, or timing-policy dispatch were added or changed.

The existing butterfly bounds still apply: inputs are below 4p, its reduced left
input and Shoup product are below 2p, and both outputs are below 4p. Since p is
less than 2^62, the sums fit in u64. The same two conditional subtractions as the
old standalone pass produce values below p. Constant-time forward and both inverse
implementations remain unchanged. Final-stage inverse normalization was prototyped
and tested, but removed after several native cases regressed by approximately
3–6%; fewer passes alone did not justify retaining it.

Native-specific regressions run even with `tfhe-ntt` enabled. They cover degrees
8–16384; 17-, 31-, and near-limit 62-bit moduli; zero, maximum and mixed lazy
residues; canonical output and lazy ranges; round trips with inverse inputs below
2p; and independent negacyclic schoolbook products (dense at small degrees,
sparse at large degrees). Debug assertions check intermediate butterfly ranges.
Release AArch64 assembly inspection confirmed that canonical reduction occurs
before final stores, with no later reduction scan. The constant-time forward and
backward instruction streams match the baseline apart from labels and panic
metadata references; this is a local compiler check, not a cross-platform timing
proof.

The NTT benchmark now uses `black_box`, covers degrees 1024–16384, and measures
each direction independently using a reused in-place buffer. The small benchmark
modulus is now 65537, which supports all these degrees. The following native
forward-VT means are from a repeated before/after comparison on an Apple M2,
macOS arm64, `rustc 1.99.0-nightly (d453bdd8f 2026-08-14)`, no optional features,
50 samples, 200 ms warmup and 600 ms measurement per case. Both revisions used
the same updated benchmark. Times exclude allocation and setup; allocation counts
were not instrumented, and the kernels remain allocation-free.

| Degree | Modulus bits | Before | After (95% confidence interval) |
| --- | --- | --- | --- |
| 1024 | 17 | 3.167 µs | 3.035 µs (3.025–3.045) |
| 1024 | 62 | 3.193 µs | 3.023 µs (3.014–3.034) |
| 4096 | 17 | 14.931 µs | 13.843 µs (13.810–13.876) |
| 4096 | 62 | 14.439 µs | 13.899 µs (13.861–13.940) |
| 8192 | 17 | 31.545 µs | 29.986 µs (29.915–30.058) |
| 8192 | 62 | 31.400 µs | 30.008 µs (29.919–30.105) |
| 16384 | 17 | 69.414 µs | 70.119 µs (69.938–70.302) |
| 16384 | 62 | 68.998 µs | 70.093 µs (69.892–70.308) |

Criterion classified degrees 1024–8192 as improvements in this repeat, the
16384/17 change within noise, and 16384/62 as an approximately 1.6% regression.
Earlier short runs were noisier and showed smaller gains. Retaining the forward
fusion trades modest gains at smaller degrees for a small measured large-degree
regression on this host; these measurements do not establish a universal speedup
or a benefit for the separate TFHE backend.

To reproduce, install the updated benchmark on the baseline before changing
production code, then run each revision with the same absolute result directory:

```sh
CRITERION_HOME=/tmp/fhe-ntt-item10 cargo bench -p fhe-math --bench ntt -- forward_vt --warm-up-time 0.2 --measurement-time 0.6 --save-baseline before
# Apply the forward fusion, then:
CRITERION_HOME=/tmp/fhe-ntt-item10 cargo bench -p fhe-math --bench ntt -- forward_vt --warm-up-time 0.2 --measurement-time 0.6 --baseline before
```

Validation passed: default and all-feature workspace tests, native and all-feature
release math tests, nightly formatting, and default/all-feature all-target Clippy
with warnings denied.

### Item 10: million-entry PIR check

Ran both commands on the same host, first with the original native NTT and then
with the retained forward fusion:

```sh
cargo run --example mulpir --release -- --database-size 1000000 --element-size 288
cargo run --example sealpir --release -- --database-size 1000000 --element-size 288
```

| Server response | Before | After | Observed change |
| --- | --- | --- | --- |
| MulPIR | 972.4 ms | 958.0 ms | 1.5% faster |
| SealPIR | 590.9 ms | 575.2 ms | 2.7% faster |

Both examples passed their database-entry equality assertions before and after.
Each timing is the example's five-response average from one invocation; random
inputs and host variability prevent treating these small changes as a statistical
performance guarantee. No example or BFV code was changed for item 10.


## MulPIR: database layout and preprocessing

Measured against `f34ec35` on macOS arm64 with
`rustc 1.99.0-nightly (d453bdd8f 2026-08-14)`, the default native backend,
and the release profile:

```sh
cargo run --example mulpir --release -- --database-size 1000000 --element-size 288
```

Temporary stage timers identified ciphertext multiplication as the largest
response cost: approximately 445 ms, compared with 273 ms for expansion and
232 ms for plaintext/ciphertext dot products. The timers were removed after
measurement.

The database packs into 14,085 plaintexts. A square 119 × 119 layout requires
119 ciphertext multiplications in the second response stage. The new layout
uses 174 × 81, reducing this to 81 multiplications. Both layouts require eight
expansion rounds and 255 Galois key switches: the implementation expands each
round's entire lower half, including the last round. The layout search minimizes
columns within the square layout's expansion budget. This preserves the BFV
parameters, encoding, and serialized query/response sizes. SealPIR's layout is
evaluated separately below.

Preprocessing now builds populated plaintexts directly and appends only the
padding, avoiding allocation and initialization of thousands of zero plaintexts
that were immediately replaced. This shared improvement also applies to SealPIR.

For the final comparison, the baseline and optimized release executables were
copied aside and run sequentially in four alternating pairs (before/after,
after/before, repeated). Each invocation computed five server responses and
verified the decrypted database entry. No builds or tests ran concurrently with
the measurements. The table reports medians across four invocations; process
wall time includes setup, preprocessing, all five responses, verification, and
cleanup, but excludes Cargo/build overhead.

| Metric | Before | After | Time reduction |
| --- | --- | --- | --- |
| Server response | 953.50 ms | 813.15 ms | 14.7% |
| Database preprocessing | 2102.65 ms | 1936.50 ms | 7.9% |
| Entire process | 7.136 s | 6.272 s | 12.1% |

The four response averages were 951.0, 956.0, 953.1, and 953.9 ms before,
and 813.0, 811.1, 814.1, and 813.3 ms after. Queries and keys are randomized;
these are local measurements rather than a cross-platform guarantee. Other
database sizes benefit according to unused capacity in their expansion round.

Regression tests check minimal column counts without extra expansion rounds
across small layouts and larger boundaries, including this workload, and decode
every database byte and padding entry for both layouts. Run the example tests
explicitly with `cargo test --example mulpir --example sealpir` (also tested
with `--release`).


## SealPIR: fewer database columns

SealPIR also selects `DatabaseLayout::FewerColumns`. For the same million-entry,
288-byte workload, its degree-4096 parameters pack 35 elements per plaintext,
requiring 28,572 plaintexts. The layout changes from 170 × 169 to 447 × 64.
Both use nine expansion rounds and 511 Galois key switches.

SealPIR does not perform MulPIR's ciphertext/ciphertext multiplications. Fewer
columns instead reduce intermediate modulus switches, ciphertext-to-plaintext
folding, and the length of the second-stage dot products. Temporary stage
measurements showed folding falling from approximately 30–32 ms to 12 ms,
and the warmed-up first dot-product/modulus-switch stage from 262–266 ms to
246–248 ms. Expansion rose by approximately 4 ms because more upper-half
outputs must be materialized; the key-switch count remains the same.

The following release comparison isolates the layout change: both versions
already include direct database construction from the MulPIR optimization.
The host, compiler, backend, and four alternating before/after pairs match the
method above. Each invocation averages five server responses and checks its
randomly selected database entry against the decrypted answer. Timings exclude
compilation, and no tests or builds ran concurrently.

```sh
cargo run --example sealpir --release -- --database-size 1000000 --element-size 288
```

| Metric (median of four invocations) | Square | Fewer columns | Time reduction |
| --- | --- | --- | --- |
| Server response | 577.45 ms | 539.75 ms | 6.5% |
| Database preprocessing | 1847.35 ms | 1846.65 ms | 0.04% |
| Entire process | 4.982 s | 4.813 s | 3.4% |

The response averages were 578.2, 577.1, 577.7, and 577.2 ms for the square
layout, and 538.8, 539.3, 549.1, and 540.2 ms for fewer columns. All eight
lookups passed. Evaluation keys remained 981.64 KiB, queries 36.05 KiB, and
responses 144.11 KiB. This is a smaller local gain than MulPIR's, consistent
with SealPIR's cheaper per-column work. The layout regression now includes
SealPIR's exact 28,572-plaintext case and verifies the expected 447 × 64 shape.


## RNS scaling: omit zero fractional corrections

Implemented after the PIR commit `0b1b2aa`, in
`crates/fhe-math/src/rns/scaler.rs`. The RNS scaler now stores only nonzero
fractional correction terms, together with their source indices. It also skips
multiplication by a zero gamma correction. These choices depend on the contexts
and scaling factor, not the input residues. Signed rounding, fixed-point
precision, modular reductions, and the public API remain the same.

For BFV downscaling from Q*P by t/Q, gamma = t*P is integral. The Garner
projections for moduli belonging to P also contain Q as a factor, so their
fractional corrections vanish. In the measured five-modulus to two-modulus case,
this removes three of five residue/correction products and the gamma product
for every coefficient. Packing the retained terms together avoids traversing
separate coefficient and index arrays. A specialization for identity scaling
keeps basis conversion free of the fractional-correction branches and scratch.

Added scalar and full-polynomial BFV scaling cases to the RNS benchmark. On the
same macOS arm64 host and nightly compiler as the PIR runs, using release mode
and the native backend, the final Criterion repeat produced the following
**means** (50 samples; BFV cases use 200 ms warmup and at least 600 ms measurement;
the existing generic cases use 3 s warmup and 5 s measurement):

| Operation | Before | After | Time reduction |
| --- | --- | --- | --- |
| Generic RNS scaling, 3 to 4 moduli | 48.20 ns | 47.04 ns | 2.4% |
| Identity basis conversion, 3 to 4 moduli | 30.06 ns | 29.93 ns | 0.4% |
| BFV RNS downscaling, 5 to 2 moduli | 47.90 ns | 34.05 ns | 28.9% |
| Degree-8192 NTT downscaling, constant-time transforms | 675.76 µs | 561.61 µs | 16.9% |
| Degree-8192 NTT downscaling, public data | 638.20 µs | 532.41 µs | 16.6% |

The final 95% confidence intervals for the BFV after-means were 33.93–34.17 ns,
559.62–563.76 µs, and 530.87–533.99 µs respectively. Earlier repeats showed
approximately 22–27% scalar and 14–16% polynomial improvements; performance
varies with host conditions. The initial implementation regressed identity
conversion by about 3%; the retained specialization removed that regression.

To reproduce, add the benchmark to the baseline before changing the scaler:

```sh
cargo bench -p fhe-math --bench rns -- --save-baseline before-rns-sparse
# Apply the scaler optimization, then:
cargo bench -p fhe-math --bench rns -- --baseline before-rns-sparse
```

Four alternating release MulPIR comparisons with the already-optimized layout,
1,000,000 entries, and 288-byte elements reduced median server-response time
from 818.05 ms to 790.75 ms (**3.3%**). The response averages were 817.6, 818.6,
815.7, and 818.5 ms before, and 789.7, 791.8, 792.9, and 787.1 ms after. All
lookups decrypted correctly. Whole-process median time fell from 6.294 s to
6.162 s (2.1%). No tests or builds ran concurrently with these measurements.

Tests exhaust small contexts, including even products, rounding ties, zero and
integer factors, identity scaling, sparse and dense corrections, reordered
source moduli, and strided input/output views. Large-context tests check exact
BigUint results on deterministic random samples and compare the sparse schedule
with the original dense schedule at centering and rounding boundaries.

### Existing fixed-point precision limitation

The new boundary tests exposed a pre-existing approximation error, reproduced
with the original scaler from `0b1b2aa`. With source moduli
`[562949954093057, 4611686018326724609, 4611686018309947393,
4611686018282684417, 4611686018257518593]`, destination moduli equal to the first
two, Q equal to their product, scaling factor 1/Q, and
`x = 1298074216154311220707323998969855 = (Q - 1)/2 - 1`, exact rounding gives
residues `[0, 0]`; both the original and optimized scalers give `[1, 1]`.
Discrepancies also occur extremely close to the source centering boundary.

The scaler's fixed-point approximations cannot distinguish every such boundary
at arbitrarily large modulus products. The optimization preserves the existing
results there; it does not fix this precision limitation. The type documentation
now states that limitation explicitly. Exact handling of these cases needs a
separate correctness change, with care to preserve constant-time processing.


## Focused core BFV benchmarks for sparse RNS scaling

Added `crates/fhe/benches/bfv_core.rs` to measure six operations without running
the full BFV benchmark suite: secret-key encryption, decryption, multiplication,
squaring, multiplication followed by relinearization, and relinearization alone.
The parameter sets are the library defaults at degree 4096 (three moduli,
109 total modulus bits) and degree 8192 (five moduli, 218 total modulus bits),
both at level zero with a 20-bit plaintext modulus.

```sh
cargo bench -p fhe --bench bfv_core
```

The harness uses seeded, fixed ciphertext inputs. Before timing, it decrypts
and checks fresh ciphertexts, products, squares, and relinearized products.
Relinearization input cloning is outside the timed region; the other operations
include their normal result allocation and destruction. Key/parameter setup is
outside every measurement, and ciphertexts are not repeatedly mutated across
iterations to avoid accumulating noise.

Both executables used the same harness and compiler configuration. The baseline
used the scaler from `0b1b2aa`; the optimized executable used the uncommitted
sparse-correction implementation. No other production source changed between
them. Four optimized-profile invocations ran sequentially in before/after,
after/before order, with no concurrent builds or tests. Each case used 30 flat
samples, 100 ms warmup, and a 600 ms measurement target. The table averages the
two invocation means for each version; all times are milliseconds.

| Operation | Degree | Before | After | Time reduction |
| --- | --- | --- | --- | --- |
| encrypt_sk | 4096 | 0.5623 ms | 0.5624 ms | -0.0% |
| encrypt_sk | 8192 | 1.8066 ms | 1.8223 ms | -0.9% |
| decrypt | 4096 | 0.4896 ms | 0.4746 ms | 3.1% |
| decrypt | 8192 | 1.4206 ms | 1.3906 ms | 2.1% |
| multiply | 4096 | 2.2222 ms | 2.0713 ms | 6.8% |
| multiply | 8192 | 8.0889 ms | 7.5647 ms | 6.5% |
| square | 4096 | 1.8225 ms | 1.6708 ms | 8.3% |
| square | 8192 | 6.6134 ms | 6.2069 ms | 6.1% |
| multiply_relinearize | 4096 | 2.5364 ms | 2.3719 ms | 6.5% |
| multiply_relinearize | 8192 | 9.6386 ms | 9.2327 ms | 4.2% |
| relinearize | 4096 | 0.3005 ms | 0.2996 ms | 0.3% |
| relinearize | 8192 | 1.5229 ms | 1.5264 ms | -0.2% |

Multiplication improves by 6.5–6.8%, and squaring by 6.1–8.3%. Their RNS scaling
work benefits, while NTTs, pointwise arithmetic, and allocation costs remain.
Adding relinearization dilutes the gain to 4.2–6.5%; relinearization itself does
not use this scaler in the tested same-level configuration. Decryption benefits
modestly (2.1–3.1%) from its final scaling step. Encryption and standalone
relinearization changed by less than 1%, which is too small to draw a useful
conclusion from these short runs.

The degree-8192 square and multiply-plus-relinearize measurements varied more:
optimized invocation means were 6.311/6.103 ms and 9.321/9.144 ms, respectively.
The table reports both runs rather than selecting the faster sample. These
results support a moderate improvement in multiplication-heavy BFV work, not
a blanket 17–29% BFV speedup inferred from the isolated scaler benchmarks.
All four runs passed the fixture assertions. The focused benchmark also passed
Clippy with warnings denied and nightly formatting checks.
