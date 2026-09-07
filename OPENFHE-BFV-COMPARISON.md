# BFV improvements from the OpenFHE comparison

Compared fhe.rs at `a5b536d` with the local OpenFHE checkout at
`92248109b288b670cd05e2fb67c7c7c3fff69dce`. This is a source comparison;
OpenFHE was not built or benchmarked. All timing comparisons below are between
Rust implementations on this machine, not performance claims against OpenFHE.

The existing library already implements HPS RNS multiplication, custom
asymmetric multiplication scalers, delayed rounding of ciphertext product
sums, and NTT-aware Galois key switching. The changes concentrate on three
remaining opportunities:

| Opportunity | OpenFHE source | Implementation here |
|---|---|---|
| Symmetric squaring | `LeveledSHEBase::EvalSquareCore`, [source](/Users/tancrede/git/openfhe-development/src/pke/lib/schemebase/base-leveledshe.cpp:666) | `Ciphertext::square` and `Multiplicator::square`; reuse extensions and compute each cross product once |
| Karatsuba multiplication | `USE_KARATSUBA` branch in `LeveledSHEBFVRNS::EvalMult`, [source](/Users/tancrede/git/openfhe-development/src/pke/lib/scheme/bfvrns/bfvrns-leveledshe.cpp:288) | Shared three-product kernel for ordinary/custom multiplication and product accumulation |
| Reusable operand preparation | `EvalFastRotationPrecompute` / `EvalFastRotation`, [rotation source](/Users/tancrede/git/openfhe-development/src/pke/lib/schemebase/base-leveledshe.cpp:455); `KeySwitchBV::EvalKeySwitchPrecomputeCore`, [decomposition source](/Users/tancrede/git/openfhe-development/src/pke/lib/keyswitch/keyswitch-bv.cpp:304) | Apply the separation of preparation from evaluation to BFV **multiplication**: cache one operand's basis extension and Shoup quotients |

The third item is an adaptation of the precomputation pattern, not a port of
OpenFHE's fast-rotation algorithm. It uses the existing fhe-math Shoup
representation and introduces no new key format.

For two-part inputs, Karatsuba computes
`c0 = a0*b0`, `c2 = a1*b1`, and
`c1 = (a0+a1)*(b0+b1)-c0-c2`. These additions occur in the extended basis,
after centered lifting. Adding before lifting would change BFV rounding.
Squaring computes `a0*a0`, `2*a0*a1`, and `a1*a1`. For larger ciphertexts,
squaring evaluates only the upper triangle of the tensor and doubles its
off-diagonal terms. General unequal-size products retain full convolution.
Both algorithms preserve the coefficients and rounding of the original
calculation; they do not consume another level or change its noise bound.

Previously the ciphertext operator detected equal inputs by comparing their
contents, and still computed both cross terms. `Multiplicator::multiply`
always extended both operands independently. Operators now recognize equal
references by identity; callers can explicitly use `square()` for cloned
values. Asymmetric custom strategies preserve separate left and right
scaling even when squaring.

The prepared API is:

```rust
let strategy = Multiplicator::default(&relinearization_key)?;
let prepared = strategy.prepare_lhs(&query)?;
let first = prepared.multiply(&encrypted_row_1)?;
let second = prepared.multiply(&encrypted_row_2)?;
```

Preparation stores the two extended parts and their sum in `NttShoup` form.
Each subsequent product skips the left basis conversion and performs three
multiplications with cached quotients. Its lifetime ties it to the strategy,
preventing accidental changes to the level, scaling, or relinearization
configuration. It owns a snapshot of the ciphertext: modifying the source
after preparation has no effect. Cached storage is zeroized on drop.

The storage tradeoff is approximately `6*N*L` u64 words, excluding shared
context tables and allocator overhead. At the benchmark parameters this is
960 KiB for degree 4096 and 3.375 MiB for degree 8192. Preparation includes
building the Shoup tables, so one-off multiplication should use the ordinary
path. Each prepared product still rounds independently. The existing
`CiphertextProductAccumulator` serves workloads that need one rounding for
a sum of products.

`Multiplicator::without_relinearization(&parameters, level)` reuses the
parameters' precomputed basis without requiring key generation.
`Ciphertext::try_mul` exposes validation errors instead of operator panics.
The checked ciphertext methods support unrelinearized inputs; multiplication
strategies and prepared operands require two parts. Checked methods reject
the empty zero sentinel, while operators retain its existing behavior.
`enable_relinearization` now also checks the key's BFV parameter instance,
preventing acceptance of a key merely because its polynomial moduli match.

Every output retains restricted timing when any input part is restricted.
The prepared representation propagates that restriction to all cached parts.
Temporary storage clears restricted polynomials on drop; explicitly public
ciphertext intermediates skip that additional memory pass.
No coefficient equality checks or value-dependent squaring shortcuts are
introduced. These changes do not alter parameter selection or ciphertext
serialization.

Other substantial OpenFHE algorithms remain candidates for separate work:

| Candidate | Evidence and reason to defer |
|---|---|
| HPS P-over-Q / leveled P-over-Q | BFV `EvalMult` branches [at lines 196 and 223](/Users/tancrede/git/openfhe-development/src/pke/lib/scheme/bfvrns/bfvrns-leveledshe.cpp:196) use specialized CRT expansion and depth-driven level dropping. Rust already exposes the asymmetric scaling strategy, but matching these fast kernels requires new CRT precomputations, rounding tests, and a noise/depth policy. |
| Hybrid key switching | [Hybrid precomputation](/Users/tancrede/git/openfhe-development/src/pke/lib/keyswitch/keyswitch-hybrid.cpp:412) partitions Q and introduces an auxiliary P basis. This affects key generation, decomposition, mod-down, parameter selection, and serialization together. |
| Hoisted rotations | OpenFHE's fast path switches before applying the automorphism. Current Rust Galois keys switch after substitution, and its unsigned residue lifts do not commute with sign-changing automorphisms. Reusing its present digits without changing the key convention or handling the lift correction would be incorrect. |

Validation includes full-convolution references for two through five parts,
the existing exact BigInt negacyclic product-sum tests, decrypted SIMD and
polynomial results, every level of a three-prime chain, relinearization with
larger key bases, modulus switching, asymmetric scalers, all input timing
permission combinations, snapshot reuse, and invalid inputs. The default
workspace tests and focused release tests with all features pass. Nightly
formatting, Clippy across all targets with warnings denied, and rustdoc with
warnings denied pass.

Reproduce the focused measurements with:

```sh
cargo bench -p fhe --bench bfv_multiplication
cargo bench -p fhe --bench bfv_core -- 'bfv_core/(multiply|square|product_sum_8)'
```

The first benchmark validates its fixtures before timing, compares distinct
encrypted right operands in batches of eight, and measures both preparation
alone and batches including preparation/destruction. Its square control uses
a separate clone to force the generic multiplication path. The measurements
use the default NTT backend, release optimization, 30 flat Criterion samples,
and Rust `1.99.0-nightly (d453bdd8f 2026-08-14)` on aarch64 macOS. They do not
include relinearization or modulus switching unless the benchmark name says
so. Feature-enabled builds are correctness checks only.

Final Criterion point estimates, in milliseconds:

| Work | Degree 4096, log Q = 109 | Degree 8192, log Q = 218 |
|---|---:|---:|
| Ordinary strategy multiplication | 2.015 | 7.349 |
| Prepared multiplication, setup excluded | 1.560 | 5.772 |
| Preparation alone, including destruction | 1.363 | 4.576 |
| Eight ordinary products | 16.372 | 59.562 |
| Eight prepared products, including setup and destruction | 13.980 | 50.849 |

Preparation saves approximately 22.6% / 21.5% per subsequent product and
14.6% for the measured batch of eight at either degree. These measurements
suggest amortization after about four products; the exact crossover depends
on the workload and allocator. Stored table memory and preparation overhead
make it inappropriate to prepare every operand for just one multiplication.

The unchanged `bfv_core` benchmark also permits a before/after comparison
against the saved pre-change baseline:

| Work | Degree 4096, before → after | Degree 8192, before → after |
|---|---:|---:|
| Ordinary multiplication | 2.044 → 2.008 | 7.473 → 7.386 |
| Ciphertext operator square | 1.654 → 1.584 | 6.019 → 5.836 |
| Fused sum of eight products | 9.457 → 9.412 | 33.896 → 33.517 |

Ordinary multiplication improves about 1–2%; basis conversion dominates, so
removing one tensor product does not produce a 25% end-to-end improvement.
The operator square improves about 3–4%. These comparisons include the
validation and storage changes, not an isolated measurement of Karatsuba.
The small fused-accumulator differences are around the benchmark's noise
threshold and should not be treated as a demonstrated speedup.

`Multiplicator::square` additionally avoids the repeated basis extension that
its old generic multiplication performed. The current generic-square control
versus explicit square measured 2.017 → 1.572 ms at degree 4096. The degree-8192
square sample was noisy (95% interval 5.788–6.612 ms; generic control
7.334–7.358 ms), so a precise speedup percentage for that case is unwarranted.

Raw logs are in `/tmp/fhe-openfhe-core-final.log` and
`/tmp/fhe-openfhe-prepared-final.log`. The original core baseline is saved as
`before-openfhe` in Criterion's output. These are local measurements and may
vary on other architectures or with the optional NTT backend.

Validation commands:

```sh
cargo test
cargo +nightly fmt --all
cargo clippy --all-targets -- -D warnings
cargo test -p fhe --release --all-features bfv::ops
RUSTDOCFLAGS='-D warnings' cargo doc -p fhe --no-deps
```
