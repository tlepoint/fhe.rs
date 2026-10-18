# Verifying `fhe-math` with Verus

The verified kernels in
[`crates/fhe-math/src/zq/verified.rs`](../crates/fhe-math/src/zq/verified.rs)
form a compositional proof: `const_time_cond_select`, then `reduce1`, then
`add_mod`. Each caller uses the preceding helper's proved contract.
`Modulus::add`, `Modulus::reduce1`, and the other conditional-selection callers
use these same functions in ordinary builds. Verus checks the executable
implementations and their contracts; the proof is not a separate copy of the
algorithms. Normal Cargo builds erase the proof code through `vstd`.

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

The expected result is `5 verified, 0 errors`: three functions and two local
proof obligations. `--no-cheating` rejects assumptions, admits, and unchecked
external bodies in our proof source. Verus, its bundled standard-library
specifications, Z3, Rust, and the compiler remain trusted infrastructure.

This checks the scalar module directly, without asking Verus to analyze the
crate's unrelated dependencies or algorithms. It needs neither Cargo dependency
compilation nor `protoc`. It is not a `cargo verus verify -p fhe-math` integration
or a verification of the whole crate. The pinned `vstd` dependency is used by
ordinary Rust builds; update it together with the verifier release and CI job.
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

All arithmetic in contracts is interpreted over mathematical integers. The
proofs cover the full stated input domains, including even and composite
moduli, and preserve `const fn` use.

The reduction and addition proofs assume valid inputs. `Modulus::new` checks
the modulus; the `Modulus::add` and `Modulus::reduce1` wrappers retain their
debug assertions. Verus does not verify these wrappers or their callers, or
insert runtime checks. Although other operations now reuse the verified
selector and reducer, their contracts have not been proved. These proofs do
not establish constant-time machine code, NTT/RNS correctness, or FHE security.

## Tests and CI

Ordinary Rust tests exercise both selection outcomes on boundary and
alternating-bit words, exhaust reduction for moduli 1 through 64 and addition
for moduli 2 through 64, and check reduction/addition boundaries at the maximum
supported widths. They also check all three kernels remain usable in `const`
expressions. Existing randomized arithmetic tests exercise their production
callers:

```sh
cargo test -p fhe-math zq::tests
```

The [Verus workflow](../.github/workflows/verus.yml) installs the pinned
checksum-checked Linux release and runs the proof on pull requests and pushes.
The ordinary test, formatting, Clippy, and MSRV jobs continue to check Rust
builds independently. When changing the kernel, rerun both the proof and these
tests. To sanity-check the contracts, removing the mask's negation must fail
the selector proof; changing the reduction comparison `x < p` to `x <= p`
must fail the reduction proof; dropping `b` from `reduce1(a + b, p)` must fail
the addition proof.
