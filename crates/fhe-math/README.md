# fhe-math [![crate version](https://img.shields.io/crates/v/fhe-math.svg)](https://crates.io/crates/fhe-math) [![documentation](https://docs.rs/fhe-math/badge.svg)](https://docs.rs/fhe-math)

Core mathematical primitives for the [`fhe.rs`](https://github.com/tlepoint/fhe.rs) ecosystem.

This crate exposes building blocks such as number theoretic transforms (NTT), residue number system (RNS) arithmetic, and ring arithmetic over `Z_q` that are used by higher level crates like [`fhe`](https://crates.io/crates/fhe).

## Features

* `ntt`, `rns`, `rq`, and `zq` modules for modular arithmetic over large rings.
* Optional `tfhe-ntt` features to enable hardware accelerated NTTs via the [`tfhe-ntt`](https://crates.io/crates/tfhe-ntt) crate.

## Installation

Add the following to your `Cargo.toml`:

```toml
[dependencies]
fhe-math = "0.2.0"
```

## Polynomial construction and dot products

Representations are sealed to `PowerBasis`, `Ntt`, and `NttShoup`. Start with
`Poly::<PowerBasis>::from_coefficients`, `from_signed_coefficients`, or
`from_biguint_coefficients`, then use consuming `into_ntt` / `into_ntt_shoup`
transitions. `from_rns_residues` and `from_rns_slice` instead import raw residues
in the selected representation without transforming them. Constructors validate
shape and reduce residues before optimized kernels see them. The old
`rq::traits::TryConvertFrom` extension point is removed.

`dot_product` and `DotProductWorkspace` accept slices. Workspace `_refs` methods
avoid allocating lists for borrowed operands; `_iter` methods consume each
iterator once before validation and arithmetic. `dot_product_into` reuses both
workspace and output storage, leaving them unchanged on validation errors.

## Checked boundaries

`NttOperator::try_forward` and `try_backward` validate length and canonical
residues before changing the slice. Their `_public` counterparts require an
explicit `VariableTime` token. The existing safe slice transforms reject wrong
lengths in release builds too; raw-pointer kernels remain `unsafe`.

Polynomial `from_bytes` uses `DecodeLimits::default()` to bound input and expanded
residue storage; `from_bytes_with_limits` accepts explicit bounds. Decode errors
retain their protobuf source. Detailed classifications are available through
`fhe_math::error`; `Error` and `Result` remain root exports.

## Testing

```bash
cargo test -p fhe-math
```

## Formal verification

`verified/modulus.rs` checks the modulus range and constructs the arithmetic
metadata used by `Modulus::new`: the exact Barrett reciprocal limbs,
leading-zero count, and optimization flag. The reciprocal uses verified
`u128` arithmetic. `zq/modulus.rs` proves the constructor's error mapping and
field storage, carries a persistent type invariant, and verifies cloning,
equality, dereferencing, and all production public wrappers. Private fields
prevent unchecked mutation. The internal runtime backend retains explicit
trusted contracts for distribution setup, sampling, CPU dispatch, and primality;
its implementation and those dependencies are not formally verified.

`zq/primes.rs` proves `generate_prime` returns the greatest eligible prime,
with the requested bit width, strict upper bound, and congruence. It proves
search termination and exact `None` cases under the same trusted primality
contract. Invalid inputs return `None` consistently, and arithmetic alignment
replaces the initial linear scan.

Verus proves conditional selection, then one-subtraction reduction, then
modular addition using those helpers, including the variable-time `reduce1_vt`
and `add_vt` paths. It also proves all four scalar `lazy_reduce*` kernels and
the in-place `lazy_reduce_vec` loop: each result is below twice the modulus
and congruent to its input. The public `primes::supports_opt` check is proved
equivalent to the optimization condition for every `u64`, including the
precondition required by the optimized reducers for valid moduli.
All canonical `reduce*` kernels are also proved, including signed inputs,
optimized and variable-time variants, all vector loops, and the NTT's
`reduce3` helpers: results are the canonical residues below the modulus.
Modular addition, subtraction, and negation have the same canonical-result
guarantee, including their variable-time scalar paths and all six in-place
vector loops. Each operation's scalar and vector kernels live together in
`verified/add.rs`, `verified/sub.rs`, or `verified/neg.rs`.
`verified/mul.rs` proves ordinary, optimized, and Shoup multiplication,
including variable-time variants, all multiplication vector loops, and Shoup
precomputation. Shoup proofs cover any `u64` first operand; lazy products are
below twice the modulus, and canonical products equal `(a*b) mod p`.
`ntt/native/butterflies.rs` proves the four production butterfly methods by
composing reduction and Shoup multiplication. Forward butterflies map inputs
below `4p` to `(x + wy, x - wy)` modulo `p`, with both outputs below `4p`.
Inverse butterflies map inputs below `2p` to `(x + y, z(x - y))` modulo `p`,
with both outputs below `2p`. Both ordinary and variable-time paths are checked,
including overflow/underflow safety. These contracts require `p_twice = 2p`
and valid Shoup twiddles; NTT construction, tables, and transform loops remain
outside the proof.
`ntt/predicates.rs` proves `supports_ntt` is exactly the prime/congruence
eligibility check for power-of-two sizes at least 8. Invalid and oversized
sizes return false before doubling, so both NTT constructors return `None`
for them. It also proves `is_primitive_root` detects exact multiplicative
order `n`: no smaller positive exponent returns one. The power-of-two root
criterion works for composite moduli too. Primality detection retains its
existing trusted contract; root search and root existence are not proved here.
`verified/center.rs` proves scalar and vector centering: each result is
congruent to its input and satisfies `-p <= 2*result < p`, including the
midpoint behavior for both odd and even moduli.
`verified/pow.rs` proves variable-time exponentiation yields `a^n mod p` below
the modulus for canonical bases and the existing `n < p` exponent bound,
including `0^0 = 1`.
`verified/inv.rs` proves inversion by the same Fermat exponentiation, including
the inverse identity and exact `None` cases. Fermat's little theorem is proved
from mathematical primality; the external primality check remains a boundary.
`verified/random.rs` proves the `random_vec` and `random_into` loops preserve
length, canonical bounds, and sample order, with exactly one sampler call per
output entry. These proofs depend on an explicit sampler contract: the
backend sampling implementation, uniformity, and RNG security remain outside
the proof; the adapter that records sample history is verified.
The ordinary Rust implementation shares these verified functions.

`verified/serialize.rs` proves the modulus encoding width, checked byte lengths,
serialization/deserialization composition, and round trips. It imports the
verified packing contracts from `fhe-util` through Cargo Verus. The byte-length
calculation divides before multiplying and checks overflow in release builds
too. These proofs cover the shared `Modulus` serialization implementation;
the wrappers use the checked object invariant. Polynomial/protobuf serialization
remains outside the proof.

Both crates enable Cargo Verus metadata, so proof dependencies use ordinary
crate imports and their specifications ship with the source packages:

```sh
cargo verus verify --fwd-verus-args-to roots -p fhe-util -p fhe-math --locked -- --no-cheating
```

See the [verification guide](../../verification/README.md)
for their exact contracts, pinned toolchain, and proof command. Other operations
are outside the scope of these proofs.

## License

This project is licensed under the [MIT license](https://opensource.org/licenses/MIT).

## Security / Stability

The code in this crate has not undergone an independent security audit.
Use at your own risk.
