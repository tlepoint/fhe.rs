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

## Testing

```bash
cargo test -p fhe-math
```

## License

This project is licensed under the [MIT license](https://opensource.org/licenses/MIT).

## Security / Stability

The code in this crate has not undergone an independent security audit.
Use at your own risk.
