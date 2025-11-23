# fhe-math [![crate version](https://img.shields.io/crates/v/fhe-math.svg)](https://crates.io/crates/fhe-math) [![documentation](https://docs.rs/fhe-math/badge.svg)](https://docs.rs/fhe-math)

Core mathematical primitives for the [`fhe.rs`](https://github.com/tlepoint/fhe.rs) ecosystem.

This crate exposes building blocks such as number theoretic transforms (NTT), residue number system (RNS) arithmetic, and ring arithmetic over `Z_q` that are used by higher level crates like [`fhe`](https://crates.io/crates/fhe).

## Features

* `ntt`, `rns`, `rq`, and `zq` modules for modular arithmetic over large rings.
* Optional `tfhe-ntt` and `tfhe-ntt-nightly` features to enable hardware accelerated NTTs via the [`tfhe-ntt`](https://crates.io/crates/tfhe-ntt) crate.

## Installation

Add the following to your `Cargo.toml`:

```toml
[dependencies]
fhe-math = "0.1.1"
```

## Testing

```bash
cargo test -p fhe-math
```

## License

This project is licensed under the [MIT license](https://opensource.org/licenses/MIT).

## Security / Stability

The code in this crate has not undergone an independent security audit.
Use at your own risk.
