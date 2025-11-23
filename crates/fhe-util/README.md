# fhe-util [![crate version](https://img.shields.io/crates/v/fhe-util.svg)](https://crates.io/crates/fhe-util) [![documentation](https://docs.rs/fhe-util/badge.svg)](https://docs.rs/fhe-util)

Utility functions for the [`fhe.rs`](https://github.com/tlepoint/fhe.rs) ecosystem.

The crate contains helper routines such as primality testing, centered binomial sampling, modular arithmetic helpers and other small utilities relied upon by the [`fhe`](https://crates.io/crates/fhe) and `fhe-math` crates.

## Installation

```toml
[dependencies]
fhe-util = "0.1.1"
```

## Testing

```bash
cargo test -p fhe-util
```

## License

This project is licensed under the [MIT license](https://opensource.org/licenses/MIT).

## Security / Stability

The code in this crate has not undergone an independent security audit.
Use at your own risk.
