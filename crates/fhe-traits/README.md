# fhe-traits [![crate version](https://img.shields.io/crates/v/fhe-traits.svg)](https://crates.io/crates/fhe-traits) [![documentation](https://docs.rs/fhe-traits/badge.svg)](https://docs.rs/fhe-traits)

Traits defining the interface for fully homomorphic encryption types and operations.

This crate provides common abstractions for parameters, plaintext and ciphertext representations, encoding, encryption, decryption and serialization used throughout the [`fhe.rs`](https://github.com/tlepoint/fhe.rs) crates.

## Installation

```toml
[dependencies]
fhe-traits = "0.1.0-beta.8"
```

## Testing

```bash
cargo test -p fhe-traits
```

## License

This project is licensed under the [MIT license](https://opensource.org/licenses/MIT).

## Security / Stability

The code in this crate has not undergone an independent security audit and the API may change before version `1.0.0`.
Use at your own risk.
