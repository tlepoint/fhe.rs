# fhe.rs: Fully Homomorphic Encryption in Rust

[![continuous integration](https://github.com/tlepoint/fhe.rs/actions/workflows/rust.yml/badge.svg?branch=main)](https://github.com/tlepoint/fhe.rs/actions/workflows/rust.yml) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![Code coverage](https://codecov.io/gh/tlepoint/fhe.rs/branch/main/graph/badge.svg?token=LCBSDMB5NS)](https://codecov.io/gh/tlepoint/fhe.rs)

This repository contains the `fhe.rs` library, an experimental cryptographic library in Rust for Ring-LWE-based homomorphic encryption, developed by [Tancrède Lepoint](https://tancre.de). For more information about the library, see [fhe.rs](https://fhe.rs).

The library features:

* An implementation of a RNS-variant of the Brakerski-Fan-Vercauteren (BFV) homomorphic encryption scheme;
* Performances comparable or better than state-of-the-art libraries in C++ and Go.

> **Warning**
> `fhe.rs` is a beta library, and **should be considered unstable with potential breaking API changes until version 1.0.0 is released!**

> **Note**
> This library is not related to the `concrete` ecosystems (Zama's fully homomorphic encryption in Rust), available at [concrete.rs](https://concrete.rs).

## fhe.rs crates

`fhe.rs` is implemented using the Rust programming language. The ecosystem is composed of four public crates (packages):

* [`fhe`](https://github.com/tlepoint/fhe.rs/tree/main/crates/fhe): This crate contains the implementations of the homomorphic encryption schemes;
* [`fhe-math`](https://github.com/tlepoint/fhe.rs/tree/main/crates/fhe-math): This crate contains the core mathematical operations for the `fhe` crate;
* [`fhe-traits`](https://github.com/tlepoint/fhe.rs/tree/main/crates/fhe-traits): This crate contains traits for homomorphic encryption schemes.
* [`fhe-util`](https://github.com/tlepoint/fhe.rs/tree/main/crates/fhe-util): This crate contains utility functions for the `fhe` crate.

## Minimum supported version / toolchain

The `fhe` crate requires the `nightly` toolchain as it uses multiple unstable features. The minimal supported version will be changed to a stable version in a future update.

## Installation

To use the latest published crate, add one or both of the following to your `Cargo.toml` file:

```toml
[dependencies]
fhe = "0.1.0-beta.1"
fhe-traits = "0.1.0-beta.1"
```

## Getting started

Below is a simple example of an homomorphic multiplication modulo `1024 (= 1 << 10)`, where one ciphertext is encrypted using the secret key, and one ciphertext encrypted using the public key. The `poly()` encoding means that the vector being encoded corresponds to the coefficients of a polynomial in `(ZZ / (1024))[x] / (x^2048+1)`.

```rust
use fhe::bfv;
use fhe_traits::*;
use std::sync::Arc;

fn main() {
    let parameters = Arc::new(
        bfv::BfvParametersBuilder::new()
            .set_degree(2048)
            .set_moduli(&[0x3fffffff000001])
            .set_plaintext_modulus(1 << 10)
            .build()
            .unwrap(),
    );

    let secret_key = bfv::SecretKey::random(&parameters);
    let public_key = bfv::PublicKey::new(&secret_key).unwrap();

    let plaintext_1 =
        bfv::Plaintext::try_encode(&[20_u64] as &[u64], bfv::Encoding::poly(), &parameters)
            .unwrap();
    let plaintext_2 =
        bfv::Plaintext::try_encode(&[-7_i64] as &[i64], bfv::Encoding::poly(), &parameters)
            .unwrap();

    let ciphertext_1: bfv::Ciphertext = secret_key.try_encrypt(&plaintext_1).unwrap();
    let ciphertext_2: bfv::Ciphertext = public_key.try_encrypt(&plaintext_2).unwrap();

    let result = &ciphertext_1 * &ciphertext_2;

    let decrypted_plaintext = secret_key.try_decrypt(&result).unwrap();
    let decrypted_vector =
        Vec::<u64>::try_decode(&decrypted_plaintext, bfv::Encoding::poly()).unwrap();

    assert_eq!(decrypted_vector[0], ((1 << 10) - 20 * 7));
}
```

## ⚠️ Security warning

The implementations contained in the `fhe.rs` ecosystem have never been independently audited for security.
Additionally, no promise on the API and ABI stability will be made until version `1.0.0` of the crates.

Use at your own risk.
