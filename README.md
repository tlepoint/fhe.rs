# fhe.rs: Fully Homomorphic Encryption in Rust

[![continuous integration](https://github.com/tlepoint/fhe.rs/actions/workflows/rust.yml/badge.svg?branch=main)](https://github.com/tlepoint/fhe.rs/actions/workflows/rust.yml) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![Code coverage](https://codecov.io/gh/tlepoint/fhe.rs/branch/main/graph/badge.svg?token=LCBSDMB5NS)](https://codecov.io/gh/tlepoint/fhe.rs)

This repository contains the `fhe.rs` library, an experimental cryptographic library in Rust for Ring-LWE-based homomorphic encryption, developed by [Tancrède Lepoint](https://tancre.de).
For more information about the library, see [fhe.rs](https://fhe.rs).

The library features:

* An implementation of a RNS-variant of the Brakerski-Fan-Vercauteren (BFV) homomorphic encryption scheme;
* Performances comparable or better than state-of-the-art libraries in C++ and Go.

> **Warning**
> `fhe.rs` is a beta library, and **should be considered unstable with potential breaking API changes until version 1.0.0 is released!**

> **Note**
> This library is **not** related to the `tfhe-rs` library (a.k.a. `concrete`), Zama's fully homomorphic encryption in Rust, available at [tfhe.rs](https://github.com/zama-ai/tfhe-rs).

## fhe.rs crates

`fhe.rs` is implemented using the Rust programming language. The ecosystem is composed of four public crates (packages):

* [![fhe crate version](https://img.shields.io/crates/v/fhe.svg)](https://crates.io/crates/fhe) [`fhe`](https://crates.io/crates/fhe): This crate contains the implementations of the homomorphic encryption schemes;
* [![fhe-math crate version](https://img.shields.io/crates/v/fhe-math.svg)](https://crates.io/crates/fhe-math) [`fhe-math`](https://crates.io/crates/fhe-math): This crate contains the core mathematical operations for the `fhe` crate;
* [![fhe-traits crate version](https://img.shields.io/crates/v/fhe-traits.svg)](https://crates.io/crates/fhe-traits) [`fhe-traits`](https://crates.io/crates/fhe-traits): This crate contains traits for homomorphic encryption schemes;
* [![fhe-util crate version](https://img.shields.io/crates/v/fhe-util.svg)](https://crates.io/crates/fhe-util) [`fhe-util`](https://crates.io/crates/fhe-util): This crate contains utility functions for the `fhe` crate.

### Installation

To install, add the following to your project's `Cargo.toml` file:

```toml
[dependencies]
fhe = "0.1.0-beta.5"
fhe-traits = "0.1.0-beta.5"
```

## Minimum supported version / toolchain

Rust **1.73** or newer.

## ⚠️ Security / Stability

The implementations contained in the `fhe.rs` ecosystem have never been independently audited for security.
Additionally, no promise on the API and ABI stability will be made until version `1.0.0` of the crates.

Use at your own risk.
