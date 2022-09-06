<h1 align="center">fhe.rs: Fully Homomorphic Encryption in Rust</h1>
<p align="center">
<a href="https://github.com/tlepoint/fhe.rs/actions/workflows/rust.yml"><img src="https://github.com/tlepoint/fhe.rs/actions/workflows/rust.yml/badge.svg?branch=main"/></a>
<a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-yellow.svg"/></a>
<a href="https://codecov.io/gh/tlepoint/fhe.rs"><img src="https://codecov.io/gh/tlepoint/fhe.rs/branch/main/graph/badge.svg?token=LCBSDMB5NS"/></a>
</p>

This repository contains the `fhe.rs` library, a cryptographic library in Rust for Ring-LWE-based homomorphic encryption, developed by [Tancrède Lepoint](https://tancre.de). For more information about the library, see [fhe.rs](https://fhe.rs).

The library features:

* An implementation of a RNS-variant of the Brakerski-Fan-Vercauteren (BFV) homomorphic encryption scheme;
* Performances comparable or better than state-of-the-art libraries in C++ and Go.

> **Warning**
> `fhe.rs` is a beta library, and **should be considered unstable with potential breaking API changes until version 1.0.0 is released!**

> **Note**
> This library is not related to the `concrete` ecosystems (Zama's fully homomorphic encryption in Rust), available at [concrete.rs](https://concrete.rs).

## fhe.rs crates

`fhe.rs` is implemented using the Rust programming language. The ecosystem is composed of two public crates (packages):

* [`fhe`](https://github.com/tlepoint/fhe.rs/tree/main/crates/fhe): This crate contains the implementations of the homomorphic encryption schemes;
* [`fhe-traits`](https://github.com/tlepoint/fhe.rs/tree/main/crates/fhe-traits): This crate contains traits for homomorphic encryption schemes.

The repository also contain several internal (private) crates, used by the public crates, implementing the core cryptographic and mathematical operations.

## Minimum supported version

The `fhe` crate requires the `nightly` toolchain as it enables multiple unstable features. The minimal supported version will be changed to a stable version in a future update.

## Installation

To use the latest published crate, add the following to your `Cargo.toml` file:

```toml
[dependencies]
fhe = "0.0.1-alpha"
```

## ⚠️ Security warning

The implementations contained in the `fhe.rs` ecosystem have never been independently audited for security. Additionally, no promise on the API and ABI stability will be made until version `1.0.0` of the crates.

Use at your own risk.
