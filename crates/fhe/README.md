# fhe [![fhe crate version](https://img.shields.io/crates/v/fhe.svg)](https://crates.io/crates/fhe) [![documentation](https://docs.rs/fhe/badge.svg)](https://docs.rs/fhe)

**A pure-Rust implementation of fully homomorphic encryption schemes based on Ring-LWE.**

This library implements [Fully Homomorphic Encryption](https://en.wikipedia.org/wiki/Homomorphic_encryption#Fully_homomorphic_encryption) schemes, i.e., encryption schemes which perform implicit additions and multiplications on plaintext values while exclusively manipulating encrypted data.

This library provides implementations of:

* BFV, the Brakerski-Fan-Vercauteren (BFV) homomorphic encryption scheme.
  More precisely, this library implements a leveled variant of the [HPS](https://eprint.iacr.org/2018/117) (Halevi--Polyakov--Shoup) RNS-variant of the scheme.

## Installation

Add the following to your `Cargo.toml`:

```toml
[dependencies]
fhe = "0.2.0"
```

## Cargo features

* `tfhe-ntt` enables the accelerated NTT implementation from `tfhe-ntt`.
* `experimental-mbfv` exposes the incomplete multiparty BFV APIs. These APIs
  have additional unresolved security requirements and must not be used in
  production or to protect sensitive data.

## Example

Below is a simple example using BFV of an homomorphic multiplication.
One ciphertext encrypts the value `20` using the secret key, and one ciphertext encrypts the value `-7` using the public key. The ciphertexts are then multiplied, and after decryption, the program checks that the decrypted value has `20 * (-7) = -140` in the first coefficient.

```rust
use fhe::bfv::{Encoding, Parameters, Plaintext, PublicKey, SecretKey};

fn main() -> fhe::Result<()> {
    let parameters = Parameters::builder()
        .degree(2048)
        .ciphertext_moduli([0x3fffffff000001])
        .plaintext_modulus(1024_u64)
        .build()?;
    let mut rng = rand::rng();
    let secret_key = SecretKey::generate(&parameters, &mut rng);
    let public_key = PublicKey::from_secret_key(&secret_key, &mut rng);

    let a = Plaintext::encode(&parameters, &[20_u64], Encoding::Polynomial)?;
    let b = Plaintext::encode_signed(&parameters, &[-7_i64], Encoding::Polynomial)?;
    let encrypted_a = secret_key.encrypt(&a, &mut rng)?;
    let encrypted_b = public_key.encrypt(&b, &mut rng)?;
    let product = encrypted_a.multiply(&encrypted_b)?;
    let decoded = secret_key.decrypt(&product)?.decode_signed(Encoding::Polynomial)?;
    assert_eq!(decoded[0], -140);
    Ok(())
}
```

Note that operations actually happen modulo the `plaintext_modulus`, here set to `1024 (= 1 << 10)`; for example, we would have had that the homomorphic multiplication of `805` and `-7` is `509 = (805 * (-7)) mod 1024`. Additionally, the `Encoding::Polynomial` encoding means that the vector being encoded corresponds to the coefficients of a polynomial in `(ZZ / (1024))[x] / (x^2048+1)` (and homomorphic multiplication happens in that ring); here since only one coefficient is provided, the value is placed in the constant coefficient. The library also contains a `Encoding::Simd` encoding, which enables component-wise operation on the values of the vector, provided the technical limitation that the plaintext modulus is congruent to `1` modulo twice the polynomial degree.

## Advanced APIs

The basic flow uses `bfv::{Parameters, Encoding, Plaintext, Ciphertext, SecretKey,
PublicKey}`. Evaluation keys, immutable multiplication plans, prepared operands,
and workspaces live in `bfv::evaluation`. Compact plaintext storage lives in
`bfv::packing`; immutable level descriptors live in `bfv::context`.

`Plaintext::encode_chunks` returns `Vec<Plaintext>`; empty input produces one
zero plaintext. `PackedPlaintextBatch` keeps contiguous storage and validates
`try_extend` before appending. Dot-product methods accept slices, with `_refs`
variants for reference slices and `_iter` variants that collect references once.
Workspaces remain caller-owned mutable values; no implicit thread pool is created.

## Examples

More examples exercizing multiple functions from the API are provided in the repository [`examples/`](./examples/). For example, this library implements [SealPIR](https://eprint.iacr.org/2017/1142) and [MulPIR](https://eprint.iacr.org/2019/1483), which can be run as follows:

```bash
cargo run --release --example sealpir
```

and

```bash
cargo run --release --example mulpir
```

## Key configuration and imports

Use `RelinearizationKey::builder(&secret_key).ciphertext_level(level).key_level(key_level).build(&mut rng)`
for named levels; its type lives in `bfv::evaluation`. Both levels default to zero.
The key level must not exceed the ciphertext level.

Byte imports apply `DecodeLimits::default()` before allocating wire objects and
expanded polynomial storage. Use `from_bytes_with_limits` to supply explicit
bounds for a larger workload. Limits cover input bytes, degree, modulus count,
big plaintext modulus bytes, polynomial slots, and charged residue storage;
they are not a process memory quota. Detailed errors live in `fhe::error`, while
`Error` and `Result` remain at the root. Decode errors retain their underlying
protobuf cause through `std::error::Error::source()`.

## Performance

Micro benchmarks can be obtained by running `cargo bench`. This crate uses [criterion.rs](https://criterion.rs) for benchmarks.

For repeated ciphertext multiplication, prepare the operand that stays fixed.
This caches its basis extension and Shoup multiplication tables. The prepared
value owns a snapshot and borrows the multiplication strategy, so its level
and scaling settings cannot accidentally change underneath it.

```rust
use fhe::bfv::{Ciphertext, evaluation::{MultiplicationPlan, RelinearizationKey}};

fn multiply_query(
    query: &Ciphertext,
    encrypted_rows: &[Ciphertext],
    key: &RelinearizationKey,
) -> fhe::Result<Vec<Ciphertext>> {
    let strategy = MultiplicationPlan::with_relinearization(key)?;
    let prepared = strategy.prepare_lhs(query)?;
    encrypted_rows.iter().map(|row| prepared.multiply(row)).collect()
}
```

Preparation trades memory and setup time for faster subsequent products. It
works with custom asymmetric scaling strategies and optional modulus switching;
every right operand must use compatible parameter settings and the strategy's input
level. Use `MultiplicationPlan::builder(&parameters).level(level).build()?` to reuse
the parameters' multiplication basis without a relinearization key.

`Ciphertext::multiply` provides checked multiplication, and `Ciphertext::square`
uses symmetry to avoid duplicate cross products. These return unrelinearized
ciphertexts. `MultiplicationPlan::square` also applies the strategy's relinearization
and modulus switching. Explicit `square()` is preferable when squaring a cloned
ciphertext: multiplication methods recognize identical references without scanning encrypted
coefficients for equality.

Run `cargo bench -p fhe --bench bfv_multiplication` to compare ordinary and
prepared products, including preparation cost for a batch of eight products.

## Unit tests

Run tests with `cargo test`.

## ⚠️ Security / Stability

The implementations in this crate have never been independently audited for security.

Use at your own risk.
