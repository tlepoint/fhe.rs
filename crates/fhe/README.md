# fhe [![fhe crate version](https://img.shields.io/crates/v/fhe.svg)](https://crates.io/crates/fhe)

**A pure-Rust implementation of fully homomorphic encryption schemes based on Ring-LWE.**

This library implements [Fully Homomorphic Encryption](https://en.wikipedia.org/wiki/Homomorphic_encryption#Fully_homomorphic_encryption) schemes, i.e., encryption schemes which perform implicit additions and multiplications on plaintext values while exclusively manipulating encrypted data.

This library provides implementations of:

* BFV, the Brakerski-Fan-Vercauteren (BFV) homomorphic encryption scheme.
More precisely, this library implements a leveled variant of the [HPS](https://eprint.iacr.org/2018/117) (Halevi--Polyakov--Shoup) RNS-variant of the scheme.

## Example

Below is a simple example using BFV of an homomorphic multiplication.
One ciphertext encrypts the value `20` using the secret key, and one ciphertext encrypts the value `-7` using the public key. The ciphertexts are then multiplied, and after decryption, the program checks that the decrypted value has `20 * (-7) = -140` in the first coefficient.

```rust
use fhe::bfv::{BfvParametersBuilder, Ciphertext, Encoding, Plaintext, PublicKey, SecretKey};
use fhe_traits::*;
use rand::{rngs::OsRng, thread_rng};
use std::sync::Arc;

fn main() {
    let parameters = Arc::new(
        BfvParametersBuilder::new()
            .set_degree(2048)
            .set_moduli(&[0x3fffffff000001])
            .set_plaintext_modulus(1 << 10)
            .build()
            .unwrap(),
    );
    let mut rng = thread_rng();

    let secret_key = SecretKey::random(&parameters, &mut OsRng);
    let public_key = PublicKey::new(&secret_key, &mut rng);

    let plaintext_1 = Plaintext::try_encode(&[20_u64], Encoding::poly(), &parameters).unwrap();
    let plaintext_2 = Plaintext::try_encode(&[-7_i64], Encoding::poly(), &parameters).unwrap();

    let ciphertext_1: Ciphertext = secret_key.try_encrypt(&plaintext_1, &mut rng).unwrap();
    let ciphertext_2: Ciphertext = public_key.try_encrypt(&plaintext_2, &mut rng).unwrap();

    let result = &ciphertext_1 * &ciphertext_2;

    let decrypted_plaintext = secret_key.try_decrypt(&result).unwrap();
    let decrypted_vector = Vec::<i64>::try_decode(&decrypted_plaintext, Encoding::poly()).unwrap();

    assert_eq!(decrypted_vector[0], -140);
}
```

Note that operations actually happen modulo the `plaintext_modulus`, here set to `1024 (= 1 << 10)`; for example, we would have had that the homomorphic multiplication of `805` and `-7` is `509 = (805 * (-7)) mod 1024`. Additionally, the `poly()` encoding means that the vector being encoded corresponds to the coefficients of a polynomial in `(ZZ / (1024))[x] / (x^2048+1)` (and homomorphic multiplication happens in that ring); here since only one coefficient is provided, the value is placed in the constant coefficient. The library also contains a `simd()` encoding, which enables component-wise operation on the values of the vector, provided the technical limitation that the plaintext modulus is congruent to `1` modulo twice the polynomial degree.

## Examples

More examples exercizing multiple functions from the API are provided in the repository [`examples/`](./examples/). For example, this library implements [SealPIR](https://eprint.iacr.org/2017/1142) and [MulPIR](https://eprint.iacr.org/2019/1483), which can be run as follows:

```bash
cargo run --release --example sealpir
```

and

```bash
cargo run --release --example mulpir
```

## Performance

Micro benchmarks can be obtained by running `cargo bench`. This crate uses [criterion.rs](https://criterion.rs) for benchmarks.

## Unit tests

Run tests with `cargo test`.

## ⚠️ Security / Stability

The implementations in this crate have never been independently audited for security.
Additionally, no promise on the API and ABI stability will be made until version `1.0.0` of the crate.

Use at your own risk.
