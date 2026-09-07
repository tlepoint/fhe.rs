# P0 protobuf fixtures

Generated from commit `0eb11ce` before the P1 API migration. The secret key is
public, deterministic test data, generated with `ChaCha8Rng::seed_from_u64(0xf100)`.

Parameters: degree 16, plaintext modulus 1153, ciphertext modulus bit lengths
`[50, 50, 50]`, default variance 10. Using that one RNG, generation order is:
secret key, public key, symmetric ciphertext of SIMD `[2, 3]`, relinearization
key, RGSW ciphertext of the same plaintext, evaluation key with inner sum and
four expansion steps. All objects are at level zero.

Files contain the original `fhe_traits::Serialize::to_bytes` output. The P1 tests
check import, re-export, and evaluation. Evaluation-key protobuf map iteration
was not deterministic, so its byte ordering is not used as an equality contract.

The root fixtures use the native NTT backend. The `tfhe/` fixtures use
`--features tfhe-ntt` on the same P0 commit and the same generator inputs.
The existing format does not ensure interoperability across those backends:
the P0 `tfhe-ntt` build also fails to recover SIMD `[2, 3]` from the native
ciphertext fixture. Tests select the fixture set matching the build's backend.
