# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

Fix a few bugs, remove the need of using nightly, and make some backward-incompatible changes:
- public key generation doesn't return a `Result` anymore
- modify the API to take as input the random number generator.

Additionally, we added more documentation.

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 8 commits contributed to the release over the course of 4 calendar days.
 - 4 days passed between releases.
 - 0 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 7 unique issues were worked on: [#123](https://github.com/tlepoint/fhe.rs/issues/123), [#127](https://github.com/tlepoint/fhe.rs/issues/127), [#130](https://github.com/tlepoint/fhe.rs/issues/130), [#132](https://github.com/tlepoint/fhe.rs/issues/132), [#133](https://github.com/tlepoint/fhe.rs/issues/133), [#134](https://github.com/tlepoint/fhe.rs/issues/134), [#135](https://github.com/tlepoint/fhe.rs/issues/135)

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **[#123](https://github.com/tlepoint/fhe.rs/issues/123)**
    - RGSW mistakenly appeared to depend on a feature; fixes #122. ([`739d4ce`](https://github.com/tlepoint/fhe.rs/commit/739d4ced784ee4aea20c57f3e042361aab7d5517))
 * **[#127](https://github.com/tlepoint/fhe.rs/issues/127)**
    - Computes correctly the number of bits in the plaintext; fixes #126 ([`432586c`](https://github.com/tlepoint/fhe.rs/commit/432586cecf83a0808cf987882c472acbf1330a36))
 * **[#130](https://github.com/tlepoint/fhe.rs/issues/130)**
    - Remove some nightly features, see #117 ([`6361fa3`](https://github.com/tlepoint/fhe.rs/commit/6361fa3ce322b16551cfe4856a49e3933d85c872))
 * **[#132](https://github.com/tlepoint/fhe.rs/issues/132)**
    - Remove the nightly features, except for code coverage and formatting ([`b573138`](https://github.com/tlepoint/fhe.rs/commit/b573138d682e69c3553c2e4ae4a1b7f7a65dbe5d))
 * **[#133](https://github.com/tlepoint/fhe.rs/issues/133)**
    - Explicitely specify the RNG everytime randomness is involved. Fixes #128 ([`8aafe43`](https://github.com/tlepoint/fhe.rs/commit/8aafe4396d0b771e6aa25257c7daa61c109eb367))
 * **[#134](https://github.com/tlepoint/fhe.rs/issues/134)**
    - Remove unnecessary casting by defining more conversions ([`f7cddb3`](https://github.com/tlepoint/fhe.rs/commit/f7cddb358f2ce28483944f99e223c07ae41b0c1c))
 * **[#135](https://github.com/tlepoint/fhe.rs/issues/135)**
    - Starting better documentation ([`13a633c`](https://github.com/tlepoint/fhe.rs/commit/13a633c0f288d27da15548942a061540365aec10))
 * **Uncategorized**
    - Bump all version to beta.3 ([`fc63e4e`](https://github.com/tlepoint/fhe.rs/commit/fc63e4ea6acbb3e9dda83a65cafdf63a081836f2))
</details>

## 0.1.0-beta.2 (2022-09-07)

This release fixes a bug that did not allow to decrypt a modulo-switched ciphertext correctly.

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 2 commits contributed to the release.
 - 0 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' where seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Release fhe v0.1.0-beta.2 ([`d13c33c`](https://github.com/tlepoint/fhe.rs/commit/d13c33caf2850753ed9ef556c41cfaf73700ecd1))
    - Remove forgotten cfg(not(feature ([`2e247f2`](https://github.com/tlepoint/fhe.rs/commit/2e247f235bbe632459259f6ca74a637a2f765187))
</details>

## 0.1.0-beta.1 (2022-09-07)

First version of the `fhe` crate, that includes the BFV cryptosystem.

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 13 commits contributed to the release over the course of 1 calendar day.
 - 0 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 6 unique issues were worked on: [#114](https://github.com/tlepoint/fhe.rs/issues/114), [#115](https://github.com/tlepoint/fhe.rs/issues/115), [#116](https://github.com/tlepoint/fhe.rs/issues/116), [#118](https://github.com/tlepoint/fhe.rs/issues/118), [#120](https://github.com/tlepoint/fhe.rs/issues/120), [#121](https://github.com/tlepoint/fhe.rs/issues/121)

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **[#114](https://github.com/tlepoint/fhe.rs/issues/114)**
    - Rename crates to fhe and fhe-traits ([`9a3d608`](https://github.com/tlepoint/fhe.rs/commit/9a3d6082976a7e0b6f3cec93c096bfaa4a07ebd6))
 * **[#115](https://github.com/tlepoint/fhe.rs/issues/115)**
    - Bump thiserror from 1.0.33 to 1.0.34 ([`e724edf`](https://github.com/tlepoint/fhe.rs/commit/e724edfec78809593e99b21ba5c9eeaaca1a191c))
 * **[#116](https://github.com/tlepoint/fhe.rs/issues/116)**
    - Use zeroizing instead of manual calls to zeroize ([`1d7bc50`](https://github.com/tlepoint/fhe.rs/commit/1d7bc50c58e8807d696d02f3d64e19f34a4ad0c3))
 * **[#118](https://github.com/tlepoint/fhe.rs/issues/118)**
    - Update the README with minimal example and fix compilation error ([`ecba998`](https://github.com/tlepoint/fhe.rs/commit/ecba99898c86a7908a7e9360a6e62826e2ccc5c6))
 * **[#120](https://github.com/tlepoint/fhe.rs/issues/120)**
    - Move internal to crates as they would be published, add changelog ([`cd3ba02`](https://github.com/tlepoint/fhe.rs/commit/cd3ba026d01275672e0c3f5e1d32aa473cde7978))
 * **[#121](https://github.com/tlepoint/fhe.rs/issues/121)**
    - Remove features, remove utilities crate, bump versions ([`570943a`](https://github.com/tlepoint/fhe.rs/commit/570943ae1822888a2ccb27412619ab3355b3ea3a))
 * **Uncategorized**
    - Release fhe v0.1.0-beta.1 ([`718f0cd`](https://github.com/tlepoint/fhe.rs/commit/718f0cdc1e5b75eaf52a5ea1078c1ed9c2bf46f5))
    - First version fhe crate ([`3f9e80c`](https://github.com/tlepoint/fhe.rs/commit/3f9e80c9bc91b068d00ec6b03ccafb07f150185a))
    - Release fhe-traits v0.1.0-beta.0, fhe-util v0.1.0-beta.0, fhe-math v0.1.0-beta.0, fhe v0.1.0-beta.0 ([`e81e1c6`](https://github.com/tlepoint/fhe.rs/commit/e81e1c60769e63c52ad3885d16249161074ca293))
    - Adjusting changelogs prior to release of fhe-traits v0.1.0-beta.0, fhe-util v0.1.0-beta.0, fhe-math v0.1.0-beta.0, fhe v0.1.0-beta.0 ([`4c9ed5b`](https://github.com/tlepoint/fhe.rs/commit/4c9ed5bc57ccaa4a9d9ac98e4883f6c5c2136b5b))
    - Add space to test ([`f5e82f3`](https://github.com/tlepoint/fhe.rs/commit/f5e82f3708bc15a7f517f19bb482fce0044cf091))
    - Update changelog ([`85a00a1`](https://github.com/tlepoint/fhe.rs/commit/85a00a1b8113e4dc8b1d4e9d19fc6c354fb6ae0e))
    - Switch version to a pre-release number ([`cd8d3b2`](https://github.com/tlepoint/fhe.rs/commit/cd8d3b2d383367239436adcc2508bdbe816b9981))
</details>

## 0.0.0 (2023-09-06)

Reserve the name for the `fhe` project.

