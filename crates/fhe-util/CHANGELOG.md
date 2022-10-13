# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 0.1.0-beta.4 (2023-10-13)

Bump dependencies versions and use workspace dependencies introduced in Rust 1.64.

## 0.1.0-beta.3 (2022-09-11)

Fix a few bugs, remove the need of using nightly, and make some backward-incompatible changes by modifying the API to take as input the random number generator.

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 4 commits contributed to the release over the course of 3 calendar days.
 - 4 days passed between releases.
 - 0 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 3 unique issues were worked on: [#130](https://github.com/tlepoint/fhe.rs/issues/130), [#132](https://github.com/tlepoint/fhe.rs/issues/132), [#133](https://github.com/tlepoint/fhe.rs/issues/133)

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **[#130](https://github.com/tlepoint/fhe.rs/issues/130)**
    - Remove some nightly features, see #117 ([`6361fa3`](https://github.com/tlepoint/fhe.rs/commit/6361fa3ce322b16551cfe4856a49e3933d85c872))
 * **[#132](https://github.com/tlepoint/fhe.rs/issues/132)**
    - Remove the nightly features, except for code coverage and formatting ([`b573138`](https://github.com/tlepoint/fhe.rs/commit/b573138d682e69c3553c2e4ae4a1b7f7a65dbe5d))
 * **[#133](https://github.com/tlepoint/fhe.rs/issues/133)**
    - Explicitely specify the RNG everytime randomness is involved. Fixes #128 ([`8aafe43`](https://github.com/tlepoint/fhe.rs/commit/8aafe4396d0b771e6aa25257c7daa61c109eb367))
 * **Uncategorized**
    - Bump all version to beta.3 ([`b300590`](https://github.com/tlepoint/fhe.rs/commit/b3005904a62d1e39c1dde908054f24d1d96e8547))
</details>

## 0.1.0-beta.1 (2022-09-07)

Bump pre-release version to match that of `fhe`.

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 3 commits contributed to the release.
 - 0 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 1 unique issue was worked on: [#121](https://github.com/tlepoint/fhe.rs/issues/121)

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **[#121](https://github.com/tlepoint/fhe.rs/issues/121)**
    - Remove features, remove utilities crate, bump versions ([`570943a`](https://github.com/tlepoint/fhe.rs/commit/570943ae1822888a2ccb27412619ab3355b3ea3a))
 * **Uncategorized**
    - Release fhe-util v0.1.0-beta.1 ([`49d32b7`](https://github.com/tlepoint/fhe.rs/commit/49d32b737bf3e943aab7861c375e269cb7971740))
    - Bump version fhe-util ([`136134c`](https://github.com/tlepoint/fhe.rs/commit/136134ccefb563780a34c356c2a646f18285630b))
</details>

## 0.1.0-beta.0 (2022-09-06)

First pre-release version of the `fhe-util` crate for fhe.rs.

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 4 commits contributed to the release.
 - 0 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 1 unique issue was worked on: [#120](https://github.com/tlepoint/fhe.rs/issues/120)

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **[#120](https://github.com/tlepoint/fhe.rs/issues/120)**
    - Move internal to crates as they would be published, add changelog ([`cd3ba02`](https://github.com/tlepoint/fhe.rs/commit/cd3ba026d01275672e0c3f5e1d32aa473cde7978))
 * **Uncategorized**
    - Release fhe-traits v0.1.0-beta.0, fhe-util v0.1.0-beta.0, fhe-math v0.1.0-beta.0, fhe v0.1.0-beta.0 ([`e81e1c6`](https://github.com/tlepoint/fhe.rs/commit/e81e1c60769e63c52ad3885d16249161074ca293))
    - Update changelog ([`85a00a1`](https://github.com/tlepoint/fhe.rs/commit/85a00a1b8113e4dc8b1d4e9d19fc6c354fb6ae0e))
    - Switch version to a pre-release number ([`cd8d3b2`](https://github.com/tlepoint/fhe.rs/commit/cd8d3b2d383367239436adcc2508bdbe816b9981))
</details>

## 0.0.0 (2023-09-06)

Reserve the name for a necessary crate of the `fhe` project.

