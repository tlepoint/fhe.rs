# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

Fix a few bugs, remove the need of using nightly, and make some backward-incompatible changes by modifying the API to take as input the random number generator.

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 5 commits contributed to the release over the course of 3 calendar days.
 - 4 days passed between releases.
 - 0 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 4 unique issues were worked on: [#130](https://github.com/tlepoint/fhe.rs/issues/130), [#132](https://github.com/tlepoint/fhe.rs/issues/132), [#133](https://github.com/tlepoint/fhe.rs/issues/133), [#134](https://github.com/tlepoint/fhe.rs/issues/134)

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **[#130](https://github.com/tlepoint/fhe.rs/issues/130)**
    - Remove some nightly features, see #117 ([`6361fa3`](https://github.com/tlepoint/fhe.rs/commit/6361fa3ce322b16551cfe4856a49e3933d85c872))
 * **[#132](https://github.com/tlepoint/fhe.rs/issues/132)**
    - Remove the nightly features, except for code coverage and formatting ([`b573138`](https://github.com/tlepoint/fhe.rs/commit/b573138d682e69c3553c2e4ae4a1b7f7a65dbe5d))
 * **[#133](https://github.com/tlepoint/fhe.rs/issues/133)**
    - Explicitely specify the RNG everytime randomness is involved. Fixes #128 ([`8aafe43`](https://github.com/tlepoint/fhe.rs/commit/8aafe4396d0b771e6aa25257c7daa61c109eb367))
 * **[#134](https://github.com/tlepoint/fhe.rs/issues/134)**
    - Remove unnecessary casting by defining more conversions ([`f7cddb3`](https://github.com/tlepoint/fhe.rs/commit/f7cddb358f2ce28483944f99e223c07ae41b0c1c))
 * **Uncategorized**
    - Bump all version to beta.3 ([`fc63e4e`](https://github.com/tlepoint/fhe.rs/commit/fc63e4ea6acbb3e9dda83a65cafdf63a081836f2))
</details>

## 0.1.0-beta.1 (2022-09-07)

Bump pre-release version to match that of `fhe`.

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 9 commits contributed to the release.
 - 0 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 2 unique issues were worked on: [#120](https://github.com/tlepoint/fhe.rs/issues/120), [#121](https://github.com/tlepoint/fhe.rs/issues/121)

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **[#120](https://github.com/tlepoint/fhe.rs/issues/120)**
    - Move internal to crates as they would be published, add changelog ([`cd3ba02`](https://github.com/tlepoint/fhe.rs/commit/cd3ba026d01275672e0c3f5e1d32aa473cde7978))
 * **[#121](https://github.com/tlepoint/fhe.rs/issues/121)**
    - Remove features, remove utilities crate, bump versions ([`570943a`](https://github.com/tlepoint/fhe.rs/commit/570943ae1822888a2ccb27412619ab3355b3ea3a))
 * **Uncategorized**
    - Release fhe-math v0.1.0-beta.1 ([`1b35a2e`](https://github.com/tlepoint/fhe.rs/commit/1b35a2ebd5e2c4d3821e6c967684fdc6e0a77441))
    - Add changelog entry for fhe-math ([`3abb768`](https://github.com/tlepoint/fhe.rs/commit/3abb768ecd236e854bc1c1baa28f2646fb81ecd6))
    - Release fhe-traits v0.1.0-beta.0, fhe-util v0.1.0-beta.0, fhe-math v0.1.0-beta.0, fhe v0.1.0-beta.0 ([`e81e1c6`](https://github.com/tlepoint/fhe.rs/commit/e81e1c60769e63c52ad3885d16249161074ca293))
    - Write changelog ([`ef65eb4`](https://github.com/tlepoint/fhe.rs/commit/ef65eb4b14fd52dfe3796d6c782127d38e551f69))
    - Adjusting changelogs prior to release of fhe-traits v0.1.0-beta.0, fhe-util v0.1.0-beta.0, fhe-math v0.1.0-beta.0, fhe v0.1.0-beta.0 ([`4c9ed5b`](https://github.com/tlepoint/fhe.rs/commit/4c9ed5bc57ccaa4a9d9ac98e4883f6c5c2136b5b))
    - Update changelog ([`85a00a1`](https://github.com/tlepoint/fhe.rs/commit/85a00a1b8113e4dc8b1d4e9d19fc6c354fb6ae0e))
    - Switch version to a pre-release number ([`cd8d3b2`](https://github.com/tlepoint/fhe.rs/commit/cd8d3b2d383367239436adcc2508bdbe816b9981))
</details>

## 0.1.0-beta.0 (2022-09-06)

First pre-release version of the `fhe-math` crate for fhe.rs.

## 0.0.0 (2023-09-06)

Reserve the name for a necessary crate of the `fhe` project.

