

## v0.1.0-beta.5 (2023-09-05)

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 21 commits contributed to the release over the course of 362 calendar days.
 - 363 days passed between releases.
 - 0 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 17 unique issues were worked on: [#130](https://github.com/tlepoint/fhe.rs/issues/130), [#132](https://github.com/tlepoint/fhe.rs/issues/132), [#133](https://github.com/tlepoint/fhe.rs/issues/133), [#134](https://github.com/tlepoint/fhe.rs/issues/134), [#138](https://github.com/tlepoint/fhe.rs/issues/138), [#139](https://github.com/tlepoint/fhe.rs/issues/139), [#140](https://github.com/tlepoint/fhe.rs/issues/140), [#142](https://github.com/tlepoint/fhe.rs/issues/142), [#143](https://github.com/tlepoint/fhe.rs/issues/143), [#144](https://github.com/tlepoint/fhe.rs/issues/144), [#145](https://github.com/tlepoint/fhe.rs/issues/145), [#146](https://github.com/tlepoint/fhe.rs/issues/146), [#148](https://github.com/tlepoint/fhe.rs/issues/148), [#149](https://github.com/tlepoint/fhe.rs/issues/149), [#157](https://github.com/tlepoint/fhe.rs/issues/157), [#168](https://github.com/tlepoint/fhe.rs/issues/168), [#170](https://github.com/tlepoint/fhe.rs/issues/170)

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
 * **[#138](https://github.com/tlepoint/fhe.rs/issues/138)**
    - Bump criterion from 0.3.6 to 0.4.0 ([`db4daf2`](https://github.com/tlepoint/fhe.rs/commit/db4daf29375497bb4331443da86fd520fd71cac8))
 * **[#139](https://github.com/tlepoint/fhe.rs/issues/139)**
    - Bump itertools from 0.10.3 to 0.10.4 ([`75d51b3`](https://github.com/tlepoint/fhe.rs/commit/75d51b3cb1c86fa1603f5e1cd9dc1f3f5859554d))
 * **[#140](https://github.com/tlepoint/fhe.rs/issues/140)**
    - Bump thiserror from 1.0.34 to 1.0.35 ([`195303e`](https://github.com/tlepoint/fhe.rs/commit/195303e93689c20ea6a282e6855762b74145fc59))
 * **[#142](https://github.com/tlepoint/fhe.rs/issues/142)**
    - Bump sha2 from 0.10.5 to 0.10.6 ([`15c3cfe`](https://github.com/tlepoint/fhe.rs/commit/15c3cfead4eea5f8a4610a2d9b571c436e85779b))
 * **[#143](https://github.com/tlepoint/fhe.rs/issues/143)**
    - Bump itertools from 0.10.4 to 0.10.5 ([`72ea4f3`](https://github.com/tlepoint/fhe.rs/commit/72ea4f36256e37e435611e4c28918baeb6e23eae))
 * **[#144](https://github.com/tlepoint/fhe.rs/issues/144)**
    - Bump protobuf from 3.1.0 to 3.2.0 ([`992e8f7`](https://github.com/tlepoint/fhe.rs/commit/992e8f72b42779e96c9cb7828f510a771524aa46))
 * **[#145](https://github.com/tlepoint/fhe.rs/issues/145)**
    - Bump thiserror from 1.0.35 to 1.0.36 ([`2750800`](https://github.com/tlepoint/fhe.rs/commit/27508002ea516d9ba41d0aa756bec7347f8404b2))
 * **[#146](https://github.com/tlepoint/fhe.rs/issues/146)**
    - Bump thiserror from 1.0.36 to 1.0.37 ([`c0022fb`](https://github.com/tlepoint/fhe.rs/commit/c0022fb24bdb191f5ebfe6a2dce31fe6d5b34523))
 * **[#148](https://github.com/tlepoint/fhe.rs/issues/148)**
    - Bump crypto-bigint from 0.4.8 to 0.4.9 ([`e5c94a2`](https://github.com/tlepoint/fhe.rs/commit/e5c94a22d0e0114869d390dc4b2e1e7d0c9d9dcd))
 * **[#149](https://github.com/tlepoint/fhe.rs/issues/149)**
    - Use workspace dependencies et versions + Release 0.1.0-beta4 ([`a0287ba`](https://github.com/tlepoint/fhe.rs/commit/a0287ba3842fcf19b45fd380c56ba7b5e52a387b))
 * **[#157](https://github.com/tlepoint/fhe.rs/issues/157)**
    - Fix clippy warnings. Fixes #154 ([`3d1d9e8`](https://github.com/tlepoint/fhe.rs/commit/3d1d9e8e28853e468dfe974253834d9c242fd0a5))
 * **[#168](https://github.com/tlepoint/fhe.rs/issues/168)**
    - Fix clippy error ([`ee67419`](https://github.com/tlepoint/fhe.rs/commit/ee6741906cfc90ab21c71953b8a1d35e90c56003))
 * **[#170](https://github.com/tlepoint/fhe.rs/issues/170)**
    - Change tabs into space, optimize ntt operator constructor ([`393316f`](https://github.com/tlepoint/fhe.rs/commit/393316ffe1d02efe70e26310ff04318b2e185e87))
 * **Uncategorized**
    - Change concrete to tfhe ([`33146f7`](https://github.com/tlepoint/fhe.rs/commit/33146f77acdc64b7c5a32494d1bd575b6bc9910f))
    - Update dependencies ([`0b4cc99`](https://github.com/tlepoint/fhe.rs/commit/0b4cc9922e6f3bfe0a16f8bf4c7f7b7be548a464))
    - Release fhe-traits v0.1.0-beta.3, fhe-util v0.1.0-beta.3, fhe-math v0.1.0-beta.3, fhe v0.1.0-beta.3 ([`c031c0e`](https://github.com/tlepoint/fhe.rs/commit/c031c0eca3a354e7d1e016dc7da2fba27f061f08))
    - Bump all version to beta.3 ([`913f84d`](https://github.com/tlepoint/fhe.rs/commit/913f84d9f510602283716a5ff310215734337956))
</details>

## v0.1.0-beta.1 (2022-09-07)

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

