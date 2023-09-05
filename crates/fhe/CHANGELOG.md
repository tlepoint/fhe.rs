

## v0.1.0-beta.5 (2023-09-05)

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 3 commits contributed to the release over the course of 250 calendar days.
 - 327 days passed between releases.
 - 0 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 2 unique issues were worked on: [#157](https://github.com/tlepoint/fhe.rs/issues/157), [#170](https://github.com/tlepoint/fhe.rs/issues/170)

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **[#157](https://github.com/tlepoint/fhe.rs/issues/157)**
    - Fix clippy warnings. Fixes #154 ([`3d1d9e8`](https://github.com/tlepoint/fhe.rs/commit/3d1d9e8e28853e468dfe974253834d9c242fd0a5))
 * **[#170](https://github.com/tlepoint/fhe.rs/issues/170)**
    - Change tabs into space, optimize ntt operator constructor ([`393316f`](https://github.com/tlepoint/fhe.rs/commit/393316ffe1d02efe70e26310ff04318b2e185e87))
 * **Uncategorized**
    - Change concrete to tfhe ([`33146f7`](https://github.com/tlepoint/fhe.rs/commit/33146f77acdc64b7c5a32494d1bd575b6bc9910f))
</details>

## v0.1.0-beta.4 (2022-10-13)

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 20 commits contributed to the release over the course of 36 calendar days.
 - 36 days passed between releases.
 - 0 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 18 unique issues were worked on: [#123](https://github.com/tlepoint/fhe.rs/issues/123), [#127](https://github.com/tlepoint/fhe.rs/issues/127), [#130](https://github.com/tlepoint/fhe.rs/issues/130), [#132](https://github.com/tlepoint/fhe.rs/issues/132), [#133](https://github.com/tlepoint/fhe.rs/issues/133), [#134](https://github.com/tlepoint/fhe.rs/issues/134), [#135](https://github.com/tlepoint/fhe.rs/issues/135), [#136](https://github.com/tlepoint/fhe.rs/issues/136), [#138](https://github.com/tlepoint/fhe.rs/issues/138), [#139](https://github.com/tlepoint/fhe.rs/issues/139), [#140](https://github.com/tlepoint/fhe.rs/issues/140), [#141](https://github.com/tlepoint/fhe.rs/issues/141), [#143](https://github.com/tlepoint/fhe.rs/issues/143), [#144](https://github.com/tlepoint/fhe.rs/issues/144), [#145](https://github.com/tlepoint/fhe.rs/issues/145), [#146](https://github.com/tlepoint/fhe.rs/issues/146), [#147](https://github.com/tlepoint/fhe.rs/issues/147), [#149](https://github.com/tlepoint/fhe.rs/issues/149)

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
 * **[#136](https://github.com/tlepoint/fhe.rs/issues/136)**
    - Add comments to SealPIR and MulPIR ([`f374841`](https://github.com/tlepoint/fhe.rs/commit/f374841f5d9d0bf2ce43fb0c4043d341edf68564))
 * **[#138](https://github.com/tlepoint/fhe.rs/issues/138)**
    - Bump criterion from 0.3.6 to 0.4.0 ([`db4daf2`](https://github.com/tlepoint/fhe.rs/commit/db4daf29375497bb4331443da86fd520fd71cac8))
 * **[#139](https://github.com/tlepoint/fhe.rs/issues/139)**
    - Bump itertools from 0.10.3 to 0.10.4 ([`75d51b3`](https://github.com/tlepoint/fhe.rs/commit/75d51b3cb1c86fa1603f5e1cd9dc1f3f5859554d))
 * **[#140](https://github.com/tlepoint/fhe.rs/issues/140)**
    - Bump thiserror from 1.0.34 to 1.0.35 ([`195303e`](https://github.com/tlepoint/fhe.rs/commit/195303e93689c20ea6a282e6855762b74145fc59))
 * **[#141](https://github.com/tlepoint/fhe.rs/issues/141)**
    - Bump indicatif from 0.17.0 to 0.17.1 ([`02f0674`](https://github.com/tlepoint/fhe.rs/commit/02f0674fb24ee8f3c6129bddbe91edae0a9d7808))
 * **[#143](https://github.com/tlepoint/fhe.rs/issues/143)**
    - Bump itertools from 0.10.4 to 0.10.5 ([`72ea4f3`](https://github.com/tlepoint/fhe.rs/commit/72ea4f36256e37e435611e4c28918baeb6e23eae))
 * **[#144](https://github.com/tlepoint/fhe.rs/issues/144)**
    - Bump protobuf from 3.1.0 to 3.2.0 ([`992e8f7`](https://github.com/tlepoint/fhe.rs/commit/992e8f72b42779e96c9cb7828f510a771524aa46))
 * **[#145](https://github.com/tlepoint/fhe.rs/issues/145)**
    - Bump thiserror from 1.0.35 to 1.0.36 ([`2750800`](https://github.com/tlepoint/fhe.rs/commit/27508002ea516d9ba41d0aa756bec7347f8404b2))
 * **[#146](https://github.com/tlepoint/fhe.rs/issues/146)**
    - Bump thiserror from 1.0.36 to 1.0.37 ([`c0022fb`](https://github.com/tlepoint/fhe.rs/commit/c0022fb24bdb191f5ebfe6a2dce31fe6d5b34523))
 * **[#147](https://github.com/tlepoint/fhe.rs/issues/147)**
    - Bump console from 0.15.1 to 0.15.2 ([`a344fab`](https://github.com/tlepoint/fhe.rs/commit/a344fabf2382dbabca7368cb32b44907735481fc))
 * **[#149](https://github.com/tlepoint/fhe.rs/issues/149)**
    - Use workspace dependencies et versions + Release 0.1.0-beta4 ([`a0287ba`](https://github.com/tlepoint/fhe.rs/commit/a0287ba3842fcf19b45fd380c56ba7b5e52a387b))
 * **Uncategorized**
    - Release fhe-traits v0.1.0-beta.3, fhe-util v0.1.0-beta.3, fhe-math v0.1.0-beta.3, fhe v0.1.0-beta.3 ([`c031c0e`](https://github.com/tlepoint/fhe.rs/commit/c031c0eca3a354e7d1e016dc7da2fba27f061f08))
    - Bump all version to beta.3 ([`913f84d`](https://github.com/tlepoint/fhe.rs/commit/913f84d9f510602283716a5ff310215734337956))
</details>

## v0.1.0-beta.2 (2022-09-07)

### Commit Statistics

<csr-read-only-do-not-edit/>

 - 2 commits contributed to the release.
 - 0 commits were understood as [conventional](https://www.conventionalcommits.org).
 - 0 issues like '(#ID)' were seen in commit messages

### Commit Details

<csr-read-only-do-not-edit/>

<details><summary>view details</summary>

 * **Uncategorized**
    - Release fhe v0.1.0-beta.2 ([`d13c33c`](https://github.com/tlepoint/fhe.rs/commit/d13c33caf2850753ed9ef556c41cfaf73700ecd1))
    - Remove forgotten cfg(not(feature ([`2e247f2`](https://github.com/tlepoint/fhe.rs/commit/2e247f235bbe632459259f6ca74a637a2f765187))
</details>

## v0.1.0-beta.1 (2022-09-07)

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

