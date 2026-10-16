# GitHub Actions maintenance

The workflows run on pull requests targeting `main`, pushes to `main`, and
manual dispatch. They deliberately have no path filters so documentation-only
pull requests also complete any required checks. Existing job names are retained
for branch protection; MSRV and workflow lint checks can additionally be required
in repository settings.

- `rust.yml`: stable checks for every target and tests (including documentation
  tests), with default and all features. Rust 1.91.1 also checks every target in
  both feature configurations. Keep the MSRV aligned with `Cargo.toml` and the
  root README.
- `lint-fmt.yml`: pinned nightly rustfmt, stable Clippy for every target with default and
  all features, and actionlint to validate the workflows themselves. Nightly is
  required for the `wrap_comments` option in `rustfmt.toml`.
  CI uses `nightly-2026-08-15` to avoid formatting drift as nightly changes. Use
  that dated toolchain locally for matching results, and review formatting when
  updating the date. actionlint uses a prebuilt release binary whose SHA256 is
  verified before extraction; Go is not needed. Update its versioned download
  URL and checksum together in `lint-fmt.yml`. Dependabot does not update this
  inline binary download.
- `security.yml`: cargo-audit scans the committed `Cargo.lock`, including a weekly
  scheduled run to detect new advisories without a code change. Findings fail the
  job directly; the workflow needs no secrets or permission to create issues.
- `coverage.yml`: combines default- and all-feature coverage using stable
  cargo-llvm-cov, saves an LCOV artifact for 14 days, and uploads to Codecov.
- `typos.yml`: spelling checks using the root `typos.toml` configuration.

Cargo builds use `--locked` so dependency changes must include an updated
`Cargo.lock`. Rust caches distinguish jobs, toolchains, and feature configurations;
only pushes to `main` save caches. Superseded pull request runs are cancelled and
all jobs have timeouts.

All actions are pinned to full commit hashes. Dependabot checks them daily and
groups action updates into one pull request. Keep version comments alongside the
hashes. The Rust toolchain action must use a commit from its `master` history with
an explicit `toolchain` input, as required by its upstream documentation.

To validate workflow edits locally, download the binary for your platform from
the [actionlint v1.7.12 release](https://github.com/rhysd/actionlint/releases/tag/v1.7.12),
verify it against the release's checksums, and put it on `PATH`. Run `actionlint`
from the repository root. CI uses the Linux x86_64 archive for `ubuntu-latest`.
Install the formatter with
`rustup toolchain install nightly-2026-08-15 --profile minimal --component rustfmt`.
To reproduce Rust checks (with `protoc` installed):

```sh
cargo +nightly-2026-08-15 fmt --all -- --check
cargo +stable check --workspace --all-targets --locked
cargo +stable check --workspace --all-targets --locked --all-features
cargo +stable test --workspace --locked
cargo +stable test --workspace --locked --all-features
cargo +stable clippy --workspace --all-targets --locked -- -D warnings
cargo +stable clippy --workspace --all-targets --locked --all-features -- -D warnings
cargo +1.91.1 check --workspace --all-targets --locked
cargo +1.91.1 check --workspace --all-targets --locked --all-features
cargo audit --file Cargo.lock
```

The validation workflows do not publish crates. `publish.yml` provides the
manual release process described below; no workflow deploys a website.

## Coverage setup

Ensure `tlepoint/fhe.rs` is activated in Codecov and its GitHub App has access to
this repository. Uploads authenticate with GitHub OIDC, so no `CODECOV_TOKEN`
secret is required. Only the coverage job requests `id-token: write`. Public fork
PRs use the action's tokenless upload support. Dependabot runs generate and save
coverage but skip the external upload because their token is read-only.

Open a PR to run the workflow automatically. Once the workflow is on `main`, it
can also be run manually as **Coverage** in the Actions tab. Check the Codecov
upload step and the repository's Codecov page. The first
successful upload on `main` refreshes the README badge. Its public URL does not
need a badge token. Upload failures fail the job; the LCOV artifact is saved first
so it remains available if Codecov is unavailable. No minimum coverage threshold
has been imposed before establishing a baseline.

The report excludes generated protobuf code, standalone tests, benchmarks,
examples, and build scripts. Stable coverage measures unit and integration test
execution; documentation tests still run in `rust.yml` but are not instrumented
here (cargo-llvm-cov requires nightly for documentation-test coverage).

To reproduce the combined report locally, install cargo-llvm-cov and run:

```sh
rustup component add llvm-tools-preview --toolchain stable
cargo +stable llvm-cov clean --workspace
cargo +stable llvm-cov --workspace --locked --no-report
cargo +stable llvm-cov --workspace --locked --all-features --no-report
mkdir -p target/coverage
cargo +stable llvm-cov report --lcov --output-path target/coverage/lcov.info --ignore-filename-regex '/(proto|tests|benches|examples)/|/build\.rs$'
```

See the [cargo-llvm-cov documentation](https://github.com/taiki-e/cargo-llvm-cov)
and [Codecov OIDC setup](https://github.com/codecov/codecov-action#using-oidc).


## Publishing to crates.io

`publish.yml` is manually dispatched from **main** and defaults to a dry run.
It accepts an existing stable tag `vMAJOR.MINOR.PATCH`, whose version must match
`crates/fhe/Cargo.toml`. The tag must resolve to a commit reachable from `main`.
Both jobs check out the validated commit SHA, so moving the tag during a run
cannot change the code being published. The tagged commit must include the
publishing workflow and scripts; merge this setup before creating a release tag.
Prerelease tags are intentionally not supported.

### One-time account setup

1. Create a GitHub environment named `crates-io`. Restrict its deployment branches
   to `main` (the workflow is dispatched from main, then checks out the release
   commit). Configure required reviewers if you want an approval before upload.
2. In the crates.io settings for **each** of `fhe-util`, `fhe-math`, and `fhe`, add
   a GitHub Trusted Publisher with owner `tlepoint`, repository `fhe.rs`, workflow
   filename `publish.yml`, and environment `crates-io`. Existing crate ownership
   is required to configure this.
3. No `CARGO_REGISTRY_TOKEN` repository secret is needed. Only the publishing job
   receives `id-token: write`; the official Rust authentication action exchanges
   GitHub OIDC credentials for a temporary token immediately before publishing.

See [Trusted Publishing](https://crates.io/docs/trusted-publishing) and the
[official authentication action](https://github.com/rust-lang/crates-io-auth-action).

### Each release

1. Prepare and merge a release PR. Bump versions only for crates being released,
   update exact internal dependency versions to match the workspace manifests,
   refresh `Cargo.lock`, and update release notes. A utility version bump requires
   updating and releasing its dependents. `fhe-util` may retain a different
   version from `fhe-math` and `fhe`.
2. Create and push a tag matching the new `fhe` version on the reviewed main
   commit, for example `v0.3.0`. Tags alone do not trigger publication.
3. In Actions, select **Publish crates → Run workflow**, use branch **main**,
   enter the tag, and leave **dry_run** checked. Review the package plan and
   verification results.
4. Run again with the same tag and **dry_run** unchecked to publish. Any configured
   environment approval occurs after validation. A successful dry run does not
   upload anything or require crates.io credentials.

Validation tests the tagged commit with default and all features, checks rustfmt
and Clippy, tests the release script, and runs `cargo publish --dry-run --locked`
for all pending packages together. Modern stable Cargo can verify unpublished
workspace dependencies in that single invocation. The publishing job repeats
package verification before authentication, then publishes `fhe-util`, `fhe-math`,
and `fhe` sequentially with normal Cargo build verification enabled.

The script queries the crates.io sparse index for exact manifest versions. It
skips versions already published, rejects yanked versions, and fails on registry
errors. Rerun the same tag after a partial failure; successful uploads cannot be
rolled back, and Cargo may report a timeout after an upload has succeeded. Check
the registry before retrying. Skipping an existing version checks version
availability, not byte-for-byte equality with the local source; unchanged
versions must represent the already released code. If every version exists, the
run is a no-op. Release runs are serialized and never cancel an active upload.

Run the release-tool regression tests locally with:

```sh
python3 -m unittest discover -s .github/scripts -p 'test_*.py'
```

Python 3.11 or newer is required (`tomllib`); the Ubuntu runner provides it.
