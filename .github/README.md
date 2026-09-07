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

These workflows validate the library; they do not publish crates or deploy a
website.

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
