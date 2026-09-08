"""Validate and publish the three public crates. Uses only Python's standard library."""

import argparse
import json
import os
import re
import subprocess
import tomllib
import urllib.error
import urllib.request
from pathlib import Path

CRATES = ("fhe-util", "fhe-math", "fhe")
ROOT = Path(__file__).resolve().parents[2]
SEMVER = r"(?:0|[1-9][0-9]*)\.(?:0|[1-9][0-9]*)\.(?:0|[1-9][0-9]*)"


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def validate_tag(tag):
    # Deliberately accept stable releases only; never interpret input as a git option.
    if not re.fullmatch("v" + SEMVER, tag):
        raise ValueError("Release tag must be vMAJOR.MINOR.PATCH (no prereleases).")
    commit = git("rev-parse", "--verify", f"refs/tags/{tag}^{{commit}}")
    subprocess.run(
        ["git", "merge-base", "--is-ancestor", commit, "refs/remotes/origin/main"],
        cwd=ROOT, check=True,
    )
    manifest = tomllib.loads(git("show", f"{commit}:crates/fhe/Cargo.toml"))
    if tag != "v" + manifest["package"]["version"]:
        raise ValueError("Release tag does not match the tagged fhe package version.")
    return commit


def packages():
    result = []
    for name in CRATES:
        with (ROOT / "crates" / name / "Cargo.toml").open("rb") as source:
            manifest = tomllib.load(source)
        package = manifest["package"]
        if package["name"] != name:
            raise ValueError(f"Unexpected package name in {name}")
        result.append({"name": name, "version": package["version"], "manifest": manifest})
    versions = {p["name"]: p["version"] for p in result}
    for package in result:
        for name, spec in package["manifest"].get("dependencies", {}).items():
            if name in versions and spec.get("version") != "=" + versions[name]:
                raise ValueError(f"{package['name']} must depend on {name} ={versions[name]}")
    return result


def published(name, version):
    # Read the sparse index: missing crate is distinct from a registry/network failure.
    request = urllib.request.Request(
        "https://index.crates.io/" + (f"3/{name[0]}/{name}" if len(name) == 3
                                    else f"{name[:2]}/{name[2:4]}/{name}"),
        headers={"User-Agent": "fhe.rs-release-workflow"},
    )
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            entries = [json.loads(line) for line in response.read().decode().splitlines()]
    except urllib.error.HTTPError as error:
        if error.code == 404:
            return False
        raise
    for entry in entries:
        if entry["vers"] == version:
            if entry["yanked"]:
                raise ValueError(f"{name} {version} is yanked; refusing to skip it")
            return True
    return False


def plan():
    pending = []
    for package in packages():
        name, version = package["name"], package["version"]
        exists = published(name, version)
        print(f"{'Already published' if exists else 'Will publish'}: {name} {version}", flush=True)
        if not exists:
            pending.append(name)
    return pending


def cargo_publish(names, dry_run):
    if not names:
        print("All manifest versions are already published.", flush=True)
        return
    command = ["cargo", "+stable", "publish", "--locked", "--registry", "crates-io"]
    for name in names:
        command.extend(["--package", name])
    if dry_run:
        command.append("--dry-run")
    subprocess.run(command, cwd=ROOT, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("validate-tag", "dry-run", "publish"))
    args = parser.parse_args()
    if args.mode == "validate-tag":
        commit = validate_tag(os.environ["RELEASE_TAG"])
        with open(os.environ["GITHUB_OUTPUT"], "a") as output:
            output.write(f"commit={commit}\n")
    elif args.mode == "dry-run":
        # A single invocation lets Cargo verify unpublished workspace dependencies.
        cargo_publish(plan(), dry_run=True)
    else:
        # Sequential publishing permits a rerun to resume after a partial release.
        for name in plan():
            cargo_publish([name], dry_run=False)


if __name__ == "__main__":
    main()
