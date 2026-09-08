"""Release safeguards, registry errors, and partial-publication regression tests."""

import io
import json
import subprocess
import unittest
import urllib.error
from unittest.mock import patch

import publish


class ReleaseTests(unittest.TestCase):
    def test_rejects_invalid_tags_before_running_git(self):
        for tag in ("main", "--help", "v1.2", "v01.2.3", "v1.2.3-rc.1", "v1.2.3\n"):
            with self.subTest(tag=tag), patch.object(publish, "git") as git:
                with self.assertRaises(ValueError):
                    publish.validate_tag(tag)
                git.assert_not_called()

    def test_tag_resolves_commit_and_checks_main_ancestry(self):
        with patch.object(publish, "git", side_effect=["abc123", '[package]\nversion="1.2.3"']), \
                patch.object(publish.subprocess, "run") as run:
            self.assertEqual(publish.validate_tag("v1.2.3"), "abc123")
            self.assertEqual(run.call_args.args[0],
                             ["git", "merge-base", "--is-ancestor", "abc123", "refs/remotes/origin/main"])

    def test_rejects_tag_version_mismatch(self):
        with patch.object(publish, "git", side_effect=["abc123", '[package]\nversion="1.2.4"']), \
                patch.object(publish.subprocess, "run"):
            with self.assertRaises(ValueError):
                publish.validate_tag("v1.2.3")

    def test_rejects_commit_outside_main(self):
        with patch.object(publish, "git", return_value="abc123"), \
                patch.object(publish.subprocess, "run", side_effect=subprocess.CalledProcessError(1, "git")):
            with self.assertRaises(subprocess.CalledProcessError):
                publish.validate_tag("v1.2.3")

    def test_workspace_dependency_versions_match(self):
        self.assertEqual([p["name"] for p in publish.packages()], list(publish.CRATES))

    def test_rejects_stale_internal_dependency_versions(self):
        manifests = [p["manifest"] for p in publish.packages()]
        manifests[1]["dependencies"]["fhe-util"]["version"] = "=0.0.0"
        with patch.object(publish.tomllib, "load", side_effect=manifests):
            with self.assertRaisesRegex(ValueError, "must depend on fhe-util"):
                publish.packages()

    def test_registry_exact_version_and_sparse_paths(self):
        for name, path in (("fhe", "3/f/fhe"), ("fhe-util", "fh/e-/fhe-util")):
            data = json.dumps({"vers": "1.2.3", "yanked": False}).encode()
            with patch.object(publish.urllib.request, "urlopen", return_value=io.BytesIO(data)) as fetch:
                self.assertTrue(publish.published(name, "1.2.3"))
                self.assertEqual(fetch.call_args.args[0].full_url, "https://index.crates.io/" + path)
            with patch.object(publish.urllib.request, "urlopen", return_value=io.BytesIO(data)):
                self.assertFalse(publish.published(name, "1.2.4"))

    def test_missing_crate_is_unpublished(self):
        error = urllib.error.HTTPError("url", 404, "missing", {}, io.BytesIO())
        self.addCleanup(error.close)
        with patch.object(publish.urllib.request, "urlopen", side_effect=error):
            self.assertFalse(publish.published("fhe", "1.2.3"))

    def test_registry_errors_fail_closed(self):
        for error in (urllib.error.HTTPError("url", 503, "down", {}, io.BytesIO()),
                      urllib.error.URLError("offline")):
            if isinstance(error, urllib.error.HTTPError):
                self.addCleanup(error.close)
            with patch.object(publish.urllib.request, "urlopen", side_effect=error):
                with self.assertRaises(urllib.error.URLError):
                    publish.published("fhe", "1.2.3")

    def test_yanked_version_is_not_skipped(self):
        data = b'{"vers":"1.2.3","yanked":true}'
        with patch.object(publish.urllib.request, "urlopen", return_value=io.BytesIO(data)):
            with self.assertRaises(ValueError):
                publish.published("fhe", "1.2.3")

    def test_partial_release_keeps_dependency_order(self):
        with patch.object(publish, "published", side_effect=[True, False, False]):
            self.assertEqual(publish.plan(), ["fhe-math", "fhe"])

    def test_dry_run_verifies_pending_packages_together(self):
        with patch.object(publish.subprocess, "run") as run:
            publish.cargo_publish(["fhe-math", "fhe"], dry_run=True)
            self.assertEqual(run.call_args.args[0], [
                "cargo", "+stable", "publish", "--locked", "--registry", "crates-io",
                "--package", "fhe-math", "--package", "fhe", "--dry-run",
            ])

    def test_completed_release_does_not_publish(self):
        with patch.object(publish.subprocess, "run") as run:
            publish.cargo_publish([], dry_run=False)
            run.assert_not_called()

    def test_publish_stops_on_first_failure(self):
        with patch("sys.argv", ["publish.py", "publish"]), \
                patch.object(publish, "plan", return_value=list(publish.CRATES)), \
                patch.object(publish, "cargo_publish", side_effect=subprocess.CalledProcessError(1, "cargo")) as cargo:
            with self.assertRaises(subprocess.CalledProcessError):
                publish.main()
            cargo.assert_called_once_with(["fhe-util"], dry_run=False)


if __name__ == "__main__":
    unittest.main()
