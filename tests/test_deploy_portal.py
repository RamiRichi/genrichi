"""
Tests for deploy_portal.sh (Phase 2 hardening).

The script is exercised against a synthetic world under /tmp: a throw-away Git
repository, a fake "production" home directory, fake systemctl/pgrep commands
and a fake cgroup tree. Production is never touched. Test mode is switched on
with GENRICHI_DEPLOY_TEST_* variables, which the script only honours for paths
under /tmp.

Also checks that the tracked portal files on the deployment allowlist contain no
hardcoded credential material (the publicly exposed historical values are
derived from Git history at test time, so this file contains no secret).

Run with (POSIX/WSL):
    python -m unittest tests.test_deploy_portal -v
"""

import os
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "deploy_portal.sh"
sys.path.insert(0, str(Path(__file__).resolve().parent))

from test_launchers_safe import scan_text  # noqa: E402

POSIX = os.name == "posix" and shutil.which("bash") and shutil.which("git") and shutil.which("sha256sum")

EXPECTED_ALLOWLIST = [
    "app.py", "config.py", "mailer.py", "models.py", "runner.py",
    "templates/base.html", "templates/dashboard.html", "templates/invoice.html", "templates/login.html",
    "templates/new_order.html", "templates/order.html", "templates/report_view.html",
    "templates/settings.html", "templates/stats.html", "templates/users.html",
    "static/img/logo.png",
]
PNG_BYTES = b"\x89PNG\r\n\x1a\n" + bytes(range(64))
SENTINEL_ENV = b"SENTINEL-not-a-real-secret\n"
SENTINEL_DB = b"SQLite format 3\x00" + b"sentinel" * 20

SYSTEMCTL_SHIM = r"""#!/bin/sh
echo "$*" >> "$SHIM_STATE/systemctl.log"
state() { if [ -f "$SHIM_STATE/$1" ]; then cat "$SHIM_STATE/$1"; else echo "$2"; fi; }
case "$1 $2" in
  "is-active genrichi.service")         s=$(state dup_state inactive); echo "$s"; [ "$s" = active ] || exit 3 ;;
  "is-active genrichi-portal.service")  s=$(state portal_state active); echo "$s"; [ "$s" = active ] || exit 3 ;;
  "is-enabled genrichi.service")        echo disabled; exit 1 ;;
  "show genrichi-portal.service")       echo /fake.slice/genrichi-portal.service ;;
  "cat genrichi-portal.service")        cat "$SHIM_STATE/unit.txt" ;;
  *) echo "unexpected systemctl call: $*" >&2; echo "$*" >> "$SHIM_STATE/forbidden.log"; exit 99 ;;
esac
"""
PGREP_SHIM = r"""#!/bin/sh
if [ -f "$SHIM_STATE/pgrep_error" ]; then echo "pgrep: simulated failure" >&2; exit "$(cat "$SHIM_STATE/pgrep_error")"; fi
if [ -f "$SHIM_STATE/pgrep_garbage" ]; then echo "not-a-number"; exit 0; fi
n=0; [ -f "$SHIM_STATE/snakemake_count" ] && n=$(cat "$SHIM_STATE/snakemake_count")
echo "$n"; [ "$n" -gt 0 ]
"""
FORBIDDEN_SHIM = r"""#!/bin/sh
echo "$0 $*" >> "$SHIM_STATE/forbidden.log"
exit 99
"""


def git(repo, *args, check=True):
    env = dict(os.environ, GIT_CONFIG_GLOBAL="/dev/null", GIT_CONFIG_SYSTEM="/dev/null",
               GIT_AUTHOR_NAME="t", GIT_AUTHOR_EMAIL="t@example.invalid",
               GIT_COMMITTER_NAME="t", GIT_COMMITTER_EMAIL="t@example.invalid")
    return subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True, env=env, check=check).stdout.strip()


def sha256(path):
    import hashlib
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def content_for(rel, version):
    if rel.endswith(".png"):
        return PNG_BYTES + version.encode()
    if rel.endswith(".py"):
        return f"# {rel} {version}\nVALUE = '{version}'\n".encode()
    return f"<p>{rel} {version}</p>\n".encode()


@unittest.skipUnless(POSIX, "needs a POSIX shell with git and sha256sum (run in WSL/Linux)")
class DeployTestCase(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp(prefix="deploytest-", dir="/tmp"))
        self.addCleanup(shutil.rmtree, self.root, ignore_errors=True)
        self.home = self.root / "home"
        self.repo = self.root / "repo"
        self.state = self.root / "state"
        self.shims = self.root / "shims"
        self.cgroup = self.root / "cgroup"
        (self.root / "tmp").mkdir()
        for d in (self.home / "portal", self.state, self.shims):
            d.mkdir(parents=True)
        self.dst = self.home / "portal"
        self._build_repo()
        self._build_production()
        self._build_shims()
        procs = self.cgroup / "fake.slice" / "genrichi-portal.service"
        procs.mkdir(parents=True)
        (procs / "cgroup.procs").write_text("4242\n")
        (self.state / "unit.txt").write_text(
            "[Service]\nUser=rami\nWorkingDirectory=%s\nExecStart=%s app.py\nRestart=on-failure\n[Install]\nWantedBy=multi-user.target\n"
            % (self.dst, self.pywrap))

    # ── fixtures ────────────────────────────────────────────────────────
    def _build_repo(self):
        self.repo.mkdir()
        git(self.repo, "init", "-q")
        for rel in EXPECTED_ALLOWLIST:
            p = self.repo / "portal" / rel
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_bytes(content_for(rel, "v1"))
        (self.repo / "portal" / ".env.example").write_text("PLACEHOLDER=<set-me>\n")
        self.commit("initial")

    def commit(self, message):
        git(self.repo, "add", "-A")
        git(self.repo, "commit", "-q", "-m", message)
        self.head = git(self.repo, "rev-parse", "HEAD")

    def _build_production(self):
        # older content for some files, identical for others, some missing (NEW)
        for rel in EXPECTED_ALLOWLIST:
            if rel in ("mailer.py", "templates/users.html") or rel.startswith("static/"):
                continue  # NEW files (static/img directory does not exist yet)
            p = self.dst / rel
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_bytes(content_for(rel, "v0" if rel in ("app.py", "models.py", "templates/base.html") else "v1"))
        (self.dst / ".env").write_bytes(SENTINEL_ENV)
        os.chmod(self.dst / ".env", 0o600)
        (self.dst / "genrichi_orders.db").write_bytes(SENTINEL_DB)
        (self.dst / "uploads").mkdir()
        (self.dst / "uploads" / "sample.fastq").write_text("uploaded data\n")
        (self.dst / "portal_logs").mkdir()
        (self.dst / "portal_logs" / "x.log").write_text("log\n")
        (self.home / "results").mkdir()
        (self.home / "results" / "r.txt").write_text("result\n")
        (self.home / ".snakemake").mkdir()
        (self.home / ".snakemake" / "meta").write_text("state\n")

    def _build_shims(self):
        for name, body in (("systemctl", SYSTEMCTL_SHIM), ("pgrep", PGREP_SHIM)):
            (self.shims / name).write_text(body)
        for name in ("sudo", "pip", "pip3", "nohup", "pkill", "killall", "cloudflared"):
            (self.shims / name).write_text(FORBIDDEN_SHIM)
        for p in self.shims.iterdir():
            p.chmod(0o755)
        # A python wrapper that records its arguments, so tests can prove every
        # validation subprocess runs with -B (no .pyc may be written).
        self.pywrap = self.root / "python-wrapper"
        self.pywrap.write_text('#!/bin/sh\necho "$*" >> "$SHIM_STATE/python.log"\nexec "%s" "$@"\n' % sys.executable)
        self.pywrap.chmod(0o755)

    def path_without(self, *names):
        """A PATH holding the fake systemctl/pgrep shims plus every system command EXCEPT `names`."""
        d = self.root / "path-without"
        d.mkdir(exist_ok=True)
        for src in ("/usr/local/bin", "/usr/bin", "/bin"):
            if not os.path.isdir(src):
                continue
            for entry in os.listdir(src):
                link = d / entry
                if entry in names or link.exists() or link.is_symlink():
                    continue
                try:
                    os.symlink(os.path.join(src, entry), link)
                except OSError:
                    pass
        for shim in self.shims.iterdir():
            if shim.name not in names:
                target = d / shim.name
                if target.is_symlink() or target.exists():
                    target.unlink()
                os.symlink(shim, target)
        return str(d)

    # ── running ─────────────────────────────────────────────────────────
    def run_script(self, *args, modules="json sqlite3", extra_env=None, path_override=None):
        env = dict(os.environ)
        env.update({
            "PATH": path_override or f"{self.shims}:{os.environ['PATH']}",
            "SHIM_STATE": str(self.state),
            "TMPDIR": str(self.root / "tmp"),
            "GENRICHI_DEPLOY_TEST_HOME": str(self.home),
            "GENRICHI_DEPLOY_TEST_REPO": str(self.repo),
            "GENRICHI_DEPLOY_TEST_PYTHON": str(self.pywrap),
            "GENRICHI_DEPLOY_TEST_CGROUP_ROOT": str(self.cgroup),
            "GENRICHI_DEPLOY_TEST_MODULES": modules,
            "GIT_CONFIG_GLOBAL": "/dev/null",
        })
        env.update(extra_env or {})
        result = subprocess.run(["bash", str(SCRIPT), *args], capture_output=True, text=True, env=env, timeout=120)
        self.output = result.stdout + result.stderr
        self.assertFalse((self.state / "forbidden.log").exists(),
                         "the script invoked a forbidden command: " + (self.state / "forbidden.log").read_text()
                         if (self.state / "forbidden.log").exists() else "")
        return result

    def snapshot(self, include_allowlisted=True):
        """path -> (sha256, mtime_ns, mode) for every file under home."""
        allow = {str(self.dst / rel) for rel in EXPECTED_ALLOWLIST}
        snap = {}
        for base, dirs, files in os.walk(self.home):
            for name in files + dirs:
                p = Path(base) / name
                if not include_allowlisted and str(p) in allow:
                    continue
                st = p.lstat()
                snap[str(p)] = (sha256(p) if p.is_file() else "dir", st.st_mtime_ns if p.is_file() else 0, stat.S_IMODE(st.st_mode))
        return snap

    def backup_dirs(self):
        parent = self.home / "deploy_backups" / "portal"
        return sorted(parent.iterdir()) if parent.exists() else []

    def assert_nothing_written(self, before):
        self.assertEqual(self.snapshot(), before, "production files changed")
        self.assertEqual(self.backup_dirs(), [], "a backup directory was created")
        self.assertEqual(list(self.dst.rglob(".deploy.*")), [], "temporary files left behind")

    def assert_refused(self, result, needle=None):
        self.assertEqual(result.returncode, 1, self.output)
        if needle:
            self.assertIn(needle, self.output)

    # ── dry-run / basics ────────────────────────────────────────────────
    def test_bash_syntax_is_valid(self):
        self.assertEqual(subprocess.run(["bash", "-n", str(SCRIPT)], capture_output=True).returncode, 0)

    def test_help_exits_zero_and_writes_nothing(self):
        before = self.snapshot()
        result = self.run_script("--help")
        self.assertEqual(result.returncode, 0)
        self.assertIn("dry-run", result.stdout.lower())
        self.assert_nothing_written(before)

    def test_unknown_option_is_a_usage_error(self):
        self.assertEqual(self.run_script("--frobnicate").returncode, 2)

    def test_default_is_dry_run_and_writes_nothing(self):
        before = self.snapshot()
        result = self.run_script()
        self.assertEqual(result.returncode, 0, self.output)
        self.assertIn("DRY-RUN", self.output)
        self.assertIn("nothing was written", self.output)
        self.assertIn(self.head, self.output)
        self.assert_nothing_written(before)

    def test_dry_run_reports_the_plan(self):
        self.run_script()
        self.assertIn("CHANGED", self.output)
        self.assertIn("NEW", self.output)
        self.assertRegex(self.output, r"3 changed, 3 new, 10 unchanged, 16 allowlisted|\d+ changed, \d+ new, \d+ unchanged, 16 allowlisted")

    def test_dry_run_exits_nonzero_when_execute_would_be_refused(self):
        (self.state / "dup_state").write_text("active")
        before = self.snapshot()
        result = self.run_script()
        self.assert_refused(result, "would be refused")
        self.assert_nothing_written(before)

    # ── execute: what gets written ──────────────────────────────────────
    def test_execute_copies_exactly_the_allowlist_from_head(self):
        before_protected = self.snapshot(include_allowlisted=False)
        files_before = {str(p) for p in self.dst.rglob("*") if p.is_file()}
        result = self.run_script("--execute", "--expect-commit", self.head)
        self.assertEqual(result.returncode, 0, self.output)
        for rel in EXPECTED_ALLOWLIST:
            self.assertEqual((self.dst / rel).read_bytes(), content_for(rel, "v1"), rel)
        new_files = {str(p) for p in self.dst.rglob("*") if p.is_file()} - files_before
        self.assertEqual(new_files, {str(self.dst / r) for r in ("mailer.py", "templates/users.html", "static/img/logo.png")})
        self.assertFalse((self.dst / ".env.example").exists(), ".env.example must not be deployed")
        allow = {str(self.dst / rel) for rel in EXPECTED_ALLOWLIST}
        created = set(self.snapshot(include_allowlisted=False)) - set(before_protected)
        stray = {p for p in created if p not in allow and "deploy_backups" not in p}
        self.assertEqual(stray, {str(self.dst / "static"), str(self.dst / "static" / "img")}, "only the new static dirs may appear")

    def test_deployed_files_are_mode_644_and_owned_by_the_operator(self):
        self.assertEqual(self.run_script("--execute").returncode, 0, self.output)
        for rel in EXPECTED_ALLOWLIST:
            st = (self.dst / rel).stat()
            self.assertEqual(stat.S_IMODE(st.st_mode), 0o644, rel)
            self.assertEqual(st.st_uid, os.getuid(), rel)

    def test_env_db_uploads_results_and_snakemake_are_untouched(self):
        protected = self.snapshot(include_allowlisted=False)
        self.assertEqual(self.run_script("--execute").returncode, 0, self.output)
        after = self.snapshot(include_allowlisted=False)
        for path, meta in protected.items():
            self.assertEqual(after.get(path), meta, f"{path} was modified")
        self.assertEqual((self.dst / ".env").read_bytes(), SENTINEL_ENV)
        self.assertEqual(stat.S_IMODE((self.dst / ".env").stat().st_mode), 0o600)
        self.assertEqual((self.dst / "genrichi_orders.db").read_bytes(), SENTINEL_DB)

    def test_no_temporary_files_remain_after_a_successful_deploy(self):
        self.assertEqual(self.run_script("--execute").returncode, 0, self.output)
        self.assertEqual(list(self.dst.rglob(".deploy.*")), [])

    def test_second_run_is_a_noop_without_a_new_backup(self):
        self.assertEqual(self.run_script("--execute").returncode, 0, self.output)
        backups = self.backup_dirs()
        result = self.run_script("--execute")
        self.assertEqual(result.returncode, 0, self.output)
        self.assertIn("Nothing to deploy", self.output)
        self.assertEqual(self.backup_dirs(), backups)

    # ── backups ─────────────────────────────────────────────────────────
    def test_backup_contains_only_overwritten_files_byte_for_byte(self):
        old = {rel: (self.dst / rel).read_bytes() for rel in ("app.py", "models.py", "templates/base.html")}
        self.assertEqual(self.run_script("--execute").returncode, 0, self.output)
        backups = self.backup_dirs()
        self.assertEqual(len(backups), 1)
        b = backups[0]
        self.assertEqual(stat.S_IMODE(b.stat().st_mode), 0o700)
        self.assertEqual(stat.S_IMODE((self.home / "deploy_backups").stat().st_mode), 0o700)
        saved = {str(p.relative_to(b / "files")): p.read_bytes() for p in (b / "files").rglob("*") if p.is_file()}
        self.assertEqual(saved, old)                                   # only the 3 CHANGED files, unchanged bytes
        for name in (".env", "genrichi_orders.db", "sample.fastq", "r.txt", "meta", "x.log"):
            self.assertEqual(list(b.rglob(name)), [], f"{name} must never be backed up")

    def test_backup_manifest_has_commit_and_sha256(self):
        old_sha = sha256(self.dst / "app.py")
        self.assertEqual(self.run_script("--execute").returncode, 0, self.output)
        manifest = self.backup_dirs()[0] / "MANIFEST.txt"
        self.assertEqual(stat.S_IMODE(manifest.stat().st_mode), 0o600)
        text = manifest.read_text()
        self.assertIn(f"source_commit={self.head}", text)
        self.assertIn("source_mode=head", text)
        self.assertIn(f"source_bytes=exact blobs of commit {self.head}", text)
        new_sha = sha256(self.dst / "app.py")
        self.assertIn(f"CHANGED app.py {old_sha} {new_sha}", text)
        self.assertRegex(text, r"NEW mailer\.py - [0-9a-f]{64}")

    def test_backup_failure_aborts_before_any_production_write(self):
        (self.home / "deploy_backups").write_text("this is a file, not a directory")
        before = {k: v for k, v in self.snapshot().items() if "deploy_backups" not in k}
        result = self.run_script("--execute")
        self.assertNotEqual(result.returncode, 0)
        after = {k: v for k, v in self.snapshot().items() if "deploy_backups" not in k}
        self.assertEqual(after, before, "production changed although the backup failed")

    # ── refusals ────────────────────────────────────────────────────────
    def refuse_case(self, prepare, needle, *args):
        prepare()
        before = self.snapshot()
        result = self.run_script("--execute", *args)
        self.assert_refused(result, needle)
        self.assert_nothing_written(before)

    def test_missing_env_is_refused(self):
        self.refuse_case(lambda: (self.dst / ".env").unlink(), ".env")

    def test_env_wider_than_640_is_refused(self):
        self.refuse_case(lambda: os.chmod(self.dst / ".env", 0o644), "640")

    def test_active_genrichi_service_is_refused(self):
        self.refuse_case(lambda: (self.state / "dup_state").write_text("active"), "genrichi.service is active")

    def test_activating_genrichi_service_is_refused(self):
        self.refuse_case(lambda: (self.state / "dup_state").write_text("activating"), "activating")

    def test_running_job_seen_in_the_portal_cgroup_is_refused(self):
        procs = self.cgroup / "fake.slice" / "genrichi-portal.service" / "cgroup.procs"
        self.refuse_case(lambda: procs.write_text("4242\n5151\n"), "analysis appears to be running")

    def test_running_snakemake_process_is_refused(self):
        self.refuse_case(lambda: (self.state / "snakemake_count").write_text("1"), "Snakemake process")

    def test_snakemake_locks_directory_is_refused(self):
        self.refuse_case(lambda: (self.home / ".snakemake" / "locks").mkdir(), "locks")

    def test_unreadable_cgroup_fails_closed(self):
        self.refuse_case(lambda: (self.cgroup / "fake.slice" / "genrichi-portal.service" / "cgroup.procs").unlink(), "cgroup")

    def test_unknown_portal_state_fails_closed(self):
        self.refuse_case(lambda: (self.state / "portal_state").write_text("weird"), "cannot determine")

    def test_unknown_tracked_portal_file_is_refused(self):
        def prepare():
            (self.repo / "portal" / "evil.py").write_text("x = 1\n")
            self.commit("add unlisted file")
        self.refuse_case(prepare, "unknown tracked file under portal/: evil.py")

    def test_tracked_dot_env_is_refused_as_unknown(self):
        def prepare():
            (self.repo / "portal" / ".env").write_text("A=B\n")
            self.commit("oops committed .env")
        self.refuse_case(prepare, "unknown tracked file under portal/: .env")

    def test_missing_allowlisted_file_is_refused(self):
        def prepare():
            git(self.repo, "rm", "-q", "portal/templates/users.html")
            self.commit("remove template")
        self.refuse_case(prepare, "missing in commit")

    def test_credential_like_literal_in_source_is_refused_and_not_echoed(self):
        secret = "Zq7-fake-not-real-1234"
        def prepare():
            (self.repo / "portal" / "app.py").write_text(f'# app\nSECRET_KEY = "{secret}"\n')
            self.commit("hardcode a secret")
        prepare()
        before = self.snapshot()
        result = self.run_script("--execute")
        self.assert_refused(result, "credential-like literal in app.py")
        self.assertNotIn(secret, self.output)
        self.assert_nothing_written(before)

    def test_crlf_in_source_is_refused(self):
        def prepare():
            (self.repo / "portal" / "app.py").write_bytes(b"# app\r\nX = 1\r\n")
            self.commit("crlf")
        self.refuse_case(prepare, "carriage returns")

    def test_python_syntax_error_in_source_is_refused(self):
        def prepare():
            (self.repo / "portal" / "runner.py").write_text("def broken(:\n")
            self.commit("syntax error")
        self.refuse_case(prepare, "syntax check failed")

    def test_missing_dependency_is_refused_and_nothing_is_installed(self):
        before = self.snapshot()
        result = self.run_script("--execute", modules="definitely_missing_module_xyz")
        self.assert_refused(result, "definitely_missing_module_xyz")
        self.assertIn("nothing is installed", self.output)
        self.assert_nothing_written(before)

    def test_destination_symlink_is_refused(self):
        other = self.root / "elsewhere"
        shutil.copytree(self.dst, other, symlinks=True)
        shutil.rmtree(self.dst)
        os.symlink(other, self.dst)
        before = {str(p): sha256(p) for p in other.rglob("*") if p.is_file()}
        result = self.run_script("--execute")
        self.assert_refused(result, "symlink")
        self.assertEqual({str(p): sha256(p) for p in other.rglob("*") if p.is_file()}, before)

    def test_expect_commit_mismatch_is_refused(self):
        self.refuse_case(lambda: None, "does not match --expect-commit", "--expect-commit", "0" * 40)

    def test_expect_commit_accepts_the_full_sha_in_any_case(self):
        self.assertEqual(self.run_script("--expect-commit", self.head.upper()).returncode, 0, self.output)
        self.assertIn("HEAD matches --expect-commit", self.output)

    def test_expect_commit_rejects_short_prefixes_and_bad_values_as_usage_errors(self):
        before = self.snapshot()
        for bad in (self.head[:1], self.head[:8], self.head[:12], self.head[:39], self.head + "0", "z" * 40, "not-a-sha"):
            with self.subTest(value=bad[:12]):
                result = self.run_script("--execute", "--expect-commit", bad)
                self.assertEqual(result.returncode, 2, self.output)
                self.assertNotIn("Backup", self.output)
        self.assert_nothing_written(before)

    def test_expect_commit_needs_a_value(self):
        self.assertEqual(self.run_script("--expect-commit").returncode, 2)
        self.assertEqual(self.run_script("--expect-commit", "").returncode, 2)

    def test_expect_commit_is_an_exact_match_not_a_prefix_match(self):
        other = self.head[:39] + ("0" if self.head[39] != "0" else "1")
        before = self.snapshot()
        result = self.run_script("--execute", "--expect-commit", other)
        self.assert_refused(result, "does not match --expect-commit")
        self.assert_nothing_written(before)

    # ── source selection: HEAD by default, working tree only on request ─
    def dirty_app(self):
        (self.repo / "portal" / "app.py").write_text("# app dirty\nVALUE = 'uncommitted'\n")

    def test_head_mode_deploys_the_committed_blob_not_the_dirty_worktree(self):
        self.dirty_app()
        self.assertEqual(self.run_script("--execute").returncode, 0, self.output)
        self.assertEqual((self.dst / "app.py").read_bytes(), content_for("app.py", "v1"))
        self.assertIn("IGNORED", self.output)

    def test_worktree_mode_refuses_a_dirty_scope_without_allow_dirty(self):
        self.dirty_app()
        before = self.snapshot()
        result = self.run_script("--execute", "--from-worktree")
        self.assert_refused(result, "--allow-dirty")
        self.assert_nothing_written(before)

    def test_worktree_mode_with_allow_dirty_deploys_working_tree_content(self):
        self.dirty_app()
        result = self.run_script("--execute", "--from-worktree", "--allow-dirty")
        self.assertEqual(result.returncode, 0, self.output)
        self.assertIn(b"uncommitted", (self.dst / "app.py").read_bytes())
        manifest = (self.backup_dirs()[0] / "MANIFEST.txt").read_text()
        self.assertIn("source_mode=worktree", manifest)
        self.assertIn("WORKING-TREE files (not guaranteed to equal any commit)", manifest)
        self.assertIn("source_commit=none", manifest)
        self.assertIn(f"worktree_base_head={self.head}", manifest)
        self.assertIn("worktree_uncommitted_entries=1", manifest)
        self.assertNotIn(f"source_commit={self.head}", manifest)
        self.assertNotIn("source_mode=head", manifest)
        self.assertIn("from the WORKING TREE", self.output)
        self.assertNotIn(f"file(s) from commit {self.head}", self.output)

    def test_worktree_mode_on_a_clean_tree_needs_no_allow_dirty(self):
        self.assertEqual(self.run_script("--from-worktree").returncode, 0, self.output)

    def test_dirty_check_is_limited_to_the_deployment_scope(self):
        (self.repo / "README.md").write_text("changed outside portal/\n")
        self.assertEqual(self.run_script("--from-worktree").returncode, 0, self.output)

    # ── services and unit ───────────────────────────────────────────────
    def test_only_read_only_systemctl_verbs_are_used(self):
        self.assertEqual(self.run_script("--execute").returncode, 0, self.output)
        verbs = {line.split()[0] for line in (self.state / "systemctl.log").read_text().splitlines() if line.strip()}
        self.assertTrue(verbs <= {"is-active", "is-enabled", "show", "cat"}, verbs)

    def test_unit_drift_is_reported_but_never_fatal_or_modified(self):
        (self.state / "unit.txt").write_text("[Service]\nUser=someone-else\nRestart=always\n")
        result = self.run_script("--execute")
        self.assertEqual(result.returncode, 0, self.output)
        self.assertIn("UNIT DRIFT", self.output)
        self.assertEqual((self.state / "unit.txt").read_text(), "[Service]\nUser=someone-else\nRestart=always\n")

    def test_the_portal_is_not_restarted(self):
        self.assertEqual(self.run_script("--execute").returncode, 0, self.output)
        self.assertIn("NOT restarted", self.output)
        self.assertFalse(any(w in (self.state / "systemctl.log").read_text() for w in ("restart", "start", "stop", "reload")))

    # ── D3: job detection must fail closed when pgrep is unusable ───────
    def test_pgrep_missing_fails_closed(self):
        before = self.snapshot()
        result = self.run_script("--execute", path_override=self.path_without("pgrep"))
        self.assert_refused(result, "pgrep is not available")
        self.assertNotIn("ok    no Snakemake process is running", self.output)
        self.assert_nothing_written(before)

    def test_pgrep_error_exit_fails_closed(self):
        for rc in ("2", "3", "127"):
            with self.subTest(pgrep_exit=rc):
                (self.state / "pgrep_error").write_text(rc)
                before = self.snapshot()
                result = self.run_script("--execute")
                self.assert_refused(result, "pgrep failed")
                self.assertNotIn("ok    no Snakemake process is running", self.output)
                self.assert_nothing_written(before)

    def test_pgrep_unusable_output_fails_closed(self):
        (self.state / "pgrep_garbage").write_text("1")
        before = self.snapshot()
        result = self.run_script("--execute")
        self.assert_refused(result, "pgrep gave unusable output")
        self.assert_nothing_written(before)

    def test_pgrep_error_is_also_reported_by_the_dry_run(self):
        (self.state / "pgrep_error").write_text("2")
        result = self.run_script()
        self.assert_refused(result, "would be refused")

    def test_pgrep_no_match_is_the_only_idle_result(self):
        result = self.run_script()          # shim: prints 0, exits 1 like the real pgrep -c
        self.assertEqual(result.returncode, 0, self.output)
        self.assertIn("no Snakemake process is running", self.output)

    # ── D7: read-only git queries must not touch .git/index ─────────────
    def test_dry_run_and_execute_do_not_rewrite_the_git_index(self):
        # make every tracked file stat-dirty, so a plain `git status` WOULD refresh (rewrite) the index
        for path in (self.repo / "portal").rglob("*"):
            if path.is_file():
                os.utime(path, (path.stat().st_atime + 1000, path.stat().st_mtime + 1000))
        index = self.repo / ".git" / "index"
        before = (index.read_bytes(), index.stat().st_mtime_ns)
        self.assertEqual(self.run_script().returncode, 0, self.output)
        self.assertEqual((index.read_bytes(), index.stat().st_mtime_ns), before, "dry-run rewrote .git/index")
        self.assertEqual(self.run_script("--from-worktree").returncode, 0, self.output)
        self.assertEqual((index.read_bytes(), index.stat().st_mtime_ns), before, "worktree dry-run rewrote .git/index")

    def test_every_git_invocation_uses_no_optional_locks(self):
        text = SCRIPT.read_text(encoding="utf-8")
        calls = [m.group(1) for line in text.splitlines() if not line.strip().startswith("#")
                 for m in re.finditer(r"(?:^|[|&;({]|\$\()\s*git\s+(\S+)", line)]
        self.assertGreaterEqual(len(calls), 1)
        self.assertEqual(set(calls), {"--no-optional-locks"}, calls)
        for line in text.splitlines():
            if re.search(r"git\s+(status|ls-files|ls-tree|rev-parse|log|branch|config|cat-file)", line) and not line.strip().startswith("#"):
                self.fail("git subcommand invoked outside repo_git(): " + line.strip())

    # ── D8: python validation subprocesses never write bytecode ─────────
    def test_every_python_subprocess_runs_with_dash_B(self):
        self.assertEqual(self.run_script().returncode, 0, self.output)
        calls = (self.state / "python.log").read_text().splitlines()
        self.assertGreaterEqual(len(calls), 1 + sum(1 for r in EXPECTED_ALLOWLIST if r.endswith(".py")))  # deps + one parse per .py
        self.assertTrue(all(c.split()[0] == "-B" for c in calls), calls)

    def test_no_pyc_files_are_created_by_validation(self):
        before = {str(p) for p in self.root.rglob("*.pyc")} | {str(p) for p in self.root.rglob("__pycache__")}
        self.assertEqual(self.run_script("--execute").returncode, 0, self.output)
        after = {str(p) for p in self.root.rglob("*.pyc")} | {str(p) for p in self.root.rglob("__pycache__")}
        self.assertEqual(after, before)

    def test_every_python_invocation_in_the_script_uses_dash_B(self):
        text = SCRIPT.read_text(encoding="utf-8")
        calls = [c for c in re.findall(r'"\$PY"\s+(\S+)', text) if c != "];"]   # ignore the [ -x "$PY" ] guard
        self.assertGreaterEqual(len(calls), 2)
        self.assertEqual(set(calls), {"-B"}, calls)
        code_lines = [l for l in text.splitlines()
                      if not l.lstrip().startswith("#") and not re.match(r"\s*(echo|info|warn|fail|pass|die)\b", l)]
        offenders = [l.strip() for l in code_lines if re.search(r"(?<![\w./-])python3?\s", l)]
        self.assertEqual(offenders, [], "python must only be invoked as \"$PY\" -B ...")

    # ── D2: worktree vs HEAD must be described truthfully ───────────────
    def test_head_mode_output_says_bytes_come_from_the_commit(self):
        self.assertEqual(self.run_script("--expect-commit", self.head).returncode, 0, self.output)
        self.assertIn("file bytes are read from this commit", self.output)
        self.assertNotIn("WORKING TREE", self.output)

    def test_worktree_mode_output_never_claims_the_bytes_are_head(self):
        self.dirty_app()
        result = self.run_script("--from-worktree", "--allow-dirty", "--expect-commit", self.head)
        self.assertEqual(result.returncode, 0, self.output)
        self.assertIn("file bytes come from the WORKING TREE, not from this commit", self.output)
        self.assertIn("pins HEAD only", self.output)
        self.assertNotIn("file bytes are read from this commit", self.output)

    # ── static properties of the script ─────────────────────────────────
    def test_script_has_no_forbidden_commands_or_credentials(self):
        text = SCRIPT.read_text(encoding="utf-8")
        self.assertEqual(scan_text(text, include_comments=True), [])
        for pattern in (r"\bsudo\b", r"\bpip3?\s+install\b", r"\bnohup\b", r"\bpkill\b", r"\bkillall\b",
                        r"python3?\s+app\.py", r"cloudflared", r"conda activate", r"config\.py before"):
            self.assertIsNone(re.search(pattern, text), f"forbidden pattern in deploy_portal.sh: {pattern}")
        self.assertIsNone(re.search(r"systemctl\s+(start|stop|restart|enable|disable|mask|unmask|daemon-reload|kill)", text))
        self.assertIn("set -euo pipefail", text)
        self.assertNotIn(b"\r", SCRIPT.read_bytes())

    def test_script_allowlist_is_exactly_the_reviewed_16_files(self):
        text = SCRIPT.read_text(encoding="utf-8")
        block = re.search(r"readonly ALLOWLIST=\(\n(.*?)\n\)", text, re.S).group(1)
        self.assertEqual([line.strip() for line in block.splitlines() if line.strip()], EXPECTED_ALLOWLIST)
        forbidden_fragments = (".env", ".db", "uploads", "results", "workflow", "resources", "scripts", "tests", "logs")
        for rel in EXPECTED_ALLOWLIST:
            self.assertFalse(any(f in rel for f in forbidden_fragments), rel)

    def test_script_does_not_open_env_or_touch_state_paths(self):
        text = SCRIPT.read_text(encoding="utf-8")
        self.assertIsNone(re.search(r"(cat|source|\.|grep|head|tail|less|cp|mv)\s+[^\n|;]*\$ENV_FILE", text))
        for line in text.splitlines():
            if re.search(r"\b(rm|mv|cp|chmod|mkdir)\b", line) and not line.strip().startswith(("#", "info", "warn", "fail", "echo", "pass")):
                self.assertNotRegex(line, r"\.env|\.db|uploads|portal_logs|results|snakemake", line)


@unittest.skipUnless(shutil.which("git"), "git not available")
class TestNoCredentialsInProductionCode(unittest.TestCase):
    """Tracked production files on the allowlist must hold no credential material."""

    QUOTED_SECRET = re.compile(
        r"""(?ix)\b[a-z0-9_]*(?:pass(?:word|wd)?|secret(?:_key)?|token|api_?key)[a-z0-9_]*\s*=\s*["'][^"'\n]{6,}["']""")
    ENV_LOOKUP = re.compile(r"os\.environ|getenv|_require_env|request\.form|placeholder=|<input|type=")

    @staticmethod
    def exposed_values():
        """The publicly exposed historical values, read from Git history in memory only."""
        first = subprocess.run(["git", "-C", str(ROOT), "log", "--all", "--diff-filter=A", "--format=%H", "--", "portal/config.py"],
                               capture_output=True, text=True).stdout.split()
        if not first:
            return None
        old = subprocess.run(["git", "-C", str(ROOT), "show", f"{first[-1]}:portal/config.py"], capture_output=True, text=True).stdout
        values = []
        for name in ("PORTAL_PASS", "SECRET_KEY", "SMTP_PASS"):
            m = re.search(rf'^{name}\s*=\s*(?:os\.environ\.get\([^,]+,\s*)?["\']([^"\']+)["\']', old, re.M)
            if m:
                values.append(m.group(1))
        return values or None

    def allowlisted_files(self):
        return [ROOT / "portal" / rel for rel in EXPECTED_ALLOWLIST]

    def test_allowlisted_files_exist(self):
        for path in self.allowlisted_files():
            self.assertTrue(path.is_file(), f"{path.relative_to(ROOT)} missing")

    def test_no_credential_like_quoted_literals(self):
        offenders = []
        for path in self.allowlisted_files():
            if path.suffix == ".png":
                continue
            for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
                if self.QUOTED_SECRET.search(line) and not self.ENV_LOOKUP.search(line):
                    offenders.append(f"{path.relative_to(ROOT)}:{lineno}")
        self.assertEqual(offenders, [], "credential-like literal(s) in deployable production code (values not shown)")

    def test_no_publicly_exposed_historical_value_in_any_deployable_file(self):
        values = self.exposed_values()
        if values is None:
            self.skipTest("Git history not available: cannot derive the exposed values")
        offenders = []
        for path in self.allowlisted_files():
            raw = path.read_bytes()
            if any(v.encode() in raw for v in values):
                offenders.append(str(path.relative_to(ROOT)))
        self.assertEqual(offenders, [], "a publicly exposed historical secret is present in a deployable file (value not shown)")

    def test_head_blobs_of_the_allowlist_are_clean_too(self):
        values = self.exposed_values()
        head = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "--verify", "HEAD"], capture_output=True, text=True)
        if values is None or head.returncode != 0:
            self.skipTest("Git history not available")
        offenders = []
        for rel in EXPECTED_ALLOWLIST:
            blob = subprocess.run(["git", "-C", str(ROOT), "show", f"HEAD:portal/{rel}"], capture_output=True).stdout
            if any(v.encode() in blob for v in values):
                offenders.append(rel)
        self.assertEqual(offenders, [], "HEAD contains an exposed value in a deployable file (value not shown)")

    def test_config_reads_secrets_from_the_environment_only(self):
        text = (ROOT / "portal" / "config.py").read_text(encoding="utf-8")
        self.assertIn('_require_env("SECRET_KEY")', text)
        self.assertIn('_require_env("PORTAL_PASS")', text)
        self.assertIsNone(re.search(r"(?im)^(?:SECRET_KEY|PORTAL_PASS|SMTP_PASS)\s*=\s*[\"']", text))

    def test_env_file_is_not_tracked_or_deployable(self):
        tracked = subprocess.run(["git", "-C", str(ROOT), "ls-files"], capture_output=True, text=True).stdout.split()
        if not tracked:
            self.skipTest("no Git metadata")
        self.assertNotIn("portal/.env", tracked)
        self.assertNotIn(".env", EXPECTED_ALLOWLIST)


if __name__ == "__main__":
    unittest.main()
