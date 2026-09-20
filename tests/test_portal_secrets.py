"""
GenRichi Phase 5.5 Step 14, Stage 1 -- portal secret-loading regression tests.

portal/config.py no longer contains plaintext SECRET_KEY / PORTAL_PASS /
SMTP_PASS literals. It reads them from the process environment (optionally
pre-populated from the gitignored portal/.env via a small dependency-free
loader), and raises RuntimeError -- clearly, with no fallback placeholder --
if a required one is missing.

Each `import config` here happens in a **fresh subprocess** (never the test
runner's own process). This is deliberate, not incidental: `portal/config.py`
is imported as the bare module name "config" (portal/ has no __init__.py),
so a second `import config` anywhere else in the same process would just
return the first cached module object without re-running its top-level
code -- silently hiding exactly the fail-fast/env-loading behavior this file
exists to prove. A subprocess sidesteps that entirely and, as a side benefit,
lets each test control the process environment precisely without ever
touching this test process's own os.environ.

As of Phase 5.5 Step 14, `_run_import_config` additionally runs from a
temporary COPY of config.py rather than the real portal/ directory, so
these tests stay correct and fully isolated now that a genuine, rotated
portal/.env legitimately exists on disk there (see that function's own
docstring for the full reasoning). Neither portal/.env nor portal/config.py
is read, written, or deleted by anything in this file.

No real secret value is used, printed, or asserted anywhere in this file --
only clearly-fake placeholder strings that exist solely to prove the
loading mechanism works, plus assertions that those placeholders (and,
separately, the historical real values) are absent from error output.
"""

import contextlib
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
PORTAL_DIR = REPO_ROOT / "portal"
CONFIG_PY = PORTAL_DIR / "config.py"
APP_PY = PORTAL_DIR / "app.py"

FAKE_SECRET_KEY = "unit-test-fake-secret-key-0001"
FAKE_PORTAL_PASS = "unit-test-fake-portal-pass-0002"
FAKE_SMTP_PASS = "unit-test-fake-smtp-pass-0003"

# The historical plaintext values found by the Step 13 security audit are NOT
# stored anywhere in the working tree. The detection targets are derived from
# Git history at test run time (those values were once assigned in
# portal/config.py), held only in memory, and never printed or placed in an
# assertion message. If Git history is unavailable the affected tests skip.
_SECRET_NAMES = ("SECRET_KEY", "PORTAL_PASS", "SMTP_PASS")
_FRAGMENT_LEN = 12


def _git_output(*args):
    try:
        done = subprocess.run(
            ["git", "--no-optional-locks", "-C", str(REPO_ROOT), *args],
            capture_output=True, timeout=120,
        )
    except (OSError, subprocess.SubprocessError):
        return ""
    return done.stdout.decode("utf-8", "replace") if done.returncode == 0 else ""


def _historical_secret_targets():
    """Every value once assigned as a quoted literal to a secret-named
    variable in any version of portal/config.py, plus a short leading
    fragment of each. Returns a set of str (possibly empty)."""
    targets = set()
    for sha in _git_output("log", "--all", "--format=%H", "--", "portal/config.py").split():
        source = _git_output("show", f"{sha}:portal/config.py")
        for name in _SECRET_NAMES:
            pattern = (rf'{name}["\']?\s*,\s*["\']([^"\']{{6,}})["\']'
                       rf'|^{name}\s*=\s*["\']([^"\']{{6,}})["\']')
            for match in re.finditer(pattern, source, re.M):
                value = match.group(1) or match.group(2)
                if not value or re.fullmatch(r"[A-Z0-9_]+", value):
                    continue  # an identifier such as a variable name, not a secret
                targets.add(value)
                if len(value) > _FRAGMENT_LEN:
                    targets.add(value[:_FRAGMENT_LEN])
    return targets


def _publishable_files():
    """Tracked files plus untracked-but-not-ignored files (what a commit or
    push could publish). The gitignored runtime portal/.env is excluded."""
    names = set(_git_output("ls-files").splitlines())
    names |= set(_git_output("ls-files", "--others", "--exclude-standard").splitlines())
    skip_suffixes = (".png", ".jpg", ".jpeg", ".gif", ".pdf", ".docx", ".gz", ".bam", ".bai", ".fa", ".fasta", ".zip")
    for name in sorted(names):
        path = REPO_ROOT / name
        if name.lower().endswith(skip_suffixes) or not path.is_file():
            continue
        try:
            if path.stat().st_size > 2_000_000:
                continue
        except OSError:
            continue
        yield name, path


@contextlib.contextmanager
def _isolated_process_state(env_var_names):
    """Snapshot the given environment variable names, sys.modules["config"],
    and sys.path around the wrapped block, then restore each EXACTLY to its
    prior state (present-with-original-value, or absent) once the block
    exits -- never a blind delete/leave-mutated. Order-independent: this
    works correctly no matter what was already there before the test ran,
    including nothing at all."""
    env_existed = {name: name in os.environ for name in env_var_names}
    env_snapshot = {name: os.environ.get(name) for name in env_var_names}
    had_config_module = "config" in sys.modules
    orig_config_module = sys.modules.get("config")
    orig_sys_path = list(sys.path)
    try:
        yield
    finally:
        for name in env_var_names:
            if env_existed[name]:
                os.environ[name] = env_snapshot[name]
            else:
                os.environ.pop(name, None)
        if had_config_module:
            sys.modules["config"] = orig_config_module
        else:
            sys.modules.pop("config", None)
        sys.path[:] = orig_sys_path


def _run_import_config(env_overrides, code="import config"):
    """Run `code` in a fresh Python subprocess importing an isolated
    temporary COPY of the real, unmodified portal/config.py -- never the
    real portal/ directory itself.

    This is deliberate: portal/config.py derives BASE_DIR from its own
    __file__, and a genuine, gitignored portal/.env now legitimately
    exists there (rotated real secrets -- see Phase 5.5 Step 14). Running
    an unmodified copy of config.py from a fresh temp directory with no
    .env file next to it means _load_dotenv() finds nothing to load --
    a real absence, not a faked one -- without ever touching, reading
    from, or deleting the real portal/.env. Every test in this module
    that needs to prove "missing secret" behavior therefore still proves
    it genuinely, regardless of what is or isn't configured on this
    machine. Tests that supply every required value via env_overrides are
    unaffected either way, since no .env is needed once all values are
    already present. config.py's own source is copied byte-for-byte,
    never edited -- production code is untouched.

    PATH (and, on Windows, a few more inherited vars) is preserved so the
    interpreter itself can start; nothing else from this test process's
    real os.environ leaks in."""
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        shutil.copy2(str(CONFIG_PY), str(tmp_path / "config.py"))

        env = {"PATH": os.environ.get("PATH", "")}
        if sys.platform == "win32":
            # Python on Windows needs a few more inherited vars to start cleanly.
            for k in ("SYSTEMROOT", "TEMP", "TMP"):
                if k in os.environ:
                    env[k] = os.environ[k]
        env.update(env_overrides)
        return subprocess.run(
            [sys.executable, "-c", code],
            cwd=str(tmp_path),
            env=env,
            capture_output=True,
            text=True,
            timeout=30,
        )


class TestMissingSecretsFailSafely(unittest.TestCase):
    def test_missing_secret_key_raises_clear_runtime_error(self):
        result = _run_import_config({"PORTAL_PASS": FAKE_PORTAL_PASS, "SMTP_ENABLED": "false"})
        self.assertNotEqual(result.returncode, 0, "import must fail when SECRET_KEY is unset")
        self.assertIn("SECRET_KEY", result.stderr)
        self.assertIn("RuntimeError", result.stderr)

    def test_missing_portal_pass_raises_clear_runtime_error(self):
        result = _run_import_config({"SECRET_KEY": FAKE_SECRET_KEY, "SMTP_ENABLED": "false"})
        self.assertNotEqual(result.returncode, 0, "import must fail when PORTAL_PASS is unset")
        self.assertIn("PORTAL_PASS", result.stderr)
        self.assertIn("RuntimeError", result.stderr)

    def test_smtp_pass_required_only_when_smtp_enabled(self):
        # SMTP disabled: must succeed with no SMTP_PASS at all.
        ok = _run_import_config({
            "SECRET_KEY": FAKE_SECRET_KEY,
            "PORTAL_PASS": FAKE_PORTAL_PASS,
            "SMTP_ENABLED": "false",
        })
        self.assertEqual(ok.returncode, 0, ok.stderr)

        # SMTP enabled, SMTP_PASS missing: must fail clearly.
        fail = _run_import_config({
            "SECRET_KEY": FAKE_SECRET_KEY,
            "PORTAL_PASS": FAKE_PORTAL_PASS,
            "SMTP_ENABLED": "true",
        })
        self.assertNotEqual(fail.returncode, 0)
        self.assertIn("SMTP_PASS", fail.stderr)
        self.assertIn("RuntimeError", fail.stderr)

    def test_no_stray_default_credential_is_used(self):
        """Fail-fast must never fall through to a placeholder value."""
        result = _run_import_config(
            {"PORTAL_PASS": FAKE_PORTAL_PASS, "SMTP_ENABLED": "false"},
            code="import config; print('SHOULD NOT REACH HERE')",
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("SHOULD NOT REACH HERE", result.stdout)


class TestEnvironmentVariablesLoadCorrectly(unittest.TestCase):
    def test_provided_env_values_are_used_verbatim(self):
        result = _run_import_config(
            {
                "SECRET_KEY": FAKE_SECRET_KEY,
                "PORTAL_PASS": FAKE_PORTAL_PASS,
                "SMTP_ENABLED": "true",
                "SMTP_PASS": FAKE_SMTP_PASS,
            },
            code=(
                "import config; "
                "print(config.SECRET_KEY); "
                "print(config.PORTAL_PASS); "
                "print(config.SMTP_PASS); "
                "print(config.SMTP_ENABLED)"
            ),
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        lines = result.stdout.splitlines()
        self.assertEqual(lines[0], FAKE_SECRET_KEY)
        self.assertEqual(lines[1], FAKE_PORTAL_PASS)
        self.assertEqual(lines[2], FAKE_SMTP_PASS)
        self.assertEqual(lines[3], "True")

    def test_dotenv_file_is_loaded_without_overriding_real_env(self):
        """Unit-tests the _load_dotenv() loader function itself, directly,
        against a throwaway temp file. (The genuine end-to-end path -- a
        real .env file driving a real config import with no other env vars
        set -- is covered separately by TestRealEndToEndDotenvPath below;
        that one doesn't need this method's real-os.environ/sys.modules
        mutation at all, since it runs entirely in a subprocess.)

        This test mutates the real process's os.environ, sys.modules, and
        sys.path for the duration of the `with` block below, but
        _isolated_process_state snapshots and restores all of it exactly --
        present-with-original-value or absent, never a blind delete -- so
        this test is safe to run regardless of what, if anything, was set
        in the environment before it started, and regardless of test
        discovery/execution order. The assertions after the `with` block
        prove that restoration actually happened.
        """
        env_names = ["SECRET_KEY", "PORTAL_PASS", "SMTP_ENABLED", "GENRICHI_TEST_DOTENV_PROBE"]
        pre_env = {name: (name in os.environ, os.environ.get(name)) for name in env_names}
        pre_had_config_module = "config" in sys.modules
        pre_config_module = sys.modules.get("config")
        pre_sys_path = list(sys.path)

        with _isolated_process_state(env_names):
            with tempfile.TemporaryDirectory() as tmp:
                dotenv_path = Path(tmp) / ".env"
                sys.path.insert(0, str(PORTAL_DIR))
                if "config" in sys.modules:
                    del sys.modules["config"]
                # config.py only ever reads portal/.env (a fixed, hardcoded
                # path) at import time -- it has no notion of this test's
                # temp file. Supply fake env vars so the module itself
                # loads successfully; the loader function is then exercised
                # directly, against our own temp file, below.
                os.environ["SECRET_KEY"] = FAKE_SECRET_KEY
                os.environ["PORTAL_PASS"] = FAKE_PORTAL_PASS
                os.environ["SMTP_ENABLED"] = "false"
                import config as portal_config  # noqa: F401  (only need _load_dotenv)

                probe_env_key = "GENRICHI_TEST_DOTENV_PROBE"
                os.environ.pop(probe_env_key, None)
                dotenv_path.write_text(f'{probe_env_key}=from-dotenv\n')
                portal_config._load_dotenv(str(dotenv_path))
                self.assertEqual(os.environ.get(probe_env_key), "from-dotenv")

                # A key already set in the real environment must NOT be
                # overridden by the .env file.
                os.environ[probe_env_key] = "from-real-environment"
                dotenv_path.write_text(f'{probe_env_key}=from-dotenv-should-be-ignored\n')
                portal_config._load_dotenv(str(dotenv_path))
                self.assertEqual(os.environ.get(probe_env_key), "from-real-environment")

        # Process state must be back to exactly what it was before this
        # test ran -- whatever that was, including "unset".
        for name in env_names:
            had, val = pre_env[name]
            with self.subTest(env_var=name):
                self.assertEqual(name in os.environ, had)
                if had:
                    self.assertEqual(os.environ.get(name), val)
        self.assertEqual("config" in sys.modules, pre_had_config_module)
        if pre_had_config_module:
            self.assertIs(sys.modules.get("config"), pre_config_module)
        self.assertEqual(sys.path, pre_sys_path)


class TestRealEndToEndDotenvPath(unittest.TestCase):
    """
    Exercises the genuine production sequence end to end:

        portal/.env -> _load_dotenv() -> os.environ -> _require_env() -> successful config import

    without ever creating or touching the real repository's portal/.env.

    portal/config.py derives BASE_DIR from its own __file__, so the only
    way to point its *own, unmodified* _load_dotenv() call at a different,
    throwaway .env file -- without weakening production code just to make
    it testable -- is to run a temporary copy of config.py from a temp
    directory: BASE_DIR then naturally resolves to that temp directory,
    and _load_dotenv() naturally looks for <tempdir>/.env, which this test
    fully controls and which is deleted with the temp directory afterward.
    config.py's own source is copied byte-for-byte, never edited.

    Every import happens in a fresh subprocess with an explicit, minimal
    environment (see _run_with_temp_env_file), so this test never touches
    this test process's own os.environ, sys.modules, or sys.path either.
    """

    def _run_with_temp_env_file(self, dotenv_contents, proc_env_overrides, code):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            shutil.copy2(str(CONFIG_PY), str(tmp_path / "config.py"))
            if dotenv_contents is not None:
                (tmp_path / ".env").write_text(dotenv_contents)

            env = {"PATH": os.environ.get("PATH", "")}
            if sys.platform == "win32":
                for k in ("SYSTEMROOT", "TEMP", "TMP"):
                    if k in os.environ:
                        env[k] = os.environ[k]
            env.update(proc_env_overrides)
            return subprocess.run(
                [sys.executable, "-c", code],
                cwd=str(tmp_path),
                env=env,
                capture_output=True,
                text=True,
                timeout=30,
            )

    def test_secret_key_and_portal_pass_loaded_from_real_dotenv_file(self):
        # Deliberately NO SECRET_KEY / PORTAL_PASS / SMTP_ENABLED in the
        # subprocess's own environment (proc_env_overrides is {}) -- they
        # must come only from the temp .env file written below.
        dotenv_contents = (
            f"SECRET_KEY={FAKE_SECRET_KEY}\n"
            f"PORTAL_PASS={FAKE_PORTAL_PASS}\n"
            f"SMTP_ENABLED=false\n"
        )
        result = self._run_with_temp_env_file(
            dotenv_contents,
            {},
            "import config; "
            "print(config.SECRET_KEY); "
            "print(config.PORTAL_PASS); "
            "print(config.SMTP_ENABLED)",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        lines = result.stdout.splitlines()
        self.assertEqual(lines[0], FAKE_SECRET_KEY)
        self.assertEqual(lines[1], FAKE_PORTAL_PASS)
        self.assertEqual(lines[2], "False")

    def test_process_env_var_takes_precedence_over_dotenv_value_end_to_end(self):
        dotenv_only_value = "unit-test-value-from-dotenv-must-be-overridden"
        process_env_value = "unit-test-value-from-real-process-environment"
        dotenv_contents = (
            f"SECRET_KEY={dotenv_only_value}\n"
            f"PORTAL_PASS={FAKE_PORTAL_PASS}\n"
            f"SMTP_ENABLED=false\n"
        )
        result = self._run_with_temp_env_file(
            dotenv_contents,
            {"SECRET_KEY": process_env_value},
            "import config; print(config.SECRET_KEY)",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.strip(), process_env_value)


class TestNoPlaintextSecretsInTrackedSource(unittest.TestCase):
    ASSIGNMENT_RE = staticmethod(
        lambda name: re.compile(rf'^\s*{name}\s*=\s*["\']', re.MULTILINE)
    )

    def test_config_py_has_no_literal_secret_assignment(self):
        text = CONFIG_PY.read_text(encoding="utf-8")
        for name in ("SECRET_KEY", "PORTAL_PASS", "SMTP_PASS"):
            with self.subTest(name=name):
                self.assertIsNone(
                    self.ASSIGNMENT_RE(name).search(text),
                    f"{name} must not be assigned a literal string in portal/config.py",
                )

    def test_historical_leaked_values_are_absent_from_tracked_source(self):
        targets = _historical_secret_targets()
        if not targets:
            self.skipTest("Git history unavailable: cannot derive the historical secret values")
        for path in (CONFIG_PY, APP_PY):
            data = path.read_bytes()
            with self.subTest(file=path.name):
                self.assertFalse(
                    any(t.encode("utf-8") in data for t in targets),
                    "a historical secret value is present (value not shown)",
                )

    def test_historical_leaked_values_are_absent_from_every_publishable_file(self):
        """No tracked or untracked-but-not-ignored file (tests, scripts, docs,
        fixtures) may carry a historical secret value or its leading fragment."""
        targets = _historical_secret_targets()
        if not targets:
            self.skipTest("Git history unavailable: cannot derive the historical secret values")
        needles = [t.encode("utf-8") for t in targets]
        offenders = []
        for name, path in _publishable_files():
            try:
                data = path.read_bytes()
            except OSError:
                continue
            if any(n in data for n in needles):
                offenders.append(name)
        self.assertEqual(offenders, [], "historical secret value(s) found in publishable file(s) (values not shown)")

    def test_this_file_declares_no_historical_secret_literals(self):
        source = Path(__file__).read_text(encoding="utf-8")
        self.assertNotRegex(source, r"(?m)^HISTORICAL_[A-Z_]+\s*=")


if __name__ == "__main__":
    unittest.main()
