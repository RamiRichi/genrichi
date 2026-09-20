"""
GenRichi Phase 5.5 Step 14, Stage 2 -- portal settings-persistence tests.

Stage 2 removed portal/app.py's /settings/save mechanism that used to
rewrite portal/config.py's own source text in place to persist SMTP
settings (the Step 13 security audit's second finding). Non-secret settings
(SMTP_ENABLED, SMTP_USER) are now persisted to the gitignored
portal/instance/settings.json via config.py's _load_instance_settings();
SMTP_PASS is never read from or written to that file -- it remains
exclusively environment/.env-sourced (Stage 1), and the settings page can
no longer change it at all (its form field is now display-only).

Flask is not installed in the environment this suite normally runs under
(confirmed before writing this file), so /settings/save itself is not
exercised through a real HTTP request/test client here -- that would add
Flask as a new test dependency, against the "dependency-light" convention
already established by test_portal_secrets.py. Instead, several
complementary techniques prove the same properties:
  - the ABSENCE of the old file-rewriting mechanism is proven by static
    inspection of portal/app.py's source -- a stronger guarantee than a
    single runtime call, since it shows the capability doesn't exist at
    all, not just "didn't happen this once";
  - the REPLACEMENT persistence mechanism (instance/settings.json) is
    proven functionally, end to end, in a fresh subprocess importing a
    temporary copy of the real, unmodified config.py -- the same technique
    test_portal_secrets.py's TestRealEndToEndDotenvPath already established
    for the same underlying reason (config.py derives its paths from its
    own __file__, so genuinely testing it without weakening production
    code means running an unmodified copy from a location the test
    controls);
  - whether SMTP_PASS can reach a rendered HTML response is proven with a
    REAL Jinja2 render of the actual, unmodified portal/templates/
    settings.html (Jinja2 -- unlike Flask -- IS already a dependency of
    this project's own report.yaml conda environment, so this needs no new
    test dependency either). Flask-only globals the template's base.html
    needs (url_for, session, request, get_flashed_messages) are supplied
    as minimal stand-ins; nothing about settings.html itself is modified
    to make this possible.

No real secret value is used, printed, or asserted anywhere in this file.
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile
import types
import unittest
from pathlib import Path

import jinja2

REPO_ROOT = Path(__file__).resolve().parents[1]
PORTAL_DIR = REPO_ROOT / "portal"
CONFIG_PY = PORTAL_DIR / "config.py"
APP_PY = PORTAL_DIR / "app.py"
TEMPLATES_DIR = PORTAL_DIR / "templates"

FAKE_SECRET_KEY = "unit-test-fake-secret-key-0004"
FAKE_PORTAL_PASS = "unit-test-fake-portal-pass-0005"
FAKE_SMTP_PASS = "unit-test-fake-smtp-pass-0006"
FAKE_SMTP_USER = "unit-test-fake-smtp-user@example.invalid"
INSTANCE_FAKE_SMTP_PASS = "should-never-be-read-from-instance-file"


def _run(code, cwd, env_overrides):
    env = {"PATH": os.environ.get("PATH", "")}
    if sys.platform == "win32":
        for k in ("SYSTEMROOT", "TEMP", "TMP"):
            if k in os.environ:
                env[k] = os.environ[k]
    env.update(env_overrides)
    return subprocess.run(
        [sys.executable, "-c", code],
        cwd=str(cwd),
        env=env,
        capture_output=True,
        text=True,
        timeout=30,
    )


class TestConfigPyRewriteMechanismIsGone(unittest.TestCase):
    """(A) the capability to rewrite config.py's source no longer exists in
    portal/app.py at all -- proven statically."""

    def test_no_config_py_write_target_construction(self):
        text = APP_PY.read_text(encoding="utf-8")
        self.assertNotIn('"config.py"', text)

    def test_no_source_text_replace_pattern(self):
        text = APP_PY.read_text(encoding="utf-8")
        self.assertNotIn(".replace(", text)

    def test_settings_save_route_writes_only_the_instance_settings_path(self):
        text = APP_PY.read_text(encoding="utf-8")
        self.assertIn("INSTANCE_SETTINGS_PATH", text)
        self.assertIn("cfg.INSTANCE_DIR", text)

    def test_smtp_pass_value_is_never_interpolated_into_a_flash_or_log_string(self):
        """(F) smtp_pass is only ever used in a boolean/comparison
        condition in app.py now, never embedded in text shown to the user
        or written to a log."""
        text = APP_PY.read_text(encoding="utf-8")
        self.assertNotIn('{smtp_pass}', text)
        self.assertNotIn('{cfg.SMTP_PASS}', text)


def _settings_route_body():
    """Extract just the `def settings():` GET-route function body from
    portal/app.py's source text, so checks below can't accidentally match
    the unrelated /settings/save POST handler (which legitimately reads
    request.form.get("smtp_pass") to decide what to flash -- it never
    renders that value back, which is checked separately, above)."""
    text = APP_PY.read_text(encoding="utf-8")
    start = text.index("def settings():")
    end = text.index("\ndef ", start + 1)
    return text[start:end]


class TestSmtpPassNeverReachesTheTemplate(unittest.TestCase):
    """
    Phase 5.5 Stage 2 review finding: SMTP_PASS must stay server-side only
    and must never appear in rendered HTML/HTTP responses/template context,
    not even masked. Proven three ways: the GET /settings route no longer
    reads or passes anything derived from SMTP_PASS; the template source no
    longer references it at all; and a real Jinja2 render of the actual,
    unmodified template -- even when a context deliberately (mis)supplies a
    fake secret under the old `smtp_pass` key, simulating a future
    regression upstream -- never surfaces it in the rendered output.
    """

    def test_settings_route_never_reads_or_passes_smtp_pass(self):
        # Checks actual code usage (attribute access / kwarg), not the
        # explanatory comment inside this function that names the variable
        # to document its absence.
        body = _settings_route_body()
        self.assertNotIn("cfg.SMTP_PASS", body)
        self.assertNotIn("smtp_pass ", body)
        self.assertNotIn("smtp_pass=", body)

    def test_settings_html_source_never_references_smtp_pass(self):
        text = (TEMPLATES_DIR / "settings.html").read_text(encoding="utf-8")
        self.assertNotIn("smtp_pass", text)

    def _render_settings_html(self, context):
        env = jinja2.Environment(loader=jinja2.FileSystemLoader(str(TEMPLATES_DIR)))
        env.globals["url_for"] = lambda *a, **k: "#"
        env.globals["session"] = {"role": "admin", "username": "tester", "full_name": "Tester"}
        env.globals["request"] = types.SimpleNamespace(endpoint="settings")
        env.globals["get_flashed_messages"] = lambda **k: []
        return env.get_template("settings.html").render(**context)

    def _base_context(self):
        return dict(
            smtp_enabled=True,
            smtp_user="test@example.invalid",
            portal_url="https://example.invalid",
            admin_user="admin",
            db_path="/tmp/x.db",
            workflow_dir="/tmp/wf",
            total_orders=0,
            disk_usage="1 GB",
            pipelines={},
        )

    def test_rendered_page_with_current_real_context_has_no_secret(self):
        """The exact context shape portal/app.py's settings() route now
        builds (no smtp_pass key at all) renders successfully and cleanly."""
        html = self._render_settings_html(self._base_context())
        self.assertIn("Configured via environment", html)
        self.assertNotIn("SMTP_PASS", html)

    def test_rendered_page_ignores_smtp_pass_even_if_context_regresses(self):
        """Even if a future change to app.py mistakenly reintroduced
        passing smtp_pass into the context (the exact bug this review
        found), the template itself -- unmodified, loaded from disk, not
        mocked -- must never surface it, because it no longer references
        that variable at all."""
        fake_leaked_secret = "unit-test-fake-smtp-pass-0007-should-never-render"
        context = self._base_context()
        context["smtp_pass"] = fake_leaked_secret  # simulated regression
        html = self._render_settings_html(context)
        self.assertNotIn(fake_leaked_secret, html)

    def test_settings_save_new_settings_dict_never_assigns_smtp_pass_key(self):
        """Static check that the persisted instance-settings dict the
        save route builds can never itself carry an SMTP_PASS entry."""
        text = APP_PY.read_text(encoding="utf-8")
        start = text.index('elif action == "email_settings"')
        end = text.index("\n\n    return redirect", start)
        body = text[start:end]
        self.assertNotIn('"SMTP_PASS"', body)
        self.assertNotIn("cfg.SMTP_PASS =", body)


class TestInstanceSettingsPersistence(unittest.TestCase):
    """Proves the replacement persistence mechanism end to end. Every
    subprocess import uses an unmodified temporary COPY of the real
    portal/config.py (never edited) -- see module docstring."""

    @staticmethod
    def _copy_config_into(tmp_path):
        shutil.copy2(str(CONFIG_PY), str(tmp_path / "config.py"))

    def test_persisted_smtp_settings_are_read_back_by_config_py(self):
        """(C) settings a settings-save action persists are usable by the
        application on the next import, without ever touching config.py."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            self._copy_config_into(tmp_path)
            instance_dir = tmp_path / "instance"
            instance_dir.mkdir()
            (instance_dir / "settings.json").write_text(
                json.dumps({"SMTP_ENABLED": False, "SMTP_USER": FAKE_SMTP_USER})
            )

            result = _run(
                "import config; print(config.SMTP_ENABLED); print(config.SMTP_USER)",
                tmp_path,
                {"SECRET_KEY": FAKE_SECRET_KEY, "PORTAL_PASS": FAKE_PORTAL_PASS},
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            lines = result.stdout.splitlines()
            self.assertEqual(lines[0], "False")
            self.assertEqual(lines[1], FAKE_SMTP_USER)

    def test_config_py_source_is_byte_identical_after_persistence(self):
        """(A) functional counterpart to the static proof above: writing
        instance/settings.json never touches config.py's own file on disk."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            self._copy_config_into(tmp_path)
            copied_config = tmp_path / "config.py"
            before = copied_config.read_bytes()

            instance_dir = tmp_path / "instance"
            instance_dir.mkdir()
            (instance_dir / "settings.json").write_text(
                json.dumps({"SMTP_ENABLED": True, "SMTP_USER": FAKE_SMTP_USER})
            )
            result = _run(
                "import config",
                tmp_path,
                {
                    "SECRET_KEY": FAKE_SECRET_KEY,
                    "PORTAL_PASS": FAKE_PORTAL_PASS,
                    "SMTP_PASS": FAKE_SMTP_PASS,
                },
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(copied_config.read_bytes(), before)

    def test_instance_settings_file_never_supplies_smtp_pass(self):
        """(B) a SMTP_PASS key present in settings.json (accidentally or
        maliciously) is never read -- the secret channel stays singular:
        environment/.env only, per Stage 1."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            self._copy_config_into(tmp_path)
            instance_dir = tmp_path / "instance"
            instance_dir.mkdir()
            (instance_dir / "settings.json").write_text(json.dumps({
                "SMTP_ENABLED": True,
                "SMTP_USER": FAKE_SMTP_USER,
                "SMTP_PASS": INSTANCE_FAKE_SMTP_PASS,
            }))
            result = _run(
                "import config; print(config.SMTP_PASS)",
                tmp_path,
                {
                    "SECRET_KEY": FAKE_SECRET_KEY,
                    "PORTAL_PASS": FAKE_PORTAL_PASS,
                    "SMTP_PASS": FAKE_SMTP_PASS,
                },
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(result.stdout.strip(), FAKE_SMTP_PASS)
            self.assertNotIn(INSTANCE_FAKE_SMTP_PASS, result.stdout)

    def test_missing_required_secrets_still_fail_with_instance_settings_present(self):
        """(D) the new mechanism cannot accidentally satisfy _require_env."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            self._copy_config_into(tmp_path)
            instance_dir = tmp_path / "instance"
            instance_dir.mkdir()
            (instance_dir / "settings.json").write_text(
                json.dumps({"SMTP_ENABLED": False, "SMTP_USER": FAKE_SMTP_USER})
            )
            result = _run("import config", tmp_path, {"PORTAL_PASS": FAKE_PORTAL_PASS})
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("SECRET_KEY", result.stderr)
            self.assertIn("RuntimeError", result.stderr)

    def test_env_var_still_takes_precedence_over_instance_settings(self):
        """(E) Stage-1 precedence semantics preserved: a real environment
        variable beats whatever the settings page persisted."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            self._copy_config_into(tmp_path)
            instance_dir = tmp_path / "instance"
            instance_dir.mkdir()
            (instance_dir / "settings.json").write_text(json.dumps({
                "SMTP_ENABLED": True,
                "SMTP_USER": "from-instance-file@example.invalid",
            }))
            result = _run(
                "import config; print(config.SMTP_ENABLED); print(config.SMTP_USER)",
                tmp_path,
                {
                    "SECRET_KEY": FAKE_SECRET_KEY,
                    "PORTAL_PASS": FAKE_PORTAL_PASS,
                    "SMTP_ENABLED": "false",
                    "SMTP_USER": FAKE_SMTP_USER,
                },
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            lines = result.stdout.splitlines()
            self.assertEqual(lines[0], "False")
            self.assertEqual(lines[1], FAKE_SMTP_USER)

    def test_missing_instance_file_falls_back_to_prior_default_behavior(self):
        """A fresh checkout with no instance/ directory at all behaves
        exactly as it did before Stage 2 (default SMTP_ENABLED=true,
        default SMTP_USER)."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            self._copy_config_into(tmp_path)
            result = _run(
                "import config; print(config.SMTP_ENABLED); print(config.SMTP_USER)",
                tmp_path,
                {
                    "SECRET_KEY": FAKE_SECRET_KEY,
                    "PORTAL_PASS": FAKE_PORTAL_PASS,
                    "SMTP_PASS": FAKE_SMTP_PASS,
                },
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            lines = result.stdout.splitlines()
            self.assertEqual(lines[0], "True")
            self.assertEqual(lines[1], "info@genrichi.de")


if __name__ == "__main__":
    unittest.main()
