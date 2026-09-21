"""
Guard against unmanaged launch paths and embedded credentials in launcher scripts.

systemd is the only supervisor of the production portal and tunnel
(genrichi-portal.service, cloudflared.service; genrichi.service stays disabled).
Tracked shell/batch files must therefore not:

  * kill processes (pkill / killall / kill <pid>)
  * detach processes with nohup
  * start the portal by hand (python app.py, flask run)
  * start the tunnel by hand (cloudflared ... tunnel ... run)
  * start/stop/disable/mask units, or enable/restart any unit other than
    genrichi-portal
  * contain credential-like literals (passwords, secrets, tokens, "Login: user / pw")

scripts/portal_status.sh is held to a stricter, read-only standard.
Findings are reported as file:line + rule name only; matched text is never printed.

Run with:
    python -m unittest tests.test_launchers_safe -v
"""

import re
import shutil
import subprocess
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STATUS_SCRIPT = ROOT / "scripts" / "portal_status.sh"
LOCAL_BAT = ROOT / "GenRichi_Start.bat"          # gitignored, local only
REMOVED_LAUNCHER = ROOT / "start_genrichi.sh"

LAUNCHER_SUFFIXES = (".sh", ".bat", ".cmd")

# Findings that exist today and are scheduled to be fixed. Every entry must still
# be a real finding (test_pending_allowlist_is_not_stale), so an entry has to be
# deleted as soon as its fix lands. Nothing may be added here for a new finding.
# (Empty: the deploy_portal.sh findings were fixed when it was hardened.)
KNOWN_PENDING = {}

# ── rules ────────────────────────────────────────────────────────────────
UNMANAGED_RULES = {
    "pkill-killall": re.compile(r"\b(pkill|killall|skill)\b"),
    "kill-pid": re.compile(r"(?<![\w./-])kill\s+(?:-[A-Za-z0-9]+\s+)?[\$%\d]"),
    "nohup": re.compile(r"\bnohup\b"),
    "python-app": re.compile(r"\bpython3?(?:\.\d+)?\s+(?:\S*[/\\])?app\.py\b"),
    "flask-run": re.compile(r"\bflask\b[^\n#]*\brun\b"),
    "cloudflared-tunnel-run": re.compile(r"\bcloudflared\b[^\n#]*\btunnel\b[^\n#]*\brun\b"),
}

PLACEHOLDER = re.compile(
    r"(?i)^(<.*>|\*+|x{3,}|change_?me|your[_-]?\w*|example\w*|placeholder|todo|none|null|true|false|\.{3})$"
)
CREDENTIAL_RULES = {
    "login-echo": re.compile(r"(?i)\b(?:login|user(?:name)?)\s*:\s*\S+\s*/\s*(?P<v>\S+)"),
    "password-literal": re.compile(r"""(?i)\b(?:password|passwd)\b\s*[:=]\s*["']?(?P<v>[^\s"'<>$`]{3,})"""),
    "secret-assignment": re.compile(
        r"""(?ix)\b[A-Z0-9_]*(?:PASS(?:WORD)?|SECRET|TOKEN|API_?KEY|PRIVATE_?KEY)[A-Z0-9_]*\s*=\s*["']?(?P<v>[^\s"'$<>{}()`]{4,})"""),
    "secret-cli-arg": re.compile(r"""(?i)--(?:token|password|passwd|secret)[= ]+(?P<v>[^\s"'$<>]{6,})"""),
    "token-shaped": re.compile(r"(?P<v>\beyJ[A-Za-z0-9_-]{20,}|\b(?:ghp_|github_pat_|xox[baprs]-|AKIA)[A-Za-z0-9_-]{10,})"),
}

MUTATING_SYSTEMCTL = {"start", "stop", "restart", "try-restart", "reload", "reload-or-restart", "disable",
                      "mask", "unmask", "kill", "isolate", "enable", "edit", "set-property", "reset-failed"}
PORTAL_ONLY_VERBS = {"enable", "restart"}        # allowed, but only for genrichi-portal
PORTAL_UNITS = {"genrichi-portal", "genrichi-portal.service"}
STATUS_READ_ONLY_VERBS = {"is-active", "is-enabled", "is-failed", "show", "status", "cat",
                          "list-units", "list-unit-files", "list-timers", "list-dependencies", "list-sockets"}

# `systemctl` counts only where a command can start: line start, after ; & | ( $( `
# or `sudo`, or (so commented-out commands are caught in strict mode) after `#`.
SYSTEMCTL = re.compile(r"(?:^\s*|[;&|(`]\s*|\$\(\s*|#\s*|\bsudo\s+(?:-\S+\s+)*)systemctl\b(?P<rest>[^\n;&|#)`]*)")


def is_comment(line: str) -> bool:
    s = line.strip().lower()
    return s.startswith("#") or s.startswith("rem ") or s == "rem" or s.startswith("::")


def systemctl_findings(line: str, strict: bool) -> set:
    found = set()
    for m in SYSTEMCTL.finditer(line):
        tokens = [t for t in m.group("rest").split() if t and not t.startswith("-")]
        if not tokens:
            continue
        verb, units = tokens[0], set(tokens[1:])
        if strict:
            if verb not in STATUS_READ_ONLY_VERBS:
                found.add("systemctl-non-readonly")
            continue
        if verb not in MUTATING_SYSTEMCTL:
            continue
        if verb in PORTAL_ONLY_VERBS and units and units <= PORTAL_UNITS:
            continue
        found.add("systemctl-mutation")
    return found


def scan_text(text: str, *, include_comments: bool = False, strict: bool = False) -> list:
    """Return sorted (line_number, rule) findings. Matched text is never returned."""
    findings = set()
    for lineno, line in enumerate(text.splitlines(), 1):
        if not include_comments and is_comment(line):
            continue
        for rule, rx in UNMANAGED_RULES.items():
            if rx.search(line):
                findings.add((lineno, rule))
        for rule, rx in CREDENTIAL_RULES.items():
            for m in rx.finditer(line):
                value = m.groupdict().get("v") or ""
                if value and PLACEHOLDER.match(value):
                    continue
                findings.add((lineno, rule))
        for rule in systemctl_findings(line, strict):
            findings.add((lineno, rule))
        if strict:
            if re.search(r"\bsudo\b", line) and not is_comment(line):
                findings.add((lineno, "sudo"))
            if not is_comment(line) and re.search(r"(?:^|[;&|(]\s*|^\s*)(?:rm|mv|cp|chmod|chown|tee|dd|truncate|ln|mkdir|touch)\s", line):
                findings.add((lineno, "filesystem-write-command"))
            if not is_comment(line) and re.search(r"\bsed\s+-[A-Za-z]*i\b", line):
                findings.add((lineno, "in-place-edit"))
    return sorted(findings)


def tracked_launchers() -> list:
    try:
        out = subprocess.run(["git", "ls-files"], cwd=str(ROOT), capture_output=True, text=True, timeout=20)
        names = out.stdout.splitlines() if out.returncode == 0 else []
    except (OSError, subprocess.SubprocessError):
        names = []
    if not names:  # no git metadata (e.g. an exported tree): walk the tree instead
        skip = {".git", ".snakemake", "results", "resources", "__pycache__", "node_modules"}
        names = [str(p.relative_to(ROOT)).replace("\\", "/") for p in ROOT.rglob("*")
                 if p.is_file() and not (set(p.relative_to(ROOT).parts) & skip)]
    paths = [ROOT / n for n in names if n.lower().endswith(LAUNCHER_SUFFIXES)]
    return sorted(p for p in paths if p.is_file())


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT)).replace("\\", "/")


def pending_for(path: Path) -> set:
    return KNOWN_PENDING.get(rel(path), set())


class TestScannerBehaviour(unittest.TestCase):
    """The scanner itself: dangerous patterns are caught, legitimate commands are not."""

    def rules(self, text, **kw):
        return {rule for _line, rule in scan_text(text, **kw)}

    def test_detects_unmanaged_launch_patterns(self):
        cases = {
            "pkill -f app": "pkill-killall",
            "killall python": "pkill-killall",
            "kill -9 $PID": "kill-pid",
            "nohup ./run.sh &": "nohup",
            "cd portal && python app.py": "python-app",
            "/opt/env/bin/python3 app.py": "python-app",
            "flask --app app run": "flask-run",
            "cloudflared tunnel run genrichi-portal": "cloudflared-tunnel-run",
            "cloudflared --config x.yml tunnel run": "cloudflared-tunnel-run",
        }
        for text, rule in cases.items():
            with self.subTest(text=text):
                self.assertIn(rule, self.rules(text))

    def test_detects_unit_mutation_except_portal_enable_restart(self):
        for text in ("sudo systemctl start genrichi.service", "systemctl stop cloudflared", "systemctl disable genrichi",
                     "systemctl mask genrichi.service", "systemctl restart cloudflared.service",
                     "sudo systemctl enable genrichi.service", "systemctl restart genrichi-portal genrichi.service"):
            with self.subTest(text=text):
                self.assertIn("systemctl-mutation", self.rules(text))

    def test_allows_legitimate_deploy_and_status_commands(self):
        for text in ("sudo systemctl enable genrichi-portal --quiet", "sudo systemctl daemon-reload",
                     "sudo systemctl restart genrichi-portal.service", "systemctl is-active genrichi-portal.service",
                     "systemctl show cloudflared.service -p MainPID --value", "systemctl status genrichi.service --no-pager",
                     "conda activate snakemake", "ss -ltnp | grep :5000", "curl -s http://127.0.0.1:5000/login",
                     "snakemake --snakefile workflow/Snakefile --cores 4", "wget -q -c https://example.org/x.gz"):
            with self.subTest(text=text):
                self.assertEqual(scan_text(text), [])

    def test_detects_credential_like_literals_without_returning_them(self):
        secret_value = "Zx9-not-a-real-value"
        cases = {
            f'echo "  Password : {secret_value}"': "password-literal",
            f"echo Login: admin / {secret_value}": "login-echo",
            f"export PORTAL_PASS={secret_value}": "secret-assignment",
            f"set SMTP_PASS={secret_value}": "secret-assignment",
            f"tool --token {secret_value}": "secret-cli-arg",
            "curl -H x eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9abcdefgh": "token-shaped",
        }
        for text, rule in cases.items():
            with self.subTest(rule=rule):
                findings = scan_text(text)
                self.assertIn(rule, {r for _l, r in findings})
                self.assertNotIn(secret_value, repr(findings))

    def test_ignores_placeholders_variables_and_comments(self):
        for text in ('PORTAL_PASS=$PORTAL_PASS', 'export SECRET_KEY="${SECRET_KEY}"', "SMTP_PASS=<set-in-env>",
                     'echo "Password : <your password>"', "PASS=changeme",
                     "# pkill -f app  (documentation only)", "REM nohup python app.py", ":: cloudflared tunnel run x"):
            with self.subTest(text=text):
                self.assertEqual(scan_text(text), [])

    def test_strict_mode_rejects_any_state_changing_command_even_in_comments(self):
        for text in ("# systemctl restart genrichi-portal", "systemctl enable genrichi-portal", "systemctl daemon-reload",
                     "# pkill python", "sudo true", "rm -f /tmp/x", "sed -i s/a/b/ file"):
            with self.subTest(text=text):
                self.assertTrue(scan_text(text, include_comments=True, strict=True))


class TestTrackedLaunchers(unittest.TestCase):
    def test_launcher_discovery_finds_the_deployment_scripts(self):
        names = {rel(p) for p in tracked_launchers()}
        self.assertIn("deploy_portal.sh", names)

    def test_no_unmanaged_launch_patterns_in_tracked_launchers(self):
        unmanaged = set(UNMANAGED_RULES) | {"systemctl-mutation"}
        offenders = []
        for path in tracked_launchers():
            text = path.read_text(encoding="utf-8", errors="replace")
            for lineno, rule in scan_text(text):
                if rule in unmanaged and rule not in pending_for(path):
                    offenders.append(f"{rel(path)}:{lineno}  {rule}")
        self.assertEqual(offenders, [], "unmanaged launch path in tracked launcher(s):\n  " + "\n  ".join(offenders))

    def test_no_credential_like_literals_in_tracked_launchers(self):
        offenders = []
        for path in tracked_launchers():
            text = path.read_text(encoding="utf-8", errors="replace")
            for lineno, rule in scan_text(text):
                if rule in CREDENTIAL_RULES and rule not in pending_for(path):
                    offenders.append(f"{rel(path)}:{lineno}  {rule}")
        self.assertEqual(offenders, [], "credential-like literal in tracked launcher(s):\n  " + "\n  ".join(offenders))

    def test_pending_allowlist_is_not_stale(self):
        stale = []
        for name, rules in KNOWN_PENDING.items():
            path = ROOT / name
            self.assertTrue(path.is_file(), f"allowlisted file {name} no longer exists; remove it from KNOWN_PENDING")
            present = {rule for _l, rule in scan_text(path.read_text(encoding="utf-8", errors="replace"))}
            stale.extend(f"{name}:{rule}" for rule in sorted(rules - present))
        self.assertEqual(stale, [], "these KNOWN_PENDING entries are fixed; delete them from the allowlist:\n  " + "\n  ".join(stale))


class TestRetiredLauncher(unittest.TestCase):
    def test_start_genrichi_script_is_gone(self):
        self.assertFalse(REMOVED_LAUNCHER.exists(), "start_genrichi.sh must not exist (unmanaged launch path)")

    def test_nothing_tracked_references_the_retired_script(self):
        offenders = []
        for path in tracked_launchers() + [p for p in (ROOT / "tests").glob("*.py") if p.name != Path(__file__).name]:
            if re.search(r"start_genrichi", path.read_text(encoding="utf-8", errors="replace")):
                offenders.append(rel(path))
        self.assertEqual(offenders, [], "still references the retired launcher")


class TestPortalStatusScript(unittest.TestCase):
    def setUp(self):
        self.assertTrue(STATUS_SCRIPT.is_file(), "scripts/portal_status.sh is missing")
        self.text = STATUS_SCRIPT.read_text(encoding="utf-8")

    def test_is_strictly_read_only_including_comments(self):
        findings = scan_text(self.text, include_comments=True, strict=True)
        self.assertEqual(findings, [], f"portal_status.sh must be read-only: {findings}")

    def test_only_queries_the_expected_units_and_port(self):
        for expected in ("genrichi-portal.service", "genrichi.service", "cloudflared.service", "5000", "/login"):
            self.assertIn(expected, self.text)
        used = {m.group("rest").split()[0] for m in SYSTEMCTL.finditer(self.text) if m.group("rest").split()}
        self.assertTrue(used <= STATUS_READ_ONLY_VERBS, f"unexpected systemctl verbs: {used - STATUS_READ_ONLY_VERBS}")

    @unittest.skipUnless(shutil.which("bash"), "bash not available")
    def test_bash_syntax_is_valid(self):
        # feed the raw bytes on stdin: path-independent (Windows paths confuse Git Bash / WSL)
        # and no newline translation on Windows
        result = subprocess.run(["bash", "-n"], input=STATUS_SCRIPT.read_bytes(), capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr.decode("utf-8", "replace"))

    def test_uses_lf_line_endings(self):
        # a CRLF checkout (Windows autocrlf) would make bash fail with "$'\r': command not found"
        self.assertNotIn(b"\r", STATUS_SCRIPT.read_bytes(), "portal_status.sh must use LF line endings")

    def run_with_args(self, *args):
        """Run the script with args (fed on stdin: path independent). A hang is a test failure."""
        try:
            return subprocess.run(["bash", "-s", "--", *args], input=STATUS_SCRIPT.read_bytes(),
                                  capture_output=True, timeout=15)
        except subprocess.TimeoutExpired:
            self.fail(f"portal_status.sh {' '.join(args)!r} did not return: argument handling loops")

    @unittest.skipUnless(shutil.which("bash"), "bash not available")
    def test_wait_without_a_value_fails_immediately(self):
        result = self.run_with_args("--wait")
        self.assertEqual(result.returncode, 2, result.stderr)
        self.assertIn(b"--wait needs a whole number of seconds", result.stderr)
        self.assertEqual(result.stdout, b"", "must fail before printing any status")

    @unittest.skipUnless(shutil.which("bash"), "bash not available")
    def test_wait_with_a_non_numeric_value_fails_immediately(self):
        for bad in ("abc", "", "-5", "1.5", "10s", "--no-public"):
            with self.subTest(value=bad):
                result = self.run_with_args("--wait", bad)
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertIn(b"--wait needs a whole number of seconds", result.stderr)
                self.assertEqual(result.stdout, b"")

    @unittest.skipUnless(shutil.which("bash"), "bash not available")
    def test_wait_error_also_applies_after_other_options(self):
        result = self.run_with_args("--no-public", "--wait")
        self.assertEqual(result.returncode, 2, result.stderr)
        self.assertIn(b"--wait needs a whole number of seconds", result.stderr)

    @unittest.skipUnless(shutil.which("bash"), "bash not available")
    def test_unknown_option_is_a_usage_error(self):
        self.assertEqual(self.run_with_args("--bogus").returncode, 2)

    def test_wait_is_validated_before_any_query_runs(self):
        parse_end = self.text.index('case "$WAIT_SECONDS" in')
        self.assertLess(parse_end, self.text.index("systemctl is-active"))
        self.assertLess(parse_end, self.text.index("curl "))


class TestLocalWindowsLauncher(unittest.TestCase):
    """GenRichi_Start.bat is gitignored and local; checked only when present."""

    def setUp(self):
        if not LOCAL_BAT.is_file():
            self.skipTest("GenRichi_Start.bat not present on this machine")
        self.text = LOCAL_BAT.read_text(encoding="utf-8", errors="replace")

    def test_calls_only_the_read_only_status_script(self):
        self.assertIn("portal_status.sh", self.text)
        self.assertNotIn("start_genrichi", self.text)

    def test_has_no_unmanaged_launch_pattern_or_credential(self):
        findings = scan_text(self.text, include_comments=True)
        self.assertEqual(findings, [], f"GenRichi_Start.bat findings (line, rule): {findings}")

    def test_does_not_manage_services(self):
        self.assertNotRegex(self.text, r"(?i)\bsystemctl\b")


if __name__ == "__main__":
    unittest.main()
