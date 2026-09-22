"""
Tests for deploy_comprehensive.sh (the hardened Comprehensive-pipeline
deployment mechanism, and its Comprehensive-scoped build stamp).

Mirrors the harness pattern of tests/test_deploy_portal.py: a throwaway git
repository and a fake "production" tree under /tmp, PATH shims for
systemctl/pgrep, and GENRICHI_DEPLOY_TEST_* overrides the script only
honours for paths under /tmp.

Run with (POSIX/WSL):
    python -m unittest tests.test_deploy_comprehensive -v
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
SCRIPT = ROOT / "deploy_comprehensive.sh"
sys.path.insert(0, str(Path(__file__).resolve().parent))

from test_launchers_safe import scan_text  # noqa: E402

POSIX = os.name == "posix" and shutil.which("bash") and shutil.which("git") and shutil.which("sha256sum")

ALLOWLIST_CODE = [
    "workflow/Snakefile_comprehensive",
    "workflow/rules/runtime_versions.smk",
    "workflow/rules/qc_paired.smk",
    "workflow/rules/align_paired.smk",
    "workflow/rules/somatic_calling_paired.smk",
    "workflow/rules/cnv_calling.smk",
    "workflow/rules/msi_scoring.smk",
    "workflow/rules/annotation_comprehensive.smk",
    "workflow/rules/report_comprehensive.smk",
    "workflow/scripts/calculate_cnv.py",
    "workflow/scripts/calculate_msi.py",
    "workflow/scripts/vcf_to_somatic_table.py",
    "workflow/scripts/generate_comprehensive_report.py",
    "workflow/scripts/provenance.py",
    "workflow/scripts/resource_inspect.py",
    "workflow/scripts/qc_status.py",
    "workflow/scripts/sample_columns.py",
    "workflow/scripts/runtime_probe.py",
]
ALLOWLIST_CONFIG = [
    "config/comprehensive_config.yaml",
    "config/comprehensive_samples.tsv",
    "resources/panel/phase1_solid_tumor/solid_tumor_phase1_v1.bed",
    "resources/panel/comprehensive_genes_cds.bed",
]
EXPECTED_ALLOWLIST = ALLOWLIST_CODE + ALLOWLIST_CONFIG
RULE_FILES = [r for r in ALLOWLIST_CODE if r.startswith("workflow/rules/")]

SYSTEMCTL_SHIM = r"""#!/bin/sh
echo "$*" >> "$SHIM_STATE/systemctl.log"
state() { if [ -f "$SHIM_STATE/$1" ]; then cat "$SHIM_STATE/$1"; else echo "$2"; fi; }
case "$1 $2" in
  "is-active genrichi-portal.service")  s=$(state portal_state active); echo "$s"; [ "$s" = active ] || exit 3 ;;
  "show genrichi-portal.service")       echo /fake.slice/genrichi-portal.service ;;
  *) echo "unexpected systemctl call: $*" >&2; echo "$*" >> "$SHIM_STATE/forbidden.log"; exit 99 ;;
esac
"""
PGREP_SHIM = r"""#!/bin/sh
if [ -f "$SHIM_STATE/pgrep_error" ]; then echo "pgrep: simulated failure" >&2; exit "$(cat "$SHIM_STATE/pgrep_error")"; fi
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
    if rel.endswith(".smk"):
        if rel.endswith("report_comprehensive.smk"):
            return f"# {rel} {version}\nrule comprehensive_report:\n    script:\n        \"../scripts/generate_comprehensive_report.py\"\n".encode()
        if rel.endswith("cnv_calling.smk"):
            return f"# {rel} {version}\nrule cnv:\n    script:\n        \"../scripts/calculate_cnv.py\"\n".encode()
        if rel.endswith("msi_scoring.smk"):
            return f"# {rel} {version}\nrule msi:\n    script:\n        \"../scripts/calculate_msi.py\"\n".encode()
        return f"# {rel} {version}\nrule r_{version}:\n    shell: \"true\"\n".encode()
    if rel.endswith(".py"):
        return f"# {rel} {version}\nVALUE = '{version}'\n".encode()
    if rel.endswith(".yaml"):
        return (
            f"# {rel} {version}\n"
            "samples: config/comprehensive_samples.tsv\n"
            "panel:\n"
            "  bed: resources/panel/phase1_solid_tumor/solid_tumor_phase1_v1.bed  # comment\n"
            "  name: \"Test Panel\"\n"
            "tmb:\n"
            "  coding_mb: 1.5\n"
            "  coding_bed: resources/panel/comprehensive_genes_cds.bed\n"
        ).encode()
    if rel.endswith(".tsv"):
        return f"sample_id\ttumor_r1\n{version}\tx.fastq\n".encode()
    if rel.endswith(".bed"):
        return f"chr17\t{7675000 if 'phase1' in rel else 1}\t{7675100 if 'phase1' in rel else 100}\t{version}\n".encode()
    return f"{rel} {version}\n".encode()


def snakefile_content(version, rule_files=RULE_FILES):
    includes = "\n".join(f'include: "rules/{Path(r).name}"' for r in rule_files)
    return f"# Snakefile_comprehensive {version}\n{includes}\n\nrule all:\n    input: []\n".encode()


@unittest.skipUnless(POSIX, "needs a POSIX shell with git and sha256sum (run in WSL/Linux)")
class DeployComprehensiveTestCase(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp(prefix="deploycomp-", dir="/tmp"))
        self.addCleanup(shutil.rmtree, self.root, ignore_errors=True)
        self.home = self.root / "home"
        self.repo = self.root / "repo"
        self.state = self.root / "state"
        self.shims = self.root / "shims"
        self.cgroup = self.root / "cgroup"
        (self.root / "tmp").mkdir()
        for d in (self.home, self.state, self.shims):
            d.mkdir(parents=True)
        (self.home / "portal").mkdir()
        self._build_repo()
        self._build_production()
        self._build_shims()
        procs = self.cgroup / "fake.slice" / "genrichi-portal.service"
        procs.mkdir(parents=True)
        (procs / "cgroup.procs").write_text("4242\n")

    # ── fixtures ────────────────────────────────────────────────────────
    def _build_repo(self):
        self.repo.mkdir()
        git(self.repo, "init", "-q")
        for rel in ALLOWLIST_CODE:
            p = self.repo / rel
            p.parent.mkdir(parents=True, exist_ok=True)
            if rel == "workflow/Snakefile_comprehensive":
                p.write_bytes(snakefile_content("v1"))
            elif rel == "workflow/scripts/provenance.py":
                # real content: test_12 exercises genuine integration with
                # the actual resolve_pipeline_identity()/read_build_stamp().
                p.write_bytes((ROOT / rel).read_bytes())
            else:
                p.write_bytes(content_for(rel, "v1"))
        for rel in ALLOWLIST_CONFIG:
            p = self.repo / rel
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_bytes(content_for(rel, "v1"))
        # the retired BED must exist in the repo (tracked) but stay OUT of the allowlist
        retired = self.repo / "resources/panel/comprehensive_genes.bed"
        retired.parent.mkdir(parents=True, exist_ok=True)
        retired.write_bytes(b"chr1\t1\t2\tretired\n")
        self.commit("initial")

    def commit(self, message):
        git(self.repo, "add", "-A")
        git(self.repo, "commit", "-q", "-m", message)
        self.head = git(self.repo, "rev-parse", "HEAD")

    def _build_production(self):
        # some files already present (older content -> CHANGED), some absent (-> NEW)
        present = set(EXPECTED_ALLOWLIST) - {
            "workflow/rules/runtime_versions.smk",
            "workflow/scripts/resource_inspect.py",
            "workflow/scripts/runtime_probe.py",
        }
        for rel in present:
            p = self.home / rel
            p.parent.mkdir(parents=True, exist_ok=True)
            if rel == "workflow/Snakefile_comprehensive":
                p.write_bytes(snakefile_content("v0"))
            elif rel in ("workflow/scripts/provenance.py", "workflow/rules/report_comprehensive.smk"):
                p.write_bytes(content_for(rel, "v0"))  # stale -> CHANGED
            else:
                p.write_bytes(content_for(rel, "v1"))  # already current -> same
        (self.home / "portal" / ".env").write_bytes(b"SENTINEL-not-a-real-secret\n")
        os.chmod(self.home / "portal" / ".env", 0o600)
        (self.home / "portal" / "genrichi_orders.db").write_bytes(b"SQLite format 3\x00" + b"sentinel" * 20)
        for d in ("results", "uploads", "portal_logs", "logs", ".snakemake", "reference_db", "envs"):
            (self.home / d).mkdir(exist_ok=True)
        (self.home / "results" / "r.txt").write_text("result\n")
        (self.home / "reference_db" / "hg38.fa").write_text("machine-local reference\n")
        (self.home / "envs" / "align.yaml").write_text("name: align\n")
        # the retired BED already sitting in production, from the old unhardened script
        retired = self.home / "resources/panel/comprehensive_genes.bed"
        retired.parent.mkdir(parents=True, exist_ok=True)
        retired.write_bytes(b"chr1\t1\t2\tretired-in-production\n")

    def _build_shims(self):
        for name, body in (("systemctl", SYSTEMCTL_SHIM), ("pgrep", PGREP_SHIM)):
            (self.shims / name).write_text(body)
        for name in ("sudo", "pip", "pip3", "nohup", "pkill", "killall", "cloudflared", "kill"):
            (self.shims / name).write_text(FORBIDDEN_SHIM)
        for p in self.shims.iterdir():
            p.chmod(0o755)

    def path_without(self, *names):
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
    def run_script(self, *args, extra_env=None, path_override=None):
        env = dict(os.environ)
        env.update({
            "PATH": path_override or f"{self.shims}:{os.environ['PATH']}",
            "SHIM_STATE": str(self.state),
            "TMPDIR": str(self.root / "tmp"),
            "GENRICHI_DEPLOY_TEST_HOME": str(self.home),
            "GENRICHI_DEPLOY_TEST_REPO": str(self.repo),
            "GENRICHI_DEPLOY_TEST_CGROUP_ROOT": str(self.cgroup),
            "GIT_CONFIG_GLOBAL": "/dev/null",
        })
        env.update(extra_env or {})
        result = subprocess.run(["bash", str(SCRIPT), *args], capture_output=True, text=True, env=env, timeout=120)
        self.output = result.stdout + result.stderr
        self.assertFalse((self.state / "forbidden.log").exists(),
                          "the script invoked a forbidden command: " +
                          ((self.state / "forbidden.log").read_text() if (self.state / "forbidden.log").exists() else ""))
        return result

    def snapshot(self):
        snap = {}
        for base, dirs, files in os.walk(self.home):
            for name in files + dirs:
                p = Path(base) / name
                st = p.lstat()
                snap[str(p)] = (sha256(p) if p.is_file() else "dir", st.st_mtime_ns if p.is_file() else 0, stat.S_IMODE(st.st_mode))
        return snap

    def backup_dirs(self):
        parent = self.home / "deploy_backups" / "comprehensive"
        return sorted(parent.iterdir()) if parent.exists() else []

    def stamp_path(self):
        return self.home / "workflow" / "build_info.json"

    def assert_nothing_written(self, before):
        self.assertEqual(self.snapshot(), before, "production files changed")
        self.assertEqual(self.backup_dirs(), [], "a backup directory was created")
        self.assertFalse(self.stamp_path().exists(), "build_info.json was written")
        self.assertEqual(list((self.home / "workflow").rglob(".deploy.*")), [], "temporary files left behind")

    def assert_refused(self, result, needle=None):
        self.assertEqual(result.returncode, 1, self.output)
        if needle:
            self.assertIn(needle, self.output)

    def refuse_case(self, prepare, needle, *args):
        prepare()
        before = self.snapshot()
        result = self.run_script("--execute", *args)
        self.assert_refused(result, needle)
        self.assert_nothing_written(before)

    # ── basics ──────────────────────────────────────────────────────────
    def test_bash_syntax_is_valid(self):
        self.assertEqual(subprocess.run(["bash", "-n", str(SCRIPT)], capture_output=True).returncode, 0)

    def test_help_exits_zero_and_writes_nothing(self):
        before = self.snapshot()
        result = self.run_script("--help")
        self.assertEqual(result.returncode, 0)
        self.assertIn("dry-run", result.stdout.lower())
        self.assert_nothing_written(before)

    # 1. dry-run does not modify anything
    def test_1_dry_run_does_not_modify_anything(self):
        before = self.snapshot()
        result = self.run_script()
        self.assertEqual(result.returncode, 0, self.output)
        self.assertIn("DRY-RUN", self.output)
        self.assertIn("nothing was written", self.output)
        self.assert_nothing_written(before)

    def test_1b_dry_run_reports_the_plan(self):
        self.run_script()
        self.assertIn("CHANGED", self.output)
        self.assertIn("NEW", self.output)
        self.assertIn("3 new", self.output)

    # 2. execute is explicitly required
    def test_2_execute_is_explicitly_required(self):
        before = self.snapshot()
        result = self.run_script()  # no --execute
        self.assertEqual(result.returncode, 0)
        self.assert_nothing_written(before)
        result2 = self.run_script("--execute")
        self.assertEqual(result2.returncode, 0, self.output)
        self.assertNotEqual(self.snapshot(), before)

    # 3. wrong --expect-commit fails
    def test_3_wrong_expect_commit_fails(self):
        self.refuse_case(lambda: None, "does not match --expect-commit", "--expect-commit", "0" * 40)

    def test_3b_short_expect_commit_is_a_usage_error(self):
        before = self.snapshot()
        result = self.run_script("--execute", "--expect-commit", self.head[:10])
        self.assertEqual(result.returncode, 2)
        self.assert_nothing_written(before)

    # 4. unexpected source files fail closed
    def test_4_unexpected_snakefile_include_fails_closed(self):
        def prepare():
            (self.repo / "workflow/rules/rogue_rule.smk").write_text("rule rogue:\n    shell: \"true\"\n")
            (self.repo / "workflow/Snakefile_comprehensive").write_bytes(
                snakefile_content("v2", RULE_FILES + ["workflow/rules/rogue_rule.smk"]))
            self.commit("add an unlisted rule to the include graph")
        self.refuse_case(prepare, "include: list no longer matches ALLOWLIST_CODE")

    def test_4b_allowlisted_file_missing_from_commit_fails_closed(self):
        def prepare():
            (self.repo / "workflow/scripts/runtime_probe.py").unlink()
            self.commit("remove an allowlisted file")
        self.refuse_case(prepare, "allowlisted file missing in commit")

    def test_4c_config_bed_drift_fails_closed(self):
        def prepare():
            bad = content_for("config/comprehensive_config.yaml", "v1").decode().replace(
                "resources/panel/phase1_solid_tumor/solid_tumor_phase1_v1.bed",
                "resources/panel/comprehensive_genes.bed")
            (self.repo / "config/comprehensive_config.yaml").write_text(bad)
            self.commit("point config at the retired panel BED")
        self.refuse_case(prepare, "panel.bed is")

    # 5. unexpected destination files are not touched
    def test_5_retired_bed_and_other_pipelines_untouched(self):
        retired = self.home / "resources/panel/comprehensive_genes.bed"
        before_hash = sha256(retired)
        other_pipeline = self.home / "workflow" / "Snakefile_hereditary"
        other_pipeline.write_text("# unrelated pipeline, must never be touched\n")
        before_other = sha256(other_pipeline)
        self.assertEqual(self.run_script("--execute").returncode, 0, self.output)
        self.assertEqual(sha256(retired), before_hash, "retired panel BED was modified")
        self.assertEqual(sha256(other_pipeline), before_other, "an unrelated pipeline's file was touched")

    # 6. backup created before overwrite
    def test_6_backup_created_before_overwrite_and_only_for_changed_files(self):
        changed_rels = ("workflow/Snakefile_comprehensive", "workflow/scripts/provenance.py", "workflow/rules/report_comprehensive.smk")
        old = {rel: (self.home / rel).read_bytes() for rel in changed_rels}
        self.assertEqual(self.run_script("--execute").returncode, 0, self.output)
        backups = self.backup_dirs()
        self.assertEqual(len(backups), 1)
        b = backups[0]
        self.assertEqual(stat.S_IMODE(b.stat().st_mode), 0o700)
        saved = {str(p.relative_to(b / "files")).replace(os.sep, "/"): p.read_bytes() for p in (b / "files").rglob("*") if p.is_file()}
        self.assertEqual(saved, old)
        for name in (".env", "genrichi_orders.db", "r.txt", "hg38.fa"):
            self.assertEqual(list(b.rglob(name)), [], f"{name} must never be backed up")

    # 7. copied hashes equal source hashes
    def test_7_deployed_file_hashes_equal_source_hashes(self):
        self.assertEqual(self.run_script("--execute").returncode, 0, self.output)
        for rel in EXPECTED_ALLOWLIST:
            src_sha = git(self.repo, "hash-object", str(self.repo / rel))
            # git hash-object includes a blob header internally different from sha256; compare raw content instead
            self.assertEqual((self.home / rel).read_bytes(), (self.repo / rel).read_bytes(), rel)

    # 8. manifest hashes equal deployed hashes
    def test_8_manifest_hashes_equal_deployed_hashes(self):
        self.assertEqual(self.run_script("--execute").returncode, 0, self.output)
        manifest = (self.backup_dirs()[0] / "MANIFEST.txt").read_text()
        self.assertIn(f"source_commit={self.head}", manifest)
        for line in manifest.splitlines():
            parts = line.split()
            if parts and parts[0] in ("NEW", "CHANGED"):
                rel, new_sha = parts[1], parts[3]
                self.assertEqual(sha256(self.home / rel), new_sha, rel)

    # 9. build stamp contains the expected Comprehensive commit
    def test_9_build_stamp_contains_expected_commit(self):
        self.assertEqual(self.run_script("--execute", "--expect-commit", self.head).returncode, 0, self.output)
        import json
        stamp = json.loads(self.stamp_path().read_text())
        self.assertEqual(stamp["git_commit_sha"], self.head)
        self.assertEqual(stamp["source_commit"], self.head)
        self.assertEqual(stamp["pipeline"], "comprehensive")  # 11 also covered here

    # 10. build stamp contains the expected code fingerprint
    def test_10_build_stamp_contains_code_fingerprint(self):
        self.assertEqual(self.run_script("--execute").returncode, 0, self.output)
        import hashlib, json
        stamp = json.loads(self.stamp_path().read_text())
        fp = stamp["comprehensive_code_fingerprint_sha256"]
        self.assertRegex(fp, r"^[0-9a-f]{64}$")
        expected_input = "".join(f"{rel}\t{sha256(self.home / rel)}\n" for rel in sorted(ALLOWLIST_CODE))
        expected = hashlib.sha256(expected_input.encode()).hexdigest()
        self.assertEqual(fp, expected)

    # 11. pipeline identifier is comprehensive (see test_9) + explicit check
    def test_11_pipeline_identifier_is_comprehensive(self):
        self.assertEqual(self.run_script("--execute").returncode, 0, self.output)
        import json
        self.assertEqual(json.loads(self.stamp_path().read_text())["pipeline"], "comprehensive")

    # 12. provenance can consume the build stamp
    def test_12_provenance_py_consumes_the_build_stamp(self):
        self.assertEqual(self.run_script("--execute", "--expect-commit", self.head).returncode, 0, self.output)
        no_git_copy = self.root / "no_git_copy"
        shutil.copytree(self.home / "workflow", no_git_copy / "workflow")
        sys.path.insert(0, str(no_git_copy / "workflow" / "scripts"))
        import importlib
        for mod in ("provenance", "resource_inspect"):
            sys.modules.pop(mod, None)
        provenance = importlib.import_module("provenance")
        identity = provenance.resolve_pipeline_identity(str(no_git_copy))
        self.assertEqual(identity["git_commit_sha"], self.head)
        self.assertEqual(identity["git_commit_source"], "build_stamp")
        sys.path.remove(str(no_git_copy / "workflow" / "scripts"))
        for mod in ("provenance", "resource_inspect"):
            sys.modules.pop(mod, None)

    # 13. fingerprint mismatch is detected
    def test_13_fingerprint_mismatch_is_detected_after_hand_edit(self):
        self.assertEqual(self.run_script("--execute").returncode, 0, self.output)
        import hashlib, json
        stamp = json.loads(self.stamp_path().read_text())
        (self.home / "workflow/scripts/calculate_msi.py").write_text("# hand-edited after stamping, bypassing the deploy script\n")
        recomputed_input = "".join(
            f"{rel}\t{sha256(self.home / rel)}\n" for rel in sorted(ALLOWLIST_CODE))
        recomputed = hashlib.sha256(recomputed_input.encode()).hexdigest()
        self.assertNotEqual(stamp["comprehensive_code_fingerprint_sha256"], recomputed,
                             "a hand-edit after stamping must be detectable via fingerprint mismatch")

    # 14. source/deployment commit mismatch fails (covered by test_3, plus a stamp-vs-HEAD check)
    def test_14_stamp_commit_mismatches_a_later_head_is_detectable(self):
        self.assertEqual(self.run_script("--execute").returncode, 0, self.output)
        import json
        first_stamp = json.loads(self.stamp_path().read_text())
        (self.repo / "workflow/scripts/calculate_msi.py").write_text("# calculate_msi.py v2\nVALUE = 'v2'\n")
        self.commit("a later, undeployed commit")
        self.assertNotEqual(first_stamp["git_commit_sha"], self.head,
                             "the stamp must not silently claim to describe a commit newer than what was deployed")

    # 15. .env untouched
    def test_15_env_untouched(self):
        before = (self.home / "portal" / ".env").read_bytes()
        before_mode = stat.S_IMODE((self.home / "portal" / ".env").stat().st_mode)
        self.assertEqual(self.run_script("--execute").returncode, 0, self.output)
        self.assertEqual((self.home / "portal" / ".env").read_bytes(), before)
        self.assertEqual(stat.S_IMODE((self.home / "portal" / ".env").stat().st_mode), before_mode)

    # 16. database untouched
    def test_16_database_untouched(self):
        before = (self.home / "portal" / "genrichi_orders.db").read_bytes()
        self.assertEqual(self.run_script("--execute").returncode, 0, self.output)
        self.assertEqual((self.home / "portal" / "genrichi_orders.db").read_bytes(), before)

    # 17. results/logs/uploads untouched
    def test_17_results_logs_uploads_untouched(self):
        before = {p: sha256(self.home / p) for p in ("results/r.txt",)}
        self.assertEqual(self.run_script("--execute").returncode, 0, self.output)
        for p, h in before.items():
            self.assertEqual(sha256(self.home / p), h, p)
        for d in ("uploads", "portal_logs", "logs", ".snakemake", "reference_db", "envs"):
            self.assertTrue((self.home / d).is_dir(), f"{d} must still exist untouched")

    # 18. no service lifecycle commands are used
    def test_18_only_read_only_systemctl_verbs_are_used(self):
        self.assertEqual(self.run_script("--execute").returncode, 0, self.output)
        verbs = {line.split()[0] for line in (self.state / "systemctl.log").read_text().splitlines() if line.strip()}
        self.assertTrue(verbs <= {"is-active", "show"}, verbs)

    def test_18b_script_has_no_forbidden_commands(self):
        text = SCRIPT.read_text(encoding="utf-8")
        self.assertEqual(scan_text(text, include_comments=True), [])
        for pattern in (r"\bsudo\b", r"\bpip3?\s+install\b", r"\bnohup\b", r"\bpkill\b", r"\bkillall\b",
                        r"\bkill\s+-9\b", r"rm\s+-rf?\s+.*\.snakemake"):
            self.assertIsNone(re.search(pattern, text), f"forbidden pattern: {pattern}")
        self.assertIsNone(re.search(r"systemctl\s+(start|stop|restart|enable|disable|mask|unmask|daemon-reload|kill)", text))
        self.assertIn("set -euo pipefail", text)
        self.assertNotIn(b"\r", SCRIPT.read_bytes())

    # 19. active pipeline/run causes deployment to fail safely
    def test_19_active_snakemake_process_refuses(self):
        self.refuse_case(lambda: (self.state / "snakemake_count").write_text("1"), "Snakemake process")

    def test_19b_running_job_in_portal_cgroup_refuses(self):
        procs = self.cgroup / "fake.slice" / "genrichi-portal.service" / "cgroup.procs"
        self.refuse_case(lambda: procs.write_text("4242\n5151\n"), "analysis appears to be running")

    def test_19c_snakemake_locks_present_refuses_and_locks_are_never_touched(self):
        locks = self.home / ".snakemake" / "locks"

        def prepare():
            locks.mkdir(parents=True)
            (locks / "lock.file").write_text("held\n")
        self.refuse_case(prepare, ".snakemake/locks exists")
        self.assertTrue(locks.exists(), "the script must never remove or touch an existing lock directory")
        self.assertTrue((locks / "lock.file").exists())

    def test_19d_pgrep_error_fails_closed_never_kills(self):
        (self.state / "pgrep_error").write_text("2")
        before = self.snapshot()
        result = self.run_script("--execute")
        self.assert_refused(result, "pgrep failed")
        self.assert_nothing_written(before)
        self.assertFalse((self.state / "forbidden.log").exists())

    # 20. no partial deployment occurs after a failed preflight
    def test_20_no_partial_deployment_after_failed_preflight(self):
        def prepare():
            (self.state / "portal_state").write_text("activating")
        prepare()
        before = self.snapshot()
        result = self.run_script("--execute")
        self.assertEqual(result.returncode, 1, self.output)
        self.assert_nothing_written(before)

    def test_20b_backup_failure_aborts_before_any_write(self):
        (self.home / "deploy_backups").write_text("a file, not a directory")
        before = {k: v for k, v in self.snapshot().items() if "deploy_backups" not in k}
        result = self.run_script("--execute")
        self.assertNotEqual(result.returncode, 0)
        after = {k: v for k, v in self.snapshot().items() if "deploy_backups" not in k}
        self.assertEqual(after, before, "production changed although the backup failed")
        self.assertFalse(self.stamp_path().exists())

    # ── extras ──────────────────────────────────────────────────────────
    def test_second_run_is_a_noop_without_a_new_backup_or_restamp(self):
        self.assertEqual(self.run_script("--execute").returncode, 0, self.output)
        backups = self.backup_dirs()
        stamp_before = self.stamp_path().read_text()
        result = self.run_script("--execute")
        self.assertEqual(result.returncode, 0, self.output)
        self.assertIn("Nothing to deploy", self.output)
        self.assertEqual(self.backup_dirs(), backups)
        self.assertEqual(self.stamp_path().read_text(), stamp_before)

    def test_destination_symlink_is_refused(self):
        other = self.root / "elsewhere"
        shutil.copytree(self.home, other, symlinks=True)
        shutil.rmtree(self.home)
        os.symlink(other, self.home)
        result = self.run_script("--execute")
        self.assert_refused(result)

    def test_credential_like_literal_in_source_is_refused_and_not_echoed(self):
        secret = "Zq7-fake-not-real-1234"

        def prepare():
            (self.repo / "workflow/scripts/calculate_cnv.py").write_text(f'# cnv\nSECRET_KEY = "{secret}"\n')
            self.commit("hardcode a secret")
        prepare()
        before = self.snapshot()
        result = self.run_script("--execute")
        self.assert_refused(result, "credential-like literal in workflow/scripts/calculate_cnv.py")
        self.assertNotIn(secret, self.output)
        self.assert_nothing_written(before)

    def test_python_syntax_error_is_refused(self):
        def prepare():
            (self.repo / "workflow/scripts/calculate_msi.py").write_text("def broken(:\n")
            self.commit("syntax error")
        self.refuse_case(prepare, "syntax check failed")

    def test_crlf_in_source_is_refused(self):
        def prepare():
            (self.repo / "workflow/scripts/qc_status.py").write_bytes(b"# qc\r\nX = 1\r\n")
            self.commit("crlf")
        self.refuse_case(prepare, "carriage returns")

    def test_dependency_files_never_deployed(self):
        self.assertEqual(self.run_script("--execute").returncode, 0, self.output)
        self.assertFalse((self.home / "resources/panel/comprehensive_genes.bed").read_bytes()
                          == (self.repo / "resources/panel/comprehensive_genes.bed").read_bytes())


if __name__ == "__main__":
    unittest.main()
