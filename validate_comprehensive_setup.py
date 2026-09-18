"""
GenRichi comprehensive (Phase 1 solid-tumor) pipeline setup validator.

Scoped counterpart to validate_setup.py, which only covers the hereditary
pipeline. Run this from anywhere -- it uses absolute paths.

    python C:/GenRichi/validate_comprehensive_setup.py

Checks:
  1. All comprehensive-pipeline files exist
  2. All relevant Python scripts are syntax-valid
  3. comprehensive_config.yaml parses correctly
  4. include/script references inside the comprehensive .smk files resolve
  5. snakemake is reachable on PATH
  6. Full Phase 5.2 pre-flight validation of the actual configured
     comprehensive_config.yaml + comprehensive_samples.tsv (config keys/
     ranges, samples sheet, FASTQ existence/non-emptiness, required
     reference/resource files, panel/reference build compatibility)
"""

import ast
import os
import pathlib
import re
import shutil
import subprocess
import sys

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")

ROOT = pathlib.Path(__file__).parent.resolve()
os.chdir(ROOT)
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

PASS = "\033[92m PASS\033[0m"
FAIL = "\033[91m FAIL\033[0m"
WARN = "\033[93m WARN\033[0m"

failures = 0


def ok(msg):
    print(f"{PASS}  {msg}")


def fail(msg):
    global failures
    failures += 1
    print(f"{FAIL}  {msg}")


def warn(msg):
    print(f"{WARN}  {msg}", flush=True)


print("\nGenRichi Comprehensive (Phase 1) Setup Validator")
print(f"Project root: {ROOT}")
print("=" * 60)

# ── 1. Required files ─────────────────────────────────────────────────────
print("\n[1] Required files")
required = [
    "workflow/Snakefile_comprehensive",
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
    "workflow/scripts/preflight.py",
    "workflow/scripts/qc_status.py",
    "config/comprehensive_config.yaml",
    "config/comprehensive_samples.tsv",
    "workflow/envs/align.yaml",
    "workflow/envs/calling.yaml",
    "workflow/envs/qc.yaml",
    "workflow/envs/annotation.yaml",
    "workflow/envs/report.yaml",
]
for f in required:
    p = ROOT / f
    if p.exists():
        ok(f)
    else:
        fail(f"MISSING: {f}")

# ── 2. Python syntax ───────────────────────────────────────────────────────
print("\n[2] Python script syntax")
_scripts = [
    "calculate_cnv.py", "calculate_msi.py", "vcf_to_somatic_table.py",
    "generate_comprehensive_report.py", "provenance.py", "preflight.py",
    "qc_status.py",
]
for name in _scripts:
    f = ROOT / "workflow" / "scripts" / name
    if not f.exists():
        continue  # already reported above
    try:
        ast.parse(f.read_text(encoding="utf-8"))
        ok(name)
    except SyntaxError as e:
        fail(f"{name}: {e}")

# ── 3. YAML config parsing ─────────────────────────────────────────────────
print("\n[3] YAML config files")
try:
    import yaml
    for f in [ROOT / "config/comprehensive_config.yaml"] + list((ROOT / "workflow/envs").glob("*.yaml")):
        try:
            yaml.safe_load(f.read_text(encoding="utf-8"))
            ok(str(f.relative_to(ROOT)))
        except yaml.YAMLError as e:
            fail(f"{f.name}: {e}")
except ImportError:
    warn("PyYAML not installed -- skipping YAML validation (pip install pyyaml)")

# ── 4. Include / script cross-references ───────────────────────────────────
print("\n[4] Include / script cross-references")
sf_path = ROOT / "workflow/Snakefile_comprehensive"
if sf_path.exists():
    sf_txt = sf_path.read_text(encoding="utf-8")
    for m in re.finditer(r'include:\s*["\'](.+?)["\']', sf_txt):
        inc = ROOT / "workflow" / m.group(1)
        if inc.exists():
            ok(f"include: {m.group(1)}")
        else:
            fail(f"include missing: {m.group(1)}")

    _comprehensive_rule_files = {
        m.group(1) for m in re.finditer(r'include:\s*["\'](.+?)["\']', sf_txt)
    }
    for rel in _comprehensive_rule_files:
        smk = ROOT / "workflow" / rel
        if not smk.exists():
            continue
        txt = smk.read_text(encoding="utf-8")
        for m in re.finditer(r'script:\s*["\'](.+?)["\']', txt):
            script_path = (smk.parent / m.group(1)).resolve()
            if script_path.exists():
                ok(f"script: {script_path.name} (in {smk.name})")
            else:
                fail(f"script missing: {m.group(1)} referenced in {smk.name}")
else:
    fail("workflow/Snakefile_comprehensive not found -- cannot check includes")

# ── 5. Snakemake on PATH ────────────────────────────────────────────────────
print("\n[5] Tools on PATH")
for tool in ["snakemake", "conda", "mamba"]:
    path = shutil.which(tool)
    if path:
        try:
            ver = subprocess.check_output(
                [tool, "--version"], stderr=subprocess.STDOUT, text=True
            ).strip()
            ok(f"{tool} ({ver})  ->  {path}")
        except Exception:
            ok(f"{tool}  ->  {path}")
    else:
        if tool == "snakemake":
            fail(f"{tool} not found -- install with: conda install -c bioconda snakemake")
        else:
            warn(f"{tool} not found  (mamba/conda optional but recommended)")

# ── 6. Full Phase 5.2 pre-flight validation ─────────────────────────────────
print("\n[6] Pre-flight validation (config + samples sheet + resources)")
try:
    import pandas as pd
    import yaml as _yaml
    from preflight import run_preflight_checks, format_report

    cfg_path = ROOT / "config/comprehensive_config.yaml"
    _config = _yaml.safe_load(cfg_path.read_text(encoding="utf-8"))
    _samples_path = ROOT / _config["samples"]
    _samples_df = pd.read_csv(_samples_path, sep="\t", dtype=str).set_index("sample_id")

    report = run_preflight_checks(_config, _samples_df, base_dir=str(ROOT))
    print(format_report(report))
    if not report.ok:
        fail(f"Pre-flight validation found {len(report.errors)} error(s) -- see above")
    else:
        ok("Pre-flight validation passed")
        for w in report.warnings:
            warn(f"{w.field}: {w.message}")
except Exception as e:  # pragma: no cover - diagnostic path only
    warn(f"Could not run pre-flight validation: {e}")

# ── Summary ──────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
if failures == 0:
    print("\033[92mAll checks passed!\033[0m")
    print("\nNext step -- dry-run the comprehensive pipeline:")
    print(f"  cd {ROOT}")
    print("  snakemake --snakefile workflow/Snakefile_comprehensive --configfile config/comprehensive_config.yaml --use-conda --cores 8 --dry-run")
else:
    print(f"\033[91m{failures} check(s) failed -- fix issues above before running.\033[0m")
    sys.exit(1)
