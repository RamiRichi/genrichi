"""
GenRichi pipeline setup validator.
Run this from anywhere — it uses absolute paths.

    python C:/GenRichi/validate_setup.py

Checks:
  1. All pipeline files exist
  2. All Python scripts are syntax-valid
  3. All YAML configs parse correctly
  4. All include/script references inside .smk files resolve
  5. snakemake is reachable on PATH
  6. Resource files (reference, annotation) exist (warnings only)
"""

import ast
import io
import os
import pathlib
import re
import shutil
import subprocess
import sys

# Force UTF-8 output on Windows (avoids CP1252 encode errors)
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")

ROOT = pathlib.Path(__file__).parent.resolve()
os.chdir(ROOT)

PASS  = "\033[92m PASS\033[0m"
FAIL  = "\033[91m FAIL\033[0m"
WARN  = "\033[93m WARN\033[0m"

failures = 0

def ok(msg):   print(f"{PASS}  {msg}")
def fail(msg):
    global failures
    failures += 1
    print(f"{FAIL}  {msg}")
def warn(msg): print(f"{WARN}  {msg}", flush=True)

print(f"\nGenRichi Setup Validator")
print(f"Project root: {ROOT}")
print("=" * 60)

# ── 1. Required files ─────────────────────────────────────────────────────────
print("\n[1] Required files")
required = [
    "workflow/Snakefile",
    "workflow/Snakefile_hereditary",
    "workflow/rules/qc.smk",
    "workflow/rules/align.smk",
    "workflow/rules/germline_calling.smk",
    "workflow/rules/annotation_germline.smk",
    "workflow/rules/acmg_classify.smk",
    "workflow/rules/report_hereditary.smk",
    "workflow/scripts/vcf_to_germline_table.py",
    "workflow/scripts/acmg_classifier.py",
    "workflow/scripts/generate_hereditary_report.py",
    "config/hereditary_config.yaml",
    "config/hereditary_samples.tsv",
    "resources/panel/hereditary_genes.bed",
    "envs/germline.yaml",
    "envs/annotation.yaml",
    "envs/report.yaml",
]
for f in required:
    p = ROOT / f
    if p.exists():
        ok(f)
    else:
        fail(f"MISSING: {f}")

# ── 2. Python syntax ──────────────────────────────────────────────────────────
print("\n[2] Python script syntax")
for f in (ROOT / "workflow/scripts").glob("*.py"):
    try:
        ast.parse(f.read_text(encoding="utf-8"))
        ok(f.name)
    except SyntaxError as e:
        fail(f"{f.name}: {e}")

# ── 3. YAML config parsing ────────────────────────────────────────────────────
print("\n[3] YAML config files")
try:
    import yaml
    for f in list((ROOT / "config").glob("*.yaml")) + list((ROOT / "envs").glob("*.yaml")):
        try:
            yaml.safe_load(f.read_text(encoding="utf-8"))
            ok(f.relative_to(ROOT))
        except yaml.YAMLError as e:
            fail(f"{f.name}: {e}")
except ImportError:
    warn("PyYAML not installed — skipping YAML validation (pip install pyyaml)")

# ── 4. Include / script references ───────────────────────────────────────────
print("\n[4] Include / script cross-references")
all_smk_text = ""
for smk in (ROOT / "workflow/rules").glob("*.smk"):
    all_smk_text += smk.read_text(encoding="utf-8")

# Includes in Snakefile_hereditary
sf_txt = (ROOT / "workflow/Snakefile_hereditary").read_text(encoding="utf-8")
for m in re.finditer(r'include:\s*["\'](.+?)["\']', sf_txt):
    inc = ROOT / "workflow" / m.group(1)
    if inc.exists():
        ok(f"include: {m.group(1)}")
    else:
        fail(f"include missing: {m.group(1)}")

# script: refs
for smk in (ROOT / "workflow/rules").glob("*.smk"):
    txt = smk.read_text(encoding="utf-8")
    for m in re.finditer(r'script:\s*["\'](.+?)["\']', txt):
        script_path = (smk.parent / m.group(1)).resolve()
        if script_path.exists():
            ok(f"script: {script_path.name} (in {smk.name})")
        else:
            fail(f"script missing: {m.group(1)} referenced in {smk.name}")

# ── 5. Snakemake on PATH ─────────────────────────────────────────────────────
print("\n[5] Tools on PATH")
for tool in ["snakemake", "conda", "mamba"]:
    path = shutil.which(tool)
    if path:
        # get version
        try:
            ver = subprocess.check_output(
                [tool, "--version"], stderr=subprocess.STDOUT, text=True
            ).strip()
            ok(f"{tool} ({ver})  →  {path}")
        except Exception:
            ok(f"{tool}  →  {path}")
    else:
        if tool == "snakemake":
            fail(f"{tool} not found — install with: conda install -c bioconda snakemake")
        else:
            warn(f"{tool} not found  (mamba/conda optional but recommended)")

# ── 6. Resource files (warnings only) ────────────────────────────────────────
print("\n[6] Resource files (warnings — not required until real run)")
resources = {
    "resources/reference/hg38.fa":                     "Reference genome (run setup_resources.sh)",
    "resources/reference/dbsnp_146.hg38.vcf.gz":      "dbSNP (run setup_resources.sh)",
    "resources/vep_cache":                             "VEP cache dir (run setup_resources.sh)",
    "resources/annotation/clinvar_20240101.vcf.gz":   "ClinVar VCF (run setup_resources.sh)",
}
for f, desc in resources.items():
    p = ROOT / f
    if p.exists():
        ok(f"{f}")
    else:
        warn(f"NOT FOUND: {f}  ({desc})")

# ── Summary ───────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
if failures == 0:
    print("\033[92mAll checks passed!\033[0m")
    print("\nNext step — dry-run the hereditary pipeline:")
    print(f"  cd {ROOT}")
    print(f"  snakemake --snakefile workflow/Snakefile_hereditary --configfile config/hereditary_config.yaml --use-conda --cores 8 --dry-run")
    print(f"\n  OR on Windows: double-click run_hereditary.bat")
else:
    print(f"\033[91m{failures} check(s) failed — fix issues above before running.\033[0m")
    sys.exit(1)
