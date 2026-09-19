"""
GenRichi Phase 5.3 — observe the tool versions actually in use.

Runs inside a rule's conda environment (so the executables it finds are the
ones that rule's siblings run), records the versions the tools themselves
report, and — for the annotation environment — verifies the VEP executable
matches the VEP cache it will be pointed at. Writes a small JSON file that the
report/provenance step reads.

CLI:
    python runtime_probe.py --label annotation --tools vep,bcftools \
        --out results/runtime/annotation_tool_versions.json \
        --vep-cache-dir DIR --vep-cache-version 113 --vep-assembly GRCh38

Exit status is non-zero (and the JSON still written) when a requested tool is
missing or the VEP executable does not match the VEP cache — a hard failure
before annotation, not a silent mismatch.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)

import resource_inspect  # noqa: E402
from resource_inspect import vep_cache_version_from_extra  # noqa: E402,F401  (re-exported for callers/tests)


# ── Parsers (pure; unit-tested against real captured output) ─────────────
def _first(regex, text, flags=re.MULTILINE):
    m = re.search(regex, text or "", flags)
    return m.group(1) if m else None


def parse_vep(text):
    versions = {k: v for k, v in re.findall(r"^\s+(ensembl[\w-]*)\s*:\s*(\S+)\s*$", text or "", re.MULTILINE)}
    return {"version": versions.get("ensembl-vep"), "detail": versions}


def parse_bwa(text):
    return {"version": _first(r"^Version:\s*(\S+)", text)}


def parse_htslib_tool(name):
    def _parse(text):
        return {"version": _first(rf"^{name}\s+(\S+)", text), "htslib": _first(r"^Using htslib\s+(\S+)", text)}
    return _parse


def parse_multiqc(text):
    return {"version": _first(r"version\s+(\S+)", text)}


def parse_fastp(text):
    return {"version": _first(r"^fastp\s+(\S+)", text)}


def parse_mosdepth(text):
    return {"version": _first(r"^mosdepth\s+(\S+)", text)}


def parse_gatk(text):
    return {"version": _first(r"GATK\)\s+v?(\S+)", text)}


TOOLS = {
    "vep": (["vep", "--help"], parse_vep),
    "bwa": (["bwa"], parse_bwa),
    "samtools": (["samtools", "--version"], parse_htslib_tool("samtools")),
    "bcftools": (["bcftools", "--version"], parse_htslib_tool("bcftools")),
    "multiqc": (["multiqc", "--version"], parse_multiqc),
    "fastp": (["fastp", "--version"], parse_fastp),
    "mosdepth": (["mosdepth", "--version"], parse_mosdepth),
    "gatk": (["gatk", "--version"], parse_gatk),
}


# ── Running / probing ────────────────────────────────────────────────────
def run_command(cmd, timeout: int = 120):
    """(returncode, combined stdout+stderr). Never raises; (None, message) if it cannot run."""
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except (OSError, subprocess.SubprocessError) as exc:
        return None, f"{type(exc).__name__}: {exc}"
    return proc.returncode, (proc.stdout or "") + (proc.stderr or "")


def probe_tool(name, runner=run_command, which=shutil.which):
    """Observed facts for one tool: found?, resolved path, version reported by the tool itself."""
    if name not in TOOLS:
        raise KeyError(f"unknown tool {name!r}")
    cmd, parser = TOOLS[name]
    path = which(cmd[0])
    record = {"found": bool(path), "path": path, "version": None, "command": " ".join(cmd)}
    if not path:
        return record
    rc, text = runner(cmd)
    parsed = parser(text)
    record.update({k: v for k, v in parsed.items() if v is not None})
    record.setdefault("version", None)
    record["first_output_line"] = next((ln.strip() for ln in (text or "").splitlines() if ln.strip()), None)
    return record


def check_vep_matches_cache(vep_version, cache_dir, cache_version, assembly):
    """
    The VEP executable's release must equal the VEP cache release it will read.
    Returns {'status': 'PASS'|'FAIL'|'SKIPPED', 'message': ..., ...}.
    """
    result = {"status": "SKIPPED", "vep_version": vep_version, "vep_release": None,
              "cache_version_pinned": cache_version or None, "cache_version_used": None,
              "assembly": assembly, "cache_dir": cache_dir, "message": ""}
    m = re.match(r"(\d+)", str(vep_version or ""))
    if not m:
        result.update(status="FAIL", message=f"could not read the VEP executable version (got {vep_version!r})")
        return result
    release = int(m.group(1))
    result["vep_release"] = release
    if not cache_dir:
        result["message"] = "no VEP cache directory supplied; compatibility not checked"
        return result
    expected_asm = resource_inspect.normalize_build(assembly)
    used = str(cache_version) if cache_version else str(release)
    result["cache_version_used"] = used
    cache_path = os.path.join(cache_dir, "homo_sapiens", f"{used}_{expected_asm}") if expected_asm else None
    if str(release) != used:
        result.update(
            status="FAIL",
            message=(f"VEP executable is release {release} ({vep_version}) but the cache version in use is {used}; "
                     "VEP requires the executable and cache to be the same release"),
        )
    elif cache_path and not os.path.isdir(cache_path):
        result.update(status="FAIL", message=f"VEP cache directory not found: {cache_path}")
    else:
        result.update(status="PASS", message=f"VEP {vep_version} matches cache {used}_{expected_asm}")
    return result


def build_record(label, tool_names, *, vep_cache=None, runner=run_command, which=shutil.which):
    """Assemble the JSON record for one environment."""
    record = {
        "schema_version": "1.0",
        "label": label,
        "observed_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "conda_prefix": os.environ.get("CONDA_PREFIX"),
        "python_version": sys.version.split()[0],
        "tools": {},
    }
    for name in tool_names:
        record["tools"][name] = probe_tool(name, runner=runner, which=which)
    if "vep" in record["tools"] and vep_cache is not None:
        record["vep_cache_check"] = check_vep_matches_cache(
            record["tools"]["vep"].get("version"),
            vep_cache.get("cache_dir"), vep_cache.get("cache_version"), vep_cache.get("assembly"),
        )
    return record


def record_problems(record) -> list:
    problems = []
    for name, info in (record.get("tools") or {}).items():
        if not info.get("found"):
            problems.append(f"required tool '{name}' not found in the {record.get('label')} environment")
        elif not info.get("version"):
            problems.append(f"could not determine the version of '{name}' (output not recognised)")
    check = record.get("vep_cache_check")
    if check and check.get("status") == "FAIL":
        problems.append(check.get("message"))
    return problems


def summary_line(record) -> str:
    parts = []
    for name, info in (record.get("tools") or {}).items():
        v = info.get("version") or "UNKNOWN"
        extra = f" (htslib {info['htslib']})" if info.get("htslib") else ""
        parts.append(f"{name} {v}{extra}")
    check = record.get("vep_cache_check")
    if check:
        parts.append(f"VEP-cache check: {check['status']}")
    return f"[runtime-probe:{record.get('label')}] " + "; ".join(parts)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--label", required=True)
    ap.add_argument("--tools", required=True, help="comma-separated: " + ",".join(TOOLS))
    ap.add_argument("--out", required=True)
    ap.add_argument("--vep-cache-dir")
    ap.add_argument("--vep-cache-version", default="")
    ap.add_argument("--vep-assembly")
    args = ap.parse_args(argv)

    names = [t.strip() for t in args.tools.split(",") if t.strip()]
    unknown = [t for t in names if t not in TOOLS]
    if unknown:
        print(f"unknown tool(s): {', '.join(unknown)}", file=sys.stderr)
        return 2
    vep_cache = None
    if "vep" in names and args.vep_cache_dir:
        vep_cache = {"cache_dir": args.vep_cache_dir, "cache_version": args.vep_cache_version or None,
                     "assembly": args.vep_assembly}
    record = build_record(args.label, names, vep_cache=vep_cache)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as fh:
        json.dump(record, fh, indent=2, sort_keys=True)
    print(summary_line(record))
    problems = record_problems(record)
    for p in problems:
        print(f"[runtime-probe:{args.label}] ERROR: {p}", file=sys.stderr)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
