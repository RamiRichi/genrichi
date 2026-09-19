"""
GenRichi Phase 5.1 — Reproducibility & Provenance.

Pure, dependency-light helpers that assemble a deterministic, machine-readable
provenance record for a single analysis run (pipeline/software version, git
commit, sample/run/panel identifiers, reference and annotation resource
versions, execution environment).

Every field degrades gracefully to `None` when the underlying information is
unavailable (e.g. no git repo, resource path not configured) — provenance
collection must never raise or block report generation.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)

import resource_inspect  # noqa: E402

# workflow/scripts/provenance.py -> workflow/scripts -> workflow -> repo root
REPO_ROOT = Path(__file__).resolve().parents[2]

SCHEMA_VERSION = "1.1"
# Tools whose observed version the provenance record is expected to carry.
EXPECTED_OBSERVED_TOOLS = ("vep", "bwa", "samtools", "bcftools", "multiqc")
BUILD_STAMP_RELPATH = os.path.join("workflow", "build_info.json")


# ── Git ───────────────────────────────────────────────────────────────────
def _run_git(args, repo_root):
    """Run a git command; return stripped stdout on success, None on any failure."""
    try:
        result = subprocess.run(
            ["git", *args],
            cwd=str(repo_root),
            capture_output=True,
            text=True,
            timeout=5,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if result.returncode != 0:
        return None
    return result.stdout.strip()


def get_git_commit(repo_root=None):
    """Full HEAD commit SHA, or None if git/repo is unavailable."""
    out = _run_git(["rev-parse", "HEAD"], repo_root or REPO_ROOT)
    return out or None


def get_git_describe(repo_root=None):
    """Human-readable pipeline version (tag-based if tags exist, else abbreviated commit); None if unavailable."""
    out = _run_git(["describe", "--tags", "--always", "--dirty"], repo_root or REPO_ROOT)
    return out or None


def get_git_dirty(repo_root=None):
    """True/False for uncommitted changes, None if git/repo is unavailable (distinct from a clean tree)."""
    out = _run_git(["status", "--porcelain"], repo_root or REPO_ROOT)
    if out is None:
        return None
    return bool(out)


# ── Pipeline identity: which code actually ran ────────────────────────────
def read_build_stamp(repo_root=None):
    """Contents of workflow/build_info.json (written at deploy time), or None if absent/unreadable."""
    path = os.path.join(str(repo_root or REPO_ROOT), BUILD_STAMP_RELPATH)
    try:
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) else None


def write_build_stamp(repo_root, out_path):
    """
    Record the git identity of a checkout into a small JSON file so a deployed
    copy of the pipeline that has no .git directory can still report the
    commit it was built from. Returns the stamp dict.
    """
    stamp = {
        "git_commit_sha": get_git_commit(repo_root),
        "git_describe": get_git_describe(repo_root),
        "git_dirty": get_git_dirty(repo_root),
        "stamped_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    with open(out_path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(stamp, fh, indent=2, sort_keys=True)
        fh.write("\n")
    return stamp


def compute_code_fingerprint(repo_root=None):
    """
    Content fingerprint of the workflow code that is deployed: SHA-256 over the
    sorted (relative path, file SHA-256) pairs of workflow/Snakefile*,
    workflow/rules/*.smk, workflow/scripts/*.py and workflow/envs/*.yaml.
    Independent of git, so it verifies "the code that actually ran" even for a
    deployed copy. None when there is no workflow/ directory.
    """
    root = Path(repo_root or REPO_ROOT) / "workflow"
    if not root.is_dir():
        return None
    entries = []
    for sub, patterns in (("", ("Snakefile*",)), ("rules", ("*.smk",)),
                          ("scripts", ("*.py",)), ("envs", ("*.yaml",))):
        d = root / sub if sub else root
        if not d.is_dir():
            continue
        for pattern in patterns:
            for f in d.glob(pattern):
                if f.is_file():
                    entries.append((f.relative_to(root).as_posix(), sha256_of_file(str(f))))
    if not entries:
        return None
    h = hashlib.sha256()
    for rel, digest in sorted(entries):
        h.update(f"{rel}\t{digest}\n".encode("utf-8"))
    return h.hexdigest()


def resolve_pipeline_identity(repo_root=None):
    """
    Commit SHA of the code actually used: live git if available, otherwise the
    deploy-time build stamp, otherwise explicitly unavailable — never a guess.
    """
    repo_root = repo_root or REPO_ROOT
    commit = get_git_commit(repo_root)
    if commit:
        return {
            "version": get_git_describe(repo_root),
            "git_commit_sha": commit,
            "git_dirty": get_git_dirty(repo_root),
            "git_commit_source": "git",
            "build_stamp": None,
        }
    stamp = read_build_stamp(repo_root)
    if stamp and stamp.get("git_commit_sha"):
        return {
            "version": stamp.get("git_describe"),
            "git_commit_sha": stamp.get("git_commit_sha"),
            "git_dirty": stamp.get("git_dirty"),
            "git_commit_source": "build_stamp",
            "build_stamp": stamp,
        }
    return {"version": None, "git_commit_sha": None, "git_dirty": None,
            "git_commit_source": "unavailable", "build_stamp": None}


# ── Resource version / hashing ───────────────────────────────────────────
_VERSION_PATTERNS = (
    (re.compile(r"[_-]v(\d+(?:\.\d+)*)"), lambda m: f"v{m.group(1)}"),
    (re.compile(r"_(\d{8})(?=[._]|$)"), lambda m: m.group(1)),
    (re.compile(r"_(\d+)\."), lambda m: m.group(1)),
)


def extract_version_token(path):
    """
    Best-effort, non-fabricating extraction of a version/build token from a
    reference-resource filename, e.g. 'dbsnp_146.hg38.vcf.gz' -> '146',
    'CosmicCodingMuts_v99_GRCh38.vcf.gz' -> 'v99',
    'clinvar_20240101.vcf.gz' -> '20240101'.

    Falls back to the bare filename when no numeric/version token is present
    in it, so this never invents a version that isn't literally in the
    configured file reference. Returns None only when no path was given.
    """
    if not path:
        return None
    name = os.path.basename(str(path))
    for pattern, extract in _VERSION_PATTERNS:
        m = pattern.search(name)
        if m:
            return extract(m)
    return name


def sha256_of_file(path, chunk_size=1 << 20):
    """SHA-256 of a (small) file's contents, or None if missing/unreadable."""
    if not path:
        return None
    try:
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(chunk_size), b""):
                h.update(chunk)
        return h.hexdigest()
    except OSError:
        return None


def sha256_of_config(config):
    """Deterministic (key-order-independent) hash of the resolved run configuration."""
    try:
        blob = json.dumps(config, sort_keys=True, default=str).encode("utf-8")
    except (TypeError, ValueError):
        return None
    return hashlib.sha256(blob).hexdigest()


def _resource_record(path):
    """path + best-effort version token for a reference/annotation resource.

    Deliberately does not hash resource contents — these files (dbSNP,
    gnomAD, COSMIC, ClinVar, genome FASTA) can be multi-gigabyte, and hashing
    them on every report would add unacceptable, unrequested runtime cost.
    """
    return {
        "path": str(path) if path else None,
        "version": extract_version_token(path),
    }


def _extract_vep_cache_version(vep_extra):
    if not vep_extra:
        return None
    m = re.search(r"--cache_version[= ](\S+)", str(vep_extra))
    return m.group(1) if m else None


# ── Execution environment ────────────────────────────────────────────────
def get_execution_environment():
    return {
        "python_version": sys.version.split()[0],
        "platform": platform.platform(),
        "hostname": platform.node() or None,
        "conda_env": os.environ.get("CONDA_DEFAULT_ENV"),
    }


# ── Observed (not merely configured) facts ────────────────────────────────
def _safe(fn, *args, **kwargs):
    """Provenance collection must never raise or block report generation."""
    try:
        return fn(*args, **kwargs)
    except Exception as exc:  # noqa: BLE001 - deliberate: observation is best-effort
        return {"error": f"{type(exc).__name__}: {exc}"}


def _safe_value(fn, *args):
    try:
        return fn(*args)
    except Exception:  # noqa: BLE001
        return None


def _resource_record_observed(path):
    """Configured record (path + filename token) plus what the file itself says."""
    record = _resource_record(path)
    if not path:
        record["observed"] = None
        return record
    observed = {"file": resource_inspect.file_facts(path)}
    if observed["file"]["exists"]:
        observed["vcf_header"] = resource_inspect.read_vcf_header(path)
    record["observed"] = observed
    return record


def _observed_reference(config, ref_cfg, vep_cfg):
    genome = ref_cfg.get("genome")
    if not isinstance(genome, str) or not os.path.isfile(genome):
        return None
    errors, warnings, facts = resource_inspect.check_fasta_indexes(genome)
    out = {
        "fasta_file": resource_inspect.file_facts(genome),
        "contig_count": facts.get("contig_count"),
        "total_bases": facts.get("total_bases"),
        "first_contig": facts.get("first_contig"),
        "fai_sha256": facts.get("fai_sha256"),
        "dict_sha256": facts.get("dict_sha256"),
        "fasta_fai_dict_consistent": (not errors) if facts else None,
        "consistency_errors": errors,
        "consistency_warnings": warnings,
    }
    if facts and not errors:
        b_err, b_warn, b_facts = resource_inspect.check_bwa_index(
            genome, expected_total_bases=facts.get("total_bases"), expected_contigs=facts.get("contig_count"))
        out["bwa_index"] = {"consistent_with_fasta": not b_err, "errors": b_err, "warnings": b_warn, **b_facts}
    manifest_path = (ref_cfg.get("manifest") if isinstance(ref_cfg.get("manifest"), str) else None)         or os.path.join("config", "reference_manifest.json")
    if os.path.isfile(manifest_path):
        matches = None
        try:
            with open(manifest_path, encoding="utf-8") as fh:
                matches = not resource_inspect.compare_to_manifest(json.load(fh), facts, vep_cfg.get("genome_build"))
        except (OSError, ValueError):
            pass
        out["manifest"] = {"path": manifest_path, "sha256": sha256_of_file(manifest_path), "matches": matches}
    else:
        out["manifest"] = {"path": manifest_path, "sha256": None, "matches": None}
    return out


def _observed_vep_cache(vep_cfg):
    cache_dir = vep_cfg.get("cache_dir")
    if not isinstance(cache_dir, str) or not os.path.isdir(cache_dir):
        return None
    errors, warnings, facts = resource_inspect.inspect_vep_cache(
        cache_dir, _extract_vep_cache_version(vep_cfg.get("extra")), vep_cfg.get("genome_build"))
    return {
        "cache_path": facts.get("cache_path"),
        "version_dir": facts.get("version_dir"),
        "cache_version": facts.get("cache_version"),
        "assembly": facts.get("assembly"),
        "info_txt_sha256": facts.get("info_txt_sha256"),
        "info": facts.get("info"),
        "consistency_errors": errors,
        "consistency_warnings": warnings,
    }


def load_observed_software(paths):
    """
    Merge the runtime-probe JSON files written inside each conda environment
    into one tool -> observed version map. Unreadable files are reported, not fatal.
    """
    tools: dict = {}
    problems: list = []
    read_any = False
    for p in paths or []:
        try:
            with open(p, encoding="utf-8") as fh:
                rec = json.load(fh)
        except (OSError, ValueError):
            problems.append(f"unreadable tool-version file: {p}")
            continue
        read_any = True
        for name, info in (rec.get("tools") or {}).items():
            entry = {
                "version": info.get("version"),
                "path": info.get("path"),
                "observed_in_environment": rec.get("label"),
                "conda_prefix": rec.get("conda_prefix"),
                "observed_at_utc": rec.get("observed_at_utc"),
            }
            if info.get("htslib"):
                entry["htslib_version"] = info["htslib"]
            if name == "vep":
                entry["component_versions"] = info.get("detail")
                entry["cache_compatibility_check"] = rec.get("vep_cache_check")
            tools[name] = entry
    if not read_any:
        status = "not_collected"
    else:
        missing = [t for t in EXPECTED_OBSERVED_TOOLS if not (tools.get(t) or {}).get("version")]
        status = "complete" if not missing else "partial: missing " + ", ".join(missing)
    return {"status": status, "tools": tools, "problems": problems}


def _clinical_annotation_sources(clinvar_record, cache_observed):
    """
    Which ClinVar feeds the report, and which ClinVar is merely an annotation
    dependency. The report's ClinVar column comes from the custom VEP
    annotation (ClinVar_CLNSIG / ClinVar_CLNDN), NOT the cache's CLIN_SIG.
    """
    obs = (clinvar_record or {}).get("observed") or {}
    header = obs.get("vcf_header") or {}
    file_facts_ = obs.get("file") or {}
    file_date = header.get("file_date")
    cache_info = (cache_observed or {}).get("info") or {}
    return {
        "report_clinvar": {
            "role": "report_clinical_classification_source",
            "used_for": "ClinVar_CLNSIG / ClinVar_CLNDN fields added by VEP --custom (exact match); "
                        "they populate the report's ClinVar column",
            "path": (clinvar_record or {}).get("path"),
            "version": file_date,
            "version_basis": "vcf_header_fileDate" if file_date else "unavailable (no ##fileDate in the VCF header)",
            "vcf_header": {k: header.get(k) for k in ("fileformat", "file_date", "source", "reference")},
            "size_bytes": file_facts_.get("size_bytes"),
            "sha256": file_facts_.get("sha256"),
        },
        "vep_cache_clinvar": {
            "role": "annotation_dependency_not_report_classification_source",
            "used_for": "VEP cache CLIN_SIG from --everything/--check_existing; not read by the report",
            "version": cache_info.get("source_ClinVar"),
            "version_basis": "vep_cache_info.txt:source_ClinVar" if cache_info.get("source_ClinVar") else "unavailable",
            "cache_path": (cache_observed or {}).get("cache_path"),
        },
    }


def _key_file_checksums(config, ref_cfg, clinvar_path, panel_bed_path):
    """SHA-256 (with size/mtime) of the important small/medium files a result depends on."""
    genome = ref_cfg.get("genome")
    tmb_bed = (config.get("tmb") or {}).get("coding_bed")
    candidates = {
        "panel_bed": panel_bed_path,
        "tmb_coding_bed": tmb_bed,
        "reference_fai": genome + ".fai" if isinstance(genome, str) else None,
        "reference_dict": resource_inspect.dict_path_for(genome) if isinstance(genome, str) else None,
        "clinvar_vcf": clinvar_path,
        "clinvar_vcf_index": clinvar_path + ".tbi" if isinstance(clinvar_path, str) else None,
    }
    return {name: resource_inspect.file_facts(path) for name, path in candidates.items() if path}


# ── Assembly ──────────────────────────────────────────────────────────────
def build_provenance(sample_id, config, run_id=None, repo_root=None, panel_bed_path=None,
                     tool_version_files=None):
    """
    Assemble a deterministic, machine-readable provenance record for one
    sample's analysis run.

    `config` is the fully-resolved Snakemake run configuration (a plain
    dict). Optional/unavailable fields (git info outside a repo, an
    unconfigured PoN, a missing run_id, ...) degrade to `None` rather than
    raising or altering pipeline behavior.
    """
    repo_root = repo_root or REPO_ROOT
    config = config if isinstance(config, dict) else {}
    ref_cfg = config.get("ref") or {}
    ann_cfg = config.get("annotation") or {}
    vep_cfg = ann_cfg.get("vep") or {}
    cosmic_cfg = ann_cfg.get("cosmic") or {}
    clinvar_cfg = ann_cfg.get("clinvar") or {}
    panel_cfg = config.get("panel") or {}

    panel_bed_path = panel_bed_path or panel_cfg.get("bed")

    identity = _safe(resolve_pipeline_identity, repo_root)
    if "error" in identity:
        identity = {"version": None, "git_commit_sha": None, "git_dirty": None,
                    "git_commit_source": "unavailable", "build_stamp": None}
    clinvar_record = _safe(_resource_record_observed, clinvar_cfg.get("vcf"))
    cache_observed = _safe(_observed_vep_cache, vep_cfg)
    software = _safe(load_observed_software, tool_version_files)

    checksums = _safe(_key_file_checksums, config, ref_cfg, clinvar_cfg.get("vcf"), panel_bed_path)
    if isinstance(checksums, dict) and isinstance(cache_observed, dict) and cache_observed.get("cache_path"):
        checksums["vep_cache_info_txt"] = _safe(
            resource_inspect.file_facts, os.path.join(cache_observed["cache_path"], "info.txt"))

    return {
        "schema_version": SCHEMA_VERSION,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sample_id": sample_id,
        "run_id": run_id,
        "pipeline": {
            "name": "GenRichi",
            "version": identity["version"],
            "git_commit_sha": identity["git_commit_sha"],
            "git_dirty": identity["git_dirty"],
            "git_commit_source": identity["git_commit_source"],
            "build_stamp": identity["build_stamp"],
            "code_fingerprint_sha256": _safe_value(compute_code_fingerprint, repo_root),
        },
        "panel": {
            "name": panel_cfg.get("name"),
            "bed_path": panel_bed_path,
            "bed_sha256": sha256_of_file(panel_bed_path),
        },
        # Configured values stay under their original keys; what the files/tools
        # themselves report is under "observed" / "software_observed".
        "reference_genome": {
            "build": vep_cfg.get("genome_build"),
            "path": ref_cfg.get("genome"),
            "observed": _safe(_observed_reference, config, ref_cfg, vep_cfg),
        },
        "reference_resources": {
            "dbsnp": _safe(_resource_record_observed, ref_cfg.get("dbsnp")),
            "gnomad": _safe(_resource_record_observed, ref_cfg.get("gnomad")),
            "panel_of_normals": _resource_record(ref_cfg.get("pon")),
        },
        "annotation_databases": {
            "vep_cache_version": _extract_vep_cache_version(vep_cfg.get("extra")),
            "vep_genome_build": vep_cfg.get("genome_build"),
            "vep_cache_observed": cache_observed,
            "cosmic": _resource_record(cosmic_cfg.get("vcf")),
            "clinvar": clinvar_record,
        },
        "clinical_annotation_sources": _safe(_clinical_annotation_sources, clinvar_record, cache_observed),
        "software_observed": software,
        "checksums": checksums,
        "execution_environment": get_execution_environment(),
        "config_sha256": sha256_of_config(config),
    }


# ── HTML embedding helpers ───────────────────────────────────────────────
def provenance_script_tag(provenance, element_id="genrichi-provenance"):
    """Deterministic, machine-readable <script type="application/json"> block for embedding in a report."""
    payload = json.dumps(provenance, indent=2, sort_keys=True, default=str)
    return f'<script type="application/json" id="{element_id}">\n{payload}\n</script>'


def provenance_footer_line(provenance):
    """One-line, human-readable provenance summary for a report footer."""
    pipeline = provenance.get("pipeline", {}) or {}
    version = pipeline.get("version") or "unknown"
    commit = pipeline.get("git_commit_sha")
    commit_short = commit[:12] if commit else "unknown"
    generated = provenance.get("generated_at_utc", "unknown")
    return f"Pipeline {version} (commit {commit_short}) &bull; Provenance recorded {generated} UTC"


def _main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="GenRichi provenance helpers")
    ap.add_argument("--write-build-stamp", metavar="OUT_JSON",
                    help="record this checkout's git commit/describe/dirty state (for deployed copies without .git)")
    ap.add_argument("--repo", default=None, help="repository root (default: this checkout)")
    args = ap.parse_args(argv)
    if not args.write_build_stamp:
        ap.print_help()
        return 2
    stamp = write_build_stamp(args.repo or REPO_ROOT, args.write_build_stamp)
    print(f"wrote {args.write_build_stamp}: commit {stamp['git_commit_sha']} dirty={stamp['git_dirty']}")
    return 0 if stamp["git_commit_sha"] else 1


if __name__ == "__main__":
    sys.exit(_main())
