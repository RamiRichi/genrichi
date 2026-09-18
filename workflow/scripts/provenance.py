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

# workflow/scripts/provenance.py -> workflow/scripts -> workflow -> repo root
REPO_ROOT = Path(__file__).resolve().parents[2]


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


# ── Assembly ──────────────────────────────────────────────────────────────
def build_provenance(sample_id, config, run_id=None, repo_root=None, panel_bed_path=None):
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

    return {
        "schema_version": "1.0",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sample_id": sample_id,
        "run_id": run_id,
        "pipeline": {
            "name": "GenRichi",
            "version": get_git_describe(repo_root),
            "git_commit_sha": get_git_commit(repo_root),
            "git_dirty": get_git_dirty(repo_root),
        },
        "panel": {
            "name": panel_cfg.get("name"),
            "bed_path": panel_bed_path,
            "bed_sha256": sha256_of_file(panel_bed_path),
        },
        "reference_genome": {
            "build": vep_cfg.get("genome_build"),
            "path": ref_cfg.get("genome"),
        },
        "reference_resources": {
            "dbsnp": _resource_record(ref_cfg.get("dbsnp")),
            "gnomad": _resource_record(ref_cfg.get("gnomad")),
            "panel_of_normals": _resource_record(ref_cfg.get("pon")),
        },
        "annotation_databases": {
            "vep_cache_version": _extract_vep_cache_version(vep_cfg.get("extra")),
            "vep_genome_build": vep_cfg.get("genome_build"),
            "cosmic": _resource_record(cosmic_cfg.get("vcf")),
            "clinvar": _resource_record(clinvar_cfg.get("vcf")),
        },
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
