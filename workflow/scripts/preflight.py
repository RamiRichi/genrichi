"""
GenRichi Phase 5.2 — Pre-flight validation for the comprehensive (Phase 1)
solid-tumor pipeline.

Pure, dependency-light (stdlib + pandas, which the Snakefile already
requires) validation of the resolved run configuration and samples sheet,
run once at Snakefile parse time, before any compute-heavy rule (alignment,
variant calling, VEP annotation) is scheduled.

Design goals:
  * One clean validation gate, not duplicated per-rule.
  * Every failure names the exact field/sample/file at fault.
  * Never a raw pandas/KeyError/Snakemake traceback as the primary
    diagnostic.
  * Does not change, threshold, or reinterpret any scientific value — it
    only checks that the inputs required to run the existing logic are
    present, well-formed, and self-consistent.
"""

from __future__ import annotations

import json
import os
import sys
from dataclasses import dataclass, field
from typing import Any, Iterable

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)

import resource_inspect  # noqa: E402

try:
    import pandas as pd
except ImportError:  # pragma: no cover - pandas is already a hard Snakefile dependency
    pd = None


# ── Data model ────────────────────────────────────────────────────────────
@dataclass
class Issue:
    field: str
    message: str


@dataclass
class PreflightReport:
    errors: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    info: list = field(default_factory=list)  # observed facts (Phase 5.3); never affects `ok`

    @property
    def ok(self) -> bool:
        return not self.errors


_MISSING = object()

REQUIRED_SAMPLE_COLUMNS = ("tumor_r1", "tumor_r2", "normal_r1", "normal_r2")

_REQUIRED_CONFIG_KEYS = (
    "samples",
    "ref.genome", "ref.dbsnp", "ref.gnomad",
    "panel.bed", "panel.name",
    "calling.filter.min_af", "calling.filter.min_depth", "calling.filter.min_alt_reads",
    "annotation.vep.cache_dir", "annotation.vep.genome_build",
    "annotation.clinvar.vcf",
    "cnv.amp_threshold", "cnv.del_threshold", "cnv.min_probes",
    "msi.threshold",
    "tmb.coding_mb", "tmb.coding_bed", "tmb.partial_threshold", "tmb.high_threshold",
    "report.company",
)
# Intentionally NOT required (existing, deliberate graceful degradation,
# left untouched by this phase): ref.pon, annotation.cosmic.vcf, report.logo.

_RANGE_CHECKS = {
    "calling.filter.min_af": lambda v: 0 <= float(v) <= 1,
    "calling.filter.min_depth": lambda v: float(v) > 0,
    "calling.filter.min_alt_reads": lambda v: float(v) > 0,
    "cnv.amp_threshold": lambda v: float(v) > 0,
    "cnv.del_threshold": lambda v: float(v) < 0,
    "cnv.min_probes": lambda v: float(v) > 0,
    "msi.threshold": lambda v: float(v) > 0,
    "tmb.coding_mb": lambda v: float(v) > 0,
    "tmb.partial_threshold": lambda v: 0 < float(v) <= 1,
    "tmb.high_threshold": lambda v: float(v) > 0,
    "annotation.vep.genome_build": lambda v: str(v) in {"GRCh37", "GRCh38", "hg19", "hg38"},
}


def _get_dotted(config: dict, dotted_path: str) -> Any:
    node = config
    for part in dotted_path.split("."):
        if not isinstance(node, dict) or part not in node:
            return _MISSING
        node = node[part]
    return node


def _is_blank(value) -> bool:
    if value is None:
        return True
    if pd is not None:
        try:
            if pd.isna(value):
                return True
        except (TypeError, ValueError):
            pass
    return not str(value).strip()


# ── Config validation ────────────────────────────────────────────────────
def validate_config_keys(config: dict) -> list:
    if not isinstance(config, dict):
        return [Issue("config", "Configuration is missing or is not a mapping.")]
    issues = []
    for key in _REQUIRED_CONFIG_KEYS:
        if _get_dotted(config, key) is _MISSING:
            issues.append(Issue(key, f"Required config key is missing: {key}"))
    return issues


def validate_config_ranges(config: dict) -> list:
    if not isinstance(config, dict):
        return []
    issues = []
    for key, validator in _RANGE_CHECKS.items():
        val = _get_dotted(config, key)
        if val is _MISSING or val is None:
            continue  # absence already reported by validate_config_keys
        try:
            ok = bool(validator(val))
        except (TypeError, ValueError):
            ok = False
        if not ok:
            issues.append(Issue(key, f"Config value out of expected range/type for {key}: {val!r}"))
    return issues


# ── Samples-sheet validation ─────────────────────────────────────────────
def validate_samples_sheet(samples_df) -> list:
    issues: list = []
    if samples_df is None or pd is None:
        return issues

    # Blank / duplicate sample_id
    dup_mask = samples_df.index.duplicated(keep=False)
    if dup_mask.any():
        dups = sorted({str(s) for s in samples_df.index[dup_mask]})
        issues.append(Issue("samples.sample_id", f"Duplicate sample_id value(s): {', '.join(dups)}"))
    for sample_id in samples_df.index:
        if _is_blank(sample_id):
            issues.append(Issue("samples.sample_id", "Blank/NaN sample_id in samples sheet."))

    missing_cols = [c for c in REQUIRED_SAMPLE_COLUMNS if c not in samples_df.columns]
    if missing_cols:
        issues.append(Issue(
            "samples.columns",
            f"Samples sheet is missing required column(s): {', '.join(missing_cols)}",
        ))
        return issues  # per-row checks below would be meaningless without these columns

    for sample_id, row in samples_df.iterrows():
        for col in REQUIRED_SAMPLE_COLUMNS:
            if _is_blank(row.get(col)):
                issues.append(Issue(
                    f"samples.{sample_id}.{col}",
                    f"Sample '{sample_id}': required field '{col}' is blank/missing.",
                ))
    return issues


def validate_fastq_files(samples_df, base_dir: str | None = None) -> list:
    issues: list = []
    if samples_df is None or pd is None:
        return issues
    if any(c not in samples_df.columns for c in REQUIRED_SAMPLE_COLUMNS):
        return issues  # already reported by validate_samples_sheet

    for sample_id, row in samples_df.iterrows():
        for col in REQUIRED_SAMPLE_COLUMNS:
            val = row.get(col)
            if _is_blank(val):
                continue  # already reported by validate_samples_sheet; avoid duplicate noise
            path = str(val)
            resolved = path if os.path.isabs(path) or not base_dir else os.path.join(base_dir, path)
            if not os.path.isfile(resolved):
                issues.append(Issue(
                    f"samples.{sample_id}.{col}",
                    f"Sample '{sample_id}': {col} file not found: {path}",
                ))
            elif os.path.getsize(resolved) == 0:
                issues.append(Issue(
                    f"samples.{sample_id}.{col}",
                    f"Sample '{sample_id}': {col} file is empty (0 bytes): {path}",
                ))
    return issues


# ── Resource / reference validation ──────────────────────────────────────
def _require_file(issues: list, path, field_name: str) -> bool:
    if not path:
        issues.append(Issue(field_name, f"Required resource not configured: {field_name}"))
        return False
    if not os.path.isfile(path):
        issues.append(Issue(field_name, f"Required resource file not found: {path}"))
        return False
    return True


def _require_index(issues: list, path: str, field_name: str, suffixes=(".tbi", ".idx")) -> None:
    if any(os.path.isfile(path + suffix) for suffix in suffixes):
        return
    issues.append(Issue(
        field_name,
        f"Required index ({'/'.join(suffixes)}) not found alongside: {path}",
    ))


def validate_resources(config: dict) -> list:
    issues: list = []
    if not isinstance(config, dict):
        return issues

    ref = config.get("ref") or {}
    ann = config.get("annotation") or {}
    vep = ann.get("vep") or {}
    panel = config.get("panel") or {}

    genome = ref.get("genome")
    if _require_file(issues, genome, "ref.genome"):
        _require_index(issues, genome, "ref.genome (.fai index)", suffixes=(".fai",))
        dict_path = os.path.splitext(genome)[0] + ".dict"
        if not os.path.isfile(dict_path):
            issues.append(Issue(
                "ref.genome (.dict)",
                f"Required GATK sequence dictionary not found: {dict_path}",
            ))

    dbsnp = ref.get("dbsnp")
    if _require_file(issues, dbsnp, "ref.dbsnp"):
        _require_index(issues, dbsnp, "ref.dbsnp (index)")

    gnomad = ref.get("gnomad")
    if _require_file(issues, gnomad, "ref.gnomad"):
        _require_index(issues, gnomad, "ref.gnomad (index)")

    pon = ref.get("pon")
    if pon:  # explicitly optional — only validated if the user configured one
        _require_file(issues, pon, "ref.pon")

    _require_file(issues, panel.get("bed"), "panel.bed")

    cache_dir = vep.get("cache_dir")
    if not cache_dir:
        issues.append(Issue("annotation.vep.cache_dir", "Required resource not configured: annotation.vep.cache_dir"))
    elif not os.path.isdir(cache_dir):
        issues.append(Issue("annotation.vep.cache_dir", f"Required VEP cache directory not found: {cache_dir}"))

    clinvar = (ann.get("clinvar") or {}).get("vcf")
    if _require_file(issues, clinvar, "annotation.clinvar.vcf"):
        _require_index(issues, clinvar, "annotation.clinvar.vcf (index)")

    # NOTE: tmb.coding_bed is deliberately NOT hard-required here. It is
    # already consumed defensively by generate_comprehensive_report.py
    # (_read_cds_bed) which degrades to an explicit "TMB: N/A" rather than a
    # misleading result when absent — unlike the MSI/CNV silent-negative
    # gaps this phase targets, that existing path is not a false-negative
    # risk, so hard-failing it here would be a behavior change, not a fix.

    return issues


# ── Runtime-resource validation (Phase 5.3): observe the files, don't trust the config ──
def _str_path(value):
    return value if isinstance(value, str) and value else None


def _resolve_manifest_path(config: dict, base_dir):
    """(path, explicitly_configured). Default: <base_dir>/config/reference_manifest.json."""
    explicit = _str_path((config.get("ref") or {}).get("manifest"))
    if explicit:
        return (explicit if os.path.isabs(explicit) or not base_dir else os.path.join(base_dir, explicit)), True
    return os.path.join(base_dir or os.getcwd(), "config", "reference_manifest.json"), False


def validate_runtime_resources(config: dict, base_dir=None):
    """
    Deep checks on the resources the pipeline will actually read.

    Returns (errors, warnings, info) as lists of Issue. Files that are simply
    missing are NOT reported here (validate_resources already does), so a
    missing resource never produces two messages.

    Errors are limited to conditions that make the run invalid or unsafe
    (inconsistent FASTA/.fai/.dict, BWA index built from another FASTA,
    unusable VCF index, VEP cache/assembly mismatch, reference differing from
    the configured manifest). Anything softer is a warning.
    """
    errors: list = []
    warnings: list = []
    info: list = []
    if not isinstance(config, dict):
        return errors, warnings, info

    ref = config.get("ref") or {}
    ann = config.get("annotation") or {}
    vep = ann.get("vep") or {}
    build = vep.get("genome_build")

    def add(target, field_name, messages):
        for m in messages:
            target.append(Issue(field_name, m))

    # ── Reference FASTA, .fai, .dict, BWA index, manifest ──
    genome = _str_path(ref.get("genome"))
    if genome and os.path.isfile(genome):
        e, w, facts = resource_inspect.check_fasta_indexes(genome)
        add(errors, "ref.genome (FASTA/.fai/.dict consistency)", e)
        add(warnings, "ref.genome (FASTA/.fai/.dict consistency)", w)
        if facts and not e:
            e2, w2, bfacts = resource_inspect.check_bwa_index(
                genome, expected_total_bases=facts.get("total_bases"), expected_contigs=facts.get("contig_count"))
            add(errors, "ref.genome (BWA index)", e2)
            add(warnings, "ref.genome (BWA index)", w2)

            manifest_path, explicit = _resolve_manifest_path(config, base_dir)
            manifest_note = ""
            if os.path.isfile(manifest_path):
                try:
                    with open(manifest_path, encoding="utf-8") as fh:
                        manifest = json.load(fh)
                except (OSError, ValueError) as exc:
                    manifest = None
                    errors.append(Issue("ref.manifest", f"Reference manifest unreadable: {manifest_path}: {exc}"))
                if manifest is not None:
                    diffs = resource_inspect.compare_to_manifest(manifest, facts, build)
                    if diffs:
                        errors.append(Issue(
                            "ref.genome (expected resource set)",
                            "The configured FASTA is not the reference recorded in "
                            f"{manifest_path}: " + "; ".join(diffs) + ". If the reference was changed "
                            "deliberately it must be re-validated and the manifest regenerated.",
                        ))
                    else:
                        manifest_note = f"; matches manifest {os.path.basename(manifest_path)}"
            elif explicit:
                errors.append(Issue("ref.manifest", f"Configured reference manifest not found: {manifest_path}"))
            else:
                warnings.append(Issue(
                    "ref.genome (expected resource set)",
                    f"No reference manifest found ({manifest_path}); cannot confirm the FASTA is the expected "
                    "resource set. Create one with: python workflow/scripts/preflight.py "
                    "--write-reference-manifest <genome.fa> <manifest.json> --genome-build GRCh38",
                ))
            info.append(Issue(
                "ref.genome",
                f"observed {facts['contig_count']} contigs, {facts['total_bases']:,} bp, "
                f".fai sha256 {facts['fai_sha256'][:16]}…; FASTA/.fai/.dict consistent; "
                f"BWA index verified ({bfacts.get('bwa_index_n_seqs', '?')} contigs){manifest_note}",
            ))

    # ── bgzip / tabix indexes usable, header build agrees with declared build ──
    vcf_resources = [("ref.dbsnp", ref.get("dbsnp")), ("ref.gnomad", ref.get("gnomad")),
                     ("annotation.clinvar.vcf", (ann.get("clinvar") or {}).get("vcf")),
                     ("ref.pon", ref.get("pon"))]
    cosmic = _str_path((ann.get("cosmic") or {}).get("vcf"))
    if cosmic and os.path.isfile(cosmic):
        vcf_resources.append(("annotation.cosmic.vcf", cosmic))
    declared = resource_inspect.normalize_build(build)
    for field_name, path in vcf_resources:
        path = _str_path(path)
        if not path or not os.path.isfile(path):
            continue
        e, w = resource_inspect.check_vcf_index(path)
        add(errors, f"{field_name} (index)", e)
        add(warnings, f"{field_name} (index)", w)
        header = resource_inspect.read_vcf_header(path)
        header_build = resource_inspect.build_family_from_text(header.get("reference"))
        if declared and header_build and header_build != declared:
            errors.append(Issue(
                field_name,
                f"VCF header says reference '{header.get('reference')}' ({header_build}) but the declared "
                f"genome build is {declared}",
            ))
        if field_name in ("annotation.clinvar.vcf", "ref.dbsnp") and not header.get("error"):
            bits = [f"{k}={header[k]}" for k in ("file_date", "dbsnp_build_id", "reference", "source") if header.get(k)]
            note = " (custom ClinVar file; source of the report's ClinVar column)" if field_name.endswith("clinvar.vcf") else ""
            info.append(Issue(field_name, "header " + ", ".join(bits) + note))

    # ── VEP cache observed from the cache itself ──
    cache_dir = _str_path(vep.get("cache_dir"))
    if cache_dir and os.path.isdir(cache_dir):
        pinned = resource_inspect.vep_cache_version_from_extra(vep.get("extra"))
        e, w, cfacts = resource_inspect.inspect_vep_cache(cache_dir, pinned, build)
        add(errors, "annotation.vep.cache_dir (cache/assembly)", e)
        add(warnings, "annotation.vep.cache_dir (cache/assembly)", w)
        cinfo = cfacts.get("info") or {}
        if cinfo and not e:
            srcs = ", ".join(f"{k[len('source_'):]} {cinfo[k]}" for k in
                             ("source_gencode", "source_ClinVar", "source_COSMIC", "source_dbSNP") if k in cinfo)
            info.append(Issue(
                "annotation.vep.cache_dir",
                f"observed cache {cfacts['version_dir']}: info.txt assembly {cinfo.get('assembly')}; sources: {srcs} "
                "(the cache's ClinVar is an annotation dependency, not the report's ClinVar source)",
            ))

    return errors, warnings, info


def check_cosmic_availability(config: dict):
    """
    COSMIC is optional. Returns (available: bool, configured_path: str|None)
    without raising or failing — callers decide whether to warn/report.
    """
    cosmic_path = ((config.get("annotation") or {}).get("cosmic") or {}).get("vcf")
    available = bool(cosmic_path) and os.path.isfile(cosmic_path)
    return available, cosmic_path


# ── Panel / reference build compatibility (deterministic, best-effort) ──
def _first_field(path: str, skip_prefixes=("#",)):
    try:
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith(skip_prefixes):
                    continue
                return line.split("\t")[0]
    except OSError:
        return None
    return None


def check_panel_reference_compatibility(panel_bed_path, genome_path) -> list:
    """
    Deterministic, cheap compatibility check: compares chr-prefix naming
    convention between the panel BED and the reference genome's .fai index.
    Never blocks a run just because compatibility *can't* be determined
    (e.g. missing .fai, which is separately reported by validate_resources)
    — only when it CAN be determined and finds a mismatch.
    """
    if not panel_bed_path or not genome_path:
        return []
    if not isinstance(panel_bed_path, (str, os.PathLike)) or not isinstance(genome_path, (str, os.PathLike)):
        return []  # malformed config values are reported by validate_config_ranges, not here
    panel_first = _first_field(panel_bed_path)
    fai_path = genome_path + ".fai"
    genome_first = _first_field(fai_path) if os.path.isfile(fai_path) else None
    if panel_first is None or genome_first is None:
        return []
    panel_chr = panel_first.lower().startswith("chr")
    genome_chr = genome_first.lower().startswith("chr")
    if panel_chr == genome_chr:
        return []
    return [Issue(
        "panel.bed / ref.genome",
        "Panel BED and reference genome use incompatible chromosome naming conventions "
        f"(panel contig '{panel_first}' vs. genome contig '{genome_first}' per {fai_path}). "
        "This would silently produce zero/near-zero coverage rather than a clear error.",
    )]


def check_genome_build_naming(config: dict) -> list:
    """
    Deterministic heuristic: flags an obvious contradiction between the
    declared annotation.vep.genome_build and the reference genome file path
    (e.g. genome_build=GRCh38 but ref.genome path clearly says hg19/GRCh37).
    Does not attempt to *infer* a build when the path gives no clear hint —
    only fires on an unambiguous, literal contradiction.
    """
    if not isinstance(config, dict):
        return []
    vep = (config.get("annotation") or {}).get("vep") or {}
    build = str(vep.get("genome_build") or "").strip().upper()
    genome_path = str((config.get("ref") or {}).get("genome") or "").lower()
    if not build or not genome_path:
        return []
    grch38_hint = any(t in genome_path for t in ("hg38", "grch38"))
    grch37_hint = any(t in genome_path for t in ("hg19", "grch37"))
    if build in ("GRCH38", "HG38") and grch37_hint and not grch38_hint:
        return [Issue(
            "annotation.vep.genome_build",
            f"genome_build='{vep.get('genome_build')}' but ref.genome path suggests GRCh37/hg19: {genome_path}",
        )]
    if build in ("GRCH37", "HG19") and grch38_hint and not grch37_hint:
        return [Issue(
            "annotation.vep.genome_build",
            f"genome_build='{vep.get('genome_build')}' but ref.genome path suggests GRCh38/hg38: {genome_path}",
        )]
    return []


# ── Top-level entry point ────────────────────────────────────────────────
def run_preflight_checks(config: dict, samples_df, base_dir: str | None = None,
                         deep_resource_checks: bool = True) -> PreflightReport:
    errors: list = []
    warnings: list = []
    info: list = []

    errors.extend(validate_config_keys(config))
    errors.extend(validate_config_ranges(config))
    errors.extend(validate_samples_sheet(samples_df))
    errors.extend(validate_fastq_files(samples_df, base_dir=base_dir))
    errors.extend(validate_resources(config))
    _panel_bed = _get_dotted(config, "panel.bed") if isinstance(config, dict) else None
    _ref_genome = _get_dotted(config, "ref.genome") if isinstance(config, dict) else None
    errors.extend(check_panel_reference_compatibility(
        None if _panel_bed is _MISSING else _panel_bed,
        None if _ref_genome is _MISSING else _ref_genome,
    ))
    errors.extend(check_genome_build_naming(config))

    if deep_resource_checks:
        deep_errors, deep_warnings, deep_info = validate_runtime_resources(config, base_dir=base_dir)
        errors.extend(deep_errors)
        warnings.extend(deep_warnings)
        info.extend(deep_info)

    if isinstance(config, dict):
        cosmic_available, cosmic_path = check_cosmic_availability(config)
        if cosmic_path and not cosmic_available:
            warnings.append(Issue(
                "annotation.cosmic.vcf",
                f"COSMIC resource configured but not found on disk: {cosmic_path} — "
                "COSMIC annotation will be skipped for this run.",
            ))
        elif not cosmic_path:
            warnings.append(Issue(
                "annotation.cosmic.vcf",
                "COSMIC not configured — COSMIC annotation will be skipped for this run (optional).",
            ))

    return PreflightReport(errors=errors, warnings=warnings, info=info)


def format_report(report: PreflightReport, title: str = "GenRichi Phase 1 Pre-flight Validation") -> str:
    lines = [title, "=" * len(title)]
    info = getattr(report, "info", None) or []
    if not report.errors and not report.warnings:
        lines.append("All pre-flight checks passed.")
    if report.errors:
        lines.append(f"\n{len(report.errors)} error(s) — pipeline will NOT run:")
        for i, issue in enumerate(report.errors, 1):
            lines.append(f"  [{i}] {issue.field}: {issue.message}")
    if report.warnings:
        lines.append(f"\n{len(report.warnings)} warning(s) — pipeline will continue:")
        for i, issue in enumerate(report.warnings, 1):
            lines.append(f"  [{i}] {issue.field}: {issue.message}")
    if info:
        lines.append("\nObserved resources / runtime facts:")
        for issue in info:
            lines.append(f"  - {issue.field}: {issue.message}")
    return "\n".join(lines) + "\n"


def _main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="GenRichi pre-flight helpers")
    ap.add_argument("--write-reference-manifest", nargs=2, metavar=("GENOME_FA", "OUT_JSON"),
                    help="fingerprint the reference (contig set, .fai/.dict hashes) for later comparison")
    ap.add_argument("--genome-build", default="GRCh38")
    ap.add_argument("--description", default=None)
    args = ap.parse_args(argv)
    if not args.write_reference_manifest:
        ap.print_help()
        return 2
    genome, out = args.write_reference_manifest
    try:
        manifest = resource_inspect.build_reference_manifest(genome, args.genome_build, args.description)
    except resource_inspect.ResourceError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    with open(out, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print(f"wrote {out}: {manifest['reference']['contig_count']} contigs, "
          f"{manifest['reference']['total_bases']:,} bp")
    return 0


if __name__ == "__main__":
    sys.exit(_main())
