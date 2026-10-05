"""
GenRichi MSI Methodology Redesign -- DEVELOPMENT PROTOTYPE. NOT PRODUCTION CODE.

Coverage audit ONLY -- does not run any classification. Determines whether
an already-existing local tumor/normal BAM pair provides adequate coverage
across the combined (Tier 1 + Tier 2) marker set, without downloading or
sourcing any new data.

Two different measurement methods are used deliberately, for performance
reasons, and this is documented rather than hidden:
  - Tier 1 (5 markers): EXACT per-read measurement via bam_extract.py's
    already-tested locus_lengths() (real samtools view + CIGAR-aware
    filtering, MAPQ>=20, flank>=5bp) -- identical logic to run_msi_dev.py's
    actual extraction, since there are only 5 loci this is cheap.
  - Tier 2 (potentially hundreds of thousands to millions of markers):
    APPROXIMATE coverage via `samtools bedcov -Q 20`, which sums per-base
    pileup depth over each region in one efficient pass per chromosome
    (MAPQ>=20 filtered, but WITHOUT the exact CIGAR-based flank/duplicate
    filtering bam_extract.py applies). Average depth = summed depth /
    region width, thresholded at the same MIN_LOCUS_DEPTH=20 used by
    msi_locus_model.py. This is explicitly an approximation used only to
    scope the coverage audit at this scale -- it is NOT used for any
    classification, and any future real classification run must use the
    exact bam_extract.py path, not this audit's numbers.

Only chromosomes with mapped-read counts (samtools idxstats) at or above
--coverage-chrom-threshold are treated as "covered" -- this distinguishes
real target coverage from off-target/spillover reads (empirically: in the
one local dataset, the covered chromosome has ~8M mapped reads vs ~1-6K
spillover elsewhere).

Usage:
    python3 coverage_audit.py --markers dev/msi_redesign/output/msi_markers.dev_combined.v1.tsv \
        --tumor-bam <path> --normal-bam <path> --out-dir dev/msi_redesign/output \
        --sample-id <id> --samtools samtools
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from datetime import datetime, timezone

DEV_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DEV_DIR)

import bam_extract as be  # noqa: E402

MIN_LOCUS_DEPTH = 20
DEFAULT_COVERAGE_CHROM_THRESHOLD = 100_000


def read_markers(path: str):
    markers = []
    with open(path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            markers.append({
                "marker_id": row["marker_id"], "chrom": row["chrom"],
                "start": int(row["start"]), "end": int(row["end"]),
                "tier": int(row["tier"]),
            })
    return markers


def idxstats_mapped_counts(samtools: str, bam_path: str) -> dict:
    out = subprocess.run([samtools, "idxstats", bam_path], capture_output=True, text=True, timeout=120)
    if out.returncode != 0:
        raise RuntimeError(f"samtools idxstats failed for {bam_path}: {out.stderr.strip()}")
    counts = {}
    for line in out.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 3:
            counts[parts[0]] = int(parts[2])
    return counts


def covered_chromosomes(samtools: str, bam_path: str, threshold: int) -> set:
    counts = idxstats_mapped_counts(samtools, bam_path)
    return {chrom for chrom, n in counts.items() if n >= threshold}


def bedcov_average_depths(samtools: str, bam_path: str, chrom: str, regions, min_mapq: int = 20):
    """regions: list of (marker_id, start, end) 0-based half-open, all on
    the same chromosome. Returns {marker_id: approx_average_depth}."""
    if not regions:
        return {}
    bed_lines = "\n".join(f"{chrom}\t{start}\t{end}\t{marker_id}" for marker_id, start, end in regions)
    proc = subprocess.run(
        [samtools, "bedcov", "-Q", str(min_mapq), "/dev/stdin", bam_path],
        input=bed_lines, capture_output=True, text=True, timeout=600,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"samtools bedcov failed on {chrom}: {proc.stderr.strip()}")
    result = {}
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        marker_id, start, end, depth_sum = parts[3], int(parts[1]), int(parts[2]), int(parts[4])
        width = max(1, end - start)
        result[marker_id] = depth_sum / width
    return result


def audit(samtools, tumor_bam, normal_bam, markers, coverage_chrom_threshold, min_mapq, min_flank):
    tumor_covered = covered_chromosomes(samtools, tumor_bam, coverage_chrom_threshold)
    normal_covered = covered_chromosomes(samtools, normal_bam, coverage_chrom_threshold)
    jointly_covered = tumor_covered & normal_covered

    tier1 = [m for m in markers if m["tier"] == 1]
    tier2 = [m for m in markers if m["tier"] == 2]

    per_marker = []

    # Tier 1: exact measurement (cheap, only 5 loci)
    for m in tier1:
        on_covered_chrom = m["chrom"] in jointly_covered
        tumor_lengths = be.locus_lengths(samtools, tumor_bam, m["chrom"], m["start"], m["end"],
                                          min_mapq=min_mapq, min_flank=min_flank)
        normal_lengths = be.locus_lengths(samtools, normal_bam, m["chrom"], m["start"], m["end"],
                                           min_mapq=min_mapq, min_flank=min_flank)
        per_marker.append({
            "marker_id": m["marker_id"], "tier": 1, "chrom": m["chrom"],
            "on_covered_chromosome": on_covered_chrom,
            "measurement_method": "exact",
            "tumor_depth": len(tumor_lengths), "normal_depth": len(normal_lengths),
        })

    # Tier 2: approximate measurement, restricted to jointly-covered chromosomes
    tier2_by_chrom = {}
    for m in tier2:
        tier2_by_chrom.setdefault(m["chrom"], []).append(m)

    for chrom, chrom_markers in tier2_by_chrom.items():
        if chrom not in jointly_covered:
            for m in chrom_markers:
                per_marker.append({
                    "marker_id": m["marker_id"], "tier": 2, "chrom": chrom,
                    "on_covered_chromosome": False, "measurement_method": "approximate",
                    "tumor_depth": 0, "normal_depth": 0,
                })
            continue
        regions = [(m["marker_id"], m["start"], m["end"]) for m in chrom_markers]
        tumor_depths = bedcov_average_depths(samtools, tumor_bam, chrom, regions, min_mapq)
        normal_depths = bedcov_average_depths(samtools, normal_bam, chrom, regions, min_mapq)
        for m in chrom_markers:
            per_marker.append({
                "marker_id": m["marker_id"], "tier": 2, "chrom": chrom,
                "on_covered_chromosome": True, "measurement_method": "approximate",
                "tumor_depth": round(tumor_depths.get(m["marker_id"], 0.0), 2),
                "normal_depth": round(normal_depths.get(m["marker_id"], 0.0), 2),
            })

    n_total = len(markers)
    n_on_covered_chrom = sum(1 for r in per_marker if r["on_covered_chromosome"])
    n_tumor_ok = sum(1 for r in per_marker if r["tumor_depth"] >= MIN_LOCUS_DEPTH)
    n_normal_ok = sum(1 for r in per_marker if r["normal_depth"] >= MIN_LOCUS_DEPTH)
    n_callable = sum(1 for r in per_marker if r["tumor_depth"] >= MIN_LOCUS_DEPTH and r["normal_depth"] >= MIN_LOCUS_DEPTH)
    tier1_records = [r for r in per_marker if r["tier"] == 1]
    n_callable_tier1 = sum(1 for r in tier1_records
                            if r["tumor_depth"] >= MIN_LOCUS_DEPTH and r["normal_depth"] >= MIN_LOCUS_DEPTH)

    summary = {
        "coverage_chrom_threshold_mapped_reads": coverage_chrom_threshold,
        "min_locus_depth_threshold": MIN_LOCUS_DEPTH,
        "tumor_covered_chromosomes": sorted(tumor_covered),
        "normal_covered_chromosomes": sorted(normal_covered),
        "jointly_covered_chromosomes": sorted(jointly_covered),
        "n_markers_total": n_total,
        "n_markers_tier1_promega": len(tier1),
        "n_markers_tier2_genome_scan": len(tier2),
        "n_markers_on_covered_chromosomes": n_on_covered_chrom,
        "n_markers_tumor_depth_ge_20": n_tumor_ok,
        "n_markers_normal_depth_ge_20": n_normal_ok,
        "n_markers_callable_both_ge_20": n_callable,
        "n_markers_tier1_callable": n_callable_tier1,
        "n_markers_tier1_total": len(tier1),
        "dataset_suitable_for_meaningful_dev_run": n_callable >= 5,  # matches MIN_INFORMATIVE_LOCI in msi_locus_model.py
    }
    return summary, per_marker


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--markers", required=True)
    ap.add_argument("--tumor-bam", required=True)
    ap.add_argument("--normal-bam", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--coverage-chrom-threshold", type=int, default=DEFAULT_COVERAGE_CHROM_THRESHOLD)
    ap.add_argument("--min-mapq", type=int, default=be.MIN_MAPQ)
    ap.add_argument("--min-flank", type=int, default=be.MIN_FLANK)
    args = ap.parse_args(argv)

    os.makedirs(args.out_dir, exist_ok=True)
    markers = read_markers(args.markers)
    summary, per_marker = audit(
        args.samtools, args.tumor_bam, args.normal_bam, markers,
        args.coverage_chrom_threshold, args.min_mapq, args.min_flank,
    )

    summary_out = dict(summary)
    summary_out.update({
        "sample_id": args.sample_id,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "analysis_type": "COVERAGE_AUDIT_ONLY_NOT_CLASSIFICATION",
        "markers_file": args.markers,
        "tumor_bam": args.tumor_bam,
        "normal_bam": args.normal_bam,
    })
    summary_path = os.path.join(args.out_dir, f"{args.sample_id}.coverage_audit.summary.json")
    with open(summary_path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(summary_out, fh, indent=2, sort_keys=True)
        fh.write("\n")

    detail_path = os.path.join(args.out_dir, f"{args.sample_id}.coverage_audit.per_marker.tsv")
    with open(detail_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("marker_id\ttier\tchrom\ton_covered_chromosome\tmeasurement_method\ttumor_depth\tnormal_depth\n")
        for r in per_marker:
            fh.write(f"{r['marker_id']}\t{r['tier']}\t{r['chrom']}\t{r['on_covered_chromosome']}\t"
                      f"{r['measurement_method']}\t{r['tumor_depth']}\t{r['normal_depth']}\n")

    print(json.dumps(summary_out, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
