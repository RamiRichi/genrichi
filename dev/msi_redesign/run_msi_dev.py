"""
GenRichi MSI Methodology Redesign -- DEVELOPMENT PROTOTYPE. NOT PRODUCTION CODE.

Phase C/D orchestrator: real marker file + real tumor/normal BAMs ->
LocusObservation objects -> the already-tested classify_sample() core.

This module does NOT implement any classification logic itself -- it only
wires scan_markers.py's output and bam_extract.py's real-BAM extraction
into msi_locus_model.classify_sample(), and writes an auditable per-locus
TSV/JSON so every number in the final call can be traced back to specific
reads at a specific locus.

This is a DEVELOPMENT ANALYSIS, not a clinical validation run. See
dev/msi_redesign/README.md and the Phase F decision-gate report for what
this run does and does not establish.

Usage:
    python3 run_msi_dev.py --markers dev/msi_redesign/output/msi_markers.v1.tsv \
        --tumor-bam /path/tumor.final.bam --normal-bam /path/normal.final.bam \
        --out-dir dev/msi_redesign/output --sample-id SAMPLE_001 --samtools samtools
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from dataclasses import asdict
from datetime import datetime, timezone

DEV_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DEV_DIR)

import bam_extract as be  # noqa: E402
import msi_locus_model as mlm  # noqa: E402


def read_markers(path: str):
    markers = []
    with open(path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            markers.append({
                "marker_id": row["marker_id"], "chrom": row["chrom"],
                "start": int(row["start"]), "end": int(row["end"]),
                "gene": row["gene"], "kind": row["kind"],
            })
    return markers


def build_observations(samtools, tumor_bam, normal_bam, markers, min_mapq, min_flank):
    """Real I/O: for every marker, extract real repeat-length lists from both
    BAMs via bam_extract.locus_lengths(). Returns (observations, per_locus_meta)
    where per_locus_meta records raw depth/coordinates for the audit trail
    even for loci later marked NOT_CALLABLE by the classification core."""
    observations = []
    per_locus_meta = []
    for m in markers:
        tumor_lengths = be.locus_lengths(
            samtools, tumor_bam, m["chrom"], m["start"], m["end"],
            min_mapq=min_mapq, min_flank=min_flank,
        )
        normal_lengths = be.locus_lengths(
            samtools, normal_bam, m["chrom"], m["start"], m["end"],
            min_mapq=min_mapq, min_flank=min_flank,
        )
        observations.append(mlm.LocusObservation(
            marker_id=m["marker_id"], tumor_lengths=tumor_lengths, normal_lengths=normal_lengths,
        ))
        per_locus_meta.append({
            "marker_id": m["marker_id"], "chrom": m["chrom"], "start": m["start"], "end": m["end"],
            "gene": m["gene"], "kind": m["kind"],
            "raw_tumor_depth": len(tumor_lengths), "raw_normal_depth": len(normal_lengths),
        })
    return observations, per_locus_meta


def write_per_locus_tsv(out_path, per_locus_meta, classifications):
    meta_by_id = {m["marker_id"]: m for m in per_locus_meta}
    with open(out_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("marker_id\tchrom\tstart\tend\tgene\tkind\tstatus\t"
                  "raw_tumor_depth\traw_normal_depth\ttumor_depth\tnormal_depth\t"
                  "tumor_mode\tnormal_mode\tshift_units\tdistance\n")
        for lc in classifications:
            meta = meta_by_id[lc.marker_id]
            fh.write(
                f"{lc.marker_id}\t{meta['chrom']}\t{meta['start']}\t{meta['end']}\t"
                f"{meta['gene']}\t{meta['kind']}\t{lc.status}\t"
                f"{meta['raw_tumor_depth']}\t{meta['raw_normal_depth']}\t"
                f"{lc.tumor_depth}\t{lc.normal_depth}\t"
                f"{lc.tumor_mode if lc.tumor_mode is not None else ''}\t"
                f"{lc.normal_mode if lc.normal_mode is not None else ''}\t"
                f"{lc.shift_units if lc.shift_units is not None else ''}\t"
                f"{lc.distance if lc.distance is not None else ''}\n"
            )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--markers", required=True)
    ap.add_argument("--tumor-bam", required=True)
    ap.add_argument("--normal-bam", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--min-mapq", type=int, default=be.MIN_MAPQ)
    ap.add_argument("--min-flank", type=int, default=be.MIN_FLANK)
    args = ap.parse_args(argv)

    os.makedirs(args.out_dir, exist_ok=True)
    markers = read_markers(args.markers)

    observations, per_locus_meta = build_observations(
        args.samtools, args.tumor_bam, args.normal_bam, markers,
        args.min_mapq, args.min_flank,
    )

    result = mlm.classify_sample(observations)

    per_locus_path = os.path.join(args.out_dir, f"{args.sample_id}.msi_dev.per_locus.tsv")
    write_per_locus_tsv(per_locus_path, per_locus_meta, result.loci)

    summary = {
        "sample_id": args.sample_id,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "analysis_type": "DEVELOPMENT_ANALYSIS_NOT_CLINICAL_VALIDATION",
        "markers_file": args.markers,
        "tumor_bam": args.tumor_bam,
        "normal_bam": args.normal_bam,
        "read_filter_params": {"min_mapq": args.min_mapq, "min_flank": args.min_flank},
        "classification_params": {
            "min_locus_depth": mlm.MIN_LOCUS_DEPTH,
            "instability_distance_cutoff": mlm.INSTABILITY_DISTANCE_CUTOFF,
            "min_shift_units": mlm.MIN_SHIFT_UNITS,
            "min_informative_loci": mlm.MIN_INFORMATIVE_LOCI,
            "msi_h_fraction_threshold": mlm.MSI_H_FRACTION_THRESHOLD,
        },
        "n_markers_in_panel": len(markers),
        "sample_status": result.sample_status,
        "n_loci_total": result.n_loci_total,
        "n_callable": result.n_callable,
        "n_unstable": result.n_unstable,
        "n_not_callable": result.n_not_callable,
        "fraction_unstable": result.fraction_unstable,
        "per_locus_tsv": per_locus_path,
    }
    summary_path = os.path.join(args.out_dir, f"{args.sample_id}.msi_dev.summary.json")
    with open(summary_path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)
        fh.write("\n")

    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
