"""
GenRichi MSI Methodology Redesign -- DEVELOPMENT PROTOTYPE. NOT PRODUCTION CODE.

READ-ONLY analysis of the Tier 2 genome-wide candidate universe
(msi_markers.v2.tsv). Does NOT filter, select, or write a final marker
panel -- it only computes descriptive statistics used to ground the
filtering-strategy design report in real numbers from the real candidate
set, instead of assumed/estimated figures.

Usage:
    python3 analyze_candidate_universe.py --markers dev/msi_redesign/output/msi_markers.v2.tsv \
        --out-dir dev/msi_redesign/output
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone


def analyze(markers_path: str):
    mono_length_hist = Counter()
    str_motif_len_hist = Counter()
    str_units_hist = Counter()
    chrom_counts = Counter()
    kind_counts = Counter()

    # nearest-neighbor distance on the same chromosome (proxy for "isolated
    # vs. clustered in a repetitive context" without a mappability track)
    by_chrom_positions = defaultdict(list)

    with open(markers_path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            chrom = row["chrom"]
            start, end = int(row["start"]), int(row["end"])
            kind = row["kind"]
            motif = row["motif"]
            length = int(row["ref_repeat_length"])

            chrom_counts[chrom] += 1
            kind_counts[kind] += 1
            by_chrom_positions[chrom].append((start, end))

            if kind == "mononucleotide":
                mono_length_hist[length] += 1
            else:
                str_motif_len_hist[len(motif)] += 1
                str_units_hist[length] += 1

    # nearest-neighbor gap distribution (sampled per chromosome, sorted once)
    gap_hist = Counter()
    n_within_50bp_of_neighbor = 0
    n_total = 0
    for chrom, positions in by_chrom_positions.items():
        positions.sort()
        for i in range(len(positions)):
            n_total += 1
            gaps = []
            if i > 0:
                gaps.append(positions[i][0] - positions[i - 1][1])
            if i < len(positions) - 1:
                gaps.append(positions[i + 1][0] - positions[i][1])
            if gaps:
                min_gap = min(gaps)
                bucket = ("<0_overlapping" if min_gap < 0 else
                          "0-10" if min_gap <= 10 else
                          "11-50" if min_gap <= 50 else
                          "51-200" if min_gap <= 200 else
                          "201-1000" if min_gap <= 1000 else ">1000")
                gap_hist[bucket] += 1
                if min_gap <= 50:
                    n_within_50bp_of_neighbor += 1

    # cumulative counts for candidate homopolymer-length caps (sensitivity
    # analysis for where to set a max-length threshold)
    mono_cumulative_by_cap = {}
    for cap in (10, 12, 15, 20, 25, 30, 40, 50):
        mono_cumulative_by_cap[cap] = sum(n for length, n in mono_length_hist.items() if length <= cap)

    str_cumulative_units_by_cap = {}
    for cap in (5, 6, 8, 10, 15, 20):
        str_cumulative_units_by_cap[cap] = sum(n for length, n in str_units_hist.items() if length <= cap)

    return {
        "n_total": sum(kind_counts.values()),
        "kind_counts": dict(kind_counts),
        "chrom_counts": dict(chrom_counts),
        "mononucleotide_length_histogram": dict(sorted(mono_length_hist.items())),
        "mononucleotide_cumulative_count_by_max_length_cap": mono_cumulative_by_cap,
        "str_motif_length_histogram": dict(sorted(str_motif_len_hist.items())),
        "str_units_histogram": dict(sorted(str_units_hist.items())),
        "str_cumulative_count_by_max_units_cap": str_cumulative_units_by_cap,
        "nearest_neighbor_gap_bp_histogram": dict(gap_hist),
        "n_loci_within_50bp_of_another_locus": n_within_50bp_of_neighbor,
        "fraction_loci_within_50bp_of_another_locus": round(n_within_50bp_of_neighbor / max(1, n_total), 4),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--markers", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args(argv)

    os.makedirs(args.out_dir, exist_ok=True)
    stats = analyze(args.markers)
    stats["generated_at_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    stats["source_markers_file"] = args.markers
    stats["analysis_type"] = "READ_ONLY_DESCRIPTIVE_STATISTICS_NOT_A_FILTERED_PANEL"

    out_path = os.path.join(args.out_dir, "candidate_universe_analysis.json")
    with open(out_path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(stats, fh, indent=2, sort_keys=True)
        fh.write("\n")

    print(json.dumps(stats, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
