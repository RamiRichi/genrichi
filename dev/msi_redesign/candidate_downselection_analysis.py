"""
GenRichi MSI Methodology Redesign -- DEVELOPMENT PROTOTYPE. NOT PRODUCTION CODE.

EXPLORATORY down-selection analysis for View A (WES-compatible, 93,687
loci) and View B (WGS/master filtered universe, 1,208,169 loci), run in
parallel. Produces distributions and scenario counts only -- it does NOT
select a final threshold, does NOT write a new filtered panel, and the
result is explicitly NOT to be called a "validated MSI panel" (see the
module-level report this script feeds).

Reuses already-computed, already-documented local annotations read-only:
  - tier2_annotated.v1.tsv (per-locus chrom/motif/repeat_length/kind, plus
    Umap continuous mappability mean/min already computed in this session)
  - view_A_wes_compatible.v1.tsv (View A membership)
  - assembly_gaps.hg38.bed (NEW this step: a fresh reference-only scan for
    N-runs in the real production hg38.fa -- no sample/BAM data, no new
    external download, see scan_assembly_gaps.py)

No HCC1395 or any other sample/BAM data is read anywhere in this script.

Usage:
    python3 candidate_downselection_analysis.py \
        --annotated dev/msi_redesign/output/tier2_annotated.v1.tsv \
        --view-a dev/msi_redesign/output/view_A_wes_compatible.v1.tsv \
        --gap-distance-tsv dev/msi_redesign/reference_annotations/view_{A,B}_gap_distance.tsv \
        --out-dir dev/msi_redesign/output
"""

from __future__ import annotations

import csv
import json
import sys
from collections import Counter, defaultdict

SPACING_THRESHOLDS_BP = [20, 50, 100, 200]
GAP_PROXIMITY_THRESHOLDS_BP = [0, 10, 50, 100]
REPEAT_LENGTH_BANDS = [
    (0, 9, "5-9bp"), (10, 14, "10-14bp"), (15, 20, "15-20bp"),
    (21, 30, "21-30bp"), (31, 10**9, "31bp+"),
]
UMAP_MEAN_BANDS = [
    (1.0, 1.0, "1.0_exactly"), (0.9, 1.0, "0.9-1.0_exclusive"),
    (0.7, 0.9, "0.7-0.9"), (0.0, 0.7, "below_0.7"),
]


def repeat_length_band(length: int) -> str:
    for lo, hi, label in REPEAT_LENGTH_BANDS:
        if lo <= length <= hi:
            return label
    return "unbanded"


def umap_mean_band(mean_score) -> str:
    if mean_score is None or mean_score == "":
        return "no_data"
    mean_score = float(mean_score)
    if mean_score >= 1.0:
        return "1.0_exactly"
    for lo, hi, label in UMAP_MEAN_BANDS[1:]:
        if lo <= mean_score < hi:
            return label
    return "unbanded"


def nearest_neighbor_min_gaps(sorted_positions):
    """sorted_positions: list of (start, end) sorted by start, same
    chromosome. Returns a list of min-gap-to-neighbor (bp) per locus,
    None for a locus with no neighbor on that chromosome (gaps can be
    negative if intervals overlap)."""
    n = len(sorted_positions)
    gaps = []
    for i in range(n):
        candidates = []
        if i > 0:
            candidates.append(sorted_positions[i][0] - sorted_positions[i - 1][1])
        if i < n - 1:
            candidates.append(sorted_positions[i + 1][0] - sorted_positions[i][1])
        gaps.append(min(candidates) if candidates else None)
    return gaps


def spacing_scenario_counts(min_gaps_by_chrom: dict, thresholds=SPACING_THRESHOLDS_BP):
    """For each threshold T, counts how many loci have a neighbor within
    T bp (i.e. would be part of a 'redundant cluster' under that
    definition) -- NOT a decision to remove them, just a count."""
    all_gaps = [g for gaps in min_gaps_by_chrom.values() for g in gaps if g is not None]
    counts = {}
    for t in thresholds:
        counts[f"n_loci_with_neighbor_within_{t}bp"] = sum(1 for g in all_gaps if g <= t)
    counts["n_loci_with_no_neighbor_on_chromosome"] = sum(
        1 for gaps in min_gaps_by_chrom.values() for g in gaps if g is None
    )
    counts["n_total"] = sum(len(gaps) for gaps in min_gaps_by_chrom.values())
    return counts


def gap_proximity_scenario_counts(distances: list, thresholds=GAP_PROXIMITY_THRESHOLDS_BP):
    counts = {}
    for t in thresholds:
        counts[f"n_loci_within_{t}bp_of_assembly_gap"] = sum(1 for d in distances if d is not None and d <= t)
    counts["n_total_with_distance_data"] = sum(1 for d in distances if d is not None)
    return counts


def load_view_marker_ids(path: str) -> set:
    ids = set()
    with open(path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            ids.add(row["marker_id"])
    return ids


def load_annotated_rows(path: str):
    with open(path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            yield row


def analyze_view(marker_ids: set, annotated_rows_iter, gap_distances: dict = None):
    kind_counts = Counter()
    motif_len_counts = Counter()
    repeat_length_band_counts = Counter()
    repeat_length_hist = Counter()
    chrom_counts = Counter()
    umap_band_counts = Counter()
    positions_by_chrom = defaultdict(list)
    marker_by_position = {}

    for row in annotated_rows_iter:
        if row["marker_id"] not in marker_ids:
            continue
        kind_counts[row["kind"]] += 1
        motif_len_counts[len(row["motif"])] += 1
        length = int(row["ref_repeat_length"])
        repeat_length_hist[length] += 1
        repeat_length_band_counts[repeat_length_band(length)] += 1
        chrom_counts[row["chrom"]] += 1
        umap_band_counts[umap_mean_band(row.get("umap_multiread_mean", ""))] += 1
        start, end = int(row["start"]), int(row["end"])
        positions_by_chrom[row["chrom"]].append((start, end))
        marker_by_position[(row["chrom"], start, end)] = row["marker_id"]

    min_gaps_by_chrom = {}
    for chrom, positions in positions_by_chrom.items():
        positions.sort()
        min_gaps_by_chrom[chrom] = nearest_neighbor_min_gaps(positions)

    result = {
        "n_total": sum(kind_counts.values()),
        "kind_counts": dict(kind_counts),
        "motif_length_counts": dict(sorted(motif_len_counts.items())),
        "repeat_length_band_counts": dict(repeat_length_band_counts),
        "repeat_length_histogram_summary": {
            "min": min(repeat_length_hist) if repeat_length_hist else None,
            "max": max(repeat_length_hist) if repeat_length_hist else None,
            "n_distinct_lengths": len(repeat_length_hist),
        },
        "chrom_counts": dict(chrom_counts),
        "umap_multiread_mean_band_counts": dict(umap_band_counts),
        "spacing_scenario_counts": spacing_scenario_counts(min_gaps_by_chrom),
    }
    if gap_distances is not None:
        distances = [gap_distances.get(mid) for mid in marker_ids if mid in gap_distances]
        result["gap_proximity_scenario_counts"] = gap_proximity_scenario_counts(distances)
        result["n_loci_with_gap_distance_data"] = len(distances)
    return result


def load_gap_distances(path: str) -> dict:
    """bedtools closest -d output -> {marker_id: distance_bp}."""
    result = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            marker_id = parts[3]
            try:
                distance = int(parts[-1])
            except ValueError:
                continue
            if distance >= 0:
                result[marker_id] = distance
    return result


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--annotated", required=True)
    ap.add_argument("--view-a-ids", required=True)
    ap.add_argument("--filtered-panel-ids", required=True, help="view B membership (= master filtered panel)")
    ap.add_argument("--gap-distance-view-a", required=True)
    ap.add_argument("--gap-distance-view-b", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    view_a_ids = load_view_marker_ids(args.view_a_ids)
    view_b_ids = load_view_marker_ids(args.filtered_panel_ids)

    gap_dist_a = load_gap_distances(args.gap_distance_view_a)
    gap_dist_b = load_gap_distances(args.gap_distance_view_b)

    # single pass over the (large) annotated file, reused for both views
    rows_a, rows_b = [], []
    for row in load_annotated_rows(args.annotated):
        mid = row["marker_id"]
        if mid in view_a_ids:
            rows_a.append(row)
        if mid in view_b_ids:
            rows_b.append(row)

    result = {
        "analysis_type": "EXPLORATORY_DOWNSELECTION_NOT_A_VALIDATED_PANEL",
        "view_A_wes_compatible": analyze_view(view_a_ids, rows_a, gap_dist_a),
        "view_B_wgs": analyze_view(view_b_ids, rows_b, gap_dist_b),
        "hcc1395_or_any_sample_data_used": False,
        "master_panel_or_views_modified": False,
        "final_panel_selected": False,
    }

    with open(args.out, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(result, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
