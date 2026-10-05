"""
GenRichi MSI Methodology Redesign -- DEVELOPMENT PROTOTYPE. NOT PRODUCTION CODE.

READ-ONLY feasibility analysis of the Tier 2 genome-wide candidate universe
against three independently sourced GRCh38 reference annotation tracks:
Umap mappability (single-read BED + continuous multi-read bedGraph),
UCSC RepeatMasker, and UCSC segmental duplications.

This module does NOT select a final threshold, does NOT exclude any locus,
and does NOT write a filtered marker panel. It only aggregates pre-computed
bedtools overlap outputs (produced by shell/bedtools, not by this module --
bedtools is the appropriate tool at this scale, ~2.68M loci) into the
descriptive statistics needed to make an evidence-based filtering decision
later.

Terminology:
  - "locus" = the microsatellite candidate itself (chrom,start,end from the
    Tier 2 scan).
  - "flank" = the 100bp region immediately outside the locus on each side
    (matching the real production sequencing read length, 100bp, confirmed
    against the actual HCC1395 BAMs).
  - "window" = locus + both flanks combined into one interval, used for the
    mappability/segdup overlap (a locus is only usable if a full read
    anchored in its flank can be placed uniquely).

RepeatMasker classes overlapping the LOCUS ITSELF are expected to be
Simple_repeat/Low_complexity (that is what a microsatellite structurally
is, in RepeatMasker's own taxonomy) -- this is NOT a defect and is reported
separately from RepeatMasker classes overlapping the FLANK, which indicate
a different, potentially problematic repetitive element sitting next to
(not making up) the microsatellite.
"""

from __future__ import annotations

import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

MICROSATELLITE_SELF_CLASSES = {"Simple_repeat", "Low_complexity"}

# Umap continuous-score bins for reporting the distribution before any
# threshold is chosen (per explicit instruction: do not invent a cutoff).
SCORE_BINS = [
    (0.0, 0.0, "0.0_exactly"),
    (0.0, 0.2, "0.0-0.2"),
    (0.2, 0.4, "0.2-0.4"),
    (0.4, 0.6, "0.4-0.6"),
    (0.6, 0.8, "0.6-0.8"),
    (0.8, 1.0, "0.8-1.0_exclusive"),
    (1.0, 1.0, "1.0_exactly"),
]


def bin_score(score: float) -> str:
    if score == 0.0:
        return "0.0_exactly"
    if score == 1.0:
        return "1.0_exactly"
    for lo, hi, label in SCORE_BINS[1:-1]:
        if lo <= score < hi:
            return label
    return "unbinned"  # should not happen for scores in [0,1]


def load_tier2_markers(markers_path: str):
    """Returns {marker_id: {"chrom":..., "tier":..., "kind":...}}."""
    markers = {}
    with open(markers_path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if int(row["tier"]) != 2:
                continue
            markers[row["marker_id"]] = {"chrom": row["chrom"], "kind": row["kind"]}
    return markers


def load_pair_tsv(path: str):
    """(id, value) TSV -> {id: set(values)}."""
    result = defaultdict(set)
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            result[parts[0]].add(parts[1])
    return result


def load_segdup_hits(path: str):
    """marker_id -> list of fracMatch floats (segdup 'fracMatch' column)."""
    result = defaultdict(list)
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            try:
                result[parts[0]].append(float(parts[1]))
            except ValueError:
                continue
    return result


def load_umap_singleread_coverage(path: str):
    """bedtools coverage output -> {marker_id: fraction_covered}."""
    result = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            marker_id = parts[3]
            fraction_covered = float(parts[-1])
            result[marker_id] = fraction_covered
    return result


def load_umap_multiread_scores(path: str):
    """bedtools map (-o mean,min,max) output -> {marker_id: (mean, min, max)}."""
    result = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            marker_id = parts[3]
            mean_s, min_s, max_s = parts[4], parts[5], parts[6]
            if mean_s == ".":
                continue  # no bedgraph coverage at all for this window
            result[marker_id] = (float(mean_s), float(min_s), float(max_s))
    return result


def classify_rmsk_classes(classes: set) -> str:
    """For a locus's SELF overlap: are all hit classes the expected
    microsatellite classes, a mix, or entirely something else?"""
    if not classes:
        return "no_rmsk_hit"
    non_self = classes - MICROSATELLITE_SELF_CLASSES
    if not non_self:
        return "microsatellite_class_only"
    if non_self and (classes & MICROSATELLITE_SELF_CLASSES):
        return "mixed_microsatellite_and_other"
    return "unexpected_non_microsatellite_class"


def analyze(
    markers_path: str,
    locus_self_rmsk_path: str,
    flank_rmsk_path: str,
    segdup_hits_path: str,
    umap_single_path: str,
    umap_multi_path: str,
):
    markers = load_tier2_markers(markers_path)
    n_total = len(markers)

    locus_self_rmsk = load_pair_tsv(locus_self_rmsk_path)
    flank_rmsk = load_pair_tsv(flank_rmsk_path)
    segdup_hits = load_segdup_hits(segdup_hits_path)
    umap_single = load_umap_singleread_coverage(umap_single_path)
    umap_multi = load_umap_multiread_scores(umap_multi_path)

    # --- RepeatMasker: locus-self classification ---
    self_classification_counts = Counter()
    self_class_hit_counts = Counter()  # which classes appear at the locus itself
    for marker_id in markers:
        classes = locus_self_rmsk.get(marker_id, set())
        self_classification_counts[classify_rmsk_classes(classes)] += 1
        for c in classes:
            self_class_hit_counts[c] += 1

    # --- RepeatMasker: flank classification (problematic-context candidates) ---
    flank_class_hit_counts = Counter()
    n_loci_with_any_flank_repeat = 0
    n_loci_with_flank_only_microsatellite_class = 0
    for marker_id in markers:
        classes = flank_rmsk.get(marker_id, set())
        if classes:
            n_loci_with_any_flank_repeat += 1
            for c in classes:
                flank_class_hit_counts[c] += 1
            if classes <= MICROSATELLITE_SELF_CLASSES:
                n_loci_with_flank_only_microsatellite_class += 1

    # --- Segmental duplications ---
    n_loci_in_segdup = len(segdup_hits)
    segdup_fracmatch_values = [v for vals in segdup_hits.values() for v in vals]

    # --- Umap single-read (binary) coverage of the window ---
    n_fully_covered = sum(1 for f in umap_single.values() if f >= 0.999)
    n_zero_covered = sum(1 for f in umap_single.values() if f <= 0.001)
    n_partial_covered = n_total - n_fully_covered - n_zero_covered
    n_no_umap_single_data = n_total - len(umap_single)

    # --- Umap multi-read (continuous) score distribution ---
    mean_score_bins = Counter()
    min_score_bins = Counter()
    for marker_id, (mean_s, min_s, max_s) in umap_multi.items():
        mean_score_bins[bin_score(round(mean_s, 6))] += 1
        min_score_bins[bin_score(round(min_s, 6))] += 1
    n_no_umap_multi_data = n_total - len(umap_multi)

    # --- Pairwise overlaps between the three filter dimensions ---
    low_mappability_ids = {mid for mid, f in umap_single.items() if f < 0.5}
    flank_other_repeat_ids = {
        mid for mid in markers
        if (flank_rmsk.get(mid, set()) - MICROSATELLITE_SELF_CLASSES)
    }
    segdup_ids = set(segdup_hits.keys())

    pairwise = {
        "low_mappability_AND_flank_other_repeat": len(low_mappability_ids & flank_other_repeat_ids),
        "low_mappability_AND_segdup": len(low_mappability_ids & segdup_ids),
        "flank_other_repeat_AND_segdup": len(flank_other_repeat_ids & segdup_ids),
        "all_three": len(low_mappability_ids & flank_other_repeat_ids & segdup_ids),
        "any_of_three": len(low_mappability_ids | flank_other_repeat_ids | segdup_ids),
    }

    # --- Clearly labeled filter scenarios (NOT a final panel; counts only) ---
    scenarios = {
        "raw_candidate_universe": n_total,
        "exclude_segdup_overlap": n_total - len(segdup_ids),
        "exclude_flank_other_repeat_class": n_total - len(flank_other_repeat_ids),
        "exclude_umap_single_read_zero_coverage": n_total - n_zero_covered,
        "exclude_umap_single_read_not_fully_covered": n_total - n_partial_covered - n_zero_covered,
        "exclude_any_of_(segdup_or_flank_other_repeat_or_umap_zero_coverage)": (
            n_total - len(segdup_ids | flank_other_repeat_ids
                          | {mid for mid, f in umap_single.items() if f <= 0.001})
        ),
    }

    return {
        "n_total_tier2_candidates": n_total,
        "repeatmasker_locus_self_classification": dict(self_classification_counts),
        "repeatmasker_locus_self_class_hit_counts": dict(self_class_hit_counts.most_common()),
        "repeatmasker_flank_class_hit_counts": dict(flank_class_hit_counts.most_common()),
        "n_loci_with_any_flank_repeatmasker_hit": n_loci_with_any_flank_repeat,
        "n_loci_with_flank_hit_but_only_microsatellite_class": n_loci_with_flank_only_microsatellite_class,
        "n_loci_with_flank_other_repeat_class": len(flank_other_repeat_ids),
        "segmental_duplication": {
            "n_loci_window_overlapping_segdup": n_loci_in_segdup,
            "n_segdup_fracmatch_values": len(segdup_fracmatch_values),
            "fracmatch_min": min(segdup_fracmatch_values) if segdup_fracmatch_values else None,
            "fracmatch_max": max(segdup_fracmatch_values) if segdup_fracmatch_values else None,
            "fracmatch_mean": (sum(segdup_fracmatch_values) / len(segdup_fracmatch_values)
                                if segdup_fracmatch_values else None),
            "provenance_note": "genomicSuperDups.txt last updated 2014 on UCSC hg38 -- "
                                "not re-analyzed since; treat as a QC flag informed by an "
                                "aging annotation, not an authoritative up-to-date map of "
                                "all GRCh38 segmental duplications",
        },
        "umap_single_read_binary_window_coverage": {
            "n_loci_with_data": len(umap_single),
            "n_loci_no_data": n_no_umap_single_data,
            "n_fully_covered_window": n_fully_covered,
            "n_zero_covered_window": n_zero_covered,
            "n_partially_covered_window": n_partial_covered,
        },
        "umap_multi_read_continuous_score_distribution": {
            "n_loci_with_data": len(umap_multi),
            "n_loci_no_data": n_no_umap_multi_data,
            "mean_score_histogram": dict(mean_score_bins),
            "min_score_histogram": dict(min_score_bins),
        },
        "pairwise_overlaps": pairwise,
        "filter_scenarios_not_a_final_panel": scenarios,
    }


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--markers", required=True)
    ap.add_argument("--locus-self-rmsk", required=True)
    ap.add_argument("--flank-rmsk", required=True)
    ap.add_argument("--segdup-hits", required=True)
    ap.add_argument("--umap-single", required=True)
    ap.add_argument("--umap-multi", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    result = analyze(
        args.markers, args.locus_self_rmsk, args.flank_rmsk,
        args.segdup_hits, args.umap_single, args.umap_multi,
    )
    result["analysis_type"] = "READ_ONLY_FEASIBILITY_ANALYSIS_NOT_A_FILTERED_PANEL"

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(result, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
