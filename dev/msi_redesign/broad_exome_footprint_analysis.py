"""
GenRichi MSI Methodology Redesign -- DEVELOPMENT PROTOTYPE. NOT PRODUCTION CODE.

READ-ONLY "Broad exome footprint compatibility" analysis of the filtered
Tier 2 development panel (tier2_filtered.v1.tsv, 1,208,169 loci).

Deliberately named "Broad exome footprint compatibility", NOT "clinical
capture compatibility" -- this checks overlap against the Broad Institute's
generic public GATK exome_calling_regions.v1 interval list (a broad,
GATK-ecosystem-standard exonic footprint definition), not against any
specific vendor's clinical hybridization-capture bait design. A locus
falling inside this footprint is NOT thereby shown to be capturable by any
particular real-world exome kit -- vendors' actual bait designs differ
(coverage padding, probe density, GC-content tolerance, etc.), and this
script makes no claim about that. It only reports overlap with a generic,
widely-used exonic coordinate reference.

Also performs a DOCUMENTATION-ONLY overlap against gnomAD's disease-
associated STR catalog (broadinstitute/str-analysis, GRCh38) -- this is
explicitly NOT used as a population-variability filter (see prior
decision: population-variability filtering was skipped this round; the
BCGSC genome-wide tandem-repeat catalog was not used due to its
CC BY-NC-ND 4.0 license being unsuitable for a potentially commercial/
clinical project without legal clearance).

This script does not modify tier2_raw.v1.tsv or tier2_filtered.v1.tsv, and
does not select a final footprint-overlap threshold -- it only reports
scenario counts for you to decide among.

Coordinate conversion note: the Broad interval_list uses Picard's 1-based
inclusive convention; converted to 0-based half-open (BED) as
bed_start = interval_start - 1, bed_end = interval_end (unchanged), matching
our own marker file's convention.

Reference-identity note: chromosome LENGTHS in our production hg38.fa match
the Broad reference's @SQ headers exactly for every one of the 455 shared
contigs (ruling out any indel-scale coordinate drift). However M5 (MD5)
checksums differ for chr1, chr17, and chrX between the two (chrM matches
exactly) -- the two files are not byte-identical copies of "GRCh38",
likely reflecting a benign, common divergence between GRCh38 distributions
(e.g. differing patch/masking treatment at a fixed set of positions, which
cannot shift coordinates given the identical lengths). Disclosed, not
silently assumed away.

Usage:
    python3 broad_exome_footprint_analysis.py \
        --filtered-panel dev/msi_redesign/output/tier2_filtered.v1.tsv \
        --coverage-tsv dev/msi_redesign/reference_annotations/filtered_vs_footprint_coverage.tsv \
        --gnomad-overlap-tsv dev/msi_redesign/reference_annotations/filtered_vs_gnomad_disease_overlap.tsv \
        --out dev/msi_redesign/output/broad_exome_footprint_feasibility.json
"""

from __future__ import annotations

import csv
import json
import sys
from collections import Counter


def load_filtered_panel(path: str):
    markers = {}
    with open(path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            markers[row["marker_id"]] = {
                "chrom": row["chrom"], "motif": row["motif"],
                "ref_repeat_length": int(row["ref_repeat_length"]), "kind": row["kind"],
            }
    return markers


def classify_coverage(fraction: float) -> str:
    if fraction >= 0.999999:
        return "fully_inside_footprint"
    if fraction <= 0.0:
        return "outside_footprint"
    return "partially_overlapping_footprint"


def load_coverage(path: str):
    """bedtools coverage output: [chrom,start,end,marker_id, n_b_overlaps,
    n_bases_covered, length, fraction_covered]."""
    result = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            result[parts[3]] = float(parts[-1])
    return result


def load_gnomad_overlap_ids(path: str):
    ids = set()
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4:
                ids.add(parts[3])
    return ids


def breakdown(marker_ids, markers):
    kind_counts = Counter()
    chrom_counts = Counter()
    motif_counts = Counter()
    length_hist = Counter()
    for mid in marker_ids:
        m = markers.get(mid)
        if m is None:
            continue
        kind_counts[m["kind"]] += 1
        chrom_counts[m["chrom"]] += 1
        motif_counts[m["motif"]] += 1
        length_hist[m["ref_repeat_length"]] += 1
    return {
        "n": len(marker_ids),
        "kind_counts": dict(kind_counts),
        "chrom_counts": dict(chrom_counts),
        "motif_counts_top20": dict(motif_counts.most_common(20)),
        "repeat_length_histogram": dict(sorted(length_hist.items())),
    }


def analyze(filtered_panel_path: str, coverage_path: str, gnomad_overlap_path: str):
    markers = load_filtered_panel(filtered_panel_path)
    coverage = load_coverage(coverage_path)
    gnomad_overlap_ids = load_gnomad_overlap_ids(gnomad_overlap_path)

    categories = {"fully_inside_footprint": set(), "partially_overlapping_footprint": set(),
                  "outside_footprint": set()}
    for marker_id in markers:
        fraction = coverage.get(marker_id, 0.0)
        categories[classify_coverage(fraction)].add(marker_id)

    scenarios = {
        "current_filtered_panel_no_footprint_filter": len(markers),
        "require_fully_inside_footprint": len(categories["fully_inside_footprint"]),
        "require_any_overlap_with_footprint": (
            len(categories["fully_inside_footprint"]) + len(categories["partially_overlapping_footprint"])
        ),
        "require_outside_footprint_only": len(categories["outside_footprint"]),
    }

    return {
        "analysis_type": "READ_ONLY_BROAD_EXOME_FOOTPRINT_COMPATIBILITY_ANALYSIS_NOT_A_FINAL_PANEL",
        "terminology_note": "This reports 'Broad exome footprint compatibility' (overlap with a generic "
                             "GATK-ecosystem exonic coordinate reference), NOT 'clinical capture "
                             "compatibility' with any specific vendor's hybridization-capture kit.",
        "n_filtered_panel_input": len(markers),
        "n_gnomad_disease_str_catalog_overlap_documentation_only": len(gnomad_overlap_ids & set(markers.keys())),
        "gnomad_overlap_used_as_filter": False,
        "breakdown_fully_inside_footprint": breakdown(categories["fully_inside_footprint"], markers),
        "breakdown_partially_overlapping_footprint": breakdown(categories["partially_overlapping_footprint"], markers),
        "breakdown_outside_footprint": breakdown(categories["outside_footprint"], markers),
        "filter_scenarios_not_a_final_panel": scenarios,
        "raw_and_filtered_universes_modified": False,
    }


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--filtered-panel", required=True)
    ap.add_argument("--coverage-tsv", required=True)
    ap.add_argument("--gnomad-overlap-tsv", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    result = analyze(args.filtered_panel, args.coverage_tsv, args.gnomad_overlap_tsv)
    with open(args.out, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(result, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
