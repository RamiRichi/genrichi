"""
GenRichi MSI Methodology Redesign -- DEVELOPMENT PROTOTYPE. NOT PRODUCTION CODE.

Development-only "Candidate MSI Panel v1" spacing scenario. NOT a
clinically validated panel, NOT a final panel, and it has no fixed panel
size: the number of representatives is simply the outcome of applying the
pre-declared rules below to each view.

Spacing (default 250 bp) is a GenRichi DEVELOPMENT HEURISTIC informed by the
real-BAM sequencing-geometry analysis (sequencing_geometry_analysis.v1.json:
same-read redundancy ends near 75 bp; fragment-level redundancy is roughly
5-7% near 250 bp). It is not a clinically validated threshold. A 75 bp
scenario is computed only as a documented sensitivity comparison.

Inputs are read-only: tier2_annotated.v1.tsv (per-locus annotations, all
already computed) plus View A membership. tier2_raw / tier2_filtered /
View A / View B are never modified. No HCC1395 variant/MSI result, no BAM,
and no classification threshold is used or touched here.

Procedure, per view (View A and View B are processed independently):
  1. CLUSTERS: single-linkage over all loci on a chromosome; a locus joins
     the current cluster if gap = start - (max end so far) <= spacing.
     Every locus gets a cluster_id (traceability).
  2. ELIGIBILITY: loci whose span exceeds READ_LENGTH - 2*MIN_FLANK (=90 bp)
     cannot be measured by a single 100 bp read under bam_extract.py's
     min_flank rule; they are labelled INELIGIBLE_EXTRACTION_SPAN and are
     never selected (they neither become representatives nor exclude anyone).
  3. PRIORITY (lexicographic, lower = better; all GenRichi development
     choices, listed in PRIORITY_CRITERIA below):
       umap_tier, rmsk_tier, kind_tier, length_tier, extraction_span_tier,
       flank_repeat_flag, -span, sha256(marker_id) (deterministic tie-break,
       no positional/chromosomal bias).
  4. SELECTION: greedy maximal independent set -- take the best remaining
     eligible locus, then mark every eligible locus within `spacing` bp of
     it (gap <= spacing) EXCLUDED_BY_SPACING with representative_id set.
     Result: all selected loci are > spacing apart and every excluded locus
     is within `spacing` of the representative that excluded it.
Chromosome/genomic distribution is monitored and reported (selected vs
pre-selection composition), not enforced by quotas.

Usage:
    python3 build_candidate_panel_spacing.py --annotated tier2_annotated.v1.tsv \
        --view-a view_A_wes_compatible.v1.tsv --out-dir OUT [--spacing 250]
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
from collections import Counter, defaultdict, namedtuple
from datetime import datetime, timezone

PRIMARY_SPACING_BP = 250
SENSITIVITY_SPACING_BP = 75
READ_LENGTH_BP = 100          # confirmed from the real BAMs (geometry analysis)
MIN_FLANK_BP = 5              # bam_extract.py MIN_FLANK
MAX_MEASURABLE_SPAN_BP = READ_LENGTH_BP - 2 * MIN_FLANK_BP

Locus = namedtuple(
    "Locus",
    "marker_id chrom start end kind motif rep_len umap_mean rmsk_cat flank_flag",
)

PRIORITY_CRITERIA = {
    "1_umap_tier": "0: Umap multi-read mean == 1.0; 1: >=0.9; 2: >=0.7; 3: <0.7; 4: no Umap data (e.g. chrY)",
    "2_rmsk_tier": "0: locus-self RepeatMasker = Simple_repeat/Low_complexity only (independently recognised); "
                   "1: no RepeatMasker hit; 2: mixed microsatellite+other class",
    "3_kind_tier": "0: mononucleotide; 1: STR (literature: mononucleotide loci are more sensitive/specific for MSI)",
    "4_length_tier": "span 15-40 bp: 0; 10-14 or 41-60 bp: 1; otherwise: 2 (development heuristic; validated Promega "
                     "mononucleotide markers are 21-27 bp, very long runs carry more sequencing noise)",
    "5_extraction_span_tier": "0: span <= 60 bp; 1: 61-90 bp (fraction of fragments able to measure it falls with span)",
    "6_flank_repeat_flag": "0: no other-class repeat in the 100 bp flanks; 1: present (QC flag used only as a tie-break)",
    "7_longer_span_first": "within all of the above, longer span preferred",
    "8_tie_break": "sha256(marker_id): deterministic, no positional or chromosomal bias",
}


def span_of(locus) -> int:
    return locus.end - locus.start


def umap_tier(mean) -> int:
    if mean is None:
        return 4
    if mean >= 1.0:
        return 0
    if mean >= 0.9:
        return 1
    if mean >= 0.7:
        return 2
    return 3


def rmsk_tier(category: str) -> int:
    if category == "microsatellite_class_only":
        return 0
    if category == "no_rmsk_hit":
        return 1
    return 2  # mixed_microsatellite_and_other (unexpected_non_microsatellite_class is already excluded upstream)


def kind_tier(kind: str) -> int:
    return 0 if kind == "mononucleotide" else 1


def length_tier(span: int) -> int:
    if 15 <= span <= 40:
        return 0
    if 10 <= span <= 14 or 41 <= span <= 60:
        return 1
    return 2


def extraction_span_tier(span: int) -> int:
    return 0 if span <= 60 else 1


def is_eligible(locus) -> bool:
    return span_of(locus) <= MAX_MEASURABLE_SPAN_BP


def priority_key(locus):
    sp = span_of(locus)
    digest = hashlib.sha256(locus.marker_id.encode("utf-8")).hexdigest()
    return (
        umap_tier(locus.umap_mean), rmsk_tier(locus.rmsk_cat), kind_tier(locus.kind),
        length_tier(sp), extraction_span_tier(sp), 1 if locus.flank_flag else 0,
        -sp, digest,
    )


def assign_clusters(sorted_loci, spacing):
    """sorted_loci: list of Locus on ONE chromosome sorted by (start, end).
    Returns a list of cluster indices (0-based, per chromosome)."""
    ids = []
    cid = -1
    max_end = None
    for locus in sorted_loci:
        if max_end is None or locus.start - max_end > spacing:
            cid += 1
            max_end = locus.end
        else:
            max_end = max(max_end, locus.end)
        ids.append(cid)
    return ids


def select_representatives(sorted_loci, spacing):
    """Returns (status, representative_index) lists aligned with sorted_loci.
    status in {SELECTED, EXCLUDED_BY_SPACING, INELIGIBLE_EXTRACTION_SPAN}."""
    n = len(sorted_loci)
    status = [None] * n
    rep = [None] * n
    for i, locus in enumerate(sorted_loci):
        if not is_eligible(locus):
            status[i] = "INELIGIBLE_EXTRACTION_SPAN"
    max_len = max((span_of(l) for l in sorted_loci), default=0)
    order = sorted((i for i in range(n) if status[i] is None),
                   key=lambda i: priority_key(sorted_loci[i]))
    for i in order:
        if status[i] is not None:
            continue
        status[i] = "SELECTED"
        li = sorted_loci[i]
        j = i + 1
        while j < n and sorted_loci[j].start - li.end <= spacing:
            if status[j] is None:
                status[j] = "EXCLUDED_BY_SPACING"
                rep[j] = i
            j += 1
        j = i - 1
        while j >= 0 and sorted_loci[j].start >= li.start - spacing - max_len:
            if status[j] is None and li.start - sorted_loci[j].end <= spacing:
                status[j] = "EXCLUDED_BY_SPACING"
                rep[j] = i
            j -= 1
    return status, rep


def process_view(loci, spacing):
    """loci: iterable of Locus. Returns per-locus records (dicts) plus
    summary counts. Nothing is removed from the input."""
    by_chrom = defaultdict(list)
    for l in loci:
        by_chrom[l.chrom].append(l)
    records = []
    n_clusters = 0
    n_multi_clusters = 0
    n_clusters_without_rep = 0
    for chrom in sorted(by_chrom):
        chrom_loci = sorted(by_chrom[chrom], key=lambda l: (l.start, l.end))
        cids = assign_clusters(chrom_loci, spacing)
        status, rep = select_representatives(chrom_loci, spacing)
        sizes = Counter(cids)
        sel_per_cluster = Counter(c for c, s in zip(cids, status) if s == "SELECTED")
        n_clusters += len(sizes)
        n_multi_clusters += sum(1 for c, k in sizes.items() if k >= 2)
        n_clusters_without_rep += sum(1 for c in sizes if sel_per_cluster[c] == 0)
        for locus, cid, st, r in zip(chrom_loci, cids, status, rep):
            rep_locus = chrom_loci[r] if r is not None else None
            records.append({
                "locus": locus, "cluster_id": f"{chrom}:C{cid}", "cluster_size": sizes[cid],
                "status": st,
                "representative_id": rep_locus.marker_id if rep_locus else "",
                "gap_to_representative_bp": (
                    (max(locus.start - rep_locus.end, rep_locus.start - locus.end)) if rep_locus else ""),
            })
    counts = Counter(r["status"] for r in records)
    return records, {
        "n_loci_in_view": len(records),
        "n_clusters_total_including_singletons": n_clusters,
        "n_clusters_with_2plus_loci": n_multi_clusters,
        "n_clusters_without_any_representative": n_clusters_without_rep,
        "n_selected_representatives": counts["SELECTED"],
        "n_excluded_by_spacing": counts["EXCLUDED_BY_SPACING"],
        "n_ineligible_extraction_span": counts["INELIGIBLE_EXTRACTION_SPAN"],
    }


def composition(loci):
    loci = list(loci)
    band = Counter()
    for l in loci:
        s = span_of(l)
        band["<10bp" if s < 10 else "10-14bp" if s <= 14 else "15-40bp" if s <= 40
             else "41-60bp" if s <= 60 else "61-90bp" if s <= 90 else ">90bp"] += 1
    return {
        "n": len(loci),
        "kind": dict(Counter(l.kind for l in loci)),
        "motif_length": dict(sorted(Counter(len(l.motif) for l in loci).items())),
        "motif_top15": dict(Counter(l.motif for l in loci).most_common(15)),
        "span_band": dict(band),
        "chromosome": dict(Counter(l.chrom for l in loci)),
        "umap_tier": dict(sorted(Counter(umap_tier(l.umap_mean) for l in loci).items())),
        "rmsk_self_category": dict(Counter(l.rmsk_cat for l in loci)),
        "length_tier": dict(sorted(Counter(length_tier(span_of(l)) for l in loci).items())),
        "flank_other_repeat_flag": dict(Counter(bool(l.flank_flag) for l in loci)),
    }


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_loci(annotated_path, view_a_ids_path):
    view_a_ids = set()
    with open(view_a_ids_path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            view_a_ids.add(row["marker_id"])
    view_b, view_a = [], []
    with open(annotated_path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["kept"] != "True":
                continue
            mean = row["umap_multiread_mean"]
            locus = Locus(
                row["marker_id"], row["chrom"], int(row["start"]), int(row["end"]),
                row["kind"], row["motif"], int(row["ref_repeat_length"]),
                float(mean) if mean != "" else None,
                row["repeatmasker_self_category"],
                row["repeatmasker_flank_other_repeat_flag"] == "True",
            )
            view_b.append(locus)
            if locus.marker_id in view_a_ids:
                view_a.append(locus)
    return view_a, view_b


TRACE_HEADER = ["marker_id", "chrom", "start", "end", "kind", "motif", "ref_repeat_length", "span",
                "umap_multiread_mean", "rmsk_self_category", "flank_other_repeat_flag",
                "umap_tier", "rmsk_tier", "kind_tier", "length_tier", "extraction_span_tier",
                "cluster_id", "cluster_size", "status", "representative_id", "gap_to_representative_bp"]


def _row(rec):
    l = rec["locus"]
    sp = span_of(l)
    return [l.marker_id, l.chrom, l.start, l.end, l.kind, l.motif, l.rep_len, sp,
            "" if l.umap_mean is None else l.umap_mean, l.rmsk_cat, l.flank_flag,
            umap_tier(l.umap_mean), rmsk_tier(l.rmsk_cat), kind_tier(l.kind), length_tier(sp),
            extraction_span_tier(sp), rec["cluster_id"], rec["cluster_size"], rec["status"],
            rec["representative_id"], rec["gap_to_representative_bp"]]


def write_tsv(path, header, rows):
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--annotated", required=True)
    ap.add_argument("--view-a", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--spacing", type=int, default=PRIMARY_SPACING_BP)
    ap.add_argument("--sensitivity-spacing", type=int, default=SENSITIVITY_SPACING_BP)
    ap.add_argument("--geometry-json", default=None)
    args = ap.parse_args(argv)

    os.makedirs(args.out_dir, exist_ok=True)
    view_a, view_b = load_loci(args.annotated, args.view_a)
    views = {"viewA_wes_compatible": view_a, "viewB_wgs": view_b}

    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "panel_name": "Candidate MSI Panel v1 (development candidate)",
        "status": "DEVELOPMENT_CANDIDATE_NOT_CLINICALLY_VALIDATED",
        "final_panel_size_fixed": False,
        "primary_spacing_bp": args.spacing,
        "spacing_note": "GenRichi development heuristic informed by real-BAM sequencing geometry; not a validated threshold",
        "sensitivity_spacing_bp": args.sensitivity_spacing,
        "read_length_bp": READ_LENGTH_BP, "min_flank_bp": MIN_FLANK_BP,
        "max_measurable_span_bp": MAX_MEASURABLE_SPAN_BP,
        "priority_criteria_all_development_choices": PRIORITY_CRITERIA,
        "inputs": {"annotated": {"path": args.annotated, "sha256": sha256_file(args.annotated)},
                   "view_A_ids": {"path": args.view_a, "sha256": sha256_file(args.view_a)}},
        "hcc1395_variant_or_msi_results_used": False,
        "classification_performed": False,
        "master_panel_or_views_modified": False,
        "production_changes": False,
        "views": {},
        "outputs": {},
    }
    if args.geometry_json:
        manifest["inputs"]["sequencing_geometry_analysis"] = {
            "path": args.geometry_json, "sha256": sha256_file(args.geometry_json)}

    selected_sets = {}
    for name, loci in views.items():
        records, summary = process_view(loci, args.spacing)
        selected = [r["locus"] for r in records if r["status"] == "SELECTED"]
        selected_sets[name] = {l.marker_id for l in selected}
        sens_records, sens_summary = process_view(loci, args.sensitivity_spacing)
        sens_selected = {r["locus"].marker_id for r in sens_records if r["status"] == "SELECTED"}

        panel_path = os.path.join(args.out_dir, f"candidate_msi_panel.v1.{name}.spacing{args.spacing}.tsv")
        trace_path = os.path.join(args.out_dir, f"candidate_msi_panel.v1.{name}.spacing{args.spacing}.traceability.tsv")
        write_tsv(panel_path, TRACE_HEADER, (_row(r) for r in records if r["status"] == "SELECTED"))
        write_tsv(trace_path, TRACE_HEADER, (_row(r) for r in records))

        manifest["views"][name] = {
            "spacing_summary": summary,
            "pre_selection_composition": composition(loci),
            "selected_composition": composition(selected),
            "sensitivity_75bp_documented_only": {
                **sens_summary,
                "n_selected_in_both_scenarios": len(selected_sets[name] & sens_selected),
                "n_selected_only_at_primary_spacing": len(selected_sets[name] - sens_selected),
                "n_selected_only_at_sensitivity_spacing": len(sens_selected - selected_sets[name]),
            },
        }
        manifest["outputs"][name] = {
            "selected_panel": {"path": panel_path, "sha256": sha256_file(panel_path)},
            "traceability_all_loci": {"path": trace_path, "sha256": sha256_file(trace_path)},
        }

    a, b = selected_sets["viewA_wes_compatible"], selected_sets["viewB_wgs"]
    view_a_ids = {l.marker_id for l in view_a}
    manifest["view_A_vs_view_B"] = {
        "n_viewA_representatives": len(a), "n_viewB_representatives": len(b),
        "n_representatives_in_both": len(a & b),
        "n_viewA_reps_not_selected_in_viewB": len(a - b),
        "n_viewB_reps_inside_viewA_footprint_set": len(b & view_a_ids),
        "n_viewB_reps_outside_viewA": len(b - view_a_ids),
    }

    manifest_path = os.path.join(args.out_dir, f"candidate_msi_panel.v1.spacing{args.spacing}.manifest.json")
    with open(manifest_path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print(f"wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
