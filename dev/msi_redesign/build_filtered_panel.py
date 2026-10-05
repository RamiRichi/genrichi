"""
GenRichi MSI Methodology Redesign -- DEVELOPMENT PROTOTYPE. NOT PRODUCTION CODE.

Filtered MSI Development Panel v1 -- PANEL CONSTRUCTION ONLY. Does not run
any classification, does not touch msi_locus_model.py's thresholds, does
not use HCC1395 or any BAM data in any way. Scoped to the Tier 2
(genome-wide, cancer-panel-independent) candidate universe only; the 5
Tier 1 Promega markers already live in their own file
(msi_markers.promega_tier1.tsv) and are untouched by this script.

Exclusion criteria applied (each independently evidence-based; see the
mappability/repeat feasibility report delivered earlier this session for
the full distributions each of these is drawn from -- none is an arbitrary
mid-range cutoff):

  UMAP_UNMAPPABLE:
    single-read window coverage fraction == 0  (no unique 100-mer anywhere
    in the locus+-100bp window), OR
    continuous multi-read MIN score == 0  (at least one position in the
    window has zero probability of unique placement).
    Rationale: both are the natural, non-arbitrary zero/non-zero boundary
    -- not an invented intermediate threshold. A locus with either
    condition is unmappable somewhere in the region a real read would need
    to anchor.

  SEGDUP_OVERLAP:
    locus+-100bp window overlaps any UCSC genomicSuperDups entry.
    Rationale: all observed fracMatch values in our data were 0.90-1.00
    (uniformly high-identity), so presence/absence is the only meaningful
    split -- there is no low-identity population to threshold against.
    Provenance caveat: this UCSC track was last updated 2014 and is used
    as an informative flag, not asserted to be a fully current segdup map.

  REPEATMASKER_SELF_NON_MICROSATELLITE:
    the locus ITSELF overlaps ONLY RepeatMasker classes other than
    Simple_repeat/Low_complexity (the "unexpected_non_microsatellite_class"
    bucket from the feasibility analysis -- dominated by SINE/LINE, i.e.
    very likely retrotransposon-internal poly-A tails, not independent
    unique-context microsatellites).
    Loci classified "mixed_microsatellite_and_other" are KEPT, not
    excluded -- the locus itself does still contain a genuine
    Simple_repeat/Low_complexity component, unlike the purely-other case.

RepeatMasker FLANK annotation is recorded but NEVER used to exclude --
it stays a QC flag (see feasibility report: 74% prevalence makes it
non-discriminating as a hard filter; Umap mappability is the direct
measure of the same underlying concern).

No MSI classification threshold (MIN_LOCUS_DEPTH, INSTABILITY_DISTANCE_
CUTOFF, etc. in msi_locus_model.py) is read, used, or implied here --
this script only selects which reference-genome LOCI are in the panel.

Usage:
    python3 build_filtered_panel.py --tier2-raw dev/msi_redesign/output/msi_markers.v2.tsv \
        --annotations-dir dev/msi_redesign/reference_annotations \
        --out-dir dev/msi_redesign/output
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
from datetime import datetime, timezone

MICROSATELLITE_SELF_CLASSES = {"Simple_repeat", "Low_complexity"}


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_pair_tsv(path: str):
    result = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            result.setdefault(parts[0], set()).add(parts[1])
    return result


def load_segdup_hits(path: str):
    result = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            try:
                result.setdefault(parts[0], []).append(float(parts[1]))
            except ValueError:
                continue
    return result


def load_umap_single(path: str):
    result = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            result[parts[3]] = float(parts[-1])
    return result


def load_umap_multi(path: str):
    result = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7 or parts[4] == ".":
                continue
            result[parts[3]] = (float(parts[4]), float(parts[5]), float(parts[6]))
    return result


def classify_self_rmsk(classes: set) -> str:
    if not classes:
        return "no_rmsk_hit"
    non_self = classes - MICROSATELLITE_SELF_CLASSES
    if not non_self:
        return "microsatellite_class_only"
    if classes & MICROSATELLITE_SELF_CLASSES:
        return "mixed_microsatellite_and_other"
    return "unexpected_non_microsatellite_class"


def determine_exclusion_reasons(
    marker_id: str, locus_self_rmsk: dict, flank_rmsk: dict, segdup: dict,
    umap_single: dict, umap_multi: dict,
) -> list:
    """Returns a list of exclusion reason strings (empty list = kept)."""
    reasons = []

    single_cov = umap_single.get(marker_id)
    multi = umap_multi.get(marker_id)
    multi_min = multi[1] if multi is not None else None
    if (single_cov is not None and single_cov == 0.0) or (multi_min is not None and multi_min == 0.0):
        reasons.append("UMAP_UNMAPPABLE")

    if marker_id in segdup:
        reasons.append("SEGDUP_OVERLAP")

    self_classes = locus_self_rmsk.get(marker_id, set())
    if classify_self_rmsk(self_classes) == "unexpected_non_microsatellite_class":
        reasons.append("REPEATMASKER_SELF_NON_MICROSATELLITE")

    return reasons


def build(tier2_raw_path: str, annotations_dir: str, out_dir: str):
    locus_self_rmsk = load_pair_tsv(os.path.join(annotations_dir, "locus_self_rmsk_classes.tsv"))
    flank_rmsk = load_pair_tsv(os.path.join(annotations_dir, "flank_rmsk_classes.tsv"))
    segdup = load_segdup_hits(os.path.join(annotations_dir, "window_segdup_hits.tsv"))
    umap_single = load_umap_single(os.path.join(annotations_dir, "window_umap_singleread_coverage.tsv"))
    umap_multi = load_umap_multi(os.path.join(annotations_dir, "window_umap_multiread_scores.tsv"))

    annotated_rows = []
    cumulative = {
        "raw_candidate_universe": 0,
        "after_excluding_umap_unmappable": 0,
        "after_excluding_umap_unmappable_and_segdup": 0,
        "after_excluding_umap_unmappable_and_segdup_and_rmsk_self": 0,
    }
    exclusion_reason_counts = {}

    with open(tier2_raw_path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            marker_id = row["marker_id"]
            cumulative["raw_candidate_universe"] += 1

            reasons = determine_exclusion_reasons(
                marker_id, locus_self_rmsk, flank_rmsk, segdup, umap_single, umap_multi
            )
            has_umap = "UMAP_UNMAPPABLE" in reasons
            has_segdup = "SEGDUP_OVERLAP" in reasons
            has_rmsk = "REPEATMASKER_SELF_NON_MICROSATELLITE" in reasons

            if not has_umap:
                cumulative["after_excluding_umap_unmappable"] += 1
            if not has_umap and not has_segdup:
                cumulative["after_excluding_umap_unmappable_and_segdup"] += 1
            if not has_umap and not has_segdup and not has_rmsk:
                cumulative["after_excluding_umap_unmappable_and_segdup_and_rmsk_self"] += 1

            for r in reasons:
                exclusion_reason_counts[r] = exclusion_reason_counts.get(r, 0) + 1
            if not reasons:
                exclusion_reason_counts["KEPT"] = exclusion_reason_counts.get("KEPT", 0) + 1

            multi = umap_multi.get(marker_id)
            segdup_vals = segdup.get(marker_id)
            self_classes = locus_self_rmsk.get(marker_id, set())
            flank_classes = flank_rmsk.get(marker_id, set())
            flank_other = bool(flank_classes - MICROSATELLITE_SELF_CLASSES)

            annotated_rows.append({
                "marker_id": marker_id,
                "chrom": row["chrom"],
                "start": row["start"],
                "end": row["end"],
                "motif": row["motif"],
                "ref_repeat_length": row["ref_repeat_length"],
                "kind": row["kind"],
                "umap_singleread_window_coverage": (
                    "" if marker_id not in umap_single else round(umap_single[marker_id], 6)
                ),
                "umap_multiread_mean": "" if multi is None else round(multi[0], 6),
                "umap_multiread_min": "" if multi is None else round(multi[1], 6),
                "umap_multiread_max": "" if multi is None else round(multi[2], 6),
                "segdup_overlap": bool(segdup_vals),
                "segdup_fracmatch_max": "" if not segdup_vals else round(max(segdup_vals), 6),
                "repeatmasker_self_class": ";".join(sorted(self_classes)) if self_classes else "",
                "repeatmasker_self_category": classify_self_rmsk(self_classes),
                "repeatmasker_flank_other_repeat_flag": flank_other,
                "repeatmasker_flank_classes": ";".join(sorted(flank_classes)) if flank_classes else "",
                "exclusion_reasons": ";".join(reasons),
                "kept": len(reasons) == 0,
            })

    os.makedirs(out_dir, exist_ok=True)

    annotated_path = os.path.join(out_dir, "tier2_annotated.v1.tsv")
    fieldnames = list(annotated_rows[0].keys()) if annotated_rows else []
    with open(annotated_path, "w", encoding="utf-8", newline="\n") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for r in annotated_rows:
            writer.writerow(r)

    filtered_path = os.path.join(out_dir, "tier2_filtered.v1.tsv")
    filtered_fieldnames = ["marker_id", "chrom", "start", "end", "motif", "ref_repeat_length", "kind"]
    with open(filtered_path, "w", encoding="utf-8", newline="\n") as fh:
        writer = csv.DictWriter(fh, fieldnames=filtered_fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for r in annotated_rows:
            if r["kept"]:
                writer.writerow(r)

    return annotated_path, filtered_path, cumulative, exclusion_reason_counts


def audit_filtered_panel(filtered_path: str):
    """Descriptive audit of the surviving filtered panel: motif, chromosome,
    repeat-length distributions -- read-only, no further filtering."""
    from collections import Counter
    kind_counts = Counter()
    motif_counts = Counter()
    chrom_counts = Counter()
    length_hist = Counter()

    with open(filtered_path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            kind_counts[row["kind"]] += 1
            motif_counts[row["motif"]] += 1
            chrom_counts[row["chrom"]] += 1
            length_hist[int(row["ref_repeat_length"])] += 1

    return {
        "n_markers": sum(kind_counts.values()),
        "kind_counts": dict(kind_counts),
        "motif_counts_top30": dict(motif_counts.most_common(30)),
        "chrom_counts": dict(chrom_counts),
        "repeat_length_histogram": dict(sorted(length_hist.items())),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--tier2-raw", required=True)
    ap.add_argument("--annotations-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args(argv)

    tier2_raw_immutable = os.path.join(args.out_dir, "tier2_raw.v1.tsv")
    os.makedirs(args.out_dir, exist_ok=True)
    with open(args.tier2_raw, "rb") as src, open(tier2_raw_immutable, "wb") as dst:
        dst.write(src.read())

    annotated_path, filtered_path, cumulative, exclusion_counts = build(
        tier2_raw_immutable, args.annotations_dir, args.out_dir
    )
    audit = audit_filtered_panel(filtered_path)

    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "panel_version": "tier2_filtered.v1",
        "scope": "Tier 2 (genome-wide, cancer-panel-independent) only -- Tier 1 Promega markers untouched, live separately in msi_markers.promega_tier1.tsv",
        "hcc1395_used_in_selection": False,
        "classification_thresholds_modified": False,
        "production_files_modified": False,
        "exclusion_criteria": {
            "UMAP_UNMAPPABLE": "single-read window coverage == 0 OR multi-read window min score == 0",
            "SEGDUP_OVERLAP": "locus+-100bp window overlaps any UCSC genomicSuperDups entry (track dated 2014, used as a flag, not asserted current)",
            "REPEATMASKER_SELF_NON_MICROSATELLITE": "locus itself overlaps RepeatMasker classes that are ALL non-Simple_repeat/Low_complexity (dominated by SINE/LINE -- likely retrotransposon poly-A tails)",
            "REPEATMASKER_FLANK": "recorded as a QC flag column only; NEVER used to exclude (74% prevalence makes it non-discriminating; Umap mappability is the direct measure)",
        },
        "cumulative_counts": cumulative,
        "exclusion_reason_counts_non_exclusive": exclusion_counts,
        "n_raw": cumulative["raw_candidate_universe"],
        "n_filtered_kept": cumulative["after_excluding_umap_unmappable_and_segdup_and_rmsk_self"],
        "files": {
            "tier2_raw_v1_tsv": {"path": tier2_raw_immutable, "sha256": sha256_file(tier2_raw_immutable)},
            "tier2_annotated_v1_tsv": {"path": annotated_path, "sha256": sha256_file(annotated_path)},
            "tier2_filtered_v1_tsv": {"path": filtered_path, "sha256": sha256_file(filtered_path)},
        },
        "filtered_panel_audit": audit,
    }
    manifest_path = os.path.join(args.out_dir, "tier2_filtered.v1.manifest.json")
    with open(manifest_path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
        fh.write("\n")

    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
