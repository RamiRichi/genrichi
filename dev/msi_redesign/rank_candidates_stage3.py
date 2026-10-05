"""
GenRichi MSI Methodology Redesign -- DEVELOPMENT / RUO CANDIDATE RANKING.
NOT PRODUCTION CODE. NOT A CLINICALLY VALIDATED PANEL. No MSI classification.

Stage 3: deterministic, fully traceable RANKING of the Stage-2 strong
mononucleotide candidates (View A 6,074; View B 60,756; shared 6,038; union
60,792). It selects nothing: no panel size is chosen, no locus is removed,
and no clinical-performance claim is made. The ranking exists so that
different panel-size scenarios can be evaluated later.

Immutable inputs (verified before anything runs; the script refuses to
proceed on any mismatch):
  - candidate_selection_stage2.v1.manifest.json and every output hash it
    records (strong_mono_primary, strong_str_comparison, traceability);
  - tier2_annotated.v1.tsv, verified against the hash recorded in the
    250 bp manifest that stage 2 itself pins (read-only; used only to pull
    reference-derived Umap uniqueness fields for the strong loci);
  - msi_markers.promega_tier1.tsv, verified against its own manifest (used
    only to derive the length anchor below).

RANKING HIERARCHY (declared here, applied lexicographically; lower key =
better rank; every key is an attribute already present in validated data or
objectively derived from it; all ordering choices are GenRichi development
heuristics unless noted):
  K1 shared_tier            0 = present in BOTH View A and View B reps
                            (portable across WES and WGS strategies),
                            1 = A-only or B-only. A-only / B-only loci are
                            ranked lower, never discarded.
  K2 promega_span_distance  |span - ANCHOR| where ANCHOR = median span of the
                            five verified Promega mononucleotide markers
                            (derived from msi_markers.promega_tier1.tsv, =25).
                            Ranks WITHIN the preserved 15-40 bp range instead
                            of treating it as flat; 15-20 / 21-30 / 31-40 bp
                            are reported separately. Literature link: the
                            clinically established mononucleotide markers are
                            21-27 bp; the anchoring is a development choice.
  K3 umap_evidence_tier     0 = multi-read mean == 1.0 (assessed, strongest
                            evidence); 1 = NOT_ASSESSED (no multi-read data;
                            never imputed, always identifiable).
  K4 geometry               higher read_measurable_fraction first (existing
                            column = (R-2f-L+1)/(R+L-1) from the 100 bp read
                            + 5 bp flank rule); no new sequencing threshold.
  K5 repeatmasker           rmsk_tier (0 microsatellite_class_only,
                            1 no_rmsk_hit, 2 mixed) then flank_other_repeat
                            flag (False before True) -- the categories
                            already defined at the 250 bp stage; no new
                            biological interpretation.
  K6 reference_uniqueness   objective Umap fields from the existing annotation
                            table: single-read window coverage (higher first),
                            then multi-read window MIN (higher first; missing
                            sorts after present, never imputed). No other
                            uniqueness data exist in the validated inputs
                            (e.g. flank k-mer multiplicity), so none is
                            proxied -- reported as UNAVAILABLE.
  K7 tie-break              sha256(marker_id): deterministic, no positional or
                            chromosomal bias.
Chromosome is NOT a ranking key: no quotas, chrY is not removed; chromosome
distribution is reported for every ranking and every relative prefix.

Two rankings are produced per view: PRIMARY (K1-K7) and a SENSITIVITY
ranking without K1 (rank_no_shared_key), so the effect of preferring shared
loci is visible. Prefix slices are reported as PERCENTAGES of each ranked
list (5/10/25/50/100%), not as fixed counts.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import statistics
import sys
from collections import Counter
from datetime import datetime, timezone

SCRIPT_VERSION = "candidate_ranking_stage3.v1"
PREFIX_FRACTIONS = [0.05, 0.10, 0.25, 0.50, 1.00]
EXPECTED_STRONG_COUNTS = {"viewA": 6074, "viewB": 60756, "shared": 6038}   # from the stage-2 report

RANK_HEADER = [
    "marker_id", "source_view_class", "in_viewA_strong", "in_viewB_strong",
    "chrom", "start", "end", "span", "span_subband", "motif", "ref_repeat_length",
    "umap_multiread_mean", "umap_status", "umap_singleread_window_coverage",
    "umap_multiread_min", "umap_multiread_max", "rmsk_self_category", "rmsk_tier",
    "flank_other_repeat_flag", "read_measurable_fraction", "is_chrY",
    "K1_shared_tier", "K2_promega_span_distance", "K3_umap_evidence_tier",
    "K4_read_measurable_fraction", "K5_rmsk_tier", "K5_flank_flag",
    "K6_singleread_cov", "K6_multiread_min", "K7_tie_hash",
    "rank_viewA", "rank_viewB", "rank_viewA_no_shared_key", "rank_viewB_no_shared_key",
]


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path):
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def span_subband(span):
    if 15 <= span <= 20:
        return "15_20"
    if 21 <= span <= 30:
        return "21_30"
    if 31 <= span <= 40:
        return "31_40"
    return "outside_15_40"


def promega_anchor(spans):
    """Median span of the verified Promega markers (objective anchor)."""
    return statistics.median(spans)


def parse_float(text):
    return float(text) if text not in ("", None) else None


def rank_key(rec, anchor, use_shared=True):
    """Lexicographic ranking key (lower = better). rec keys:
    shared(bool) span umap_status read_measurable_fraction rmsk_tier flank(bool)
    singleread_cov multiread_min marker_id."""
    sc = rec["singleread_cov"]
    mm = rec["multiread_min"]
    return (
        (0 if rec["shared"] else 1) if use_shared else 0,                 # K1
        abs(rec["span"] - anchor),                                        # K2
        0 if rec["umap_status"] == "TIER0_MEAN_1.0" else 1,               # K3
        -rec["read_measurable_fraction"],                                 # K4
        rec["rmsk_tier"], 1 if rec["flank"] else 0,                       # K5
        (0, -sc) if sc is not None else (1, 0.0),                         # K6a
        (0, -mm) if mm is not None else (1, 0.0),                         # K6b
        hashlib.sha256(rec["marker_id"].encode("utf-8")).hexdigest(),     # K7
    )


def rank_records(recs, anchor, use_shared=True):
    """Returns {marker_id: rank (1-based)} for recs; deterministic and unique."""
    ordered = sorted(recs, key=lambda r: rank_key(r, anchor, use_shared))
    return {r["marker_id"]: i + 1 for i, r in enumerate(ordered)}, ordered


def prefix_size(n, fraction):
    return max(1, math.ceil(n * fraction)) if n else 0


def slice_composition(recs):
    return {
        "n": len(recs),
        "chromosome": dict(sorted(Counter(r["chrom"] for r in recs).items())),
        "span_subband": dict(Counter(span_subband(r["span"]) for r in recs)),
        "umap_status": dict(Counter(r["umap_status"] for r in recs)),
        "rmsk_self_category": dict(Counter(r["rmsk_category"] for r in recs)),
        "shared_vs_only": dict(Counter("shared" if r["shared"] else "view_only" for r in recs)),
        "chrY": sum(1 for r in recs if r["chrom"] == "chrY"),
    }


def prefix_report(ordered):
    n = len(ordered)
    return {f"top_{int(f * 100)}pct": slice_composition(ordered[:prefix_size(n, f)]) for f in PREFIX_FRACTIONS}


def build_records(strong_rows, annotated_lookup):
    recs = []
    for r in strong_rows:
        mid = r["marker_id"]
        ann = annotated_lookup[mid]
        recs.append({
            "marker_id": mid, "source_view_class": r["source_view_class"],
            "in_a": r["in_viewA_reps"] == "True", "in_b": r["in_viewB_reps"] == "True",
            "shared": r["source_view_class"] == "A_and_B",
            "chrom": r["chrom"], "start": int(r["start"]), "end": int(r["end"]),
            "span": int(r["span"]), "motif": r["motif"], "ref_repeat_length": r["ref_repeat_length"],
            "umap_mean": r["umap_multiread_mean"], "umap_status": r["umap_status"],
            "rmsk_category": r["rmsk_self_category"], "rmsk_tier": int(r["rmsk_tier"]),
            "flank": r["flank_other_repeat_flag"] == "True",
            "read_measurable_fraction": float(r["read_measurable_fraction"]),
            "is_chrY": r["is_chrY"] == "True",
            "singleread_cov": parse_float(ann["umap_singleread_window_coverage"]),
            "multiread_min": parse_float(ann["umap_multiread_min"]),
            "multiread_max": parse_float(ann["umap_multiread_max"]),
        })
    return recs


def verify_inputs(stage2_manifest_path, promega_manifest_path, promega_tsv, annotated_path,
                  expect_stage2_manifest_sha256=None):
    """Returns (stage2_manifest dict, input_hashes dict). Raises SystemExit on any mismatch."""
    hashes = {}
    with open(stage2_manifest_path, encoding="utf-8") as fh:
        m = json.load(fh)
    hashes["stage2_manifest"] = sha256_file(stage2_manifest_path)
    if expect_stage2_manifest_sha256 and hashes["stage2_manifest"] != expect_stage2_manifest_sha256:
        raise SystemExit("stage-2 manifest hash differs from the expected value")
    for name, info in m["outputs"].items():
        got = sha256_file(info["path"])
        if got != info["sha256"]:
            raise SystemExit(f"stage-2 output {name} does not match the manifest hash")
        hashes[f"stage2_{name}"] = got
    for flag in ("hcc1395_results_used", "classification_performed", "production_changes",
                 "final_panel_size_fixed", "previous_250bp_outputs_modified"):
        if m[flag] is not False:
            raise SystemExit(f"stage-2 manifest safety flag {flag} is not False")
    with open(m["prior_manifest"]["path"], encoding="utf-8") as fh:
        prior = json.load(fh)
    want_annot = prior["inputs"]["annotated"]["sha256"]
    got_annot = sha256_file(annotated_path)
    if got_annot != want_annot:
        raise SystemExit("tier2_annotated.v1.tsv does not match the hash pinned by the 250 bp manifest")
    hashes["tier2_annotated"] = got_annot
    with open(promega_manifest_path, encoding="utf-8") as fh:
        pm = json.load(fh)
    got_tier1 = sha256_file(promega_tsv)
    if got_tier1 != pm["marker_file_sha256"]:
        raise SystemExit("Promega tier-1 file does not match its manifest hash")
    hashes["promega_tier1"] = got_tier1
    return m, hashes


def write_tsv(path, header, rows):
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--stage2-manifest", required=True)
    ap.add_argument("--annotated", required=True)
    ap.add_argument("--promega-tier1", required=True)
    ap.add_argument("--promega-manifest", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--expect-stage2-manifest-sha256", default=None)
    ap.add_argument("--skip-expected-count-check", action="store_true")
    args = ap.parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)

    m2, input_hashes = verify_inputs(args.stage2_manifest, args.promega_manifest, args.promega_tier1,
                                     args.annotated, args.expect_stage2_manifest_sha256)

    strong_rows = read_tsv(m2["outputs"]["strong_mono_primary"]["path"])
    ids = {r["marker_id"] for r in strong_rows}
    wanted = {}
    with open(args.annotated, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["marker_id"] in ids:
                wanted[row["marker_id"]] = row
    if len(wanted) != len(ids):
        raise SystemExit("some strong loci are missing from the annotation table")

    recs = build_records(strong_rows, wanted)
    a_recs = [r for r in recs if r["in_a"]]
    b_recs = [r for r in recs if r["in_b"]]
    shared = [r for r in recs if r["shared"]]
    counts = {"viewA": len(a_recs), "viewB": len(b_recs), "shared": len(shared), "union": len(recs)}
    if not args.skip_expected_count_check and \
            {k: counts[k] for k in EXPECTED_STRONG_COUNTS} != EXPECTED_STRONG_COUNTS:
        raise SystemExit(f"strong-set counts {counts} differ from the validated stage-2 counts")

    tier1_spans = [int(r["end"]) - int(r["start"]) for r in read_tsv(args.promega_tier1)]
    anchor = promega_anchor(tier1_spans)

    ranks = {}
    ordered = {}
    for view, subset in (("viewA", a_recs), ("viewB", b_recs)):
        for label, use_shared in (("primary", True), ("no_shared_key", False)):
            r, o = rank_records(subset, anchor, use_shared)
            ranks[(view, label)] = r
            ordered[(view, label)] = o

    def row_out(r):
        mid = r["marker_id"]
        cov = "" if r["singleread_cov"] is None else r["singleread_cov"]
        mmin = "" if r["multiread_min"] is None else r["multiread_min"]
        mmax = "" if r["multiread_max"] is None else r["multiread_max"]
        return [
            mid, r["source_view_class"], r["in_a"], r["in_b"], r["chrom"], r["start"], r["end"], r["span"],
            span_subband(r["span"]), r["motif"], r["ref_repeat_length"], r["umap_mean"], r["umap_status"],
            cov, mmin, mmax, r["rmsk_category"], r["rmsk_tier"], r["flank"],
            r["read_measurable_fraction"], r["is_chrY"],
            0 if r["shared"] else 1, abs(r["span"] - anchor),
            0 if r["umap_status"] == "TIER0_MEAN_1.0" else 1,
            r["read_measurable_fraction"], r["rmsk_tier"], 1 if r["flank"] else 0,
            cov, mmin, hashlib.sha256(mid.encode("utf-8")).hexdigest()[:16],
            ranks[("viewA", "primary")].get(mid, ""), ranks[("viewB", "primary")].get(mid, ""),
            ranks[("viewA", "no_shared_key")].get(mid, ""), ranks[("viewB", "no_shared_key")].get(mid, ""),
        ]

    union_sorted = sorted(recs, key=lambda r: rank_key(r, anchor, True))
    union_path = os.path.join(args.out_dir, "candidate_ranking_stage3.v1.ranked_union.tsv")
    write_tsv(union_path, RANK_HEADER, (row_out(r) for r in union_sorted))
    view_paths = {}
    for view in ("viewA", "viewB"):
        p = os.path.join(args.out_dir, f"candidate_ranking_stage3.v1.ranked_{view}.tsv")
        write_tsv(p, RANK_HEADER, (row_out(r) for r in ordered[(view, "primary")]))
        view_paths[view] = p

    report = {}
    for view, subset in (("viewA", a_recs), ("viewB", b_recs)):
        prim, sens = ordered[(view, "primary")], ordered[(view, "no_shared_key")]
        n = len(prim)
        overlap = {}
        for f in PREFIX_FRACTIONS:
            k = prefix_size(n, f)
            overlap[f"top_{int(f * 100)}pct"] = len({r["marker_id"] for r in prim[:k]} &
                                                    {r["marker_id"] for r in sens[:k]})
        report[view] = {
            "n_ranked": n,
            "full_set_composition": slice_composition(prim),
            "primary_ranking_prefixes": prefix_report(prim),
            "sensitivity_no_shared_key_prefixes": prefix_report(sens),
            "primary_vs_sensitivity_prefix_overlap_counts": overlap,
        }

    unavailable = {
        "flank_kmer_multiplicity_or_other_uniqueness": "UNAVAILABLE in validated inputs; not proxied",
        "umap_multiread_for_chrY": "UNAVAILABLE (no multi-read bedGraph coverage); kept as NOT_ASSESSED, never imputed",
    }
    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "stage": "deterministic ranking of stage-2 strong candidates (development / RUO)",
        "status": "DEVELOPMENT_RANKING_NOT_CLINICALLY_VALIDATED",
        "final_panel_size_fixed": False, "panel_size_selected": False,
        "script": {"version": SCRIPT_VERSION, "path": os.path.abspath(__file__),
                   "sha256": sha256_file(os.path.abspath(__file__)), "python": platform.python_version()},
        "inputs_verified_hashes": input_hashes,
        "stage2_manifest_path": args.stage2_manifest,
        "counts": counts,
        "parameters": {"promega_tier1_path": args.promega_tier1,
                       "promega_tier1_spans": tier1_spans, "promega_span_anchor_bp": anchor,
                       "prefix_fractions": PREFIX_FRACTIONS,
                       "hierarchy": ["K1 shared_tier", "K2 promega_span_distance", "K3 umap_evidence_tier",
                                     "K4 read_measurable_fraction (desc)", "K5 rmsk_tier then flank flag",
                                     "K6 single-read cov (desc) then multi-read min (desc)", "K7 sha256(marker_id)"],
                       "chromosome_is_a_ranking_key": False, "chromosome_quotas": False},
        "unavailable_information": unavailable,
        "hcc1395_results_used": False, "classification_performed": False,
        "production_changes": False, "stage2_outputs_modified": False,
        "report": report,
        "outputs": {"ranked_union": {"path": union_path, "sha256": sha256_file(union_path)},
                    "ranked_viewA": {"path": view_paths["viewA"], "sha256": sha256_file(view_paths["viewA"])},
                    "ranked_viewB": {"path": view_paths["viewB"], "sha256": sha256_file(view_paths["viewB"])}},
    }
    manifest_path = os.path.join(args.out_dir, "candidate_ranking_stage3.v1.manifest.json")
    with open(manifest_path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print(f"wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
