"""
GenRichi MSI Methodology Redesign -- DEVELOPMENT / RUO SCENARIO PREVIEW.
NOT PRODUCTION CODE. NOT A CLINICALLY VALIDATED PANEL. No MSI classification.

Stage 4: scenario previews (top-N slices) over the Stage-3 rankings. It does
NOT rewrite or re-run any ranking algorithm and does NOT select a panel or a
panel size; it reads the Stage-3 rank columns (immutable, hash-verified) and
reports what each candidate size WOULD contain, so the two policies can be
compared before any decision.

Two ranking concepts, kept SEPARATE and never combined into one score.
Naming note: Stage 3 called the shared-first ranking "primary" and the
shared-free ranking "no_shared_key"; from this stage on the labels are:

  PRIMARY QUALITY RANKING  = Stage-3 columns rank_viewA_no_shared_key /
      rank_viewB_no_shared_key. Shared status does NOT influence rank.
      Hierarchy (documented Stage-3 keys K2-K7 with K1 removed):
        1. |span - 25| (25 = median span of the five verified Promega markers)
        2. Umap multi-read mean == 1.0 before NOT_ASSESSED (never imputed)
        3. higher read_measurable_fraction (100 bp read, 5 bp flank geometry)
        4. RepeatMasker tier, then flank flag (categories already defined)
        5. objective Umap uniqueness fields: single-read coverage (higher
           first), multi-read min (higher first, missing after present)
        6. sha256(marker_id) deterministic tie-break
  SECONDARY COMMON-CORE RANKING = Stage-3 columns rank_viewA / rank_viewB
      (the existing Shared-first ranking: K1 shared tier, then K2-K7).
      Kept as a sensitivity/comparison for a future common-core strategy.

Chromosome is not a ranking criterion; no quotas; chrY is never removed.
Shared / A-only / B-only status is reported in every scenario.
Scenario sizes are previews only -- none is chosen as a final panel. No
membership files are written (membership is simply rank <= N in the Stage-3
union table), so the previews cannot be mistaken for selected panels.

No HCC1395 result, no BAM, no MSI classification, no production file.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import sys
from collections import Counter
from datetime import datetime, timezone

SCRIPT_VERSION = "scenario_preview_stage4.v1"
DEFAULT_SIZES = [500, 1000, 1500, 2000, 5000]
VIEWS = {"viewA": ("in_viewA_strong", "rank_viewA_no_shared_key", "rank_viewA"),
         "viewB": ("in_viewB_strong", "rank_viewB_no_shared_key", "rank_viewB")}
RANKING_LABELS = {"primary_quality_ranking": 1, "secondary_common_core_ranking": 2}   # index into VIEWS tuple
SUBBANDS = ["15_20", "21_30", "31_40"]


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path):
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def top_n(rows, rank_col, n):
    """Rows of the view sorted by rank_col, first n. Rows lacking the rank
    (not in this view) are ignored."""
    ranked = sorted((r for r in rows if r[rank_col] != ""), key=lambda r: int(r[rank_col]))
    return ranked[:n], len(ranked)


def scenario_stats(top_rows, n_available):
    n = len(top_rows)
    classes = Counter(r["source_view_class"] for r in top_rows)
    return {
        "n": n,
        "pct_of_available_strong_set": round(100.0 * n / n_available, 3) if n_available else None,
        "shared": classes.get("A_and_B", 0),
        "A_only": classes.get("A_only", 0),
        "B_only": classes.get("B_only", 0),
        "span_subband": {b: sum(1 for r in top_rows if r["span_subband"] == b) for b in SUBBANDS},
        "umap_status": dict(Counter(r["umap_status"] for r in top_rows)),
        "rmsk_self_category": dict(Counter(r["rmsk_self_category"] for r in top_rows)),
        "flank_other_repeat_flag": dict(Counter(r["flank_other_repeat_flag"] for r in top_rows)),
        "chromosome": dict(sorted(Counter(r["chrom"] for r in top_rows).items())),
        "chrY": sum(1 for r in top_rows if r["chrom"] == "chrY"),
    }


def verify_chain(stage3_manifest_path, expect_sha256=None):
    """Refuses (SystemExit) unless the Stage-3 manifest, its outputs, the
    Stage-2 manifest + outputs and the 250 bp outputs all still match their
    recorded hashes and every safety flag is False. Returns (manifest, hashes)."""
    hashes = {}
    with open(stage3_manifest_path, encoding="utf-8") as fh:
        m3 = json.load(fh)
    hashes["stage3_manifest"] = sha256_file(stage3_manifest_path)
    if expect_sha256 and hashes["stage3_manifest"] != expect_sha256:
        raise SystemExit("stage-3 manifest hash differs from the expected value")
    for flag in ("hcc1395_results_used", "classification_performed", "production_changes",
                 "final_panel_size_fixed", "panel_size_selected", "stage2_outputs_modified"):
        if m3[flag] is not False:
            raise SystemExit(f"stage-3 safety flag {flag} is not False")
    for name, info in m3["outputs"].items():
        got = sha256_file(info["path"])
        if got != info["sha256"]:
            raise SystemExit(f"stage-3 output {name} does not match its recorded hash")
        hashes[f"stage3_{name}"] = got
    s2_path = m3["stage2_manifest_path"]
    got = sha256_file(s2_path)
    if got != m3["inputs_verified_hashes"]["stage2_manifest"]:
        raise SystemExit("stage-2 manifest changed since stage 3")
    hashes["stage2_manifest"] = got
    with open(s2_path, encoding="utf-8") as fh:
        m2 = json.load(fh)
    for name, info in m2["outputs"].items():
        got = sha256_file(info["path"])
        if got != info["sha256"]:
            raise SystemExit(f"stage-2 output {name} changed")
        hashes[f"stage2_{name}"] = got
    prior_path = m2["prior_manifest"]["path"]
    got = sha256_file(prior_path)
    if got != m2["prior_manifest"]["sha256"]:
        raise SystemExit("250 bp manifest changed")
    hashes["spacing250_manifest"] = got
    with open(prior_path, encoding="utf-8") as fh:
        prior = json.load(fh)
    for view, outs in prior["outputs"].items():
        for key, info in outs.items():
            got = sha256_file(info["path"])
            if got != info["sha256"]:
                raise SystemExit(f"250 bp output changed: {view}/{key}")
            hashes[f"spacing250_{view}_{key}"] = got
    return m3, hashes


def build_previews(union_rows, sizes):
    previews = {}
    for view, (in_col, primary_col, shared_first_col) in VIEWS.items():
        view_rows = [r for r in union_rows if r[in_col] == "True"]
        previews[view] = {}
        for size in sizes:
            entry = {}
            sets = {}
            for label, col in (("primary_quality_ranking", primary_col),
                               ("secondary_common_core_ranking", shared_first_col)):
                top, available = top_n(view_rows, col, size)
                entry[label] = scenario_stats(top, available)
                entry["n_available_strong_set"] = available
                sets[label] = {r["marker_id"] for r in top}
            entry["loci_in_both_rankings_top_n"] = len(sets["primary_quality_ranking"] &
                                                        sets["secondary_common_core_ranking"])
            previews[view][str(size)] = entry
    return previews


def write_comparison_tsv(path, previews):
    header = ["view", "size", "ranking", "n", "pct_of_available", "shared", "A_only", "B_only",
              "sub_15_20", "sub_21_30", "sub_31_40", "umap_TIER0", "umap_NOT_ASSESSED",
              "rmsk_microsatellite_class_only", "rmsk_no_rmsk_hit", "rmsk_mixed", "flank_flag_True",
              "chrY", "loci_in_both_rankings"]
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for view, sizes in previews.items():
            for size, entry in sizes.items():
                for label in ("primary_quality_ranking", "secondary_common_core_ranking"):
                    s = entry[label]
                    w.writerow([view, size, label, s["n"], s["pct_of_available_strong_set"], s["shared"],
                                s["A_only"], s["B_only"], s["span_subband"]["15_20"], s["span_subband"]["21_30"],
                                s["span_subband"]["31_40"], s["umap_status"].get("TIER0_MEAN_1.0", 0),
                                s["umap_status"].get("NOT_ASSESSED", 0),
                                s["rmsk_self_category"].get("microsatellite_class_only", 0),
                                s["rmsk_self_category"].get("no_rmsk_hit", 0),
                                s["rmsk_self_category"].get("mixed_microsatellite_and_other", 0),
                                s["flank_other_repeat_flag"].get("True", 0), s["chrY"],
                                entry["loci_in_both_rankings_top_n"]])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--stage3-manifest", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--expect-stage3-manifest-sha256", default=None)
    ap.add_argument("--sizes", type=int, nargs="*", default=DEFAULT_SIZES)
    args = ap.parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)

    m3, hashes = verify_chain(args.stage3_manifest, args.expect_stage3_manifest_sha256)
    union_rows = read_tsv(m3["outputs"]["ranked_union"]["path"])
    previews = build_previews(union_rows, args.sizes)

    comp_path = os.path.join(args.out_dir, "scenario_preview_stage4.v1.comparison.tsv")
    write_comparison_tsv(comp_path, previews)

    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "stage": "scenario previews over the stage-3 rankings (development / RUO)",
        "status": "DEVELOPMENT_SCENARIO_PREVIEW_NOT_A_SELECTED_PANEL",
        "final_panel_size_fixed": False, "panel_size_selected": False,
        "script": {"version": SCRIPT_VERSION, "path": os.path.abspath(__file__),
                   "sha256": sha256_file(os.path.abspath(__file__)), "python": platform.python_version()},
        "stage3_manifest_path": args.stage3_manifest,
        "verified_upstream_hashes": hashes,
        "ranking_label_mapping": {
            "primary_quality_ranking": "stage-3 columns rank_viewA_no_shared_key / rank_viewB_no_shared_key (Shared status does not influence rank)",
            "secondary_common_core_ranking": "stage-3 columns rank_viewA / rank_viewB (Shared-first; sensitivity/comparison only)",
            "not_combined_into_one_score": True,
        },
        "primary_quality_ranking_hierarchy": [
            "1. |span - 25| (25 = median span of verified Promega markers)",
            "2. Umap multi-read mean == 1.0 before NOT_ASSESSED (never imputed)",
            "3. higher read_measurable_fraction (100 bp read + 5 bp flank geometry)",
            "4. RepeatMasker tier, then flank flag (existing categories)",
            "5. objective Umap uniqueness fields (single-read coverage, multi-read min)",
            "6. sha256(marker_id) deterministic tie-break"],
        "chromosome_is_a_ranking_key": False, "chromosome_quotas": False, "chrY_removed": False,
        "sizes_previewed": args.sizes,
        "hcc1395_results_used": False, "classification_performed": False, "production_changes": False,
        "ranking_algorithm_rewritten": False, "stage3_outputs_modified": False,
        "scenarios": previews,
        "outputs": {"comparison_table": {"path": comp_path, "sha256": sha256_file(comp_path)}},
    }
    manifest_path = os.path.join(args.out_dir, "scenario_preview_stage4.v1.manifest.json")
    with open(manifest_path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print(f"wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
