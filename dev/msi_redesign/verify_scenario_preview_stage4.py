"""
INDEPENDENT verification of scenario_preview_stage4.v1 outputs.

Does NOT import scenario_preview_stage4.py (nor any ranking code). It
  1. re-hashes the whole upstream chain (stage-3, stage-2, 250 bp outputs),
  2. re-derives every scenario statistic from the raw Stage-3 union table
     with its own code and compares it with the manifest,
  3. checks, from the attribute columns alone, that both rankings are
     genuinely ordered by their declared hierarchy (adjacent-pair
     monotonicity) -- in particular that the PRIMARY QUALITY ranking is
     independent of Shared status, and that the common-core ranking is
     Shared-first,
  4. checks the safety flags.

Usage: python3 verify_scenario_preview_stage4.py --manifest scenario_preview_stage4.v1.manifest.json
Exit code 0 = all checks passed.
"""

import csv
import hashlib
import json
import sys
from collections import Counter


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def tsv(path):
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def load(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def flt(x):
    return float(x) if x != "" else None


def order_key(r, anchor, with_shared):
    cov, mn = flt(r["umap_singleread_window_coverage"]), flt(r["umap_multiread_min"])
    key = []
    if with_shared:
        key.append(0 if r["source_view_class"] == "A_and_B" else 1)
    key += [abs(int(r["span"]) - anchor),
            0 if r["umap_status"] == "TIER0_MEAN_1.0" else 1,
            -float(r["read_measurable_fraction"]),
            int(r["rmsk_tier"]), 1 if r["flank_other_repeat_flag"] == "True" else 0,
            (0, -cov) if cov is not None else (1, 0.0),
            (0, -mn) if mn is not None else (1, 0.0),
            r["K7_tie_hash"]]
    return tuple(key)


def verify(manifest_path):
    fails = []

    def check(cond, msg):
        if not cond:
            fails.append(msg)

    m4 = load(manifest_path)
    m3 = load(m4["stage3_manifest_path"])
    check(sha(m4["stage3_manifest_path"]) == m4["verified_upstream_hashes"]["stage3_manifest"], "stage-3 manifest changed")
    for name, info in m3["outputs"].items():
        check(sha(info["path"]) == info["sha256"], f"stage-3 output changed: {name}")
    m2 = load(m3["stage2_manifest_path"])
    check(sha(m3["stage2_manifest_path"]) == m3["inputs_verified_hashes"]["stage2_manifest"], "stage-2 manifest changed")
    for name, info in m2["outputs"].items():
        check(sha(info["path"]) == info["sha256"], f"stage-2 output changed: {name}")
    check(sha(m2["prior_manifest"]["path"]) == m2["prior_manifest"]["sha256"], "250 bp manifest changed")
    prior = load(m2["prior_manifest"]["path"])
    for view, outs in prior["outputs"].items():
        for key, info in outs.items():
            check(sha(info["path"]) == info["sha256"], f"250 bp output changed: {view}/{key}")
    for name, info in m4["outputs"].items():
        check(sha(info["path"]) == info["sha256"], f"stage-4 output hash mismatch: {name}")

    for flag in ("final_panel_size_fixed", "panel_size_selected", "hcc1395_results_used", "classification_performed",
                 "production_changes", "ranking_algorithm_rewritten", "stage3_outputs_modified",
                 "chromosome_is_a_ranking_key", "chromosome_quotas", "chrY_removed"):
        check(m4[flag] is False, f"safety flag not False: {flag}")
    check(m4["ranking_label_mapping"]["not_combined_into_one_score"] is True, "rankings were combined")

    anchor = m3["parameters"]["promega_span_anchor_bp"]
    rows = tsv(m3["outputs"]["ranked_union"]["path"])

    # K1 column consistent with the view class (input sanity)
    check(all((r["K1_shared_tier"] == "0") == (r["source_view_class"] == "A_and_B") for r in rows),
          "K1_shared_tier inconsistent with source_view_class")

    views = {"viewA": ("in_viewA_strong", "rank_viewA_no_shared_key", "rank_viewA"),
             "viewB": ("in_viewB_strong", "rank_viewB_no_shared_key", "rank_viewB")}
    for view, (in_col, primary_col, shared_col) in views.items():
        vrows = [r for r in rows if r[in_col] == "True"]
        n_av = len(vrows)
        by_primary = sorted(vrows, key=lambda r: int(r[primary_col]))
        by_shared = sorted(vrows, key=lambda r: int(r[shared_col]))
        check([int(r[primary_col]) for r in by_primary] == list(range(1, n_av + 1)), f"{view}: primary ranks not 1..N")
        check([int(r[shared_col]) for r in by_shared] == list(range(1, n_av + 1)), f"{view}: shared-first ranks not 1..N")
        # adjacent-pair monotonicity from attributes alone
        bad_p = sum(1 for a, b in zip(by_primary, by_primary[1:])
                    if order_key(a, anchor, False) > order_key(b, anchor, False))
        bad_s = sum(1 for a, b in zip(by_shared, by_shared[1:])
                    if order_key(a, anchor, True) > order_key(b, anchor, True))
        check(bad_p == 0, f"{view}: primary quality ranking violates its hierarchy at {bad_p} adjacent pairs")
        check(bad_s == 0, f"{view}: shared-first ranking violates its hierarchy at {bad_s} adjacent pairs")
        # the primary ranking must NOT be shared-first (else the policy change did not happen)
        first_nonshared = next((i for i, r in enumerate(by_primary) if r["source_view_class"] != "A_and_B"), None)
        if view == "viewB":
            check(first_nonshared is not None and first_nonshared < len([r for r in vrows if r["source_view_class"] == "A_and_B"]),
                  "viewB primary quality ranking still looks shared-first")

        for size_text, entry in m4["scenarios"][view].items():
            size = int(size_text)
            check(entry["n_available_strong_set"] == n_av, f"{view}/{size}: available count differs")
            sets = {}
            for label, ordered in (("primary_quality_ranking", by_primary), ("secondary_common_core_ranking", by_shared)):
                top = ordered[:size]
                sets[label] = {r["marker_id"] for r in top}
                got = entry[label]
                cls = Counter(r["source_view_class"] for r in top)
                want = {
                    "n": len(top), "shared": cls.get("A_and_B", 0), "A_only": cls.get("A_only", 0),
                    "B_only": cls.get("B_only", 0),
                    "span_subband": {b: sum(1 for r in top if r["span_subband"] == b) for b in ("15_20", "21_30", "31_40")},
                    "umap_status": dict(Counter(r["umap_status"] for r in top)),
                    "rmsk_self_category": dict(Counter(r["rmsk_self_category"] for r in top)),
                    "flank_other_repeat_flag": dict(Counter(r["flank_other_repeat_flag"] for r in top)),
                    "chromosome": dict(sorted(Counter(r["chrom"] for r in top).items())),
                    "chrY": sum(1 for r in top if r["chrom"] == "chrY"),
                }
                for k, v in want.items():
                    check(got[k] == v, f"{view}/{size}/{label}: {k} differs from independent recount")
                check(abs(got["pct_of_available_strong_set"] - round(100.0 * len(top) / n_av, 3)) < 1e-9,
                      f"{view}/{size}/{label}: percentage differs")
            check(entry["loci_in_both_rankings_top_n"] ==
                  len(sets["primary_quality_ranking"] & sets["secondary_common_core_ranking"]),
                  f"{view}/{size}: overlap differs")
    return fails


def main(argv):
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    args = ap.parse_args(argv)
    fails = verify(args.manifest)
    if fails:
        print("VERIFICATION FAILED")
        for f in fails:
            print(" -", f)
        return 1
    print("VERIFICATION PASSED: upstream hashes, independent scenario recounts and ranking-order checks agree")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
