"""
INDEPENDENT verification of candidate_ranking_stage3.v1 outputs.

Does NOT import rank_candidates_stage3.py. It re-reads the raw immutable
inputs (stage-2 strong file, annotation table, Promega tier-1 file), rebuilds
the ranking with its own comparator (functools.cmp_to_key, not the
implementation's tuple keys), and compares every rank, count, hash and
prefix composition with what the stage wrote.

Usage: python3 verify_ranking_stage3.py --manifest candidate_ranking_stage3.v1.manifest.json
Exit code 0 = all checks passed.
"""

import csv
import functools
import hashlib
import json
import statistics
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


def verify(manifest_path):
    fails = []

    def check(cond, msg):
        if not cond:
            fails.append(msg)

    with open(manifest_path, encoding="utf-8") as f:
        m = json.load(f)
    with open(m["stage2_manifest_path"], encoding="utf-8") as f:
        s2 = json.load(f)

    # immutable inputs and stage-2 outputs unchanged
    check(sha(m["stage2_manifest_path"]) == m["inputs_verified_hashes"]["stage2_manifest"], "stage-2 manifest changed")
    for name, info in s2["outputs"].items():
        check(sha(info["path"]) == info["sha256"], f"stage-2 output changed: {name}")
    for name, info in m["outputs"].items():
        check(sha(info["path"]) == info["sha256"], f"ranking output hash mismatch: {name}")
    for flag in ("hcc1395_results_used", "classification_performed", "production_changes",
                 "final_panel_size_fixed", "panel_size_selected", "stage2_outputs_modified"):
        check(m[flag] is False, f"safety flag not False: {flag}")
    check(m["parameters"]["chromosome_is_a_ranking_key"] is False and m["parameters"]["chromosome_quotas"] is False,
          "chromosome used as key/quota")

    # raw inputs
    strong = tsv(s2["outputs"]["strong_mono_primary"]["path"])
    ids = {r["marker_id"] for r in strong}
    with open(s2["prior_manifest"]["path"], encoding="utf-8") as f:
        prior = json.load(f)
    annotated_path = prior["inputs"]["annotated"]["path"]
    check(sha(annotated_path) == prior["inputs"]["annotated"]["sha256"], "annotation table changed")
    ann = {}
    with open(annotated_path, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["marker_id"] in ids:
                ann[row["marker_id"]] = row
    check(len(ann) == len(ids), "strong loci missing from annotation table")
    tier1_spans = list(int(r["end"]) - int(r["start"]) for r in tsv(m["parameters"].get("promega_tier1_path",
                       "dev/msi_redesign/output/msi_markers.promega_tier1.tsv")))
    anchor = statistics.median(tier1_spans)
    check(anchor == m["parameters"]["promega_span_anchor_bp"], "Promega anchor differs")

    rows = []
    for r in strong:
        a = ann[r["marker_id"]]
        span = int(r["end"]) - int(r["start"])
        sc = a["umap_singleread_window_coverage"]
        mm = a["umap_multiread_min"]
        rows.append({
            "id": r["marker_id"], "chrom": r["chrom"], "span": span,
            "shared": r["in_viewA_reps"] == "True" and r["in_viewB_reps"] == "True",
            "inA": r["in_viewA_reps"] == "True", "inB": r["in_viewB_reps"] == "True",
            "assessed": r["umap_multiread_mean"] != "" and float(r["umap_multiread_mean"]) >= 1.0,
            "geom": float(r["read_measurable_fraction"]), "rmsk": int(r["rmsk_tier"]),
            "flank": r["flank_other_repeat_flag"] == "True",
            "sc": float(sc) if sc != "" else None, "mm": float(mm) if mm != "" else None,
            "hash": hashlib.sha256(r["marker_id"].encode("utf-8")).hexdigest(),
        })

    def lower(x, y):
        return -1 if x < y else (1 if x > y else 0)

    def cmp(a, b, use_shared):
        steps = []
        if use_shared:
            steps.append((0 if a["shared"] else 1, 0 if b["shared"] else 1))
        steps.append((abs(a["span"] - anchor), abs(b["span"] - anchor)))
        steps.append((0 if a["assessed"] else 1, 0 if b["assessed"] else 1))
        steps.append((b["geom"], a["geom"]))                        # higher geometry first
        steps.append((a["rmsk"], b["rmsk"]))
        steps.append((1 if a["flank"] else 0, 1 if b["flank"] else 0))
        for f in ("sc", "mm"):
            steps.append((1 if a[f] is None else 0, 1 if b[f] is None else 0))   # missing sorts after present
            if a[f] is not None and b[f] is not None:
                steps.append((b[f], a[f]))                          # higher first
        steps.append((a["hash"], b["hash"]))
        for x, y in steps:
            c = lower(x, y)
            if c:
                return c
        return 0

    expected_counts = {"viewA": 0, "viewB": 0}
    for view, flag in (("viewA", "inA"), ("viewB", "inB")):
        subset = [r for r in rows if r[flag]]
        expected_counts[view] = len(subset)
        for label, use_shared, col in (("primary", True, f"rank_{view}"),
                                       ("no_shared", False, f"rank_{view}_no_shared_key")):
            ordered = sorted(subset, key=functools.cmp_to_key(lambda a, b: cmp(a, b, use_shared)))
            want = {r["id"]: i + 1 for i, r in enumerate(ordered)}
            got_rows = tsv(m["outputs"]["ranked_union"]["path"])
            got = {r["marker_id"]: r[col] for r in got_rows if r[col] != ""}
            check(len(got) == len(want), f"{view}/{label}: ranked count differs")
            check(all(str(want[k]) == got.get(k) for k in want), f"{view}/{label}: ranks differ from independent ranking")
            check(sorted(int(v) for v in got.values()) == list(range(1, len(want) + 1)),
                  f"{view}/{label}: ranks are not a permutation 1..N")
            if label == "primary":
                file_rows = tsv(m["outputs"][f"ranked_{view}"]["path"])
                check([r["marker_id"] for r in file_rows] == [r["id"] for r in ordered],
                      f"{view}: per-view ranked file order differs")
                check(all(int(r[col]) == i + 1 for i, r in enumerate(file_rows)), f"{view}: file rank column not sequential")
                # prefix composition (top 10%)
                import math
                k = max(1, math.ceil(len(ordered) * 0.10))
                chrom = dict(sorted(Counter(r["chrom"] for r in ordered[:k]).items()))
                check(chrom == m["report"][view]["primary_ranking_prefixes"]["top_10pct"]["chromosome"],
                      f"{view}: top-10% chromosome composition differs")
                check(sum(1 for r in ordered if r["chrom"] == "chrY") == m["report"][view]["full_set_composition"]["chrY"],
                      f"{view}: chrY count differs")
    check(expected_counts["viewA"] == m["counts"]["viewA"] and expected_counts["viewB"] == m["counts"]["viewB"],
          "view counts differ")
    check(sum(1 for r in rows if r["shared"]) == m["counts"]["shared"], "shared count differs")
    check(len(rows) == m["counts"]["union"], "union count differs")
    # not-assessed loci stay identifiable and never imputed
    got_rows = tsv(m["outputs"]["ranked_union"]["path"])
    for r in got_rows:
        if r["umap_status"] == "NOT_ASSESSED":
            check(r["umap_multiread_min"] == "" and r["umap_multiread_mean"] == "", "NOT_ASSESSED locus has an imputed Umap value")
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
    print("VERIFICATION PASSED: independent ranking, counts, hashes and prefix compositions agree")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
