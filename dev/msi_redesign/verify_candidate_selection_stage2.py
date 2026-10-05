"""
INDEPENDENT verification of candidate_selection_stage2.v1 outputs.

Deliberately does NOT import select_candidates_stage2.py. It re-derives every
rule from the raw 250 bp panel/traceability files with its own plain-csv
code and compares against (a) the traceability table the stage wrote and
(b) the summary numbers in the manifest, and re-hashes inputs/outputs.

Usage:
    python3 verify_candidate_selection_stage2.py --manifest MANIFEST.json \
        [--master-hashes tier2_raw=7bb1b9ed,tier2_filtered=d038ec32,view_A=f59e0a73 --output-dir DIR]
Exit code 0 = all checks passed.
"""

import csv
import hashlib
import json
import os
import sys
from collections import Counter


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def read_tsv(path):
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def verify(manifest_path, master_hashes=None, output_dir=None):
    failures = []

    def check(cond, msg):
        if not cond:
            failures.append(msg)

    with open(manifest_path, encoding="utf-8") as fh:
        m = json.load(fh)
    read_len, flank = 100, 5

    # 1. hashes of inputs and outputs
    for path, info in m["inputs"].items():
        check(sha(path) == info["sha256"], f"input hash mismatch: {path}")
    for name, info in m["outputs"].items():
        check(sha(info["path"]) == info["sha256"], f"output hash mismatch: {name}")
    with open(m["prior_manifest"]["path"], encoding="utf-8") as fh:
        prior = json.load(fh)
    check(sha(m["prior_manifest"]["path"]) == m["prior_manifest"]["sha256"], "prior manifest hash changed")
    for view in ("viewA_wes_compatible", "viewB_wgs"):
        for key in ("selected_panel", "traceability_all_loci"):
            info = prior["outputs"][view][key]
            check(sha(info["path"]) == info["sha256"], f"250bp output changed: {view}/{key}")
    if master_hashes and output_dir:
        for fname, prefix in master_hashes.items():
            check(sha(os.path.join(output_dir, fname)).startswith(prefix), f"master file changed: {fname}")

    # 2. re-derive membership and rules from the RAW 250bp panels
    panel_paths = [p for p in m["inputs"] if p.endswith(".tsv") and "traceability" not in p]
    a_path = [p for p in panel_paths if "viewA" in os.path.basename(p)][0]
    b_path = [p for p in panel_paths if "viewB" in os.path.basename(p)][0]
    raw_a = {r["marker_id"]: r for r in read_tsv(a_path)}
    raw_b = {r["marker_id"]: r for r in read_tsv(b_path)}
    union = set(raw_a) | set(raw_b)

    def rules(r):
        span = int(r["end"]) - int(r["start"])
        mono = r["kind"] == "mononucleotide"
        s3 = 15 <= span <= 40
        mean = r["umap_multiread_mean"]
        umap_ok = (mean == "") or float(mean) >= 1.0
        if span > read_len - 2 * flank:
            outcome = "FAIL_S1_SPAN_GT_90"
        elif not s3:
            band = "lt10" if span < 10 else "10_14" if span < 15 else "41_60" if span <= 60 else "gt60"
            outcome = f"FAIL_S3_SPAN_BAND_{band}"
        elif not umap_ok:
            outcome = "FAIL_S4_UMAP_BELOW_TOP_TIER"
        else:
            outcome = "PASS_ALL_CRITERIA"
        prefix = "MONO" if mono else "STR_COMPARISON"
        final = ("STRONG_" + ("MONO_PRIMARY" if mono else "STR_COMPARISON")) if outcome == "PASS_ALL_CRITERIA" \
            else ("NOT_STRONG_MONO" if mono else "NOT_STRONG_STR")
        return final, f"{prefix}:{outcome}", span

    # 3. compare with the traceability table the stage wrote
    trace = read_tsv(m["outputs"]["traceability_all_representatives"]["path"])
    check(len(trace) == len(union), f"row count {len(trace)} != union {len(union)}")
    check({t["marker_id"] for t in trace} == union, "traceability marker set != union of raw panels")
    check(len({t["marker_id"] for t in trace}) == len(trace), "duplicate marker rows in traceability")
    bad = 0
    for t in trace:
        mid = t["marker_id"]
        raw = raw_a.get(mid) or raw_b.get(mid)
        final, reason, span = rules(raw)
        cls = "A_and_B" if mid in raw_a and mid in raw_b else "A_only" if mid in raw_a else "B_only"
        if (t["final_class"], t["reason"], int(t["span"]), t["source_view_class"], t["representative_id"]) != \
                (final, reason, span, cls, mid):
            bad += 1
    check(bad == 0, f"{bad} traceability rows disagree with independently derived rules")

    strong_rows = read_tsv(m["outputs"]["strong_mono_primary"]["path"])
    comp_rows = read_tsv(m["outputs"]["strong_str_comparison"]["path"])
    check({r["marker_id"] for r in strong_rows} == {t["marker_id"] for t in trace if t["final_class"] == "STRONG_MONO_PRIMARY"},
          "strong_mono_primary file != rows classed STRONG_MONO_PRIMARY")
    check({r["marker_id"] for r in comp_rows} == {t["marker_id"] for t in trace if t["final_class"] == "STRONG_STR_COMPARISON"},
          "strong_str_comparison file != rows classed STRONG_STR_COMPARISON")

    # 4. recompute the manifest summary numbers independently
    groups = {
        "viewA": lambda mid: mid in raw_a, "viewB": lambda mid: mid in raw_b,
        "shared_A_and_B": lambda mid: mid in raw_a and mid in raw_b,
        "A_only": lambda mid: mid in raw_a and mid not in raw_b,
        "B_only": lambda mid: mid in raw_b and mid not in raw_a,
    }
    for g, pred in groups.items():
        ids = [mid for mid in union if pred(mid)]
        derived = Counter()
        chry = Counter()
        for mid in ids:
            raw = raw_a.get(mid) or raw_b.get(mid)
            final, _, span = rules(raw)
            derived[final] += 1
            if final == "STRONG_MONO_PRIMARY" and raw["chrom"] == "chrY":
                chry["mono"] += 1
        s = m["summary"][g]
        check(s["stage_attrition"]["S0_input_250bp_representatives"] == len(ids), f"{g}: S0 count differs")
        check(s["stage_attrition"]["MONO_PRIMARY"]["S4_pass_umap_top_tier_or_not_assessed_=_STRONG"]
              == derived["STRONG_MONO_PRIMARY"], f"{g}: strong mono count differs")
        check(s["stage_attrition"]["STR_COMPARISON"]["S4_pass_umap_top_tier_or_not_assessed_=_STRONG"]
              == derived["STRONG_STR_COMPARISON"], f"{g}: strong STR count differs")
        check(s["chrY_impact"]["strong_mono_chrY"] == chry["mono"], f"{g}: chrY strong count differs")
        check(s["chrY_impact"]["strong_mono_with_chrY"] - s["chrY_impact"]["strong_mono_without_chrY"] == chry["mono"],
              f"{g}: chrY with/without arithmetic differs")
    check(m["view_relationship"].get("A_and_B", 0) == len(set(raw_a) & set(raw_b)), "shared count differs")
    check(m["view_relationship"].get("A_only", 0) == len(set(raw_a) - set(raw_b)), "A-only count differs")
    check(m["view_relationship"].get("B_only", 0) == len(set(raw_b) - set(raw_a)), "B-only count differs")
    check(m["production_changes"] is False and m["hcc1395_results_used"] is False
          and m["classification_performed"] is False and m["final_panel_size_fixed"] is False,
          "manifest safety flags not all false")
    return failures


def main(argv):
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--master-hashes", default="")
    ap.add_argument("--output-dir", default=None)
    args = ap.parse_args(argv)
    masters = dict(kv.split("=") for kv in args.master_hashes.split(",") if kv)
    files = {"tier2_raw": "tier2_raw.v1.tsv", "tier2_filtered": "tier2_filtered.v1.tsv",
             "view_A": "view_A_wes_compatible.v1.tsv"}
    master_hashes = {files[k]: v for k, v in masters.items()} if masters else None
    failures = verify(args.manifest, master_hashes, args.output_dir)
    if failures:
        print("VERIFICATION FAILED")
        for f in failures:
            print(" -", f)
        return 1
    print("VERIFICATION PASSED: all independent checks agree")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
