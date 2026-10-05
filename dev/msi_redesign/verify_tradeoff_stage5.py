"""
INDEPENDENT verification of tradeoff_stage5.v1 outputs.

Imports NO Stage 1-5 implementation code. From the raw hash-verified Stage-3
union table it independently
  1. re-hashes the whole chain (stage-4, stage-3, stage-2, 250 bp),
  2. rebuilds every quality-first / shared-first scenario membership from the
     Stage-3 rank columns and recounts every reported statistic,
  3. re-derives the shared-core order from the ATTRIBUTE columns (not from the
     rank columns) and checks it equals the rank-based order,
  4. rebuilds every hybrid (core k + extension to T) with its own code and
     recounts composition, fractions and overlaps,
  5. recomputes the Pareto dominance fronts (both size conventions and within
     fixed size) and the 15-20 bp band positions,
  6. checks safety flags.

Usage: python3 verify_tradeoff_stage5.py --manifest tradeoff_stage5.v1.manifest.json
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


def load(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def tsv(path):
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def flt(x):
    return float(x) if x != "" else None


def count_stats(members, n_avail):
    n = len(members)
    cls = Counter(r["source_view_class"] for r in members)
    sub = {b: sum(1 for r in members if r["span_subband"] == b) for b in ("15_20", "21_30", "31_40")}
    return {
        "n": n, "shared": cls.get("A_and_B", 0), "A_only": cls.get("A_only", 0), "B_only": cls.get("B_only", 0),
        "span_subband": sub,
        "umap_status": dict(Counter(r["umap_status"] for r in members)),
        "rmsk_self_category": dict(Counter(r["rmsk_self_category"] for r in members)),
        "flank_other_repeat_flag": dict(Counter(r["flank_other_repeat_flag"] for r in members)),
        "chromosome": dict(sorted(Counter(r["chrom"] for r in members).items())),
        "chrY": sum(1 for r in members if r["chrom"] == "chrY"),
        "shared_fraction": round(cls.get("A_and_B", 0) / n, 6),
        "fraction_21_30": round(sub["21_30"] / n, 6), "fraction_15_20": round(sub["15_20"] / n, 6),
        "pct_of_available_strong_set": round(100.0 * n / n_avail, 3),
    }


def attribute_key(r):
    cov, mn = flt(r["K6_singleread_cov"]), flt(r["K6_multiread_min"])
    return (float(r["K2_promega_span_distance"]), 0 if r["umap_status"] == "TIER0_MEAN_1.0" else 1,
            -float(r["K4_read_measurable_fraction"]), int(r["K5_rmsk_tier"]), int(r["K5_flank_flag"]),
            (0, -cov) if cov is not None else (1, 0.0), (0, -mn) if mn is not None else (1, 0.0),
            r["K7_tie_hash"])


def dominated_map(configs, mode):
    def vec(c):
        v = [c["shared_fraction"], c["fraction_21_30"], -c["quality_percentile"]]
        if mode == "benefit":
            v.append(c["size"])
        elif mode == "cost":
            v.append(-c["size"])
        return v
    vecs = {c["id"]: vec(c) for c in configs}
    front, dom = [], {}
    for c in configs:
        doms = sorted(o["id"] for o in configs if o["id"] != c["id"]
                      and all(x >= y for x, y in zip(vecs[o["id"]], vecs[c["id"]]))
                      and any(x > y for x, y in zip(vecs[o["id"]], vecs[c["id"]])))
        if doms:
            dom[c["id"]] = doms[0]
        else:
            front.append(c["id"])
    return front, dom


def band_positions(order):
    first = next((i + 1 for i, r in enumerate(order) if r["span_subband"] == "15_20"), None)
    last = max((i + 1 for i, r in enumerate(order) if r["span_subband"] == "21_30"), default=None)
    cross, cum = {}, 0
    for i, r in enumerate(order):
        cum += r["span_subband"] == "15_20"
        for lvl in (0.01, 0.05, 0.10, 0.25, 0.50):
            if str(lvl) not in cross and cum / (i + 1) >= lvl:
                cross[str(lvl)] = i + 1
    return first, last, sum(1 for r in order if r["span_subband"] == "21_30"), cross


def verify(manifest_path):
    fails = []

    def check(cond, msg):
        if not cond:
            fails.append(msg)

    m5 = load(manifest_path)
    m4 = load(m5["stage4_manifest_path"])
    m3 = load(m4["stage3_manifest_path"])
    m2 = load(m3["stage2_manifest_path"])
    prior = load(m2["prior_manifest"]["path"])
    check(sha(m5["stage4_manifest_path"]) == m5["verified_upstream_hashes"]["stage4_manifest"], "stage-4 manifest changed")
    for name, info in m4["outputs"].items():
        check(sha(info["path"]) == info["sha256"], f"stage-4 output changed: {name}")
    check(sha(m4["stage3_manifest_path"]) == m4["verified_upstream_hashes"]["stage3_manifest"], "stage-3 manifest changed")
    for name, info in m3["outputs"].items():
        check(sha(info["path"]) == info["sha256"], f"stage-3 output changed: {name}")
    check(sha(m3["stage2_manifest_path"]) == m3["inputs_verified_hashes"]["stage2_manifest"], "stage-2 manifest changed")
    for name, info in m2["outputs"].items():
        check(sha(info["path"]) == info["sha256"], f"stage-2 output changed: {name}")
    check(sha(m2["prior_manifest"]["path"]) == m2["prior_manifest"]["sha256"], "250 bp manifest changed")
    for view, outs in prior["outputs"].items():
        for key, info in outs.items():
            check(sha(info["path"]) == info["sha256"], f"250 bp output changed: {view}/{key}")
    for name, info in m5["outputs"].items():
        check(sha(info["path"]) == info["sha256"], f"stage-5 output hash mismatch: {name}")
    for flag in ("final_panel_size_fixed", "panel_size_selected", "membership_files_written", "hcc1395_results_used",
                 "classification_performed", "production_changes", "upstream_outputs_modified"):
        check(m5[flag] is False, f"safety flag not False: {flag}")

    rows = tsv(m3["outputs"]["ranked_union"]["path"])
    for r in rows:
        r["inA"], r["inB"] = r["in_viewA_strong"] == "True", r["in_viewB_strong"] == "True"
        r["shared"] = r["source_view_class"] == "A_and_B"
        for c in ("rank_viewA", "rank_viewB", "rank_viewA_no_shared_key", "rank_viewB_no_shared_key"):
            r[c] = int(r[c]) if r[c] != "" else None

    def order(flag, col):
        return sorted((r for r in rows if r[flag]), key=lambda r: r[col])

    aq, as_, bq, bs = (order("inA", "rank_viewA_no_shared_key"), order("inA", "rank_viewA"),
                       order("inB", "rank_viewB_no_shared_key"), order("inB", "rank_viewB"))
    n_a, n_b = len(aq), len(bq)
    p = m5["parameters"]
    check((p["n_viewA_strong"], p["n_viewB_strong"]) == (n_a, n_b), "view sizes differ")

    # 2. scenarios
    for view, (q, s, n_av) in (("viewA", (aq, as_, n_a)), ("viewB", (bq, bs, n_b))):
        for size in p["sizes"]:
            got = m5["scenarios"][view][str(size)]
            for label, lst in (("quality_first", q), ("shared_first", s)):
                want = count_stats(lst[:size], n_av)
                for k, v in want.items():
                    check(got[label][k] == v, f"{view}/{size}/{label}: {k} differs")
            both = len({r["marker_id"] for r in q[:size]} & {r["marker_id"] for r in s[:size]})
            check(got["loci_in_both_strategies"] == both, f"{view}/{size}: overlap differs")

    # 3. shared core order from ATTRIBUTES (independent of rank columns)
    shared_rows = [r for r in rows if r["shared"]]
    core_attr = sorted(shared_rows, key=attribute_key)
    core_rank = sorted(shared_rows, key=lambda r: r["rank_viewB_no_shared_key"])
    check([r["marker_id"] for r in core_attr] == [r["marker_id"] for r in core_rank],
          "shared core order from attributes differs from the rank-based order")
    check(len(shared_rows) == p["n_shared_strong"], "shared count differs")
    kmax = max(p["hybrid_cores"])
    check([r["marker_id"] for r in core_rank[:kmax]] == [r["marker_id"] for r in bs[:kmax]],
          "core is not the shared-first prefix")

    # 4. hybrids
    configs = []
    for size in p["sizes"]:
        for tag, lst in (("QF", bq), ("SF", bs)):
            st = count_stats(lst[:size], n_b)
            configs.append({"id": f"{tag}_{size}", "size": size, "shared_fraction": st["shared_fraction"],
                            "fraction_21_30": st["fraction_21_30"],
                            "quality_percentile": round(sum(r["rank_viewB_no_shared_key"] for r in lst[:size]) / size / n_b, 6)})
    for k in p["hybrid_cores"]:
        for total in p["hybrid_totals"]:
            if total <= k:
                continue
            core = core_rank[:k]
            core_ids = {r["marker_id"] for r in core}
            ext = [r for r in bq if r["marker_id"] not in core_ids][: total - k]
            members = core + ext
            ids = [r["marker_id"] for r in members]
            check(len(ids) == total and len(set(ids)) == total, f"hybrid c{k}/T{total}: duplicates or wrong size")
            got = m5["hybrid"][f"core_{k}"][f"total_{total}"]
            st = count_stats(members, n_b)
            for key, v in st.items():
                check(got["stats"][key] == v, f"hybrid c{k}/T{total}: {key} differs")
            check((got["common_core_size"], got["wgs_extension_size"], got["total_wgs_size"], got["total_unique_loci"])
                  == (k, total - k, total, total), f"hybrid c{k}/T{total}: size fields differ")
            shared_members = [r for r in members if r["shared"]]
            check(got["shared_loci_in_panel"] == len(shared_members), f"hybrid c{k}/T{total}: shared count differs")
            check(got["shared_loci_in_extension_beyond_core"] == sum(1 for r in ext if r["shared"]),
                  f"hybrid c{k}/T{total}: shared-in-extension differs")
            check(got["wgs_only_loci_in_panel"] == sum(1 for r in members if not r["shared"]),
                  f"hybrid c{k}/T{total}: wgs-only differs")
            for label, part in (("core_span_subband", core), ("extension_span_subband", ext),
                                ("shared_subset_span_subband", shared_members)):
                want = {b: sum(1 for r in part if r["span_subband"] == b) for b in ("15_20", "21_30", "31_40")}
                check(got[label] == want, f"hybrid c{k}/T{total}: {label} differs")
            check(got["overlap_with_quality_first_top_total"] == len(set(ids) & {r["marker_id"] for r in bq[:total]}),
                  f"hybrid c{k}/T{total}: overlap with quality-first differs")
            check(got["overlap_with_shared_first_top_total"] == len(set(ids) & {r["marker_id"] for r in bs[:total]}),
                  f"hybrid c{k}/T{total}: overlap with shared-first differs")
            qp = round(sum(r["rank_viewB_no_shared_key"] for r in members) / total / n_b, 6)
            check(got["mean_quality_rank_percentile_viewB"] == qp, f"hybrid c{k}/T{total}: quality percentile differs")
            configs.append({"id": f"HY_c{k}_T{total}", "size": total, "shared_fraction": st["shared_fraction"],
                            "fraction_21_30": st["fraction_21_30"], "quality_percentile": qp})

    # 5. Pareto
    front_b, dom_b = dominated_map(configs, "benefit")
    front_c, dom_c = dominated_map(configs, "cost")
    P = m5["pareto"]
    check(sorted(P["front_size_as_benefit"]) == sorted(front_b) and P["dominated_by_size_as_benefit"] == dom_b,
          "Pareto front (size as benefit) differs")
    check(sorted(P["front_size_as_cost"]) == sorted(front_c) and P["dominated_by_size_as_cost"] == dom_c,
          "Pareto front (size as cost) differs")
    for size in sorted({c["size"] for c in configs}):
        f, d = dominated_map([c for c in configs if c["size"] == size], "fixed")
        got = P["within_fixed_size"][str(size)]
        check(sorted(got["front"]) == sorted(f) and got["dominated_by"] == d, f"within-size front differs at {size}")
    check({c["id"] for c in P["configs"]} == {c["id"] for c in configs}, "Pareto config set differs")

    # 6. band positions
    for key, lst in (("viewA_shared_first", as_), ("viewB_shared_first", bs), ("viewB_quality_first", bq)):
        first, last, n2130, cross = band_positions(lst)
        got = m5["shared_first_band_analysis"][key]
        check((got["first_15_20_rank"], got["last_21_30_rank"], got["n_21_30_in_list"]) == (first, last, n2130),
              f"{key}: band positions differ")
        check(got["first_rank_where_cumulative_15_20_share_reaches"] == cross, f"{key}: share crossings differ")
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
    print("VERIFICATION PASSED: independent memberships, statistics, hybrids, Pareto fronts and hashes agree")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
