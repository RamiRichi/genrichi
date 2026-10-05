"""
GenRichi MSI Methodology Redesign -- DEVELOPMENT / RUO TRADE-OFF ANALYSIS.
NOT PRODUCTION CODE. NOT A CLINICALLY VALIDATED PANEL. No MSI classification.

Stage 5: trade-off / Pareto analysis between
  A. QUALITY-FIRST   = Stage-3 columns rank_*_no_shared_key (unchanged)
  B. SHARED-FIRST    = Stage-3 columns rank_* (unchanged; common-core)
  C. HYBRID          = a ranked shared core, then WGS supplementation.
It selects no panel and no size, writes no membership files, and makes no
performance claim. Repeat-length categories are DEVELOPMENT PREFERENCES
(the 25 bp anchor is the median span of the five verified Promega markers),
not validated MSI performance thresholds. Terminology is deliberately
neutral: ranking priority, development preference, candidate composition,
technical suitability, repeat-length category.

Immutable inputs: the Stage-4 manifest (pinned by hash) and everything
upstream of it are re-hashed before any work; the Stage-3 union table is
read only.

HYBRID DEFINITION (explicit, deterministic, duplicate-free):
  core order   = the shared strong loci (in both View A and View B strong
                 sets; 6,038) ordered by the Stage-3 quality-first rank
                 (rank_viewB_no_shared_key). Shared status is a MEMBERSHIP
                 filter for the core only -- it is not a ranking key; the
                 order comes purely from the quality attributes. (For any
                 core size k <= 6,038 this equals the Stage-3 shared-first
                 top-k; that identity is asserted.)
  core(k)      = first k loci of the core order.
  extension    = the first (T - k) loci of the View-B quality-first ranking
                 that are NOT already in core(k). Extension loci may be
                 B-only or additional shared loci that rank highly on
                 quality; no locus appears twice.
  WGS panel    = core(k) + extension, total exactly T unique loci.

PARETO. Configurations: quality-first(T), shared-first(T) for every scenario
size, and hybrid(k, T). All live in the View-B strong set, so one common
reference (View-B quality rank) is valid. Objectives (a config is
DOMINATED if another is >= on all and > on at least one):
  - shared fraction (shared loci / size)      : higher is preferred
  - 21-30 bp fraction                         : higher is preferred
  - quality percentile = mean(rank_viewB_no_shared_key) / N_viewB
                                              : lower is preferred
  - panel size: direction is a matter of design, so BOTH conventions are
    reported: 'size_as_benefit' (larger preferred, more loci) and
    'size_as_cost' (smaller preferred, simpler panel). Fronts are also
    reported within each fixed size. No winner is chosen.

No HCC1395 result, no BAM, no classification, no production file.
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

import scenario_preview_stage4 as s4   # reused unchanged (stats + hash-chain verification)

SCRIPT_VERSION = "tradeoff_pareto_stage5.v1"
SCENARIO_SIZES = [500, 750, 1000, 1250, 1500, 2000, 3000, 5000]
HYBRID_CORES = [500, 750, 1000, 1250, 1500]
HYBRID_TOTALS = [2000, 3000, 5000]
BAND_SHARE_LEVELS = [0.01, 0.05, 0.10, 0.25, 0.50]   # descriptive reference levels, not thresholds
SUBBANDS = ["15_20", "21_30", "31_40"]


def sha256_file(path):
    return s4.sha256_file(path)


def verify_upstream(stage4_manifest_path, expect_sha256=None):
    """Refuses (SystemExit) unless the Stage-4 manifest, its outputs and the
    whole chain below it still match their recorded hashes."""
    with open(stage4_manifest_path, encoding="utf-8") as fh:
        m4 = json.load(fh)
    hashes = {"stage4_manifest": sha256_file(stage4_manifest_path)}
    if expect_sha256 and hashes["stage4_manifest"] != expect_sha256:
        raise SystemExit("stage-4 manifest hash differs from the expected value")
    for flag in ("final_panel_size_fixed", "panel_size_selected", "hcc1395_results_used",
                 "classification_performed", "production_changes", "ranking_algorithm_rewritten",
                 "stage3_outputs_modified", "chromosome_is_a_ranking_key", "chromosome_quotas", "chrY_removed"):
        if m4[flag] is not False:
            raise SystemExit(f"stage-4 safety flag {flag} is not False")
    for name, info in m4["outputs"].items():
        got = sha256_file(info["path"])
        if got != info["sha256"]:
            raise SystemExit(f"stage-4 output {name} does not match its recorded hash")
        hashes[f"stage4_{name}"] = got
    m3, chain = s4.verify_chain(m4["stage3_manifest_path"])
    if chain["stage3_manifest"] != m4["verified_upstream_hashes"]["stage3_manifest"]:
        raise SystemExit("stage-3 manifest differs from the one stage 4 verified")
    hashes.update(chain)
    return m4, m3, hashes


def to_int(x):
    return int(x) if x != "" else None


def load_union(path):
    rows = s4.read_tsv(path)
    for r in rows:
        for c in ("rank_viewA", "rank_viewB", "rank_viewA_no_shared_key", "rank_viewB_no_shared_key"):
            r[c] = to_int(r[c])
        r["_shared"] = r["source_view_class"] == "A_and_B"
        r["_inA"] = r["in_viewA_strong"] == "True"
        r["_inB"] = r["in_viewB_strong"] == "True"
    return rows


def ordered(rows, in_key, rank_col):
    return sorted((r for r in rows if r[in_key]), key=lambda r: r[rank_col])


def stats_for(members, n_available):
    s = s4.scenario_stats(members, n_available)
    n = s["n"]
    s["shared_fraction"] = round(s["shared"] / n, 6) if n else None
    s["wgs_only_or_view_only_fraction"] = round((s["A_only"] + s["B_only"]) / n, 6) if n else None
    s["fraction_21_30"] = round(s["span_subband"]["21_30"] / n, 6) if n else None
    s["fraction_15_20"] = round(s["span_subband"]["15_20"] / n, 6) if n else None
    return s


def core_order(rows):
    """Shared strong loci ordered by the Stage-3 quality-first rank. Shared
    status only filters membership; it never enters the ordering."""
    core = sorted((r for r in rows if r["_shared"]), key=lambda r: r["rank_viewB_no_shared_key"])
    core_a = sorted((r for r in rows if r["_shared"]), key=lambda r: r["rank_viewA_no_shared_key"])
    if [r["marker_id"] for r in core] != [r["marker_id"] for r in core_a]:
        raise SystemExit("shared-locus order differs between the View A and View B quality rankings")
    return core


def hybrid_members(core_sorted, b_quality_order, k, total):
    core = core_sorted[:k]
    core_ids = {r["marker_id"] for r in core}
    if len(core) < k:
        raise ValueError("core larger than the shared set")
    extension = []
    for r in b_quality_order:
        if r["marker_id"] in core_ids:
            continue
        extension.append(r)
        if len(extension) == total - k:
            break
    members = core + extension
    if len(members) != total or len({r["marker_id"] for r in members}) != total:
        raise ValueError("hybrid membership is not exactly `total` unique loci")
    return core, extension, members


def quality_percentile(members, n_b):
    return round(sum(r["rank_viewB_no_shared_key"] for r in members) / len(members) / n_b, 6)


def dominates(a, b, keys):
    ge = all(a[k] >= b[k] for k in keys)
    gt = any(a[k] > b[k] for k in keys)
    return ge and gt


def pareto_front(configs, mode):
    """mode: 'size_as_benefit' | 'size_as_cost' | 'fixed_size'."""
    def vec(c):
        v = {"shared_fraction": c["shared_fraction"], "fraction_21_30": c["fraction_21_30"],
             "neg_quality_percentile": -c["quality_percentile"]}
        if mode == "size_as_benefit":
            v["size_score"] = c["size"]
        elif mode == "size_as_cost":
            v["size_score"] = -c["size"]
        return v
    vecs = {c["id"]: vec(c) for c in configs}
    keys = list(next(iter(vecs.values())).keys())
    front, dominated_by = [], {}
    for c in configs:
        doms = sorted(o["id"] for o in configs if o["id"] != c["id"] and dominates(vecs[o["id"]], vecs[c["id"]], keys))
        if doms:
            dominated_by[c["id"]] = doms[0]
        else:
            front.append(c["id"])
    return front, dominated_by


def band_curve(order):
    """Exact positions in a ranked list where 15-20 bp content starts."""
    flags = [r["span_subband"] == "15_20" for r in order]
    first = next((i + 1 for i, f in enumerate(flags) if f), None)
    last_2130 = max((i + 1 for i, r in enumerate(order) if r["span_subband"] == "21_30"), default=None)
    cum, crossings = 0, {}
    for i, f in enumerate(flags):
        cum += f
        share = cum / (i + 1)
        for lvl in BAND_SHARE_LEVELS:
            if str(lvl) not in crossings and share >= lvl:
                crossings[str(lvl)] = i + 1
    return {"first_15_20_rank": first, "last_21_30_rank": last_2130,
            "n_21_30_in_list": sum(1 for r in order if r["span_subband"] == "21_30"),
            "first_rank_where_cumulative_15_20_share_reaches": crossings}


def write_tsv(path, header, rows):
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--stage4-manifest", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--expect-stage4-manifest-sha256", default=None)
    ap.add_argument("--sizes", type=int, nargs="*", default=SCENARIO_SIZES)
    ap.add_argument("--hybrid-cores", type=int, nargs="*", default=HYBRID_CORES)
    ap.add_argument("--hybrid-totals", type=int, nargs="*", default=HYBRID_TOTALS)
    args = ap.parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)

    m4, m3, hashes = verify_upstream(args.stage4_manifest, args.expect_stage4_manifest_sha256)
    rows = load_union(m3["outputs"]["ranked_union"]["path"])

    a_q = ordered(rows, "_inA", "rank_viewA_no_shared_key")
    a_s = ordered(rows, "_inA", "rank_viewA")
    b_q = ordered(rows, "_inB", "rank_viewB_no_shared_key")
    b_s = ordered(rows, "_inB", "rank_viewB")
    n_a, n_b = len(a_q), len(b_q)
    orders = {"viewA": (a_q, a_s, n_a), "viewB": (b_q, b_s, n_b)}

    scenarios = {}
    for view, (q, s, n_av) in orders.items():
        scenarios[view] = {}
        for size in args.sizes:
            qs, ss = q[:size], s[:size]
            scenarios[view][str(size)] = {
                "n_available_strong_set": n_av,
                "quality_first": stats_for(qs, n_av), "shared_first": stats_for(ss, n_av),
                "loci_in_both_strategies": len({r["marker_id"] for r in qs} & {r["marker_id"] for r in ss}),
            }

    core_sorted = core_order(rows)
    for k in (max(args.hybrid_cores),):
        if [r["marker_id"] for r in core_sorted[:k]] != [r["marker_id"] for r in b_s[:k]]:
            raise SystemExit("core order differs from the Stage-3 shared-first top-k (unexpected)")

    configs = []
    for size in args.sizes:
        for tag, order in (("QF", b_q), ("SF", b_s)):
            members = order[:size]
            st = stats_for(members, n_b)
            configs.append({"id": f"{tag}_{size}", "strategy": "quality_first" if tag == "QF" else "shared_first",
                            "size": size, "core": None, "shared_fraction": st["shared_fraction"],
                            "fraction_21_30": st["fraction_21_30"],
                            "quality_percentile": quality_percentile(members, n_b)})

    qf_sets = {size: {r["marker_id"] for r in b_q[:size]} for size in args.hybrid_totals}
    sf_sets = {size: {r["marker_id"] for r in b_s[:size]} for size in args.hybrid_totals}
    hybrid = {}
    for k in args.hybrid_cores:
        hybrid[f"core_{k}"] = {}
        for total in args.hybrid_totals:
            if total <= k:
                continue
            core, ext, members = hybrid_members(core_sorted, b_q, k, total)
            st = stats_for(members, n_b)
            ids = {r["marker_id"] for r in members}
            shared_members = [r for r in members if r["_shared"]]
            entry = {
                "common_core_size": k, "wgs_extension_size": total - k, "total_wgs_size": total,
                "total_unique_loci": len(ids),
                "stats": st,
                "shared_loci_in_panel": len(shared_members),
                "shared_loci_in_extension_beyond_core": sum(1 for r in ext if r["_shared"]),
                "wgs_only_loci_in_panel": sum(1 for r in members if not r["_shared"]),
                "shared_subset_span_subband": {b: sum(1 for r in shared_members if r["span_subband"] == b) for b in SUBBANDS},
                "core_span_subband": {b: sum(1 for r in core if r["span_subband"] == b) for b in SUBBANDS},
                "extension_span_subband": {b: sum(1 for r in ext if r["span_subband"] == b) for b in SUBBANDS},
                "overlap_with_quality_first_top_total": len(ids & qf_sets[total]),
                "overlap_with_shared_first_top_total": len(ids & sf_sets[total]),
                "mean_quality_rank_percentile_viewB": quality_percentile(members, n_b),
            }
            hybrid[f"core_{k}"][f"total_{total}"] = entry
            configs.append({"id": f"HY_c{k}_T{total}", "strategy": "hybrid", "size": total, "core": k,
                            "shared_fraction": st["shared_fraction"], "fraction_21_30": st["fraction_21_30"],
                            "quality_percentile": entry["mean_quality_rank_percentile_viewB"]})

    front_benefit, dom_benefit = pareto_front(configs, "size_as_benefit")
    front_cost, dom_cost = pareto_front(configs, "size_as_cost")
    within = {}
    for size in sorted({c["size"] for c in configs}):
        group = [c for c in configs if c["size"] == size]
        f, d = pareto_front(group, "fixed_size")
        within[str(size)] = {"front": f, "dominated_by": d}

    band_analysis = {"viewA_shared_first": band_curve(a_s), "viewB_shared_first": band_curve(b_s),
                     "viewB_quality_first": band_curve(b_q)}
    overlaps_b = {str(size): scenarios["viewB"][str(size)]["loci_in_both_strategies"] for size in args.sizes}

    # summary tables (TSV; no membership lists)
    stat_cols = ["n", "pct_of_available_strong_set", "shared", "A_only", "B_only", "shared_fraction",
                 "fraction_15_20", "fraction_21_30"]

    def stat_row(view, size, label, s, overlap):
        return [view, size, label, *[s[c] for c in stat_cols], *[s["span_subband"][b] for b in SUBBANDS],
                s["umap_status"].get("NOT_ASSESSED", 0), s["rmsk_self_category"].get("microsatellite_class_only", 0),
                s["rmsk_self_category"].get("no_rmsk_hit", 0), s["rmsk_self_category"].get("mixed_microsatellite_and_other", 0),
                s["flank_other_repeat_flag"].get("True", 0), s["chrY"], overlap]
    head = ["view", "size", "ranking", *stat_cols, "sub_15_20", "sub_21_30", "sub_31_40", "umap_NOT_ASSESSED",
            "rmsk_microsatellite_class_only", "rmsk_no_rmsk_hit", "rmsk_mixed", "flank_flag_True", "chrY",
            "loci_in_both_strategies"]
    paths = {}
    for label, key in (("quality_first", "quality_first"), ("shared_first", "shared_first")):
        p = os.path.join(args.out_dir, f"tradeoff_stage5.v1.{label}_scenarios.tsv")
        write_tsv(p, head, [stat_row(v, sz, label, sc[key], sc["loci_in_both_strategies"])
                            for v in scenarios for sz, sc in scenarios[v].items()])
        paths[f"{label}_scenarios"] = p
    p = os.path.join(args.out_dir, "tradeoff_stage5.v1.hybrid.tsv")
    write_tsv(p, ["core", "extension", "total", "unique", "shared_in_panel", "shared_beyond_core", "wgs_only",
                  "shared_fraction", "wgs_only_fraction", "sub_15_20", "sub_21_30", "sub_31_40", "chrY",
                  "quality_percentile", "overlap_quality_first", "overlap_shared_first"],
              [[e["common_core_size"], e["wgs_extension_size"], e["total_wgs_size"], e["total_unique_loci"],
                e["shared_loci_in_panel"], e["shared_loci_in_extension_beyond_core"], e["wgs_only_loci_in_panel"],
                e["stats"]["shared_fraction"], e["stats"]["wgs_only_or_view_only_fraction"],
                *[e["stats"]["span_subband"][b] for b in SUBBANDS], e["stats"]["chrY"],
                e["mean_quality_rank_percentile_viewB"], e["overlap_with_quality_first_top_total"],
                e["overlap_with_shared_first_top_total"]]
               for ce in hybrid.values() for e in ce.values()])
    paths["hybrid"] = p
    p = os.path.join(args.out_dir, "tradeoff_stage5.v1.pareto.tsv")
    write_tsv(p, ["id", "strategy", "size", "core", "shared_fraction", "fraction_21_30", "quality_percentile",
                  "on_front_size_as_benefit", "on_front_size_as_cost", "on_front_within_size"],
              [[c["id"], c["strategy"], c["size"], c["core"] if c["core"] else "", c["shared_fraction"],
                c["fraction_21_30"], c["quality_percentile"], c["id"] in front_benefit, c["id"] in front_cost,
                c["id"] in within[str(c["size"])]["front"]] for c in configs])
    paths["pareto"] = p

    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "stage": "trade-off / Pareto analysis (development / RUO)",
        "status": "DEVELOPMENT_ANALYSIS_NOT_A_SELECTED_PANEL",
        "final_panel_size_fixed": False, "panel_size_selected": False, "membership_files_written": False,
        "script": {"version": SCRIPT_VERSION, "path": os.path.abspath(__file__),
                   "sha256": sha256_file(os.path.abspath(__file__)), "python": platform.python_version()},
        "stage4_manifest_path": args.stage4_manifest, "verified_upstream_hashes": hashes,
        "parameters": {"sizes": args.sizes, "hybrid_cores": args.hybrid_cores, "hybrid_totals": args.hybrid_totals,
                       "band_share_reference_levels_not_thresholds": BAND_SHARE_LEVELS,
                       "anchor_bp_development_preference": m3["parameters"]["promega_span_anchor_bp"],
                       "n_viewA_strong": n_a, "n_viewB_strong": n_b, "n_shared_strong": len(core_sorted),
                       "quality_percentile_definition": "mean(rank_viewB_no_shared_key over members) / N_viewB; lower = higher ranking priority",
                       "pareto_objectives": ["shared_fraction (max)", "fraction_21_30 (max)",
                                             "quality_percentile (min)", "size: both conventions reported"]},
        "definitions": {
            "quality_first": "stage-3 rank_*_no_shared_key, unchanged",
            "shared_first": "stage-3 rank_*, unchanged",
            "hybrid": "core = first k shared loci by quality rank (shared status is a membership filter only); "
                      "extension = next (T-k) loci of the View-B quality ranking not already in the core; no duplicates",
            "core_equals_shared_first_prefix_asserted": True},
        "scenarios": scenarios, "viewB_strategy_overlaps": overlaps_b, "hybrid": hybrid,
        "pareto": {"configs": configs,
                   "front_size_as_benefit": front_benefit, "dominated_by_size_as_benefit": dom_benefit,
                   "front_size_as_cost": front_cost, "dominated_by_size_as_cost": dom_cost,
                   "within_fixed_size": within},
        "shared_first_band_analysis": band_analysis,
        "hcc1395_results_used": False, "classification_performed": False, "production_changes": False,
        "upstream_outputs_modified": False,
        "outputs": {k: {"path": v, "sha256": sha256_file(v)} for k, v in paths.items()},
    }
    manifest_path = os.path.join(args.out_dir, "tradeoff_stage5.v1.manifest.json")
    with open(manifest_path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print(f"wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
