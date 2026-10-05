"""
SYNTHETIC, DEVELOPMENT-ONLY tests for dev/msi_redesign/tradeoff_pareto_stage5.py
and its independent verifier verify_tradeoff_stage5.py.

*** Fabricated loci and manifests only. Real results come from a separate ***
*** run over the hash-verified Stage-3 rankings.                          ***

Run with:
    python -m unittest tests.test_msi_tradeoff_pareto_stage5_dev -v
"""

import csv
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path

DEV_DIR = Path(__file__).resolve().parents[1] / "dev" / "msi_redesign"
sys.path.insert(0, str(DEV_DIR))

import tradeoff_pareto_stage5 as t5  # noqa: E402
import verify_tradeoff_stage5 as v5  # noqa: E402

HEADER = ["marker_id", "source_view_class", "in_viewA_strong", "in_viewB_strong", "chrom", "span", "span_subband",
          "umap_status", "rmsk_self_category", "flank_other_repeat_flag", "K2_promega_span_distance",
          "K4_read_measurable_fraction", "K5_rmsk_tier", "K5_flank_flag", "K6_singleread_cov", "K6_multiread_min",
          "K7_tie_hash", "rank_viewA", "rank_viewB", "rank_viewA_no_shared_key", "rank_viewB_no_shared_key"]


def geom(span):
    return round((100 - 10 - span + 1) / (100 + span - 1), 4)


def subband(span):
    return "15_20" if span <= 20 else "21_30" if span <= 30 else "31_40"


def make_loci():
    """8 shared, 4 A-only, 12 B-only with varied spans; one chrY not-assessed shared locus."""
    spans_shared = [25, 24, 26, 22, 28, 18, 33, 21]
    spans_a = [23, 27, 17, 36]
    spans_b = [25, 25, 24, 30, 20, 19, 26, 23, 29, 16, 22, 35]
    loci = []
    for i, s in enumerate(spans_shared):
        loci.append(("s%02d" % i, "A_and_B", s, "chrY" if i == 7 else "chr%d" % (1 + i % 5), i != 7))
    for i, s in enumerate(spans_a):
        loci.append(("a%02d" % i, "A_only", s, "chr%d" % (1 + i), True))
    for i, s in enumerate(spans_b):
        loci.append(("b%02d" % i, "B_only", s, "chr%d" % (1 + i % 6), True))
    return loci


def qkey(l, with_shared, idx):
    mid, cls, span, chrom, assessed = l
    head = [0 if cls == "A_and_B" else 1] if with_shared else []
    return tuple(head + [abs(span - 25), 0 if assessed else 1, -geom(span), 1, 0, (0, -1.0),
                         (0, -1.0) if assessed else (1, 0.0), "h%03d" % idx])


def build_rows():
    loci = make_loci()
    idx = {l[0]: i for i, l in enumerate(loci)}
    ranks = {}
    for view, classes in (("A", ("A_and_B", "A_only")), ("B", ("A_and_B", "B_only"))):
        subset = [l for l in loci if l[1] in classes]
        for tag, ws in (("", True), ("_nsk", False)):
            for r, l in enumerate(sorted(subset, key=lambda x: qkey(x, ws, idx[x[0]])), start=1):
                ranks[(view, tag, l[0])] = r
    rows = []
    for l in loci:
        mid, cls, span, chrom, assessed = l
        rows.append([mid, cls, cls != "B_only", cls != "A_only", chrom, span, subband(span),
                     "TIER0_MEAN_1.0" if assessed else "NOT_ASSESSED", "no_rmsk_hit", False, abs(span - 25),
                     geom(span), 1, 0, 1.0, 1.0 if assessed else "", "h%03d" % idx[mid],
                     ranks.get(("A", "", mid), ""), ranks.get(("B", "", mid), ""),
                     ranks.get(("A", "_nsk", mid), ""), ranks.get(("B", "_nsk", mid), "")])
    return rows


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def build_chain(tmp):
    def tiny(name):
        p = tmp / name
        p.write_text(name)
        return p

    def entry(p):
        return {"path": str(p), "sha256": sha(p)}

    union = tmp / "ranked_union.tsv"
    with open(union, "w", newline="\n") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(HEADER)
        w.writerows(build_rows())
    f250 = {v: {k: entry(tiny(f"{v}_{k}.tsv")) for k in ("selected_panel", "traceability_all_loci")}
            for v in ("viewA_wes_compatible", "viewB_wgs")}
    m250 = tmp / "m250.json"
    m250.write_text(json.dumps({"outputs": f250}))
    m2 = tmp / "m2.json"
    m2.write_text(json.dumps({"outputs": {"strong_mono_primary": entry(tiny("strong.tsv"))}, "prior_manifest": entry(m250)}))
    m3 = tmp / "m3.json"
    m3.write_text(json.dumps({
        "outputs": {"ranked_union": entry(union), "ranked_viewA": entry(tiny("ra.tsv")), "ranked_viewB": entry(tiny("rb.tsv"))},
        "stage2_manifest_path": str(m2), "inputs_verified_hashes": {"stage2_manifest": sha(m2)},
        "parameters": {"promega_span_anchor_bp": 25},
        "hcc1395_results_used": False, "classification_performed": False, "production_changes": False,
        "final_panel_size_fixed": False, "panel_size_selected": False, "stage2_outputs_modified": False}))
    m4 = tmp / "m4.json"
    m4.write_text(json.dumps({
        "outputs": {"comparison_table": entry(tiny("comparison.tsv"))},
        "stage3_manifest_path": str(m3), "verified_upstream_hashes": {"stage3_manifest": sha(m3)},
        "final_panel_size_fixed": False, "panel_size_selected": False, "hcc1395_results_used": False,
        "classification_performed": False, "production_changes": False, "ranking_algorithm_rewritten": False,
        "stage3_outputs_modified": False, "chromosome_is_a_ranking_key": False, "chromosome_quotas": False,
        "chrY_removed": False}))
    return m4, m3, union


ARGS = ["--sizes", "3", "5", "8", "--hybrid-cores", "2", "3", "--hybrid-totals", "6", "10"]


class TestPure(unittest.TestCase):
    def _rows(self):
        rows = []
        for r in build_rows():
            d = dict(zip(HEADER, r))
            for c in ("rank_viewA", "rank_viewB", "rank_viewA_no_shared_key", "rank_viewB_no_shared_key"):
                d[c] = int(d[c]) if d[c] != "" else None
            d["_shared"] = d["source_view_class"] == "A_and_B"
            d["_inA"], d["_inB"] = d["in_viewA_strong"] is True, d["in_viewB_strong"] is True
            d["chrom"] = d["chrom"]
            rows.append(d)
        return rows

    def test_dominance_and_fronts(self):
        cfgs = [{"id": "x", "size": 10, "shared_fraction": 1.0, "fraction_21_30": 0.5, "quality_percentile": 0.2},
                {"id": "y", "size": 10, "shared_fraction": 0.5, "fraction_21_30": 1.0, "quality_percentile": 0.1},
                {"id": "z", "size": 10, "shared_fraction": 0.4, "fraction_21_30": 0.9, "quality_percentile": 0.3}]
        front, dom = t5.pareto_front(cfgs, "fixed_size")
        self.assertEqual(sorted(front), ["x", "y"])          # z is dominated by y
        self.assertEqual(dom, {"z": "y"})

    def test_size_convention_changes_the_front(self):
        small = {"id": "s", "size": 5, "shared_fraction": 0.5, "fraction_21_30": 0.5, "quality_percentile": 0.5}
        big = {"id": "b", "size": 50, "shared_fraction": 0.5, "fraction_21_30": 0.5, "quality_percentile": 0.5}
        self.assertEqual(t5.pareto_front([small, big], "size_as_benefit")[0], ["b"])
        self.assertEqual(t5.pareto_front([small, big], "size_as_cost")[0], ["s"])
        self.assertEqual(sorted(t5.pareto_front([small, big], "fixed_size")[0]), ["b", "s"])   # identical -> neither dominates

    def test_core_order_is_quality_order_of_shared_loci_and_equals_shared_first_prefix(self):
        rows = self._rows()
        core = t5.core_order(rows)
        self.assertEqual(len(core), 8)
        b_s = t5.ordered(rows, "_inB", "rank_viewB")
        self.assertEqual([r["marker_id"] for r in core], [r["marker_id"] for r in b_s[:8]])

    def test_core_order_rejects_inconsistent_view_orders(self):
        rows = self._rows()
        shared = [r for r in rows if r["_shared"]]
        shared[0]["rank_viewA_no_shared_key"], shared[1]["rank_viewA_no_shared_key"] = \
            shared[1]["rank_viewA_no_shared_key"], shared[0]["rank_viewA_no_shared_key"]
        with self.assertRaises(SystemExit):
            t5.core_order(rows)

    def test_hybrid_is_duplicate_free_core_prefix_plus_quality_extension(self):
        rows = self._rows()
        core_sorted = t5.core_order(rows)
        b_q = t5.ordered(rows, "_inB", "rank_viewB_no_shared_key")
        core, ext, members = t5.hybrid_members(core_sorted, b_q, 3, 10)
        ids = [r["marker_id"] for r in members]
        self.assertEqual(len(ids), 10)
        self.assertEqual(len(set(ids)), 10)
        self.assertEqual(ids[:3], [r["marker_id"] for r in core_sorted[:3]])
        core_ids = {r["marker_id"] for r in core}
        self.assertFalse(core_ids & {r["marker_id"] for r in ext})
        # extension follows the View-B quality order, skipping only core loci
        expected = [r["marker_id"] for r in b_q if r["marker_id"] not in core_ids][:7]
        self.assertEqual([r["marker_id"] for r in ext], expected)
        # extension may legitimately contain shared loci beyond the core
        self.assertTrue(any(r["_shared"] for r in ext))

    def test_hybrid_rejects_impossible_sizes(self):
        rows = self._rows()
        core_sorted = t5.core_order(rows)
        b_q = t5.ordered(rows, "_inB", "rank_viewB_no_shared_key")
        with self.assertRaises(ValueError):
            t5.hybrid_members(core_sorted, b_q, 3, 500)
        with self.assertRaises(ValueError):
            t5.hybrid_members(core_sorted, b_q, 99, 100)

    def test_quality_percentile_arithmetic(self):
        members = [{"rank_viewB_no_shared_key": 1}, {"rank_viewB_no_shared_key": 3}]
        self.assertEqual(t5.quality_percentile(members, 10), 0.2)

    def test_band_curve_positions_and_crossings(self):
        order = [{"span_subband": b} for b in ["21_30", "21_30", "15_20", "21_30", "15_20", "31_40"]]
        c = t5.band_curve(order)
        self.assertEqual(c["first_15_20_rank"], 3)
        self.assertEqual(c["last_21_30_rank"], 4)
        self.assertEqual(c["n_21_30_in_list"], 3)
        cross = c["first_rank_where_cumulative_15_20_share_reaches"]
        self.assertEqual(cross["0.25"], 3)          # 1/3 >= 0.25 at position 3
        self.assertNotIn("0.5", cross)              # the 15-20 share never reaches 50% in this list
        self.assertEqual(cross["0.01"], 3)

    def test_stats_add_fractions(self):
        rows = self._rows()
        b_q = t5.ordered(rows, "_inB", "rank_viewB_no_shared_key")
        st = t5.stats_for(b_q[:4], len(b_q))
        self.assertEqual(st["n"], 4)
        self.assertAlmostEqual(st["fraction_21_30"], st["span_subband"]["21_30"] / 4)
        self.assertAlmostEqual(st["shared_fraction"] + st["wgs_only_or_view_only_fraction"], 1.0)


class TestEndToEnd(unittest.TestCase):
    def _run(self, tmp):
        m4, m3, union = build_chain(tmp)
        out = tmp / "out"
        rc = t5.main(["--stage4-manifest", str(m4), "--out-dir", str(out), *ARGS])
        self.assertEqual(rc, 0)
        return out, m4, m3, union

    def test_outputs_manifest_contents_and_no_membership_files(self):
        with tempfile.TemporaryDirectory() as t:
            out, *_ = self._run(Path(t))
            m = json.loads((out / "tradeoff_stage5.v1.manifest.json").read_text())
            self.assertFalse(m["panel_size_selected"])
            self.assertFalse(m["membership_files_written"])
            self.assertEqual(sorted(p.name for p in out.iterdir()), sorted([
                "tradeoff_stage5.v1.manifest.json", "tradeoff_stage5.v1.quality_first_scenarios.tsv",
                "tradeoff_stage5.v1.shared_first_scenarios.tsv", "tradeoff_stage5.v1.hybrid.tsv",
                "tradeoff_stage5.v1.pareto.tsv"]))
            e = m["hybrid"]["core_3"]["total_10"]
            self.assertEqual((e["common_core_size"], e["wgs_extension_size"], e["total_wgs_size"],
                              e["total_unique_loci"]), (3, 7, 10, 10))
            self.assertEqual(e["stats"]["shared"], e["shared_loci_in_panel"])
            self.assertEqual(e["shared_loci_in_panel"] + e["wgs_only_loci_in_panel"], 10)
            self.assertNotIn("total_6", m["hybrid"].get("core_7", {}))
            ids = {c["id"] for c in m["pareto"]["configs"]}
            self.assertIn("QF_5", ids)
            self.assertIn("SF_5", ids)
            self.assertIn("HY_c2_T6", ids)
            self.assertTrue(set(m["pareto"]["front_size_as_benefit"]) <= ids)
            self.assertIn("5", m["pareto"]["within_fixed_size"])
            self.assertTrue(m["definitions"]["core_equals_shared_first_prefix_asserted"])

    def test_independent_verifier_passes_and_detects_tampering(self):
        with tempfile.TemporaryDirectory() as t:
            out, *_ = self._run(Path(t))
            manifest = out / "tradeoff_stage5.v1.manifest.json"
            self.assertEqual(v5.verify(str(manifest)), [])
            m = json.loads(manifest.read_text())
            m["hybrid"]["core_3"]["total_10"]["stats"]["shared"] += 1
            manifest.write_text(json.dumps(m))
            self.assertTrue(any("differs" in f for f in v5.verify(str(manifest))))

    def test_verifier_detects_pareto_front_tampering(self):
        with tempfile.TemporaryDirectory() as t:
            out, *_ = self._run(Path(t))
            manifest = out / "tradeoff_stage5.v1.manifest.json"
            m = json.loads(manifest.read_text())
            m["pareto"]["front_size_as_cost"] = m["pareto"]["front_size_as_cost"][:-1] or ["QF_3"]
            manifest.write_text(json.dumps(m))
            self.assertTrue(any("Pareto" in f for f in v5.verify(str(manifest))))

    def test_hybrid_overlaps_reported_against_both_strategies(self):
        with tempfile.TemporaryDirectory() as t:
            out, *_ = self._run(Path(t))
            m = json.loads((out / "tradeoff_stage5.v1.manifest.json").read_text())
            for core in m["hybrid"].values():
                for e in core.values():
                    self.assertLessEqual(e["overlap_with_quality_first_top_total"], e["total_wgs_size"])
                    self.assertLessEqual(e["overlap_with_shared_first_top_total"], e["total_wgs_size"])

    def test_refuses_on_wrong_hash_unsafe_flag_or_changed_upstream(self):
        with tempfile.TemporaryDirectory() as t:
            tmp = Path(t)
            m4, m3, union = build_chain(tmp)
            with self.assertRaises(SystemExit):
                t5.main(["--stage4-manifest", str(m4), "--out-dir", str(tmp / "o"),
                         "--expect-stage4-manifest-sha256", "0" * 64, *ARGS])
            meta = json.loads(m4.read_text())
            meta["panel_size_selected"] = True
            m4.write_text(json.dumps(meta))
            with self.assertRaises(SystemExit):
                t5.main(["--stage4-manifest", str(m4), "--out-dir", str(tmp / "o"), *ARGS])
            meta["panel_size_selected"] = False
            m4.write_text(json.dumps(meta))
            union.write_text(union.read_text() + "\n")
            with self.assertRaises(SystemExit):
                t5.main(["--stage4-manifest", str(m4), "--out-dir", str(tmp / "o"), *ARGS])


if __name__ == "__main__":
    unittest.main()
