"""
SYNTHETIC, DEVELOPMENT-ONLY tests for dev/msi_redesign/scenario_preview_stage4.py
and its independent verifier verify_scenario_preview_stage4.py.

*** Fabricated loci and manifests only. Real scenario previews come from a ***
*** separate run over the hash-verified Stage-3 rankings.                  ***

Run with:
    python -m unittest tests.test_msi_scenario_preview_stage4_dev -v
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

import scenario_preview_stage4 as s4  # noqa: E402
import verify_scenario_preview_stage4 as v4  # noqa: E402

UNION_HEADER = ["marker_id", "source_view_class", "in_viewA_strong", "in_viewB_strong", "chrom", "span",
                "span_subband", "umap_status", "umap_singleread_window_coverage", "umap_multiread_min",
                "rmsk_self_category", "rmsk_tier", "flank_other_repeat_flag", "read_measurable_fraction",
                "K1_shared_tier", "K7_tie_hash", "rank_viewA", "rank_viewB",
                "rank_viewA_no_shared_key", "rank_viewB_no_shared_key"]

# id, class, span, chrom, umap_assessed
LOCI = [("p1", "B_only", 25, "chr2", True), ("p2", "A_and_B", 24, "chr1", True),
        ("p3", "A_and_B", 26, "chr3", True), ("p4", "A_and_B", 21, "chrY", False),
        ("p5", "A_only", 22, "chr4", True), ("p6", "B_only", 30, "chr5", True),
        ("p7", "A_and_B", 35, "chr6", True), ("p8", "B_only", 16, "chr7", True)]


def geom(span):
    return round((100 - 10 - span + 1) / (100 + span - 1), 4)


def subband(span):
    return "15_20" if span <= 20 else "21_30" if span <= 30 else "31_40"


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def key(loc, with_shared):
    mid, cls, span, chrom, assessed = loc
    k = [0 if cls == "A_and_B" else 1] if with_shared else []
    return tuple(k + [abs(span - 25), 0 if assessed else 1, -geom(span), 1, 0, (0, -1.0) if assessed else (0, -1.0),
                      (0, -1.0) if assessed else (1, 0.0), mid])


def build_union_rows():
    in_a = [l for l in LOCI if l[1] in ("A_and_B", "A_only")]
    in_b = [l for l in LOCI if l[1] in ("A_and_B", "B_only")]
    ranks = {}
    for view, subset in (("A", in_a), ("B", in_b)):
        for tag, ws in (("", True), ("_no_shared_key", False)):
            for i, l in enumerate(sorted(subset, key=lambda x: key(x, ws))):
                ranks[(view, tag, l[0])] = i + 1
    rows = []
    for l in LOCI:
        mid, cls, span, chrom, assessed = l
        rows.append([mid, cls, cls in ("A_and_B", "A_only"), cls in ("A_and_B", "B_only"), chrom, span,
                     subband(span), "TIER0_MEAN_1.0" if assessed else "NOT_ASSESSED", 1.0,
                     1.0 if assessed else "", "no_rmsk_hit", 1, False, geom(span),
                     0 if cls == "A_and_B" else 1, mid,
                     ranks.get(("A", "", mid), ""), ranks.get(("B", "", mid), ""),
                     ranks.get(("A", "_no_shared_key", mid), ""), ranks.get(("B", "_no_shared_key", mid), "")])
    return rows


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
        w.writerow(UNION_HEADER)
        w.writerows(build_union_rows())
    rank_a, rank_b = tiny("ranked_viewA.tsv"), tiny("ranked_viewB.tsv")
    f250 = {v: {k: entry(tiny(f"{v}_{k}.tsv")) for k in ("selected_panel", "traceability_all_loci")}
            for v in ("viewA_wes_compatible", "viewB_wgs")}
    m250 = tmp / "m250.json"
    m250.write_text(json.dumps({"outputs": f250}))
    m2 = tmp / "m2.json"
    m2.write_text(json.dumps({"outputs": {"strong_mono_primary": entry(tiny("strong.tsv"))},
                              "prior_manifest": entry(m250)}))
    m3 = tmp / "m3.json"
    m3.write_text(json.dumps({
        "outputs": {"ranked_union": entry(union), "ranked_viewA": entry(rank_a), "ranked_viewB": entry(rank_b)},
        "stage2_manifest_path": str(m2), "inputs_verified_hashes": {"stage2_manifest": sha(m2)},
        "parameters": {"promega_span_anchor_bp": 25},
        "hcc1395_results_used": False, "classification_performed": False, "production_changes": False,
        "final_panel_size_fixed": False, "panel_size_selected": False, "stage2_outputs_modified": False}))
    return m3, union


class TestScenarioStats(unittest.TestCase):
    def test_top_n_orders_by_rank_and_ignores_rows_not_in_view(self):
        rows = [{"rank": "2", "id": "b"}, {"rank": "", "id": "x"}, {"rank": "1", "id": "a"}]
        top, available = s4.top_n(rows, "rank", 5)
        self.assertEqual([r["id"] for r in top], ["a", "b"])
        self.assertEqual(available, 2)

    def test_scenario_stats_fields(self):
        top = [{"source_view_class": "A_and_B", "span_subband": "21_30", "umap_status": "TIER0_MEAN_1.0",
                "rmsk_self_category": "no_rmsk_hit", "flank_other_repeat_flag": "False", "chrom": "chr1"},
               {"source_view_class": "B_only", "span_subband": "15_20", "umap_status": "NOT_ASSESSED",
                "rmsk_self_category": "microsatellite_class_only", "flank_other_repeat_flag": "True", "chrom": "chrY"}]
        s = s4.scenario_stats(top, 8)
        self.assertEqual((s["n"], s["shared"], s["A_only"], s["B_only"], s["chrY"]), (2, 1, 0, 1, 1))
        self.assertEqual(s["pct_of_available_strong_set"], 25.0)
        self.assertEqual(s["span_subband"], {"15_20": 1, "21_30": 1, "31_40": 0})
        self.assertEqual(s["umap_status"], {"TIER0_MEAN_1.0": 1, "NOT_ASSESSED": 1})

    def test_labels_are_the_two_distinct_stage3_columns(self):
        self.assertEqual(s4.VIEWS["viewB"][1], "rank_viewB_no_shared_key")   # primary quality = shared-free
        self.assertEqual(s4.VIEWS["viewB"][2], "rank_viewB")                 # secondary = shared-first
        self.assertEqual(s4.DEFAULT_SIZES, [500, 1000, 1500, 2000, 5000])


class TestEndToEnd(unittest.TestCase):
    def _run(self, tmp, m3):
        out = tmp / "out"
        rc = s4.main(["--stage3-manifest", str(m3), "--out-dir", str(out), "--sizes", "2", "3", "10"])
        self.assertEqual(rc, 0)
        return out

    def test_previews_show_policy_difference_and_no_membership_files(self):
        with tempfile.TemporaryDirectory() as t:
            tmp = Path(t)
            m3, _ = build_chain(tmp)
            out = self._run(tmp, m3)
            m = json.loads((out / "scenario_preview_stage4.v1.manifest.json").read_text())
            b3 = m["scenarios"]["viewB"]["3"]
            self.assertEqual(b3["n_available_strong_set"], 7)
            prim, sec = b3["primary_quality_ranking"], b3["secondary_common_core_ranking"]
            self.assertEqual((prim["shared"], prim["B_only"]), (2, 1))       # p1 (B-only, perfect) enters
            self.assertEqual((sec["shared"], sec["B_only"]), (3, 0))         # shared-first keeps view-only out
            self.assertEqual(sec["chrY"], 1)                                 # chrY p4 present, never removed
            self.assertEqual(b3["loci_in_both_rankings_top_n"], 2)
            a3 = m["scenarios"]["viewA"]["3"]
            self.assertEqual(a3["primary_quality_ranking"]["A_only"], 1)     # A-only p5 outranks shared p4/p7
            self.assertEqual(a3["secondary_common_core_ranking"]["A_only"], 0)
            big = m["scenarios"]["viewB"]["10"]["primary_quality_ranking"]   # size above available -> capped, not padded
            self.assertEqual(big["n"], 7)
            self.assertEqual(big["pct_of_available_strong_set"], 100.0)
            self.assertFalse(m["panel_size_selected"])
            self.assertTrue(m["ranking_label_mapping"]["not_combined_into_one_score"])
            self.assertEqual(sorted(p.name for p in out.iterdir()),
                             ["scenario_preview_stage4.v1.comparison.tsv", "scenario_preview_stage4.v1.manifest.json"])

    def test_independent_verifier_passes_detects_tampering_and_policy_swap(self):
        with tempfile.TemporaryDirectory() as t:
            tmp = Path(t)
            m3, union = build_chain(tmp)
            out = self._run(tmp, m3)
            manifest = str(out / "scenario_preview_stage4.v1.manifest.json")
            self.assertEqual(v4.verify(manifest), [])
            comp = out / "scenario_preview_stage4.v1.comparison.tsv"
            comp.write_text(comp.read_text() + "x")
            self.assertTrue(any("hash mismatch" in f for f in v4.verify(manifest)))

    def test_verifier_flags_a_primary_ranking_that_is_secretly_shared_first(self):
        with tempfile.TemporaryDirectory() as t:
            tmp = Path(t)
            m3, union = build_chain(tmp)
            # corrupt the Stage-3 union so the "primary" columns equal the shared-first ones
            rows = list(csv.DictReader(open(union), delimiter="\t"))
            for r in rows:
                r["rank_viewB_no_shared_key"] = r["rank_viewB"]
                r["rank_viewA_no_shared_key"] = r["rank_viewA"]
            with open(union, "w", newline="\n") as fh:
                w = csv.DictWriter(fh, fieldnames=UNION_HEADER, delimiter="\t", lineterminator="\n")
                w.writeheader()
                w.writerows(rows)
            # refresh recorded hashes so only the ordering logic can catch it
            meta = json.loads(m3.read_text())
            meta["outputs"]["ranked_union"]["sha256"] = sha(union)
            m3.write_text(json.dumps(meta))
            out = tmp / "out"
            s4.main(["--stage3-manifest", str(m3), "--out-dir", str(out), "--sizes", "3"])
            fails = v4.verify(str(out / "scenario_preview_stage4.v1.manifest.json"))
            self.assertTrue(any("primary quality ranking" in f for f in fails))

    def test_refuses_on_tampered_stage3_output_wrong_hash_or_unsafe_flag(self):
        with tempfile.TemporaryDirectory() as t:
            tmp = Path(t)
            m3, union = build_chain(tmp)
            with self.assertRaises(SystemExit):
                s4.main(["--stage3-manifest", str(m3), "--out-dir", str(tmp / "o"),
                         "--expect-stage3-manifest-sha256", "0" * 64])
            meta = json.loads(m3.read_text())
            meta["panel_size_selected"] = True
            m3.write_text(json.dumps(meta))
            with self.assertRaises(SystemExit):
                s4.main(["--stage3-manifest", str(m3), "--out-dir", str(tmp / "o")])
            meta["panel_size_selected"] = False
            m3.write_text(json.dumps(meta))
            union.write_text(union.read_text() + "\n")
            with self.assertRaises(SystemExit):
                s4.main(["--stage3-manifest", str(m3), "--out-dir", str(tmp / "o")])


if __name__ == "__main__":
    unittest.main()
