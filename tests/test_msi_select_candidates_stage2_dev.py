"""
SYNTHETIC, DEVELOPMENT-ONLY tests for dev/msi_redesign/select_candidates_stage2.py
and its independent verifier verify_candidate_selection_stage2.py.

*** Fabricated loci only. Real results come from a separate run against ***
*** the verified 250 bp candidate panels, reported with manifest/hashes. ***

Run with:
    python -m unittest tests.test_msi_select_candidates_stage2_dev -v
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

import select_candidates_stage2 as s2  # noqa: E402
import verify_candidate_selection_stage2 as ver  # noqa: E402

PANEL_HEADER = ["marker_id", "chrom", "start", "end", "kind", "motif", "ref_repeat_length", "span",
                "umap_multiread_mean", "rmsk_self_category", "flank_other_repeat_flag", "umap_tier",
                "rmsk_tier", "kind_tier", "length_tier", "extraction_span_tier", "cluster_id",
                "cluster_size", "status", "representative_id", "gap_to_representative_bp"]


def prow(mid, chrom="chr1", start=100, span=20, kind="mononucleotide", motif="A", mean="1.0",
         status="SELECTED", rep="", cluster="chr1:C0", size=1):
    return [mid, chrom, start, start + span, kind, motif, span, span, mean, "no_rmsk_hit", "False",
            0, 1, 0, 0, 0, cluster, size, status, rep, ""]


def write_rows(path, rows):
    with open(path, "w", newline="\n") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(PANEL_HEADER)
        w.writerows(rows)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class TestPureRules(unittest.TestCase):
    def test_span_band_boundaries(self):
        got = {s: s2.span_band(s) for s in (5, 9, 10, 14, 15, 40, 41, 60, 61, 90)}
        self.assertEqual(got, {5: "lt10", 9: "lt10", 10: "10_14", 14: "10_14", 15: "15_40", 40: "15_40",
                                41: "41_60", 60: "41_60", 61: "gt60", 90: "gt60"})

    def test_strong_subbands(self):
        self.assertEqual([s2.strong_subband(s) for s in (14, 15, 20, 21, 30, 31, 40, 41)],
                          [None, "15_20", "15_20", "21_30", "21_30", "31_40", "31_40", None])

    def test_geometry_constants_are_explicit(self):
        self.assertEqual(s2.READ_LENGTH_BP, 100)
        self.assertEqual(s2.MIN_FLANK_BP, 5)
        self.assertEqual(s2.MAX_MEASURABLE_SPAN_BP, 90)

    def test_read_measurable_fraction_formula(self):
        self.assertAlmostEqual(s2.read_measurable_fraction(10), 81 / 109)
        self.assertAlmostEqual(s2.read_measurable_fraction(90), 1 / 189)
        self.assertEqual(s2.read_measurable_fraction(91), 0.0)
        vals = [s2.read_measurable_fraction(s) for s in range(5, 91)]
        self.assertEqual(vals, sorted(vals, reverse=True))          # monotone decreasing in span

    def test_umap_status(self):
        self.assertEqual(s2.umap_status("1.0"), "TIER0_MEAN_1.0")
        self.assertEqual(s2.umap_status("0.99"), "BELOW_TOP_TIER")
        self.assertEqual(s2.umap_status(""), "NOT_ASSESSED")

    def test_view_class(self):
        self.assertEqual(s2.view_class(True, True), "A_and_B")
        self.assertEqual(s2.view_class(True, False), "A_only")
        self.assertEqual(s2.view_class(False, True), "B_only")

    def test_evaluate_locus_paths(self):
        ev = s2.evaluate_locus
        self.assertEqual(ev("mononucleotide", 20, "1.0")[4:], ("STRONG_MONO_PRIMARY", "MONO:PASS_ALL_CRITERIA"))
        self.assertEqual(ev("mononucleotide", 20, "0.9")[4:], ("NOT_STRONG_MONO", "MONO:FAIL_S4_UMAP_BELOW_TOP_TIER"))
        self.assertEqual(ev("mononucleotide", 12, "1.0")[4:], ("NOT_STRONG_MONO", "MONO:FAIL_S3_SPAN_BAND_10_14"))
        self.assertEqual(ev("mononucleotide", 8, "1.0")[5], "MONO:FAIL_S3_SPAN_BAND_lt10")
        self.assertEqual(ev("mononucleotide", 45, "1.0")[5], "MONO:FAIL_S3_SPAN_BAND_41_60")
        self.assertEqual(ev("mononucleotide", 70, "1.0")[5], "MONO:FAIL_S3_SPAN_BAND_gt60")
        self.assertEqual(ev("mononucleotide", 95, "1.0")[5], "MONO:FAIL_S1_SPAN_GT_90")
        self.assertEqual(ev("str", 20, "1.0")[4:], ("STRONG_STR_COMPARISON", "STR_COMPARISON:PASS_ALL_CRITERIA"))
        self.assertEqual(ev("str", 12, "1.0")[4], "NOT_STRONG_STR")

    def test_str_is_retained_in_comparison_track_not_dropped(self):
        s1, track, s3, s4, final, reason = s2.evaluate_locus("str", 20, "1.0")
        self.assertEqual(track, "STR_COMPARISON")
        self.assertTrue(s1 and s3 and s4)

    def test_missing_umap_does_not_fail_stage_so_chry_is_not_removed_implicitly(self):
        self.assertEqual(s2.evaluate_locus("mononucleotide", 20, "")[4], "STRONG_MONO_PRIMARY")

    def test_rules_are_deterministic(self):
        self.assertEqual(s2.evaluate_locus("mononucleotide", 33, "1.0"), s2.evaluate_locus("mononucleotide", 33, "1.0"))


class TestEndToEnd(unittest.TestCase):
    def _build(self, tmp):
        panel_a = tmp / "candidate_msi_panel.v1.viewA_wes_compatible.spacing250.tsv"
        trace_a = tmp / "candidate_msi_panel.v1.viewA_wes_compatible.spacing250.traceability.tsv"
        panel_b = tmp / "candidate_msi_panel.v1.viewB_wgs.spacing250.tsv"
        trace_b = tmp / "candidate_msi_panel.v1.viewB_wgs.spacing250.traceability.tsv"
        m1 = prow("m1", span=20)                               # strong mono, shared
        m2 = prow("m2", start=1000, span=20, kind="str", motif="AC")       # A only, STR strong comparison
        m3 = prow("m3", start=2000, span=9)                    # A only, too short
        m4 = prow("m4", chrom="chrY", start=500, span=25, mean="")         # B only, chrY no Umap data
        m5 = prow("m5", start=3000, span=22, mean="0.8")       # B only, below top Umap tier
        write_rows(panel_a, [m1, m2, m3])
        write_rows(panel_b, [m1, m4, m5])
        excluded = prow("x1", start=150, status="EXCLUDED_BY_SPACING", rep="m1")
        write_rows(trace_a, [m1, m2, m3, excluded])
        write_rows(trace_b, [m1, m4, m5, excluded])
        prior = {"outputs": {
            "viewA_wes_compatible": {"selected_panel": {"sha256": sha(panel_a)}, "traceability_all_loci": {"sha256": sha(trace_a)}},
            "viewB_wgs": {"selected_panel": {"sha256": sha(panel_b)}, "traceability_all_loci": {"sha256": sha(trace_b)}}}}
        for view in ("viewA_wes_compatible", "viewB_wgs"):
            for key, path in (("selected_panel", panel_a if view.startswith("viewA") else panel_b),
                              ("traceability_all_loci", trace_a if view.startswith("viewA") else trace_b)):
                prior["outputs"][view][key]["path"] = str(path)
        prior_path = tmp / "prior.json"
        prior_path.write_text(json.dumps(prior))
        return panel_a, panel_b, trace_a, trace_b, prior_path

    def _run(self, tmp):
        panel_a, panel_b, trace_a, trace_b, prior_path = self._build(tmp)
        out = tmp / "out"
        rc = s2.main(["--panel-a", str(panel_a), "--panel-b", str(panel_b), "--trace-a", str(trace_a),
                      "--trace-b", str(trace_b), "--prior-manifest", str(prior_path), "--out-dir", str(out)])
        self.assertEqual(rc, 0)
        return out, (panel_a, panel_b, trace_a, trace_b, prior_path)

    def test_counts_classification_and_chry_reporting(self):
        with tempfile.TemporaryDirectory() as t:
            out, _ = self._run(Path(t))
            m = json.loads((out / "candidate_selection_stage2.v1.manifest.json").read_text())
            self.assertEqual(m["view_relationship"], {"A_and_B": 1, "A_only": 2, "B_only": 2})
            self.assertEqual(m["n_union_representatives"], 5)
            a, b, sh = m["summary"]["viewA"], m["summary"]["viewB"], m["summary"]["shared_A_and_B"]
            key = "S4_pass_umap_top_tier_or_not_assessed_=_STRONG"
            self.assertEqual(a["stage_attrition"]["MONO_PRIMARY"][key], 1)         # m1
            self.assertEqual(a["stage_attrition"]["STR_COMPARISON"][key], 1)       # m2
            self.assertEqual(b["stage_attrition"]["MONO_PRIMARY"][key], 2)         # m1 + chrY m4 (not removed)
            self.assertEqual(b["chrY_impact"]["strong_mono_chrY"], 1)
            self.assertEqual(b["chrY_impact"]["strong_mono_without_chrY"], 1)
            self.assertEqual(sh["stage_attrition"]["S0_input_250bp_representatives"], 1)
            self.assertEqual(m["summary"]["B_only"]["umap_not_assessed_by_chromosome"], {"chrY": 1})
            self.assertFalse(m["final_panel_size_fixed"])
            self.assertFalse(m["hcc1395_results_used"])
            rows = list(csv.DictReader(open(out / "candidate_selection_stage2.v1.traceability.tsv"), delimiter="\t"))
            by_id = {r["marker_id"]: r for r in rows}
            self.assertEqual(len(rows), 5)                                        # nothing dropped
            self.assertEqual(by_id["m1"]["source_view_class"], "A_and_B")
            self.assertEqual(by_id["m1"]["n_spacing_excluded_represented_viewA"], "1")
            self.assertEqual(by_id["m5"]["reason"], "MONO:FAIL_S4_UMAP_BELOW_TOP_TIER")
            self.assertEqual(by_id["m3"]["reason"], "MONO:FAIL_S3_SPAN_BAND_lt10")
            self.assertEqual(by_id["m4"]["final_class"], "STRONG_MONO_PRIMARY")

    def test_independent_verifier_passes_and_detects_tampering(self):
        with tempfile.TemporaryDirectory() as t:
            out, (panel_a, panel_b, trace_a, trace_b, prior_path) = self._run(Path(t))
            manifest = str(out / "candidate_selection_stage2.v1.manifest.json")
            m = json.loads(Path(manifest).read_text())
            m["prior_manifest"]["path"] = str(prior_path)
            self.assertEqual(ver.verify(manifest), [])
            trace = out / "candidate_selection_stage2.v1.traceability.tsv"
            original = trace.read_text()
            trace.write_text(original.replace("STRONG_MONO_PRIMARY", "NOT_STRONG_MONO", 1))
            failures = ver.verify(manifest)
            self.assertTrue(any("hash mismatch" in f or "disagree" in f for f in failures))

    def test_wrong_input_hash_aborts(self):
        with tempfile.TemporaryDirectory() as t:
            tmp = Path(t)
            panel_a, panel_b, trace_a, trace_b, prior_path = self._build(tmp)
            panel_a.write_text(panel_a.read_text() + "\n")           # alters a 250 bp input
            with self.assertRaises(SystemExit):
                s2.main(["--panel-a", str(panel_a), "--panel-b", str(panel_b), "--trace-a", str(trace_a),
                         "--trace-b", str(trace_b), "--prior-manifest", str(prior_path),
                         "--out-dir", str(tmp / "out")])


if __name__ == "__main__":
    unittest.main()
