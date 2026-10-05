"""
SYNTHETIC, DEVELOPMENT-ONLY tests for dev/msi_redesign/rank_candidates_stage3.py
and its independent verifier verify_ranking_stage3.py.

*** Fabricated loci/manifests only. Real results come from a separate run ***
*** against the validated stage-2 strong candidate outputs.               ***

Run with:
    python -m unittest tests.test_msi_rank_candidates_stage3_dev -v
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

import rank_candidates_stage3 as r3  # noqa: E402
import verify_ranking_stage3 as ver3  # noqa: E402

STRONG_HEADER = ["marker_id", "source_view_class", "in_viewA_reps", "in_viewB_reps", "chrom", "start", "end",
                 "span", "motif", "ref_repeat_length", "umap_multiread_mean", "umap_status",
                 "rmsk_self_category", "rmsk_tier", "flank_other_repeat_flag", "read_measurable_fraction", "is_chrY"]


def geom(span):
    return round((100 - 10 - span + 1) / (100 + span - 1), 4)


def srow(mid, span=25, view="A_and_B", chrom="chr1", umap="1.0", rmsk=1, flank=False):
    in_a, in_b = view in ("A_and_B", "A_only"), view in ("A_and_B", "B_only")
    status = "TIER0_MEAN_1.0" if umap != "" else "NOT_ASSESSED"
    return [mid, view, in_a, in_b, chrom, 1000, 1000 + span, span, "A", span, umap, status,
            "no_rmsk_hit", rmsk, flank, geom(span), chrom == "chrY"]


def rec(mid, span=25, shared=True, status="TIER0_MEAN_1.0", rmsk=1, flank=False, sc=1.0, mm=1.0, frac=None):
    return {"marker_id": mid, "span": span, "shared": shared, "umap_status": status,
            "read_measurable_fraction": geom(span) if frac is None else frac, "rmsk_tier": rmsk,
            "flank": flank, "singleread_cov": sc, "multiread_min": mm, "chrom": "chr1", "rmsk_category": "no_rmsk_hit"}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


class TestPureRules(unittest.TestCase):
    def test_span_subbands(self):
        self.assertEqual([r3.span_subband(s) for s in (14, 15, 20, 21, 30, 31, 40, 41)],
                          ["outside_15_40", "15_20", "15_20", "21_30", "21_30", "31_40", "31_40", "outside_15_40"])

    def test_promega_anchor_is_median_of_verified_spans(self):
        self.assertEqual(r3.promega_anchor([25, 27, 21, 23, 26]), 25)

    def test_prefix_size_uses_percentages_and_ceil(self):
        self.assertEqual(r3.prefix_size(6074, 0.10), 608)
        self.assertEqual(r3.prefix_size(10, 0.05), 1)
        self.assertEqual(r3.prefix_size(0, 0.5), 0)
        self.assertEqual(r3.prefix_size(100, 1.0), 100)

    def test_k1_shared_first_but_only_when_enabled(self):
        shared, only = rec("a", span=30, shared=True), rec("b", span=25, shared=False)
        self.assertLess(r3.rank_key(shared, 25, True), r3.rank_key(only, 25, True))
        self.assertLess(r3.rank_key(only, 25, False), r3.rank_key(shared, 25, False))

    def test_k2_distance_to_anchor_ranks_within_15_40(self):
        self.assertLess(r3.rank_key(rec("a", span=25), 25), r3.rank_key(rec("b", span=27), 25))
        self.assertLess(r3.rank_key(rec("a", span=27), 25), r3.rank_key(rec("b", span=18), 25))
        self.assertEqual(r3.rank_key(rec("a", span=24), 25)[1], r3.rank_key(rec("b", span=26), 25)[1])

    def test_k3_not_assessed_stays_identifiable_and_ranks_after_assessed(self):
        assessed, na = rec("a"), rec("b", status="NOT_ASSESSED", sc=1.0, mm=None)
        self.assertLess(r3.rank_key(assessed, 25), r3.rank_key(na, 25))

    def test_k4_higher_read_measurable_fraction_wins_ties_on_distance(self):
        self.assertLess(r3.rank_key(rec("a", span=24), 25), r3.rank_key(rec("b", span=26), 25))

    def test_k5_rmsk_tier_then_flank(self):
        base = dict(span=25)
        self.assertLess(r3.rank_key(rec("a", rmsk=0, **base), 25), r3.rank_key(rec("b", rmsk=1, **base), 25))
        self.assertLess(r3.rank_key(rec("a", flank=False, **base), 25), r3.rank_key(rec("b", flank=True, **base), 25))

    def test_k6_uniqueness_higher_first_and_missing_never_imputed(self):
        hi, lo = rec("a", sc=1.0), rec("b", sc=0.5)
        self.assertLess(r3.rank_key(hi, 25), r3.rank_key(lo, 25))
        present, missing = rec("a", mm=0.9), rec("b", mm=None)
        self.assertLess(r3.rank_key(present, 25), r3.rank_key(missing, 25))

    def test_k7_tie_break_deterministic_and_total(self):
        a, b = rec("id_a"), rec("id_b")
        self.assertNotEqual(r3.rank_key(a, 25), r3.rank_key(b, 25))
        self.assertEqual(r3.rank_key(a, 25), r3.rank_key(rec("id_a"), 25))

    def test_chromosome_is_not_a_ranking_key(self):
        x, y = rec("same"), rec("same")
        y["chrom"] = "chrY"
        self.assertEqual(r3.rank_key(x, 25), r3.rank_key(y, 25))

    def test_rank_records_permutation_and_input_not_reordered_away(self):
        recs = [rec(f"m{i}", span=15 + i) for i in range(10)]
        ranks, ordered = r3.rank_records(recs, 25)
        self.assertEqual(sorted(ranks.values()), list(range(1, 11)))
        self.assertEqual(len(ordered), 10)


class TestEndToEnd(unittest.TestCase):
    def _build(self, tmp):
        rows = [srow("s1", 25), srow("s8", 25, rmsk=2), srow("s7", 25, chrom="chrY", umap=""),
                srow("s2", 24), srow("s3", 26), srow("s6", 40),
                srow("s4", 25, view="A_only"), srow("s5", 25, view="B_only")]
        strong = tmp / "strong.tsv"
        with open(strong, "w", newline="\n") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(STRONG_HEADER)
            w.writerows(rows)
        trace = tmp / "trace.tsv"
        trace.write_text("marker_id\nx\n")
        comp = tmp / "comp.tsv"
        comp.write_text("marker_id\ny\n")
        annotated = tmp / "annotated.tsv"
        with open(annotated, "w", newline="\n") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["marker_id", "umap_singleread_window_coverage", "umap_multiread_min", "umap_multiread_max"])
            for r in rows:
                na = r[10] == ""
                w.writerow([r[0], 1.0, "" if na else 1.0, "" if na else 1.0])
        prior = tmp / "prior.json"
        prior.write_text(json.dumps({"inputs": {"annotated": {"path": str(annotated), "sha256": sha(annotated)}}}))
        stage2 = tmp / "stage2.json"
        stage2.write_text(json.dumps({
            "outputs": {"traceability_all_representatives": {"path": str(trace), "sha256": sha(trace)},
                        "strong_mono_primary": {"path": str(strong), "sha256": sha(strong)},
                        "strong_str_comparison": {"path": str(comp), "sha256": sha(comp)}},
            "hcc1395_results_used": False, "classification_performed": False, "production_changes": False,
            "final_panel_size_fixed": False, "previous_250bp_outputs_modified": False,
            "prior_manifest": {"path": str(prior)}}))
        tier1 = tmp / "tier1.tsv"
        tier1.write_text("marker_id\tstart\tend\n" + "\n".join(
            f"p{i}\t{100 * i}\t{100 * i + s}" for i, s in enumerate([25, 27, 21, 23, 26])) + "\n")
        pman = tmp / "pman.json"
        pman.write_text(json.dumps({"marker_file_sha256": sha(tier1)}))
        return stage2, annotated, tier1, pman, strong

    def _run(self, tmp, extra=()):
        stage2, annotated, tier1, pman, strong = self._build(tmp)
        out = tmp / "out"
        rc = r3.main(["--stage2-manifest", str(stage2), "--annotated", str(annotated), "--promega-tier1", str(tier1),
                      "--promega-manifest", str(pman), "--out-dir", str(out), "--skip-expected-count-check", *extra])
        self.assertEqual(rc, 0)
        return out, (stage2, annotated, tier1, pman, strong)

    def test_ranking_order_hierarchy_and_sensitivity(self):
        with tempfile.TemporaryDirectory() as t:
            out, _ = self._run(Path(t))
            b = list(csv.DictReader(open(out / "candidate_ranking_stage3.v1.ranked_viewB.tsv"), delimiter="\t"))
            self.assertEqual([r["marker_id"] for r in b], ["s1", "s8", "s7", "s2", "s3", "s6", "s5"])
            a = list(csv.DictReader(open(out / "candidate_ranking_stage3.v1.ranked_viewA.tsv"), delimiter="\t"))
            self.assertEqual([r["marker_id"] for r in a], ["s1", "s8", "s7", "s2", "s3", "s6", "s4"])
            union = {r["marker_id"]: r for r in csv.DictReader(
                open(out / "candidate_ranking_stage3.v1.ranked_union.tsv"), delimiter="\t")}
            self.assertEqual(len(union), 8)
            # no-shared-key sensitivity ranking lets the view-only locus s5 rise above distance-1 shared loci
            self.assertLess(int(union["s5"]["rank_viewB_no_shared_key"]), int(union["s2"]["rank_viewB_no_shared_key"]))
            self.assertGreater(int(union["s5"]["rank_viewB"]), int(union["s2"]["rank_viewB"]))
            self.assertEqual(union["s4"]["rank_viewB"], "")                 # A-only is not in View B
            self.assertEqual(union["s7"]["umap_status"], "NOT_ASSESSED")    # never dropped, never imputed
            self.assertEqual(union["s7"]["umap_multiread_min"], "")
            self.assertEqual(union["s7"]["is_chrY"], "True")                # chrY kept

    def test_manifest_contents_and_no_size_selected(self):
        with tempfile.TemporaryDirectory() as t:
            out, _ = self._run(Path(t))
            m = json.loads((out / "candidate_ranking_stage3.v1.manifest.json").read_text())
            self.assertFalse(m["panel_size_selected"])
            self.assertFalse(m["parameters"]["chromosome_quotas"])
            self.assertEqual(m["parameters"]["promega_span_anchor_bp"], 25)
            self.assertEqual(m["counts"], {"viewA": 7, "viewB": 7, "shared": 6, "union": 8})
            self.assertIn("top_100pct", m["report"]["viewB"]["primary_ranking_prefixes"])
            self.assertEqual(m["report"]["viewB"]["primary_ranking_prefixes"]["top_100pct"]["chrY"], 1)
            self.assertIn("UNAVAILABLE", m["unavailable_information"]["flank_kmer_multiplicity_or_other_uniqueness"])

    def test_reproducible_identical_outputs_on_rerun(self):
        with tempfile.TemporaryDirectory() as t:
            tmp = Path(t)
            out1, inputs = self._run(tmp)
            first = {p.name: sha(p) for p in out1.glob("*.tsv")}
            stage2, annotated, tier1, pman, _ = inputs
            out2 = tmp / "out2"
            r3.main(["--stage2-manifest", str(stage2), "--annotated", str(annotated), "--promega-tier1", str(tier1),
                     "--promega-manifest", str(pman), "--out-dir", str(out2), "--skip-expected-count-check"])
            self.assertEqual(first, {p.name: sha(p) for p in out2.glob("*.tsv")})

    def test_independent_verifier_passes_and_detects_tampering(self):
        with tempfile.TemporaryDirectory() as t:
            out, _ = self._run(Path(t))
            manifest = str(out / "candidate_ranking_stage3.v1.manifest.json")
            self.assertEqual(ver3.verify(manifest), [])
            path = out / "candidate_ranking_stage3.v1.ranked_viewB.tsv"
            path.write_text(path.read_text().replace("s1\tA_and_B", "sX\tA_and_B", 1))
            self.assertTrue(ver3.verify(manifest))

    def test_refuses_when_stage2_output_hash_mismatches(self):
        with tempfile.TemporaryDirectory() as t:
            tmp = Path(t)
            stage2, annotated, tier1, pman, strong = self._build(tmp)
            strong.write_text(strong.read_text() + "\n")
            with self.assertRaises(SystemExit):
                r3.main(["--stage2-manifest", str(stage2), "--annotated", str(annotated), "--promega-tier1", str(tier1),
                         "--promega-manifest", str(pman), "--out-dir", str(tmp / "o"), "--skip-expected-count-check"])

    def test_refuses_on_wrong_expected_manifest_hash_and_wrong_counts(self):
        with tempfile.TemporaryDirectory() as t:
            tmp = Path(t)
            stage2, annotated, tier1, pman, _ = self._build(tmp)
            base = ["--stage2-manifest", str(stage2), "--annotated", str(annotated), "--promega-tier1", str(tier1),
                    "--promega-manifest", str(pman), "--out-dir", str(tmp / "o")]
            with self.assertRaises(SystemExit):
                r3.main(base + ["--skip-expected-count-check", "--expect-stage2-manifest-sha256", "0" * 64])
            with self.assertRaises(SystemExit):
                r3.main(base)             # synthetic counts differ from the validated 6,074 / 60,756 / 6,038

    def test_refuses_when_annotation_or_promega_input_changed(self):
        with tempfile.TemporaryDirectory() as t:
            tmp = Path(t)
            stage2, annotated, tier1, pman, _ = self._build(tmp)
            tier1.write_text(tier1.read_text() + "\n")
            with self.assertRaises(SystemExit):
                r3.main(["--stage2-manifest", str(stage2), "--annotated", str(annotated), "--promega-tier1", str(tier1),
                         "--promega-manifest", str(pman), "--out-dir", str(tmp / "o"), "--skip-expected-count-check"])


if __name__ == "__main__":
    unittest.main()
