"""
SYNTHETIC, DEVELOPMENT-ONLY tests for
dev/msi_redesign/build_candidate_panel_spacing.py.

*** Fabricated loci only. The real Candidate MSI Panel v1 (development ***
*** candidate, not clinically validated) comes from a separate real run ***
*** against View A / View B, reported with its own manifest and SHA-256. ***

Run with:
    python -m unittest tests.test_msi_build_candidate_panel_spacing_dev -v
"""

import csv
import json
import random
import sys
import tempfile
import unittest
from pathlib import Path

DEV_DIR = Path(__file__).resolve().parents[1] / "dev" / "msi_redesign"
sys.path.insert(0, str(DEV_DIR))

import build_candidate_panel_spacing as bcp  # noqa: E402


def L(mid, start, end, kind="mononucleotide", motif="A", rep_len=None, umap=1.0,
      rmsk="no_rmsk_hit", flank=False, chrom="chr1"):
    return bcp.Locus(mid, chrom, start, end, kind, motif, rep_len or (end - start), umap, rmsk, flank)


class TestTiers(unittest.TestCase):
    def test_umap_tier(self):
        self.assertEqual([bcp.umap_tier(x) for x in (1.0, 0.95, 0.8, 0.5, None)], [0, 1, 2, 3, 4])

    def test_rmsk_tier(self):
        self.assertEqual(bcp.rmsk_tier("microsatellite_class_only"), 0)
        self.assertEqual(bcp.rmsk_tier("no_rmsk_hit"), 1)
        self.assertEqual(bcp.rmsk_tier("mixed_microsatellite_and_other"), 2)

    def test_kind_and_length_tiers(self):
        self.assertEqual(bcp.kind_tier("mononucleotide"), 0)
        self.assertEqual(bcp.kind_tier("str"), 1)
        self.assertEqual([bcp.length_tier(s) for s in (8, 12, 20, 45, 70)], [2, 1, 0, 1, 2])
        self.assertEqual([bcp.extraction_span_tier(s) for s in (30, 60, 61, 90)], [0, 0, 1, 1])

    def test_eligibility_follows_read_length_and_flank(self):
        self.assertEqual(bcp.MAX_MEASURABLE_SPAN_BP, 90)
        self.assertTrue(bcp.is_eligible(L("a", 0, 90)))
        self.assertFalse(bcp.is_eligible(L("a", 0, 91)))


class TestPriority(unittest.TestCase):
    def test_umap_dominates(self):
        good = L("g", 0, 20, umap=1.0, kind="str")
        bad = L("b", 0, 20, umap=0.8)
        self.assertLess(bcp.priority_key(good), bcp.priority_key(bad))

    def test_mononucleotide_beats_str_when_umap_and_rmsk_equal(self):
        self.assertLess(bcp.priority_key(L("m", 0, 20)), bcp.priority_key(L("s", 0, 20, kind="str")))

    def test_length_tier_prefers_informative_span(self):
        self.assertLess(bcp.priority_key(L("a", 0, 20)), bcp.priority_key(L("b", 0, 9)))

    def test_tie_break_is_deterministic(self):
        a, b = L("id_a", 0, 20), L("id_b", 0, 20)
        self.assertNotEqual(bcp.priority_key(a), bcp.priority_key(b))     # hash separates otherwise-equal loci
        self.assertEqual(bcp.priority_key(a), bcp.priority_key(L("id_a", 0, 20)))


class TestClusters(unittest.TestCase):
    def test_gap_boundary_and_chaining(self):
        loci = [L("a", 0, 10), L("b", 100, 110), L("c", 400, 410)]
        self.assertEqual(bcp.assign_clusters(loci, 250), [0, 0, 1])
        chain = [L("a", 0, 10), L("b", 200, 210), L("c", 400, 410)]
        self.assertEqual(bcp.assign_clusters(chain, 250), [0, 0, 0])   # single linkage chains

    def test_exact_spacing_is_same_cluster(self):
        self.assertEqual(bcp.assign_clusters([L("a", 0, 10), L("b", 260, 270)], 250), [0, 0])
        self.assertEqual(bcp.assign_clusters([L("a", 0, 10), L("b", 261, 271)], 250), [0, 1])


class TestSelection(unittest.TestCase):
    def test_best_locus_in_cluster_wins_and_excludes_neighbours(self):
        loci = [L("a", 0, 10, umap=0.8), L("b", 200, 210), L("c", 400, 410, umap=0.8)]
        status, rep = bcp.select_representatives(loci, 250)
        self.assertEqual(status, ["EXCLUDED_BY_SPACING", "SELECTED", "EXCLUDED_BY_SPACING"])
        self.assertEqual(rep, [1, None, 1])

    def test_chain_can_yield_multiple_representatives(self):
        loci = [L("a", 0, 10), L("b", 200, 210, umap=0.8), L("c", 400, 410, umap=0.9)]
        status, rep = bcp.select_representatives(loci, 250)
        self.assertEqual(status, ["SELECTED", "EXCLUDED_BY_SPACING", "SELECTED"])
        self.assertEqual(rep[1], 0)

    def test_gap_exactly_at_spacing_is_excluded_one_more_is_kept(self):
        s1, _ = bcp.select_representatives([L("a", 0, 10), L("b", 260, 270, umap=0.8)], 250)
        self.assertEqual(s1, ["SELECTED", "EXCLUDED_BY_SPACING"])
        s2, _ = bcp.select_representatives([L("a", 0, 10), L("b", 261, 271, umap=0.8)], 250)
        self.assertEqual(s2, ["SELECTED", "SELECTED"])

    def test_ineligible_locus_is_never_selected_and_does_not_exclude(self):
        loci = [L("long", 0, 95, kind="str", motif="AC"), L("near", 140, 160, umap=0.8)]
        status, rep = bcp.select_representatives(loci, 250)
        self.assertEqual(status, ["INELIGIBLE_EXTRACTION_SPAN", "SELECTED"])

    def test_overlapping_loci_conflict(self):
        status, _ = bcp.select_representatives([L("a", 0, 20), L("b", 10, 30, umap=0.8)], 250)
        self.assertEqual(status, ["SELECTED", "EXCLUDED_BY_SPACING"])

    def test_left_scan_finds_long_earlier_locus(self):
        # long earlier locus (span 80) must be found when a later, better locus is selected first
        loci = [L("long", 0, 80, umap=0.8), L("mid", 250, 255, umap=0.8), L("best", 300, 310)]
        status, rep = bcp.select_representatives(loci, 250)
        self.assertEqual(status[2], "SELECTED")
        self.assertEqual(status[0], "EXCLUDED_BY_SPACING")   # gap 300-80 = 220 <= 250
        self.assertEqual(rep[0], 2)

    def test_random_property_min_gap_and_traceability(self):
        rng = random.Random(7)
        for spacing in (75, 250):
            pos = sorted(rng.sample(range(0, 20000, 7), 400))
            loci = [L(f"m{i}", p, p + rng.choice([8, 10, 15, 30, 95]),
                      umap=rng.choice([1.0, 0.95, 0.8]),
                      kind=rng.choice(["mononucleotide", "str"])) for i, p in enumerate(pos)]
            status, rep = bcp.select_representatives(loci, spacing)
            sel = [i for i, s in enumerate(status) if s == "SELECTED"]
            for x, y in zip(sel, sel[1:]):
                self.assertGreater(loci[y].start - loci[x].end, spacing)
            for i, s in enumerate(status):
                if s == "EXCLUDED_BY_SPACING":
                    r = rep[i]
                    self.assertEqual(status[r], "SELECTED")
                    gap = max(loci[i].start - loci[r].end, loci[r].start - loci[i].end)
                    self.assertLessEqual(gap, spacing)
                if s == "INELIGIBLE_EXTRACTION_SPAN":
                    self.assertGreater(loci[i].end - loci[i].start, 90)
                self.assertIsNotNone(s)


class TestProcessView(unittest.TestCase):
    def test_no_locus_lost_and_counts_add_up(self):
        loci = [L("a", 0, 10, chrom="chr1"), L("b", 100, 110, umap=0.8, chrom="chr1"),
                L("c", 5000, 5010, chrom="chr1"), L("d", 0, 10, chrom="chr2"),
                L("e", 0, 95, kind="str", motif="AC", chrom="chr2")]
        records, s = bcp.process_view(loci, 250)
        self.assertEqual(len(records), 5)
        self.assertEqual(s["n_selected_representatives"] + s["n_excluded_by_spacing"]
                         + s["n_ineligible_extraction_span"], 5)
        self.assertEqual(s["n_selected_representatives"], 3)
        self.assertEqual(s["n_excluded_by_spacing"], 1)
        self.assertEqual(s["n_clusters_total_including_singletons"], 3)
        self.assertEqual(s["n_clusters_with_2plus_loci"], 2)
        ids = {r["locus"].marker_id: r for r in records}
        self.assertEqual(ids["b"]["representative_id"], "a")
        self.assertEqual(ids["a"]["cluster_id"], ids["b"]["cluster_id"])
        self.assertTrue(ids["a"]["cluster_id"].startswith("chr1:C"))


class TestComposition(unittest.TestCase):
    def test_counts(self):
        c = bcp.composition([L("a", 0, 20), L("b", 0, 12, kind="str", motif="AC", chrom="chr2", umap=None)])
        self.assertEqual(c["n"], 2)
        self.assertEqual(c["kind"], {"mononucleotide": 1, "str": 1})
        self.assertEqual(c["motif_length"], {1: 1, 2: 1})
        self.assertEqual(c["chromosome"], {"chr1": 1, "chr2": 1})
        self.assertEqual(c["umap_tier"], {0: 1, 4: 1})


ANNOT_HEADER = ["marker_id", "chrom", "start", "end", "motif", "ref_repeat_length", "kind",
                "umap_multiread_mean", "repeatmasker_self_category",
                "repeatmasker_flank_other_repeat_flag", "kept"]


def annot_row(mid, chrom, start, end, kept="True", mean="1.0", kind="mononucleotide", motif="A", cat="no_rmsk_hit"):
    return [mid, chrom, start, end, motif, end - start, kind, mean, cat, "False", kept]


class TestEndToEnd(unittest.TestCase):
    def test_main_writes_outputs_and_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            annotated = tmp / "annot.tsv"
            with open(annotated, "w", newline="\n") as fh:
                w = csv.writer(fh, delimiter="\t", lineterminator="\n")
                w.writerow(ANNOT_HEADER)
                w.writerows([
                    annot_row("m1", "chr1", 100, 120),
                    annot_row("m2", "chr1", 200, 215, mean="0.8"),
                    annot_row("m3", "chr1", 9000, 9020),
                    annot_row("m4", "chr1", 9100, 9120, kept="False"),   # excluded upstream: must not appear
                ])
            view_a = tmp / "view_a.tsv"
            view_a.write_text("marker_id\tchrom\nm1\tchr1\nm3\tchr1\n")
            out = tmp / "out"
            self.assertEqual(bcp.main(["--annotated", str(annotated), "--view-a", str(view_a),
                                        "--out-dir", str(out)]), 0)
            manifest = json.loads((out / "candidate_msi_panel.v1.spacing250.manifest.json").read_text())
            self.assertEqual(manifest["status"], "DEVELOPMENT_CANDIDATE_NOT_CLINICALLY_VALIDATED")
            self.assertFalse(manifest["final_panel_size_fixed"])
            self.assertFalse(manifest["hcc1395_variant_or_msi_results_used"])
            b = manifest["views"]["viewB_wgs"]["spacing_summary"]
            a = manifest["views"]["viewA_wes_compatible"]["spacing_summary"]
            self.assertEqual(b["n_loci_in_view"], 3)          # m4 (kept=False) never enters
            self.assertEqual(b["n_selected_representatives"], 2)
            self.assertEqual(a["n_loci_in_view"], 2)
            panel = out / "candidate_msi_panel.v1.viewB_wgs.spacing250.tsv"
            with open(panel, encoding="utf-8") as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
            self.assertEqual(sorted(r["marker_id"] for r in rows), ["m1", "m3"])
            trace = out / "candidate_msi_panel.v1.viewB_wgs.spacing250.traceability.tsv"
            with open(trace, encoding="utf-8") as fh:
                trows = {r["marker_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
            self.assertEqual(trows["m2"]["status"], "EXCLUDED_BY_SPACING")
            self.assertEqual(trows["m2"]["representative_id"], "m1")
            self.assertEqual(len(manifest["outputs"]["viewB_wgs"]["selected_panel"]["sha256"]), 64)


if __name__ == "__main__":
    unittest.main()
