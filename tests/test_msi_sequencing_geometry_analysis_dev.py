"""
SYNTHETIC, DEVELOPMENT-ONLY tests for
dev/msi_redesign/sequencing_geometry_analysis.py.

*** Hand-computed fragment geometries and fabricated SAM lines only. The ***
*** real geometry numbers come from a separate run against the real ***
*** HCC1395 BAMs, used as sequencing-geometry references only.

Run with:
    python -m unittest tests.test_msi_sequencing_geometry_analysis_dev -v
"""

import sys
import unittest
from collections import Counter
from pathlib import Path

DEV_DIR = Path(__file__).resolve().parents[1] / "dev" / "msi_redesign"
sys.path.insert(0, str(DEV_DIR))

import sequencing_geometry_analysis as sga  # noqa: E402


class TestBasics(unittest.TestCase):
    def test_percentile(self):
        vals = [1, 2, 3, 4, 5]
        self.assertEqual(sga.percentile(vals, 50), 3)
        self.assertEqual(sga.percentile(vals, 25), 2)
        self.assertEqual(sga.percentile(vals, 0), 1)
        self.assertEqual(sga.percentile(vals, 100), 5)
        self.assertIsNone(sga.percentile([], 50))

    def test_ref_span_from_cigar(self):
        self.assertEqual(sga.ref_span_from_cigar("100M"), 100)
        self.assertEqual(sga.ref_span_from_cigar("76M2I10M"), 86)      # I does not consume ref
        self.assertEqual(sga.ref_span_from_cigar("10S60M2D30M"), 92)   # S excluded, D included

    def test_merge_and_intersect(self):
        self.assertEqual(sga.merge_intervals([(5, 10), (11, 15), (20, 25)]), [(5, 15), (20, 25)])
        self.assertEqual(sga.intersect_intervals([(0, 10), (20, 30)], [(5, 25)]), [(5, 10), (20, 25)])
        self.assertEqual(sga.total_length([(5, 15), (20, 25)]), 17)

    def test_summarize_values(self):
        s = sga.summarize_values([100, 200, 300])
        self.assertEqual(s["n"], 3)
        self.assertEqual(s["percentiles"]["P50"], 200)


class TestFragmentGeometry(unittest.TestCase):
    def test_adjacent_loci_two_non_overlapping_mates(self):
        # T=200, R=100, L=10, flank=5, gap=0 (hand-computed: 162 / 142 / 142)
        self.assertEqual(sga.fragment_joint_counts(200, 100, 10, 10, 0), (162, 142, 142))

    def test_cross_mate_joint_is_not_same_read(self):
        # gap=100: locus 1 in mate 1 and locus 2 in mate 2 -> fragment-joint only
        self.assertEqual(sga.fragment_joint_counts(200, 100, 10, 10, 100), (162, 0, 71))

    def test_far_apart_loci_are_independent(self):
        self.assertEqual(sga.fragment_joint_counts(200, 100, 10, 10, 200), (162, 0, 0))

    def test_gap_20(self):
        self.assertEqual(sga.fragment_joint_counts(200, 100, 10, 10, 20), (162, 102, 113))

    def test_fragment_too_short_to_measure_anything(self):
        self.assertEqual(sga.fragment_joint_counts(15, 100, 10, 10, 0), (0, 0, 0))

    def test_same_read_never_exceeds_fragment_joint(self):
        for t in (120, 200, 350):
            for gap in (0, 20, 50, 100, 150):
                n1, ns, nf = sga.fragment_joint_counts(t, 100, 10, 10, gap)
                self.assertLessEqual(ns, nf)
                self.assertLessEqual(nf, n1)


class TestPooledRedundancy(unittest.TestCase):
    def test_pooling_weights_by_locus1_offsets(self):
        counter = Counter({(200, 100): 3})
        r = sga.pooled_redundancy(counter, 10, 10, 0)
        self.assertAlmostEqual(r["fragment_joint"], 142 / 162)
        self.assertAlmostEqual(r["same_read_joint"], 142 / 162)

    def test_empty_counter(self):
        self.assertEqual(sga.pooled_redundancy(Counter(), 10, 10, 0),
                          {"same_read_joint": None, "fragment_joint": None})

    def test_redundancy_decreases_with_spacing_for_fixed_geometry(self):
        counter = Counter({(200, 100): 5, (300, 100): 5})
        vals = [sga.pooled_redundancy(counter, 10, 10, g)["fragment_joint"]
                for g in (0, 50, 100, 200, 400)]
        self.assertEqual(vals, sorted(vals, reverse=True))
        self.assertEqual(vals[-1], 0.0)

    def test_spacing_crossing(self):
        curve = {0: {"k": 0.9}, 10: {"k": 0.4}, 20: {"k": 0.04}, 30: {"k": 0.0}}
        self.assertEqual(sga.spacing_crossing(curve, "k", 0.5), 10)
        self.assertEqual(sga.spacing_crossing(curve, "k", 0.05), 20)
        self.assertIsNone(sga.spacing_crossing({0: {"k": 0.9}}, "k", 0.1))


class TestParsing(unittest.TestCase):
    def _line(self, **kw):
        f = ["read1", "99", "chr17", "1000", "60", "100M", "=", "1150", "-250", "A" * 100]
        return "\t".join(f)

    def test_records_parsed_and_tlen_made_absolute(self):
        out = list(sga.read_geometry_records([self._line()]))
        self.assertEqual(out, [(250, 100, 100)])

    def test_skips_other_reference_and_zero_tlen_and_short_lines(self):
        bad_ref = self._line().replace("\t=\t", "\tchr2\t")
        zero_tlen = self._line().replace("-250", "0")
        short = "a\tb\tc"
        self.assertEqual(list(sga.read_geometry_records([bad_ref, zero_tlen, short])), [])


class TestLociAffected(unittest.TestCase):
    def test_counts_and_greedy_scenario(self):
        by_chrom = {"chr1": [(0, 10), (30, 40), (500, 510)]}
        r = sga.loci_affected_by_spacing(by_chrom, thresholds=[10, 20, 500])
        s = r["by_spacing_bp"]
        self.assertEqual(r["n_total"], 3)
        self.assertEqual(s["10"]["n_loci_with_neighbor_within"], 0)
        self.assertEqual(s["20"]["n_loci_with_neighbor_within"], 2)
        self.assertEqual(s["500"]["n_loci_with_neighbor_within"], 3)
        self.assertEqual(s["10"]["n_loci_remaining_if_one_per_cluster_scenario_only"], 3)
        self.assertEqual(s["20"]["n_loci_remaining_if_one_per_cluster_scenario_only"], 2)
        self.assertEqual(s["500"]["n_loci_remaining_if_one_per_cluster_scenario_only"], 1)

    def test_no_locus_is_ever_removed_from_input(self):
        by_chrom = {"chr1": [(0, 10), (30, 40)]}
        before = [list(v) for v in by_chrom.values()]
        sga.loci_affected_by_spacing(by_chrom, thresholds=[50])
        self.assertEqual([list(v) for v in by_chrom.values()], before)


class TestDescribeSample(unittest.TestCase):
    def test_describe(self):
        tlens = [200, 250, 300]
        halves = [Counter({(200, 100): 1, (300, 100): 1}), Counter({(250, 100): 1})]
        desc, full = sga.describe_sample(tlens, [100, 100, 100], halves)
        self.assertEqual(desc["n_fragments_sampled"], 3)
        self.assertEqual(desc["read_length_distribution"], {100: 3})
        self.assertEqual(desc["fraction_mates_overlapping_(T<2R)"], 0.0)
        self.assertEqual(sum(full.values()), 3)


if __name__ == "__main__":
    unittest.main()
