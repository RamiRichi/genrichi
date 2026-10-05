"""
SYNTHETIC, DEVELOPMENT-ONLY tests for
dev/msi_redesign/candidate_downselection_analysis.py.

*** Fabricated marker/annotation data only. The real down-selection ***
*** report comes from a separate real run against View A/View B and the ***
*** real reference-derived annotations, with its own provenance.

Run with:
    python -m unittest tests.test_msi_candidate_downselection_analysis_dev -v
"""

import sys
import unittest
from pathlib import Path

DEV_DIR = Path(__file__).resolve().parents[1] / "dev" / "msi_redesign"
sys.path.insert(0, str(DEV_DIR))

import candidate_downselection_analysis as cda  # noqa: E402


class TestRepeatLengthBand(unittest.TestCase):
    def test_bands(self):
        self.assertEqual(cda.repeat_length_band(5), "5-9bp")
        self.assertEqual(cda.repeat_length_band(9), "5-9bp")
        self.assertEqual(cda.repeat_length_band(10), "10-14bp")
        self.assertEqual(cda.repeat_length_band(20), "15-20bp")
        self.assertEqual(cda.repeat_length_band(21), "21-30bp")
        self.assertEqual(cda.repeat_length_band(31), "31bp+")
        self.assertEqual(cda.repeat_length_band(300), "31bp+")


class TestUmapMeanBand(unittest.TestCase):
    def test_no_data(self):
        self.assertEqual(cda.umap_mean_band(""), "no_data")
        self.assertEqual(cda.umap_mean_band(None), "no_data")

    def test_exact_one(self):
        self.assertEqual(cda.umap_mean_band(1.0), "1.0_exactly")
        self.assertEqual(cda.umap_mean_band("1.0"), "1.0_exactly")

    def test_bands(self):
        self.assertEqual(cda.umap_mean_band(0.95), "0.9-1.0_exclusive")
        self.assertEqual(cda.umap_mean_band(0.8), "0.7-0.9")
        self.assertEqual(cda.umap_mean_band(0.5), "below_0.7")


class TestNearestNeighborMinGaps(unittest.TestCase):
    def test_single_locus_has_no_neighbor(self):
        self.assertEqual(cda.nearest_neighbor_min_gaps([(100, 110)]), [None])

    def test_two_loci_gap_computed_both_directions(self):
        positions = [(100, 110), (130, 140)]
        gaps = cda.nearest_neighbor_min_gaps(positions)
        self.assertEqual(gaps, [20, 20])

    def test_three_loci_picks_minimum_of_both_sides(self):
        positions = [(0, 10), (30, 40), (200, 210)]
        gaps = cda.nearest_neighbor_min_gaps(positions)
        # locus 2 (30,40): left gap=20, right gap=160 -> min=20
        self.assertEqual(gaps, [20, 20, 160])

    def test_overlapping_intervals_give_negative_gap(self):
        positions = [(0, 20), (10, 30)]
        gaps = cda.nearest_neighbor_min_gaps(positions)
        self.assertEqual(gaps, [-10, -10])


class TestSpacingScenarioCounts(unittest.TestCase):
    def test_counts_within_thresholds(self):
        min_gaps_by_chrom = {"chr1": [10, 60, None, 5]}
        counts = cda.spacing_scenario_counts(min_gaps_by_chrom, thresholds=[20, 100])
        self.assertEqual(counts["n_loci_with_neighbor_within_20bp"], 2)  # 10, 5
        self.assertEqual(counts["n_loci_with_neighbor_within_100bp"], 3)  # 10, 60, 5
        self.assertEqual(counts["n_loci_with_no_neighbor_on_chromosome"], 1)
        self.assertEqual(counts["n_total"], 4)


class TestGapProximityScenarioCounts(unittest.TestCase):
    def test_counts_within_thresholds(self):
        distances = [0, 5, 50, 200, None]
        counts = cda.gap_proximity_scenario_counts(distances, thresholds=[0, 10, 100])
        self.assertEqual(counts["n_loci_within_0bp_of_assembly_gap"], 1)
        self.assertEqual(counts["n_loci_within_10bp_of_assembly_gap"], 2)
        self.assertEqual(counts["n_loci_within_100bp_of_assembly_gap"], 3)
        self.assertEqual(counts["n_total_with_distance_data"], 4)


class TestAnalyzeView(unittest.TestCase):
    def test_basic_aggregation(self):
        marker_ids = {"m1", "m2", "m3"}
        rows = [
            {"marker_id": "m1", "chrom": "chr1", "start": "100", "end": "110",
             "motif": "A", "ref_repeat_length": "10", "kind": "mononucleotide",
             "umap_multiread_mean": "1.0"},
            {"marker_id": "m2", "chrom": "chr1", "start": "500", "end": "510",
             "motif": "AC", "ref_repeat_length": "5", "kind": "str",
             "umap_multiread_mean": "0.5"},
            {"marker_id": "m3", "chrom": "chr2", "start": "900", "end": "910",
             "motif": "A", "ref_repeat_length": "10", "kind": "mononucleotide",
             "umap_multiread_mean": ""},
        ]
        result = cda.analyze_view(marker_ids, rows)
        self.assertEqual(result["n_total"], 3)
        self.assertEqual(result["kind_counts"], {"mononucleotide": 2, "str": 1})
        self.assertEqual(result["motif_length_counts"], {1: 2, 2: 1})
        self.assertEqual(result["chrom_counts"], {"chr1": 2, "chr2": 1})
        self.assertEqual(result["umap_multiread_mean_band_counts"],
                          {"1.0_exactly": 1, "below_0.7": 1, "no_data": 1})


if __name__ == "__main__":
    unittest.main()
