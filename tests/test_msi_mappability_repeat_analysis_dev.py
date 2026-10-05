"""
SYNTHETIC, DEVELOPMENT-ONLY tests for
dev/msi_redesign/mappability_repeat_analysis.py's pure aggregation logic.

*** These use fabricated (id, class/value) pairs, not real RepeatMasker/ ***
*** Umap/segdup data. The real feasibility numbers reported to the user ***
*** come from a separate real run against the real downloaded tracks, ***
*** with its own provenance (source URLs, release dates, SHA-256).

Run with:
    python -m unittest tests.test_msi_mappability_repeat_analysis_dev -v
"""

import sys
import tempfile
import unittest
from pathlib import Path

DEV_DIR = Path(__file__).resolve().parents[1] / "dev" / "msi_redesign"
sys.path.insert(0, str(DEV_DIR))

import mappability_repeat_analysis as mra  # noqa: E402


class TestBinScore(unittest.TestCase):
    def test_exact_zero(self):
        self.assertEqual(mra.bin_score(0.0), "0.0_exactly")

    def test_exact_one(self):
        self.assertEqual(mra.bin_score(1.0), "1.0_exactly")

    def test_mid_range_bins(self):
        self.assertEqual(mra.bin_score(0.05), "0.0-0.2")
        self.assertEqual(mra.bin_score(0.35), "0.2-0.4")
        self.assertEqual(mra.bin_score(0.55), "0.4-0.6")
        self.assertEqual(mra.bin_score(0.75), "0.6-0.8")
        self.assertEqual(mra.bin_score(0.95), "0.8-1.0_exclusive")

    def test_bin_boundaries_are_half_open_low_inclusive(self):
        self.assertEqual(mra.bin_score(0.2), "0.2-0.4")
        self.assertEqual(mra.bin_score(0.4), "0.4-0.6")


class TestClassifyRmskClasses(unittest.TestCase):
    def test_no_hit(self):
        self.assertEqual(mra.classify_rmsk_classes(set()), "no_rmsk_hit")

    def test_only_microsatellite_classes(self):
        self.assertEqual(mra.classify_rmsk_classes({"Simple_repeat"}), "microsatellite_class_only")
        self.assertEqual(
            mra.classify_rmsk_classes({"Simple_repeat", "Low_complexity"}),
            "microsatellite_class_only",
        )

    def test_mixed_classes(self):
        self.assertEqual(
            mra.classify_rmsk_classes({"Simple_repeat", "SINE"}),
            "mixed_microsatellite_and_other",
        )

    def test_unexpected_non_microsatellite_class(self):
        self.assertEqual(mra.classify_rmsk_classes({"SINE"}), "unexpected_non_microsatellite_class")
        self.assertEqual(
            mra.classify_rmsk_classes({"LINE", "LTR"}),
            "unexpected_non_microsatellite_class",
        )


class TestLoaders(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def _write(self, name, lines):
        path = Path(self.tmp.name) / name
        path.write_text("\n".join(lines) + "\n")
        return str(path)

    def test_load_pair_tsv_groups_by_id(self):
        path = self._write("pairs.tsv", ["m1\tSimple_repeat", "m1\tLow_complexity", "m2\tSINE"])
        result = mra.load_pair_tsv(path)
        self.assertEqual(result["m1"], {"Simple_repeat", "Low_complexity"})
        self.assertEqual(result["m2"], {"SINE"})

    def test_load_segdup_hits_parses_floats(self):
        path = self._write("segdup.tsv", ["m1\t0.98", "m1\t0.95", "m2\t0.99"])
        result = mra.load_segdup_hits(path)
        self.assertEqual(result["m1"], [0.98, 0.95])
        self.assertEqual(result["m2"], [0.99])

    def test_load_umap_singleread_coverage_takes_fraction_column(self):
        # bedtools coverage: chrom start end marker_id ... fraction(last col)
        path = self._write("cov.tsv", ["chr1\t0\t100\tm1\t1\t100\t100\t1.0000000",
                                        "chr1\t200\t300\tm2\t0\t0\t100\t0.0000000"])
        result = mra.load_umap_singleread_coverage(path)
        self.assertEqual(result["m1"], 1.0)
        self.assertEqual(result["m2"], 0.0)

    def test_load_umap_multiread_scores_skips_dot_rows(self):
        path = self._write("scores.tsv", ["chr1\t0\t100\tm1\t0.5\t0.1\t0.9",
                                           "chr1\t200\t300\tm2\t.\t.\t."])
        result = mra.load_umap_multiread_scores(path)
        self.assertEqual(result["m1"], (0.5, 0.1, 0.9))
        self.assertNotIn("m2", result)


class TestAnalyzeEndToEnd(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def _write(self, name, lines):
        path = Path(self.tmp.name) / name
        path.write_text("\n".join(lines) + "\n")
        return str(path)

    def test_full_pipeline_on_three_synthetic_loci(self):
        markers_path = self._write("markers.tsv", [
            "marker_id\tchrom\tstart\tend\tmotif\tref_repeat_length\tkind\ttier\ttier_label",
            "m1\tchr1\t100\t110\tA\t10\tmononucleotide\t2\tgenome_scan_quantitative",
            "m2\tchr1\t500\t510\tA\t10\tmononucleotide\t2\tgenome_scan_quantitative",
            "m3\tchr1\t900\t910\tA\t10\tmononucleotide\t2\tgenome_scan_quantitative",
        ])
        # m1: clean microsatellite-only self hit, no flank issue, full mappability, no segdup
        # m2: flank overlaps a SINE (problematic), partial mappability
        # m3: window overlaps a segmental duplication
        locus_self = self._write("self.tsv", ["m1\tSimple_repeat", "m2\tSimple_repeat", "m3\tSimple_repeat"])
        flank = self._write("flank.tsv", ["m2\tSINE"])
        segdup = self._write("segdup.tsv", ["m3\t0.97"])
        umap_single = self._write("umap_single.tsv", [
            "chr1\t0\t210\tm1\t1\t210\t210\t1.0000000",
            "chr1\t400\t610\tm2\t1\t100\t210\t0.4761905",
            "chr1\t800\t1010\tm3\t1\t210\t210\t1.0000000",
        ])
        umap_multi = self._write("umap_multi.tsv", [
            "chr1\t0\t210\tm1\t0.95\t0.90\t1.00",
            "chr1\t400\t610\tm2\t0.30\t0.00\t0.60",
            "chr1\t800\t1010\tm3\t0.99\t0.98\t1.00",
        ])

        result = mra.analyze(markers_path, locus_self, flank, segdup, umap_single, umap_multi)

        self.assertEqual(result["n_total_tier2_candidates"], 3)
        self.assertEqual(result["repeatmasker_locus_self_classification"],
                          {"microsatellite_class_only": 3})
        self.assertEqual(result["n_loci_with_flank_other_repeat_class"], 1)
        self.assertEqual(result["segmental_duplication"]["n_loci_window_overlapping_segdup"], 1)
        self.assertEqual(result["umap_single_read_binary_window_coverage"]["n_fully_covered_window"], 2)
        self.assertEqual(result["filter_scenarios_not_a_final_panel"]["raw_candidate_universe"], 3)
        self.assertEqual(result["filter_scenarios_not_a_final_panel"]["exclude_segdup_overlap"], 2)
        self.assertEqual(
            result["filter_scenarios_not_a_final_panel"]["exclude_flank_other_repeat_class"], 2
        )


if __name__ == "__main__":
    unittest.main()
