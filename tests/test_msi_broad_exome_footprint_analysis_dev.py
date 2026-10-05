"""
SYNTHETIC, DEVELOPMENT-ONLY tests for
dev/msi_redesign/broad_exome_footprint_analysis.py.

*** Fabricated coverage/overlap data only. The real footprint-compatibility ***
*** numbers reported to the user come from a separate real run against the ***
*** real filtered panel and the real Broad exome_calling_regions.v1 file, ***
*** with its own provenance reported alongside it.

Run with:
    python -m unittest tests.test_msi_broad_exome_footprint_analysis_dev -v
"""

import sys
import tempfile
import unittest
from pathlib import Path

DEV_DIR = Path(__file__).resolve().parents[1] / "dev" / "msi_redesign"
sys.path.insert(0, str(DEV_DIR))

import broad_exome_footprint_analysis as bfa  # noqa: E402


class TestClassifyCoverage(unittest.TestCase):
    def test_fully_inside(self):
        self.assertEqual(bfa.classify_coverage(1.0), "fully_inside_footprint")

    def test_outside(self):
        self.assertEqual(bfa.classify_coverage(0.0), "outside_footprint")

    def test_partial(self):
        self.assertEqual(bfa.classify_coverage(0.5), "partially_overlapping_footprint")

    def test_near_one_rounds_to_fully_inside(self):
        self.assertEqual(bfa.classify_coverage(0.9999999), "fully_inside_footprint")


class TestAnalyzeEndToEnd(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def _write(self, name, lines):
        path = Path(self.tmp.name) / name
        path.write_text("\n".join(lines) + "\n")
        return str(path)

    def test_full_pipeline_on_three_synthetic_loci(self):
        panel_path = self._write("filtered.tsv", [
            "marker_id\tchrom\tstart\tend\tmotif\tref_repeat_length\tkind",
            "m1\tchr1\t100\t110\tA\t10\tmononucleotide",
            "m2\tchr1\t500\t510\tAC\t5\tstr",
            "m3\tchr1\t900\t910\tA\t10\tmononucleotide",
        ])
        coverage_path = self._write("coverage.tsv", [
            "chr1\t100\t110\tm1\t1\t10\t10\t1.0000000",
            "chr1\t500\t510\tm2\t1\t3\t10\t0.3000000",
            "chr1\t900\t910\tm3\t0\t0\t10\t0.0000000",
        ])
        gnomad_path = self._write("gnomad.tsv", ["chr1\t100\t110\tm1\tchr1\t95\t115\tHTT"])

        result = bfa.analyze(panel_path, coverage_path, gnomad_path)

        self.assertEqual(result["n_filtered_panel_input"], 3)
        self.assertEqual(result["breakdown_fully_inside_footprint"]["n"], 1)
        self.assertEqual(result["breakdown_partially_overlapping_footprint"]["n"], 1)
        self.assertEqual(result["breakdown_outside_footprint"]["n"], 1)
        self.assertEqual(result["n_gnomad_disease_str_catalog_overlap_documentation_only"], 1)
        self.assertFalse(result["gnomad_overlap_used_as_filter"])
        self.assertFalse(result["raw_and_filtered_universes_modified"])
        self.assertEqual(
            result["filter_scenarios_not_a_final_panel"]["require_fully_inside_footprint"], 1
        )
        self.assertEqual(
            result["filter_scenarios_not_a_final_panel"]["require_any_overlap_with_footprint"], 2
        )
        self.assertEqual(
            result["filter_scenarios_not_a_final_panel"]["require_outside_footprint_only"], 1
        )
        self.assertEqual(
            result["filter_scenarios_not_a_final_panel"]["current_filtered_panel_no_footprint_filter"], 3
        )

    def test_missing_coverage_data_defaults_to_outside(self):
        panel_path = self._write("filtered.tsv", [
            "marker_id\tchrom\tstart\tend\tmotif\tref_repeat_length\tkind",
            "m1\tchr1\t100\t110\tA\t10\tmononucleotide",
        ])
        coverage_path = self._write("coverage.tsv", [])
        gnomad_path = self._write("gnomad.tsv", [])
        result = bfa.analyze(panel_path, coverage_path, gnomad_path)
        self.assertEqual(result["breakdown_outside_footprint"]["n"], 1)


if __name__ == "__main__":
    unittest.main()
