"""
SYNTHETIC, DEVELOPMENT-ONLY tests for
dev/msi_redesign/build_sequencing_strategy_views.py.

*** Fabricated marker/coverage data only. The real view files and manifest ***
*** come from a separate real run against the real master filtered panel ***
*** and the real footprint coverage data, with their own provenance.

Run with:
    python -m unittest tests.test_msi_build_sequencing_strategy_views_dev -v
"""

import sys
import tempfile
import unittest
from pathlib import Path

DEV_DIR = Path(__file__).resolve().parents[1] / "dev" / "msi_redesign"
sys.path.insert(0, str(DEV_DIR))

import build_sequencing_strategy_views as bsv  # noqa: E402


class TestSelectFullyInside(unittest.TestCase):
    def test_selects_only_full_coverage(self):
        coverage = {"m1": 1.0, "m2": 0.5, "m3": 0.0, "m4": 0.9999999}
        self.assertEqual(bsv.select_fully_inside(coverage), {"m1", "m4"})

    def test_empty_coverage_gives_empty_set(self):
        self.assertEqual(bsv.select_fully_inside({}), set())

    def test_custom_threshold(self):
        coverage = {"m1": 0.95, "m2": 0.5}
        self.assertEqual(bsv.select_fully_inside(coverage, min_fraction=0.9), {"m1"})


class TestBuildView(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def test_writes_only_member_rows_preserving_columns(self):
        master_path = Path(self.tmp.name) / "master.tsv"
        master_path.write_text(
            "marker_id\tchrom\tstart\tend\n"
            "m1\tchr1\t100\t110\n"
            "m2\tchr1\t500\t510\n"
            "m3\tchr1\t900\t910\n"
        )
        out_path = Path(self.tmp.name) / "view.tsv"
        n_written = bsv.build_view(str(master_path), {"m1", "m3"}, str(out_path))

        self.assertEqual(n_written, 2)
        lines = out_path.read_text().splitlines()
        self.assertEqual(lines[0], "marker_id\tchrom\tstart\tend")
        self.assertEqual(lines[1], "m1\tchr1\t100\t110")
        self.assertEqual(lines[2], "m3\tchr1\t900\t910")
        self.assertEqual(len(lines), 3)

    def test_empty_membership_writes_header_only(self):
        master_path = Path(self.tmp.name) / "master.tsv"
        master_path.write_text("marker_id\tchrom\n" "m1\tchr1\n")
        out_path = Path(self.tmp.name) / "view.tsv"
        n_written = bsv.build_view(str(master_path), set(), str(out_path))
        self.assertEqual(n_written, 0)
        self.assertEqual(out_path.read_text().splitlines(), ["marker_id\tchrom"])

    def test_master_file_is_not_modified(self):
        master_path = Path(self.tmp.name) / "master.tsv"
        original = "marker_id\tchrom\nm1\tchr1\nm2\tchr1\n"
        master_path.write_text(original)
        out_path = Path(self.tmp.name) / "view.tsv"
        bsv.build_view(str(master_path), {"m1"}, str(out_path))
        self.assertEqual(master_path.read_text(), original)


if __name__ == "__main__":
    unittest.main()
