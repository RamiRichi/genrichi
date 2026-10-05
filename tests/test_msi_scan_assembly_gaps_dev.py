"""
SYNTHETIC, DEVELOPMENT-ONLY tests for dev/msi_redesign/scan_assembly_gaps.py.

Run with:
    python -m unittest tests.test_msi_scan_assembly_gaps_dev -v
"""

import sys
import unittest
from pathlib import Path

DEV_DIR = Path(__file__).resolve().parents[1] / "dev" / "msi_redesign"
sys.path.insert(0, str(DEV_DIR))

import scan_assembly_gaps as sag  # noqa: E402


class TestFindNRuns(unittest.TestCase):
    def test_single_run(self):
        seq = "ACGT" + "N" * 5 + "ACGT"
        self.assertEqual(list(sag.find_n_runs(seq)), [(4, 9)])

    def test_no_run(self):
        self.assertEqual(list(sag.find_n_runs("ACGTACGT")), [])

    def test_multiple_runs(self):
        seq = "N" * 3 + "ACGT" + "N" * 2
        self.assertEqual(list(sag.find_n_runs(seq)), [(0, 3), (7, 9)])

    def test_single_n_counts(self):
        self.assertEqual(list(sag.find_n_runs("ACNGT")), [(2, 3)])


if __name__ == "__main__":
    unittest.main()
