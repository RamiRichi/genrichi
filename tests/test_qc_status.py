"""Unit tests for workflow/scripts/qc_status.py (GenRichi Phase 5.2)."""

import sys
import unittest
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import qc_status  # noqa: E402


class TestMsiQcStatus(unittest.TestCase):
    def test_zero_variants_is_insufficient(self):
        self.assertEqual(qc_status.msi_qc_status(0), qc_status.INSUFFICIENT_DATA)

    def test_none_is_insufficient(self):
        self.assertEqual(qc_status.msi_qc_status(None), qc_status.INSUFFICIENT_DATA)

    def test_one_variant_is_valid(self):
        self.assertEqual(qc_status.msi_qc_status(1), qc_status.VALID)

    def test_many_variants_is_valid(self):
        self.assertEqual(qc_status.msi_qc_status(500), qc_status.VALID)


class TestCnvQcStatus(unittest.TestCase):
    def test_zero_overlap_is_insufficient(self):
        self.assertEqual(qc_status.cnv_qc_status(0), qc_status.INSUFFICIENT_DATA)

    def test_none_is_insufficient(self):
        self.assertEqual(qc_status.cnv_qc_status(None), qc_status.INSUFFICIENT_DATA)

    def test_one_region_is_valid(self):
        self.assertEqual(qc_status.cnv_qc_status(1), qc_status.VALID)

    def test_many_regions_is_valid(self):
        self.assertEqual(qc_status.cnv_qc_status(901), qc_status.VALID)


if __name__ == "__main__":
    unittest.main()
