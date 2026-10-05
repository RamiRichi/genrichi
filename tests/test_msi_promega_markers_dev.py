"""
SYNTHETIC, DEVELOPMENT-ONLY tests for dev/msi_redesign/promega_markers.py's
pure helper functions (revcomp, find_homopolymer_runs, locate_primer).

These do NOT re-verify the real GRCh38 coordinates (that verification run
against the real hg38.fa is reported separately, with its own manifest/hash
output) -- they only prove the primer-location and homopolymer-finding
logic behaves correctly on fabricated sequences.

Run with:
    python -m unittest tests.test_msi_promega_markers_dev -v
"""

import sys
import unittest
from pathlib import Path

DEV_DIR = Path(__file__).resolve().parents[1] / "dev" / "msi_redesign"
sys.path.insert(0, str(DEV_DIR))

import promega_markers as pm  # noqa: E402


class TestRevComp(unittest.TestCase):
    def test_simple_sequence(self):
        self.assertEqual(pm.revcomp("ACGT"), "ACGT")

    def test_asymmetric_sequence(self):
        self.assertEqual(pm.revcomp("AACCATGCTTGCAAACCACT"), "AGTGGTTTGCAAGCATGGTT")


class TestFindHomopolymerRuns(unittest.TestCase):
    def test_finds_run_at_min_length(self):
        seq = "GC" + "A" * 6 + "GC"
        hits = pm.find_homopolymer_runs(seq, min_len=6)
        self.assertEqual(hits, [(2, 8, "A", 6)])

    def test_below_min_length_not_reported(self):
        seq = "GC" + "A" * 5 + "GC"
        hits = pm.find_homopolymer_runs(seq, min_len=6)
        self.assertEqual(hits, [])

    def test_picks_longest_when_multiple_runs_present(self):
        seq = "A" * 6 + "GG" + "T" * 20
        hits = pm.find_homopolymer_runs(seq, min_len=6)
        longest = max(hits, key=lambda h: h[3])
        self.assertEqual(longest, (8, 28, "T", 20))


class TestLocatePrimer(unittest.TestCase):
    def test_forward_orientation_found(self):
        primer = "ACGTACGT"
        seq = "TTTT" + primer + "TTTT"
        hit = pm.locate_primer(seq, primer)
        self.assertEqual(hit, (4, 12, "forward_as_given"))

    def test_reverse_complement_orientation_found(self):
        primer = "AAACCCTTT"  # not a revcomp palindrome, unlike ACGTACGT/AACCGGTT
        seq = "GGGG" + pm.revcomp(primer) + "GGGG"
        hit = pm.locate_primer(seq, primer)
        self.assertEqual(hit, (4, 13, "reverse_complement"))

    def test_absent_primer_returns_none(self):
        hit = pm.locate_primer("TTTTTTTTTTTTTTTT", "ACGTACGT")
        self.assertIsNone(hit)


class TestVerificationHashDeterminism(unittest.TestCase):
    def test_same_payload_same_hash(self):
        self.assertEqual(pm.sha256_text("chr1\t10\t20\tA\t10"), pm.sha256_text("chr1\t10\t20\tA\t10"))

    def test_different_payload_different_hash(self):
        self.assertNotEqual(pm.sha256_text("chr1\t10\t20\tA\t10"), pm.sha256_text("chr1\t10\t20\tA\t11"))


if __name__ == "__main__":
    unittest.main()
