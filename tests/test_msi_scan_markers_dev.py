"""
SYNTHETIC, DEVELOPMENT-ONLY tests for dev/msi_redesign/scan_markers.py.

*** All sequences here are fabricated strings, not real reference genome ***
*** data. They prove the scanning/hashing LOGIC is correct and ***
*** deterministic -- they say nothing about the real hg38-derived marker ***
*** panel's biological suitability (that is reported separately from the ***
*** real scan in the Phase A development-analysis report).

Run with:
    python -m unittest tests.test_msi_scan_markers_dev -v
"""

import sys
import tempfile
import unittest
from pathlib import Path

DEV_DIR = Path(__file__).resolve().parents[1] / "dev" / "msi_redesign"
sys.path.insert(0, str(DEV_DIR))

import scan_markers as sm  # noqa: E402


class TestMononucleotideScan(unittest.TestCase):
    def test_finds_a_clean_run_at_the_minimum_length(self):
        seq = "GC" + "A" * 10 + "GC"
        hits = list(sm.find_mononucleotide_runs(seq, min_len=10))
        self.assertEqual(hits, [(2, 12, "A", 10)])

    def test_run_shorter_than_minimum_is_not_reported(self):
        seq = "GC" + "A" * 9 + "GC"
        hits = list(sm.find_mononucleotide_runs(seq, min_len=10))
        self.assertEqual(hits, [])

    def test_two_separate_runs_are_both_found(self):
        seq = "A" * 10 + "GGGG" + "T" * 12
        hits = list(sm.find_mononucleotide_runs(seq, min_len=10))
        self.assertEqual(hits, [(0, 10, "A", 10), (14, 26, "T", 12)])

    def test_ambiguous_base_breaks_a_run(self):
        seq = "A" * 6 + "N" + "A" * 6
        hits = list(sm.find_mononucleotide_runs(seq, min_len=10))
        self.assertEqual(hits, [])  # neither half reaches min_len alone


class TestShortTandemRepeatScan(unittest.TestCase):
    def test_finds_a_dinucleotide_repeat(self):
        seq = "GC" + "AT" * 6 + "GC"  # 6 units of "AT", well above min_units[2]=5
        hits = sm.find_short_tandem_repeats(seq, min_units={2: 5, 3: 5, 4: 5})
        self.assertEqual(hits, [(2, 14, "AT", 6)])

    def test_below_minimum_units_is_not_reported(self):
        seq = "GC" + "AT" * 3 + "GC"
        hits = sm.find_short_tandem_repeats(seq, min_units={2: 5, 3: 5, 4: 5})
        self.assertEqual(hits, [])

    def test_mononucleotide_span_is_excluded_from_str_scan(self):
        # A pure mononucleotide run must not also be reported as a
        # "repeat motif of length 1" by the STR scanner (STR scanner only
        # considers 2-4bp motifs and explicitly skips single-base motifs).
        seq = "A" * 12
        str_hits = sm.find_short_tandem_repeats(seq, min_units={2: 5, 3: 5, 4: 5})
        self.assertEqual(str_hits, [])

    def test_longer_motif_is_preferred_over_shorter_at_same_position(self):
        # "ACAC" repeated is a dinucleotide (AC)x6; must not be reported as
        # a spurious 4bp motif match instead, and must not double-count.
        seq = "AC" * 6
        hits = sm.find_short_tandem_repeats(seq, min_units={2: 5, 3: 5, 4: 5})
        self.assertEqual(hits, [(0, 12, "AC", 6)])


class TestDeterminism(unittest.TestCase):
    def test_scanning_the_same_sequence_twice_gives_identical_results(self):
        seq = "GC" + "A" * 10 + "TAG" * 6 + "N" + "CAGCAG" * 5
        first = (list(sm.find_mononucleotide_runs(seq)), sm.find_short_tandem_repeats(seq))
        second = (list(sm.find_mononucleotide_runs(seq)), sm.find_short_tandem_repeats(seq))
        self.assertEqual(first, second)

    def test_marker_file_hash_is_reproducible_for_identical_content(self):
        with tempfile.TemporaryDirectory() as tmp:
            p1 = Path(tmp) / "a.tsv"
            p2 = Path(tmp) / "b.tsv"
            markers = [{"marker_id": "MS_chr1_10_20", "chrom": "chr1", "start": 10, "end": 20,
                        "motif": "A", "ref_repeat_length": 10, "gene": "TESTGENE", "kind": "mononucleotide"}]
            sm.write_marker_tsv(markers, str(p1), "GRCh38")
            sm.write_marker_tsv(markers, str(p2), "GRCh38")
            self.assertEqual(sm.sha256_file(str(p1)), sm.sha256_file(str(p2)))

    def test_marker_file_hash_changes_if_a_marker_changes(self):
        with tempfile.TemporaryDirectory() as tmp:
            p1 = Path(tmp) / "a.tsv"
            p2 = Path(tmp) / "b.tsv"
            base = {"marker_id": "MS_chr1_10_20", "chrom": "chr1", "start": 10, "end": 20,
                    "motif": "A", "ref_repeat_length": 10, "gene": "TESTGENE", "kind": "mononucleotide"}
            changed = dict(base, ref_repeat_length=11)
            sm.write_marker_tsv([base], str(p1), "GRCh38")
            sm.write_marker_tsv([changed], str(p2), "GRCh38")
            self.assertNotEqual(sm.sha256_file(str(p1)), sm.sha256_file(str(p2)))


class TestPanelBedParsing(unittest.TestCase):
    def test_reads_bed4_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            bed = Path(tmp) / "panel.bed"
            bed.write_text("chr17\t100\t200\tTP53\nchr17\t300\t400\tBRCA1\n")
            rows = sm.read_panel_bed(str(bed))
            self.assertEqual(rows, [("chr17", 100, 200, "TP53"), ("chr17", 300, 400, "BRCA1")])

    def test_skips_blank_and_comment_lines(self):
        with tempfile.TemporaryDirectory() as tmp:
            bed = Path(tmp) / "panel.bed"
            bed.write_text("# header\n\nchr17\t100\t200\tTP53\n")
            rows = sm.read_panel_bed(str(bed))
            self.assertEqual(rows, [("chr17", 100, 200, "TP53")])


class TestChromSortOrder(unittest.TestCase):
    def test_karyotype_order_not_lexicographic(self):
        chroms = ["chr2", "chr10", "chr1", "chrX", "chrM"]
        ordered = sorted(chroms, key=sm.chrom_sort_key)
        self.assertEqual(ordered, ["chr1", "chr2", "chr10", "chrX", "chrM"])


if __name__ == "__main__":
    unittest.main()
