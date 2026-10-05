"""
SYNTHETIC, DEVELOPMENT-ONLY tests for dev/msi_redesign/scan_markers_v2.py
(the genome-wide, cancer-panel-independent marker scanner).

*** All sequences here are fabricated strings, not real reference genome ***
*** data. They prove the regex-based scanning/normalization LOGIC matches ***
*** the intended semantics (maximal runs, minimal-period motif reduction, ***
*** homopolymer-STR overlap exclusion) -- they say nothing about the real ***
*** hg38 whole-genome scan's biological suitability, which is reported ***
*** separately from the real scan output.

Run with:
    python -m unittest tests.test_msi_scan_markers_v2_dev -v
"""

import sys
import tempfile
import unittest
from pathlib import Path

DEV_DIR = Path(__file__).resolve().parents[1] / "dev" / "msi_redesign"
sys.path.insert(0, str(DEV_DIR))

import scan_markers_v2 as sm2  # noqa: E402


class TestMinimalPeriod(unittest.TestCase):
    def test_already_primitive_motif_unchanged(self):
        self.assertEqual(sm2.minimal_period("AC"), "AC")
        self.assertEqual(sm2.minimal_period("ACG"), "ACG")

    def test_reduces_dinucleotide_disguised_as_tetranucleotide(self):
        self.assertEqual(sm2.minimal_period("ATAT"), "AT")

    def test_reduces_trinucleotide_disguised_as_hexanucleotide(self):
        self.assertEqual(sm2.minimal_period("ACGACG"), "ACG")

    def test_reduces_to_period_one_for_pure_homopolymer_text(self):
        self.assertEqual(sm2.minimal_period("AAAA"), "A")

    def test_five_base_motif_with_no_smaller_period_unchanged(self):
        self.assertEqual(sm2.minimal_period("ACGTA"), "ACGTA")


class TestHomopolymerScan(unittest.TestCase):
    def test_finds_run_at_minimum_length(self):
        seq = "GC" + "A" * 8 + "GC"
        hits = list(sm2.find_homopolymer_runs(seq, min_len=8, max_len=50))
        self.assertEqual(hits, [(2, 10, "A", 8)])

    def test_run_below_minimum_not_reported(self):
        seq = "GC" + "A" * 7 + "GC"
        hits = list(sm2.find_homopolymer_runs(seq, min_len=8, max_len=50))
        self.assertEqual(hits, [])

    def test_run_above_maximum_excluded(self):
        seq = "GC" + "A" * 60 + "GC"
        hits = list(sm2.find_homopolymer_runs(seq, min_len=8, max_len=50))
        self.assertEqual(hits, [])

    def test_run_at_maximum_length_included(self):
        seq = "GC" + "T" * 50 + "GC"
        hits = list(sm2.find_homopolymer_runs(seq, min_len=8, max_len=50))
        self.assertEqual(hits, [(2, 52, "T", 50)])

    def test_two_separate_runs_both_found(self):
        seq = "A" * 8 + "GGGG" + "T" * 10
        hits = list(sm2.find_homopolymer_runs(seq, min_len=8, max_len=50))
        self.assertEqual(hits, [(0, 8, "A", 8), (12, 22, "T", 10)])


class TestShortTandemRepeatScan(unittest.TestCase):
    def test_finds_dinucleotide_repeat(self):
        seq = "GC" + "AT" * 6 + "GC"
        hits = list(sm2.find_short_tandem_repeats(seq, min_units=5))
        self.assertEqual(hits, [(2, 14, "AT", 6)])

    def test_below_minimum_units_not_reported(self):
        seq = "GC" + "AT" * 4 + "GC"
        hits = list(sm2.find_short_tandem_repeats(seq, min_units=5))
        self.assertEqual(hits, [])

    def test_finds_hexanucleotide_repeat_new_in_v2(self):
        # v1 capped motif length at 4bp; v2 extends to 6bp per the
        # MSIsensor-pro precedent -- this motif length was unreachable in v1.
        # Flanking bases are chosen so they do NOT coincide with a rotation
        # of the repeat motif (which would ambiguously shift the boundary).
        seq = "TT" + "ACGTAC" * 5 + "TT"
        hits = list(sm2.find_short_tandem_repeats(seq, min_units=5))
        self.assertEqual(hits, [(2, 32, "ACGTAC", 5)])

    def test_tetranucleotide_disguised_dinucleotide_is_normalized(self):
        # "ATAT" x N is really "AT" x 2N -- must be reported with the
        # minimal period, not the raw captured group.
        seq = "GC" + "AT" * 10 + "GC"
        hits = list(sm2.find_short_tandem_repeats(seq, min_units=5))
        self.assertEqual(hits, [(2, 22, "AT", 10)])

    def test_pure_homopolymer_text_is_not_reported_as_str(self):
        seq = "A" * 12
        hits = list(sm2.find_short_tandem_repeats(seq, min_units=5))
        self.assertEqual(hits, [])


class TestScanChromosomeOverlapExclusion(unittest.TestCase):
    def test_str_overlapping_homopolymer_span_is_excluded(self):
        # A homopolymer run of "A"x10 could also register as a spurious
        # 2bp "AA" STR match at the same span -- scan_chromosome must keep
        # only the homopolymer classification, not both.
        seq = "GC" + "A" * 10 + "GC"
        markers = sm2.scan_chromosome("chrTest", seq)
        kinds = [m["kind"] for m in markers]
        self.assertEqual(kinds.count("mononucleotide"), 1)
        self.assertNotIn("str", kinds)

    def test_independent_homopolymer_and_str_both_reported(self):
        seq = "GC" + "A" * 10 + "GC" + "TG" * 6 + "GC"
        markers = sm2.scan_chromosome("chrTest", seq)
        kinds = sorted(m["kind"] for m in markers)
        self.assertEqual(kinds, ["mononucleotide", "str"])

    def test_marker_ids_are_unique_and_positionally_derived(self):
        seq = "GC" + "A" * 10 + "GC" + "TG" * 6 + "GC"
        markers = sm2.scan_chromosome("chrTest", seq)
        ids = [m["marker_id"] for m in markers]
        self.assertEqual(len(ids), len(set(ids)))
        for m in markers:
            self.assertIn(str(m["start"]), m["marker_id"])
            self.assertIn(str(m["end"]), m["marker_id"])


class TestDeterminism(unittest.TestCase):
    def test_scanning_same_sequence_twice_gives_identical_results(self):
        seq = "GC" + "A" * 10 + "TAG" * 6 + "N" + "CAGCAG" * 5 + "TG" * 7
        first = sm2.scan_chromosome("chrTest", seq)
        second = sm2.scan_chromosome("chrTest", seq)
        self.assertEqual(first, second)

    def test_marker_file_hash_reproducible(self):
        with tempfile.TemporaryDirectory() as tmp:
            p1 = Path(tmp) / "a.tsv"
            p2 = Path(tmp) / "b.tsv"
            markers = [{"marker_id": "MS2_chr1_10_20", "chrom": "chr1", "start": 10, "end": 20,
                        "motif": "A", "ref_repeat_length": 10, "kind": "mononucleotide"}]
            sm2.write_marker_tsv(markers, str(p1), "GRCh38")
            sm2.write_marker_tsv(markers, str(p2), "GRCh38")
            self.assertEqual(sm2.sha256_file(str(p1)), sm2.sha256_file(str(p2)))

    def test_marker_file_hash_changes_if_marker_changes(self):
        with tempfile.TemporaryDirectory() as tmp:
            p1 = Path(tmp) / "a.tsv"
            p2 = Path(tmp) / "b.tsv"
            base = {"marker_id": "MS2_chr1_10_20", "chrom": "chr1", "start": 10, "end": 20,
                    "motif": "A", "ref_repeat_length": 10, "kind": "mononucleotide"}
            changed = dict(base, ref_repeat_length=11)
            sm2.write_marker_tsv([base], str(p1), "GRCh38")
            sm2.write_marker_tsv([changed], str(p2), "GRCh38")
            self.assertNotEqual(sm2.sha256_file(str(p1)), sm2.sha256_file(str(p2)))


class TestOverlapsAnyHelper(unittest.TestCase):
    def test_detects_overlap(self):
        self.assertTrue(sm2._overlaps_any(5, 15, [(10, 20)]))

    def test_detects_no_overlap(self):
        self.assertFalse(sm2._overlaps_any(5, 10, [(10, 20)]))  # half-open, touching not overlapping

    def test_works_with_multiple_intervals(self):
        intervals = [(0, 5), (10, 20), (30, 40)]
        self.assertTrue(sm2._overlaps_any(18, 25, intervals))
        self.assertFalse(sm2._overlaps_any(21, 29, intervals))


if __name__ == "__main__":
    unittest.main()
