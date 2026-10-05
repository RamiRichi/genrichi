"""
SYNTHETIC, DEVELOPMENT-ONLY tests for dev/msi_redesign/bam_extract.py.

*** All reads/CIGAR strings/BAM records here are fabricated for these ***
*** tests. None of this is real sequencing data. These tests prove the ***
*** read-filtering and CIGAR-length LOGIC is correct on deliberately ***
*** constructed examples, including one small synthetic BAM built with ***
*** real `samtools` to exercise the actual I/O path -- still entirely ***
*** synthetic data, not biological evidence.

Run with:
    python -m unittest tests.test_msi_bam_extract_dev -v
"""

import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

DEV_DIR = Path(__file__).resolve().parents[1] / "dev" / "msi_redesign"
sys.path.insert(0, str(DEV_DIR))

import bam_extract as be  # noqa: E402

SAMTOOLS = shutil.which("samtools")


def rec(flag=0, pos=1, mapq=60, cigar="30M"):
    return be.SamRecord(flag=flag, pos=pos, mapq=mapq, cigar=cigar)


class TestCigarParsing(unittest.TestCase):
    def test_simple_cigar(self):
        self.assertEqual(be.parse_cigar("76M"), [(76, "M")])

    def test_mixed_cigar(self):
        self.assertEqual(be.parse_cigar("10S60M2I10M"), [(10, "S"), (60, "M"), (2, "I"), (10, "M")])

    def test_star_cigar_is_empty(self):
        self.assertEqual(be.parse_cigar("*"), [])


class TestReadFilters(unittest.TestCase):
    def test_a_clean_primary_read_passes(self):
        self.assertTrue(be.passes_read_filters(rec(flag=0, mapq=60, cigar="30M")))

    def test_unmapped_read_is_excluded(self):
        self.assertFalse(be.passes_read_filters(rec(flag=0x4)))

    def test_secondary_alignment_is_excluded(self):
        self.assertFalse(be.passes_read_filters(rec(flag=0x100)))

    def test_supplementary_alignment_is_excluded(self):
        self.assertFalse(be.passes_read_filters(rec(flag=0x800)))

    def test_duplicate_is_excluded(self):
        self.assertFalse(be.passes_read_filters(rec(flag=0x400)))

    def test_low_mapq_is_excluded(self):
        self.assertFalse(be.passes_read_filters(rec(mapq=10), min_mapq=20))

    def test_star_cigar_read_is_excluded(self):
        self.assertFalse(be.passes_read_filters(rec(cigar="*")))


class TestObservedRepeatLength(unittest.TestCase):
    """locus_start/locus_end are 0-based half-open, matching scan_markers.py."""

    def test_clean_match_no_indel_returns_locus_width(self):
        # Read covers ref [0, 30) with plain match; locus is [10, 20) ->
        # width 10, with 10bp flank on each side.
        r = rec(pos=1, cigar="30M")  # SAM POS=1 -> 0-based ref_pos starts at 0
        length = be.observed_repeat_length(r, locus_start=10, locus_end=20, min_flank=5)
        self.assertEqual(length, 10)

    def test_insertion_inside_locus_increases_observed_length(self):
        # ref [0,10) M, then 3bp insertion at ref_pos=10 (inside locus
        # [10,20) at its very start), then [10,30) M covering the rest of
        # the locus and flank.
        r = rec(pos=1, cigar="10M3I20M")
        length = be.observed_repeat_length(r, locus_start=10, locus_end=20, min_flank=5)
        # locus width 10 + 3bp insertion attributed to it = 13
        self.assertEqual(length, 13)

    def test_deletion_inside_locus_decreases_observed_length(self):
        # ref [0,15) M, 3bp deletion at ref_pos 15..18 (inside locus
        # [10,20)), then ref [18,30) M.
        r = rec(pos=1, cigar="15M3D12M")
        length = be.observed_repeat_length(r, locus_start=10, locus_end=20, min_flank=5)
        # locus width 10, minus 3bp deleted (no read bases) = 7
        self.assertEqual(length, 7)

    def test_insufficient_flank_before_locus_returns_none(self):
        # Read starts at ref_pos=8 (0-based), only 2bp before the locus
        # start of 10 -> below min_flank=5.
        r = rec(pos=9, cigar="30M")  # SAM POS=9 -> 0-based ref_pos=8
        length = be.observed_repeat_length(r, locus_start=10, locus_end=20, min_flank=5)
        self.assertIsNone(length)

    def test_insufficient_flank_after_locus_returns_none(self):
        # Read covers ref [0, 23) only -> 3bp after locus_end=20, below min_flank=5.
        r = rec(pos=1, cigar="23M")
        length = be.observed_repeat_length(r, locus_start=10, locus_end=20, min_flank=5)
        self.assertIsNone(length)

    def test_soft_clip_before_match_still_counts_aligned_flank(self):
        # 5bp soft-clip (no reference contribution), then 30M covering
        # ref [0,30) -- the clip must not be mistaken for aligned flank,
        # but the subsequent M still gives full flank on both sides.
        r = rec(pos=1, cigar="5S30M")
        length = be.observed_repeat_length(r, locus_start=10, locus_end=20, min_flank=5)
        self.assertEqual(length, 10)

    def test_read_entirely_outside_locus_returns_none(self):
        r = rec(pos=1, cigar="5M")  # ref [0,5), locus is [10,20)
        length = be.observed_repeat_length(r, locus_start=10, locus_end=20, min_flank=5)
        self.assertIsNone(length)

    def test_empty_cigar_returns_none(self):
        r = rec(cigar="*")
        length = be.observed_repeat_length(r, locus_start=10, locus_end=20, min_flank=5)
        self.assertIsNone(length)


@unittest.skipUnless(SAMTOOLS, "samtools not on PATH")
class TestSyntheticBamEndToEnd(unittest.TestCase):
    """Builds one small, fully synthetic BAM (fictional contig, fabricated
    reads) with real `samtools`, to exercise the actual samtools-view I/O
    path -- not just the pure CIGAR function. Still entirely synthetic:
    the contig name/length and every read are invented for this test."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="msi-dev-bam-"))
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def _build_bam(self, reads):
        sam_path = self.tmp / "synthetic.sam"
        bam_path = self.tmp / "synthetic.bam"
        header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:synthetic_contig\tLN:1000\n"
        lines = [header]
        for i, (pos, mapq, cigar, flag) in enumerate(reads):
            # minimal valid SAM record; SEQ/QUAL as placeholders sized to CIGAR's read-consuming length
            read_len = sum(n for n, op in be.parse_cigar(cigar) if op in be.READ_CONSUMING)
            seq = "A" * max(read_len, 1)
            qual = "I" * max(read_len, 1)
            lines.append(f"read{i}\t{flag}\tsynthetic_contig\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t{qual}\n")
        sam_path.write_text("".join(lines))
        subprocess.run([SAMTOOLS, "view", "-b", "-o", str(bam_path), str(sam_path)], check=True, capture_output=True)
        subprocess.run([SAMTOOLS, "sort", "-o", str(bam_path), str(bam_path)], check=True, capture_output=True)
        subprocess.run([SAMTOOLS, "index", str(bam_path)], check=True, capture_output=True)
        return bam_path

    def test_real_samtools_view_returns_expected_reads_and_lengths(self):
        # 12 clean reads covering ref [0,30) with sufficient flank around
        # locus [10,20); 3 reads with MAPQ below threshold (excluded); 2
        # reads flagged as duplicates (excluded).
        reads = (
            [(1, 60, "30M", 0)] * 12
            + [(1, 5, "30M", 0)] * 3          # low MAPQ
            + [(1, 60, "30M", 0x400)] * 2     # duplicate flag
        )
        bam = self._build_bam(reads)
        lengths = be.locus_lengths(SAMTOOLS, str(bam), "synthetic_contig", 10, 20, min_mapq=20, min_flank=5)
        self.assertEqual(len(lengths), 12)
        self.assertTrue(all(length == 10 for length in lengths))

    def test_real_samtools_view_reflects_a_fabricated_expansion(self):
        # 10 "normal-like" reads with no indel at the locus, and 10
        # "tumor-like" reads with a fabricated 2bp insertion inside the
        # locus -- proves the real I/O path carries the CIGAR-derived
        # length difference through correctly, on synthetic data only.
        reads = [(1, 60, "30M", 0)] * 10 + [(1, 60, "10M2I20M", 0)] * 10
        bam = self._build_bam(reads)
        lengths = be.locus_lengths(SAMTOOLS, str(bam), "synthetic_contig", 10, 20, min_mapq=20, min_flank=5)
        self.assertEqual(sorted(lengths), sorted([10] * 10 + [12] * 10))


if __name__ == "__main__":
    unittest.main()
