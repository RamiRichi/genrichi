"""
SYNTHETIC, DEVELOPMENT-ONLY tests for dev/msi_redesign/analyze_candidate_universe.py.

*** All marker rows here are fabricated for these tests, not real scan ***
*** output. They prove the descriptive-statistics aggregation logic is ***
*** correct -- the real candidate-universe numbers used in the filtering ***
*** design report were produced by a separate real run against the real ***
*** Tier 2 marker file, reported with their own provenance.

Run with:
    python -m unittest tests.test_msi_analyze_candidate_universe_dev -v
"""

import csv
import sys
import tempfile
import unittest
from pathlib import Path

DEV_DIR = Path(__file__).resolve().parents[1] / "dev" / "msi_redesign"
sys.path.insert(0, str(DEV_DIR))

import analyze_candidate_universe as acu  # noqa: E402

FIELDS = ["marker_id", "chrom", "start", "end", "motif", "ref_repeat_length", "reference_build", "kind"]


def write_tsv(path, rows):
    with open(path, "w", newline="\n", encoding="utf-8") as fh:
        fh.write("\t".join(FIELDS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[f]) for f in FIELDS) + "\n")


class TestAnalyzeCandidateUniverse(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def _row(self, marker_id, chrom, start, end, motif, length, kind):
        return {"marker_id": marker_id, "chrom": chrom, "start": start, "end": end,
                "motif": motif, "ref_repeat_length": length, "reference_build": "GRCh38", "kind": kind}

    def test_basic_counts_and_kind_split(self):
        rows = [
            self._row("m1", "chr1", 100, 108, "A", 8, "mononucleotide"),
            self._row("m2", "chr1", 200, 210, "AC", 5, "str"),
        ]
        path = Path(self.tmp.name) / "markers.tsv"
        write_tsv(path, rows)
        stats = acu.analyze(str(path))
        self.assertEqual(stats["n_total"], 2)
        self.assertEqual(stats["kind_counts"], {"mononucleotide": 1, "str": 1})
        self.assertEqual(stats["chrom_counts"], {"chr1": 2})

    def test_mononucleotide_length_histogram(self):
        rows = [
            self._row("m1", "chr1", 0, 8, "A", 8, "mononucleotide"),
            self._row("m2", "chr1", 1000, 1010, "T", 10, "mononucleotide"),
            self._row("m3", "chr1", 2000, 2010, "A", 10, "mononucleotide"),
        ]
        path = Path(self.tmp.name) / "markers.tsv"
        write_tsv(path, rows)
        stats = acu.analyze(str(path))
        self.assertEqual(stats["mononucleotide_length_histogram"], {8: 1, 10: 2})

    def test_cumulative_count_by_max_length_cap_is_monotonic(self):
        rows = [self._row(f"m{i}", "chr1", i * 100, i * 100 + 8 + i, "A", 8 + i, "mononucleotide")
                for i in range(5)]
        path = Path(self.tmp.name) / "markers.tsv"
        write_tsv(path, rows)
        stats = acu.analyze(str(path))
        caps = stats["mononucleotide_cumulative_count_by_max_length_cap"]
        values = [caps[k] for k in sorted(caps, key=int)]
        self.assertEqual(values, sorted(values))  # non-decreasing as cap widens

    def test_str_motif_length_and_units_histograms(self):
        rows = [
            self._row("m1", "chr1", 0, 10, "AC", 5, "str"),
            self._row("m2", "chr1", 100, 118, "ACG", 6, "str"),
        ]
        path = Path(self.tmp.name) / "markers.tsv"
        write_tsv(path, rows)
        stats = acu.analyze(str(path))
        self.assertEqual(stats["str_motif_length_histogram"], {2: 1, 3: 1})
        self.assertEqual(stats["str_units_histogram"], {5: 1, 6: 1})

    def test_nearest_neighbor_gap_detects_close_loci(self):
        # two loci 5bp apart on the same chromosome -> both flagged "within 50bp"
        rows = [
            self._row("m1", "chr1", 100, 110, "A", 10, "mononucleotide"),
            self._row("m2", "chr1", 115, 125, "T", 10, "mononucleotide"),
        ]
        path = Path(self.tmp.name) / "markers.tsv"
        write_tsv(path, rows)
        stats = acu.analyze(str(path))
        self.assertEqual(stats["n_loci_within_50bp_of_another_locus"], 2)

    def test_isolated_locus_not_flagged_as_close(self):
        rows = [
            self._row("m1", "chr1", 100, 110, "A", 10, "mononucleotide"),
            self._row("m2", "chr1", 10000, 10010, "T", 10, "mononucleotide"),
        ]
        path = Path(self.tmp.name) / "markers.tsv"
        write_tsv(path, rows)
        stats = acu.analyze(str(path))
        self.assertEqual(stats["n_loci_within_50bp_of_another_locus"], 0)

    def test_different_chromosomes_do_not_count_as_neighbors(self):
        rows = [
            self._row("m1", "chr1", 100, 110, "A", 10, "mononucleotide"),
            self._row("m2", "chr2", 105, 115, "T", 10, "mononucleotide"),
        ]
        path = Path(self.tmp.name) / "markers.tsv"
        write_tsv(path, rows)
        stats = acu.analyze(str(path))
        self.assertEqual(stats["n_loci_within_50bp_of_another_locus"], 0)


if __name__ == "__main__":
    unittest.main()
