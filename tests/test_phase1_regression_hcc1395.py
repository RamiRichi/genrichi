"""
GenRichi Phase 5.3 -- Phase 1 scientific/regression tests.

Distinct from Phase 5.1 (workflow/scripts/provenance.py) and Phase 5.2
(workflow/scripts/preflight.py, qc_status.py) infrastructure tests: this
file locks in the actual, documented HCC1395 Phase 1 baseline scientific
result (docs/PHASE1_BASELINE.md), so a future change upstream in the
pipeline cannot silently drift from the already-validated baseline without
an automated test failing.

WHAT THIS TEST DOES AND DOES NOT DO -- read this before changing anything:

  * TestHCC1395FrozenBaselineRegression is a regression of a KNOWN, FROZEN
    output snapshot: tests/fixtures/phase1_hcc1395/somatic_variants.tsv, a
    verbatim, byte-identical copy (verified via md5sum against the live
    execution host during this phase's implementation) of the real, already
    -generated HCC1395_demo.somatic_variants.tsv produced by an actual,
    documented pipeline run (docs/PHASE1_BASELINE.md). It is NOT raw FASTQ
    or reference data -- it is a ~900-byte table of already-computed variant
    calls. HCC1395 is a public, de-identified benchmark cell line (Griffith
    Lab / SEQC2), not a real patient, and this single TP53 hotspot row is
    already extensively published in the public literature about this cell
    line.
  * It does NOT invoke Snakemake, GATK, VEP, conda, or any external tool,
    and does NOT download or regenerate anything. It is fully deterministic
    and runnable on any machine with only the Python standard library.
  * Because of that, it CANNOT catch a regression introduced upstream of
    this fixture (alignment, calling, filtering, annotation) unless a
    human/agent deliberately re-runs the real pipeline and refreshes this
    fixture from a new, reviewed run. True end-to-end verification remains
    the manual procedure documented in docs/PHASE1_BASELINE.md.
  * TestHCC1395LiveOutputIfPresent is a SEPARATE, OPTIONAL class that
    additionally cross-checks a LIVE output file if one already happens to
    exist on disk (e.g. on the validated execution host after a real
    pipeline run). It still never invokes the pipeline itself -- it only
    reads a file if present -- and is skipped everywhere else, including
    this repository checkout.

Run with:
    python -m unittest tests.test_phase1_regression_hcc1395 -v
"""

import csv
import unittest
from pathlib import Path

FIXTURE_DIR = Path(__file__).resolve().parent / "fixtures" / "phase1_hcc1395"
FIXTURE_TSV = FIXTURE_DIR / "somatic_variants.tsv"

# ── Documented baseline (docs/PHASE1_BASELINE.md, "Validation Checkpoints") ──
# Do not change these values without updating that document and refreshing
# the fixture above from a new, reviewed pipeline run.
EXPECTED_PASS_VARIANT_COUNT = 1
EXPECTED_GENE = "TP53"
EXPECTED_CHROM = "chr17"
EXPECTED_POS = "7675088"
EXPECTED_REF = "C"
EXPECTED_ALT = "T"
EXPECTED_HGVSC = "c.524G>A"
EXPECTED_HGVSP = "p.Arg175His"
EXPECTED_CONSEQUENCE = "missense_variant"
EXPECTED_DEPTH = "83"        # DP
EXPECTED_ALT_READS = "82"    # AD (alt allele)
EXPECTED_VAF = "0.988"
EXPECTED_FILTER = "PASS"
EXCLUDED_GENE = "BRCA1"      # documented as correctly excluded by the Phase 1 BED


def _load_rows(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class TestHCC1395FrozenBaselineRegression(unittest.TestCase):
    """Regression of a KNOWN, frozen output snapshot -- see module docstring."""

    @classmethod
    def setUpClass(cls):
        cls.rows = _load_rows(FIXTURE_TSV)

    def _tp53_row(self):
        tp53_rows = [r for r in self.rows if r["gene"] == EXPECTED_GENE]
        self.assertEqual(len(tp53_rows), 1, "expected exactly one TP53 row in the frozen baseline")
        return tp53_rows[0]

    def test_fixture_is_a_small_derived_artifact_not_raw_data(self):
        self.assertTrue(FIXTURE_TSV.exists())
        self.assertLess(FIXTURE_TSV.stat().st_size, 10_000)

    def test_pass_variant_count(self):
        pass_rows = [r for r in self.rows if r["filter"] == EXPECTED_FILTER]
        self.assertEqual(len(pass_rows), EXPECTED_PASS_VARIANT_COUNT)

    def test_tp53_variant_identity(self):
        row = self._tp53_row()
        self.assertEqual(row["chrom"], EXPECTED_CHROM)
        self.assertEqual(row["pos"], EXPECTED_POS)
        self.assertEqual(row["ref"], EXPECTED_REF)
        self.assertEqual(row["alt"], EXPECTED_ALT)
        self.assertEqual(row["gene"], EXPECTED_GENE)

    def test_tp53_hgvs_and_consequence(self):
        row = self._tp53_row()
        self.assertEqual(row["hgvsc"], EXPECTED_HGVSC)
        self.assertEqual(row["hgvsp"], EXPECTED_HGVSP)
        self.assertEqual(row["consequence"], EXPECTED_CONSEQUENCE)

    def test_tp53_pass_status(self):
        self.assertEqual(self._tp53_row()["filter"], EXPECTED_FILTER)

    def test_tp53_vaf_dp_ad(self):
        row = self._tp53_row()
        self.assertEqual(row["depth"], EXPECTED_DEPTH)
        self.assertEqual(row["alt_reads"], EXPECTED_ALT_READS)
        self.assertEqual(row["vaf"], EXPECTED_VAF)
        # internal consistency: alt_reads/depth must match the reported VAF
        self.assertAlmostEqual(
            int(row["alt_reads"]) / int(row["depth"]), float(row["vaf"]), places=3
        )

    def test_brca1_intronic_variant_excluded(self):
        genes = {r["gene"] for r in self.rows}
        self.assertNotIn(EXCLUDED_GENE, genes)


class TestHCC1395LiveOutputIfPresent(unittest.TestCase):
    """
    OPTIONAL cross-check against a LIVE, already-generated pipeline output
    file, only if one happens to already exist on disk. Never invokes
    Snakemake/GATK/VEP -- only reads a file if present. Skipped everywhere
    a completed HCC1395 run isn't already materialized (including this
    repository checkout, since `results/` is gitignored).
    """

    LIVE_PATH = (
        Path(__file__).resolve().parents[1]
        / "results" / "HCC1395_demo" / "annotation" / "HCC1395_demo.somatic_variants.tsv"
    )

    @unittest.skipUnless(LIVE_PATH.exists(), "no live HCC1395 pipeline output present on this machine")
    def test_live_output_matches_documented_baseline(self):
        rows = _load_rows(self.LIVE_PATH)
        pass_rows = [r for r in rows if r["filter"] == EXPECTED_FILTER]
        self.assertEqual(len(pass_rows), EXPECTED_PASS_VARIANT_COUNT)

        tp53_rows = [r for r in rows if r["gene"] == EXPECTED_GENE]
        self.assertEqual(len(tp53_rows), 1)
        row = tp53_rows[0]
        self.assertEqual(row["hgvsc"], EXPECTED_HGVSC)
        self.assertEqual(row["hgvsp"], EXPECTED_HGVSP)
        self.assertEqual(row["depth"], EXPECTED_DEPTH)
        self.assertEqual(row["alt_reads"], EXPECTED_ALT_READS)
        self.assertEqual(row["vaf"], EXPECTED_VAF)
        self.assertEqual(row["filter"], EXPECTED_FILTER)

        genes = {r["gene"] for r in rows}
        self.assertNotIn(EXCLUDED_GENE, genes)


if __name__ == "__main__":
    unittest.main()
