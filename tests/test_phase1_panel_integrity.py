"""
GenRichi Phase 5.3 -- Phase 1 scientific/regression tests.

Verifies the authoritative Phase 1 panel BED
(resources/panel/phase1_solid_tumor/solid_tumor_phase1_v1.bed) remains
internally consistent with its own documented audit metadata
(resources/panel/phase1_solid_tumor/phase1_audit.json) and with the counts
recorded in docs/PHASE1_BASELINE.md. Read-only: does NOT modify the BED or
its audit file.

On the SHA256 check below: no pre-existing hash for this BED is documented
anywhere in the repository (checked: phase1_audit.json,
docs/PHASE1_BASELINE.md, resources/panel/phase1_solid_tumor/README.md --
none contain a hash field). The primary checks in this file therefore rely
entirely on gene count / interval count / duplicate count / per-gene
audit-row cross-validation, all of which ARE already documented in
phase1_audit.json -- nothing here is invented.

A SECONDARY, separately-labelled content-hash check pins the SHA256 already
independently computed and verified twice in this repository's history:
once by workflow/scripts/provenance.py during the Phase 5.1 HCC1395
regression run (embedded in that run's provenance JSON), and cross-checked
again via `sha256sum` during that phase's review turn. That value is
reproduced here as a content-integrity pin on the CURRENT, unmodified file
-- not invented from nothing. If the panel is ever intentionally revised,
this constant must be updated deliberately, with the reason documented here
and in docs/PHASE1_BASELINE.md, never silently.

Run with:
    python -m unittest tests.test_phase1_panel_integrity -v
"""

import hashlib
import json
import unittest
from collections import Counter
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
BED_PATH = REPO_ROOT / "resources" / "panel" / "phase1_solid_tumor" / "solid_tumor_phase1_v1.bed"
AUDIT_PATH = REPO_ROOT / "resources" / "panel" / "phase1_solid_tumor" / "phase1_audit.json"

# Documented baseline (phase1_audit.json + docs/PHASE1_BASELINE.md "Phase 1 Scope" table).
EXPECTED_GENE_COUNT = 55
EXPECTED_INTERVAL_COUNT = 901
EXPECTED_TOTAL_BP = 200785
EXPECTED_DUPLICATES = 0

# See module docstring: previously-verified, not invented.
PREVIOUSLY_VERIFIED_SHA256 = "cb0ce560f78334ea1e57974c6bf4d52643d9f7a445bc893a6eea2974d3265b36"


def _load_bed_lines():
    with open(BED_PATH, encoding="utf-8") as fh:
        return [line.rstrip("\n").split("\t") for line in fh if line.strip()]


def _load_audit():
    with open(AUDIT_PATH, encoding="utf-8") as fh:
        return json.load(fh)


class TestPhase1PanelIntegrity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bed_lines = _load_bed_lines()
        cls.audit = _load_audit()

    def test_bed_interval_count(self):
        self.assertEqual(len(self.bed_lines), EXPECTED_INTERVAL_COUNT)

    def test_bed_unique_gene_count(self):
        genes = {row[3] for row in self.bed_lines}
        self.assertEqual(len(genes), EXPECTED_GENE_COUNT)

    def test_bed_total_target_bp(self):
        total_bp = sum(int(row[2]) - int(row[1]) for row in self.bed_lines)
        self.assertEqual(total_bp, EXPECTED_TOTAL_BP)

    def test_no_duplicate_intervals(self):
        counts = Counter((row[0], row[1], row[2]) for row in self.bed_lines)
        duplicates = [k for k, v in counts.items() if v > 1]
        self.assertEqual(duplicates, [])

    def test_bed_intervals_are_well_formed(self):
        # start < end and non-negative, for every interval -- a basic
        # sanity property that should always hold for a valid BED.
        for row in self.bed_lines:
            start, end = int(row[1]), int(row[2])
            self.assertGreaterEqual(start, 0)
            self.assertLess(start, end)

    def test_audit_metadata_matches_documented_baseline(self):
        self.assertEqual(self.audit["genes_verified"], EXPECTED_GENE_COUNT)
        self.assertEqual(self.audit["total_bed_intervals"], EXPECTED_INTERVAL_COUNT)
        self.assertEqual(self.audit["total_unique_target_bp"], EXPECTED_TOTAL_BP)
        self.assertEqual(self.audit["duplicate_intervals_removed"], EXPECTED_DUPLICATES)

    def test_audit_rows_sum_to_documented_totals(self):
        sum_intervals = sum(r["final_interval_count"] for r in self.audit["audit_rows"])
        sum_bp = sum(r["final_target_bp"] for r in self.audit["audit_rows"])
        self.assertEqual(sum_intervals, EXPECTED_INTERVAL_COUNT)
        self.assertEqual(sum_bp, EXPECTED_TOTAL_BP)
        self.assertEqual(len(self.audit["audit_rows"]), EXPECTED_GENE_COUNT)

    def test_bed_genes_match_audit_genes_exactly(self):
        bed_genes = {row[3] for row in self.bed_lines}
        audit_genes = {r["gene"] for r in self.audit["audit_rows"]}
        self.assertEqual(bed_genes, audit_genes)

    def test_per_gene_interval_counts_match_audit(self):
        bed_gene_counts = Counter(row[3] for row in self.bed_lines)
        mismatches = [
            (r["gene"], bed_gene_counts.get(r["gene"], 0), r["final_interval_count"])
            for r in self.audit["audit_rows"]
            if bed_gene_counts.get(r["gene"], 0) != r["final_interval_count"]
        ]
        self.assertEqual(mismatches, [])

    def test_bed_content_hash_matches_previously_verified_value(self):
        actual = hashlib.sha256(BED_PATH.read_bytes()).hexdigest()
        self.assertEqual(
            actual, PREVIOUSLY_VERIFIED_SHA256,
            msg=(
                "Phase 1 panel BED content changed since the last verified "
                "snapshot. If this is an intentional, reviewed panel "
                "update, update PREVIOUSLY_VERIFIED_SHA256 explicitly and "
                "document why here and in docs/PHASE1_BASELINE.md -- do "
                "not update it silently."
            ),
        )


if __name__ == "__main__":
    unittest.main()
