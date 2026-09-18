"""
GenRichi Phase 5.3 -- Phase 1 scientific/regression tests.

Tests for workflow/scripts/sample_columns.py and its integration into
workflow/scripts/vcf_to_somatic_table.py: hardening the tumor/normal VCF
sample-column identification that was previously a positional assumption
("tumor is always the last column") with no validation.

Real-world naming, confirmed against the live HCC1395 pipeline run's actual
VCF header during this phase's implementation:

    #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT \
        HCC1395_demo_normal  HCC1395_demo_tumor

i.e. normal first, tumor last -- matching workflow/rules/align_paired.smk's
"{sample}_tumor"/"{sample}_normal" SM tags and
workflow/rules/somatic_calling_paired.smk's --normal-sample flag.

Does NOT change Mutect2 parameters or scientific filtering logic -- this is
purely a data-extraction correctness hardening.

Run with:
    python -m unittest tests.test_phase1_sample_columns -v
"""

import runpy
import sys
import tempfile
import types
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = REPO_ROOT / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import sample_columns  # noqa: E402

VCF_TO_TABLE_SCRIPT = SCRIPTS_DIR / "vcf_to_somatic_table.py"

FIXED_VCF_COLUMNS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]


class TestResolveTumorNormalColumns(unittest.TestCase):
    """Unit tests for the pure resolver function."""

    # (a) expected normal -> tumor ordering (real, confirmed Mutect2 layout)
    def test_expected_normal_then_tumor_ordering(self):
        col_names = FIXED_VCF_COLUMNS + ["HCC1395_demo_normal", "HCC1395_demo_tumor"]
        tumor_col, normal_col = sample_columns.resolve_tumor_normal_columns(
            col_names, expected_sample_id="HCC1395_demo"
        )
        self.assertEqual(tumor_col, "HCC1395_demo_tumor")
        self.assertEqual(normal_col, "HCC1395_demo_normal")

    def test_resolution_without_expected_sample_id(self):
        col_names = FIXED_VCF_COLUMNS + ["S1_normal", "S1_tumor"]
        tumor_col, normal_col = sample_columns.resolve_tumor_normal_columns(col_names)
        self.assertEqual(tumor_col, "S1_tumor")
        self.assertEqual(normal_col, "S1_normal")

    # (b) unexpected / reversed ordering -- must still resolve correctly by
    # NAME, proving the fix is robust to ordering (not just safe-failing)
    def test_reversed_tumor_then_normal_ordering_still_resolves_correctly(self):
        col_names = FIXED_VCF_COLUMNS + ["HCC1395_demo_tumor", "HCC1395_demo_normal"]
        tumor_col, normal_col = sample_columns.resolve_tumor_normal_columns(
            col_names, expected_sample_id="HCC1395_demo"
        )
        self.assertEqual(tumor_col, "HCC1395_demo_tumor")
        self.assertEqual(normal_col, "HCC1395_demo_normal")

    # (c) malformed / missing sample columns
    def test_only_one_sample_column_raises(self):
        col_names = FIXED_VCF_COLUMNS + ["HCC1395_demo_tumor"]
        with self.assertRaises(sample_columns.SampleColumnError):
            sample_columns.resolve_tumor_normal_columns(col_names)

    def test_three_sample_columns_raises(self):
        col_names = FIXED_VCF_COLUMNS + ["A_tumor", "A_normal", "B_tumor"]
        with self.assertRaises(sample_columns.SampleColumnError):
            sample_columns.resolve_tumor_normal_columns(col_names)

    def test_zero_sample_columns_raises(self):
        with self.assertRaises(sample_columns.SampleColumnError):
            sample_columns.resolve_tumor_normal_columns(list(FIXED_VCF_COLUMNS))

    def test_unnamed_sample_columns_raise(self):
        # e.g. an externally-supplied VCF that doesn't follow GenRichi's
        # "{sample}_tumor"/"{sample}_normal" SM-tag convention at all.
        col_names = FIXED_VCF_COLUMNS + ["SAMPLE_A", "SAMPLE_B"]
        with self.assertRaises(sample_columns.SampleColumnError):
            sample_columns.resolve_tumor_normal_columns(col_names)

    def test_both_columns_named_tumor_is_ambiguous_and_raises(self):
        col_names = FIXED_VCF_COLUMNS + ["A_tumor", "B_tumor"]
        with self.assertRaises(sample_columns.SampleColumnError):
            sample_columns.resolve_tumor_normal_columns(col_names)

    def test_both_columns_named_normal_is_ambiguous_and_raises(self):
        col_names = FIXED_VCF_COLUMNS + ["A_normal", "B_normal"]
        with self.assertRaises(sample_columns.SampleColumnError):
            sample_columns.resolve_tumor_normal_columns(col_names)

    def test_expected_sample_id_mismatch_raises(self):
        # Correctly-named columns, but for a DIFFERENT sample than the rule
        # instance expects -- e.g. a wrong/stale VCF fed into the wrong
        # sample's rule invocation.
        col_names = FIXED_VCF_COLUMNS + ["OTHER_SAMPLE_normal", "OTHER_SAMPLE_tumor"]
        with self.assertRaises(sample_columns.SampleColumnError):
            sample_columns.resolve_tumor_normal_columns(
                col_names, expected_sample_id="HCC1395_demo"
            )

    def test_expected_sample_id_match_does_not_raise(self):
        col_names = FIXED_VCF_COLUMNS + ["HCC1395_demo_normal", "HCC1395_demo_tumor"]
        # should not raise
        sample_columns.resolve_tumor_normal_columns(
            col_names, expected_sample_id="HCC1395_demo"
        )

    def test_error_messages_are_actionable_not_a_traceback(self):
        col_names = FIXED_VCF_COLUMNS + ["SAMPLE_A", "SAMPLE_B"]
        with self.assertRaises(sample_columns.SampleColumnError) as ctx:
            sample_columns.resolve_tumor_normal_columns(col_names)
        message = str(ctx.exception)
        self.assertNotIn("Traceback", message)
        self.assertIn("SAMPLE_A", message)
        self.assertIn("SAMPLE_B", message)


def _fmt_for(col_name: str) -> str:
    """
    FORMAT values keyed by column name (not position): a "*_tumor" column
    gets AD=1,82 -> total depth 83, alt_reads 82 (matching the real HCC1395
    baseline exactly); a "*_normal" column gets clearly DIFFERENT values
    (AD=39,1 -> depth 40, alt_reads 1), so a test can confirm the
    TUMOR-named column's specific values were extracted -- not merely "a"
    column's values, which would pass even with a silent swap.
    """
    if col_name.endswith("_tumor"):
        return "0/1:1,82:83:.:.:."
    return "0/1:39,1:40:.:.:."


def _minimal_vcf(tmp: Path, name: str, sample_cols) -> str:
    """A minimal, syntactically valid VCF with one PASS record."""
    path = tmp / name
    header = "\t".join(FIXED_VCF_COLUMNS + list(sample_cols))
    fmt_cols = "\t".join(_fmt_for(c) for c in sample_cols)
    lines = [
        "##fileformat=VCFv4.2\n",
        f"#{header}\n",
        f"chr17\t7675088\t.\tC\tT\t.\tPASS\t.\tGT:AD:DP:GQ:PL:AF\t{fmt_cols}\n",
    ]
    path.write_text("".join(lines), encoding="utf-8")
    return str(path)


def _load_single_row(tsv_path: Path):
    import csv
    with open(tsv_path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 1, f"expected exactly 1 output row, got {len(rows)}"
    return rows[0]


class TestVcfToSomaticTableIntegration(unittest.TestCase):
    """
    Integration tests running the REAL vcf_to_somatic_table.py script (via
    runpy, same pattern as the Phase 5.1/5.2 integration tests) to confirm
    the hardening works end-to-end, not just in the isolated helper.
    """

    def setUp(self):
        self.tmpdir_ctx = tempfile.TemporaryDirectory()
        self.tmp = Path(self.tmpdir_ctx.name)
        self.addCleanup(self.tmpdir_ctx.cleanup)

    def _run(self, vcf_path, sample_id="HCC1395_demo"):
        out_path = self.tmp / "out.tsv"
        fake_sm = types.SimpleNamespace(
            input=types.SimpleNamespace(vcf=vcf_path),
            output=types.SimpleNamespace(tsv=str(out_path)),
            wildcards=types.SimpleNamespace(sample=sample_id),
        )
        runpy.run_path(
            str(VCF_TO_TABLE_SCRIPT), init_globals={"snakemake": fake_sm}, run_name="__main__"
        )
        return out_path

    def test_correct_ordering_produces_expected_tumor_values(self):
        vcf = _minimal_vcf(self.tmp, "correct.vcf", ["HCC1395_demo_normal", "HCC1395_demo_tumor"])
        row = _load_single_row(self._run(vcf))
        self.assertEqual(row["depth"], "83")
        self.assertEqual(row["alt_reads"], "82")

    def test_reversed_ordering_still_produces_correct_tumor_values(self):
        # Proves the hardening is robust to ordering, not just safe-failing:
        # even with tumor listed FIRST, the tumor-named column's own values
        # (not the normal-named column's distinct values) must be extracted.
        vcf = _minimal_vcf(self.tmp, "reversed.vcf", ["HCC1395_demo_tumor", "HCC1395_demo_normal"])
        row = _load_single_row(self._run(vcf))
        self.assertEqual(row["depth"], "83")
        self.assertEqual(row["alt_reads"], "82")

    def test_missing_sample_column_fails_loudly_not_silently(self):
        vcf = _minimal_vcf(self.tmp, "malformed.vcf", ["HCC1395_demo_tumor"])
        with self.assertRaises(sample_columns.SampleColumnError):
            self._run(vcf)

    def test_wrong_sample_id_fails_loudly(self):
        vcf = _minimal_vcf(self.tmp, "wrong_sample.vcf", ["OTHER_normal", "OTHER_tumor"])
        with self.assertRaises(sample_columns.SampleColumnError):
            self._run(vcf, sample_id="HCC1395_demo")

    def test_missing_wildcards_degrades_gracefully_when_names_are_valid(self):
        # If snakemake.wildcards isn't present at all (defensive), the
        # script must still work as long as sample columns are unambiguous.
        vcf = _minimal_vcf(self.tmp, "no_wildcards.vcf", ["S1_normal", "S1_tumor"])
        out_path_target = self.tmp / "out2.tsv"
        fake_sm = types.SimpleNamespace(
            input=types.SimpleNamespace(vcf=vcf),
            output=types.SimpleNamespace(tsv=str(out_path_target)),
        )
        runpy.run_path(
            str(VCF_TO_TABLE_SCRIPT), init_globals={"snakemake": fake_sm}, run_name="__main__"
        )
        row = _load_single_row(out_path_target)
        self.assertEqual(row["depth"], "83")
        self.assertEqual(row["alt_reads"], "82")


if __name__ == "__main__":
    unittest.main()
