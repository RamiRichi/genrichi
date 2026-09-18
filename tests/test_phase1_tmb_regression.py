"""
GenRichi Phase 5.3 -- Phase 1 scientific/regression tests.

Focused regression tests locking in the EXISTING TMB behavior of
workflow/scripts/generate_comprehensive_report.py: the coding-Mb
denominator, partial-run detection, the full-run calculation path, and the
existing TMB-High/Low threshold comparison.

This file does NOT modify generate_comprehensive_report.py, does NOT change
the TMB methodology/thresholds/panel/denominator definition, and does NOT
invent new expected values -- every assertion below is either an exact,
hand-computed arithmetic result from a controlled fixture (so the numbers
are independently verifiable), or a literal reproduction of the report
script's own existing, unmodified string templates. The purpose is purely
to lock the current behavior against future silent regression.

All tests run the REAL script via runpy (same established pattern as
tests/test_provenance_integration.py and tests/test_phase52_integration.py)
against small, synthetic, deterministic fixtures -- no external dataset or
pipeline execution required.

Run with:
    python -m unittest tests.test_phase1_tmb_regression -v
"""

import json
import re
import runpy
import sys
import tempfile
import types
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = REPO_ROOT / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

REPORT_SCRIPT = SCRIPTS_DIR / "generate_comprehensive_report.py"

TMB_HIGH_THRESHOLD = 10.0
PARTIAL_THRESHOLD = 0.80


def _variants_tsv(tmp: Path, name: str, n_coding: int, n_noncoding: int) -> str:
    """
    n_coding rows use "missense_variant" (a member of the report script's
    own, unmodified `coding_consequences` set); n_noncoding rows use
    "intron_variant" (not a member), so the resulting TMB numerator is
    exactly n_coding, independent of n_noncoding.
    """
    rows = ["gene\thgvsc\thgvsp\tconsequence\timpact\tdepth\tvaf\tgnomad_af\tclinvar_sig\tcosmic_id\n"]
    for i in range(n_coding):
        rows.append(f"GENE{i}\tc.1A>T\tp.Lys1Ile\tmissense_variant\tMODERATE\t50\t0.3\t0\t\t\n")
    for i in range(n_noncoding):
        rows.append(f"NGENE{i}\tc.1A>T\tp.=\tintron_variant\tMODIFIER\t50\t0.3\t0\t\t\n")
    path = tmp / name
    path.write_text("".join(rows), encoding="utf-8")
    return str(path)


def _tmb_coding_bed_one_mb(tmp: Path) -> str:
    """Exactly 1,000,000 bp in a single interval -> coding-Mb denominator == 1.0 exactly."""
    path = tmp / "tmb_one_mb.bed"
    path.write_text("chr17\t0\t1000000\tTP53\n", encoding="utf-8")
    return str(path)


def _panel_bed_two_chroms(tmp: Path) -> str:
    path = tmp / "panel_two_chroms.bed"
    path.write_text("chr17\t0\t1\nchr18\t0\t1\n", encoding="utf-8")
    return str(path)


def _mosdepth_full_coverage(tmp: Path) -> str:
    """Both panel chromosomes covered at >=1x -> not a partial run (2/2 = 100% >= 80%)."""
    path = tmp / "mosdepth_full.txt"
    path.write_text(
        "chrom\tlength\tbases\tmean\tmin\tmax\n"
        "chr17_region\t1000\t50000\t50.0\t1\t100\n"
        "chr18_region\t1000\t50000\t50.0\t1\t100\n"
        "total_region\t2000\t100000\t50.0\t1\t100\n",
        encoding="utf-8",
    )
    return str(path)


def _mosdepth_partial_coverage(tmp: Path) -> str:
    """Only chr17 covered of 2 panel chromosomes -> 1/2 = 50% < 80% -> partial run."""
    path = tmp / "mosdepth_partial.txt"
    path.write_text(
        "chrom\tlength\tbases\tmean\tmin\tmax\n"
        "chr17_region\t1000\t50000\t50.0\t1\t100\n"
        "total_region\t1000\t50000\t50.0\t1\t100\n",
        encoding="utf-8",
    )
    return str(path)


def _common_fixtures(tmp: Path):
    msi = tmp / "sample.msi"
    msi.write_text("Total_Number_of_Sites\tNumber_of_Somatic_Sites\t%\tqc_status\n1.5\t0\t0.0\tVALID\n", encoding="utf-8")
    fastp = tmp / "fastp.json"
    fastp.write_text(json.dumps({"summary": {"before_filtering": {"total_reads": 100, "q30_rate": 0.95}}}), encoding="utf-8")
    flagstat = tmp / "sample.flagstat"
    flagstat.write_text("95 + 0 mapped (95.00%:N/A)\n", encoding="utf-8")
    markdup = tmp / "sample.markdup_metrics.txt"
    markdup.write_text("LIBRARY\tPERCENT_DUPLICATION\tOTHER\nlib1\t0.05\t1\n", encoding="utf-8")
    call_cns = tmp / "sample.call.cns"
    call_cns.write_text("#qc_status\tVALID\nchromosome\tstart\tend\tgene\tlog2\tcn\n", encoding="utf-8")
    return dict(msi=str(msi), fastp=str(fastp), flagstat=str(flagstat), markdup=str(markdup), call_cns=str(call_cns))


def _run_report(tmp: Path, *, n_coding, n_noncoding, mosdepth_path, panel_bed_path, tmb_bed_path):
    common = _common_fixtures(tmp)
    variants = _variants_tsv(tmp, "variants.tsv", n_coding, n_noncoding)

    html_out = tmp / "report.html"
    provenance_out = tmp / "provenance.json"

    fake_input = types.SimpleNamespace(
        variants=variants, cnr=common["call_cns"], call_cns=common["call_cns"],
        cnv_scatter=str(tmp / "no_scatter.png"), msi_score=common["msi"],
        tumor_fastp=common["fastp"], normal_fastp=common["fastp"],
        tumor_flagstat=common["flagstat"], normal_flagstat=common["flagstat"],
        tumor_mosdepth=mosdepth_path, normal_mosdepth=mosdepth_path,
        tumor_markdup=common["markdup"], normal_markdup=common["markdup"],
    )
    fake_output = types.SimpleNamespace(html=str(html_out), provenance=str(provenance_out))
    fake_params = types.SimpleNamespace(
        sample_id="TMBTEST", run_id=None, patient_id="P1", sex="Female", tumor_type="Test",
        panel_name="Test Panel", company="GenRichi", logo=str(tmp / "no_logo.png"), show_synonymous=False,
        msi_threshold=10.0, tmb_coding_bed=tmb_bed_path, tmb_partial_threshold=PARTIAL_THRESHOLD,
        tmb_high_threshold=TMB_HIGH_THRESHOLD, panel_bed=panel_bed_path,
        cnv_amp_threshold=0.58, cnv_del_threshold=-1.0,
    )
    fake_config = {
        "ref": {"genome": "/data/ref/hg38.fa", "dbsnp": "/data/dbsnp/dbsnp_146.hg38.vcf.gz",
                "gnomad": "/data/gnomad/af-only-gnomad.hg38.vcf.gz", "pon": None},
        "panel": {"name": "Test Panel", "bed": panel_bed_path},
        "annotation": {
            "vep": {"genome_build": "GRCh38", "extra": "--cache_version 113"},
            "cosmic": {"vcf": None},
            "clinvar": {"vcf": "/data/clinvar/clinvar_20240101.vcf.gz"},
        },
    }
    fake_sm = types.SimpleNamespace(input=fake_input, output=fake_output, params=fake_params, config=fake_config)
    runpy.run_path(str(REPORT_SCRIPT), init_globals={"snakemake": fake_sm}, run_name="__main__")
    return html_out.read_text(encoding="utf-8")


class TestTmbCodingDenominator(unittest.TestCase):
    """Locks the coding-Mb denominator computed from tmb.coding_bed (_read_cds_bed)."""

    def setUp(self):
        self.tmpdir_ctx = tempfile.TemporaryDirectory()
        self.tmp = Path(self.tmpdir_ctx.name)
        self.addCleanup(self.tmpdir_ctx.cleanup)

    def test_one_megabase_bed_yields_exact_1_0000_mb_denominator(self):
        html = _run_report(
            self.tmp, n_coding=5, n_noncoding=0,
            mosdepth_path=_mosdepth_full_coverage(self.tmp),
            panel_bed_path=_panel_bed_two_chroms(self.tmp),
            tmb_bed_path=_tmb_coding_bed_one_mb(self.tmp),
        )
        self.assertIn("TMB denominator: 1.0000&nbsp;Mb", html)
        self.assertIn("(1,000,000&nbsp;bp, 1&nbsp;intervals, 1&nbsp;genes)", html)


class TestTmbPartialRunDetection(unittest.TestCase):
    """Locks the existing partial-run -> TMB N/A behavior (matches the real chr17-only demo)."""

    def setUp(self):
        self.tmpdir_ctx = tempfile.TemporaryDirectory()
        self.tmp = Path(self.tmpdir_ctx.name)
        self.addCleanup(self.tmpdir_ctx.cleanup)

    def test_partial_coverage_reports_tmb_na(self):
        html = _run_report(
            self.tmp, n_coding=20, n_noncoding=0,  # would be well above threshold if NOT partial
            mosdepth_path=_mosdepth_partial_coverage(self.tmp),
            panel_bed_path=_panel_bed_two_chroms(self.tmp),
            tmb_bed_path=_tmb_coding_bed_one_mb(self.tmp),
        )
        self.assertIn("TMB: N/A (partial-region run", html)
        self.assertIn("Partial-Region Run", html)
        self.assertIn("Denominator would be: 1.0000&nbsp;Mb", html)
        self.assertNotIn("TMB-High", html)
        self.assertNotIn("TMB-Low", html)

    def test_full_coverage_of_both_panel_chroms_is_not_partial(self):
        html = _run_report(
            self.tmp, n_coding=1, n_noncoding=0,
            mosdepth_path=_mosdepth_full_coverage(self.tmp),
            panel_bed_path=_panel_bed_two_chroms(self.tmp),
            tmb_bed_path=_tmb_coding_bed_one_mb(self.tmp),
        )
        self.assertNotIn("partial-region run", html)
        self.assertIn("TMB: 1.0 mut/Mb [TMB-Low]", html)


class TestTmbFullRunCalculationAndThreshold(unittest.TestCase):
    """Locks the full-run TMB calculation path and the existing High/Low threshold comparison."""

    def setUp(self):
        self.tmpdir_ctx = tempfile.TemporaryDirectory()
        self.tmp = Path(self.tmpdir_ctx.name)
        self.addCleanup(self.tmpdir_ctx.cleanup)

    def _run(self, n_coding, n_noncoding=3):
        return _run_report(
            self.tmp, n_coding=n_coding, n_noncoding=n_noncoding,
            mosdepth_path=_mosdepth_full_coverage(self.tmp),
            panel_bed_path=_panel_bed_two_chroms(self.tmp),
            tmb_bed_path=_tmb_coding_bed_one_mb(self.tmp),
        )

    def test_twelve_coding_variants_over_one_mb_is_tmb_high(self):
        # 12 coding variants / 1.0 Mb = 12.0 mut/Mb >= 10.0 threshold
        html = self._run(n_coding=12)
        self.assertIn("TMB: 12.0 mut/Mb [TMB-High]", html)

    def test_five_coding_variants_over_one_mb_is_tmb_low(self):
        # 5 coding variants / 1.0 Mb = 5.0 mut/Mb < 10.0 threshold
        html = self._run(n_coding=5)
        self.assertIn("TMB: 5.0 mut/Mb [TMB-Low]", html)

    def test_exactly_at_threshold_is_tmb_high(self):
        # 10 / 1.0 Mb = 10.0 == threshold; existing comparison is >=
        html = self._run(n_coding=10)
        self.assertIn("TMB: 10.0 mut/Mb [TMB-High]", html)

    def test_noncoding_variants_do_not_affect_tmb_numerator(self):
        # Same coding count (5), very different non-coding count -> identical TMB value.
        html_a = self._run(n_coding=5, n_noncoding=0)
        html_b = self._run(n_coding=5, n_noncoding=50)
        self.assertIn("TMB: 5.0 mut/Mb [TMB-Low]", html_a)
        self.assertIn("TMB: 5.0 mut/Mb [TMB-Low]", html_b)


if __name__ == "__main__":
    unittest.main()
