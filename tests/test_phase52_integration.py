"""
Integration tests for GenRichi Phase 5.2 (Input Validation & Failure
Handling) -- MSI/CNV insufficient-data semantics and COSMIC-availability
surfacing, exercised by running the actual Snakemake `script:` files
(calculate_msi.py, calculate_cnv.py, generate_comprehensive_report.py)
against fixture inputs, the same way tests/test_provenance_integration.py
exercises the report script for Phase 5.1.

None of these tests change or re-derive the underlying scientific
calculations/thresholds -- they only verify that a genuine negative result
and an insufficient-data result are never rendered identically.

Run with:
    python -m unittest tests.test_phase52_integration -v
"""

import gzip
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

MSI_SCRIPT = SCRIPTS_DIR / "calculate_msi.py"
CNV_SCRIPT = SCRIPTS_DIR / "calculate_cnv.py"
REPORT_SCRIPT = SCRIPTS_DIR / "generate_comprehensive_report.py"


def _run_script(path, fake_snakemake):
    runpy.run_path(str(path), init_globals={"snakemake": fake_snakemake}, run_name="__main__")


# ── MSI ───────────────────────────────────────────────────────────────────
def _msi_vcf(tmp: Path, name: str, pass_variants: list) -> str:
    """pass_variants: list of (ref, alt) tuples, each written as a PASS record."""
    path = tmp / name
    lines = [
        "##fileformat=VCFv4.2\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    ]
    for i, (ref, alt) in enumerate(pass_variants):
        lines.append(f"chr17\t{7670000+i}\t.\t{ref}\t{alt}\t.\tPASS\t.\n")
    path.write_text("".join(lines), encoding="utf-8")
    return str(path)


def _run_msi(tmp: Path, vcf_path: str, coding_mb=1.5, threshold=10.0):
    score = tmp / "sample.msi"
    dis = tmp / "sample.msi_dis"
    somatic = tmp / "sample.msi_somatic"
    fake_sm = types.SimpleNamespace(
        input=types.SimpleNamespace(vcf=vcf_path),
        output=types.SimpleNamespace(score=str(score), dis=str(dis), somatic=str(somatic)),
        params=types.SimpleNamespace(coding_mb=coding_mb, threshold=threshold),
    )
    _run_script(MSI_SCRIPT, fake_sm)
    return score, dis, somatic


class TestMsiInsufficientData(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)
        self.tmp = Path(self.tmpdir.name)

    def test_zero_pass_variants_is_insufficient_data(self):
        vcf = _msi_vcf(self.tmp, "empty.vcf", [])
        score, dis, _ = _run_msi(self.tmp, vcf)

        score_lines = score.read_text(encoding="utf-8").strip().splitlines()
        self.assertEqual(score_lines[0].split("\t")[-1], "qc_status")
        data_fields = score_lines[1].split("\t")
        self.assertEqual(data_fields[-1], "INSUFFICIENT_DATA")

        dis_text = dis.read_text(encoding="utf-8")
        self.assertIn("qc_status\tINSUFFICIENT_DATA", dis_text)
        self.assertIn("msi_status\tINSUFFICIENT_DATA", dis_text)
        self.assertNotIn("msi_status\tMSS", dis_text)

    def test_nonzero_pass_variants_below_threshold_is_valid_mss(self):
        # A single SNV (not an indel) -> total_variants=1 (nonzero, VALID),
        # ms_indels=0 -> genuinely, validly MSS. Must NOT be flagged
        # INSUFFICIENT_DATA merely because the indel count is zero.
        vcf = _msi_vcf(self.tmp, "one_snv.vcf", [("G", "A")])
        score, dis, _ = _run_msi(self.tmp, vcf)

        data_fields = score.read_text(encoding="utf-8").strip().splitlines()[1].split("\t")
        self.assertEqual(data_fields[-1], "VALID")
        self.assertIn("msi_status\tMSS", dis.read_text(encoding="utf-8"))

    def test_many_indels_above_threshold_is_valid_msi_high(self):
        # 20 homopolymer-context indels over a small coding_mb -> clearly
        # above threshold. Confirms VALID path still classifies MSI-High
        # correctly -- Phase 5.2 must not touch this comparison.
        variants = [("GGG", "G")] * 20
        vcf = _msi_vcf(self.tmp, "many_indels.vcf", variants)
        score, dis, _ = _run_msi(self.tmp, vcf, coding_mb=1.0, threshold=10.0)

        data_fields = score.read_text(encoding="utf-8").strip().splitlines()[1].split("\t")
        self.assertEqual(data_fields[-1], "VALID")
        self.assertIn("msi_status\tMSI-H", dis.read_text(encoding="utf-8"))


# ── CNV ───────────────────────────────────────────────────────────────────
def _regions_file(tmp: Path, name: str, rows: list) -> str:
    """rows: list of (chrom, start, end, region_name, mean_depth)."""
    path = tmp / name
    with gzip.open(path, "wt") as fh:
        for chrom, start, end, region, depth in rows:
            fh.write(f"{chrom}\t{start}\t{end}\t{region}\t{depth}\n")
    return str(path)


def _run_cnv(tmp: Path, tumor_path: str, normal_path: str):
    cnr = tmp / "sample.cnr"
    cns = tmp / "sample.cns"
    call_cns = tmp / "sample.call.cns"
    scatter = tmp / "sample-scatter.png"
    fake_sm = types.SimpleNamespace(
        input=types.SimpleNamespace(tumor_regions=tumor_path, normal_regions=normal_path),
        output=types.SimpleNamespace(cnr=str(cnr), cns=str(cns), call_cns=str(call_cns), scatter=str(scatter)),
        params=types.SimpleNamespace(amp_threshold=0.58, del_threshold=-1.0, min_probes=5),
    )
    _run_script(CNV_SCRIPT, fake_sm)
    return cnr, cns, call_cns, scatter


class TestCnvInsufficientData(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)
        self.tmp = Path(self.tmpdir.name)

    def test_empty_tumor_normal_overlap_is_insufficient_data(self):
        # Tumor and normal region files use non-overlapping coordinates ->
        # the inner join is empty.
        tumor = _regions_file(self.tmp, "tumor.regions.bed.gz", [
            ("chr17", 0, 100, "GENE1", 50.0),
            ("chr17", 100, 200, "GENE2", 60.0),
        ])
        normal = _regions_file(self.tmp, "normal.regions.bed.gz", [
            ("chr17", 500, 600, "GENE1", 45.0),
            ("chr17", 600, 700, "GENE2", 55.0),
        ])
        cnr, cns, call_cns, scatter = _run_cnv(self.tmp, tumor, normal)

        # Must not crash, and must still produce a parseable, correctly
        # columned (if empty) call.cns file with an explicit status marker.
        text = call_cns.read_text(encoding="utf-8")
        self.assertTrue(text.startswith("#qc_status\tINSUFFICIENT_DATA"))

        import pandas as pd
        df = pd.read_csv(call_cns, sep="\t", comment="#")
        self.assertEqual(len(df), 0)
        self.assertIn("log2", df.columns)  # Phase 5.2 empty-columns fix
        self.assertTrue(scatter.exists())  # scatter plot still generated, not crashed

    def test_full_overlap_is_valid(self):
        rows = [("chr17", i * 100, i * 100 + 100, f"GENE{i}", 50.0) for i in range(5)]
        tumor = _regions_file(self.tmp, "tumor2.regions.bed.gz", rows)
        normal = _regions_file(self.tmp, "normal2.regions.bed.gz", rows)
        _, _, call_cns, _ = _run_cnv(self.tmp, tumor, normal)

        text = call_cns.read_text(encoding="utf-8")
        self.assertTrue(text.startswith("#qc_status\tVALID"))

        import pandas as pd
        df = pd.read_csv(call_cns, sep="\t", comment="#")
        self.assertGreater(len(df), 0)


# ── Report-level rendering (MSI/CNV insufficient-data + COSMIC availability) ─
def _report_fixtures(tmp: Path):
    variants = tmp / "variants.tsv"
    variants.write_text(
        "gene\thgvsc\thgvsp\tconsequence\timpact\tdepth\tvaf\tgnomad_af\tclinvar_sig\tcosmic_id\n"
        "TP53\tc.524G>A\tp.Arg175His\tmissense_variant\tHIGH\t83\t0.988\t0\tPathogenic\tCOSV1234\n",
        encoding="utf-8",
    )
    msi = tmp / "sample.msi"
    fastp = tmp / "fastp.json"
    fastp.write_text(json.dumps({"summary": {"before_filtering": {"total_reads": 100, "q30_rate": 0.95}}}), encoding="utf-8")
    flagstat = tmp / "sample.flagstat"
    flagstat.write_text("95 + 0 mapped (95.00%:N/A)\n", encoding="utf-8")
    mosdepth = tmp / "sample.mosdepth.summary.txt"
    mosdepth.write_text(
        "chrom\tlength\tbases\tmean\tmin\tmax\n"
        "chr17_region\t1000\t50000\t50.0\t1\t100\n"
        "total_region\t1000\t50000\t50.0\t1\t100\n",
        encoding="utf-8",
    )
    markdup = tmp / "sample.markdup_metrics.txt"
    markdup.write_text("LIBRARY\tPERCENT_DUPLICATION\tOTHER\nlib1\t0.05\t1\n", encoding="utf-8")
    panel_bed = tmp / "panel.bed"
    panel_bed.write_text("chr17\t0\t1\n", encoding="utf-8")
    tmb_bed = tmp / "tmb_cds.bed"
    tmb_bed.write_text("chr17\t100\t200\tTP53\n", encoding="utf-8")
    return dict(variants=str(variants), msi=str(msi), fastp=str(fastp), flagstat=str(flagstat),
                mosdepth=str(mosdepth), markdup=str(markdup), panel_bed=str(panel_bed), tmb_bed=str(tmb_bed))


def _run_report(tmp: Path, fixtures, msi_content, call_cns_content, cosmic_path):
    msi_path = Path(fixtures["msi"])
    msi_path.write_text(msi_content, encoding="utf-8")
    call_cns = tmp / "sample.call.cns"
    call_cns.write_text(call_cns_content, encoding="utf-8")

    html_out = tmp / "report.html"
    provenance_out = tmp / "provenance.json"

    fake_input = types.SimpleNamespace(
        variants=fixtures["variants"], cnr=str(call_cns), call_cns=str(call_cns),
        cnv_scatter=str(tmp / "no_scatter.png"), msi_score=str(msi_path),
        tumor_fastp=fixtures["fastp"], normal_fastp=fixtures["fastp"],
        tumor_flagstat=fixtures["flagstat"], normal_flagstat=fixtures["flagstat"],
        tumor_mosdepth=fixtures["mosdepth"], normal_mosdepth=fixtures["mosdepth"],
        tumor_markdup=fixtures["markdup"], normal_markdup=fixtures["markdup"],
    )
    fake_output = types.SimpleNamespace(html=str(html_out), provenance=str(provenance_out))
    fake_params = types.SimpleNamespace(
        sample_id="TESTSAMPLE", run_id=None, patient_id="P1", sex="Female", tumor_type="Breast_Cancer",
        panel_name="Test Panel", company="GenRichi", logo=str(tmp / "no_logo.png"), show_synonymous=False,
        msi_threshold=10.0, tmb_coding_bed=fixtures["tmb_bed"], tmb_partial_threshold=0.80,
        tmb_high_threshold=10.0, panel_bed=fixtures["panel_bed"], cnv_amp_threshold=0.58, cnv_del_threshold=-1.0,
    )
    fake_config = {
        "ref": {"genome": "/data/ref/hg38.fa", "dbsnp": "/data/dbsnp/dbsnp_146.hg38.vcf.gz",
                "gnomad": "/data/gnomad/af-only-gnomad.hg38.vcf.gz", "pon": None},
        "panel": {"name": "Test Panel", "bed": fixtures["panel_bed"]},
        "annotation": {
            "vep": {"genome_build": "GRCh38", "extra": "--cache_version 113"},
            "cosmic": {"vcf": cosmic_path},
            "clinvar": {"vcf": "/data/clinvar/clinvar_20240101.vcf.gz"},
        },
    }
    fake_sm = types.SimpleNamespace(input=fake_input, output=fake_output, params=fake_params, config=fake_config)
    _run_script(REPORT_SCRIPT, fake_sm)
    return html_out, provenance_out


class TestReportInsufficientDataRendering(unittest.TestCase):
    def setUp(self):
        self.tmpdir_ctx = tempfile.TemporaryDirectory()
        self.tmp = Path(self.tmpdir_ctx.name)
        self.addCleanup(self.tmpdir_ctx.cleanup)
        self.fixtures = _report_fixtures(self.tmp)

    def test_msi_insufficient_data_never_shown_as_mss(self):
        msi_content = "Total_Number_of_Sites\tNumber_of_Somatic_Sites\t%\tqc_status\n1.5\t0\t0.0\tINSUFFICIENT_DATA\n"
        call_cns = "chromosome\tstart\tend\tgene\tlog2\tcn\n"
        html_out, _ = _run_report(self.tmp, self.fixtures, msi_content, call_cns, cosmic_path=None)
        html = html_out.read_text(encoding="utf-8")
        self.assertIn("INSUFFICIENT_DATA", html)
        self.assertNotIn(">MSS<", html)
        self.assertNotIn("0.0%", html)  # must not show a percentage for insufficient data

    def test_msi_valid_mss_rendered_normally(self):
        msi_content = "Total_Number_of_Sites\tNumber_of_Somatic_Sites\t%\tqc_status\n1.5\t0\t0.0\tVALID\n"
        call_cns = "chromosome\tstart\tend\tgene\tlog2\tcn\n"
        html_out, _ = _run_report(self.tmp, self.fixtures, msi_content, call_cns, cosmic_path=None)
        html = html_out.read_text(encoding="utf-8")
        self.assertIn(">MSS<", html)
        self.assertNotIn("INSUFFICIENT_DATA", html)

    def test_legacy_3_column_msi_file_treated_as_valid(self):
        # Backward compatibility: files written before Phase 5.2 (no 4th column).
        msi_content = "Total_Number_of_Sites\tSomatic_Sites\tPercent\n1.5\t0\t0.0\n"
        call_cns = "chromosome\tstart\tend\tgene\tlog2\tcn\n"
        html_out, _ = _run_report(self.tmp, self.fixtures, msi_content, call_cns, cosmic_path=None)
        html = html_out.read_text(encoding="utf-8")
        self.assertIn(">MSS<", html)
        self.assertNotIn("INSUFFICIENT_DATA", html)

    def test_cnv_insufficient_data_never_shown_as_no_calls(self):
        msi_content = "Total_Number_of_Sites\tNumber_of_Somatic_Sites\t%\tqc_status\n1.5\t0\t0.0\tVALID\n"
        call_cns = "#qc_status\tINSUFFICIENT_DATA\nchromosome\tstart\tend\tgene\tlog2\tprobes\tweight\tcn\ttype\n"
        html_out, _ = _run_report(self.tmp, self.fixtures, msi_content, call_cns, cosmic_path=None)
        html = html_out.read_text(encoding="utf-8")
        self.assertIn("Insufficient tumor/normal coverage overlap", html)
        self.assertNotIn("No significant CNV calls", html)

    def test_cnv_valid_zero_calls_rendered_as_genuine_negative(self):
        msi_content = "Total_Number_of_Sites\tNumber_of_Somatic_Sites\t%\tqc_status\n1.5\t0\t0.0\tVALID\n"
        call_cns = "#qc_status\tVALID\nchromosome\tstart\tend\tgene\tlog2\tprobes\tweight\tcn\ttype\n"
        html_out, _ = _run_report(self.tmp, self.fixtures, msi_content, call_cns, cosmic_path=None)
        html = html_out.read_text(encoding="utf-8")
        self.assertIn("No significant CNV calls", html)
        self.assertNotIn("Insufficient tumor/normal coverage overlap", html)

    def test_legacy_call_cns_without_marker_treated_as_valid(self):
        msi_content = "Total_Number_of_Sites\tNumber_of_Somatic_Sites\t%\tqc_status\n1.5\t0\t0.0\tVALID\n"
        call_cns = "chromosome\tstart\tend\tgene\tlog2\tcn\n"  # no leading marker line
        html_out, _ = _run_report(self.tmp, self.fixtures, msi_content, call_cns, cosmic_path=None)
        html = html_out.read_text(encoding="utf-8")
        self.assertIn("No significant CNV calls", html)

    def test_cosmic_unavailable_notice_shown(self):
        msi_content = "Total_Number_of_Sites\tNumber_of_Somatic_Sites\t%\tqc_status\n1.5\t0\t0.0\tVALID\n"
        call_cns = "chromosome\tstart\tend\tgene\tlog2\tcn\n"
        html_out, provenance_out = _run_report(
            self.tmp, self.fixtures, msi_content, call_cns,
            cosmic_path=str(self.tmp / "does_not_exist_cosmic.vcf.gz"),
        )
        html = html_out.read_text(encoding="utf-8")
        self.assertIn("COSMIC annotation was not performed", html)
        provenance = json.loads(provenance_out.read_text(encoding="utf-8"))
        self.assertFalse(provenance["annotation_databases"]["cosmic_available"])

    def test_cosmic_available_no_notice_shown(self):
        cosmic_file = self.tmp / "cosmic.vcf.gz"
        cosmic_file.write_bytes(b"stub")
        msi_content = "Total_Number_of_Sites\tNumber_of_Somatic_Sites\t%\tqc_status\n1.5\t0\t0.0\tVALID\n"
        call_cns = "chromosome\tstart\tend\tgene\tlog2\tcn\n"
        html_out, provenance_out = _run_report(
            self.tmp, self.fixtures, msi_content, call_cns, cosmic_path=str(cosmic_file),
        )
        html = html_out.read_text(encoding="utf-8")
        self.assertNotIn("COSMIC annotation was not performed", html)
        provenance = json.loads(provenance_out.read_text(encoding="utf-8"))
        self.assertTrue(provenance["annotation_databases"]["cosmic_available"])


if __name__ == "__main__":
    unittest.main()
