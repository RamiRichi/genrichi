"""
Integration/end-to-end test for GenRichi Phase 5.1 (Reproducibility & Provenance).

Runs the real `workflow/scripts/generate_comprehensive_report.py` — the
script behind the HCC1395 Phase 1 baseline's `comprehensive_report` rule —
against minimal fixture inputs (standing in for a Snakemake `script:`
invocation) and confirms:

  * the report still renders normally (no regression in existing behavior)
  * a deterministic, machine-readable provenance block is embedded in the
    HTML output
  * a provenance JSON sidecar file is written and is consistent with the
    embedded block

Run with:
    python -m unittest tests.test_provenance_integration -v
"""

import json
import re
import sys
import tempfile
import types
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = REPO_ROOT / "workflow" / "scripts"
REPORT_SCRIPT = SCRIPTS_DIR / "generate_comprehensive_report.py"

sys.path.insert(0, str(SCRIPTS_DIR))


def _write_fixtures(tmpdir: Path):
    """Write minimal, valid input fixtures for generate_comprehensive_report.py."""
    variants = tmpdir / "variants.tsv"
    variants.write_text(
        "gene\thgvsc\thgvsp\tconsequence\timpact\tdepth\tvaf\tgnomad_af\tclinvar_sig\tcosmic_id\n"
        "TP53\tc.524G>A\tp.Arg175His\tmissense_variant\tHIGH\t83\t0.988\t0\tPathogenic\tCOSV1234\n",
        encoding="utf-8",
    )

    call_cns = tmpdir / "call.cns"
    call_cns.write_text(
        "chromosome\tstart\tend\tgene\tlog2\tcn\n"
        "chr17\t100\t200\tTP53\t0.1\t2\n",
        encoding="utf-8",
    )

    msi = tmpdir / "sample.msi"
    msi.write_text("Total_Number_of_Sites\tSomatic_Sites\tPercent\n100\t0\t0.0\n", encoding="utf-8")

    fastp = tmpdir / "fastp.json"
    fastp.write_text(
        json.dumps({"summary": {"before_filtering": {"total_reads": 100, "q30_rate": 0.95}}}),
        encoding="utf-8",
    )

    flagstat = tmpdir / "sample.flagstat"
    flagstat.write_text("95 + 0 mapped (95.00%:N/A)\n", encoding="utf-8")

    mosdepth = tmpdir / "sample.mosdepth.summary.txt"
    mosdepth.write_text(
        "chrom\tlength\tbases\tmean\tmin\tmax\n"
        "chr17_region\t1000\t50000\t50.0\t1\t100\n"
        "total_region\t1000\t50000\t50.0\t1\t100\n",
        encoding="utf-8",
    )

    markdup = tmpdir / "sample.markdup_metrics.txt"
    markdup.write_text(
        "LIBRARY\tPERCENT_DUPLICATION\tOTHER\nlib1\t0.05\t1\n", encoding="utf-8"
    )

    panel_bed = tmpdir / "panel.bed"
    panel_bed.write_text("chr17\t0\t1\n", encoding="utf-8")

    tmb_bed = tmpdir / "tmb_cds.bed"
    tmb_bed.write_text("chr17\t100\t200\tTP53\n", encoding="utf-8")

    return {
        "variants": str(variants),
        "call_cns": str(call_cns),
        "msi": str(msi),
        "fastp": str(fastp),
        "flagstat": str(flagstat),
        "mosdepth": str(mosdepth),
        "markdup": str(markdup),
        "panel_bed": str(panel_bed),
        "tmb_bed": str(tmb_bed),
    }


class TestComprehensiveReportProvenanceIntegration(unittest.TestCase):
    def setUp(self):
        self.tmpdir_ctx = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self.tmpdir_ctx.name)
        self.addCleanup(self.tmpdir_ctx.cleanup)
        self.fixtures = _write_fixtures(self.tmpdir)

        self.html_out = self.tmpdir / "report.html"
        self.provenance_out = self.tmpdir / "provenance.json"

        fake_input = types.SimpleNamespace(
            variants=self.fixtures["variants"],
            cnr=self.fixtures["call_cns"],  # unused by the script; any real path is fine
            call_cns=self.fixtures["call_cns"],
            cnv_scatter=str(self.tmpdir / "does_not_exist.png"),  # handled gracefully
            msi_score=self.fixtures["msi"],
            tumor_fastp=self.fixtures["fastp"],
            normal_fastp=self.fixtures["fastp"],
            tumor_flagstat=self.fixtures["flagstat"],
            normal_flagstat=self.fixtures["flagstat"],
            tumor_mosdepth=self.fixtures["mosdepth"],
            normal_mosdepth=self.fixtures["mosdepth"],
            tumor_markdup=self.fixtures["markdup"],
            normal_markdup=self.fixtures["markdup"],
        )
        fake_output = types.SimpleNamespace(
            html=str(self.html_out), provenance=str(self.provenance_out)
        )
        fake_params = types.SimpleNamespace(
            sample_id="HCC1395_demo",
            run_id=None,
            patient_id="HCC1395",
            sex="Female",
            tumor_type="Breast_Cancer",
            panel_name="GenRichi Comprehensive Cancer Panel v1",
            company="GenRichi",
            logo=str(self.tmpdir / "no_logo.png"),
            show_synonymous=False,
            msi_threshold=10.0,
            tmb_coding_bed=self.fixtures["tmb_bed"],
            tmb_partial_threshold=0.80,
            tmb_high_threshold=10.0,
            panel_bed=self.fixtures["panel_bed"],
            cnv_amp_threshold=0.58,
            cnv_del_threshold=-1.0,
        )
        fake_config = {
            "ref": {
                "genome": "/data/ref/hg38.fa",
                "dbsnp": "/data/dbsnp/dbsnp_146.hg38.vcf.gz",
                "gnomad": "/data/gnomad/af-only-gnomad.hg38.vcf.gz",
                "pon": None,
            },
            "panel": {
                "name": "GenRichi Comprehensive Cancer Panel v1",
                "bed": self.fixtures["panel_bed"],
            },
            "annotation": {
                "vep": {
                    "genome_build": "GRCh38",
                    "extra": "--everything --canonical --cache_version 113",
                },
                "cosmic": {"vcf": "/data/cosmic/CosmicCodingMuts_v99_GRCh38.vcf.gz"},
                "clinvar": {"vcf": "/data/clinvar/clinvar_20240101.vcf.gz"},
            },
        }

        self.fake_snakemake = types.SimpleNamespace(
            input=fake_input,
            output=fake_output,
            params=fake_params,
            config=fake_config,
        )

    def _run_report_script(self):
        import runpy

        runpy.run_path(
            str(REPORT_SCRIPT),
            init_globals={"snakemake": self.fake_snakemake},
            run_name="generate_comprehensive_report",
        )

    def test_report_renders_and_embeds_provenance(self):
        self._run_report_script()

        self.assertTrue(self.html_out.exists(), "report HTML was not written")
        html = self.html_out.read_text(encoding="utf-8")

        # Existing clinical content still renders (no regression).
        self.assertIn("HCC1395_demo", html)
        self.assertIn("TP53", html)
        self.assertIn("Breast_Cancer", html)

        # Provenance block is present, well-formed, and machine-readable.
        m = re.search(
            r'<script type="application/json" id="genrichi-provenance">\s*(\{.*?\})\s*</script>',
            html,
            re.DOTALL,
        )
        self.assertIsNotNone(m, "provenance <script> block not found in report HTML")
        embedded = json.loads(m.group(1))
        self.assertEqual(embedded["sample_id"], "HCC1395_demo")
        self.assertEqual(embedded["schema_version"], "1.1")
        self.assertIn("pipeline", embedded)
        self.assertIn("generated_at_utc", embedded)

        # Human-readable footer line is present too.
        self.assertIn("Provenance recorded", html)

    def test_provenance_sidecar_json_matches_embedded_block(self):
        self._run_report_script()

        self.assertTrue(self.provenance_out.exists(), "provenance sidecar JSON was not written")
        sidecar = json.loads(self.provenance_out.read_text(encoding="utf-8"))

        html = self.html_out.read_text(encoding="utf-8")
        m = re.search(
            r'<script type="application/json" id="genrichi-provenance">\s*(\{.*?\})\s*</script>',
            html,
            re.DOTALL,
        )
        embedded = json.loads(m.group(1))

        self.assertEqual(sidecar, embedded)

    def test_missing_run_id_and_pon_do_not_break_report(self):
        # run_id is None and ref.pon is None in the fixture config above —
        # confirm this optional-metadata gap doesn't raise or omit the file.
        self._run_report_script()
        sidecar = json.loads(self.provenance_out.read_text(encoding="utf-8"))
        self.assertIsNone(sidecar["run_id"])
        self.assertIsNone(sidecar["reference_resources"]["panel_of_normals"]["path"])

    def test_report_step_records_observed_tool_versions_and_clinvar_source(self):
        # The report rule passes the runtime-probe JSON files as `input.tool_versions`;
        # the sidecar must carry what those tools reported and identify the report's ClinVar source.
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from fixture_builders import write_vcf_with_index

        clinvar = self.tmpdir / "clinvar.vcf.gz"
        write_vcf_with_index(clinvar, ["##fileformat=VCFv4.1", "##fileDate=2026-05-17", "##source=ClinVar",
                                       "##reference=GRCh38"])
        self.fake_snakemake.config["annotation"]["clinvar"]["vcf"] = str(clinvar)

        probe = self.tmpdir / "annotation_tool_versions.json"
        probe.write_text(json.dumps({
            "label": "annotation", "observed_at_utc": "2026-09-19T00:00:00+00:00",
            "conda_prefix": "/envs/annotation",
            "tools": {"vep": {"found": True, "path": "/envs/annotation/bin/vep", "version": "113.0"},
                      "bcftools": {"found": True, "path": "/x/bcftools", "version": "1.19", "htslib": "1.19.1"}},
        }), encoding="utf-8")
        self.fake_snakemake.input.tool_versions = [str(probe)]

        self._run_report_script()
        sidecar = json.loads(self.provenance_out.read_text(encoding="utf-8"))
        tools = sidecar["software_observed"]["tools"]
        self.assertEqual(tools["vep"]["version"], "113.0")
        self.assertEqual(tools["bcftools"]["htslib_version"], "1.19.1")
        report_cv = sidecar["clinical_annotation_sources"]["report_clinvar"]
        self.assertEqual(report_cv["version"], "2026-05-17")
        self.assertEqual(report_cv["role"], "report_clinical_classification_source")

    def test_report_without_tool_version_inputs_still_renders(self):
        # older callers / fixtures have no `tool_versions` input at all
        self.assertFalse(hasattr(self.fake_snakemake.input, "tool_versions"))
        self._run_report_script()
        sidecar = json.loads(self.provenance_out.read_text(encoding="utf-8"))
        self.assertEqual(sidecar["software_observed"]["status"], "not_collected")


if __name__ == "__main__":
    unittest.main()
