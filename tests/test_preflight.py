"""
Unit tests for workflow/scripts/preflight.py (GenRichi Phase 5.2).

Run with:
    python -m unittest tests.test_preflight -v
"""

import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import preflight  # noqa: E402


def _valid_config(tmp: Path):
    """A minimally complete, internally-consistent config pointing at real (tiny) fixture files."""
    genome = tmp / "hg38.fa"
    genome.write_text(">chr1\nACGT\n", encoding="utf-8")
    (tmp / "hg38.fa.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    (tmp / "hg38.dict").write_text("@HD\tVN:1.6\n", encoding="utf-8")

    dbsnp = tmp / "dbsnp_146.hg38.vcf.gz"
    dbsnp.write_bytes(b"stub")
    (tmp / "dbsnp_146.hg38.vcf.gz.tbi").write_bytes(b"stub")

    gnomad = tmp / "af-only-gnomad.hg38.vcf.gz"
    gnomad.write_bytes(b"stub")
    (tmp / "af-only-gnomad.hg38.vcf.gz.tbi").write_bytes(b"stub")

    clinvar = tmp / "clinvar.vcf.gz"
    clinvar.write_bytes(b"stub")
    (tmp / "clinvar.vcf.gz.tbi").write_bytes(b"stub")

    panel_bed = tmp / "panel.bed"
    panel_bed.write_text("chr1\t0\t4\tTP53\n", encoding="utf-8")

    vep_cache = tmp / "vep_cache"
    vep_cache.mkdir()

    return {
        "samples": "config/comprehensive_samples.tsv",
        "ref": {"genome": str(genome), "dbsnp": str(dbsnp), "gnomad": str(gnomad), "pon": None},
        "panel": {"bed": str(panel_bed), "name": "Test Panel"},
        "calling": {"filter": {"min_af": 0.02, "min_depth": 20, "min_alt_reads": 3}},
        "annotation": {
            "vep": {"cache_dir": str(vep_cache), "genome_build": "GRCh38"},
            "clinvar": {"vcf": str(clinvar)},
            "cosmic": {"vcf": None},
        },
        "cnv": {"amp_threshold": 0.58, "del_threshold": -1.0, "min_probes": 5},
        "msi": {"threshold": 10.0},
        "tmb": {"coding_mb": 1.5, "coding_bed": "resources/panel/comprehensive_genes_cds.bed",
                "partial_threshold": 0.8, "high_threshold": 10.0},
        "report": {"company": "GenRichi"},
    }


def _samples_df(rows, columns=("tumor_r1", "tumor_r2", "normal_r1", "normal_r2")):
    df = pd.DataFrame(rows).set_index("sample_id")
    for c in columns:
        if c not in df.columns:
            pass
    return df


class TestConfigKeyValidation(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)
        self.config = _valid_config(Path(self.tmpdir.name))

    def test_all_required_keys_present_no_issues(self):
        issues = preflight.validate_config_keys(self.config)
        self.assertEqual(issues, [])

    def test_missing_key_reported_by_name(self):
        del self.config["panel"]["bed"]
        issues = preflight.validate_config_keys(self.config)
        fields = [i.field for i in issues]
        self.assertIn("panel.bed", fields)

    def test_missing_top_level_section_reported(self):
        del self.config["msi"]
        issues = preflight.validate_config_keys(self.config)
        fields = [i.field for i in issues]
        self.assertIn("msi.threshold", fields)

    def test_non_dict_config_does_not_raise(self):
        issues = preflight.validate_config_keys(None)
        self.assertTrue(len(issues) >= 1)


class TestConfigRangeValidation(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)
        self.config = _valid_config(Path(self.tmpdir.name))

    def test_valid_values_no_issues(self):
        self.assertEqual(preflight.validate_config_ranges(self.config), [])

    def test_min_af_out_of_range(self):
        self.config["calling"]["filter"]["min_af"] = 5.0
        issues = preflight.validate_config_ranges(self.config)
        self.assertTrue(any(i.field == "calling.filter.min_af" for i in issues))

    def test_del_threshold_wrong_sign(self):
        self.config["cnv"]["del_threshold"] = 1.0  # should be negative
        issues = preflight.validate_config_ranges(self.config)
        self.assertTrue(any(i.field == "cnv.del_threshold" for i in issues))

    def test_unknown_genome_build_flagged(self):
        self.config["annotation"]["vep"]["genome_build"] = "T2T-CHM13"
        issues = preflight.validate_config_ranges(self.config)
        self.assertTrue(any(i.field == "annotation.vep.genome_build" for i in issues))

    def test_non_numeric_value_flagged_not_raised(self):
        self.config["msi"]["threshold"] = "not-a-number"
        issues = preflight.validate_config_ranges(self.config)
        self.assertTrue(any(i.field == "msi.threshold" for i in issues))

    def test_missing_key_does_not_duplicate_range_issue(self):
        del self.config["msi"]["threshold"]
        issues = preflight.validate_config_ranges(self.config)
        self.assertFalse(any(i.field == "msi.threshold" for i in issues))


class TestSamplesSheetValidation(unittest.TestCase):
    def test_missing_required_column(self):
        df = pd.DataFrame([
            {"sample_id": "S1", "tumor_r1": "a", "tumor_r2": "b", "normal_r1": "c"}
        ]).set_index("sample_id")
        issues = preflight.validate_samples_sheet(df)
        self.assertTrue(any(i.field == "samples.columns" for i in issues))
        self.assertIn("normal_r2", issues[0].message)

    def test_blank_field_reported_with_sample_and_column(self):
        df = pd.DataFrame([
            {"sample_id": "S1", "tumor_r1": "a", "tumor_r2": "b", "normal_r1": "", "normal_r2": "d"}
        ]).set_index("sample_id")
        issues = preflight.validate_samples_sheet(df)
        self.assertTrue(any(i.field == "samples.S1.normal_r1" for i in issues))

    def test_nan_field_reported(self):
        df = pd.DataFrame([
            {"sample_id": "S1", "tumor_r1": "a", "tumor_r2": None, "normal_r1": "c", "normal_r2": "d"}
        ]).set_index("sample_id")
        issues = preflight.validate_samples_sheet(df)
        self.assertTrue(any(i.field == "samples.S1.tumor_r2" for i in issues))

    def test_duplicate_sample_id_reported(self):
        df = pd.DataFrame([
            {"sample_id": "S1", "tumor_r1": "a", "tumor_r2": "b", "normal_r1": "c", "normal_r2": "d"},
            {"sample_id": "S1", "tumor_r1": "e", "tumor_r2": "f", "normal_r1": "g", "normal_r2": "h"},
        ]).set_index("sample_id")
        issues = preflight.validate_samples_sheet(df)
        self.assertTrue(any(i.field == "samples.sample_id" and "Duplicate" in i.message for i in issues))

    def test_valid_sheet_no_issues(self):
        df = pd.DataFrame([
            {"sample_id": "S1", "tumor_r1": "a", "tumor_r2": "b", "normal_r1": "c", "normal_r2": "d"},
            {"sample_id": "S2", "tumor_r1": "e", "tumor_r2": "f", "normal_r1": "g", "normal_r2": "h"},
        ]).set_index("sample_id")
        self.assertEqual(preflight.validate_samples_sheet(df), [])


class TestFastqValidation(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)
        self.tmp = Path(self.tmpdir.name)

    def test_missing_fastq_reported(self):
        df = pd.DataFrame([{
            "sample_id": "S1",
            "tumor_r1": str(self.tmp / "nope_R1.fastq.gz"),
            "tumor_r2": str(self.tmp / "nope_R2.fastq.gz"),
            "normal_r1": str(self.tmp / "nope_N1.fastq.gz"),
            "normal_r2": str(self.tmp / "nope_N2.fastq.gz"),
        }]).set_index("sample_id")
        issues = preflight.validate_fastq_files(df)
        self.assertEqual(len(issues), 4)
        self.assertTrue(all("not found" in i.message for i in issues))

    def test_empty_fastq_reported(self):
        empty = self.tmp / "empty_R1.fastq.gz"
        empty.write_bytes(b"")
        real = self.tmp / "real_R2.fastq.gz"
        real.write_bytes(b"not empty")
        df = pd.DataFrame([{
            "sample_id": "S1", "tumor_r1": str(empty), "tumor_r2": str(real),
            "normal_r1": str(real), "normal_r2": str(real),
        }]).set_index("sample_id")
        issues = preflight.validate_fastq_files(df)
        self.assertEqual(len(issues), 1)
        self.assertIn("empty", issues[0].message)
        self.assertEqual(issues[0].field, "samples.S1.tumor_r1")

    def test_valid_fastqs_no_issues(self):
        real = self.tmp / "real.fastq.gz"
        real.write_bytes(b"not empty")
        df = pd.DataFrame([{
            "sample_id": "S1", "tumor_r1": str(real), "tumor_r2": str(real),
            "normal_r1": str(real), "normal_r2": str(real),
        }]).set_index("sample_id")
        self.assertEqual(preflight.validate_fastq_files(df), [])

    def test_blank_field_not_double_reported(self):
        real = self.tmp / "real.fastq.gz"
        real.write_bytes(b"not empty")
        df = pd.DataFrame([{
            "sample_id": "S1", "tumor_r1": "", "tumor_r2": str(real),
            "normal_r1": str(real), "normal_r2": str(real),
        }]).set_index("sample_id")
        # blank fields are validate_samples_sheet's job, not FASTQ existence's
        self.assertEqual(preflight.validate_fastq_files(df), [])


class TestResourceValidation(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)
        self.config = _valid_config(Path(self.tmpdir.name))

    def test_fully_valid_resources_no_issues(self):
        self.assertEqual(preflight.validate_resources(self.config), [])

    def test_missing_genome_reported(self):
        self.config["ref"]["genome"] = str(Path(self.tmpdir.name) / "does_not_exist.fa")
        issues = preflight.validate_resources(self.config)
        self.assertTrue(any(i.field == "ref.genome" for i in issues))

    def test_missing_fai_reported_separately_from_missing_genome(self):
        genome = Path(self.tmpdir.name) / "hg38.fa"
        (Path(self.tmpdir.name) / "hg38.fa.fai").unlink()
        self.config["ref"]["genome"] = str(genome)
        issues = preflight.validate_resources(self.config)
        self.assertTrue(any(".fai" in i.field for i in issues))
        self.assertFalse(any(i.field == "ref.genome" for i in issues))

    def test_missing_vep_cache_dir_reported(self):
        self.config["annotation"]["vep"]["cache_dir"] = str(Path(self.tmpdir.name) / "no_such_dir")
        issues = preflight.validate_resources(self.config)
        self.assertTrue(any(i.field == "annotation.vep.cache_dir" for i in issues))

    def test_unconfigured_pon_is_not_an_error(self):
        self.config["ref"]["pon"] = None
        self.assertEqual(preflight.validate_resources(self.config), [])

    def test_configured_but_missing_pon_is_an_error(self):
        self.config["ref"]["pon"] = str(Path(self.tmpdir.name) / "missing_pon.vcf.gz")
        issues = preflight.validate_resources(self.config)
        self.assertTrue(any(i.field == "ref.pon" for i in issues))


class TestPanelReferenceCompatibility(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)
        self.tmp = Path(self.tmpdir.name)

    def test_matching_chr_prefix_no_issue(self):
        bed = self.tmp / "panel.bed"
        bed.write_text("chr1\t0\t10\tGENE\n", encoding="utf-8")
        genome = self.tmp / "ref.fa"
        genome.write_text(">chr1\nACGT\n", encoding="utf-8")
        (self.tmp / "ref.fa.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
        issues = preflight.check_panel_reference_compatibility(str(bed), str(genome))
        self.assertEqual(issues, [])

    def test_mismatched_chr_prefix_reported(self):
        bed = self.tmp / "panel.bed"
        bed.write_text("chr1\t0\t10\tGENE\n", encoding="utf-8")
        genome = self.tmp / "ref.fa"
        genome.write_text(">1\nACGT\n", encoding="utf-8")
        (self.tmp / "ref.fa.fai").write_text("1\t4\t6\t4\t5\n", encoding="utf-8")
        issues = preflight.check_panel_reference_compatibility(str(bed), str(genome))
        self.assertEqual(len(issues), 1)
        self.assertIn("naming", issues[0].message)

    def test_missing_fai_skips_check_gracefully(self):
        bed = self.tmp / "panel.bed"
        bed.write_text("chr1\t0\t10\tGENE\n", encoding="utf-8")
        genome = self.tmp / "ref.fa"
        genome.write_text(">chr1\nACGT\n", encoding="utf-8")
        # no .fai written
        issues = preflight.check_panel_reference_compatibility(str(bed), str(genome))
        self.assertEqual(issues, [])


class TestGenomeBuildNaming(unittest.TestCase):
    def test_consistent_grch38_no_issue(self):
        config = {"ref": {"genome": "/data/ref/hg38.fa"}, "annotation": {"vep": {"genome_build": "GRCh38"}}}
        self.assertEqual(preflight.check_genome_build_naming(config), [])

    def test_grch38_build_with_hg19_path_flagged(self):
        config = {"ref": {"genome": "/data/ref/hg19.fa"}, "annotation": {"vep": {"genome_build": "GRCh38"}}}
        issues = preflight.check_genome_build_naming(config)
        self.assertEqual(len(issues), 1)

    def test_ambiguous_path_not_flagged(self):
        config = {"ref": {"genome": "/data/ref/genome.fa"}, "annotation": {"vep": {"genome_build": "GRCh38"}}}
        self.assertEqual(preflight.check_genome_build_naming(config), [])


class TestCosmicAvailability(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)

    def test_not_configured(self):
        available, path = preflight.check_cosmic_availability({"annotation": {"cosmic": {"vcf": None}}})
        self.assertFalse(available)
        self.assertIsNone(path)

    def test_configured_but_missing_file(self):
        missing = str(Path(self.tmpdir.name) / "cosmic.vcf.gz")
        available, path = preflight.check_cosmic_availability({"annotation": {"cosmic": {"vcf": missing}}})
        self.assertFalse(available)
        self.assertEqual(path, missing)

    def test_configured_and_present(self):
        real = Path(self.tmpdir.name) / "cosmic.vcf.gz"
        real.write_bytes(b"stub")
        available, path = preflight.check_cosmic_availability({"annotation": {"cosmic": {"vcf": str(real)}}})
        self.assertTrue(available)


class TestRunPreflightChecksEndToEnd(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)
        self.tmp = Path(self.tmpdir.name)
        self.config = _valid_config(self.tmp)

    def _valid_samples_df(self):
        real = self.tmp / "reads.fastq.gz"
        real.write_bytes(b"not empty")
        return pd.DataFrame([{
            "sample_id": "S1", "tumor_r1": str(real), "tumor_r2": str(real),
            "normal_r1": str(real), "normal_r2": str(real),
        }]).set_index("sample_id")

    def test_fully_valid_run_has_no_errors(self):
        report = preflight.run_preflight_checks(self.config, self._valid_samples_df())
        self.assertTrue(report.ok, msg=[(i.field, i.message) for i in report.errors])
        # COSMIC is unconfigured in _valid_config -> expect exactly one warning
        self.assertEqual(len(report.warnings), 1)
        self.assertEqual(report.warnings[0].field, "annotation.cosmic.vcf")

    def test_multiple_problems_all_surfaced(self):
        del self.config["msi"]["threshold"]  # missing key
        self.config["cnv"]["del_threshold"] = 1.0  # wrong sign
        samples_df = self._valid_samples_df()
        samples_df.loc["S1", "tumor_r1"] = ""  # blank field
        report = preflight.run_preflight_checks(self.config, samples_df)
        self.assertFalse(report.ok)
        fields = [i.field for i in report.errors]
        self.assertIn("msi.threshold", fields)
        self.assertIn("cnv.del_threshold", fields)
        self.assertIn("samples.S1.tumor_r1", fields)

    def test_format_report_is_actionable_not_a_traceback(self):
        del self.config["panel"]["bed"]
        report = preflight.run_preflight_checks(self.config, self._valid_samples_df())
        text = preflight.format_report(report)
        self.assertNotIn("Traceback", text)
        self.assertNotIn("KeyError", text)
        self.assertIn("panel.bed", text)
        self.assertIn("error(s)", text)

    def test_format_report_all_clear(self):
        report = preflight.run_preflight_checks(self.config, self._valid_samples_df())
        # still has the cosmic warning, so it won't say "All pre-flight checks passed"
        text = preflight.format_report(report)
        self.assertIn("warning(s)", text)
        self.assertNotIn("error(s)", text)


if __name__ == "__main__":
    unittest.main()
