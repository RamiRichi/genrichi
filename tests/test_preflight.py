"""
Unit tests for workflow/scripts/preflight.py (GenRichi Phase 5.2).

Run with:
    python -m unittest tests.test_preflight -v
"""

import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

sys.path.insert(0, str(Path(__file__).resolve().parent))

import preflight  # noqa: E402
import resource_inspect  # noqa: E402
from fixture_builders import (  # noqa: E402
    write_reference, write_vcf_with_index, write_vep_cache,
)


def _valid_config(tmp: Path):
    """A complete, internally-consistent config pointing at real (tiny) fixture files."""
    genome = write_reference(tmp)
    write_vcf_with_index(tmp / "dbsnp_146.hg38.vcf.gz",
                         ["##fileformat=VCFv4.0", "##fileDate=20151104", "##source=dbSNP",
                          "##dbSNP_BUILD_ID=146", "##reference=GRCh38.p2"])
    write_vcf_with_index(tmp / "af-only-gnomad.hg38.vcf.gz", ["##fileformat=VCFv4.2"])
    write_vcf_with_index(tmp / "clinvar.vcf.gz",
                         ["##fileformat=VCFv4.1", "##fileDate=2026-05-17", "##source=ClinVar",
                          "##reference=GRCh38"])
    dbsnp, gnomad, clinvar = (tmp / "dbsnp_146.hg38.vcf.gz", tmp / "af-only-gnomad.hg38.vcf.gz",
                              tmp / "clinvar.vcf.gz")

    panel_bed = tmp / "panel.bed"
    panel_bed.write_text("chr1\t0\t4\tTP53\n", encoding="utf-8")

    vep_cache = write_vep_cache(tmp)

    (tmp / "config").mkdir(exist_ok=True)
    manifest = resource_inspect.build_reference_manifest(str(genome), "GRCh38", "test reference")
    (tmp / "config" / "reference_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")

    return {
        "samples": "config/comprehensive_samples.tsv",
        "ref": {"genome": str(genome), "dbsnp": str(dbsnp), "gnomad": str(gnomad), "pon": None},
        "panel": {"bed": str(panel_bed), "name": "Test Panel"},
        "calling": {"filter": {"min_af": 0.02, "min_depth": 20, "min_alt_reads": 3}},
        "annotation": {
            "vep": {"cache_dir": str(vep_cache), "genome_build": "GRCh38",
                    "extra": "--everything --cache_version 113"},
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
        report = preflight.run_preflight_checks(self.config, self._valid_samples_df(), base_dir=str(self.tmp))
        self.assertTrue(report.ok, msg=[(i.field, i.message) for i in report.errors])
        # COSMIC is unconfigured in _valid_config -> expect exactly one warning
        self.assertEqual(len(report.warnings), 1)
        self.assertEqual(report.warnings[0].field, "annotation.cosmic.vcf")

    def test_multiple_problems_all_surfaced(self):
        del self.config["msi"]["threshold"]  # missing key
        self.config["cnv"]["del_threshold"] = 1.0  # wrong sign
        samples_df = self._valid_samples_df()
        samples_df.loc["S1", "tumor_r1"] = ""  # blank field
        report = preflight.run_preflight_checks(self.config, samples_df, base_dir=str(self.tmp))
        self.assertFalse(report.ok)
        fields = [i.field for i in report.errors]
        self.assertIn("msi.threshold", fields)
        self.assertIn("cnv.del_threshold", fields)
        self.assertIn("samples.S1.tumor_r1", fields)

    def test_format_report_is_actionable_not_a_traceback(self):
        del self.config["panel"]["bed"]
        report = preflight.run_preflight_checks(self.config, self._valid_samples_df(), base_dir=str(self.tmp))
        text = preflight.format_report(report)
        self.assertNotIn("Traceback", text)
        self.assertNotIn("KeyError", text)
        self.assertIn("panel.bed", text)
        self.assertIn("error(s)", text)

    def test_format_report_all_clear(self):
        report = preflight.run_preflight_checks(self.config, self._valid_samples_df(), base_dir=str(self.tmp))
        # still has the cosmic warning, so it won't say "All pre-flight checks passed"
        text = preflight.format_report(report)
        self.assertIn("warning(s)", text)
        self.assertNotIn("error(s)", text)


class _TmpMixin:
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)
        self.tmp = Path(self.tmpdir.name)


class TestFastaIndexConsistency(_TmpMixin, unittest.TestCase):
    def test_consistent_reference_has_no_errors_and_reports_facts(self):
        genome = write_reference(self.tmp)
        errors, warnings, facts = resource_inspect.check_fasta_indexes(str(genome))
        self.assertEqual((errors, warnings), ([], []))
        self.assertEqual(facts["contig_count"], 1)
        self.assertEqual(facts["total_bases"], 4)
        self.assertEqual(len(facts["fai_sha256"]), 64)
        self.assertEqual(len(facts["dict_sha256"]), 64)

    def test_dict_length_disagreeing_with_fai_is_an_error(self):
        genome = write_reference(self.tmp)
        (self.tmp / "hg38.dict").write_bytes(b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:5\n")
        errors, _w, _f = resource_inspect.check_fasta_indexes(str(genome))
        self.assertTrue(any("disagree" in e for e in errors), errors)

    def test_dict_with_different_contig_count_is_an_error(self):
        genome = write_reference(self.tmp)
        (self.tmp / "hg38.dict").write_bytes(b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:4\n@SQ\tSN:chr2\tLN:4\n")
        errors, _w, _f = resource_inspect.check_fasta_indexes(str(genome))
        self.assertTrue(any("2 dictionary contigs vs 1" in e for e in errors), errors)

    def test_header_only_dict_is_an_error(self):
        genome = write_reference(self.tmp)
        (self.tmp / "hg38.dict").write_bytes(b"@HD\tVN:1.6\n")
        errors, _w, _f = resource_inspect.check_fasta_indexes(str(genome))
        self.assertTrue(any("no @SQ" in e for e in errors), errors)

    def test_truncated_fasta_is_an_error(self):
        genome = write_reference(self.tmp)
        genome.write_bytes(b">chr1\nAC")
        errors, _w, _f = resource_inspect.check_fasta_indexes(str(genome))
        self.assertTrue(any("beyond the end of the FASTA" in e for e in errors), errors)

    def test_fai_missing_a_contig_the_fasta_has_is_an_error(self):
        genome = write_reference(self.tmp)
        genome.write_bytes(b">chr1\nACGT\n>chr2\n" + b"A" * 200 + b"\n")
        errors, _w, _f = resource_inspect.check_fasta_indexes(str(genome))
        self.assertTrue(any("stale" in e for e in errors), errors)

    def test_fai_for_a_different_fasta_is_an_error(self):
        genome = write_reference(self.tmp)
        genome.write_bytes(b">chrZ\nACGT\n")
        errors, _w, _f = resource_inspect.check_fasta_indexes(str(genome))
        self.assertTrue(any("first contig" in e for e in errors), errors)

    def test_malformed_fai_is_reported_not_raised(self):
        genome = write_reference(self.tmp)
        (self.tmp / "hg38.fa.fai").write_bytes(b"chr1\tfour\t6\t4\t5\n")
        errors, _w, _f = resource_inspect.check_fasta_indexes(str(genome))
        self.assertTrue(any("Unusable FASTA index" in e for e in errors), errors)

    def test_run_preflight_surfaces_inconsistency_with_named_field(self):
        config = _valid_config(self.tmp)
        (self.tmp / "hg38.dict").write_bytes(b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:9\n")
        errors, _w, _i = preflight.validate_runtime_resources(config, base_dir=str(self.tmp))
        self.assertTrue(any(i.field == "ref.genome (FASTA/.fai/.dict consistency)" for i in errors))


class TestBwaIndex(_TmpMixin, unittest.TestCase):
    def test_complete_index_ok(self):
        genome = write_reference(self.tmp)
        errors, _w, facts = resource_inspect.check_bwa_index(str(genome), 4, 1)
        self.assertEqual(errors, [])
        self.assertEqual(facts["bwa_index_l_pac"], 4)

    def test_missing_index_file_is_an_error(self):
        genome = write_reference(self.tmp)
        (self.tmp / "hg38.fa.sa").unlink()
        errors, _w, _f = resource_inspect.check_bwa_index(str(genome), 4, 1)
        self.assertTrue(any("missing .sa" in e for e in errors), errors)

    def test_index_built_from_a_different_fasta_is_an_error(self):
        genome = write_reference(self.tmp)
        (self.tmp / "hg38.fa.ann").write_bytes(b"3209286105 455 11\n")
        errors, _w, _f = resource_inspect.check_bwa_index(str(genome), 4, 1)
        self.assertTrue(any("different FASTA" in e for e in errors), errors)

    def test_truncated_pac_is_an_error(self):
        genome = write_reference(self.tmp)
        (self.tmp / "hg38.fa.pac").write_bytes(b"\x00")
        errors, _w, _f = resource_inspect.check_bwa_index(str(genome), 4, 1)
        self.assertTrue(any(".pac has 1 bytes" in e for e in errors), errors)

    def test_empty_index_file_is_an_error(self):
        genome = write_reference(self.tmp)
        (self.tmp / "hg38.fa.bwt").write_bytes(b"")
        errors, _w, _f = resource_inspect.check_bwa_index(str(genome), 4, 1)
        self.assertTrue(any("empty" in e for e in errors), errors)


class TestVcfIndexUsability(_TmpMixin, unittest.TestCase):
    def test_bgzf_vcf_with_valid_tbi_ok(self):
        vcf = self.tmp / "a.vcf.gz"
        write_vcf_with_index(vcf, ["##fileformat=VCFv4.2"])
        self.assertEqual(resource_inspect.check_vcf_index(str(vcf)), ([], []))

    def test_plain_gzip_is_rejected(self):
        vcf = self.tmp / "a.vcf.gz"
        vcf.write_bytes(gzip.compress(b"##fileformat=VCFv4.2\n"))
        Path(str(vcf) + ".tbi").write_bytes(gzip.compress(b"TBI\x01" + b"\x00" * 8))
        errors, _w = resource_inspect.check_vcf_index(str(vcf))
        self.assertTrue(any("not BGZF" in e for e in errors), errors)

    def test_index_with_wrong_magic_is_rejected(self):
        vcf = self.tmp / "a.vcf.gz"
        write_vcf_with_index(vcf, ["##fileformat=VCFv4.2"])
        Path(str(vcf) + ".tbi").write_bytes(b"stub")
        errors, _w = resource_inspect.check_vcf_index(str(vcf))
        self.assertTrue(any("not a valid TBI" in e for e in errors), errors)

    def test_index_only_slightly_older_than_data_is_not_flagged(self):
        import os
        vcf = self.tmp / "a.vcf.gz"
        write_vcf_with_index(vcf, ["##fileformat=VCFv4.2"])
        st = os.stat(vcf)
        os.utime(str(vcf) + ".tbi", (st.st_atime, st.st_mtime - 600))  # 10 min older: normal
        self.assertEqual(resource_inspect.check_vcf_index(str(vcf)), ([], []))

    def test_index_much_older_than_data_is_a_warning_not_an_error(self):
        import os
        vcf = self.tmp / "a.vcf.gz"
        write_vcf_with_index(vcf, ["##fileformat=VCFv4.2"])
        st = os.stat(vcf)
        os.utime(str(vcf) + ".tbi", (st.st_atime, st.st_mtime - 3 * 86400))
        errors, warnings = resource_inspect.check_vcf_index(str(vcf))
        self.assertEqual(errors, [])
        self.assertTrue(any("older than" in w for w in warnings), warnings)

    def test_header_reference_contradicting_declared_build_is_an_error(self):
        config = _valid_config(self.tmp)
        write_vcf_with_index(self.tmp / "clinvar.vcf.gz",
                             ["##fileformat=VCFv4.1", "##fileDate=2026-05-17", "##reference=GRCh37"])
        errors, _w, _i = preflight.validate_runtime_resources(config, base_dir=str(self.tmp))
        self.assertTrue(any(i.field == "annotation.clinvar.vcf" and "GRCh37" in i.message for i in errors))

    def test_vcf_header_facts_are_read_from_the_file_not_the_filename(self):
        vcf = self.tmp / "clinvar.vcf.gz"  # filename carries no version at all
        write_vcf_with_index(vcf, ["##fileformat=VCFv4.1", "##fileDate=2026-05-17",
                                   "##source=ClinVar", "##reference=GRCh38"])
        header = resource_inspect.read_vcf_header(str(vcf))
        self.assertEqual(header["file_date"], "2026-05-17")
        self.assertEqual(header["source"], "ClinVar")
        self.assertEqual(header["reference"], "GRCh38")


class TestVepCacheInspection(_TmpMixin, unittest.TestCase):
    def test_pinned_version_present_reports_observed_assembly(self):
        cache = write_vep_cache(self.tmp)
        errors, warnings, facts = resource_inspect.inspect_vep_cache(str(cache), "113", "GRCh38")
        self.assertEqual((errors, warnings), ([], []))
        self.assertEqual(facts["assembly"], "GRCh38")
        self.assertEqual(facts["version_dir"], "113_GRCh38")
        self.assertEqual(facts["info"]["source_ClinVar"], "202404")
        self.assertEqual(len(facts["info_txt_sha256"]), 64)

    def test_pinned_version_not_installed_is_an_error_listing_what_is(self):
        cache = write_vep_cache(self.tmp)
        errors, _w, _f = resource_inspect.inspect_vep_cache(str(cache), "116", "GRCh38")
        self.assertEqual(len(errors), 1)
        self.assertIn("116_GRCh38", errors[0])
        self.assertIn("113_GRCh38", errors[0])

    def test_assembly_in_info_txt_contradicting_config_is_an_error(self):
        cache = write_vep_cache(self.tmp, info_assembly="GRCh37")
        errors, _w, _f = resource_inspect.inspect_vep_cache(str(cache), "113", "GRCh38")
        self.assertTrue(any("assembly mismatch" in e for e in errors), errors)

    def test_wrong_assembly_directory_is_an_error(self):
        cache = write_vep_cache(self.tmp, assembly="GRCh37")
        errors, _w, _f = resource_inspect.inspect_vep_cache(str(cache), "113", "GRCh38")
        self.assertTrue(errors)

    def test_cache_without_info_txt_is_an_error(self):
        cache = write_vep_cache(self.tmp)
        (cache / "homo_sapiens" / "113_GRCh38" / "info.txt").unlink()
        errors, _w, _f = resource_inspect.inspect_vep_cache(str(cache), "113", "GRCh38")
        self.assertTrue(any("no info.txt" in e for e in errors), errors)

    def test_unpinned_single_cache_is_a_warning_not_an_error(self):
        cache = write_vep_cache(self.tmp)
        errors, warnings, facts = resource_inspect.inspect_vep_cache(str(cache), None, "GRCh38")
        self.assertEqual(errors, [])
        self.assertTrue(any("not pinned" in w for w in warnings), warnings)
        self.assertEqual(facts["version_dir"], "113_GRCh38")

    def test_unpinned_with_several_caches_is_an_error(self):
        write_vep_cache(self.tmp, version="113")
        d = self.tmp / "vep_cache" / "homo_sapiens" / "115_GRCh38"
        d.mkdir()
        (d / "info.txt").write_text("assembly\tGRCh38\n", encoding="utf-8")
        errors, _w, _f = resource_inspect.inspect_vep_cache(str(self.tmp / "vep_cache"), None, "GRCh38")
        self.assertTrue(any("several" in e for e in errors), errors)


class TestReferenceManifest(_TmpMixin, unittest.TestCase):
    def _run(self, config):
        return preflight.validate_runtime_resources(config, base_dir=str(self.tmp))

    def test_matching_manifest_reports_a_match_and_no_warning(self):
        config = _valid_config(self.tmp)
        errors, warnings, info = self._run(config)
        self.assertEqual(errors, [])
        self.assertEqual(warnings, [])
        self.assertTrue(any("matches manifest" in i.message for i in info), [i.message for i in info])

    def test_reference_that_differs_from_the_manifest_is_an_error(self):
        config = _valid_config(self.tmp)
        # a different but internally consistent reference (2 contigs) at the same path
        genome = Path(config["ref"]["genome"])
        genome.write_bytes(b">chr1\nACGT\n>chr2\nACGT\n")
        (self.tmp / "hg38.fa.fai").write_bytes(b"chr1\t4\t6\t4\t5\nchr2\t4\t17\t4\t5\n")
        (self.tmp / "hg38.dict").write_bytes(b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:4\n@SQ\tSN:chr2\tLN:4\n")
        (self.tmp / "hg38.fa.ann").write_bytes(b"8 2 11\n")
        (self.tmp / "hg38.fa.pac").write_bytes(b"\x00\x00\x00")
        errors, _w, _i = self._run(config)
        self.assertTrue(any(i.field == "ref.genome (expected resource set)" for i in errors),
                        [(i.field, i.message) for i in errors])

    def test_missing_default_manifest_is_a_warning_not_an_error(self):
        config = _valid_config(self.tmp)
        (self.tmp / "config" / "reference_manifest.json").unlink()
        errors, warnings, _i = self._run(config)
        self.assertEqual(errors, [])
        self.assertTrue(any("No reference manifest" in w.message for w in warnings))

    def test_explicitly_configured_but_missing_manifest_is_an_error(self):
        config = _valid_config(self.tmp)
        config["ref"]["manifest"] = str(self.tmp / "nope.json")
        errors, _w, _i = self._run(config)
        self.assertTrue(any(i.field == "ref.manifest" for i in errors))

    def test_manifest_genome_build_contradicting_config_is_an_error(self):
        config = _valid_config(self.tmp)
        path = self.tmp / "config" / "reference_manifest.json"
        manifest = json.loads(path.read_text(encoding="utf-8"))
        manifest["genome_build"] = "GRCh37"
        path.write_text(json.dumps(manifest), encoding="utf-8")
        errors, _w, _i = self._run(config)
        self.assertTrue(any("expected resource set" in i.field for i in errors))

    def test_manifest_cli_round_trip(self):
        genome = write_reference(self.tmp)
        out = self.tmp / "m.json"
        self.assertEqual(preflight._main(["--write-reference-manifest", str(genome), str(out)]), 0)
        manifest = json.loads(out.read_text(encoding="utf-8"))
        self.assertEqual(manifest["reference"]["contig_count"], 1)
        self.assertEqual(resource_inspect.compare_to_manifest(
            manifest, resource_inspect.check_fasta_indexes(str(genome))[2], "GRCh38"), [])


class TestObservedFactsInPreflightReport(_TmpMixin, unittest.TestCase):
    def _samples(self):
        real = self.tmp / "reads.fastq.gz"
        real.write_bytes(b"not empty")
        return pd.DataFrame([{"sample_id": "S1", "tumor_r1": str(real), "tumor_r2": str(real),
                              "normal_r1": str(real), "normal_r2": str(real)}]).set_index("sample_id")

    def test_observed_facts_are_info_not_warnings_or_errors(self):
        config = _valid_config(self.tmp)
        report = preflight.run_preflight_checks(config, self._samples(), base_dir=str(self.tmp))
        self.assertTrue(report.ok)
        self.assertEqual([w.field for w in report.warnings], ["annotation.cosmic.vcf"])
        fields = [i.field for i in report.info]
        for expected in ("ref.genome", "ref.dbsnp", "annotation.clinvar.vcf", "annotation.vep.cache_dir", "runtime.tools"):
            self.assertIn(expected, fields)

    def test_report_text_shows_observed_versions_from_headers_and_cache(self):
        config = _valid_config(self.tmp)
        text = preflight.format_report(
            preflight.run_preflight_checks(config, self._samples(), base_dir=str(self.tmp)))
        self.assertIn("file_date=2026-05-17", text)
        self.assertIn("dbsnp_build_id=146", text)
        self.assertIn("113_GRCh38", text)
        self.assertIn("ClinVar 202404", text)
        self.assertIn("not the report's ClinVar source", text)

    def test_known_optional_cosmic_warning_stays_explicit(self):
        config = _valid_config(self.tmp)
        config["annotation"]["cosmic"]["vcf"] = str(self.tmp / "CosmicCodingMuts.vcf.gz")
        report = preflight.run_preflight_checks(config, self._samples(), base_dir=str(self.tmp))
        self.assertTrue(report.ok)
        cosmic = [w for w in report.warnings if w.field == "annotation.cosmic.vcf"]
        self.assertEqual(len(cosmic), 1)
        self.assertIn("COSMIC annotation will be skipped", cosmic[0].message)

    def test_deep_checks_can_be_disabled(self):
        config = _valid_config(self.tmp)
        (self.tmp / "hg38.dict").write_bytes(b"@HD\tVN:1.6\n")
        self.assertFalse(preflight.run_preflight_checks(config, self._samples(), base_dir=str(self.tmp)).ok)
        self.assertTrue(preflight.run_preflight_checks(
            config, self._samples(), base_dir=str(self.tmp), deep_resource_checks=False).ok)

    def test_missing_resources_are_not_double_reported(self):
        config = _valid_config(self.tmp)
        config["ref"]["genome"] = str(self.tmp / "gone.fa")
        report = preflight.run_preflight_checks(config, self._samples(), base_dir=str(self.tmp))
        genome_msgs = [e for e in report.errors if e.field.startswith("ref.genome")]
        self.assertEqual(len(genome_msgs), 1, [(e.field, e.message) for e in genome_msgs])


if __name__ == "__main__":
    unittest.main()
