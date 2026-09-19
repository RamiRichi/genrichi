"""
Unit tests for workflow/scripts/provenance.py (GenRichi Phase 5.1).

Run with:
    python -m unittest tests.test_provenance -v
"""

import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import provenance  # noqa: E402
from fixture_builders import (  # noqa: E402
    write_reference, write_vcf_with_index, write_vep_cache,
)


GIT_AVAILABLE = shutil.which("git") is not None


class TestExtractVersionToken(unittest.TestCase):
    def test_numeric_version_with_build_suffix(self):
        self.assertEqual(
            provenance.extract_version_token("dbsnp_146.hg38.vcf.gz"), "146"
        )

    def test_v_prefixed_version(self):
        self.assertEqual(
            provenance.extract_version_token("CosmicCodingMuts_v99_GRCh38.vcf.gz"),
            "v99",
        )

    def test_date_based_version(self):
        self.assertEqual(
            provenance.extract_version_token("clinvar_20240101.vcf.gz"), "20240101"
        )

    def test_no_version_token_falls_back_to_filename(self):
        # Real GenRichi config uses these exact filenames with no embedded version —
        # must not invent one; the filename itself is the honest fallback.
        self.assertEqual(
            provenance.extract_version_token("af-only-gnomad.hg38.vcf.gz"),
            "af-only-gnomad.hg38.vcf.gz",
        )
        self.assertEqual(
            provenance.extract_version_token("clinvar.vcf.gz"), "clinvar.vcf.gz"
        )

    def test_none_path_returns_none(self):
        self.assertIsNone(provenance.extract_version_token(None))
        self.assertIsNone(provenance.extract_version_token(""))


class TestHashing(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)

    def test_sha256_of_file_matches_known_digest(self):
        p = Path(self.tmpdir.name) / "sample.bed"
        p.write_text("chr17\t100\t200\tTP53\n", encoding="utf-8")
        import hashlib

        expected = hashlib.sha256(p.read_bytes()).hexdigest()
        self.assertEqual(provenance.sha256_of_file(str(p)), expected)

    def test_sha256_of_file_deterministic_across_calls(self):
        p = Path(self.tmpdir.name) / "sample.bed"
        p.write_text("chr17\t100\t200\tTP53\n", encoding="utf-8")
        self.assertEqual(
            provenance.sha256_of_file(str(p)), provenance.sha256_of_file(str(p))
        )

    def test_sha256_of_file_missing_path_returns_none(self):
        missing = Path(self.tmpdir.name) / "does_not_exist.bed"
        self.assertIsNone(provenance.sha256_of_file(str(missing)))

    def test_sha256_of_file_none_path_returns_none(self):
        self.assertIsNone(provenance.sha256_of_file(None))

    def test_sha256_of_config_key_order_independent(self):
        a = {"panel": {"name": "X"}, "ref": {"genome": "g.fa"}}
        b = {"ref": {"genome": "g.fa"}, "panel": {"name": "X"}}
        self.assertEqual(provenance.sha256_of_config(a), provenance.sha256_of_config(b))

    def test_sha256_of_config_changes_with_content(self):
        a = {"ref": {"genome": "g.fa"}}
        b = {"ref": {"genome": "different.fa"}}
        self.assertNotEqual(provenance.sha256_of_config(a), provenance.sha256_of_config(b))


class TestGitHelpers(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)

    @unittest.skipUnless(GIT_AVAILABLE, "git not available")
    def test_get_git_commit_matches_real_repo_head(self):
        expected = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=str(provenance.REPO_ROOT),
            capture_output=True,
            text=True,
        ).stdout.strip()
        self.assertTrue(expected)
        self.assertEqual(provenance.get_git_commit(), expected)
        self.assertEqual(len(provenance.get_git_commit()), 40)

    @unittest.skipUnless(GIT_AVAILABLE, "git not available")
    def test_get_git_describe_returns_nonempty_string_in_real_repo(self):
        described = provenance.get_git_describe()
        self.assertIsInstance(described, str)
        self.assertTrue(described)

    def test_get_git_commit_none_outside_a_repo(self):
        # A fresh empty directory is not a git repository — must degrade to
        # None, not raise.
        self.assertIsNone(provenance.get_git_commit(repo_root=self.tmpdir.name))

    def test_get_git_describe_none_outside_a_repo(self):
        self.assertIsNone(provenance.get_git_describe(repo_root=self.tmpdir.name))

    def test_get_git_dirty_none_outside_a_repo(self):
        self.assertIsNone(provenance.get_git_dirty(repo_root=self.tmpdir.name))

    def test_git_helpers_survive_nonexistent_directory(self):
        ghost = str(Path(self.tmpdir.name) / "does_not_exist_at_all")
        self.assertIsNone(provenance.get_git_commit(repo_root=ghost))
        self.assertIsNone(provenance.get_git_describe(repo_root=ghost))
        self.assertIsNone(provenance.get_git_dirty(repo_root=ghost))


class TestBuildProvenance(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)

    def test_missing_optional_metadata_does_not_raise(self):
        # Empty config, no git repo, no panel BED, no run_id — every optional
        # field must degrade gracefully rather than raising or breaking output.
        record = provenance.build_provenance(
            sample_id="SAMPLE01",
            config={},
            run_id=None,
            repo_root=self.tmpdir.name,
            panel_bed_path=None,
        )
        self.assertEqual(record["sample_id"], "SAMPLE01")
        self.assertIsNone(record["run_id"])
        self.assertIsNone(record["pipeline"]["git_commit_sha"])
        self.assertIsNone(record["pipeline"]["git_dirty"])
        self.assertIsNone(record["panel"]["name"])
        self.assertIsNone(record["panel"]["bed_path"])
        self.assertIsNone(record["panel"]["bed_sha256"])
        self.assertIsNone(record["reference_genome"]["build"])
        self.assertIsNone(record["reference_resources"]["dbsnp"]["path"])
        self.assertIsNone(record["annotation_databases"]["vep_cache_version"])
        # required/deterministic fields still present
        self.assertEqual(record["schema_version"], "1.1")
        self.assertIn("generated_at_utc", record)
        self.assertIsNotNone(record["config_sha256"])

    def test_build_provenance_with_full_config_extracts_real_values(self):
        config = {
            "ref": {
                "genome": "/data/ref/hg38.fa",
                "dbsnp": "/data/dbsnp/dbsnp_146.hg38.vcf.gz",
                "gnomad": "/data/gnomad/af-only-gnomad.hg38.vcf.gz",
                "pon": None,
            },
            "panel": {"name": "GenRichi Comprehensive Cancer Panel v1", "bed": None},
            "annotation": {
                "vep": {
                    "genome_build": "GRCh38",
                    "extra": "--everything --cache_version 113",
                },
                "cosmic": {"vcf": "/data/cosmic/CosmicCodingMuts_v99_GRCh38.vcf.gz"},
                "clinvar": {"vcf": "/data/clinvar/clinvar_20240101.vcf.gz"},
            },
        }
        bed = Path(self.tmpdir.name) / "panel.bed"
        bed.write_text("chr17\t100\t200\tTP53\n", encoding="utf-8")

        record = provenance.build_provenance(
            sample_id="HCC1395_demo",
            config=config,
            run_id="RUN-42",
            repo_root=self.tmpdir.name,
            panel_bed_path=str(bed),
        )
        self.assertEqual(record["run_id"], "RUN-42")
        self.assertEqual(record["reference_genome"]["build"], "GRCh38")
        self.assertEqual(record["reference_resources"]["dbsnp"]["version"], "146")
        self.assertEqual(record["annotation_databases"]["cosmic"]["version"], "v99")
        self.assertEqual(record["annotation_databases"]["clinvar"]["version"], "20240101")
        self.assertEqual(record["annotation_databases"]["vep_cache_version"], "113")
        self.assertIsNotNone(record["panel"]["bed_sha256"])

    def test_build_provenance_is_json_serializable(self):
        record = provenance.build_provenance(
            sample_id="S1", config={"a": 1}, repo_root=self.tmpdir.name
        )
        json.dumps(record)  # must not raise


class TestHtmlEmbedding(unittest.TestCase):
    def test_provenance_script_tag_contains_valid_json(self):
        record = {"schema_version": "1.0", "sample_id": "S1", "run_id": None}
        tag = provenance.provenance_script_tag(record)
        self.assertIn('<script type="application/json" id="genrichi-provenance">', tag)
        self.assertIn("</script>", tag)
        inner = tag.split(">", 1)[1].rsplit("</script>", 1)[0].strip()
        parsed = json.loads(inner)
        self.assertEqual(parsed, record)

    def test_provenance_footer_line_handles_missing_commit(self):
        record = {
            "pipeline": {"version": None, "git_commit_sha": None},
            "generated_at_utc": "2026-09-18T00:00:00+00:00",
        }
        line = provenance.provenance_footer_line(record)
        self.assertIn("unknown", line)
        self.assertIn("2026-09-18T00:00:00+00:00", line)


def _observed_config(tmp: Path):
    """Config pointing at structurally real tiny resources; the ClinVar filename deliberately carries no version."""
    genome = write_reference(tmp)
    write_vcf_with_index(tmp / "dbsnp_146.hg38.vcf.gz",
                         ["##fileformat=VCFv4.0", "##fileDate=20151104", "##dbSNP_BUILD_ID=146",
                          "##reference=GRCh38.p2"])
    write_vcf_with_index(tmp / "gnomad.vcf.gz", ["##fileformat=VCFv4.2"])
    write_vcf_with_index(tmp / "clinvar.vcf.gz",
                         ["##fileformat=VCFv4.1", "##fileDate=2026-05-17", "##source=ClinVar",
                          "##reference=GRCh38"])
    cache = write_vep_cache(tmp)
    bed = tmp / "panel.bed"
    bed.write_text("chr1\t0\t4\tTP53\n", encoding="utf-8")
    cds = tmp / "cds.bed"
    cds.write_text("chr1\t0\t4\tTP53\n", encoding="utf-8")
    return {
        "ref": {"genome": str(genome), "dbsnp": str(tmp / "dbsnp_146.hg38.vcf.gz"),
                "gnomad": str(tmp / "gnomad.vcf.gz"), "pon": None},
        "panel": {"name": "Test Panel", "bed": str(bed)},
        "tmb": {"coding_bed": str(cds)},
        "annotation": {
            "vep": {"cache_dir": str(cache), "genome_build": "GRCh38",
                    "extra": "--everything --cache_version 113"},
            "cosmic": {"vcf": None},
            "clinvar": {"vcf": str(tmp / "clinvar.vcf.gz")},
        },
    }


def _tool_file(path: Path, label, tools, vep_check=None):
    record = {"schema_version": "1.0", "label": label, "observed_at_utc": "2026-09-19T00:00:00+00:00",
              "conda_prefix": f"/envs/{label}", "tools": tools}
    if vep_check:
        record["vep_cache_check"] = vep_check
    path.write_text(json.dumps(record), encoding="utf-8")
    return str(path)


def _all_tool_files(tmp: Path):
    return [
        _tool_file(tmp / "annotation.json", "annotation", {
            "vep": {"found": True, "path": "/envs/annotation/bin/vep", "version": "113.0",
                    "detail": {"ensembl-vep": "113.0"}},
            "bcftools": {"found": True, "path": "/envs/annotation/bin/bcftools", "version": "1.19", "htslib": "1.19.1"},
        }, vep_check={"status": "PASS", "message": "VEP 113.0 matches cache 113_GRCh38"}),
        _tool_file(tmp / "alignment.json", "alignment", {
            "bwa": {"found": True, "path": "/envs/alignment/bin/bwa", "version": "0.7.17-r1188"},
            "samtools": {"found": True, "path": "/envs/alignment/bin/samtools", "version": "1.19.2", "htslib": "1.21"},
        }),
        _tool_file(tmp / "qc.json", "qc", {
            "multiqc": {"found": True, "path": "/envs/qc/bin/multiqc", "version": "1.34"},
        }),
    ]


class TestObservedProvenance(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)
        self.tmp = Path(self.tmpdir.name)
        self.config = _observed_config(self.tmp)
        self.repo = self.tmp / "repo"  # not a git repo
        self.repo.mkdir()

    def _build(self, **kw):
        return provenance.build_provenance("S1", self.config, repo_root=str(self.repo), **kw)

    # ── ClinVar ──────────────────────────────────────────────────────────
    def test_report_clinvar_version_comes_from_the_vcf_header_not_the_filename(self):
        prov = self._build()
        report_cv = prov["clinical_annotation_sources"]["report_clinvar"]
        self.assertEqual(report_cv["role"], "report_clinical_classification_source")
        self.assertEqual(report_cv["version"], "2026-05-17")
        self.assertEqual(report_cv["version_basis"], "vcf_header_fileDate")
        self.assertEqual(report_cv["vcf_header"]["source"], "ClinVar")
        self.assertIn("ClinVar_CLNSIG", report_cv["used_for"])
        # the configured (filename-derived) token is still recorded, and is honestly just the filename
        self.assertEqual(prov["annotation_databases"]["clinvar"]["version"], "clinvar.vcf.gz")

    def test_report_clinvar_has_sha256_of_the_actual_file(self):
        import hashlib
        prov = self._build()
        expected = hashlib.sha256((self.tmp / "clinvar.vcf.gz").read_bytes()).hexdigest()
        self.assertEqual(prov["clinical_annotation_sources"]["report_clinvar"]["sha256"], expected)
        self.assertEqual(prov["checksums"]["clinvar_vcf"]["sha256"], expected)

    def test_vep_cache_clinvar_is_recorded_separately_as_a_dependency(self):
        prov = self._build()
        cache_cv = prov["clinical_annotation_sources"]["vep_cache_clinvar"]
        self.assertEqual(cache_cv["role"], "annotation_dependency_not_report_classification_source")
        self.assertEqual(cache_cv["version"], "202404")
        self.assertNotEqual(cache_cv["version"], prov["clinical_annotation_sources"]["report_clinvar"]["version"])

    def test_missing_fileDate_is_reported_as_unavailable_not_inferred(self):
        write_vcf_with_index(self.tmp / "clinvar.vcf.gz", ["##fileformat=VCFv4.1", "##source=ClinVar"])
        report_cv = self._build()["clinical_annotation_sources"]["report_clinvar"]
        self.assertIsNone(report_cv["version"])
        self.assertIn("unavailable", report_cv["version_basis"])

    # ── VEP cache observed from the cache itself ─────────────────────────
    def test_vep_cache_version_and_assembly_are_observed(self):
        prov = self._build()
        observed = prov["annotation_databases"]["vep_cache_observed"]
        self.assertEqual(observed["assembly"], "GRCh38")
        self.assertEqual(observed["version_dir"], "113_GRCh38")
        self.assertEqual(observed["info"]["source_gencode"], "GENCODE 47")
        self.assertEqual(len(observed["info_txt_sha256"]), 64)
        self.assertEqual(prov["annotation_databases"]["vep_cache_version"], "113")  # configured, kept
        self.assertEqual(len(prov["checksums"]["vep_cache_info_txt"]["sha256"]), 64)

    # ── Reference / resources ────────────────────────────────────────────
    def test_reference_facts_are_observed_configured_values_kept(self):
        prov = self._build()
        ref = prov["reference_genome"]
        self.assertEqual(ref["build"], "GRCh38")
        self.assertEqual(ref["observed"]["contig_count"], 1)
        self.assertTrue(ref["observed"]["fasta_fai_dict_consistent"])
        self.assertTrue(ref["observed"]["bwa_index"]["consistent_with_fasta"])
        self.assertEqual(len(ref["observed"]["fai_sha256"]), 64)
        self.assertEqual(len(prov["checksums"]["reference_fai"]["sha256"]), 64)
        self.assertEqual(len(prov["checksums"]["reference_dict"]["sha256"]), 64)

    def test_dbsnp_observed_from_header(self):
        observed = self._build()["reference_resources"]["dbsnp"]["observed"]
        self.assertEqual(observed["vcf_header"]["dbsnp_build_id"], "146")
        self.assertEqual(observed["vcf_header"]["reference"], "GRCh38.p2")

    def test_large_files_are_not_hashed_but_size_and_mtime_are_recorded(self):
        import resource_inspect
        facts = resource_inspect.file_facts(str(self.tmp / "clinvar.vcf.gz"), sha256_max_bytes=10)
        self.assertIsNone(facts["sha256"])
        self.assertIn("size >", facts["sha256_skipped_reason"])
        self.assertGreater(facts["size_bytes"], 10)
        self.assertIsNotNone(facts["mtime_utc"])

    def test_small_and_medium_key_files_have_sha256(self):
        checksums = self._build()["checksums"]
        for name in ("panel_bed", "tmb_coding_bed", "reference_fai", "reference_dict",
                     "clinvar_vcf", "clinvar_vcf_index", "vep_cache_info_txt"):
            self.assertEqual(len(checksums[name]["sha256"]), 64, name)

    # ── Observed software ────────────────────────────────────────────────
    def test_observed_software_versions_are_merged_from_probe_files(self):
        software = self._build(tool_version_files=_all_tool_files(self.tmp))["software_observed"]
        tools = software["tools"]
        self.assertEqual(tools["vep"]["version"], "113.0")
        self.assertEqual(tools["bcftools"]["version"], "1.19")
        self.assertEqual(tools["bcftools"]["htslib_version"], "1.19.1")
        self.assertEqual(tools["samtools"]["htslib_version"], "1.21")
        self.assertEqual(tools["bwa"]["version"], "0.7.17-r1188")
        self.assertEqual(tools["multiqc"]["version"], "1.34")
        self.assertEqual(tools["vep"]["cache_compatibility_check"]["status"], "PASS")
        self.assertEqual(tools["vep"]["conda_prefix"], "/envs/annotation")
        self.assertEqual(software["status"], "complete")

    def test_incomplete_probe_set_is_flagged_partial(self):
        software = self._build(tool_version_files=_all_tool_files(self.tmp)[:1])["software_observed"]
        self.assertTrue(software["status"].startswith("partial: missing"))
        self.assertIn("bwa", software["status"])

    def test_no_probe_files_means_not_collected(self):
        self.assertEqual(self._build()["software_observed"]["status"], "not_collected")

    def test_unreadable_probe_file_is_reported_not_fatal(self):
        bad = self.tmp / "bad.json"
        bad.write_text("{not json", encoding="utf-8")
        software = self._build(tool_version_files=[str(bad)])["software_observed"]
        self.assertEqual(software["status"], "not_collected")
        self.assertTrue(software["problems"])

    # ── Pipeline identity ────────────────────────────────────────────────
    def test_no_git_and_no_stamp_is_explicitly_unavailable(self):
        pipeline = self._build()["pipeline"]
        self.assertIsNone(pipeline["git_commit_sha"])
        self.assertEqual(pipeline["git_commit_source"], "unavailable")

    def test_build_stamp_supplies_the_commit_for_a_deployed_copy_without_git(self):
        (self.repo / "workflow").mkdir()
        (self.repo / "workflow" / "build_info.json").write_text(json.dumps({
            "git_commit_sha": "a" * 40, "git_describe": "v1-3-gaaaaaaa-dirty", "git_dirty": True,
            "stamped_at_utc": "2026-09-19T00:00:00+00:00"}), encoding="utf-8")
        pipeline = self._build()["pipeline"]
        self.assertEqual(pipeline["git_commit_sha"], "a" * 40)
        self.assertEqual(pipeline["git_commit_source"], "build_stamp")
        self.assertTrue(pipeline["git_dirty"])
        self.assertEqual(pipeline["version"], "v1-3-gaaaaaaa-dirty")

    @unittest.skipUnless(GIT_AVAILABLE, "git not available")
    def test_live_git_is_used_and_labelled_when_available(self):
        prov = provenance.build_provenance("S1", self.config)  # the real repository
        self.assertEqual(prov["pipeline"]["git_commit_source"], "git")
        self.assertEqual(len(prov["pipeline"]["git_commit_sha"]), 40)

    @unittest.skipUnless(GIT_AVAILABLE, "git not available")
    def test_write_build_stamp_round_trips(self):
        out = self.tmp / "stamp.json"
        stamp = provenance.write_build_stamp(provenance.REPO_ROOT, str(out))
        self.assertEqual(json.loads(out.read_text(encoding="utf-8"))["git_commit_sha"], stamp["git_commit_sha"])
        self.assertEqual(len(stamp["git_commit_sha"]), 40)

    def test_code_fingerprint_tracks_deployed_workflow_content(self):
        self.assertIsNone(provenance.compute_code_fingerprint(str(self.repo)))  # no workflow/ dir
        scripts = self.repo / "workflow" / "scripts"
        scripts.mkdir(parents=True)
        (self.repo / "workflow" / "Snakefile_x").write_text("rule all: pass\n", encoding="utf-8")
        (scripts / "a.py").write_text("x = 1\n", encoding="utf-8")
        first = provenance.compute_code_fingerprint(str(self.repo))
        self.assertEqual(len(first), 64)
        self.assertEqual(first, provenance.compute_code_fingerprint(str(self.repo)))  # deterministic
        (scripts / "a.py").write_text("x = 2\n", encoding="utf-8")
        self.assertNotEqual(first, provenance.compute_code_fingerprint(str(self.repo)))
        before_backup = provenance.compute_code_fingerprint(str(self.repo))
        (self.repo / "workflow" / "_old_backup").mkdir()
        (self.repo / "workflow" / "_old_backup" / "junk.py").write_text("y", encoding="utf-8")
        self.assertEqual(provenance.compute_code_fingerprint(str(self.repo)), before_backup)  # backups ignored

    # ── Robustness ───────────────────────────────────────────────────────
    def test_generated_at_is_preserved_and_record_serialisable(self):
        prov = self._build(tool_version_files=_all_tool_files(self.tmp))
        self.assertIn("generated_at_utc", prov)
        json.dumps(prov)

    def test_garbage_config_values_never_raise(self):
        config = {"ref": {"genome": 123, "dbsnp": ["x"], "manifest": 5},
                  "annotation": {"vep": {"cache_dir": 7, "extra": None}, "clinvar": {"vcf": {"a": 1}}},
                  "tmb": {"coding_bed": 0}}
        prov = provenance.build_provenance("S1", config, repo_root=str(self.repo))
        json.dumps(prov)
        self.assertEqual(prov["software_observed"]["status"], "not_collected")


@unittest.skipUnless(GIT_AVAILABLE, "git not available")
class TestGitDirtyState(unittest.TestCase):
    """git_dirty must reflect reality: true for uncommitted changes, false only for a clean tree."""

    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)
        self.repo = Path(self.tmpdir.name) / "repo"
        (self.repo / "workflow" / "scripts").mkdir(parents=True)
        (self.repo / "workflow" / "scripts" / "a.py").write_text("x = 1\n", encoding="utf-8")
        self._git("init", "-q")
        self._git("add", "-A")
        self._git("-c", "user.name=t", "-c", "user.email=t@example.invalid", "commit", "-q", "-m", "init")

    def _git(self, *args):
        subprocess.run(["git", *args], cwd=str(self.repo), check=True, capture_output=True)

    def _config(self):
        return {"ref": {}, "annotation": {"vep": {}}}

    def test_clean_checkout_reports_git_dirty_false(self):
        pipeline = provenance.build_provenance("S1", self._config(), repo_root=str(self.repo))["pipeline"]
        self.assertIs(pipeline["git_dirty"], False)
        self.assertEqual(pipeline["git_commit_source"], "git")
        self.assertEqual(len(pipeline["git_commit_sha"]), 40)
        self.assertFalse(pipeline["version"].endswith("-dirty"))

    def test_uncommitted_change_reports_git_dirty_true(self):
        (self.repo / "workflow" / "scripts" / "a.py").write_text("x = 2\n", encoding="utf-8")
        pipeline = provenance.build_provenance("S1", self._config(), repo_root=str(self.repo))["pipeline"]
        self.assertIs(pipeline["git_dirty"], True)
        self.assertTrue(pipeline["version"].endswith("-dirty"))

    def test_untracked_file_also_counts_as_dirty(self):
        (self.repo / "workflow" / "scripts" / "new.py").write_text("y = 1\n", encoding="utf-8")
        self.assertIs(provenance.get_git_dirty(str(self.repo)), True)

    def test_build_stamp_of_a_clean_checkout_yields_false_in_a_no_git_deployment(self):
        # deploy step: stamp the clean checkout, copy workflow/ (without .git) elsewhere
        stamp_path = self.repo / "workflow" / "build_info.json"
        stamp = provenance.write_build_stamp(str(self.repo), str(stamp_path))
        self.assertIs(stamp["git_dirty"], False)
        deployed = Path(self.tmpdir.name) / "deployed"
        shutil.copytree(self.repo / "workflow", deployed / "workflow")
        self.assertFalse((deployed / ".git").exists())
        pipeline = provenance.build_provenance("S1", self._config(), repo_root=str(deployed))["pipeline"]
        self.assertEqual(pipeline["git_commit_source"], "build_stamp")
        self.assertIs(pipeline["git_dirty"], False)
        self.assertEqual(pipeline["git_commit_sha"], stamp["git_commit_sha"])

    def test_build_stamp_of_a_dirty_checkout_is_not_laundered_to_clean(self):
        (self.repo / "workflow" / "scripts" / "a.py").write_text("x = 3\n", encoding="utf-8")
        stamp = provenance.write_build_stamp(str(self.repo), str(self.repo / "workflow" / "build_info.json"))
        self.assertIs(stamp["git_dirty"], True)
        deployed = Path(self.tmpdir.name) / "deployed"
        shutil.copytree(self.repo / "workflow", deployed / "workflow")
        pipeline = provenance.build_provenance("S1", self._config(), repo_root=str(deployed))["pipeline"]
        self.assertIs(pipeline["git_dirty"], True)


class TestReportUsesCustomClinVarSource(unittest.TestCase):
    """Guards the claim recorded in provenance: the report's ClinVar column is the custom VEP annotation."""

    def _read(self, name):
        return (SCRIPTS_DIR / name).read_text(encoding="utf-8")

    def test_somatic_table_reads_the_custom_clinvar_fields(self):
        text = self._read("vcf_to_somatic_table.py")
        self.assertIn('anno.get("ClinVar_CLNSIG"', text)
        self.assertIn('anno.get("ClinVar_CLNDN"', text)

    def test_neither_report_step_reads_the_cache_clin_sig(self):
        for name in ("vcf_to_somatic_table.py", "generate_comprehensive_report.py"):
            self.assertNotIn("CLIN_SIG", self._read(name), name)

    def test_vep_rule_still_annotates_with_the_custom_clinvar_vcf(self):
        rule = (SCRIPTS_DIR.parent / "rules" / "annotation_comprehensive.smk").read_text(encoding="utf-8")
        self.assertIn("--custom {params.clinvar_vcf},ClinVar,vcf,exact,0,CLNSIG,CLNDN", rule)


if __name__ == "__main__":
    unittest.main()
