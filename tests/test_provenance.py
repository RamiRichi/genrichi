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

import provenance  # noqa: E402


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
        self.assertEqual(record["schema_version"], "1.0")
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


if __name__ == "__main__":
    unittest.main()
