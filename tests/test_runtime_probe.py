"""
Unit tests for workflow/scripts/runtime_probe.py (GenRichi Phase 5.3).

The parser fixtures below are real output captured from the validated
HCC1395 conda environments (VEP 113.0, bcftools 1.19/htslib 1.19.1, bwa
0.7.17, samtools 1.19.2 on htslib 1.21, MultiQC 1.34, fastp 0.23.4,
mosdepth 0.3.10, GATK 4.6.1.0).

Run with:
    python -m unittest tests.test_runtime_probe -v
"""

import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import runtime_probe  # noqa: E402

VEP_HELP = """#----------------------------------#
# ENSEMBL VARIANT EFFECT PREDICTOR #
#----------------------------------#

Versions:
  ensembl              : 113.58650ec
  ensembl-compara      : 113.9cf749d
  ensembl-funcgen      : 113.e30608c
  ensembl-io           : 113.bee6816
  ensembl-variation    : 113.ebfce74
  ensembl-vep          : 113.0

Help: dev@ensembl.org , helpdesk@ensembl.org
"""
BCFTOOLS = "bcftools 1.19\nUsing htslib 1.19.1\nCopyright (C) 2023 Genome Research Ltd.\n"
SAMTOOLS = "samtools 1.19.2\nUsing htslib 1.21\nCopyright (C) 2024 Genome Research Ltd.\n"
BWA = "\nProgram: bwa (alignment via Burrows-Wheeler transformation)\nVersion: 0.7.17-r1188\nContact: Heng Li\n"
GATK = ("Using GATK jar /x/gatk-package-4.6.1.0-local.jar\nRunning:\n    java -jar x --version\n"
        "The Genome Analysis Toolkit (GATK) v4.6.1.0\nHTSJDK Version: 4.1.3\nPicard Version: 3.3.0\n")


class TestParsers(unittest.TestCase):
    def test_vep_reports_the_executable_version_not_the_api_hash(self):
        parsed = runtime_probe.parse_vep(VEP_HELP)
        self.assertEqual(parsed["version"], "113.0")
        self.assertEqual(parsed["detail"]["ensembl"], "113.58650ec")

    def test_bcftools_and_its_htslib(self):
        self.assertEqual(runtime_probe.parse_htslib_tool("bcftools")(BCFTOOLS),
                         {"version": "1.19", "htslib": "1.19.1"})

    def test_samtools_htslib_can_differ_from_bcftools_htslib(self):
        self.assertEqual(runtime_probe.parse_htslib_tool("samtools")(SAMTOOLS),
                         {"version": "1.19.2", "htslib": "1.21"})

    def test_bwa_version_from_usage_text(self):
        self.assertEqual(runtime_probe.parse_bwa(BWA)["version"], "0.7.17-r1188")

    def test_multiqc_fastp_mosdepth_gatk(self):
        self.assertEqual(runtime_probe.parse_multiqc("multiqc, version 1.34\n")["version"], "1.34")
        self.assertEqual(runtime_probe.parse_fastp("fastp 0.23.4\n")["version"], "0.23.4")
        self.assertEqual(runtime_probe.parse_mosdepth("mosdepth 0.3.10\n")["version"], "0.3.10")
        self.assertEqual(runtime_probe.parse_gatk(GATK)["version"], "4.6.1.0")

    def test_unrecognised_output_gives_none_not_a_guess(self):
        self.assertIsNone(runtime_probe.parse_vep("garbage")["version"])
        self.assertIsNone(runtime_probe.parse_bwa("")["version"])
        self.assertIsNone(runtime_probe.parse_multiqc(None)["version"])


class TestProbeTool(unittest.TestCase):
    def test_missing_executable_is_reported_as_not_found(self):
        record = runtime_probe.probe_tool("vep", which=lambda _n: None)
        self.assertFalse(record["found"])
        self.assertIsNone(record["version"])

    def test_found_tool_version_comes_from_the_tool_output(self):
        record = runtime_probe.probe_tool(
            "bcftools", runner=lambda cmd: (0, BCFTOOLS), which=lambda n: "/env/bin/" + n)
        self.assertTrue(record["found"])
        self.assertEqual(record["path"], "/env/bin/bcftools")
        self.assertEqual(record["version"], "1.19")
        self.assertEqual(record["htslib"], "1.19.1")

    def test_nonzero_exit_is_normal_for_bwa_usage(self):
        record = runtime_probe.probe_tool("bwa", runner=lambda cmd: (1, BWA), which=lambda n: "/x/bwa")
        self.assertEqual(record["version"], "0.7.17-r1188")

    def test_unknown_tool_name_is_rejected(self):
        with self.assertRaises(KeyError):
            runtime_probe.probe_tool("not-a-tool")

    def test_run_command_never_raises_for_a_missing_binary(self):
        rc, text = runtime_probe.run_command(["definitely-not-a-real-binary-xyz"])
        self.assertIsNone(rc)
        self.assertTrue(text)

    def test_run_command_does_not_raise_on_non_utf8_output(self):
        """Reproduces the 2026-09-22 production failure (order GR-20260922-1E88F4):
        Debian/Ubuntu's reproducible-builds gcc writes a single raw 0xAB/0xBB byte
        (not the two-byte UTF-8 encoding of U+00AB/U+00BB) into samtools' compiled-in
        -ffile-prefix-map=... build-flags string, which subprocess.run(text=True)
        (strict UTF-8) cannot decode -- crashing observe_tools_alignment before
        alignment, variant calling, annotation, report generation or provenance are
        ever reached. run_command() must degrade gracefully instead."""
        raw = (b"samtools 1.19.2\nUsing htslib 1.21\n"
               b"Compiler flags: -ffile-prefix-map=\xabBUILDPATH\xbb=. -O2\n")
        cmd = [sys.executable, "-c",
               "import sys; sys.stdout.buffer.write(" + repr(raw) + ")"]
        rc, text = runtime_probe.run_command(cmd)
        self.assertEqual(rc, 0)
        self.assertIn("samtools 1.19.2", text)
        self.assertIn("�", text)  # the bad byte, replaced, not raised
        # the parser must still extract the version from the clean line above it
        self.assertEqual(runtime_probe.parse_htslib_tool("samtools")(text)["version"], "1.19.2")

    def test_probe_tool_survives_non_utf8_output_end_to_end(self):
        raw_samtools = (b"samtools 1.19.2\nUsing htslib 1.21\n"
                         b"Compiler flags: -ffile-prefix-map=\xabBUILDPATH\xbb=. -O2\n")
        record = runtime_probe.probe_tool(
            "samtools",
            runner=lambda cmd: runtime_probe.run_command(
                [sys.executable, "-c", "import sys; sys.stdout.buffer.write(" + repr(raw_samtools) + ")"]),
            which=lambda n: "/env/bin/" + n,
        )
        self.assertTrue(record["found"])
        self.assertEqual(record["version"], "1.19.2")
        self.assertEqual(record["htslib"], "1.21")

    @unittest.skipUnless(shutil.which("samtools"), "samtools not on PATH")
    def test_real_installed_samtools_version_is_probed_without_crashing(self):
        """The actual regression: run_command() against the real samtools binary
        installed on this machine, not a mock -- this is exactly the call
        observe_tools_alignment makes in production."""
        record = runtime_probe.probe_tool("samtools")
        self.assertTrue(record["found"])
        self.assertIsNotNone(record["version"], f"first_output_line={record.get('first_output_line')!r}")

    @unittest.skipUnless(shutil.which("bwa"), "bwa not on PATH")
    def test_real_installed_bwa_version_is_probed_without_crashing(self):
        record = runtime_probe.probe_tool("bwa")
        self.assertTrue(record["found"])
        self.assertIsNotNone(record["version"])

    @unittest.skipUnless(shutil.which("mosdepth"), "mosdepth not on PATH")
    def test_real_installed_mosdepth_version_is_probed_without_crashing(self):
        record = runtime_probe.probe_tool("mosdepth")
        self.assertTrue(record["found"])
        self.assertIsNotNone(record["version"])


class TestVepCacheCompatibility(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)
        self.cache = Path(self.tmpdir.name)
        (self.cache / "homo_sapiens" / "113_GRCh38").mkdir(parents=True)

    def test_matching_release_passes(self):
        result = runtime_probe.check_vep_matches_cache("113.0", str(self.cache), "113", "GRCh38")
        self.assertEqual(result["status"], "PASS")

    def test_executable_newer_than_cache_fails(self):
        result = runtime_probe.check_vep_matches_cache("116.2", str(self.cache), "113", "GRCh38")
        self.assertEqual(result["status"], "FAIL")
        self.assertIn("release 116", result["message"])
        self.assertIn("113", result["message"])

    def test_unpinned_uses_the_executable_release_and_needs_that_cache(self):
        self.assertEqual(
            runtime_probe.check_vep_matches_cache("113.4", str(self.cache), None, "GRCh38")["status"], "PASS")
        failing = runtime_probe.check_vep_matches_cache("116.2", str(self.cache), None, "GRCh38")
        self.assertEqual(failing["status"], "FAIL")
        self.assertIn("not found", failing["message"])

    def test_unreadable_vep_version_fails_rather_than_passing(self):
        self.assertEqual(runtime_probe.check_vep_matches_cache(None, str(self.cache), "113", "GRCh38")["status"], "FAIL")

    def test_no_cache_dir_is_skipped(self):
        self.assertEqual(runtime_probe.check_vep_matches_cache("113.0", None, "113", "GRCh38")["status"], "SKIPPED")

    def test_cache_version_parsed_from_vep_extra_args(self):
        self.assertEqual(runtime_probe.vep_cache_version_from_extra("--everything --cache_version 113"), "113")
        self.assertEqual(runtime_probe.vep_cache_version_from_extra("--cache_version=115 --x"), "115")
        self.assertIsNone(runtime_probe.vep_cache_version_from_extra("--everything"))
        self.assertIsNone(runtime_probe.vep_cache_version_from_extra(None))


class TestRecordAndCli(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)
        self.tmp = Path(self.tmpdir.name)

    def _fake_tools(self, vep_version="113.0"):
        outputs = {"vep": VEP_HELP.replace("113.0", vep_version), "bcftools": BCFTOOLS}
        return (lambda cmd: (0, outputs[cmd[0]])), (lambda n: "/env/bin/" + n)

    def test_record_contains_observed_versions_and_environment(self):
        runner, which = self._fake_tools()
        cache = self.tmp / "cache"
        (cache / "homo_sapiens" / "113_GRCh38").mkdir(parents=True)
        record = runtime_probe.build_record(
            "annotation", ["vep", "bcftools"],
            vep_cache={"cache_dir": str(cache), "cache_version": "113", "assembly": "GRCh38"},
            runner=runner, which=which)
        self.assertEqual(record["tools"]["vep"]["version"], "113.0")
        self.assertEqual(record["tools"]["bcftools"]["htslib"], "1.19.1")
        self.assertEqual(record["vep_cache_check"]["status"], "PASS")
        self.assertEqual(runtime_probe.record_problems(record), [])
        self.assertIn("vep 113.0", runtime_probe.summary_line(record))
        json.dumps(record)

    def test_vep_cache_mismatch_is_a_problem(self):
        runner, which = self._fake_tools(vep_version="116.2")
        cache = self.tmp / "cache"
        (cache / "homo_sapiens" / "113_GRCh38").mkdir(parents=True)
        record = runtime_probe.build_record(
            "annotation", ["vep", "bcftools"],
            vep_cache={"cache_dir": str(cache), "cache_version": "113", "assembly": "GRCh38"},
            runner=runner, which=which)
        problems = runtime_probe.record_problems(record)
        self.assertEqual(len(problems), 1)
        self.assertIn("release 116", problems[0])

    def test_missing_required_tool_is_a_problem(self):
        record = runtime_probe.build_record("annotation", ["vep"], which=lambda n: None)
        self.assertTrue(any("not found" in p for p in runtime_probe.record_problems(record)))

    def test_cli_writes_json_and_returns_zero_for_a_healthy_tool(self):
        out = self.tmp / "sub" / "versions.json"
        fake = {"mosdepth": ([sys.executable, "-c", "print('mosdepth 9.9.9')"], runtime_probe.parse_mosdepth)}
        with mock.patch.dict(runtime_probe.TOOLS, fake):
            rc = runtime_probe.main(["--label", "alignment", "--tools", "mosdepth", "--out", str(out)])
        self.assertEqual(rc, 0)
        record = json.loads(out.read_text(encoding="utf-8"))
        self.assertEqual(record["tools"]["mosdepth"]["version"], "9.9.9")
        self.assertEqual(record["label"], "alignment")

    def test_cli_returns_nonzero_but_still_writes_json_when_a_tool_is_missing(self):
        out = self.tmp / "versions.json"
        fake = {"mosdepth": (["definitely-not-a-real-binary-xyz"], runtime_probe.parse_mosdepth)}
        with mock.patch.dict(runtime_probe.TOOLS, fake):
            rc = runtime_probe.main(["--label", "alignment", "--tools", "mosdepth", "--out", str(out)])
        self.assertEqual(rc, 1)
        self.assertFalse(json.loads(out.read_text(encoding="utf-8"))["tools"]["mosdepth"]["found"])

    def test_cli_rejects_unknown_tools(self):
        self.assertEqual(runtime_probe.main(["--label", "x", "--tools", "nope", "--out", str(self.tmp / "o.json")]), 2)


if __name__ == "__main__":
    unittest.main()
