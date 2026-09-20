"""
GenRichi Phase 5.5 -- conda environment ownership regression tests.

Guards the outcome of the Phase 5.5 environment audit/retirement:
  - the 5 active envs (align/calling/annotation/qc/report) remain present
    and carry their native Snakemake pin files;
  - germline.yaml, cnv.yaml, msi.yaml stay retired (do not silently
    reappear);
  - the hereditary/germline pathway keeps resolving to calling.yaml, not a
    reintroduced germline.yaml;
  - no `conda:` directive anywhere in the workflow ever points at a
    non-active env file (structural guard against this whole class of
    drift recurring).

Read-only: inspects file existence and rule-file text only, never builds or
solves a conda environment.
"""

import re
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
WORKFLOW_ENVS = REPO_ROOT / "workflow" / "envs"
ROOT_ENVS = REPO_ROOT / "envs"
RULES_DIR = REPO_ROOT / "workflow" / "rules"

ACTIVE_ENVS = ["align", "calling", "annotation", "qc", "report"]
RETIRED_ENVS = ["germline", "cnv", "msi"]

CONDA_DIRECTIVE_RE = re.compile(r'conda:\s*"\.\./envs/([\w.-]+\.yaml)"')


def _smk_and_snakefile_paths():
    paths = list(RULES_DIR.glob("*.smk"))
    paths += [p for p in (REPO_ROOT / "workflow").glob("Snakefile*") if p.is_file()]
    return paths


class TestFiveActiveEnvironmentsRemainRequired(unittest.TestCase):
    def test_active_env_yamls_exist(self):
        for name in ACTIVE_ENVS:
            with self.subTest(env=name):
                self.assertTrue(
                    (WORKFLOW_ENVS / f"{name}.yaml").is_file(),
                    f"workflow/envs/{name}.yaml is missing",
                )

    def test_active_env_pin_files_exist(self):
        for name in ACTIVE_ENVS:
            with self.subTest(env=name):
                self.assertTrue(
                    (WORKFLOW_ENVS / f"{name}.linux-64.pin.txt").is_file(),
                    f"workflow/envs/{name}.linux-64.pin.txt is missing",
                )


class TestRetiredEnvironmentsStayRetired(unittest.TestCase):
    def test_germline_yaml_retired(self):
        self.assertFalse((WORKFLOW_ENVS / "germline.yaml").exists())
        self.assertFalse((ROOT_ENVS / "germline.yaml").exists())

    def test_cnv_yaml_retired(self):
        self.assertFalse((WORKFLOW_ENVS / "cnv.yaml").exists())

    def test_msi_yaml_retired(self):
        self.assertFalse((WORKFLOW_ENVS / "msi.yaml").exists())


class TestGermlinePathwayResolvesToCallingYaml(unittest.TestCase):
    def test_germline_calling_rules_use_calling_yaml(self):
        text = (RULES_DIR / "germline_calling.smk").read_text(encoding="utf-8")
        directives = CONDA_DIRECTIVE_RE.findall(text)
        self.assertTrue(directives, "no conda: directives found in germline_calling.smk")
        for env_file in directives:
            self.assertEqual(env_file, "calling.yaml")


class TestNoRuleSelectsANonActiveEnvironment(unittest.TestCase):
    def test_every_conda_directive_resolves_to_an_active_env(self):
        active_yaml_names = {f"{name}.yaml" for name in ACTIVE_ENVS}
        offenders = []
        for path in _smk_and_snakefile_paths():
            text = path.read_text(encoding="utf-8")
            for env_file in CONDA_DIRECTIVE_RE.findall(text):
                if env_file not in active_yaml_names:
                    offenders.append(f"{path.relative_to(REPO_ROOT)} -> {env_file}")
        self.assertEqual(
            offenders,
            [],
            "found conda: directive(s) pointing at a non-active env: " + ", ".join(offenders),
        )


if __name__ == "__main__":
    unittest.main()
