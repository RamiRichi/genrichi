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

Phase 5.5 Step 11 additionally guards the deployment pin-propagation fix:
root envs/ must mirror the 5 active pin files (deploy_comprehensive.sh and
deploy_hereditary.sh symlink workflow/envs -> envs, so anything missing from
envs/ is invisible to a deployed Snakemake run), and that mirrored layout
must actually be discoverable by Snakemake's native pin-file lookup when
reached through that exact symlink shape.
"""

import re
import shutil
import tempfile
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


class TestRootEnvsMirrorsPinFiles(unittest.TestCase):
    """Guards the Step 11 fix: deploy_*.sh symlink workflow/envs -> envs, so
    the pin files must also live in root envs/, not just workflow/envs/."""

    def test_root_envs_pin_files_match_workflow_envs(self):
        for name in ACTIVE_ENVS:
            with self.subTest(env=name):
                src = WORKFLOW_ENVS / f"{name}.linux-64.pin.txt"
                dst = ROOT_ENVS / f"{name}.linux-64.pin.txt"
                self.assertTrue(dst.is_file(), f"envs/{name}.linux-64.pin.txt is missing")
                self.assertEqual(
                    dst.read_bytes(),
                    src.read_bytes(),
                    f"envs/{name}.linux-64.pin.txt differs from workflow/envs/{name}.linux-64.pin.txt",
                )

    def test_root_envs_no_retired_environments(self):
        for stem in RETIRED_ENVS:
            for suffix in (".yaml", ".linux-64.pin.txt"):
                with self.subTest(file=f"{stem}{suffix}"):
                    self.assertFalse((ROOT_ENVS / f"{stem}{suffix}").exists())


class TestDeployedSymlinkLayoutExposesPinFiles(unittest.TestCase):
    """
    Reproduces the exact mechanism deploy_comprehensive.sh / deploy_hereditary.sh
    use (a directory-level `workflow/envs -> envs` symlink) inside a throwaway
    tmpdir -- same tempfile-based deployed-copy pattern already used by
    tests/test_provenance.py -- and confirms Snakemake's own pin-file discovery
    algorithm (Env._get_aux_file: replace the .yaml/.yml suffix with
    .<platform>.pin.txt and check the resulting path exists -- see
    snakemake/deployment/conda.py, read and verified against the installed
    Snakemake 9.21.0 source in Phase 5.5 Steps 3-4) would find the pin file
    through that symlinked layout. Does not import snakemake itself -- no new
    test dependency -- it replicates the documented algorithm exactly rather
    than guessing at it.
    """

    @staticmethod
    def _pin_path_for(yaml_path: Path) -> Path:
        name = yaml_path.name
        if name.endswith(".yaml"):
            stem = name[: -len(".yaml")]
        elif name.endswith(".yml"):
            stem = name[: -len(".yml")]
        else:
            raise ValueError(f"not a .yaml/.yml file: {yaml_path}")
        return yaml_path.with_name(f"{stem}.linux-64.pin.txt")

    def test_pin_files_discoverable_through_deploy_style_symlink(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            root_envs = tmp_path / "envs"
            root_envs.mkdir()
            for name in ACTIVE_ENVS:
                shutil.copy2(ROOT_ENVS / f"{name}.yaml", root_envs / f"{name}.yaml")
                shutil.copy2(
                    ROOT_ENVS / f"{name}.linux-64.pin.txt",
                    root_envs / f"{name}.linux-64.pin.txt",
                )

            workflow_dir = tmp_path / "workflow"
            workflow_dir.mkdir()
            symlinked_envs = workflow_dir / "envs"
            try:
                symlinked_envs.symlink_to(root_envs, target_is_directory=True)
            except OSError as exc:
                self.skipTest(f"symlinks not supported on this platform/permissions: {exc}")

            for name in ACTIVE_ENVS:
                with self.subTest(env=name):
                    yaml_path = symlinked_envs / f"{name}.yaml"
                    self.assertTrue(yaml_path.exists())
                    pin_path = self._pin_path_for(yaml_path)
                    self.assertEqual(pin_path.name, f"{name}.linux-64.pin.txt")
                    self.assertTrue(
                        pin_path.exists(),
                        f"pin file for {name} not discoverable through the deploy-style symlink",
                    )


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
