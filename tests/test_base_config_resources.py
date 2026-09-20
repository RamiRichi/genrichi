"""
GenRichi Phase 5.5 Step 12-B -- base config.yaml resource-reference regression test.

config/config.yaml's panel.bed must resolve to a file that actually exists
and is tracked in the repository. Step 12-B fixed a stale reference
(resources/panel/hotspots.bed, which never existed) to the tracked
resources/panel/example_hotspots.bed. Also guards that README.md's Quick
Start still carries the historical/unvalidated-pathway disclaimer added in
the same fix, so the two do not silently drift apart again.

Dependency-free by design: PyYAML is not installed in the conda environment
this suite normally runs under (confirmed in Step 12-B), so panel.bed is
extracted with a small regex -- config/config.yaml has exactly one `bed:`
key, so this is unambiguous -- rather than a full YAML parse. Same
regex-over-text style already used by tests/test_environment_ownership.py
for .smk parsing.
"""

import re
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
CONFIG_YAML = REPO_ROOT / "config" / "config.yaml"
README = REPO_ROOT / "README.md"

PANEL_BED_RE = re.compile(r"^\s*bed:\s*(\S+)", re.MULTILINE)


class TestBaseConfigPanelBed(unittest.TestCase):
    def test_panel_bed_points_to_an_existing_tracked_file(self):
        text = CONFIG_YAML.read_text(encoding="utf-8")
        match = PANEL_BED_RE.search(text)
        self.assertIsNotNone(match, "no 'bed:' line found in config/config.yaml")
        bed_path = REPO_ROOT / match.group(1)
        self.assertTrue(
            bed_path.is_file(),
            f"config/config.yaml panel.bed points to a nonexistent file: {match.group(1)}",
        )


class TestReadmeQuickStartDisclaimer(unittest.TestCase):
    def test_quick_start_carries_historical_pathway_disclaimer(self):
        text = README.read_text(encoding="utf-8")
        start = text.find("## Quick Start")
        self.assertNotEqual(start, -1, "README.md has no '## Quick Start' section")
        end = text.find("\n## ", start + 1)
        section = text[start : end if end != -1 else None]
        for phrase in ("historical", "validated Phase 1 workflow", "HCC1395"):
            with self.subTest(phrase=phrase):
                self.assertIn(
                    phrase,
                    section,
                    f"Quick Start section is missing expected disclaimer text: {phrase!r}",
                )


if __name__ == "__main__":
    unittest.main()
