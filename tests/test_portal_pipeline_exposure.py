"""
GenRichi Phase 5.5 Step 12-E -- portal pipeline-exposure contract test.

A Phase 5.5 audit (Steps 12-C/12-D) found that `workflow/Snakefile` ("hotspot")
has a deterministic DAG-breaking bug (a filename mismatch between align.smk's
mosdepth output and report.smk's expected input) plus an unfixed historical
mosdepth --quantize crash risk, and that `workflow/Snakefile_hereditary` has
no regression baseline or classifier test coverage. Both were removed from
portal/config.py's PIPELINE_MAP (and PANEL_PRICES) so the portal cannot sell
either as a customer-facing order -- without touching either pipeline's own
code. This guards that restriction against silently reverting.

Dependency-light by design: portal/config.py only does `import os` and
defines plain dict/set literals, so it is imported directly (no Flask
required) the same way tests/test_qc_status.py imports workflow/scripts
modules -- by inserting the containing directory onto sys.path.

Since Phase 5.5 Step 14, portal/config.py requires SECRET_KEY and
PORTAL_PASS to be present in the environment (it raises RuntimeError
otherwise -- see tests/test_portal_secrets.py for that behavior itself).
This file only exercises PIPELINE_MAP/PANEL_PRICES, so it supplies
obviously-fake, non-secret placeholder values via os.environ.setdefault
(never overriding a real value if one is already set) purely so the
import succeeds; SMTP is left disabled so no SMTP_PASS is needed either.
"""

import os
import sys
import unittest
from pathlib import Path

PORTAL_DIR = Path(__file__).resolve().parents[1] / "portal"
sys.path.insert(0, str(PORTAL_DIR))

os.environ.setdefault("SECRET_KEY", "unit-test-placeholder-not-a-real-secret")
os.environ.setdefault("PORTAL_PASS", "unit-test-placeholder-not-a-real-secret")
os.environ.setdefault("SMTP_ENABLED", "false")

import config as portal_config  # noqa: E402


class TestPortalPipelineExposure(unittest.TestCase):
    def test_comprehensive_remains_present(self):
        self.assertIn("comprehensive", portal_config.PIPELINE_MAP)
        self.assertIn("comprehensive", portal_config.PANEL_PRICES)

    def test_hotspot_is_absent(self):
        self.assertNotIn("hotspot", portal_config.PIPELINE_MAP)
        self.assertNotIn("hotspot", portal_config.PANEL_PRICES)

    def test_hereditary_is_absent(self):
        self.assertNotIn("hereditary", portal_config.PIPELINE_MAP)
        self.assertNotIn("hereditary", portal_config.PANEL_PRICES)

    def test_pipeline_map_has_exactly_one_entry(self):
        self.assertEqual(list(portal_config.PIPELINE_MAP.keys()), ["comprehensive"])
        self.assertEqual(list(portal_config.PANEL_PRICES.keys()), ["comprehensive"])

    def test_comprehensive_entry_is_unchanged(self):
        entry = portal_config.PIPELINE_MAP["comprehensive"]
        self.assertEqual(entry["snakefile"], "Snakefile_comprehensive")
        self.assertEqual(entry["configfile"], "config/comprehensive_config.yaml")
        self.assertEqual(entry["label"], "Somatic Comprehensive Panel")
        self.assertEqual(
            entry["description"],
            "Full somatic workup — SNV/indel + CNV + MSI + TMB (55 cancer genes)",
        )
        self.assertEqual(portal_config.PANEL_PRICES["comprehensive"], 750.00)

    def test_paired_panels_unchanged(self):
        self.assertEqual(portal_config.PAIRED_PANELS, {"comprehensive"})

    def test_new_order_default_panel_type_is_still_valid(self):
        app_py = (PORTAL_DIR / "app.py").read_text(encoding="utf-8")
        self.assertIn(
            'request.form.get("panel_type", "comprehensive")',
            app_py,
            "new_order's panel_type fallback must default to a key that still "
            "exists in PIPELINE_MAP, not a removed one",
        )


if __name__ == "__main__":
    unittest.main()
