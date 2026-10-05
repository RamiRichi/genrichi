"""
SYNTHETIC, DEVELOPMENT-ONLY tests for dev/msi_redesign/build_filtered_panel.py.

*** Fabricated marker/annotation data only. The real filtered panel and ***
*** its manifest come from a separate real run against the real Tier 2 ***
*** universe and the real downloaded reference tracks, with its own ***
*** provenance reported alongside it.

Run with:
    python -m unittest tests.test_msi_build_filtered_panel_dev -v
"""

import csv
import sys
import tempfile
import unittest
from pathlib import Path

DEV_DIR = Path(__file__).resolve().parents[1] / "dev" / "msi_redesign"
sys.path.insert(0, str(DEV_DIR))

import build_filtered_panel as bfp  # noqa: E402


class TestClassifySelfRmsk(unittest.TestCase):
    def test_no_hit(self):
        self.assertEqual(bfp.classify_self_rmsk(set()), "no_rmsk_hit")

    def test_microsatellite_only(self):
        self.assertEqual(bfp.classify_self_rmsk({"Simple_repeat"}), "microsatellite_class_only")

    def test_mixed(self):
        self.assertEqual(bfp.classify_self_rmsk({"Simple_repeat", "SINE"}), "mixed_microsatellite_and_other")

    def test_unexpected_other(self):
        self.assertEqual(bfp.classify_self_rmsk({"SINE", "LINE"}), "unexpected_non_microsatellite_class")


class TestDetermineExclusionReasons(unittest.TestCase):
    def test_clean_locus_is_kept(self):
        reasons = bfp.determine_exclusion_reasons(
            "m1", {"m1": {"Simple_repeat"}}, {}, {}, {"m1": 1.0}, {"m1": (1.0, 1.0, 1.0)},
        )
        self.assertEqual(reasons, [])

    def test_umap_zero_singleread_excludes(self):
        reasons = bfp.determine_exclusion_reasons(
            "m1", {}, {}, {}, {"m1": 0.0}, {"m1": (0.5, 0.5, 0.5)},
        )
        self.assertIn("UMAP_UNMAPPABLE", reasons)

    def test_umap_zero_min_multiread_excludes(self):
        reasons = bfp.determine_exclusion_reasons(
            "m1", {}, {}, {}, {"m1": 0.8}, {"m1": (0.5, 0.0, 1.0)},
        )
        self.assertIn("UMAP_UNMAPPABLE", reasons)

    def test_nonzero_scores_do_not_trigger_umap_exclusion(self):
        reasons = bfp.determine_exclusion_reasons(
            "m1", {}, {}, {}, {"m1": 0.1}, {"m1": (0.1, 0.01, 0.5)},
        )
        self.assertNotIn("UMAP_UNMAPPABLE", reasons)

    def test_segdup_overlap_excludes(self):
        reasons = bfp.determine_exclusion_reasons(
            "m1", {}, {}, {"m1": [0.95]}, {"m1": 1.0}, {"m1": (1.0, 1.0, 1.0)},
        )
        self.assertIn("SEGDUP_OVERLAP", reasons)

    def test_rmsk_self_non_microsatellite_excludes(self):
        reasons = bfp.determine_exclusion_reasons(
            "m1", {"m1": {"SINE"}}, {}, {}, {"m1": 1.0}, {"m1": (1.0, 1.0, 1.0)},
        )
        self.assertIn("REPEATMASKER_SELF_NON_MICROSATELLITE", reasons)

    def test_mixed_self_class_is_not_excluded_by_rmsk_rule(self):
        reasons = bfp.determine_exclusion_reasons(
            "m1", {"m1": {"SINE", "Simple_repeat"}}, {}, {}, {"m1": 1.0}, {"m1": (1.0, 1.0, 1.0)},
        )
        self.assertNotIn("REPEATMASKER_SELF_NON_MICROSATELLITE", reasons)

    def test_flank_repeat_class_never_causes_exclusion(self):
        reasons = bfp.determine_exclusion_reasons(
            "m1", {"m1": {"Simple_repeat"}}, {"m1": {"SINE", "LINE"}}, {}, {"m1": 1.0}, {"m1": (1.0, 1.0, 1.0)},
        )
        self.assertEqual(reasons, [])

    def test_multiple_reasons_all_recorded(self):
        reasons = bfp.determine_exclusion_reasons(
            "m1", {"m1": {"LINE"}}, {}, {"m1": [0.9]}, {"m1": 0.0}, {"m1": (0.0, 0.0, 0.0)},
        )
        self.assertEqual(
            set(reasons),
            {"UMAP_UNMAPPABLE", "SEGDUP_OVERLAP", "REPEATMASKER_SELF_NON_MICROSATELLITE"},
        )

    def test_missing_umap_data_does_not_crash_or_falsely_exclude(self):
        reasons = bfp.determine_exclusion_reasons("m1", {}, {}, {}, {}, {})
        self.assertEqual(reasons, [])


class TestBuildEndToEnd(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.annotations_dir = Path(self.tmp.name) / "annotations"
        self.annotations_dir.mkdir()
        self.out_dir = Path(self.tmp.name) / "output"

    def _write(self, name, lines):
        path = self.annotations_dir / name
        path.write_text("\n".join(lines) + "\n" if lines else "")
        return str(path)

    def test_raw_file_is_copied_unchanged_and_filtered_is_a_subset(self):
        raw_path = Path(self.tmp.name) / "tier2_raw_input.tsv"
        raw_path.write_text(
            "marker_id\tchrom\tstart\tend\tmotif\tref_repeat_length\treference_build\tkind\n"
            "m1\tchr1\t100\t110\tA\t10\tGRCh38\tmononucleotide\n"
            "m2\tchr1\t500\t510\tA\t10\tGRCh38\tmononucleotide\n"
        )
        self._write("locus_self_rmsk_classes.tsv", ["m1\tSimple_repeat", "m2\tSINE"])
        self._write("flank_rmsk_classes.tsv", [])
        self._write("window_segdup_hits.tsv", [])
        self._write("window_umap_singleread_coverage.tsv", [
            "chr1\t0\t210\tm1\t1\t210\t210\t1.0000000",
            "chr1\t400\t610\tm2\t1\t210\t210\t1.0000000",
        ])
        self._write("window_umap_multiread_scores.tsv", [
            "chr1\t0\t210\tm1\t1.0\t1.0\t1.0",
            "chr1\t400\t610\tm2\t1.0\t1.0\t1.0",
        ])

        import build_filtered_panel as mod
        tier2_raw_immutable = str(self.out_dir / "tier2_raw.v1.tsv")
        self.out_dir.mkdir()
        with open(raw_path, "rb") as src, open(tier2_raw_immutable, "wb") as dst:
            dst.write(src.read())
        annotated_path, filtered_path, cumulative, exclusion_counts = mod.build(
            tier2_raw_immutable, str(self.annotations_dir), str(self.out_dir)
        )

        # raw file byte-identical to input
        self.assertEqual(Path(raw_path).read_bytes(), Path(tier2_raw_immutable).read_bytes())

        # m1 kept (clean), m2 excluded (SINE self-class)
        with open(filtered_path, encoding="utf-8") as fh:
            filtered_ids = [row["marker_id"] for row in csv.DictReader(fh, delimiter="\t")]
        self.assertEqual(filtered_ids, ["m1"])
        self.assertEqual(cumulative["raw_candidate_universe"], 2)
        self.assertEqual(cumulative["after_excluding_umap_unmappable_and_segdup_and_rmsk_self"], 1)
        self.assertEqual(exclusion_counts.get("REPEATMASKER_SELF_NON_MICROSATELLITE"), 1)
        self.assertEqual(exclusion_counts.get("KEPT"), 1)


if __name__ == "__main__":
    unittest.main()
