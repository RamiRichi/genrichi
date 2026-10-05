"""
SYNTHETIC, DEVELOPMENT-ONLY tests for dev/msi_redesign/msi_locus_model.py.

*** These tests use fabricated, non-biological repeat-length data. ***
*** They prove the classification LOGIC behaves correctly on deliberately ***
*** constructed examples. They are NOT biological validation, NOT a claim ***
*** of clinical accuracy, and NOT evidence that the (not-yet-built) real ***
*** BAM-parsing/marker-scan layers work against real sequencing data. ***

This module is not part of the production test suite's pipeline-facing
checks (it imports nothing from workflow/scripts and exercises no
Snakemake rule) -- it exists solely to review the redesigned classification
core before any production integration is considered.

Run with:
    python -m unittest tests.test_msi_locus_model_dev -v
"""

import sys
import unittest
from pathlib import Path

DEV_DIR = Path(__file__).resolve().parents[1] / "dev" / "msi_redesign"
sys.path.insert(0, str(DEV_DIR))

import msi_locus_model as m  # noqa: E402


def synthetic_observation(marker_id, tumor_lengths, normal_lengths):
    """Builds a fabricated LocusObservation. Values are invented for this
    test only and do not correspond to any real sample."""
    return m.LocusObservation(marker_id=marker_id, tumor_lengths=tumor_lengths, normal_lengths=normal_lengths)


class TestLocusClassification(unittest.TestCase):
    """Per-locus classification on deliberately constructed, synthetic examples."""

    def test_identical_distributions_are_stable(self):
        # Tumor and normal both show a clean length-20 peak -- no instability.
        obs = synthetic_observation("SYN-1", tumor_lengths=[20] * 30, normal_lengths=[20] * 30)
        result = m.classify_locus(obs)
        self.assertEqual(result.status, m.STABLE)
        self.assertEqual(result.tumor_mode, 20)
        self.assertEqual(result.normal_mode, 20)
        self.assertEqual(result.shift_units, 0)

    def test_clear_length_shift_is_unstable(self):
        # Normal is a clean length-20 peak; tumor has shifted to a length-18
        # dominant peak with a residual length-20 subpopulation -- the
        # textbook synthetic signature of microsatellite instability.
        obs = synthetic_observation(
            "SYN-2",
            tumor_lengths=[18] * 22 + [20] * 8,
            normal_lengths=[20] * 30,
        )
        result = m.classify_locus(obs)
        self.assertEqual(result.status, m.UNSTABLE)
        self.assertEqual(result.tumor_mode, 18)
        self.assertEqual(result.normal_mode, 20)
        self.assertEqual(result.shift_units, 2)
        self.assertGreaterEqual(result.distance, m.INSTABILITY_DISTANCE_CUTOFF)

    def test_small_shift_below_distance_cutoff_stays_stable(self):
        # Tumor mode shifts by 1 unit, but only a small minority of reads --
        # distance stays below the cutoff, so this must NOT be over-called.
        obs = synthetic_observation(
            "SYN-3",
            tumor_lengths=[20] * 27 + [19] * 3,
            normal_lengths=[20] * 30,
        )
        result = m.classify_locus(obs)
        self.assertEqual(result.status, m.STABLE)
        self.assertLess(result.distance, m.INSTABILITY_DISTANCE_CUTOFF)

    def test_low_depth_in_tumor_is_not_callable(self):
        obs = synthetic_observation("SYN-4", tumor_lengths=[20] * 5, normal_lengths=[20] * 30)
        result = m.classify_locus(obs)
        self.assertEqual(result.status, m.NOT_CALLABLE)

    def test_low_depth_in_normal_is_not_callable(self):
        obs = synthetic_observation("SYN-5", tumor_lengths=[20] * 30, normal_lengths=[20] * 4)
        result = m.classify_locus(obs)
        self.assertEqual(result.status, m.NOT_CALLABLE)

    def test_not_callable_locus_reports_depths_but_no_mode(self):
        obs = synthetic_observation("SYN-6", tumor_lengths=[20] * 2, normal_lengths=[20] * 2)
        result = m.classify_locus(obs)
        self.assertEqual(result.status, m.NOT_CALLABLE)
        self.assertEqual(result.tumor_depth, 2)
        self.assertEqual(result.normal_depth, 2)
        self.assertIsNone(result.tumor_mode)


class TestSampleAggregation(unittest.TestCase):
    """Sample-level aggregation on deliberately constructed synthetic panels."""

    def test_zero_loci_is_insufficient_data(self):
        result = m.classify_sample([])
        self.assertEqual(result.sample_status, m.INSUFFICIENT_DATA)
        self.assertEqual(result.n_loci_total, 0)

    def test_too_few_callable_loci_is_indeterminate(self):
        # 3 loci total, only 2 callable (below MIN_INFORMATIVE_LOCI=5) --
        # evidence exists but is too thin for any MSS/MSI-H call.
        obs = [
            synthetic_observation("L1", [20] * 30, [20] * 30),   # callable, stable
            synthetic_observation("L2", [20] * 30, [20] * 30),   # callable, stable
            synthetic_observation("L3", [20] * 2, [20] * 2),     # not callable
        ]
        result = m.classify_sample(obs)
        self.assertEqual(result.sample_status, m.INDETERMINATE)
        self.assertEqual(result.n_callable, 2)
        self.assertIsNone(result.fraction_unstable)

    def test_all_stable_synthetic_panel_is_valid_mss(self):
        # 10 synthetic loci, all clearly stable, all well above min depth --
        # fabricated MSS-like panel.
        obs = [synthetic_observation(f"L{i}", [20] * 30, [20] * 30) for i in range(10)]
        result = m.classify_sample(obs)
        self.assertEqual(result.sample_status, m.VALID_MSS)
        self.assertEqual(result.n_callable, 10)
        self.assertEqual(result.n_unstable, 0)
        self.assertEqual(result.fraction_unstable, 0.0)

    def test_majority_unstable_synthetic_panel_is_valid_msi_h(self):
        # 10 synthetic loci: 7 fabricated as clearly unstable (shifted modal
        # peak, majority-shifted reads), 3 stable -- fraction_unstable=0.70,
        # comfortably above MSI_H_FRACTION_THRESHOLD=0.40.
        unstable = [
            synthetic_observation(f"U{i}", tumor_lengths=[18] * 22 + [20] * 8, normal_lengths=[20] * 30)
            for i in range(7)
        ]
        stable = [
            synthetic_observation(f"S{i}", tumor_lengths=[20] * 30, normal_lengths=[20] * 30)
            for i in range(3)
        ]
        result = m.classify_sample(unstable + stable)
        self.assertEqual(result.sample_status, m.VALID_MSI_H)
        self.assertEqual(result.n_callable, 10)
        self.assertEqual(result.n_unstable, 7)
        self.assertEqual(result.fraction_unstable, 0.70)

    def test_fraction_just_below_threshold_is_valid_mss_not_msi_h(self):
        # 10 loci, exactly 3 unstable -> fraction=0.30, below the 0.40
        # threshold. Confirms the boundary is not accidentally inclusive.
        unstable = [
            synthetic_observation(f"U{i}", tumor_lengths=[18] * 22 + [20] * 8, normal_lengths=[20] * 30)
            for i in range(3)
        ]
        stable = [
            synthetic_observation(f"S{i}", tumor_lengths=[20] * 30, normal_lengths=[20] * 30)
            for i in range(7)
        ]
        result = m.classify_sample(unstable + stable)
        self.assertEqual(result.sample_status, m.VALID_MSS)
        self.assertEqual(result.fraction_unstable, 0.30)

    def test_not_callable_loci_are_excluded_not_counted_as_stable(self):
        # 8 loci total: 5 stable+callable, 3 fabricated as NOT_CALLABLE
        # (low depth). The 3 must be excluded from the denominator entirely,
        # not silently folded into "stable".
        stable = [synthetic_observation(f"S{i}", [20] * 30, [20] * 30) for i in range(5)]
        uncallable = [synthetic_observation(f"NC{i}", [20] * 3, [20] * 3) for i in range(3)]
        result = m.classify_sample(stable + uncallable)
        self.assertEqual(result.n_loci_total, 8)
        self.assertEqual(result.n_callable, 5)
        self.assertEqual(result.n_not_callable, 3)
        self.assertEqual(result.sample_status, m.VALID_MSS)  # 5 callable >= MIN_INFORMATIVE_LOCI

    def test_four_states_are_all_distinct_values(self):
        # Guards the core safety principle from the design doc: no two of
        # the four sample-level states may ever collapse into each other.
        states = {m.INSUFFICIENT_DATA, m.INDETERMINATE, m.VALID_MSS, m.VALID_MSI_H}
        self.assertEqual(len(states), 4)


if __name__ == "__main__":
    unittest.main()
