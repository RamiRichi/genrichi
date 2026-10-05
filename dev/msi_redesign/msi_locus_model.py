"""
GenRichi MSI Methodology Redesign -- DEVELOPMENT PROTOTYPE. NOT PRODUCTION CODE.

Status: design-phase prototype, not wired into any Snakemake rule, not
imported by workflow/scripts/calculate_msi.py, not deployed, not clinically
validated in any way. See dev/msi_redesign/README.md for the full design
document and its rationale.

This module implements ONLY the per-locus classification / sample-level
aggregation CORE of the proposed locus-based MSI methodology -- i.e. "given
the observed tumor and normal repeat-length distributions at each marker
locus, classify stability and aggregate to a sample-level call." It
deliberately does NOT yet implement:
  - the reference-genome marker-scan step (deriving the marker panel itself
    from hg38.fa + a panel BED), or
  - the BAM/CRAM read-parsing step that produces real repeat-length
    distributions from aligned reads.
Those are the next implementation steps, once this classification core has
been reviewed. Every number below (MIN_LOCUS_DEPTH, INSTABILITY_DISTANCE_CUTOFF,
MIN_SHIFT_UNITS, MIN_INFORMATIVE_LOCI, MSI_H_FRACTION_THRESHOLD) is a
provisional starting parameter grounded in general panel-based MSI
literature convention -- NONE of it has been calibrated against a real
positive or negative control yet. Do not treat any classification produced
by this module as clinically meaningful.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import Sequence

# ── Provisional parameters (NOT calibrated; see module docstring) ───────────
MIN_LOCUS_DEPTH = 20                # min informative reads required, per sample, per locus
INSTABILITY_DISTANCE_CUTOFF = 0.20  # total-variation-distance cutoff for "different distribution"
MIN_SHIFT_UNITS = 1                 # minimum modal-length shift (repeat units) to call UNSTABLE
MIN_INFORMATIVE_LOCI = 5            # minimum callable loci required for any sample-level call
MSI_H_FRACTION_THRESHOLD = 0.40     # fraction of unstable/callable loci -> MSI-H

# ── Locus-level result states ────────────────────────────────────────────────
STABLE = "STABLE"
UNSTABLE = "UNSTABLE"
NOT_CALLABLE = "NOT_CALLABLE"

# ── Sample-level result states (the redesigned 4-way safety guard, see design doc §5) ─
INSUFFICIENT_DATA = "INSUFFICIENT_DATA"   # zero loci evaluated at all
INDETERMINATE = "INDETERMINATE"           # evaluated, but too few callable loci
VALID_MSS = "VALID_MSS"                   # enough callable loci, fraction below threshold
VALID_MSI_H = "VALID_MSI_H"               # enough callable loci, fraction at/above threshold


@dataclass
class LocusObservation:
    """Repeat lengths observed in reads at one marker locus, one sample."""
    marker_id: str
    tumor_lengths: Sequence[int]
    normal_lengths: Sequence[int]


@dataclass
class LocusClassification:
    marker_id: str
    status: str                 # STABLE / UNSTABLE / NOT_CALLABLE
    tumor_depth: int
    normal_depth: int
    tumor_mode: "int | None" = None
    normal_mode: "int | None" = None
    shift_units: "int | None" = None
    distance: "float | None" = None


@dataclass
class SampleClassification:
    sample_status: str          # one of the 4 sample-level states above
    n_loci_total: int = 0
    n_callable: int = 0
    n_unstable: int = 0
    n_not_callable: int = 0
    fraction_unstable: "float | None" = None
    loci: list = field(default_factory=list)   # list[LocusClassification]


def _mode(lengths: Sequence[int]) -> "int | None":
    if not lengths:
        return None
    counts = Counter(lengths)
    max_count = max(counts.values())
    # Counter.most_common() resolves ties by first-seen order, which would
    # make a locus call depend on samtools' read ordering. Use the smaller
    # repeat length as a deterministic tie-break; this remains provisional.
    return min(length for length, count in counts.items() if count == max_count)


def _histogram(lengths: Sequence[int]) -> dict:
    if not lengths:
        return {}
    total = len(lengths)
    counts = Counter(lengths)
    return {length: count / total for length, count in counts.items()}


def _total_variation_distance(hist_a: dict, hist_b: dict) -> float:
    """Symmetric distance in [0, 1] between two normalized length histograms.
    Dependency-free (no scipy) -- 0.5 * sum(|p(x) - q(x)|) over all lengths."""
    all_lengths = set(hist_a) | set(hist_b)
    return 0.5 * sum(abs(hist_a.get(length, 0.0) - hist_b.get(length, 0.0)) for length in all_lengths)


def classify_locus(
    observation: LocusObservation,
    min_depth: int = MIN_LOCUS_DEPTH,
    distance_cutoff: float = INSTABILITY_DISTANCE_CUTOFF,
    min_shift_units: int = MIN_SHIFT_UNITS,
) -> LocusClassification:
    """Classify one marker locus as STABLE / UNSTABLE / NOT_CALLABLE from its
    observed tumor and normal repeat-length distributions.

    A locus is NOT_CALLABLE if either sample has fewer than `min_depth`
    informative reads -- it is then excluded from the sample-level score
    entirely (never silently counted as "stable")."""
    t_depth = len(observation.tumor_lengths)
    n_depth = len(observation.normal_lengths)
    if t_depth < min_depth or n_depth < min_depth:
        return LocusClassification(
            marker_id=observation.marker_id, status=NOT_CALLABLE,
            tumor_depth=t_depth, normal_depth=n_depth,
        )

    t_mode = _mode(observation.tumor_lengths)
    n_mode = _mode(observation.normal_lengths)
    distance = _total_variation_distance(
        _histogram(observation.tumor_lengths), _histogram(observation.normal_lengths)
    )
    shift = abs(t_mode - n_mode)

    unstable = shift >= min_shift_units and distance >= distance_cutoff
    return LocusClassification(
        marker_id=observation.marker_id,
        status=UNSTABLE if unstable else STABLE,
        tumor_depth=t_depth, normal_depth=n_depth,
        tumor_mode=t_mode, normal_mode=n_mode,
        shift_units=shift, distance=round(distance, 4),
    )


def classify_sample(
    observations: Sequence[LocusObservation],
    min_depth: int = MIN_LOCUS_DEPTH,
    distance_cutoff: float = INSTABILITY_DISTANCE_CUTOFF,
    min_shift_units: int = MIN_SHIFT_UNITS,
    min_informative_loci: int = MIN_INFORMATIVE_LOCI,
    msi_h_fraction_threshold: float = MSI_H_FRACTION_THRESHOLD,
) -> SampleClassification:
    """Aggregate per-locus classifications into one of the 4 sample-level
    states. Never lets "too little evidence" render as a real MSS/MSI-H call."""
    if not observations:
        return SampleClassification(sample_status=INSUFFICIENT_DATA)

    loci = [
        classify_locus(obs, min_depth=min_depth, distance_cutoff=distance_cutoff,
                        min_shift_units=min_shift_units)
        for obs in observations
    ]
    n_total = len(loci)
    callable_loci = [l for l in loci if l.status != NOT_CALLABLE]
    n_callable = len(callable_loci)
    n_not_callable = n_total - n_callable
    n_unstable = sum(1 for l in callable_loci if l.status == UNSTABLE)

    if n_callable < min_informative_loci:
        return SampleClassification(
            sample_status=INDETERMINATE, n_loci_total=n_total, n_callable=n_callable,
            n_unstable=n_unstable, n_not_callable=n_not_callable, loci=loci,
        )

    fraction_unstable = n_unstable / n_callable
    status = VALID_MSI_H if fraction_unstable >= msi_h_fraction_threshold else VALID_MSS
    return SampleClassification(
        sample_status=status, n_loci_total=n_total, n_callable=n_callable,
        n_unstable=n_unstable, n_not_callable=n_not_callable,
        fraction_unstable=round(fraction_unstable, 4), loci=loci,
    )
