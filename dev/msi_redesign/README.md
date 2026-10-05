# MSI Methodology Redesign — Development Phase

**Status: design + prototype classification core only. Not integrated into
the production pipeline. Not clinically validated. No production file has
been modified to create this directory.**

## Why this exists

The production MSI implementation (`workflow/scripts/calculate_msi.py`) is a
coarse proxy: it counts already-called somatic indels whose REF/ALT allele
string happens to contain a run of 3+ identical bases, normalized by a fixed
1.5 Mb borrowed from the TMB denominator. It has no defined microsatellite
marker set, never independently measures repeat lengths, and never compares
tumor against normal at the sequence-repeat level. It has never been
validated against a real MSI-H sample (none currently exists anywhere in
this repository or production), and — separately from data availability —
it is not the right kind of algorithm to validate against one even if it
did.

This directory holds the redesign: a locus-based methodology modeled on
established tumor/normal microsatellite analysis principles (the same
family of approach used by MSIsensor-pro / mSINGS / MANTIS), built without
reintroducing the `msisensor-pro` binary this project has previously had to
avoid due to htslib/libdeflate conda conflicts.

## What's here

- `msi_locus_model.py` — the per-locus classification and sample-level
  aggregation core. Pure functions, no I/O, no Snakemake dependency,
  independently unit-testable. This is the part reviewed and tested in this
  development phase.
- `tests/test_msi_locus_model_dev.py` (repo `tests/` dir, `_dev` suffix) —
  synthetic unit tests. Every input is fabricated for the test; none of it
  is biological data or a claim of clinical accuracy.

## What's deliberately NOT here yet

1. **The marker-scan step** — deriving the actual microsatellite marker
   panel (genomic coordinates, motif, reference repeat length) from
   `hg38.fa` intersected with a panel BED. The design commits to this being
   *generated deterministically and content-hashed*, not hand-copied from a
   remembered coordinate list, specifically so the marker panel itself is
   auditable and regenerable rather than a source of unverifiable claims.
2. **The BAM/CRAM read-parsing step** — turning real aligned reads into the
   per-locus `LocusObservation` (tumor/normal repeat-length lists) that
   `msi_locus_model.py` consumes. Two candidate approaches were evaluated
   (see the full design conversation): a `pysam`-based CIGAR-aware reader
   (not currently an installed dependency), or a zero-new-dependency
   `samtools mpileup`-based text parser (more fragile to implement
   correctly). Recommendation: attempt the `samtools`-based route first
   given this project's history of conda environment fragility; fall back
   to adding `pysam` only if that proves impractical.
3. **Any production integration** — no change to `msi_scoring.smk`,
   `calculate_msi.py`, `comprehensive_config.yaml`, or
   `generate_comprehensive_report.py`. Integration is a separate, later
   decision, made only after this core is reviewed and after the two
   missing pieces above are built and independently tested.
4. **Any real validation run.** Nothing in this directory has been run
   against HCC1395 or any other real sample. The validation matrix (MSS via
   HCC1395 / independently documented; MSI-H via a real, independently
   characterized positive control that does not currently exist in this
   environment) is defined in the design conversation but intentionally not
   executed here.

## Provisional parameters (not calibrated)

`MIN_LOCUS_DEPTH=20`, `INSTABILITY_DISTANCE_CUTOFF=0.20`,
`MIN_SHIFT_UNITS=1`, `MIN_INFORMATIVE_LOCI=5`,
`MSI_H_FRACTION_THRESHOLD=0.40` — starting values grounded in general
panel-based MSI literature convention, not reverse-engineered from any
known sample's expected result. All five need empirical recalibration once
real MSS and MSI-H controls are available (see the design document's
validation-strategy section).
