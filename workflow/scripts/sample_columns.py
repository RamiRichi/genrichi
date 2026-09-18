"""
GenRichi Phase 5.3 -- tumor/normal VCF sample-column resolution.

Pure, dependency-light helper used by vcf_to_somatic_table.py to identify
which VCF sample column is tumor and which is normal, replacing a
previously unvalidated assumption that the tumor sample is always the last
column.

GenRichi's paired BAMs are read-grouped with SM tags "{sample}_tumor" and
"{sample}_normal" (see workflow/rules/align_paired.smk), and Mutect2 is
invoked with --normal-sample "{sample}_normal" (see
workflow/rules/somatic_calling_paired.smk) -- these two SM tags are exactly
what Mutect2 preserves as the VCF's sample column names. This module
resolves tumor vs. normal by NAME (not position), and raises a clear,
actionable error rather than silently returning a swapped column when the
expected relationship cannot be established.

No Snakemake dependency, so this is independently unit-testable.
"""

from __future__ import annotations


class SampleColumnError(ValueError):
    """Raised when tumor/normal VCF sample columns cannot be safely resolved."""


def resolve_tumor_normal_columns(col_names, expected_sample_id=None):
    """
    Given a VCF #CHROM header's column names (the full list -- the 8 fixed
    VCF columns, FORMAT, and the sample columns), return
    (tumor_col, normal_col): the exact sample-column names to use.

    Resolution is by name (each BAM's read-group SM tag is preserved by
    Mutect2 as the VCF sample name), never by position, and this function
    raises SampleColumnError -- rather than silently guessing -- whenever
    the expected tumor/normal relationship cannot be confirmed:
      * not exactly 2 sample columns
      * not exactly one column ending in "_tumor" and one ending in "_normal"
      * (when expected_sample_id is given) the resolved names don't match
        "{expected_sample_id}_tumor" / "{expected_sample_id}_normal" exactly
    """
    sample_cols = list(col_names[9:])
    if len(sample_cols) != 2:
        raise SampleColumnError(
            f"Expected exactly 2 VCF sample columns (tumor + normal) after "
            f"the FORMAT column, found {len(sample_cols)}: {sample_cols!r}. "
            f"Refusing to guess which is tumor -- this would risk silently "
            f"swapping VAF/DP/AD."
        )

    tumor_matches = [c for c in sample_cols if c.endswith("_tumor")]
    normal_matches = [c for c in sample_cols if c.endswith("_normal")]

    if len(tumor_matches) != 1 or len(normal_matches) != 1:
        raise SampleColumnError(
            f"Could not identify tumor/normal VCF sample columns by name "
            f"(expected exactly one column ending '_tumor' and one ending "
            f"'_normal'; found columns: {sample_cols!r}). Refusing to fall "
            f"back to positional selection -- this would risk silently "
            f"swapping VAF/DP/AD."
        )

    tumor_col, normal_col = tumor_matches[0], normal_matches[0]

    if expected_sample_id is not None:
        expected_tumor = f"{expected_sample_id}_tumor"
        expected_normal = f"{expected_sample_id}_normal"
        if tumor_col != expected_tumor or normal_col != expected_normal:
            raise SampleColumnError(
                f"VCF sample columns do not match the expected sample_id "
                f"'{expected_sample_id}': expected tumor column "
                f"'{expected_tumor}' and normal column '{expected_normal}', "
                f"found tumor='{tumor_col}' normal='{normal_col}'."
            )

    return tumor_col, normal_col
