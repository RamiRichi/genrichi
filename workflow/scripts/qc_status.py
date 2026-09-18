"""
GenRichi Phase 5.2 — explicit VALID vs INSUFFICIENT_DATA status semantics for
MSI and CNV, used by calculate_msi.py and calculate_cnv.py.

Pure functions, no Snakemake dependency, so they're independently
unit-testable. This module does not change any scientific calculation,
threshold, or comparison — it only distinguishes "evidence was evaluated and
came back negative" from "there was no evidence to evaluate at all", so that
a report can never present the second case as if it were the first.
"""

from __future__ import annotations

VALID = "VALID"
INSUFFICIENT_DATA = "INSUFFICIENT_DATA"


def msi_qc_status(total_variants: int) -> str:
    """
    VALID              -> at least one PASS variant was evaluated (a
                           zero-indel result is still scientifically
                           meaningful: it means no MS-context indels were
                           found among variants that WERE assessed).
    INSUFFICIENT_DATA   -> zero PASS variants were evaluated at all; there is
                           no evidence to support any MSI/MSS call.
    """
    return VALID if total_variants and total_variants > 0 else INSUFFICIENT_DATA


def cnv_qc_status(n_overlapping_regions: int) -> str:
    """
    VALID              -> tumor/normal per-region coverage had at least one
                           overlapping region to compare (a call of zero
                           significant segments is still scientifically
                           meaningful).
    INSUFFICIENT_DATA   -> zero overlapping regions between tumor and normal;
                           there is no evidence to support any CNV call.
    """
    return VALID if n_overlapping_regions and n_overlapping_regions > 0 else INSUFFICIENT_DATA
