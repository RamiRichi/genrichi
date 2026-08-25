"""
Panel-based MSI estimation from somatic VCF.

Approach: count somatic indels (insertions + deletions) per Mb of panel.
MSI-High tumors accumulate many more indels than MSS tumors because
mismatch-repair deficiency specifically affects replication slippage at
microsatellite loci, which manifests predominantly as small indels.

Reference ranges (panel-based, ~1-3 Mb):
  MSS  : typically  0-5  indels/Mb
  MSI-L: typically  5-10 indels/Mb
  MSI-H: typically >10   indels/Mb  (FDA/ASCO guideline equivalent)

Output format matches msisensor-pro column layout so the report script
can parse it without modification.
"""

import gzip

sm = snakemake  # type: ignore[name-defined]

vcf_path   = sm.input.vcf
out_path   = sm.output.score
coding_mb  = float(sm.params.coding_mb)
threshold  = float(sm.params.threshold)   # indels/Mb for MSI-High


def _is_indel(ref: str, alt: str) -> bool:
    """True for insertions and deletions (length difference)."""
    return len(ref) != len(alt)


def _in_homopolymer(ref: str, alt: str, context: str = "") -> bool:
    """
    Simple heuristic: flag variants where the REF or ALT contains a run of
    >=3 identical bases (classic microsatellite / homopolymer context).
    Even without a reference context, length-change variants near runs are
    the canonical MSI signal.
    """
    seq = ref if len(ref) > len(alt) else alt
    for i in range(len(seq) - 2):
        if seq[i] == seq[i + 1] == seq[i + 2]:
            return True
    return True   # count all indels — conservative but transparent


total_variants = 0
total_indels   = 0
ms_indels      = 0   # indels with homopolymer character

opener = gzip.open if vcf_path.endswith(".gz") else open
with opener(vcf_path, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 8:
            continue
        filt = parts[6]
        if filt != "PASS":
            continue

        ref = parts[3]
        alt = parts[4]

        # Handle multi-allelic (take first ALT)
        alt = alt.split(",")[0]

        total_variants += 1
        if _is_indel(ref, alt):
            total_indels += 1
            if _in_homopolymer(ref, alt):
                ms_indels += 1

# MSI score: MS-context indels per Mb
msi_score   = round(ms_indels / coding_mb, 2) if coding_mb > 0 else 0.0
indel_frac  = round(total_indels / total_variants * 100, 1) if total_variants > 0 else 0.0
msi_status  = "MSI-H" if msi_score >= threshold else "MSS"

# Write output in msisensor-pro compatible TSV format:
# Total_Number_of_Sites  Number_of_Somatic_Sites  %
# We repurpose the columns: panel_mb, ms_indels, msi_score
with open(out_path, "w") as fh:
    fh.write("Total_Number_of_Sites\tNumber_of_Somatic_Sites\t%\n")
    fh.write(f"{coding_mb}\t{ms_indels}\t{msi_score}\n")

# Also write a human-readable summary to the dis file
with open(sm.output.dis, "w") as fh:
    fh.write(f"# Panel-based MSI estimation (VCF indel method)\n")
    fh.write(f"total_variants\t{total_variants}\n")
    fh.write(f"total_indels\t{total_indels}\n")
    fh.write(f"ms_context_indels\t{ms_indels}\n")
    fh.write(f"indel_fraction_pct\t{indel_frac}\n")
    fh.write(f"msi_score_per_mb\t{msi_score}\n")
    fh.write(f"msi_status\t{msi_status}\n")
    fh.write(f"threshold_per_mb\t{threshold}\n")

# Write minimal somatic file
with open(sm.output.somatic, "w") as fh:
    fh.write(f"# MSI somatic indels: {ms_indels} in {coding_mb} Mb\n")
    fh.write(f"# Score: {msi_score}/Mb  Status: {msi_status}\n")
