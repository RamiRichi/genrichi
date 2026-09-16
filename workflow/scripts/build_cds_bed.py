#!/usr/bin/env python3
"""
Build panel CDS BED by intersecting MANE Select v1.4 GTF CDS features
with the gene-body panel BED.

Usage:
    python build_cds_bed.py --mane MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz \
                            --panel resources/panel/comprehensive_genes.bed \
                            --out resources/panel/comprehensive_genes_cds.bed

Output: BED4 file of merged, non-overlapping CDS intervals clipped to panel gene bodies.
Stdlib only — no bedtools, no pandas, no external dependencies.
"""

import argparse
import gzip
import re
import sys
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(description="Panel CDS BED builder from MANE Select GTF")
    p.add_argument("--mane", required=True, help="MANE Select GTF (.gtf or .gtf.gz)")
    p.add_argument("--panel", required=True, help="Panel gene-body BED (BED3 or BED4)")
    p.add_argument("--out", required=True, help="Output CDS BED path")
    return p.parse_args()


def normalize_chrom(chrom: str) -> str:
    """Normalize Ensembl/NCBI chromosome names to UCSC-style chr-prefixed."""
    if chrom.startswith("chr"):
        return chrom
    if chrom in ("MT", "M"):
        return "chrM"
    return f"chr{chrom}"


def read_panel_bed(path: str):
    """Return list of (chrom, start, end, gene). Coordinates are 0-based half-open."""
    intervals = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom = normalize_chrom(parts[0])
            start = int(parts[1])
            end = int(parts[2])
            gene = parts[3].strip() if len(parts) > 3 else ""
            intervals.append((chrom, start, end, gene))
    return intervals


def build_panel_index(panel_intervals):
    """Index panel intervals by chromosome."""
    idx: dict = defaultdict(list)
    for chrom, start, end, gene in panel_intervals:
        idx[chrom].append((start, end, gene))
    return idx


def parse_gtf_attribute(attrs: str, key: str) -> str:
    m = re.search(rf'{key}\s+"([^"]+)"', attrs)
    return m.group(1) if m else ""


def read_mane_cds(gtf_path: str, panel_idx: dict):
    """
    Parse MANE GTF and return CDS records clipped to panel intervals.
    GTF coordinates are 1-based inclusive; converted to 0-based half-open here.
    Only CDS features from protein_coding transcripts are kept.
    """
    opener = gzip.open if gtf_path.endswith(".gz") else open
    records = []  # (chrom, start, end, gene)
    n_raw = 0

    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != "CDS":
                continue

            chrom = normalize_chrom(parts[0])
            if chrom not in panel_idx:
                continue  # chromosome not in panel at all

            # GTF is 1-based inclusive → 0-based half-open
            cds_start = int(parts[3]) - 1
            cds_end = int(parts[4])
            attrs = parts[8]

            transcript_type = parse_gtf_attribute(attrs, "transcript_type")
            if transcript_type and transcript_type != "protein_coding":
                continue

            gene = parse_gtf_attribute(attrs, "gene_name")
            if not gene:
                gene = parse_gtf_attribute(attrs, "gene_id")

            n_raw += 1

            # Clip to each overlapping panel interval
            for pstart, pend, panel_gene in panel_idx[chrom]:
                if cds_start >= pend or cds_end <= pstart:
                    continue  # no overlap
                clipped_start = max(cds_start, pstart)
                clipped_end = min(cds_end, pend)
                records.append((chrom, clipped_start, clipped_end, panel_gene or gene))

    print(f"  Raw protein-coding CDS lines processed: {n_raw}", file=sys.stderr)
    return records


def merge_intervals(records):
    """
    Sort and merge overlapping/adjacent intervals per chromosome.
    Returns sorted list of (chrom, start, end, gene).
    """
    by_chrom: dict = defaultdict(list)
    for chrom, start, end, gene in records:
        by_chrom[chrom].append((start, end, gene))

    # Chromosome sort order: chr1-22 numerically, then chrX, chrY, chrM
    def chrom_key(c):
        c = c.replace("chr", "")
        if c.isdigit():
            return (0, int(c))
        return (1, c)

    merged = []
    for chrom in sorted(by_chrom, key=chrom_key):
        segs = sorted(by_chrom[chrom], key=lambda x: x[0])
        cur_s, cur_e, cur_g = segs[0]
        for s, e, g in segs[1:]:
            if s <= cur_e:
                cur_e = max(cur_e, e)
            else:
                merged.append((chrom, cur_s, cur_e, cur_g))
                cur_s, cur_e, cur_g = s, e, g
        merged.append((chrom, cur_s, cur_e, cur_g))

    return merged


def main():
    args = parse_args()

    print("Reading panel BED...", file=sys.stderr)
    panel = read_panel_bed(args.panel)
    panel_idx = build_panel_index(panel)
    panel_chroms = set(panel_idx.keys())
    print(f"  {len(panel)} intervals on {len(panel_chroms)} chromosomes: {sorted(panel_chroms)}", file=sys.stderr)

    print("Reading MANE Select GTF CDS features...", file=sys.stderr)
    cds_records = read_mane_cds(args.mane, panel_idx)
    print(f"  CDS records overlapping panel (before merge): {len(cds_records)}", file=sys.stderr)

    print("Merging overlapping intervals...", file=sys.stderr)
    merged = merge_intervals(cds_records)

    total_bp = sum(e - s for _, s, e, _ in merged)
    genes = {g for _, _, _, g in merged}
    total_mb = total_bp / 1_000_000

    print(f"Writing {args.out}...", file=sys.stderr)
    with open(args.out, "w") as fh:
        for chrom, start, end, gene in merged:
            fh.write(f"{chrom}\t{start}\t{end}\t{gene}\n")

    # Machine-readable summary to stdout
    print(f"CDS_BP={total_bp}")
    print(f"CDS_MB={total_mb:.6f}")
    print(f"N_INTERVALS={len(merged)}")
    print(f"N_GENES={len(genes)}")

    # Human-readable summary to stderr
    print(f"\n=== CDS BED Summary ===", file=sys.stderr)
    print(f"  Intervals (merged, non-overlapping): {len(merged)}", file=sys.stderr)
    print(f"  Genes represented: {len(genes)}", file=sys.stderr)
    print(f"  Total coding territory: {total_bp:,} bp  ({total_mb:.6f} Mb)", file=sys.stderr)
    print(f"  Output: {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
