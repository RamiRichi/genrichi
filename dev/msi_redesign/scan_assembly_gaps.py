"""
GenRichi MSI Methodology Redesign -- DEVELOPMENT PROTOTYPE. NOT PRODUCTION CODE.

Scans the real production hg38.fa (primary assembly only, same scope as
scan_markers_v2.py) for runs of 'N' (assembly gaps / unresolved bases) and
writes a BED of their positions. Used read-only, reference-genome-only, to
check whether a candidate MSI locus sits near an assembly gap -- a locus
right next to a gap has physically insufficient flanking reference
sequence for any read to anchor uniquely, regardless of any sample's real
sequencing data. This uses no sample/BAM data and no new external
download; it is a pure re-scan of the reference already on disk.

Usage:
    python3 scan_assembly_gaps.py --reference /path/to/hg38.fa \
        --out dev/msi_redesign/reference_annotations/assembly_gaps.hg38.bed
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys

PRIMARY_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
_N_RUN_RE = re.compile(r"N+")


def extract_chromosome_sequence(samtools: str, reference: str, chrom: str) -> str:
    out = subprocess.run([samtools, "faidx", reference, chrom], capture_output=True, text=True, timeout=600)
    if out.returncode != 0:
        raise RuntimeError(f"samtools faidx failed for {chrom}: {out.stderr.strip()}")
    lines = out.stdout.splitlines()
    return "".join(lines[1:]).upper()


def find_n_runs(seq: str):
    """Yields (start, end) 0-based half-open for each maximal run of 'N'."""
    for m in _N_RUN_RE.finditer(seq):
        yield m.start(), m.end()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--reference", required=True)
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--out", required=True)
    ap.add_argument("--chroms", nargs="*", default=None)
    args = ap.parse_args(argv)

    chroms = args.chroms if args.chroms else PRIMARY_CHROMS
    n_total = 0
    with open(args.out, "w", encoding="utf-8", newline="\n") as fh:
        for chrom in chroms:
            seq = extract_chromosome_sequence(args.samtools, args.reference, chrom)
            count = 0
            for start, end in find_n_runs(seq):
                fh.write(f"{chrom}\t{start}\t{end}\n")
                count += 1
            n_total += count
            print(f"  {chrom}: {count} gap runs", file=sys.stderr)
    print(f"wrote {n_total} gap runs -> {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
