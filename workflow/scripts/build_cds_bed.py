#!/usr/bin/env python3
"""
Build a panel CDS BED from the MANE Select GTF and a panel BED.

For every gene named in the panel BED (column 4) the single transcript tagged
``MANE_Select`` is resolved from the GTF, its CDS exons are clipped to that
gene's own panel intervals, and overlapping/adjacent pieces are merged within
the gene (never across genes).

Usage:
    python build_cds_bed.py --mane MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz \
                            --panel resources/panel/phase1_solid_tumor/solid_tumor_phase1_v1.bed \
                            --out resources/panel/phase1_solid_tumor/solid_tumor_phase1_v1_cds.bed \
                            [--audit audit.json]

Output: BED4 (chrom, start, end, gene), 0-based half-open, sorted chr1-22, chrX,
chrY, chrM then by position. The script fails (non-zero exit, message on
stderr) if the panel has no gene names, a panel gene has no MANE_Select
transcript or more than one, a gene yields no CDS, or the output would contain
empty / out-of-panel / cross-gene-merged records.
Stdlib only.
"""

import argparse
import gzip
import json
import re
import sys
from collections import defaultdict

MANE_TAG = "MANE_Select"


class CdsBuildError(Exception):
    pass


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Panel CDS BED builder from MANE Select GTF")
    p.add_argument("--mane", required=True, help="MANE Select GTF (.gtf or .gtf.gz)")
    p.add_argument("--panel", required=True, help="Panel BED with gene names in column 4")
    p.add_argument("--out", required=True, help="Output CDS BED path")
    p.add_argument("--audit", help="Optional JSON audit path (per-gene transcript and bp)")
    return p.parse_args(argv)


def normalize_chrom(chrom: str) -> str:
    """Normalize Ensembl/NCBI chromosome names to UCSC-style chr-prefixed."""
    if chrom.startswith("chr"):
        return chrom
    if chrom in ("MT", "M"):
        return "chrM"
    return f"chr{chrom}"


def chrom_key(chrom: str):
    """chr1-22 numerically, then chrX, chrY, chrM. Anything else is an error."""
    c = chrom[3:] if chrom.startswith("chr") else chrom
    if c.isdigit() and 1 <= int(c) <= 22:
        return (0, int(c))
    order = {"X": 1, "Y": 2, "M": 3}
    if c in order:
        return (order[c], 0)
    raise CdsBuildError(f"unsupported chromosome '{chrom}' (expected chr1-22, chrX, chrY, chrM)")


def read_panel_bed(path: str):
    """Return {gene: [(chrom, start, end), ...]}; 0-based half-open."""
    panel = defaultdict(list)
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split("\t")
            if len(parts) < 4 or not parts[3].strip():
                raise CdsBuildError(
                    f"{path}:{lineno}: panel BED needs a gene name in column 4 (got {len(parts)} column(s))"
                )
            try:
                start, end = int(parts[1]), int(parts[2])
            except ValueError:
                raise CdsBuildError(f"{path}:{lineno}: non-integer coordinates")
            if start < 0 or end <= start:
                raise CdsBuildError(f"{path}:{lineno}: empty or invalid interval {start}-{end}")
            chrom = normalize_chrom(parts[0])
            chrom_key(chrom)
            panel[parts[3].strip()].append((chrom, start, end))
    if not panel:
        raise CdsBuildError(f"{path}: no panel intervals found")
    return panel


def gtf_attr(attrs: str, key: str) -> str:
    m = re.search(rf'(?:^|[\s;]){key}\s+"([^"]*)"', attrs)
    return m.group(1) if m else ""


def gtf_tags(attrs: str):
    return set(re.findall(r'(?:^|[\s;])tag\s+"([^"]*)"', attrs))


def open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def resolve_mane_transcripts(gtf_path: str, genes):
    """Return {gene: transcript_id} for exactly one MANE_Select transcript per gene.

    Uses transcript features only. Fails on a missing or non-unique transcript.
    """
    found = defaultdict(set)
    wanted = set(genes)
    with open_text(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "transcript":
                continue
            gene = gtf_attr(parts[8], "gene_name")
            if gene not in wanted or MANE_TAG not in gtf_tags(parts[8]):
                continue
            found[gene].add(gtf_attr(parts[8], "transcript_id"))
    missing = sorted(g for g in wanted if not found.get(g))
    multiple = {g: sorted(t) for g, t in found.items() if len(t) > 1}
    problems = []
    if missing:
        problems.append(f"no {MANE_TAG} transcript in GTF for: {', '.join(missing)}")
    if multiple:
        problems.append(
            f"more than one {MANE_TAG} transcript for: "
            + "; ".join(f"{g} ({', '.join(t)})" for g, t in sorted(multiple.items()))
        )
    if problems:
        raise CdsBuildError("; ".join(problems))
    return {g: next(iter(t)) for g, t in found.items()}


def read_transcript_cds(gtf_path: str, transcripts):
    """Return {gene: [(chrom, start, end)]} CDS exons of the chosen transcripts (0-based half-open)."""
    by_tx = {tx: g for g, tx in transcripts.items()}
    cds = defaultdict(list)
    with open_text(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            tx = gtf_attr(parts[8], "transcript_id")
            gene = by_tx.get(tx)
            if gene is None:
                continue
            # GTF is 1-based inclusive -> BED 0-based half-open
            cds[gene].append((normalize_chrom(parts[0]), int(parts[3]) - 1, int(parts[4])))
    return cds


def clip_to_panel(cds_exons, panel_regions):
    """Clip CDS exons to the same gene's panel intervals (same chromosome only)."""
    clipped = []
    for chrom, cs, ce in cds_exons:
        for pchrom, ps, pe in panel_regions:
            if chrom != pchrom or cs >= pe or ce <= ps:
                continue
            clipped.append((chrom, max(cs, ps), min(ce, pe)))
    return clipped


def merge_within_gene(segments):
    """Merge overlapping or adjacent (chrom, start, end) segments of ONE gene."""
    merged = []
    for chrom, s, e in sorted(segments, key=lambda x: (chrom_key(x[0]), x[1], x[2])):
        if merged and merged[-1][0] == chrom and s <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], e)
        else:
            merged.append([chrom, s, e])
    return [tuple(m) for m in merged]


def build_cds(panel, transcripts, cds_by_gene):
    """Return sorted [(chrom, start, end, gene)] and a per-gene audit dict."""
    records, audit = [], {}
    no_cds = []
    for gene in sorted(panel):
        pieces = merge_within_gene(clip_to_panel(cds_by_gene.get(gene, []), panel[gene]))
        if not pieces:
            no_cds.append(gene)
            continue
        records.extend((c, s, e, gene) for c, s, e in pieces)
        audit[gene] = {
            "transcript_id": transcripts[gene],
            "cds_intervals": len(pieces),
            "cds_bp": sum(e - s for _, s, e in pieces),
        }
    if no_cds:
        raise CdsBuildError(f"no CDS overlaps the panel intervals for: {', '.join(no_cds)}")
    records.sort(key=lambda r: (chrom_key(r[0]), r[1], r[2], r[3]))
    return records, audit


def validate_output(records, panel):
    """Every record non-empty, inside a same-gene panel interval, one gene per record."""
    for chrom, s, e, gene in records:
        if e <= s:
            raise CdsBuildError(f"empty output interval {chrom}:{s}-{e} ({gene})")
        if not any(chrom == pc and s >= ps and e <= pe for pc, ps, pe in panel[gene]):
            raise CdsBuildError(f"output interval {chrom}:{s}-{e} lies outside the panel intervals of {gene}")
    # merged output must be non-overlapping within a gene
    last = {}
    for chrom, s, e, gene in records:
        prev = last.get(gene)
        if prev and prev[0] == chrom and s < prev[2]:
            raise CdsBuildError(f"overlapping output intervals for {gene}")
        last[gene] = (chrom, s, e)


def run(mane, panel_path, out, audit_path=None):
    panel = read_panel_bed(panel_path)
    print(f"  Panel: {sum(len(v) for v in panel.values())} intervals, {len(panel)} genes", file=sys.stderr)
    transcripts = resolve_mane_transcripts(mane, panel)
    cds_by_gene = read_transcript_cds(mane, transcripts)
    records, audit = build_cds(panel, transcripts, cds_by_gene)
    validate_output(records, panel)

    with open(out, "w", newline="\n") as fh:
        for chrom, start, end, gene in records:
            fh.write(f"{chrom}\t{start}\t{end}\t{gene}\n")
    total_bp = sum(e - s for _, s, e, _ in records)
    if audit_path:
        with open(audit_path, "w", newline="\n") as fh:
            json.dump({"genes": len(audit), "intervals": len(records), "cds_bp": total_bp,
                       "per_gene": audit}, fh, indent=2, sort_keys=True)
            fh.write("\n")
    return records, total_bp


def main(argv=None):
    args = parse_args(argv)
    try:
        records, total_bp = run(args.mane, args.panel, args.out, args.audit)
    except CdsBuildError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    genes = {g for _, _, _, g in records}
    print(f"CDS_BP={total_bp}")
    print(f"CDS_MB={total_bp / 1_000_000:.6f}")
    print(f"N_INTERVALS={len(records)}")
    print(f"N_GENES={len(genes)}")
    print(f"Wrote {args.out}: {len(records)} intervals, {len(genes)} genes, {total_bp:,} bp", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
