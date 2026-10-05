"""
GenRichi MSI Methodology Redesign -- DEVELOPMENT PROTOTYPE. NOT PRODUCTION CODE.

Phase A v2: reference-genome-wide microsatellite marker scanning, INDEPENDENT
of the 55-gene cancer panel BED. This is a separate, additive module -- it
does not modify scan_markers.py (v1, panel-restricted), which stays exactly
as tested and used for the existing 39-marker v1 panel.

Rationale (see the evidence-review report): every established sequencing-
native MSI method (MSIsensor/MSIsensor-pro, MANTIS) scans candidate
microsatellite sites from the reference genome itself, not from an unrelated
gene-capture panel. Restricting discovery to the 55-gene cancer BED is why
the v1 panel only had 39 loci; this module removes that restriction.

Parameters follow the MSIsensor-pro `scan` module precedent (independently
documented in this session's evidence review, not invented here):
  - homopolymer length: 8-50 bp (MSIsensor-pro scan defaults)
  - other repeat motifs: length 2-6 bp, minimum 5 repeat units
    (MSIsensor-pro scan defaults: microsatellite max length 6, min repeat
    times 5)
These are marker-GENERATION parameters, not classification thresholds --
dev/msi_redesign/msi_locus_model.py's classification parameters
(MIN_LOCUS_DEPTH, INSTABILITY_DISTANCE_CUTOFF, etc.) are untouched.

Scans the PRIMARY GRCh38 assembly only (chr1-22, chrX, chrY, chrM) -- not
the ~430 alt/decoy/random contigs also present in reference_db/ref/hg38.fa.
This is a deliberate, documented scope choice (matches standard practice;
alt/decoy contigs are not part of the analyzable single-copy genome for this
purpose), recorded in the manifest.

Uses compiled regex (C-level) rather than the v1 module's manual
character-by-character scan, because a manual pure-Python scan over ~3.09
Gb of primary-assembly sequence is impractically slow. Equivalence with the
v1 semantics (maximal, non-overlapping, longest-motif-preferred, minimal-
period-normalized) is covered by tests/test_msi_scan_markers_v2_dev.py.

Usage:
    python3 scan_markers_v2.py --reference /path/to/hg38.fa \
        --out-dir dev/msi_redesign/output --samtools samtools --version-tag v2
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
from datetime import datetime, timezone

MIN_HOMOPOLYMER_LEN = 8
MAX_HOMOPOLYMER_LEN = 50
STR_MOTIF_MIN_LEN = 2
STR_MOTIF_MAX_LEN = 6
MIN_STR_UNITS = 5

PRIMARY_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]

_HOMOPOLYMER_RE = re.compile(r"([ACGT])\1{%d,}" % (MIN_HOMOPOLYMER_LEN - 1))
_STR_RE = re.compile(r"([ACGT]{%d,%d})\1{%d,}" % (STR_MOTIF_MIN_LEN, STR_MOTIF_MAX_LEN, MIN_STR_UNITS - 1))


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def fetch_fai_lengths(fai_path: str) -> dict:
    lengths = {}
    with open(fai_path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            lengths[parts[0]] = int(parts[1])
    return lengths


def extract_chromosome_sequence(samtools: str, reference: str, chrom: str) -> str:
    out = subprocess.run(
        [samtools, "faidx", reference, chrom],
        capture_output=True, text=True, timeout=600,
    )
    if out.returncode != 0:
        raise RuntimeError(f"samtools faidx failed for {chrom}: {out.stderr.strip()}")
    lines = out.stdout.splitlines()
    return "".join(lines[1:]).upper()


def minimal_period(motif: str) -> str:
    """Reduce a repeated unit to its smallest periodic substring, e.g.
    'ATAT' -> 'AT', 'ACGACG' -> 'ACG'. Returns motif unchanged if primitive."""
    n = len(motif)
    for period_len in range(1, n):
        if n % period_len != 0:
            continue
        candidate = motif[:period_len]
        if candidate * (n // period_len) == motif:
            return candidate
    return motif


def find_homopolymer_runs(seq: str, min_len: int = MIN_HOMOPOLYMER_LEN, max_len: int = MAX_HOMOPOLYMER_LEN):
    """Yields (start, end, base, length) for maximal runs of a single base,
    min_len <= length <= max_len (runs longer than max_len are skipped, not
    truncated -- they are not reported as MSIsensor-pro-style markers)."""
    for m in _HOMOPOLYMER_RE.finditer(seq):
        start, end = m.start(), m.end()
        length = end - start
        if length <= max_len:
            yield start, end, seq[start], length


def find_short_tandem_repeats(seq: str, min_units: int = MIN_STR_UNITS):
    """Yields (start, end, motif, n_units) for tandem repeats of a 2-6bp
    primitive motif with >= min_units repeat units. A single combined regex
    pass naturally gives maximal, non-overlapping matches with the greedy
    engine preferring the longest captured group at each position. Matches
    reducible to a period-1 motif (i.e. actually a homopolymer) are
    excluded here -- they are reported by find_homopolymer_runs instead."""
    for m in _STR_RE.finditer(seq):
        start, end = m.start(), m.end()
        raw_motif = m.group(1)
        period = minimal_period(raw_motif)
        if len(period) < STR_MOTIF_MIN_LEN:
            continue  # period-1 -> this is a homopolymer, not a STR
        units = (end - start) // len(period)
        yield start, end, period, units


def _overlaps_any(start: int, end: int, sorted_intervals) -> bool:
    """sorted_intervals: list of (start, end) sorted by start, non-overlapping."""
    import bisect
    starts = [iv[0] for iv in sorted_intervals]
    idx = bisect.bisect_right(starts, start) - 1
    for j in (idx, idx + 1):
        if 0 <= j < len(sorted_intervals):
            ivs, ive = sorted_intervals[j]
            if ivs < end and start < ive:
                return True
    return False


def scan_chromosome(chrom: str, seq: str) -> list:
    markers = []
    homopolymer_intervals = []
    for start, end, base, length in find_homopolymer_runs(seq):
        homopolymer_intervals.append((start, end))
        markers.append({
            "marker_id": f"MS2_{chrom}_{start}_{end}", "chrom": chrom,
            "start": start, "end": end, "motif": base,
            "ref_repeat_length": length, "kind": "mononucleotide",
        })
    homopolymer_intervals.sort()

    for start, end, motif, units in find_short_tandem_repeats(seq):
        if _overlaps_any(start, end, homopolymer_intervals):
            continue
        markers.append({
            "marker_id": f"MS2_{chrom}_{start}_{end}", "chrom": chrom,
            "start": start, "end": end, "motif": motif,
            "ref_repeat_length": units, "kind": "str",
        })
    markers.sort(key=lambda m: (m["start"], m["end"]))
    return markers


def write_marker_tsv(markers, out_path, reference_build):
    with open(out_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("marker_id\tchrom\tstart\tend\tmotif\tref_repeat_length\treference_build\tkind\n")
        for m in markers:
            fh.write(f"{m['marker_id']}\t{m['chrom']}\t{m['start']}\t{m['end']}\t{m['motif']}\t"
                      f"{m['ref_repeat_length']}\t{reference_build}\t{m['kind']}\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--reference", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--reference-build", default="GRCh38")
    ap.add_argument("--version-tag", default="v2")
    ap.add_argument("--chroms", nargs="*", default=None, help="override chromosome list (default: primary assembly)")
    args = ap.parse_args(argv)

    os.makedirs(args.out_dir, exist_ok=True)
    fai_path = args.reference + ".fai"
    lengths = fetch_fai_lengths(fai_path)
    chroms = args.chroms if args.chroms else PRIMARY_CHROMS
    for c in chroms:
        if c not in lengths:
            raise SystemExit(f"chromosome {c} not found in {fai_path}")

    all_markers = []
    per_chrom_counts = {}
    t0 = datetime.now(timezone.utc)
    for chrom in chroms:
        seq = extract_chromosome_sequence(args.samtools, args.reference, chrom)
        chrom_markers = scan_chromosome(chrom, seq)
        per_chrom_counts[chrom] = len(chrom_markers)
        all_markers.extend(chrom_markers)
        print(f"  {chrom}: {len(seq)} bp -> {len(chrom_markers)} markers", file=sys.stderr)
    elapsed = (datetime.now(timezone.utc) - t0).total_seconds()

    marker_path = os.path.join(args.out_dir, f"msi_markers.{args.version_tag}.tsv")
    write_marker_tsv(all_markers, marker_path, args.reference_build)
    marker_hash = sha256_file(marker_path)

    motif_dist = {}
    kind_dist = {}
    for m in all_markers:
        motif_dist[m["motif"]] = motif_dist.get(m["motif"], 0) + 1
        kind_dist[m["kind"]] = kind_dist.get(m["kind"], 0) + 1

    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "reference_build": args.reference_build,
        "source_reference_path": args.reference,
        "source_reference_fai_sha256": sha256_file(fai_path),
        "scan_scope": "primary_assembly_only",
        "chromosomes_scanned": chroms,
        "scan_parameters": {
            "min_homopolymer_length": MIN_HOMOPOLYMER_LEN,
            "max_homopolymer_length": MAX_HOMOPOLYMER_LEN,
            "str_motif_min_length": STR_MOTIF_MIN_LEN,
            "str_motif_max_length": STR_MOTIF_MAX_LEN,
            "min_str_repeat_units": MIN_STR_UNITS,
            "parameter_source": "MSIsensor-pro 'scan' module defaults (independently documented; see evidence review)",
        },
        "marker_file_path": marker_path,
        "marker_file_sha256": marker_hash,
        "n_markers": len(all_markers),
        "n_markers_per_chromosome": per_chrom_counts,
        "motif_distribution": motif_dist,
        "kind_distribution": kind_dist,
        "scan_wall_clock_seconds": elapsed,
        "not_derived_from_cancer_panel_bed": True,
    }
    manifest_path = os.path.join(args.out_dir, f"msi_markers.{args.version_tag}.manifest.json")
    with open(manifest_path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
        fh.write("\n")

    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
