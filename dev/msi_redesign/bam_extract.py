"""
GenRichi MSI Methodology Redesign -- DEVELOPMENT PROTOTYPE. NOT PRODUCTION CODE.

Phase B: tumor/normal BAM locus extraction -> repeat-length distributions.

Deliberately uses `samtools view` (already installed, no new dependency --
matches the "Option D" recommendation from the design document) rather than
adding pysam. All read-filtering and CIGAR-length logic is implemented as
PURE functions (no I/O) so it can be unit-tested against fabricated SAM
records/CIGAR strings without any real BAM file -- the I/O layer
(`fetch_reads_for_locus`, calling the real `samtools view`) is a thin
wrapper used only when real data is actually available (Phase D).

Read filtering policy (documented, not hidden in code):
  - exclude unmapped        (FLAG & 0x4)
  - exclude secondary       (FLAG & 0x100)
  - exclude supplementary   (FLAG & 0x800)
  - exclude PCR/optical duplicates (FLAG & 0x400) -- duplicates are already
    marked upstream by the production mark_duplicates_paired rule; this
    extractor reads the same duplicate-marked BAM and honours those flags
    rather than re-implementing duplicate detection
  - require MAPQ >= MIN_MAPQ
  - require >= MIN_FLANK bp of aligned (M/=/X) reference-consuming sequence
    strictly before locus_start and strictly after locus_end (a read that
    only clips into the repeat, or starts/ends inside it, is excluded --
    it cannot reliably show the repeat's true boundaries)
"""

from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass

MIN_MAPQ = 20
MIN_FLANK = 5

_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")

FLAG_UNMAPPED = 0x4
FLAG_SECONDARY = 0x100
FLAG_DUPLICATE = 0x400
FLAG_SUPPLEMENTARY = 0x800
FLAG_EXCLUDE_MASK = FLAG_UNMAPPED | FLAG_SECONDARY | FLAG_DUPLICATE | FLAG_SUPPLEMENTARY

REF_CONSUMING = set("MDN=X")
READ_CONSUMING = set("MIS=X")


@dataclass
class SamRecord:
    """The minimal subset of a SAM record this module needs."""
    flag: int
    pos: int          # 1-based leftmost mapped position (SAM convention)
    mapq: int
    cigar: str


def parse_cigar(cigar: str):
    """'76M2I10M' -> [(76,'M'), (2,'I'), (10,'M')]. '*' -> []."""
    if cigar in ("", "*"):
        return []
    return [(int(length), op) for length, op in _CIGAR_RE.findall(cigar)]


def passes_read_filters(record: SamRecord, min_mapq: int = MIN_MAPQ) -> bool:
    if record.flag & FLAG_EXCLUDE_MASK:
        return False
    if record.mapq < min_mapq:
        return False
    if not record.cigar or record.cigar == "*":
        return False
    return True


def observed_repeat_length(
    record: SamRecord, locus_start: int, locus_end: int,
    min_flank: int = MIN_FLANK,
) -> "int | None":
    """
    Length (in bases) of this read's representation of the reference
    interval [locus_start, locus_end) (0-based, half-open), accounting for
    insertions/deletions inside that interval.

    Returns None if the read does not have >= min_flank bp of aligned
    reference-consuming sequence strictly outside the interval on BOTH
    sides (insufficient flank -> this read cannot be used for this locus).

    ref_pos is 0-based and tracks the next reference base the CIGAR walk is
    about to consume; locus_start/locus_end are also 0-based, matching the
    marker file's coordinate convention (see scan_markers.py).
    """
    ops = parse_cigar(record.cigar)
    if not ops:
        return None

    ref_pos = record.pos - 1  # SAM POS is 1-based; convert to 0-based
    flank_before = 0
    flank_after = 0
    observed_len = 0
    entered_locus = False
    exited_locus = False

    for length, op in ops:
        if op in ("S", "H", "P"):
            continue  # clipping/padding: no reference or repeat-length contribution
        if op in ("M", "=", "X"):
            span_start, span_end = ref_pos, ref_pos + length
            # portion strictly before the locus
            before = max(0, min(span_end, locus_start) - span_start)
            flank_before += before if not entered_locus else 0
            # portion inside the locus
            inside = max(0, min(span_end, locus_end) - max(span_start, locus_start))
            if inside > 0:
                entered_locus = True
                observed_len += inside
            # portion strictly after the locus
            if span_end > locus_end:
                exited_locus = True
                after = span_end - max(span_start, locus_end)
                flank_after += max(0, after)
            elif exited_locus:
                flank_after += length
            ref_pos += length
        elif op == "D":
            span_start, span_end = ref_pos, ref_pos + length
            inside = max(0, min(span_end, locus_end) - max(span_start, locus_start))
            if inside > 0:
                entered_locus = True
                # a deletion inside the locus consumes reference but no read
                # bases -> contributes 0 to observed_len (already correct: add nothing)
            ref_pos += length
        elif op == "I":
            # zero-width in reference space; attribute to ref_pos if that
            # position falls inside the locus (an insertion "at" this point)
            if locus_start <= ref_pos < locus_end:
                entered_locus = True
                observed_len += length
            # I does not advance ref_pos
        elif op == "N":
            ref_pos += length  # skipped reference (rare in DNA-seq; treat like a large gap)

    if not entered_locus:
        return None
    if flank_before < min_flank or flank_after < min_flank:
        return None
    return observed_len


def fetch_reads_for_locus(samtools: str, bam_path: str, chrom: str, start: int, end: int):
    """Real I/O: `samtools view` a 1-based inclusive region and yield SamRecord.
    BED-style [start, end) 0-based -> samtools region chrom:(start+1)-end."""
    region = f"{chrom}:{start + 1}-{end}"
    proc = subprocess.run(
        [samtools, "view", bam_path, region],
        capture_output=True, text=True, timeout=120,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"samtools view failed for {region} in {bam_path}: {proc.stderr.strip()}")
    for line in proc.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 6:
            continue
        yield SamRecord(flag=int(fields[1]), pos=int(fields[3]), mapq=int(fields[4]), cigar=fields[5])


def locus_lengths(samtools: str, bam_path: str, chrom: str, start: int, end: int,
                   min_mapq: int = MIN_MAPQ, min_flank: int = MIN_FLANK):
    """Real I/O + pure logic combined: observed repeat lengths for every
    filter-passing read overlapping this locus in this BAM."""
    lengths = []
    for record in fetch_reads_for_locus(samtools, bam_path, chrom, start, end):
        if not passes_read_filters(record, min_mapq=min_mapq):
            continue
        length = observed_repeat_length(record, start, end, min_flank=min_flank)
        if length is not None:
            lengths.append(length)
    return lengths
