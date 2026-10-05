"""
GenRichi MSI Methodology Redesign -- DEVELOPMENT PROTOTYPE. NOT PRODUCTION CODE.

Phase A: deterministic microsatellite marker scanning.

Scans the EXACT existing production reference (reference_db/ref/hg38.fa)
within the EXACT existing panel BED
(resources/panel/phase1_solid_tumor/solid_tumor_phase1_v1.bed), using
`samtools faidx` to extract only the panel-region sequence (never loading
the 3.2GB genome into memory), and flags:
  - mononucleotide runs of length >= MIN_MONO_LEN
  - 2-4bp short tandem repeats with >= MIN_STR_UNITS[motif_len] units

The scan is a pure function of (reference sequence, panel BED, these
parameters): re-running it against the same inputs always produces a
byte-identical marker file. No coordinate in the output is copied from
memory or any external source -- every one is derived directly from the
reference sequence this script actually reads.

Usage:
    python3 scan_markers.py --reference /path/to/hg38.fa --panel-bed /path/to/panel.bed \
        --out-dir dev/msi_redesign/output --samtools samtools

Output (all written under --out-dir, never under workflow/ or reference_db/):
    msi_markers.v1.tsv       -- the marker panel itself
    msi_markers.v1.manifest.json  -- scan parameters, source hashes, marker-file hash
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone

MIN_MONO_LEN = 10                 # mononucleotide run, min length (bp)
MIN_STR_UNITS = {2: 5, 3: 5, 4: 5}  # motif length -> minimum repeat units
CHROM_ORDER = [f"chr{i}" for i in list(range(1, 23)) + ["X", "Y", "M"]]


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_panel_bed(path: str):
    """BED4: chrom, start, end, gene. 0-based half-open, as documented in
    resources/panel/phase1_solid_tumor/README.md."""
    rows = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            gene = parts[3] if len(parts) > 3 else ""
            rows.append((chrom, start, end, gene))
    return rows


def extract_region_sequence(samtools: str, reference: str, chrom: str, start: int, end: int) -> str:
    """1-based, inclusive region for samtools faidx; BED is 0-based
    half-open, so region = chrom:(start+1)-end."""
    region = f"{chrom}:{start + 1}-{end}"
    out = subprocess.run(
        [samtools, "faidx", reference, region],
        capture_output=True, text=True, timeout=60,
    )
    if out.returncode != 0:
        raise RuntimeError(f"samtools faidx failed for {region}: {out.stderr.strip()}")
    lines = out.stdout.splitlines()
    return "".join(lines[1:]).upper()  # drop the ">header" line


def find_mononucleotide_runs(seq: str, min_len: int = MIN_MONO_LEN):
    """Yields (local_start, local_end, base, length) for each maximal run of
    a single repeated base with length >= min_len. Non-overlapping,
    leftmost-greedy, deterministic single left-to-right pass."""
    i = 0
    n = len(seq)
    while i < n:
        base = seq[i]
        if base not in "ACGT":
            i += 1
            continue
        j = i + 1
        while j < n and seq[j] == base:
            j += 1
        length = j - i
        if length >= min_len:
            yield (i, j, base, length)
        i = j


def find_short_tandem_repeats(seq: str, min_units: dict = MIN_STR_UNITS):
    """Yields (local_start, local_end, motif, n_units) for 2-4bp motifs
    repeated >= min_units[len(motif)] times. Skips any span already
    consumed by a mononucleotide run or an earlier (longer-motif-first)
    STR match, so no base is double-counted -- deterministic, single pass
    per motif length, longest motif checked first at each position."""
    n = len(seq)
    consumed = [False] * n
    mono_spans = [(s, e) for s, e, _, _ in find_mononucleotide_runs(seq, MIN_MONO_LEN)]
    for s, e in mono_spans:
        for k in range(s, e):
            consumed[k] = True

    results = []
    for motif_len in (4, 3, 2):  # longest motif first, so e.g. AAAA isn't mis-split as AA AA
        min_u = min_units[motif_len]
        min_span = motif_len * min_u
        i = 0
        while i <= n - min_span:
            if consumed[i]:
                i += 1
                continue
            motif = seq[i:i + motif_len]
            if not all(c in "ACGT" for c in motif) or len(set(motif)) == 1:
                i += 1
                continue  # a "motif" of one repeated base is just a mononucleotide run
            j = i + motif_len
            units = 1
            while j + motif_len <= n and not any(consumed[j:j + motif_len]) and seq[j:j + motif_len] == motif:
                j += motif_len
                units += 1
            if units >= min_u:
                results.append((i, j, motif, units))
                for k in range(i, j):
                    consumed[k] = True
                i = j
            else:
                i += 1
    return sorted(results)


def chrom_sort_key(chrom: str):
    try:
        return (CHROM_ORDER.index(chrom), 0)
    except ValueError:
        return (len(CHROM_ORDER), chrom)


def scan_panel(samtools: str, reference: str, panel_bed_rows) -> list:
    markers = []
    for chrom, start, end, gene in panel_bed_rows:
        seq = extract_region_sequence(samtools, reference, chrom, start, end)
        for local_s, local_e, base, length in find_mononucleotide_runs(seq):
            g_start, g_end = start + local_s, start + local_e
            markers.append({
                "marker_id": f"MS_{chrom}_{g_start}_{g_end}",
                "chrom": chrom, "start": g_start, "end": g_end,
                "motif": base, "ref_repeat_length": length, "gene": gene,
                "kind": "mononucleotide",
            })
        for local_s, local_e, motif, units in find_short_tandem_repeats(seq):
            g_start, g_end = start + local_s, start + local_e
            markers.append({
                "marker_id": f"MS_{chrom}_{g_start}_{g_end}",
                "chrom": chrom, "start": g_start, "end": g_end,
                "motif": motif, "ref_repeat_length": units, "gene": gene,
                "kind": "str",
            })
    markers.sort(key=lambda m: (chrom_sort_key(m["chrom"]), m["start"], m["end"]))
    # de-duplicate identical (chrom,start,end) that could arise if the same
    # gene appears twice in the BED (adjacent/overlapping exon padding)
    seen = set()
    deduped = []
    for m in markers:
        key = (m["chrom"], m["start"], m["end"])
        if key in seen:
            continue
        seen.add(key)
        deduped.append(m)
    return deduped


def write_marker_tsv(markers, out_path, reference_build):
    with open(out_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("marker_id\tchrom\tstart\tend\tmotif\tref_repeat_length\treference_build\tgene\tkind\n")
        for m in markers:
            fh.write(f"{m['marker_id']}\t{m['chrom']}\t{m['start']}\t{m['end']}\t{m['motif']}\t"
                      f"{m['ref_repeat_length']}\t{reference_build}\t{m['gene']}\t{m['kind']}\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--reference", required=True)
    ap.add_argument("--panel-bed", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--reference-build", default="GRCh38")
    ap.add_argument("--version-tag", default="v1")
    args = ap.parse_args(argv)

    os.makedirs(args.out_dir, exist_ok=True)
    panel_rows = read_panel_bed(args.panel_bed)
    markers = scan_panel(args.samtools, args.reference, panel_rows)

    marker_path = os.path.join(args.out_dir, f"msi_markers.{args.version_tag}.tsv")
    write_marker_tsv(markers, marker_path, args.reference_build)
    marker_hash = sha256_file(marker_path)

    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "reference_build": args.reference_build,
        "source_reference_path": args.reference,
        "source_panel_bed_path": args.panel_bed,
        "source_panel_bed_sha256": sha256_file(args.panel_bed),
        "scan_parameters": {
            "min_mononucleotide_length": MIN_MONO_LEN,
            "min_str_units_by_motif_length": MIN_STR_UNITS,
        },
        "marker_file_path": marker_path,
        "marker_file_sha256": marker_hash,
        "n_panel_bed_regions": len(panel_rows),
        "n_markers": len(markers),
    }
    manifest_path = os.path.join(args.out_dir, f"msi_markers.{args.version_tag}.manifest.json")
    with open(manifest_path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
        fh.write("\n")

    print(f"wrote {len(markers)} markers -> {marker_path}")
    print(f"marker_file_sha256={marker_hash}")
    print(f"manifest -> {manifest_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
