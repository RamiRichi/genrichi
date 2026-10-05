"""
GenRichi MSI Methodology Redesign -- DEVELOPMENT PROTOTYPE. NOT PRODUCTION CODE.

Tier 1: independent GRCh38 verification of the five Promega MSI Analysis
System INSTABILITY markers (BAT-25, BAT-26, NR-21, NR-24, MONO-27).

Penta C and Penta D are DELIBERATELY EXCLUDED -- in the Promega system they
are pentanucleotide sample-identity/QC markers, not MSI instability markers,
and must not be used in any instability calculation (per the Promega TM255
technical manual and explicit user instruction for this development phase).

Verification method (does not trust any paper's stated coordinates blindly):
for each marker, this script (1) starts from a documented candidate genomic
window derived from an independently sourced primer pair, (2) searches for
both primer sequences (allowing either orientation) inside that window
against the REAL production hg38.fa, (3) locates the mononucleotide repeat
between the two confirmed primer positions, and (4) records the exact
GRCh38 coordinates, repeat motif/length, and a SHA-256 verification hash of
(coordinates + primers + extracted reference context). If a primer cannot
be confirmed by exact sequence match in the given window, the marker is
reported as UNVERIFIED rather than silently guessed.

Sources (see the evidence-review report delivered earlier in this session):
  - BAT-25, BAT-26, NR-21, NR-24 primer positions/sequences: Table SII,
    supplementary material of an in-house NGS MSI test validation paper
    (Spandidos Publications; hg38 primer positions as printed), cross-
    checked against the same four primer sequences independently given in
    US Patent 7951564B2 ("Primers for amplifying mononucleotide satellite
    markers") -- both sources agree on primer sequence.
  - MONO-27 / NR-27: US Patent 7951564B2 states it is a 27-adenine repeat
    in the 5' region of GenBank AF070674 ("inhibitor of apoptosis
    protein-1" mRNA, alias MIHC). MIHC is the HGNC-registered alias for
    BIRC3 (NOT BIRC2/API1, despite the superficially similar "API1" alias
    on BIRC2 -- this was cross-checked via HGNC before use). The primer
    pair and a ~30bp anchor immediately flanking the repeat were located
    by exact string match against the real AF070674 cDNA sequence (ENA),
    then that anchor was located in the real hg38.fa within the BIRC2/
    BIRC3 gene-cluster region on chr11.

Usage:
    python3 promega_markers.py --reference /path/to/hg38.fa \
        --out-dir dev/msi_redesign/output --samtools samtools
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone

MIN_HOMOPOLYMER_LEN = 6  # lenient, for locus-confirmation only (not marker generation)

# name -> (chrom, search_start, search_end, forward_primer, reverse_primer, source)
CANDIDATE_REGIONS = {
    "BAT-25": (
        "chr4", 54731000, 54733000,
        "TCGCCTCCAAGAATGTAAGT", "TCTGCATTTTAACTATGGCTC",
        "Spandidos Table SII (hg38 primer position, c-Kit/KIT intron); "
        "primer sequence cross-checked against US Patent 7951564B2",
    ),
    "BAT-26": (
        "chr2", 47413500, 47415500,
        "TGACTACTTTTGACTTCAGCC", "AACCATTCAACATTTTTAACCC",
        "Spandidos Table SII (hg38 primer position, MSH2 intron 5); "
        "primer sequence cross-checked against US Patent 7951564B2",
    ),
    "NR-21": (
        "chr14", 23182300, 23184000,
        "ATTCCTACTCCGCATTCACA", "TAAATGTATGTCTCCCCTGG",
        "Spandidos Table SII (hg38 primer position, SLC7A8 region); "
        "primer sequence cross-checked against US Patent 7951564B2",
    ),
    "NR-24": (
        "chr2", 95182800, 95184400,
        "CCATTGCTGAATTTTACCTC", "ATTGTGCCATTGCATTCCAA",
        "Spandidos Table SII (hg38 primer position, ZNF2 region); "
        "primer sequence cross-checked against US Patent 7951564B2",
    ),
    "MONO-27": (
        "chr11", 102150000, 102450000,
        "AACCATGCTTGCAAACCACT", "CGATAATACTAGCAATGACC",
        "US Patent 7951564B2 (27A repeat, GenBank AF070674 cDNA, alias MIHC "
        "= BIRC3 per HGNC lookup); primer pair located by exact match "
        "against the real AF070674 cDNA (ENA), then that flanking anchor "
        "located in hg38.fa within the BIRC2/BIRC3 gene cluster",
    ),
}


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def extract(samtools: str, reference: str, chrom: str, start_1based: int, end_1based: int) -> str:
    region = f"{chrom}:{start_1based}-{end_1based}"
    out = subprocess.run([samtools, "faidx", reference, region], capture_output=True, text=True, timeout=120)
    if out.returncode != 0:
        raise RuntimeError(f"samtools faidx failed for {region}: {out.stderr.strip()}")
    lines = out.stdout.splitlines()
    return "".join(lines[1:]).upper()


def find_homopolymer_runs(seq: str, min_len: int = MIN_HOMOPOLYMER_LEN):
    i, n = 0, len(seq)
    hits = []
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
            hits.append((i, j, base, length))
        i = j
    return hits


def locate_primer(seq: str, primer: str):
    """Returns (local_start, local_end, orientation) or None."""
    idx = seq.find(primer)
    if idx != -1:
        return idx, idx + len(primer), "forward_as_given"
    rc = revcomp(primer)
    idx = seq.find(rc)
    if idx != -1:
        return idx, idx + len(rc), "reverse_complement"
    return None


def verify_marker(samtools: str, reference: str, name: str):
    chrom, search_start, search_end, fwd, rev, source = CANDIDATE_REGIONS[name]
    window_seq = extract(samtools, reference, chrom, search_start, search_end)

    fwd_hit = locate_primer(window_seq, fwd)
    rev_hit = locate_primer(window_seq, rev)
    if fwd_hit is None or rev_hit is None:
        return {
            "marker_id": name, "status": "UNVERIFIED",
            "reason": f"primer not found in window (fwd={fwd_hit is not None}, rev={rev_hit is not None})",
            "search_window": f"{chrom}:{search_start}-{search_end}",
            "source": source,
        }

    positions = sorted([fwd_hit[0], fwd_hit[1], rev_hit[0], rev_hit[1]])
    interior_local_start, interior_local_end = positions[1], positions[2]
    interior_seq = window_seq[interior_local_start:interior_local_end]
    homopolymer_hits = find_homopolymer_runs(interior_seq, min_len=MIN_HOMOPOLYMER_LEN)

    if not homopolymer_hits:
        return {
            "marker_id": name, "status": "UNVERIFIED",
            "reason": "both primers located but no homopolymer >=6bp found between them",
            "search_window": f"{chrom}:{search_start}-{search_end}",
            "interior_sequence": interior_seq,
            "source": source,
        }

    local_s, local_e, base, length = max(homopolymer_hits, key=lambda h: h[3])
    g_start_1based = search_start + interior_local_start + local_s
    g_end_1based = search_start + interior_local_start + local_e - 1
    bed_start = g_start_1based - 1  # 0-based half-open
    bed_end = g_end_1based

    context_start = max(search_start, g_start_1based - 10)
    context_end = g_end_1based + 10
    context_seq = extract(samtools, reference, chrom, context_start, context_end)

    verification_payload = f"{chrom}\t{bed_start}\t{bed_end}\t{base}\t{length}\t{fwd}\t{rev}\t{context_seq}"
    verification_hash = sha256_text(verification_payload)

    return {
        "marker_id": name, "status": "VERIFIED",
        "chrom": chrom, "start": bed_start, "end": bed_end,
        "motif": base, "ref_repeat_length": length,
        "forward_primer": fwd, "reverse_primer": rev,
        "forward_primer_orientation": fwd_hit[2], "reverse_primer_orientation": rev_hit[2],
        "reference_build": "GRCh38",
        "context_sequence": context_seq,
        "context_region": f"{chrom}:{context_start}-{context_end}",
        "verification_hash_sha256": verification_hash,
        "source": source,
        "role": "msi_instability_marker",
    }


def write_marker_tsv(verified_markers, out_path):
    with open(out_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("marker_id\tchrom\tstart\tend\tmotif\tref_repeat_length\treference_build\trole\tverification_hash_sha256\n")
        for m in verified_markers:
            fh.write(f"{m['marker_id']}\t{m['chrom']}\t{m['start']}\t{m['end']}\t{m['motif']}\t"
                      f"{m['ref_repeat_length']}\t{m['reference_build']}\t{m['role']}\t{m['verification_hash_sha256']}\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--reference", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--samtools", default="samtools")
    args = ap.parse_args(argv)

    os.makedirs(args.out_dir, exist_ok=True)
    results = {name: verify_marker(args.samtools, args.reference, name) for name in CANDIDATE_REGIONS}

    verified = [r for r in results.values() if r["status"] == "VERIFIED"]
    unverified = [r for r in results.values() if r["status"] != "VERIFIED"]

    marker_path = os.path.join(args.out_dir, "msi_markers.promega_tier1.tsv")
    write_marker_tsv(verified, marker_path)
    marker_hash = sha256_file(marker_path) if verified else None

    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "tier": 1,
        "description": "Promega MSI Analysis System instability markers (Penta C/D sample-ID markers deliberately excluded)",
        "source_reference_path": args.reference,
        "n_markers_verified": len(verified),
        "n_markers_unverified": len(unverified),
        "marker_file_path": marker_path,
        "marker_file_sha256": marker_hash,
        "results": results,
    }
    manifest_path = os.path.join(args.out_dir, "msi_markers.promega_tier1.manifest.json")
    with open(manifest_path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
        fh.write("\n")

    print(json.dumps(manifest, indent=2, sort_keys=True))
    if unverified:
        print(f"\nWARNING: {len(unverified)} marker(s) could not be verified: "
              f"{[m['marker_id'] for m in unverified]}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
