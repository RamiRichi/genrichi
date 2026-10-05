"""
GenRichi MSI Methodology Redesign -- DEVELOPMENT PROTOTYPE. NOT PRODUCTION CODE.

Combines the two independently generated marker tiers into one versioned
development marker definition:

  Tier 1 -- the five verified Promega MSI Analysis System instability
            markers (BAT-25, BAT-26, NR-21, NR-24, MONO-27). Clinically
            established loci, independently GRCh38-verified in this
            session (see promega_markers.py). Penta C/D are NOT included
            (sample-identity markers, not instability markers).
  Tier 2 -- the genome-wide, cancer-panel-independent quantitative marker
            set from scan_markers_v2.py, generated with MSIsensor-pro
            'scan'-precedent parameters.

Neither tier is claimed to be clinically validated for GenRichi -- this
file only merges two already-generated, already-hashed inputs; it performs
no new marker discovery or classification logic itself.

Usage:
    python3 combine_markers.py \
        --tier1 dev/msi_redesign/output/msi_markers.promega_tier1.tsv \
        --tier2 dev/msi_redesign/output/msi_markers.v2.tsv \
        --out-dir dev/msi_redesign/output
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
from datetime import datetime, timezone


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tier1(path: str):
    rows = []
    with open(path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rows.append({
                "marker_id": row["marker_id"], "chrom": row["chrom"],
                "start": int(row["start"]), "end": int(row["end"]),
                "motif": row["motif"], "ref_repeat_length": int(row["ref_repeat_length"]),
                "kind": "mononucleotide", "tier": 1, "tier_label": "promega_clinical_anchor",
            })
    return rows


def read_tier2(path: str):
    rows = []
    with open(path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rows.append({
                "marker_id": row["marker_id"], "chrom": row["chrom"],
                "start": int(row["start"]), "end": int(row["end"]),
                "motif": row["motif"], "ref_repeat_length": int(row["ref_repeat_length"]),
                "kind": row["kind"], "tier": 2, "tier_label": "genome_scan_quantitative",
            })
    return rows


def write_combined_tsv(rows, out_path):
    with open(out_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("marker_id\tchrom\tstart\tend\tmotif\tref_repeat_length\tkind\ttier\ttier_label\n")
        for m in rows:
            fh.write(f"{m['marker_id']}\t{m['chrom']}\t{m['start']}\t{m['end']}\t{m['motif']}\t"
                      f"{m['ref_repeat_length']}\t{m['kind']}\t{m['tier']}\t{m['tier_label']}\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--tier1", required=True)
    ap.add_argument("--tier2", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args(argv)

    os.makedirs(args.out_dir, exist_ok=True)
    tier1_rows = read_tier1(args.tier1)
    tier2_rows = read_tier2(args.tier2)
    combined = tier1_rows + tier2_rows
    combined.sort(key=lambda m: (m["tier"], m["chrom"], m["start"]))

    out_path = os.path.join(args.out_dir, "msi_markers.dev_combined.v1.tsv")
    write_combined_tsv(combined, out_path)
    combined_hash = sha256_file(out_path)

    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "clinically_validated_for_genrichi": False,
        "note": "Neither tier is clinically validated for GenRichi. Tier 1 markers are "
                "clinically established in the field (Promega MSI Analysis System) but "
                "have not been validated on GenRichi's pipeline/data. Tier 2 markers are "
                "a development-only quantitative candidate set, not independently validated.",
        "tier1_source_file": args.tier1,
        "tier1_source_sha256": sha256_file(args.tier1),
        "tier1_n_markers": len(tier1_rows),
        "tier2_source_file": args.tier2,
        "tier2_source_sha256": sha256_file(args.tier2),
        "tier2_n_markers": len(tier2_rows),
        "combined_file_path": out_path,
        "combined_file_sha256": combined_hash,
        "combined_n_markers": len(combined),
    }
    manifest_path = os.path.join(args.out_dir, "msi_markers.dev_combined.v1.manifest.json")
    with open(manifest_path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
        fh.write("\n")

    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
