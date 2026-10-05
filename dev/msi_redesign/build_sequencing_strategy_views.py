"""
GenRichi MSI Methodology Redesign -- DEVELOPMENT PROTOTYPE. NOT PRODUCTION CODE.

Sequencing-strategy VIEWS over the master filtered panel
(tier2_filtered.v1.tsv, 1,208,169 loci) -- NOT a replacement for it, and
NOT a new filtering pass. The master panel is never modified or
duplicated destructively; each view is a membership subset derived from
genomic-footprint overlap that was already computed, so a future
WES+WGS-combined (or a different capture kit, or a real targeted-panel)
decision never requires re-deriving the expensive genome-wide Umap/
RepeatMasker/segdup annotation -- only a cheap interval-overlap pass
against a new footprint BED.

  View A -- WES-compatible: loci FULLY inside the Broad generic exome
    footprint (see broad_exome_footprint_analysis.py). 93,687 loci.
  View B -- WGS: the entire master filtered panel, unchanged. Not
    duplicated as a separate file -- documented as "= tier2_filtered.v1.tsv"
    with that file's own SHA-256 as the provenance anchor.
  View C -- targeted-panel feasibility: DEFERRED. No capture technology
    has been chosen yet. This view is a documented placeholder, not
    populated with any data, so its later construction is a drop-in reuse
    of build_view() below once a real capture BED is available.

Statistical/biological down-selection WITHIN a view (the next planned
step) is explicitly NOT performed by this script.

Provenance carried forward, not silently dropped: the Broad reference
used to build the WES footprint has identical chromosome LENGTHS to our
production hg38.fa across all 455 shared contigs, but different MD5/M5
checksums for chr1, chr17, and chrX (chrM matches exactly). This remains
an open, undiagnosed provenance question -- not treated as a problem, but
not silently resolved either.

Usage:
    python3 build_sequencing_strategy_views.py \
        --master-panel dev/msi_redesign/output/tier2_filtered.v1.tsv \
        --footprint-coverage dev/msi_redesign/reference_annotations/filtered_vs_footprint_coverage.tsv \
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


def load_footprint_coverage(path: str):
    """bedtools coverage output -> {marker_id: fraction_covered}."""
    result = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            result[parts[3]] = float(parts[-1])
    return result


def select_fully_inside(coverage: dict, min_fraction: float = 0.999999) -> set:
    return {marker_id for marker_id, frac in coverage.items() if frac >= min_fraction}


def build_view(master_panel_path: str, member_ids: set, out_path: str):
    """Writes the subset of master_panel_path whose marker_id is in
    member_ids, preserving row order and all columns -- a pure membership
    filter, no new annotation or scoring logic."""
    n_written = 0
    with open(master_panel_path, encoding="utf-8") as src, \
         open(out_path, "w", encoding="utf-8", newline="\n") as dst:
        reader = csv.reader(src, delimiter="\t")
        writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
        header = next(reader)
        writer.writerow(header)
        marker_col = header.index("marker_id")
        for row in reader:
            if row[marker_col] in member_ids:
                writer.writerow(row)
                n_written += 1
    return n_written


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--master-panel", required=True)
    ap.add_argument("--footprint-coverage", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args(argv)

    os.makedirs(args.out_dir, exist_ok=True)
    master_sha256 = sha256_file(args.master_panel)

    coverage = load_footprint_coverage(args.footprint_coverage)
    view_a_ids = select_fully_inside(coverage)

    view_a_path = os.path.join(args.out_dir, "view_A_wes_compatible.v1.tsv")
    n_view_a = build_view(args.master_panel, view_a_ids, view_a_path)
    view_a_sha256 = sha256_file(view_a_path)

    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "master_panel": {
            "path": args.master_panel,
            "sha256": master_sha256,
            "n_markers": None,  # filled below
        },
        "views": {
            "view_A_wes_compatible": {
                "definition": "loci FULLY inside the Broad generic exome footprint "
                               "(broad_exome_footprint_analysis.py fully_inside_footprint category)",
                "derived_from_master_panel": True,
                "master_panel_modified": False,
                "path": view_a_path,
                "sha256": view_a_sha256,
                "n_markers": n_view_a,
            },
            "view_B_wgs": {
                "definition": "the entire master filtered panel, unchanged",
                "derived_from_master_panel": True,
                "master_panel_modified": False,
                "path": args.master_panel,
                "sha256": master_sha256,
                "n_markers": None,  # = master panel n_markers, filled below
                "note": "not duplicated as a separate file -- this view IS tier2_filtered.v1.tsv",
            },
            "view_C_targeted_panel": {
                "definition": "DEFERRED -- no capture technology decided yet",
                "derived_from_master_panel": None,
                "master_panel_modified": False,
                "path": None,
                "sha256": None,
                "n_markers": None,
                "note": "reuse build_view() with a real capture BED's overlap-derived member_ids "
                        "once a targeted-panel capture technology is chosen; no re-derivation of "
                        "Umap/RepeatMasker/segdup annotation is needed to build this later",
            },
        },
        "next_step_not_yet_performed": "statistical/biological down-selection WITHIN each view",
        "reference_provenance_caveat": (
            "The Broad exome footprint reference has identical chromosome LENGTHS to our "
            "production hg38.fa across all 455 shared contigs, but different MD5/M5 checksums "
            "for chr1, chr17, and chrX (chrM matches exactly). Not treated as a defect; kept "
            "explicitly documented until the exact GRCh38 distribution difference is identified."
        ),
    }

    # fill n_markers for master panel / view B by counting once
    with open(args.master_panel, encoding="utf-8") as fh:
        n_master = sum(1 for _ in fh) - 1
    manifest["master_panel"]["n_markers"] = n_master
    manifest["views"]["view_B_wgs"]["n_markers"] = n_master

    manifest_path = os.path.join(args.out_dir, "sequencing_strategy_views.v1.manifest.json")
    with open(manifest_path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
        fh.write("\n")

    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
