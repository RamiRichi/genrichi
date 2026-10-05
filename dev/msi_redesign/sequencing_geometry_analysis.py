"""
GenRichi MSI Methodology Redesign -- DEVELOPMENT PROTOTYPE. NOT PRODUCTION CODE.

Sequencing-geometry analysis to inform (NOT make) the candidate-locus
spacing decision.

The real HCC1395 tumor/normal BAMs are used ONLY as references for
sequencing geometry: read length, properly-paired template (insert) size,
and mate placement. No variant calls, no MSI calls, no tumor biology and no
per-marker performance are read or used; no locus is selected or removed.

Redundancy metric (evidence from geometry, derived directly from the
sampled fragments -- not a heuristic cut-off):
  For a candidate locus 1 (length L1) and a neighbor locus 2 (length L2)
  separated by a gap of d bp (gap = start2 - end1), consider every offset at
  which a fragment could overlap locus 1. A locus is MEASURABLE from a mate
  if it lies fully inside that mate's aligned span with at least `flank` bp
  of aligned sequence on each side (the same MIN_FLANK rule bam_extract.py
  applies). For each real fragment geometry (template size T, mate span R):
      same_read_joint  = P(locus 2 is measurable from the SAME mate that
                           measures locus 1  | locus 1 measurable)
      fragment_joint   = P(locus 2 is measurable from EITHER mate of the
                           fragment                | locus 1 measurable)
  Both are pooled over all sampled fragments (weighting each fragment by
  the number of offsets at which it informs locus 1). A value near 1 means
  measurements at the two loci come from the same reads/fragments (highly
  redundant, statistically correlated); near 0 means largely independent
  evidence.

Assumptions (stated, not hidden):
  - mate 2 is assumed to have the same aligned reference span as the
    observed read (100bp reads; soft-clipping variation ignored);
  - fragments are assumed uniformly placed relative to a locus (no
    capture-bait positional bias modeled), so this is a geometry estimate,
    not a measurement of any specific capture design;
  - sampling is restricted to the one well-covered region of the available
    dataset (chr17), so it is a reference for geometry only.

Nothing here turns a distance into a filter. It reports loci-affected counts
and geometry-based redundancy per spacing so a later decision can be
justified.

Usage:
    python3 sequencing_geometry_analysis.py --tumor-bam T.bam --normal-bam N.bam \
        --view-a view_A.tsv --view-b tier2_filtered.tsv --out result.json
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone

SPACINGS_BP = [20, 50, 75, 100, 150, 200, 250, 300]
CURVE_SPACINGS_BP = list(range(0, 410, 10))
LOCUS_LENGTHS_BP = [10, 25]   # sensitivity: short homopolymer vs. longer repeat
PRIMARY_LOCUS_LENGTH_BP = 10
MIN_FLANK = 5                 # same value bam_extract.py uses
PERCENTILES = [10, 25, 50, 75, 90, 95, 99]
_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
_REF_CONSUMING = set("MDN=X")


def percentile(sorted_values, q):
    """Linear-interpolation percentile (q in 0..100) of a sorted list."""
    if not sorted_values:
        return None
    if len(sorted_values) == 1:
        return float(sorted_values[0])
    pos = (len(sorted_values) - 1) * q / 100.0
    lo = int(pos)
    hi = min(lo + 1, len(sorted_values) - 1)
    frac = pos - lo
    return sorted_values[lo] * (1 - frac) + sorted_values[hi] * frac


def ref_span_from_cigar(cigar: str) -> int:
    return sum(int(n) for n, op in _CIGAR_RE.findall(cigar) if op in _REF_CONSUMING)


def merge_intervals(intervals):
    """Merge inclusive integer intervals [(lo,hi),...]."""
    merged = []
    for lo, hi in sorted(i for i in intervals if i[1] >= i[0]):
        if merged and lo <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))
    return merged


def total_length(intervals):
    return sum(hi - lo + 1 for lo, hi in intervals)


def intersect_intervals(a, b):
    """Intersection of two merged, sorted inclusive interval lists."""
    out = []
    i = j = 0
    while i < len(a) and j < len(b):
        lo = max(a[i][0], b[j][0])
        hi = min(a[i][1], b[j][1])
        if hi >= lo:
            out.append((lo, hi))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


def shift_intervals(intervals, delta):
    return [(lo + delta, hi + delta) for lo, hi in intervals]


def mates_for_fragment(template_size: int, read_span: int):
    """Mate aligned spans as half-open [a,b) relative to fragment start."""
    r = min(read_span, template_size)
    return [(0, r), (template_size - r, template_size)]


def measurable_starts(mate, locus_len: int, flank: int):
    """Inclusive interval of locus start offsets fully inside `mate` with
    `flank` aligned bp on each side, or None."""
    a, b = mate
    lo, hi = a + flank, b - flank - locus_len
    return (lo, hi) if hi >= lo else None


def fragment_joint_counts(template_size, read_span, len1, len2, gap, flank=MIN_FLANK):
    """Returns (n_locus1, n_same_read_joint, n_fragment_joint): counts of
    locus-1 start offsets for which locus 1 is measurable, and for which
    locus 2 (start = x + len1 + gap) is also measurable from the same mate /
    from either mate."""
    mates = mates_for_fragment(template_size, read_span)
    delta = len1 + gap
    per_mate_1 = [measurable_starts(m, len1, flank) for m in mates]
    per_mate_2 = [measurable_starts(m, len2, flank) for m in mates]

    set1 = merge_intervals([i for i in per_mate_1 if i])
    n1 = total_length(set1)
    if n1 == 0:
        return 0, 0, 0

    # locus 2 measurable from either mate: x + delta in mate-2 measurable set
    set2_shifted = merge_intervals(shift_intervals([i for i in per_mate_2 if i], -delta))
    n_fragment = total_length(intersect_intervals(set1, set2_shifted))

    # same mate for both loci
    same = []
    for m1, m2 in zip(per_mate_1, per_mate_2):
        if m1 and m2:
            inter = intersect_intervals([m1], shift_intervals([m2], -delta))
            same.extend(inter)
    n_same = total_length(merge_intervals(same))
    return n1, n_same, n_fragment


def pooled_redundancy(geometry_counter, len1, len2, gap, flank=MIN_FLANK):
    """geometry_counter: Counter{(template_size, read_span): n_fragments}."""
    tot1 = tot_same = tot_frag = 0
    for (t, r), count in geometry_counter.items():
        n1, ns, nf = fragment_joint_counts(t, r, len1, len2, gap, flank)
        tot1 += n1 * count
        tot_same += ns * count
        tot_frag += nf * count
    if tot1 == 0:
        return {"same_read_joint": None, "fragment_joint": None}
    return {"same_read_joint": tot_same / tot1, "fragment_joint": tot_frag / tot1}


def spacing_crossing(curve, key, level):
    """Smallest curve spacing at which `key` falls at/below `level`."""
    for gap, vals in sorted(curve.items()):
        v = vals[key]
        if v is not None and v <= level:
            return gap
    return None


def summarize_values(values):
    s = sorted(values)
    if not s:
        return {"n": 0}
    out = {"n": len(s), "min": s[0], "max": s[-1], "mean": sum(s) / len(s)}
    out["percentiles"] = {f"P{q}": percentile(s, q) for q in PERCENTILES}
    return out


def read_geometry_records(lines):
    """Parse `samtools view` text lines into (template_size, read_span,
    read_len) tuples. Expects first-in-pair, proper-pair, same-reference
    records; malformed / non-pair lines are skipped."""
    for line in lines:
        f = line.rstrip("\n").split("\t")
        if len(f) < 10 or f[6] != "=":
            continue
        try:
            tlen = abs(int(f[8]))
            span = ref_span_from_cigar(f[5])
            read_len = len(f[9])
        except ValueError:
            continue
        if tlen == 0 or span == 0:
            continue
        yield tlen, span, read_len


def sample_bam_geometry(samtools, bam, region, seed, fraction, min_mapq=20):
    cmd = [samtools, "view", "-f", "66", "-F", "3844", "-q", str(min_mapq),
           "-s", f"{seed}.{str(fraction).split('.')[1]}", bam, region]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    tlens, read_lens = [], []
    geometry = [Counter(), Counter()]  # two disjoint halves (parity) for stability
    n = 0
    for tlen, span, read_len in read_geometry_records(proc.stdout):
        tlens.append(tlen)
        read_lens.append(read_len)
        geometry[n % 2][(tlen, span)] += 1
        n += 1
    proc.stdout.close()
    if proc.wait() != 0:
        raise RuntimeError(f"samtools view failed for {bam}")
    return tlens, read_lens, geometry


def describe_sample(tlens, read_lens, geometry_halves):
    read_len_hist = Counter(read_lens)
    t_sorted = sorted(tlens)
    n = len(t_sorted)
    full = Counter()
    for h in geometry_halves:
        full.update(h)
    overlap_mates = sum(c for (t, r), c in full.items() if t < 2 * r)
    shorter_than_read = sum(c for (t, r), c in full.items() if t <= r)
    return {
        "n_fragments_sampled": n,
        "read_length_distribution": dict(sorted(read_len_hist.items())),
        "template_size": summarize_values(tlens),
        "n_template_size_gt_1000": sum(1 for t in tlens if t > 1000),
        "fraction_mates_overlapping_(T<2R)": overlap_mates / n if n else None,
        "fraction_template_le_read_span": shorter_than_read / n if n else None,
    }, full


def redundancy_tables(full_geometry, half_geometries):
    tables = {}
    for L in LOCUS_LENGTHS_BP:
        rows = {}
        for gap in SPACINGS_BP:
            rows[str(gap)] = {
                "all_fragments": pooled_redundancy(full_geometry, L, L, gap),
                "half_1": pooled_redundancy(half_geometries[0], L, L, gap),
                "half_2": pooled_redundancy(half_geometries[1], L, L, gap),
            }
        tables[f"locus_length_{L}bp"] = rows
    curve = {gap: pooled_redundancy(full_geometry, PRIMARY_LOCUS_LENGTH_BP,
                                     PRIMARY_LOCUS_LENGTH_BP, gap)
             for gap in CURVE_SPACINGS_BP}
    crossings = {}
    for key in ("same_read_joint", "fragment_joint"):
        crossings[key] = {f"first_gap_bp_at_or_below_{lvl}": spacing_crossing(curve, key, lvl)
                          for lvl in (0.5, 0.25, 0.10, 0.05, 0.01, 0.001)}
    return tables, {"curve_locus_length_10bp": {str(k): v for k, v in curve.items()},
                    "crossings_locus_length_10bp": crossings}


def load_loci(path):
    by_chrom = defaultdict(list)
    with open(path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            by_chrom[row["chrom"]].append((int(row["start"]), int(row["end"])))
    for v in by_chrom.values():
        v.sort()
    return by_chrom


def loci_affected_by_spacing(by_chrom, thresholds=SPACINGS_BP):
    """Counts only -- no locus is removed. `affected` = has a neighboring
    locus with gap (start2 - end1) <= T. `greedy_remaining` = how many loci
    a left-to-right one-per-cluster de-duplication would keep (scenario
    count, not applied)."""
    total = 0
    affected = Counter()
    greedy_kept = Counter()
    for positions in by_chrom.values():
        n = len(positions)
        total += n
        for i in range(n):
            gaps = []
            if i > 0:
                gaps.append(positions[i][0] - positions[i - 1][1])
            if i < n - 1:
                gaps.append(positions[i + 1][0] - positions[i][1])
            g = min(gaps) if gaps else None
            for t in thresholds:
                if g is not None and g <= t:
                    affected[t] += 1
        for t in thresholds:
            last_end = None
            for start, end in positions:
                if last_end is None or start - last_end > t:
                    greedy_kept[t] += 1
                    last_end = end
    return {
        "n_total": total,
        "by_spacing_bp": {
            str(t): {
                "n_loci_with_neighbor_within": affected[t],
                "pct_loci_affected": 100.0 * affected[t] / total if total else None,
                "n_loci_remaining_if_one_per_cluster_scenario_only": greedy_kept[t],
            } for t in thresholds
        },
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--tumor-bam", required=True)
    ap.add_argument("--normal-bam", required=True)
    ap.add_argument("--view-a", required=True)
    ap.add_argument("--view-b", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--region", default="chr17")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--fraction", default="0.05")
    args = ap.parse_args(argv)

    result = {
        "analysis_type": "SEQUENCING_GEOMETRY_EVIDENCE_NOT_A_FILTER",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "parameters": {"region": args.region, "seed": args.seed, "fraction": args.fraction,
                        "min_flank": MIN_FLANK, "locus_lengths_bp": LOCUS_LENGTHS_BP,
                        "spacings_bp": SPACINGS_BP, "min_mapq": 20},
        "bams": {"tumor": args.tumor_bam, "normal": args.normal_bam},
        "variant_or_msi_calls_or_marker_performance_used": False,
        "master_panel_or_views_modified": False,
        "final_panel_selected": False,
        "samples": {},
    }

    for label, bam in (("tumor", args.tumor_bam), ("normal", args.normal_bam)):
        tlens, read_lens, halves = sample_bam_geometry(
            args.samtools, bam, args.region, args.seed, args.fraction)
        desc, full = describe_sample(tlens, read_lens, halves)
        tables, curve = redundancy_tables(full, halves)
        result["samples"][label] = {"geometry": desc, "redundancy_by_spacing": tables, **curve}

    result["loci_affected_by_spacing"] = {
        "view_A_wes_compatible": loci_affected_by_spacing(load_loci(args.view_a)),
        "view_B_wgs": loci_affected_by_spacing(load_loci(args.view_b)),
    }

    t = result["samples"]["tumor"]["redundancy_by_spacing"][f"locus_length_{PRIMARY_LOCUS_LENGTH_BP}bp"]
    n = result["samples"]["normal"]["redundancy_by_spacing"][f"locus_length_{PRIMARY_LOCUS_LENGTH_BP}bp"]
    result["tumor_vs_normal_consistency_locus_length_10bp"] = {
        str(g): {k: (abs(t[str(g)]["all_fragments"][k] - n[str(g)]["all_fragments"][k])
                      if t[str(g)]["all_fragments"][k] is not None
                      and n[str(g)]["all_fragments"][k] is not None else None)
                 for k in ("same_read_joint", "fragment_joint")}
        for g in SPACINGS_BP
    }

    with open(args.out, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(result, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
