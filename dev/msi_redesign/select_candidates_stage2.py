"""
GenRichi MSI Methodology Redesign -- DEVELOPMENT / RUO CANDIDATE SELECTION.
NOT PRODUCTION CODE. NOT A CLINICALLY VALIDATED PANEL. No MSI classification.

Stage 2: informativeness + technical-suitability down-selection, applied to
the already-verified 250 bp candidate representatives of View A and View B
(candidate_msi_panel.v1.*.spacing250.tsv). It answers one question -- "after
applying pre-declared, evidence-based criteria, how many strong candidates
remain naturally?" -- and deliberately does NOT choose or hard-code a final
panel size. Nothing is deleted: every input representative gets a row in the
traceability table with all stage flags, ranking attributes and a reason.

Input hierarchy (pre-declared; all criteria are recorded per locus):
  S0  input: the 250 bp representatives (View A reps, View B reps).
  S1  technical extraction (sequencing geometry): span <= READ_LENGTH -
      2*MIN_FLANK = 90 bp, the mathematical limit for a 100 bp read with 5 bp
      of aligned flank each side (bam_extract.py rule; read length confirmed
      from the real BAMs). Loci above it cannot be measured at all. Expected
      attrition is 0 because the 250 bp stage already excluded them.
      Ranking attribute (no cutoff): read_measurable_fraction =
      (R - 2f - L + 1) / (R + L - 1), the share of read placements
      overlapping a locus of span L that can measure it; documented, derived
      only from the read length/flank rule.
  S2  repeat type: mononucleotide -> MONO_PRIMARY track; STR/non-mono ->
      STR_COMPARISON track. STR loci are NOT discarded: they run through the
      same S3/S4 criteria so they remain available for comparison.
  S3  informative repeat length: span 15-40 bp (user-specified priority band;
      the validated Promega mononucleotide markers, 21-27 bp, lie inside it).
      All five bands are reported separately: <10, 10-14, 15-40, 41-60,
      >60 (61-90 bp).
  S4  mappability: Umap multi-read mean == 1.0 (the existing top tier). Loci
      with NO multi-read data (no coverage in the bedGraph) are labelled
      NOT_ASSESSED and are NOT failed by this stage, so that chrY is not
      removed implicitly; the number of NOT_ASSESSED loci per chromosome is
      reported. Existing Umap tier information is preserved on every row.
  --  RepeatMasker locus-self and flank information are preserved on every
      row and reported before/after; NO new RepeatMasker filter and no new
      biological interpretation is introduced (the upstream exclusion of
      the non-microsatellite class in tier2_filtered.v1 already applies).
  --  Chromosomes: no quotas; chrY is never removed here. Its effect is
      reported (strong counts with and without chrY).

View relationship: every representative is classified A_and_B (shared),
A_only, or B_only from the actual files.

"Strong" = passes S1, S3 and S4. STRONG_MONO_PRIMARY (mono) is the primary
result; STRONG_STR_COMPARISON is reported separately for comparison.
Scenario sizes are feasibility arithmetic only; no size is selected.

No HCC1395 result, no BAM and no production file is read or changed.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone

SCRIPT_VERSION = "candidate_selection_stage2.v1"
READ_LENGTH_BP = 100
MIN_FLANK_BP = 5
MAX_MEASURABLE_SPAN_BP = READ_LENGTH_BP - 2 * MIN_FLANK_BP
STRONG_SPAN_MIN_BP = 15
STRONG_SPAN_MAX_BP = 40
SCENARIO_SIZES = [500, 1000, 1500, 2000, 5000, 10000]

SPAN_BANDS = ["lt10", "10_14", "15_40", "41_60", "gt60"]
STRONG_SUBBANDS = ["15_20", "21_30", "31_40"]
GROUPS = ["viewA", "viewB", "shared_A_and_B", "A_only", "B_only"]

TRACE_HEADER = [
    "marker_id", "representative_id", "source_view_class", "in_viewA_reps", "in_viewB_reps",
    "cluster_id_viewA", "cluster_size_viewA", "n_spacing_excluded_represented_viewA",
    "cluster_id_viewB", "cluster_size_viewB", "n_spacing_excluded_represented_viewB",
    "chrom", "start", "end", "kind", "motif", "ref_repeat_length", "span", "span_band",
    "umap_multiread_mean", "umap_tier", "umap_status", "rmsk_self_category", "rmsk_tier",
    "flank_other_repeat_flag", "read_measurable_fraction", "is_chrY",
    "S1_geometry_span_le_90", "track", "S3_span_15_40", "S4_umap_top_tier_or_not_assessed",
    "final_class", "reason",
]


def span_band(span: int) -> str:
    if span < 10:
        return "lt10"
    if span <= 14:
        return "10_14"
    if span <= 40 and span >= 15:
        return "15_40"
    if span <= 60:
        return "41_60"
    return "gt60"


def strong_subband(span: int):
    if 15 <= span <= 20:
        return "15_20"
    if 21 <= span <= 30:
        return "21_30"
    if 31 <= span <= 40:
        return "31_40"
    return None


def read_measurable_fraction(span: int, read_len=READ_LENGTH_BP, flank=MIN_FLANK_BP) -> float:
    measurable = read_len - 2 * flank - span + 1
    if measurable <= 0:
        return 0.0
    return measurable / (read_len + span - 1)


def umap_status(mean_text: str) -> str:
    if mean_text == "" or mean_text is None:
        return "NOT_ASSESSED"
    return "TIER0_MEAN_1.0" if float(mean_text) >= 1.0 else "BELOW_TOP_TIER"


def view_class(in_a: bool, in_b: bool) -> str:
    if in_a and in_b:
        return "A_and_B"
    return "A_only" if in_a else "B_only"


def evaluate_locus(kind: str, span: int, umap_mean_text: str):
    """Pure per-locus rule evaluation. Returns
    (s1, track, s3, s4, final_class, reason)."""
    s1 = span <= MAX_MEASURABLE_SPAN_BP
    track = "MONO_PRIMARY" if kind == "mononucleotide" else "STR_COMPARISON"
    s3 = STRONG_SPAN_MIN_BP <= span <= STRONG_SPAN_MAX_BP
    us = umap_status(umap_mean_text)
    s4 = us in ("TIER0_MEAN_1.0", "NOT_ASSESSED")
    prefix = "MONO" if track == "MONO_PRIMARY" else "STR_COMPARISON"
    if not s1:
        outcome = "FAIL_S1_SPAN_GT_90"
    elif not s3:
        outcome = f"FAIL_S3_SPAN_BAND_{span_band(span)}"
    elif not s4:
        outcome = "FAIL_S4_UMAP_BELOW_TOP_TIER"
    else:
        outcome = "PASS_ALL_CRITERIA"
    if outcome == "PASS_ALL_CRITERIA":
        final = "STRONG_MONO_PRIMARY" if track == "MONO_PRIMARY" else "STRONG_STR_COMPARISON"
    else:
        final = "NOT_STRONG_MONO" if track == "MONO_PRIMARY" else "NOT_STRONG_STR"
    return s1, track, s3, s4, final, f"{prefix}:{outcome}"


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_panel(path):
    rows = {}
    with open(path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows[r["marker_id"]] = r
    return rows


def load_represented_counts(trace_path):
    counts = Counter()
    with open(trace_path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["status"] == "EXCLUDED_BY_SPACING":
                counts[r["representative_id"]] += 1
    return counts


def build_records(panel_a, panel_b, represented_a, represented_b):
    records = []
    for mid in sorted(set(panel_a) | set(panel_b)):
        ra, rb = panel_a.get(mid), panel_b.get(mid)
        base = ra or rb
        span = int(base["end"]) - int(base["start"])
        s1, track, s3, s4, final, reason = evaluate_locus(base["kind"], span, base["umap_multiread_mean"])
        records.append({
            "marker_id": mid, "representative_id": mid,
            "source_view_class": view_class(ra is not None, rb is not None),
            "in_viewA_reps": ra is not None, "in_viewB_reps": rb is not None,
            "cluster_id_viewA": ra["cluster_id"] if ra else "", "cluster_size_viewA": ra["cluster_size"] if ra else "",
            "n_spacing_excluded_represented_viewA": represented_a.get(mid, 0) if ra else "",
            "cluster_id_viewB": rb["cluster_id"] if rb else "", "cluster_size_viewB": rb["cluster_size"] if rb else "",
            "n_spacing_excluded_represented_viewB": represented_b.get(mid, 0) if rb else "",
            "chrom": base["chrom"], "start": int(base["start"]), "end": int(base["end"]),
            "kind": base["kind"], "motif": base["motif"], "ref_repeat_length": base["ref_repeat_length"],
            "span": span, "span_band": span_band(span),
            "umap_multiread_mean": base["umap_multiread_mean"], "umap_tier": base["umap_tier"],
            "umap_status": umap_status(base["umap_multiread_mean"]),
            "rmsk_self_category": base["rmsk_self_category"], "rmsk_tier": base["rmsk_tier"],
            "flank_other_repeat_flag": base["flank_other_repeat_flag"],
            "read_measurable_fraction": round(read_measurable_fraction(span), 4),
            "is_chrY": base["chrom"] == "chrY",
            "S1_geometry_span_le_90": s1, "track": track, "S3_span_15_40": s3,
            "S4_umap_top_tier_or_not_assessed": s4, "final_class": final, "reason": reason,
        })
    return records


def in_group(rec, group):
    a, b = rec["in_viewA_reps"], rec["in_viewB_reps"]
    return {"viewA": a, "viewB": b, "shared_A_and_B": a and b, "A_only": a and not b, "B_only": b and not a}[group]


def composition(recs):
    recs = list(recs)
    return {
        "n": len(recs),
        "kind": dict(Counter(r["kind"] for r in recs)),
        "span_band": {b: sum(1 for r in recs if r["span_band"] == b) for b in SPAN_BANDS},
        "umap_status": dict(Counter(r["umap_status"] for r in recs)),
        "umap_tier": dict(sorted(Counter(str(r["umap_tier"]) for r in recs).items())),
        "rmsk_self_category": dict(Counter(r["rmsk_self_category"] for r in recs)),
        "flank_other_repeat_flag": dict(Counter(str(r["flank_other_repeat_flag"]) for r in recs)),
        "motif_length": dict(sorted(Counter(len(r["motif"]) for r in recs).items())),
        "chromosome": dict(Counter(r["chrom"] for r in recs)),
        "strong_subbands_15_40": {b: sum(1 for r in recs if strong_subband(r["span"]) == b) for b in STRONG_SUBBANDS},
    }


def stage_counts(recs):
    """Sequential attrition, reported per track (nothing is dropped from the
    traceability table)."""
    out = {"S0_input_250bp_representatives": len(recs),
           "S1_pass_geometry_span_le_90": sum(1 for r in recs if r["S1_geometry_span_le_90"])}
    for track in ("MONO_PRIMARY", "STR_COMPARISON"):
        t = [r for r in recs if r["track"] == track and r["S1_geometry_span_le_90"]]
        s3 = [r for r in t if r["S3_span_15_40"]]
        s4 = [r for r in s3 if r["S4_umap_top_tier_or_not_assessed"]]
        out[track] = {
            "S2_entering_track": len(t),
            "S3_pass_span_15_40": len(s3),
            "S4_pass_umap_top_tier_or_not_assessed_=_STRONG": len(s4),
            "S4_of_which_umap_NOT_ASSESSED": sum(1 for r in s4 if r["umap_status"] == "NOT_ASSESSED"),
            "S4_of_which_chrY": sum(1 for r in s4 if r["is_chrY"]),
            "S4_strong_excluding_chrY": sum(1 for r in s4 if not r["is_chrY"]),
            "S3_chrY": sum(1 for r in s3 if r["is_chrY"]),
            "track_chrY_entering": sum(1 for r in t if r["is_chrY"]),
        }
    return out


def summarize(records):
    summary = {}
    for g in GROUPS:
        recs = [r for r in records if in_group(r, g)]
        strong_mono = [r for r in recs if r["final_class"] == "STRONG_MONO_PRIMARY"]
        strong_str = [r for r in recs if r["final_class"] == "STRONG_STR_COMPARISON"]
        summary[g] = {
            "stage_attrition": stage_counts(recs),
            "composition_S0_input": composition(recs),
            "composition_strong_mono_primary": composition(strong_mono),
            "composition_strong_str_comparison": composition(strong_str),
            "chrY_impact": {
                "chrY_in_S0": sum(1 for r in recs if r["is_chrY"]),
                "strong_mono_with_chrY": len(strong_mono),
                "strong_mono_without_chrY": sum(1 for r in strong_mono if not r["is_chrY"]),
                "strong_mono_chrY": sum(1 for r in strong_mono if r["is_chrY"]),
                "strong_str_with_chrY": len(strong_str),
                "strong_str_without_chrY": sum(1 for r in strong_str if not r["is_chrY"]),
            },
            "umap_not_assessed_by_chromosome": dict(Counter(r["chrom"] for r in recs if r["umap_status"] == "NOT_ASSESSED")),
            "scenario_size_feasibility_arithmetic_only": {
                str(n): {
                    "strong_mono_available": len(strong_mono),
                    "size_le_available": n <= len(strong_mono),
                    "pct_of_available": round(100.0 * n / len(strong_mono), 3) if strong_mono else None,
                } for n in SCENARIO_SIZES
            },
        }
    return summary


def write_tsv(path, header, rows):
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--panel-a", required=True)
    ap.add_argument("--panel-b", required=True)
    ap.add_argument("--trace-a", required=True)
    ap.add_argument("--trace-b", required=True)
    ap.add_argument("--prior-manifest", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)

    with open(args.prior_manifest, encoding="utf-8") as fh:
        prior = json.load(fh)
    expected = {
        args.panel_a: prior["outputs"]["viewA_wes_compatible"]["selected_panel"]["sha256"],
        args.trace_a: prior["outputs"]["viewA_wes_compatible"]["traceability_all_loci"]["sha256"],
        args.panel_b: prior["outputs"]["viewB_wgs"]["selected_panel"]["sha256"],
        args.trace_b: prior["outputs"]["viewB_wgs"]["traceability_all_loci"]["sha256"],
    }
    input_hashes = {}
    for path, want in expected.items():
        got = sha256_file(path)
        if got != want:
            raise SystemExit(f"input {path} does not match the 250 bp manifest hash ({got} != {want})")
        input_hashes[path] = got

    panel_a, panel_b = load_panel(args.panel_a), load_panel(args.panel_b)
    records = build_records(panel_a, panel_b, load_represented_counts(args.trace_a),
                            load_represented_counts(args.trace_b))
    summary = summarize(records)

    trace_path = os.path.join(args.out_dir, "candidate_selection_stage2.v1.traceability.tsv")
    write_tsv(trace_path, TRACE_HEADER, ([r[h] for h in TRACE_HEADER] for r in records))
    strong_path = os.path.join(args.out_dir, "candidate_selection_stage2.v1.strong_mono_primary.tsv")
    write_tsv(strong_path, TRACE_HEADER,
              ([r[h] for h in TRACE_HEADER] for r in records if r["final_class"] == "STRONG_MONO_PRIMARY"))
    comp_path = os.path.join(args.out_dir, "candidate_selection_stage2.v1.strong_str_comparison.tsv")
    write_tsv(comp_path, TRACE_HEADER,
              ([r[h] for h in TRACE_HEADER] for r in records if r["final_class"] == "STRONG_STR_COMPARISON"))

    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "stage": "informativeness + technical-suitability down-selection (development / RUO)",
        "status": "DEVELOPMENT_CANDIDATE_NOT_CLINICALLY_VALIDATED",
        "final_panel_size_fixed": False,
        "script": {"version": SCRIPT_VERSION, "path": os.path.abspath(__file__),
                   "sha256": sha256_file(os.path.abspath(__file__)), "python": platform.python_version()},
        "inputs": {p: {"sha256": h} for p, h in input_hashes.items()},
        "prior_manifest": {"path": args.prior_manifest, "sha256": sha256_file(args.prior_manifest)},
        "parameters": {
            "read_length_bp": READ_LENGTH_BP, "min_flank_bp": MIN_FLANK_BP,
            "max_measurable_span_bp": MAX_MEASURABLE_SPAN_BP,
            "strong_span_bp": [STRONG_SPAN_MIN_BP, STRONG_SPAN_MAX_BP],
            "span_bands": {"lt10": "<10", "10_14": "10-14", "15_40": "15-40", "41_60": "41-60", "gt60": ">60 (61-90)"},
            "umap_rule": "multi-read mean == 1.0, or NOT_ASSESSED when no multi-read data (no implicit chrY removal)",
            "scenario_sizes_arithmetic_only": SCENARIO_SIZES,
        },
        "rules": {
            "S1": "span <= 90 (100 bp read, 5 bp flank each side)",
            "S2": "mononucleotide -> MONO_PRIMARY; STR -> STR_COMPARISON (retained, not discarded)",
            "S3": "span 15-40 bp (user-specified priority band)",
            "S4": "Umap multi-read mean == 1.0 or NOT_ASSESSED",
            "RepeatMasker": "information preserved; no new filter or interpretation",
            "chromosomes": "no quotas; chrY never removed here; effect reported",
        },
        "hcc1395_results_used": False, "classification_performed": False,
        "production_changes": False, "previous_250bp_outputs_modified": False,
        "n_union_representatives": len(records),
        "view_relationship": dict(Counter(r["source_view_class"] for r in records)),
        "summary": summary,
        "outputs": {
            "traceability_all_representatives": {"path": trace_path, "sha256": sha256_file(trace_path)},
            "strong_mono_primary": {"path": strong_path, "sha256": sha256_file(strong_path)},
            "strong_str_comparison": {"path": comp_path, "sha256": sha256_file(comp_path)},
        },
    }
    manifest_path = os.path.join(args.out_dir, "candidate_selection_stage2.v1.manifest.json")
    with open(manifest_path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print(f"wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
