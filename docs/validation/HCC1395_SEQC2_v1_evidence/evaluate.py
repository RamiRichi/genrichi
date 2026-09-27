#!/usr/bin/env python3
"""Limited check of GenRichi Phase 1 PASS calls against SEQC2 somatic truth (v1.2.1).

DESIGN, fixed before any pipeline result was seen
-------------------------------------------------
Truth-evaluable universe  U = Phase-1 panel (merged) ∩ SEQC2 High-Confidence_Regions_v1.2.bed.
                          Anything outside U is NOT evaluable: a truth variant there does not exist for us,
                          a call there is not called FP.
PRIMARY stratum    P = U ∩ {tumor depth >= 20 and normal depth >= 20}   (matches the pipeline's own
                                                                          calling.filter.min_depth=20; user-directed
                                                                          primary mask, 2026-09-27)
ADDITIONAL strata 10 = U ∩ {tumor >= 10 and normal >= 10}   (more permissive, for context)
                  30 = U ∩ {tumor >= 30 and normal >= 20}
                  U  = no coverage requirement (transparency only)
Depth: samtools depth -a -Q 20 -q 20 on the pipeline's final BAMs (duplicates/secondary/QC-fail excluded by default).
A truth variant at a base below the stratum's coverage threshold is NOT counted FN in that stratum: it is listed as
"not evaluable (coverage)". A call outside the stratum's regions is not counted FP.
Matching: exact (chrom,pos,ref,alt) after `bcftools norm -m -any -f ref`. SNV = len(ref)==len(alt)==1;
indel = len(ref)!=len(alt); everything else (MNV/complex) is reported apart and excluded from both.
Calls: results/<sample>/snv/<sample>.pass.vcf.gz.  Reasons for FN are read from filtered.vcf.gz / mutect2.vcf.gz.
"""
import json
import math
import subprocess
import sys
from pathlib import Path

T = Path.home() / "seqc2_phase1_test"
EV = T / "eval"
SAMPLE = sys.argv[1]
RUN = T / "run"
REF = Path.home() / "genrichi/reference_db/ref/hg38.fa"
CONDA = Path.home() / "genrichi/.snakemake/conda"
BCF = CONDA / "5a6a3d0b5d257e1cf1abb4ba0c8e5a66_/bin/bcftools"
SAMT = CONDA / "b3e920885114bc7c28eb9e53d02501ec_/bin/samtools"
STRATA = [("P_primary_T20_N20", 20, 20), ("ADD_T10_N10", 10, 10), ("ADD_T30_N20", 30, 20)]


def sh(cmd):
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True, executable="/bin/bash").stdout


def bed_bases(path):
    n = b = 0
    for line in open(path):
        f = line.split()
        n += 1
        b += int(f[2]) - int(f[1])
    return n, b


def norm_records(vcf, regions):
    # -T/--targets-file (not -R/--regions-file) is required here: -R needs a random-access
    # index, which a piped/streamed bgzip on stdin does not have; -T streams and filters.
    out = sh(f"{BCF} norm -m -any -f {REF} {vcf} 2>/dev/null | {BCF} view -H -T {regions} -")
    recs = {}
    for line in out.splitlines():
        f = line.split("\t")
        recs[(f[0], int(f[1]), f[3], f[4])] = f
    return recs


def kind(key):
    _, _, ref, alt = key
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) != len(alt):
        return "indel"
    return "other"


def wilson(k, n, z=1.96):
    if n == 0:
        return None
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return [round((c - h) / d, 3), round((c + h) / d, 3)]


def depth_table():
    tb = RUN / f"results/{SAMPLE}/tumor/align/{SAMPLE}.tumor.final.bam"
    nb = RUN / f"results/{SAMPLE}/normal/align/{SAMPLE}.normal.final.bam"
    dep = {}
    for line in sh(f"{SAMT} depth -a -Q 20 -q 20 -b {EV}/panel_hc.bed {tb} {nb}").splitlines():
        c, p, t, n = line.split("\t")
        dep[(c, int(p))] = (int(t), int(n))
    return dep, {"tumor_bam": str(tb), "normal_bam": str(nb)}


def stratum_bed(name, dep, tmin, nmin):
    path = EV / f"stratum_{name}.bed"
    with open(path, "w") as fh:
        for (c, p), (t, n) in sorted(dep.items()):
            if t >= tmin and n >= nmin:
                fh.write(f"{c}\t{p - 1}\t{p}\n")
    merged = EV / f"stratum_{name}.merged.bed"
    sh(f"sort -k1,1 -k2,2n {path} | bedtools merge -i - > {merged}")
    return merged


def main():
    pass_vcf = RUN / f"results/{SAMPLE}/snv/{SAMPLE}.pass.vcf.gz"
    filt_vcf = RUN / f"results/{SAMPLE}/snv/{SAMPLE}.filtered.vcf.gz"
    raw_vcf = RUN / f"results/{SAMPLE}/snv/{SAMPLE}.mutect2.vcf.gz"
    truth_vcfs = {t: T / f"truth/high-confidence_{t}_in_HC_regions_v1.2.1.vcf.gz" for t in ("sSNV", "sINDEL")}
    dep, bams = depth_table()

    def truth_in(regions):
        d = {}
        for t in truth_vcfs:
            d.update(norm_records(truth_vcfs[t], regions))
        return d

    truth_panel = truth_in(EV / "panel.merged.bed")
    truth_U = truth_in(EV / "panel_hc.bed")
    calls_panel = norm_records(pass_vcf, EV / "panel.merged.bed")
    calls_U = norm_records(pass_vcf, EV / "panel_hc.bed")
    filt_U = norm_records(filt_vcf, EV / "panel_hc.bed")
    raw_U = norm_records(raw_vcf, EV / "panel_hc.bed")

    n_int, n_bp = bed_bases(EV / "panel.merged.bed")
    u_int, u_bp = bed_bases(EV / "panel_hc.bed")
    out = {
        "sample": SAMPLE, "pass_vcf": str(pass_vcf), **bams,
        "design": "see module docstring of evaluate.py (fixed before results)",
        "universe": {
            "panel_intervals": n_int, "panel_bases": n_bp,
            "panel_and_HC_intervals": u_int, "panel_and_HC_bases": u_bp,
            "panel_bases_outside_HC": n_bp - u_bp,
            "note": "chrX panel regions are outside SEQC2 HC regions (autosomes only) and are not evaluable",
        },
        "truth_variants": {
            k: {"in_panel": sum(1 for r in truth_panel if kind(r) == k),
                "in_panel_and_HC": sum(1 for r in truth_U if kind(r) == k)}
            for k in ("SNV", "indel", "other")},
        "pass_calls": {
            k: {"in_panel": sum(1 for r in calls_panel if kind(r) == k),
                "in_panel_and_HC": sum(1 for r in calls_U if kind(r) == k),
                "in_panel_outside_HC_not_evaluable": sum(1 for r in calls_panel if kind(r) == k)
                - sum(1 for r in calls_U if kind(r) == k)}
            for k in ("SNV", "indel", "other")},
        "strata": [],
    }
    per_variant = []
    fp_rows = []
    for name, tmin, nmin in STRATA + [("U_no_coverage_filter", 0, 0)]:
        regions = stratum_bed(name, dep, tmin, nmin) if tmin else EV / "panel_hc.bed"
        r_int, r_bp = bed_bases(regions)
        truth = truth_in(regions)
        calls = norm_records(pass_vcf, regions)
        res = {}
        for k in ("SNV", "indel", "other"):
            tk = {r for r in truth if kind(r) == k}
            ck = {r for r in calls if kind(r) == k}
            tp, fp, fn = tk & ck, ck - tk, tk - ck
            prec = len(tp) / (len(tp) + len(fp)) if (len(tp) + len(fp)) else None
            rec = len(tp) / (len(tp) + len(fn)) if (len(tp) + len(fn)) else None
            f1 = None if prec is None or rec is None else (0.0 if prec + rec == 0 else round(2 * prec * rec / (prec + rec), 3))
            res[k] = {"truth_variants": len(tk), "TP": len(tp), "FP": len(fp), "FN": len(fn),
                      "precision": None if prec is None else round(prec, 3),
                      "recall": None if rec is None else round(rec, 3),
                      "recall_wilson95_illustrative": wilson(len(tp), len(tk)),
                      "F1": f1,
                      "FP_per_Mb_evaluable": round(len(fp) / (r_bp / 1e6), 2) if r_bp else None,
                      "TP_sites": sorted(tp), "FP_sites": sorted(fp), "FN_sites": sorted(fn)}
            if name == "P_primary_T20_N20":
                for site in sorted(fp):
                    fp_rows.append(("FP", k, site, calls[site][6], "\t".join(calls[site][8:])))
        not_eval = sorted(set(r for r in truth_U if r not in truth))
        out["strata"].append({
            "name": name, "min_tumor_depth": tmin, "min_normal_depth": nmin,
            "intervals": r_int, "bases": r_bp,
            "truth_variants_in_panel_and_HC": {k: sum(1 for r in truth_U if kind(r) == k) for k in ("SNV", "indel", "other")},
            "truth_variants_evaluable_here": {k: sum(1 for r in truth if kind(r) == k) for k in ("SNV", "indel", "other")},
            "truth_not_evaluable_due_to_coverage": [
                {"site": f"{c}:{p} {ref}>{alt}", "type": kind((c, p, ref, alt)),
                 "tumor_depth": dep.get((c, p), (None, None))[0], "normal_depth": dep.get((c, p), (None, None))[1]}
                for (c, p, ref, alt) in not_eval],
            "results": res})
    # per-truth-variant table with diagnosis of misses
    for (c, p, ref, alt) in sorted(truth_U):
        key = (c, p, ref, alt)
        td, nd = dep.get((c, p), (None, None))
        per_variant.append({
            "site": f"{c}:{p} {ref}>{alt}", "type": kind(key), "tumor_depth": td, "normal_depth": nd,
            "in_pass_vcf": key in calls_U,
            "FILTER_in_filtered_vcf": (filt_U[key][6] if key in filt_U else "absent"),
            "in_raw_mutect2_vcf": key in raw_U})
    out["per_truth_variant"] = per_variant
    (EV / f"{SAMPLE}.metrics.json").write_text(json.dumps(out, indent=2, default=list))
    with open(EV / f"{SAMPLE}.primary_FP_calls.tsv", "w") as fh:
        fh.write("class\ttype\tsite\tFILTER\tFORMAT_and_samples\n")
        for cls, k, site, flt, rest in fp_rows:
            fh.write(f"{cls}\t{k}\t{site[0]}:{site[1]} {site[2]}>{site[3]}\t{flt}\t{rest}\n")
    print(json.dumps({k: out[k] for k in ("universe", "truth_variants", "pass_calls")}, indent=2))
    for s in out["strata"]:
        print("\n==", s["name"], "bases", s["bases"], "intervals", s["intervals"])
        print(" truth in panel∩HC:", s["truth_variants_in_panel_and_HC"], " evaluable here:", s["truth_variants_evaluable_here"])
        for k in ("SNV", "indel"):
            r = s["results"][k]
            print(f"  {k}: TP={r['TP']} FP={r['FP']} FN={r['FN']} precision={r['precision']} recall={r['recall']} F1={r['F1']} FP/Mb={r['FP_per_Mb_evaluable']}")
        print("  not evaluable (coverage):", s["truth_not_evaluable_due_to_coverage"])
    print("\nper truth variant:")
    for v in per_variant:
        print(" ", v)


main()
