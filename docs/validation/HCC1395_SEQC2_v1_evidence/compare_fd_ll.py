#!/usr/bin/env python3
"""Cross-replicate comparison of GenRichi Phase 1 PASS calls: WES_FD vs WES_LL (same
HCC1395/HCC1395BL biological pair, different SEQC2 sequencing site).

Does NOT modify evaluate.py, its masks, or any pipeline/config file. Standalone,
read-only over existing run outputs in ~/seqc2_phase1_test/run/results/.

Shared-evaluable mask = Phase-1-panel ∩ SEQC2-High-Confidence-Regions ∩
{FD tumor depth >=20 AND FD normal depth >=20 AND LL tumor depth >=20 AND LL normal depth >=20}.
Depth: samtools depth -a -Q20 -q20 on each run's own final BAMs (same method as evaluate.py).
Matching: exact (chrom,pos,ref,alt) after `bcftools norm -m -any -f ref` (same rule as evaluate.py).
"""
import json
import subprocess
from pathlib import Path

T = Path.home() / "seqc2_phase1_test"
EV = T / "eval"
RUN = T / "run"
REF = Path.home() / "genrichi/reference_db/ref/hg38.fa"
CONDA = Path.home() / "genrichi/.snakemake/conda"
BCF = CONDA / "5a6a3d0b5d257e1cf1abb4ba0c8e5a66_/bin/bcftools"
SAMT = CONDA / "b3e920885114bc7c28eb9e53d02501ec_/bin/samtools"
SAMPLES = {"FD": "SEQC2_HCC1395_FD1", "LL": "SEQC2_HCC1395_LL1"}


def sh(cmd):
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True, executable="/bin/bash").stdout


def depth_table(sample):
    tb = RUN / f"results/{sample}/tumor/align/{sample}.tumor.final.bam"
    nb = RUN / f"results/{sample}/normal/align/{sample}.normal.final.bam"
    dep = {}
    for line in sh(f"{SAMT} depth -a -Q 20 -q 20 -b {EV}/panel_hc.bed {tb} {nb}").splitlines():
        c, p, t, n = line.split("\t")
        dep[(c, int(p))] = (int(t), int(n))
    return dep


def norm_records(vcf, regions):
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


def main():
    print("Computing per-run depth over panel∩HC (%s)..." % (EV / "panel_hc.bed"))
    dep = {k: depth_table(s) for k, s in SAMPLES.items()}

    shared = []
    for (c, p) in dep["FD"]:
        ft, fn = dep["FD"][(c, p)]
        lt, ln = dep["LL"].get((c, p), (0, 0))
        if ft >= 20 and fn >= 20 and lt >= 20 and ln >= 20:
            shared.append((c, p - 1, p))
    shared_bed = EV / "shared_fd_ll_mask.bed"
    with open(shared_bed, "w") as fh:
        for c, s, e in shared:
            fh.write(f"{c}\t{s}\t{e}\n")
    shared_merged = EV / "shared_fd_ll_mask.merged.bed"
    sh(f"sort -k1,1 -k2,2n {shared_bed} | bedtools merge -i - > {shared_merged}")
    n_int, n_bp = 0, 0
    for line in open(shared_merged):
        f = line.split()
        n_int += 1
        n_bp += int(f[2]) - int(f[1])

    truth = {}
    for t in ("sSNV", "sINDEL"):
        truth.update(norm_records(T / f"truth/high-confidence_{t}_in_HC_regions_v1.2.1.vcf.gz", shared_merged))

    calls = {}
    for k, s in SAMPLES.items():
        calls[k] = norm_records(RUN / f"results/{s}/snv/{s}.pass.vcf.gz", shared_merged)

    both = set(calls["FD"]) & set(calls["LL"])
    only_fd = set(calls["FD"]) - set(calls["LL"])
    only_ll = set(calls["LL"]) - set(calls["FD"])
    all_called = set(calls["FD"]) | set(calls["LL"])
    tp = all_called & set(truth)
    fp_any = all_called - set(truth)

    def fmt(keys):
        return sorted(f"{c}:{p} {r}>{a} ({kind((c,p,r,a))})" for c, p, r, a in keys)

    report = {
        "shared_mask": {"intervals": n_int, "bases": n_bp,
                        "definition": "panel ∩ SEQC2-HC ∩ {FD & LL, tumor>=20 & normal>=20}"},
        "truth_variants_in_shared_mask": {k: sum(1 for r in truth if kind(r) == k) for k in ("SNV", "indel", "other")},
        "calls_in_shared_mask": {kk: {k: sum(1 for r in calls[kk] if kind(r) == k) for k in ("SNV", "indel", "other")}
                                  for kk in SAMPLES},
        "agreement": {
            "called_in_both_runs": fmt(both),
            "called_in_FD_only": fmt(only_fd),
            "called_in_LL_only": fmt(only_ll),
        },
        "vs_truth_within_shared_mask": {
            "TP_total_sites_matching_truth": fmt(tp),
            "TP_count": len(tp),
            "unmatched_to_truth_any_run": fmt(fp_any),
            "unmatched_count": len(fp_any),
            "unmatched_called_by": {f"{c}:{p} {r}>{a}": [k for k in SAMPLES if (c, p, r, a) in calls[k]]
                                     for (c, p, r, a) in fp_any},
        },
    }

    # --- specific investigation: FD's unmatched chr16:23635069 in the LL run ---
    site = ("chr16", 23635069, "G", "A")
    ll_dep = dep["LL"].get((site[0], site[1]), (None, None))
    fd_dep = dep["FD"].get((site[0], site[1]), (None, None))
    d = f"{RUN}/results/{SAMPLES['LL']}/snv"
    def view_site(vcf):
        out = sh(f"{BCF} view -H {vcf} {site[0]}:{site[1]}-{site[1]}").strip()
        return out if out else "absent"
    chr16_probe = {
        "site": f"{site[0]}:{site[1]} {site[2]}>{site[3]}",
        "FD_depth_tumor_normal": fd_dep,
        "LL_depth_tumor_normal": ll_dep,
        "LL_in_pass_vcf": view_site(f"{d}/{SAMPLES['LL']}.pass.vcf.gz"),
        "LL_in_filtered_vcf": view_site(f"{d}/{SAMPLES['LL']}.filtered.vcf.gz"),
        "LL_in_raw_mutect2_vcf": view_site(f"{d}/{SAMPLES['LL']}.mutect2.vcf.gz"),
    }
    report["chr16_23635069_FD_unmatched_call_investigated_in_LL"] = chr16_probe

    (EV / "fd_vs_ll_shared_mask_comparison.json").write_text(json.dumps(report, indent=2, default=list))
    print(json.dumps(report, indent=2, default=list))


main()
