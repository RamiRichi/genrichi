"""
INDEPENDENT verification of empirical_evidence_stage6a.v1 outputs.

Imports NO Stage 1-6 implementation code (not even bam_extract / msi_locus_model:
the provisional parameters are re-read from their SOURCE TEXT). It uses a
different read-access path from the implementation (one `samtools view
<bam> <region>` query per locus, versus one multi-region BED pass) and a
different CIGAR formulation (aligned-block interval overlap, versus an
operation walk). From the raw inputs it independently

  1. pins the Stage-5 manifest and re-hashes the whole upstream chain,
  2. re-hashes every input (registry, reference index, every BAM and index,
     source files) and every output,
  3. re-reads the provisional parameters from the prototype source text,
  4. re-checks reference / contig / coordinate integrity (BAM @SQ vs .fai;
     every union locus vs the FASTA read DIRECTLY by .fai offsets),
  5. recomputes the duplicate-content classes of the BAMs,
  6. recomputes every locus x sample-file row from the reads and compares
     each cell,
  7. recomputes the per-locus summary, the evidence matrix (dimension
     statuses, NOT_ASSESSED vs INSUFFICIENT, no zero-substitution) and the
     scenario comparison (memberships rebuilt from the Stage-3 rank columns),
  8. recomputes the parameter-sensitivity tables, and
  9. checks the circularity guard (no strat_* column among evidence columns),
     gating flags (Stage 6B not executed) and safety flags.

Usage (from the repository root):
    python3 dev/msi_redesign/verify_empirical_evidence_stage6a.py --manifest <manifest.json>
Exit code 0 = all checks passed.
"""

import csv
import hashlib
import json
import math
import os
import re
import subprocess
import sys
from collections import Counter

PINNED_STAGE5 = "9a5ecebf89e230cedf94e9bf25cce643dbb816829d0608190aa5ac4241442bef"
HERE = os.path.dirname(os.path.abspath(__file__))

LS_COLS = ["file_id", "biological_unit_id", "role", "order_id", "dataset_origin", "unit_label",
           "marker_id", "chrom", "start", "end", "ref_repeat_length",
           "n_overlap_alignments", "n_excluded_flag", "n_excluded_mapq", "n_pass_filters", "n_insufficient_flank",
           "n_informative", "n_informative_fragments", "callable_fraction", "informative_fraction_of_pass",
           "median_mapq", "n_informative_fwd", "n_informative_rev", "fwd_fraction", "median_rel_offset",
           "status", "length_hist", "modal_length", "modal_fraction", "ref_allele_fraction",
           "mean_delta_from_ref", "sd_length"]
STRAT_COLS = ["strat_source_view_class", "strat_span_subband", "strat_umap_status",
              "strat_rmsk_self_category", "strat_read_measurable_fraction", "strat_rank_viewB_no_shared_key"]
STRAT_SRC = ["source_view_class", "span_subband", "umap_status", "rmsk_self_category",
             "read_measurable_fraction", "rank_viewB_no_shared_key"]
ROLE_COLS = ["n_files", "n_assessed", "informative_min", "informative_median", "informative_max",
             "callable_fraction_median", "fwd_fraction_median", "median_mapq_median", "rel_offset_median",
             "modal_length_values", "modal_fraction_median", "ref_allele_fraction_median"]
LTS_COLS = (["marker_id", "chrom", "start", "end", "ref_repeat_length", "motif", "data_scope",
             "n_analysis_files", "n_missing_files", "missing_fraction"]
            + ["tumor_" + c for c in ROLE_COLS] + ["normal_" + c for c in ROLE_COLS]
            + ["tn_pairs_assessed", "tn_tvd_min", "tn_tvd_median", "tn_tvd_max", "tumor_rerun_modal_length_concordant"]
            + STRAT_COLS)
EVID = ["technical_evidence", "coverage_evidence", "stability_evidence", "MSI_separation_evidence",
        "reproducibility_evidence", "missingness", "sample_count", "ground_truth_count",
        "ground_truth_MSS_units", "ground_truth_MSI_H_units", "file_count"]
SUPP = ["tech_median_mapq_median", "tech_fwd_fraction_median", "tech_rel_offset_median",
        "cov_tumor_informative_median", "cov_normal_informative_median", "cov_callable_fraction_median",
        "stab_tn_tvd_median", "stab_modal_fraction_median"]
EM_COLS = ["candidate", "chrom", "start", "end", "data_scope"] + EVID + SUPP + STRAT_COLS
SC_COLS = ["scenario_id", "view", "strategy", "size", "core", "n_loci", "n_in_scope", "fraction_in_scope",
           "n_out_of_scope_NOT_ASSESSED", "technical_ASSESSED", "technical_INSUFFICIENT",
           "coverage_ASSESSED", "coverage_INSUFFICIENT", "stability_ASSESSED", "stability_INSUFFICIENT",
           "MSI_separation_NOT_ASSESSED", "reproducibility_INSUFFICIENT", "reproducibility_NOT_ASSESSED",
           "missingness_median_in_scope", "cov_tumor_informative_median_p10", "cov_tumor_informative_median_p50",
           "cov_tumor_informative_median_p90", "in_scope_shared", "in_scope_viewB_or_A_only",
           "core_in_scope", "extension_in_scope", "evidence_scope_caveat"]
SENS_MAPQ, SENS_FLANK, SENS_DEPTH = [0, 20, 30, 40], [1, 5, 10, 20], [1, 5, 10, 20, 50]


def sha(path):
    h = hashlib.sha256()
    with open(norm(path), "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def norm(path):
    """Stage 1-5 manifests record Windows backslashes (immutable); normalise at read time only."""
    return path.replace("\\", "/")


def load(path):
    with open(norm(path), encoding="utf-8") as f:
        return json.load(f)


def tsv(path):
    with open(norm(path), encoding="utf-8", newline="") as f:
        rd = csv.reader(f, delimiter="\t")
        head = next(rd)
        return head, [dict(zip(head, r)) for r in rd]


def source_const(fname, name):
    with open(os.path.join(HERE, fname), encoding="utf-8") as fh:
        text = fh.read()
    m = re.search(r"^%s\s*=\s*(\d+)" % re.escape(name), text, re.M)
    return int(m.group(1)) if m else None


def med(v):
    v = sorted(v)
    if not v:
        return None
    n = len(v)
    return round(float(v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2), 6)


def quant(v, q):
    v = sorted(v)
    if not v:
        return None
    pos = (len(v) - 1) * q
    lo = math.floor(pos)
    hi = min(lo + 1, len(v) - 1)
    return round(v[lo] + (v[hi] - v[lo]) * (pos - lo), 6)


def same(cell, want, tol=2e-6):
    """Compare a TSV cell to a recomputed value (None -> 'NA')."""
    if want is None:
        return cell == "NA"
    if isinstance(want, bool):
        return cell == ("True" if want else "False")
    if isinstance(want, (int, float)):
        try:
            return abs(float(cell) - float(want)) <= tol
        except ValueError:
            return False
    return cell == str(want)


# --------------------------------------------------------------------------- independent read analysis
def blocks_of(cigar):
    """Aligned blocks (ref_start, ref_end, query_start), insertions (next_ref_pos, n), ref_span, read_len."""
    aligned, ins, rp, qp, read_len = [], [], 0, 0, 0
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        n = int(n)
        if op in "M=X":
            aligned.append((rp, rp + n, qp))
            rp += n
            qp += n
            read_len += n
        elif op == "I":
            ins.append((rp, n))
            qp += n
            read_len += n
        elif op == "S":
            qp += n
            read_len += n
        elif op in "DN":
            rp += n
    return aligned, ins, rp, read_len


def analyse(fields, start, end):
    qname, flag, pos, mapq, cigar = fields[0], int(fields[1]), int(fields[3]), int(fields[4]), fields[5]
    if flag & 4 or cigar in ("*", ""):
        return None
    aligned, ins, span, read_len = blocks_of(cigar)
    r0 = pos - 1
    if not (r0 < end and r0 + span > start):
        return None
    fb = sum(max(0, min(re_ + r0, start) - (rs + r0)) for rs, re_, _ in aligned)
    fa = sum(max(0, (re_ + r0) - max(rs + r0, end)) for rs, re_, _ in aligned)
    obs = sum(max(0, min(re_ + r0, end) - max(rs + r0, start)) for rs, re_, _ in aligned)
    obs += sum(n for rp, n in ins if start <= rp + r0 < end)
    off = next((qs + max(rs + r0, start) - (rs + r0) for rs, re_, qs in aligned if re_ + r0 > start), None)
    return {"qname": qname, "excl": bool(flag & (0x100 | 0x400 | 0x800)), "mapq": mapq, "rev": bool(flag & 16),
            "fb": fb, "fa": fa, "obs": obs, "off": off, "rl": read_len}


def region_reads(samtools, bam, chrom, start, end):
    out = subprocess.run([samtools, "view", bam, f"{chrom}:{start + 1}-{end}"], capture_output=True, text=True, encoding="utf-8", errors="replace")
    if out.returncode != 0:
        raise SystemExit(f"samtools failed: {out.stderr[:200]}")
    excl, reads = 0, []
    for line in out.stdout.splitlines():
        rec = analyse(line.split("\t", 6), start, end)
        if rec is None:
            continue
        if rec["excl"]:
            excl += 1
        else:
            reads.append(rec)
    return reads, excl


def informative(r, mq, fl):
    return r["mapq"] >= mq and r["fb"] >= fl and r["fa"] >= fl


def hist_of(reads, mq, fl):
    h = {}
    for r in reads:
        if informative(r, mq, fl):
            h[r["obs"]] = h.get(r["obs"], 0) + 1
    return h


def allele(h, ref_len):
    n = sum(h.values())
    top = max(h.values())
    modal = min(k for k, v in h.items() if v == top)
    mean = sum(k * v for k, v in h.items()) / n
    ex2 = sum(k * k * v for k, v in h.items()) / n
    return {"modal": modal, "modal_fraction": h[modal] / n, "ref_frac": h.get(ref_len, 0) / n,
            "delta": mean - ref_len, "sd": math.sqrt(max(0.0, ex2 - mean * mean))}


def locus_sample(reads, excl, ref_len, mq, fl, depth):
    low = sum(1 for r in reads if r["mapq"] < mq)
    inf = [r for r in reads if informative(r, mq, fl)]
    n, n_ov = len(inf), excl + len(reads)
    h = hist_of(reads, mq, fl)
    fwd = sum(1 for r in inf if not r["rev"])
    row = {"n_overlap_alignments": n_ov, "n_excluded_flag": excl, "n_excluded_mapq": low,
           "n_pass_filters": len(reads) - low, "n_insufficient_flank": len(reads) - low - n, "n_informative": n,
           "n_informative_fragments": len({r["qname"] for r in inf}),
           "callable_fraction": n / n_ov if n_ov else None,
           "informative_fraction_of_pass": n / (len(reads) - low) if len(reads) - low else None,
           "median_mapq": med([r["mapq"] for r in reads]), "n_informative_fwd": fwd, "n_informative_rev": n - fwd,
           "fwd_fraction": fwd / n if n else None,
           "median_rel_offset": med([r["off"] / r["rl"] for r in inf if r["off"] is not None and r["rl"]]),
           "length_hist": ";".join(f"{k}:{h[k]}" for k in sorted(h)) if n else None,
           "status": "MISSING_NO_INFORMATIVE_READS" if n == 0 else
                     ("INSUFFICIENT_BELOW_PROVISIONAL_MIN_DEPTH" if n < depth else "ASSESSED"),
           "modal_length": None, "modal_fraction": None, "ref_allele_fraction": None,
           "mean_delta_from_ref": None, "sd_length": None}
    if row["status"] == "ASSESSED":
        a = allele(h, ref_len)
        row.update(modal_length=a["modal"], modal_fraction=a["modal_fraction"], ref_allele_fraction=a["ref_frac"],
                   mean_delta_from_ref=a["delta"], sd_length=a["sd"])
    return row, h


def tv(a, b):
    na, nb = sum(a.values()), sum(b.values())
    return 0.5 * sum(abs(a.get(k, 0) / na - b.get(k, 0) / nb) for k in set(a) | set(b))


# --------------------------------------------------------------------------- FASTA (direct, by .fai offsets)
class Fasta:
    def __init__(self, path):
        self.path = path
        self.fai = {}
        with open(path + ".fai", encoding="utf-8") as f:
            for line in f:
                p = line.split("\t")
                self.fai[p[0]] = tuple(int(x) for x in p[1:5])
        self.fh = open(path, "rb")   # closed with the object (short-lived verifier process)

    def close(self):
        self.fh.close()

    def fetch(self, chrom, start, end):
        length, offset, lb, lw = self.fai[chrom]
        first, last = start, end - 1
        b0 = offset + (first // lb) * lw + first % lb
        b1 = offset + (last // lb) * lw + last % lb
        self.fh.seek(b0)
        return self.fh.read(b1 - b0 + 1).replace(b"\n", b"").replace(b"\r", b"").decode().upper()


# --------------------------------------------------------------------------- verification
def verify(manifest_path, skip_upstream=False):
    fails = []

    def check(cond, msg):
        if not cond:
            fails.append(msg)

    m = load(manifest_path)
    samtools = os.environ.get("STAGE6A_SAMTOOLS", "samtools")

    # 1. pin + chain --------------------------------------------------------
    if not skip_upstream:
        check(sha(m["stage5_manifest_path"]) == PINNED_STAGE5, "stage-5 manifest differs from the pinned hash")
        check(m["expected_stage5_manifest_sha256"] == PINNED_STAGE5, "manifest records a different pin")
        m5 = load(m["stage5_manifest_path"])
        for name, info in m5["outputs"].items():
            check(sha(info["path"]) == info["sha256"], f"stage-5 output changed: {name}")
        m4 = load(m5["stage4_manifest_path"])
        check(sha(m5["stage4_manifest_path"]) == m5["verified_upstream_hashes"]["stage4_manifest"], "stage-4 manifest changed")
        for name, info in m4["outputs"].items():
            check(sha(info["path"]) == info["sha256"], f"stage-4 output changed: {name}")
        m3 = load(m4["stage3_manifest_path"])
        check(sha(m4["stage3_manifest_path"]) == m4["verified_upstream_hashes"]["stage3_manifest"], "stage-3 manifest changed")
        for name, info in m3["outputs"].items():
            check(sha(info["path"]) == info["sha256"], f"stage-3 output changed: {name}")
        m2 = load(m3["stage2_manifest_path"])
        check(sha(m3["stage2_manifest_path"]) == m3["inputs_verified_hashes"]["stage2_manifest"], "stage-2 manifest changed")
        for name, info in m2["outputs"].items():
            check(sha(info["path"]) == info["sha256"], f"stage-2 output changed: {name}")
        check(sha(m2["prior_manifest"]["path"]) == m2["prior_manifest"]["sha256"], "250 bp manifest changed")
        prior = load(m2["prior_manifest"]["path"])
        for view, outs in prior["outputs"].items():
            for key, info in outs.items():
                check(sha(info["path"]) == info["sha256"], f"250 bp output changed: {view}/{key}")
        union_path = m3["outputs"]["ranked_union"]["path"]
        check(m["verified_upstream_hashes"].get("stage3_ranked_union") == sha(union_path), "recorded union hash differs")
        p5 = m5["parameters"]
        sizes, cores, totals = p5["sizes"], p5["hybrid_cores"], p5["hybrid_totals"]
    else:
        union_path = m["inputs"]["stage3_ranked_union"]["path"]
        sizes, cores, totals = m["parameters"]["sizes"], m["parameters"]["hybrid_cores"], m["parameters"]["hybrid_totals"]

    # 2. inputs and outputs -------------------------------------------------
    inp = m["inputs"]
    for k in ("registry", "reference_fai", "stage1_scan_manifest", "stage3_ranked_union"):
        check(sha(inp[k]["path"]) == inp[k]["sha256"], f"input changed: {k}")
    for fid, b in inp["bams"].items():
        check(sha(b["path"]) == b["bam_sha256"], f"BAM changed: {fid}")
        check(sha(b["path"] + ".bai") == b["bai_sha256"], f"BAM index changed: {fid}")
    for name, info in m["outputs"].items():
        check(sha(info["path"]) == info["sha256"], f"output hash mismatch: {name}")
    for fname, h in m["source_files_sha256"].items():
        check(sha(os.path.join(HERE, fname)) == h, f"source file changed since the run: {fname}")

    # 3. provisional parameters from source text ----------------------------
    mq, fl, depth = (source_const("bam_extract.py", "MIN_MAPQ"), source_const("bam_extract.py", "MIN_FLANK"),
                     source_const("msi_locus_model.py", "MIN_LOCUS_DEPTH"))
    P = m["parameters"]
    check((mq, fl, depth) == (P["MIN_MAPQ"]["value"], P["MIN_FLANK"]["value"], P["MIN_LOCUS_DEPTH"]["value"]),
          "provisional parameters differ from the prototype source text")
    for k in ("MIN_MAPQ", "MIN_FLANK", "MIN_LOCUS_DEPTH"):
        check(P[k]["validated"] is False and P[k]["status"] == "PROVISIONAL_DEVELOPMENT_PARAMETER", f"{k} not labelled provisional")
    check(P["not_optimised_or_calibrated_on_HCC1395"] is True, "calibration flag missing")

    # flags ------------------------------------------------------------------
    for flag in ("final_panel_size_fixed", "panel_size_selected", "membership_files_written", "msi_classification_performed",
                 "msi_discrimination_performed", "stage_6b_executed", "clinical_claims", "composite_score_computed",
                 "production_changes", "upstream_outputs_modified", "downloads_performed"):
        check(m.get(flag) is False, f"flag not False: {flag}")
    check(m["stage_6B"]["executed"] is False and m["stage_6B"]["status"].startswith("SPECIFIED_IN_DESIGN_GATED"), "6B gate flag wrong")
    check(len(m["claims_not_made"]) >= 8, "claims_not_made list missing")
    check(not set(m["evidence_columns"]) & set(m["stratification_columns_descriptive_only"]), "evidence and strat columns overlap")
    check(list(m["stratification_columns_descriptive_only"]) == STRAT_COLS, "strat column list differs")
    check(all(c.startswith("strat_") for c in STRAT_COLS), "strat prefix violated")
    check(not any(c.startswith("strat_") for c in EM_COLS[:5] + EVID + SUPP), "an evidence column carries the strat prefix")

    # data ---------------------------------------------------------------------
    reg = load(inp["registry"]["path"])
    units, scope = reg["biological_units"], reg["analysis_scope_contigs"]
    fasta = Fasta(reg["reference"]["fasta"])
    fai_names = [(n, fasta.fai[n][0]) for n in fasta.fai]
    uhead, urows = tsv(union_path)
    for r in urows:
        r["_shared"] = r["source_view_class"] == "A_and_B"
        r["_inA"], r["_inB"] = r["in_viewA_strong"] == "True", r["in_viewB_strong"] == "True"

    # 4. reference / contigs / coordinates ------------------------------------
    for f in reg["files"]:
        hdr = subprocess.run([samtools, "view", "-H", f["bam"]], capture_output=True, text=True, encoding="utf-8", errors="replace").stdout
        sq = []
        for line in hdr.splitlines():
            if line.startswith("@SQ"):
                d = dict(x.split(":", 1) for x in line.split("\t")[1:])
                sq.append((d["SN"], int(d["LN"])))
        check(sq == fai_names, f"{f['file_id']}: BAM @SQ differs from the reference index")
        idx = subprocess.run([samtools, "idxstats", f["bam"]], capture_output=True, text=True,
                             encoding="utf-8", errors="replace").stdout.splitlines()
        counts = {l.split("\t")[0]: int(l.split("\t")[2]) for l in idx if l}
        prim = [c for c in counts if re.fullmatch(r"chr([0-9]+|X|Y|M)", c)]
        check(all(counts.get(c, 0) > 0 for c in scope), f"{f['file_id']}: a declared scope contig holds no reads")
        top = max(prim, key=lambda c: (counts[c], c))
        check(top in scope, f"{f['file_id']}: the contig with most reads ({top}) is not in the declared scope")
        off = sum(counts[c] for c in prim if c not in scope)
        got = m["reference_checks"][f["file_id"]]
        check(got["mapped_reads_off_scope_primary_contigs"] == off, f"{f['file_id']}: off-scope read count differs")
        check(got["mapped_reads_total"] == sum(counts.values()), f"{f['file_id']}: total mapped reads differ")
        check(abs(got["off_scope_read_fraction_of_primary"] - off / (off + sum(counts[c] for c in scope))) < 2e-6,
              f"{f['file_id']}: off-scope fraction differs")
    check(reg.get("off_scope_policy") == "NOT_ASSESSED_residual_reads_ignored", "off_scope_policy missing in registry")
    bad = 0
    for r in urows:
        s, e = int(r["start"]), int(r["end"])
        if r["chrom"] not in fasta.fai or e - s != int(r["ref_repeat_length"]) or not (0 <= s < e <= fasta.fai[r["chrom"]][0]):
            bad += 1
            continue
        seq, mo = fasta.fetch(r["chrom"], s, e), r["motif"].upper()
        if len(seq) != e - s or len(seq) % len(mo) or seq != mo * (len(seq) // len(mo)):
            bad += 1
    fasta.close()
    check(bad == 0, f"{bad} loci fail the direct FASTA coordinate/repeat check")
    check(m["locus_coordinate_check"]["n_loci_checked"] == len(urows) and m["locus_coordinate_check"]["n_mismatch"] == 0,
          "manifest coordinate check inconsistent")

    # 5. duplicate classes ------------------------------------------------------
    sigs = {}
    for f in reg["files"]:
        cmd = f"'{samtools}' view '{f['bam']}' | awk -F'\\t' '{{print $1\"\\t\"$2\"\\t\"$4\"\\t\"$6}}' | sha256sum"
        out = subprocess.run(["bash", "-c", "set -o pipefail; " + cmd], capture_output=True, text=True, encoding="utf-8", errors="replace")
        sigs[f["file_id"]] = out.stdout.split()[0] if out.returncode == 0 else None
    for f in reg["files"]:
        rec = m["duplicate_checks"][f["file_id"]]
        check(sigs[f["file_id"]] == rec["content_signature"], f"{f['file_id']}: content signature differs")
    first, analysis = {}, []
    for f in reg["files"]:
        key = (f["biological_unit_id"], f["role"], sigs[f["file_id"]])
        dup_of = first.get(key)
        first.setdefault(key, f["file_id"])
        check(m["duplicate_checks"][f["file_id"]]["exact_duplicate_of"] == dup_of, f"{f['file_id']}: duplicate assignment differs")
        if dup_of is None:
            analysis.append(f)
    check(m["analysis_file_ids"] == [f["file_id"] for f in analysis], "analysis file list differs")
    by_sig_units = {}
    for f in reg["files"]:
        by_sig_units.setdefault(sigs[f["file_id"]], set()).add(f["biological_unit_id"])
    check(all(len(v) == 1 for v in by_sig_units.values()), "identical content under different biological units")

    # 6. locus x sample rows -----------------------------------------------------
    inscope = [r for r in urows if r["chrom"] in scope]
    ls_head, ls_tsv = tsv(m["outputs"]["locus_sample_technical"]["path"])
    check(ls_head == LS_COLS, "locus_sample_technical header differs")
    tsv_index = {(r["file_id"], r["marker_id"]): r for r in ls_tsv}
    check(len(ls_tsv) == len(analysis) * len(inscope) == len(tsv_index), "locus_sample row count differs")
    ls = {}
    reads_store = {}
    for f in analysis:
        for r in inscope:
            s, e, L = int(r["start"]), int(r["end"]), int(r["ref_repeat_length"])
            reads, excl = region_reads(samtools, f["bam"], r["chrom"], s, e)
            row, h = locus_sample(reads, excl, L, mq, fl, depth)
            row.update(file_id=f["file_id"], role=f["role"], biological_unit_id=f["biological_unit_id"],
                       order_id=f["order_id"], marker_id=r["marker_id"], hist=h)
            ls[(f["file_id"], r["marker_id"])] = row
            reads_store[(f["file_id"], r["marker_id"])] = (reads, L)
            cell = tsv_index.get((f["file_id"], r["marker_id"]))
            if cell is None:
                check(False, f"missing locus_sample row {f['file_id']}/{r['marker_id']}")
                continue
            for c in LS_COLS[11:]:
                if not same(cell[c], row[c]):
                    check(False, f"locus_sample {f['file_id']}/{r['marker_id']}: {c} = {cell[c]!r}, recomputed {row[c]!r}")
                    break
            for c, want in (("dataset_origin", units[f["biological_unit_id"]]["dataset_origin"]),
                            ("unit_label", units[f["biological_unit_id"]]["label"]), ("chrom", r["chrom"]),
                            ("start", s), ("end", e), ("ref_repeat_length", L)):
                check(same(cell[c], want), f"locus_sample {f['file_id']}/{r['marker_id']}: {c} differs")

    # 7. summary + evidence matrix ------------------------------------------------
    lts_head, lts_tsv = tsv(m["outputs"]["locus_technical_summary"]["path"])
    check(lts_head == LTS_COLS, "locus_technical_summary header differs")
    lts_index = {r["marker_id"]: r for r in lts_tsv}
    check(len(lts_tsv) == len(inscope), "locus_technical_summary row count differs")
    em_head, em_tsv = tsv(m["outputs"]["evidence_matrix"]["path"])
    check(em_head == EM_COLS, "evidence_matrix header differs")
    check(len(em_tsv) == len(urows), "evidence_matrix must have one row per union locus")
    em_rows = {}
    for u, cell in zip(urows, em_tsv):
        check(cell["candidate"] == u["marker_id"], "evidence_matrix candidate order differs from the union")
        for sc, src in zip(STRAT_COLS, STRAT_SRC):
            want_cell = u[src] if u[src] != "" else "NA"     # convention: undefined -> NA (e.g. no View-B rank for A-only loci)
            check(cell[sc] == want_cell, f"evidence_matrix {u['marker_id']}: {sc} differs")
        mid = u["marker_id"]
        if u["chrom"] not in scope:
            for c in EVID:
                check(cell[c] == "NOT_ASSESSED", f"out-of-scope {mid}: {c} must be NOT_ASSESSED, found {cell[c]!r}")
            for c in SUPP:
                check(cell[c] == "NA", f"out-of-scope {mid}: {c} must be NA")
            check(cell["data_scope"] == "OUT_OF_SCOPE_NO_LOCAL_DATA", f"{mid}: data_scope wrong")
            em_rows[mid] = {"scope": False}
            continue
        rows = [ls[(f["file_id"], mid)] for f in analysis]
        ass = [r for r in rows if r["status"] == "ASSESSED"]
        for role in ("tumor", "normal"):
            sel = [r for r in rows if r["role"] == role]
            sa = [r for r in sel if r["status"] == "ASSESSED"]
            inf = [r["n_informative"] for r in sel]
            want = {"n_files": len(sel), "n_assessed": len(sa),
                    "informative_min": min(inf) if inf else None, "informative_median": med(inf),
                    "informative_max": max(inf) if inf else None,
                    "callable_fraction_median": med([r["callable_fraction"] for r in sel if r["callable_fraction"] is not None]),
                    "fwd_fraction_median": med([r["fwd_fraction"] for r in sa]),
                    "median_mapq_median": med([r["median_mapq"] for r in sel if r["median_mapq"] is not None]),
                    "rel_offset_median": med([r["median_rel_offset"] for r in sa if r["median_rel_offset"] is not None]),
                    "modal_length_values": ";".join(str(x) for x in sorted({r["modal_length"] for r in sa})) or None,
                    "modal_fraction_median": med([r["modal_fraction"] for r in sa]),
                    "ref_allele_fraction_median": med([r["ref_allele_fraction"] for r in sa])}
            for c, v in want.items():
                check(same(lts_index[mid][f"{role}_{c}"], v), f"summary {mid}: {role}_{c} = {lts_index[mid][f'{role}_{c}']!r}, recomputed {v!r}")
        orders = {}
        for r in rows:
            orders.setdefault((r["biological_unit_id"], r["order_id"]), {})[r["role"]] = r
        tvds = [round(tv(p["tumor"]["hist"], p["normal"]["hist"]), 6) for p in orders.values()
                if "tumor" in p and "normal" in p and p["tumor"]["status"] == "ASSESSED" and p["normal"]["status"] == "ASSESSED"]
        tmodal = [r["modal_length"] for r in rows if r["role"] == "tumor" and r["status"] == "ASSESSED"]
        n_missing = sum(1 for r in rows if r["status"] == "MISSING_NO_INFORMATIVE_READS")
        for c, v in (("n_analysis_files", len(rows)), ("n_missing_files", n_missing),
                     ("missing_fraction", n_missing / len(rows)), ("tn_pairs_assessed", len(tvds)),
                     ("tn_tvd_min", min(tvds) if tvds else None), ("tn_tvd_median", med(tvds)),
                     ("tn_tvd_max", max(tvds) if tvds else None),
                     ("tumor_rerun_modal_length_concordant", (len(set(tmodal)) == 1) if len(tmodal) >= 2 else None)):
            check(same(lts_index[mid][c], v), f"summary {mid}: {c} differs")
        cov_ok = (any(r["role"] == "tumor" for r in ass) and any(r["role"] == "normal" for r in ass))
        stab_ok = bool(tvds)
        au = {r["biological_unit_id"] for r in ass}
        lab = Counter(units[u_]["label"] for u_ in au)
        n_units_files = len({r["biological_unit_id"] for r in rows})
        st = lambda ok: "ASSESSED" if ok else "INSUFFICIENT"
        want = {"technical_evidence": st(bool(ass)), "coverage_evidence": st(cov_ok), "stability_evidence": st(stab_ok),
                "MSI_separation_evidence": "NOT_ASSESSED",
                "reproducibility_evidence": "INSUFFICIENT" if n_units_files < 2 else st(len(au) >= 2),
                "missingness": n_missing / len(rows), "sample_count": len(au),
                "ground_truth_count": lab.get("MSS", 0) + lab.get("MSI-H", 0), "ground_truth_MSS_units": lab.get("MSS", 0),
                "ground_truth_MSI_H_units": lab.get("MSI-H", 0), "file_count": len(rows),
                "tech_median_mapq_median": med([r["median_mapq"] for r in ass]) if ass else None,
                "tech_fwd_fraction_median": med([r["fwd_fraction"] for r in ass]) if ass else None,
                "tech_rel_offset_median": med([r["median_rel_offset"] for r in ass if r["median_rel_offset"] is not None]) if ass else None,
                "cov_tumor_informative_median": med([r["n_informative"] for r in rows if r["role"] == "tumor"]) if cov_ok else None,
                "cov_normal_informative_median": med([r["n_informative"] for r in rows if r["role"] == "normal"]) if cov_ok else None,
                "cov_callable_fraction_median": med([r["callable_fraction"] for r in ass]) if cov_ok else None,
                "stab_tn_tvd_median": med(tvds) if stab_ok else None,
                "stab_modal_fraction_median": med([r["modal_fraction"] for r in ass]) if stab_ok else None}
        check(cell["data_scope"] == "IN_SCOPE_LOCAL_BAM", f"{mid}: data_scope wrong")
        for c, v in want.items():
            check(same(cell[c], v), f"evidence_matrix {mid}: {c} = {cell[c]!r}, recomputed {v!r}")
            if isinstance(v, str) and v == "INSUFFICIENT":
                pass
        # INSUFFICIENT dimensions must not carry numeric support values
        if not cov_ok:
            check(all(cell[c] == "NA" for c in ("cov_tumor_informative_median", "cov_normal_informative_median", "cov_callable_fraction_median")),
                  f"{mid}: coverage numeric values present although coverage_evidence is INSUFFICIENT")
        if not stab_ok:
            check(cell["stab_tn_tvd_median"] == "NA" and cell["stab_modal_fraction_median"] == "NA",
                  f"{mid}: stability numeric values present although stability_evidence is INSUFFICIENT")
        em_rows[mid] = {"scope": True, **want}
    check(all(c["MSI_separation_evidence"] == "NOT_ASSESSED" for c in em_tsv), "MSI_separation_evidence must be NOT_ASSESSED everywhere")

    # 8. scenario comparison --------------------------------------------------------
    sc_head, sc_tsv = tsv(m["outputs"]["scenario_comparison"]["path"])
    check(sc_head == SC_COLS, "scenario_comparison header differs")
    by = {r["scenario_id"]: r for r in sc_tsv}

    def order(flag, col):
        return sorted((r for r in urows if r[flag]), key=lambda r: int(r[col]))

    aq, as_, bq, bs = (order("_inA", "rank_viewA_no_shared_key"), order("_inA", "rank_viewA"),
                       order("_inB", "rank_viewB_no_shared_key"), order("_inB", "rank_viewB"))
    core_all = sorted((r for r in urows if r["_shared"]), key=lambda r: int(r["rank_viewB_no_shared_key"]))
    expect = []
    for view, q, s in (("viewA", aq, as_), ("viewB", bq, bs)):
        for size in sizes:
            expect.append((f"{view}_quality_first_{size}", view, "quality_first", size, None, q[:size], None))
            expect.append((f"{view}_shared_first_{size}", view, "shared_first", size, None, s[:size], None))
    for k in cores:
        for total in totals:
            if total <= k:
                continue
            core = core_all[:k]
            ids = {r["marker_id"] for r in core}
            members = core + [r for r in bq if r["marker_id"] not in ids][: total - k]
            check(len({r["marker_id"] for r in members}) == total, f"hybrid c{k}/T{total}: not unique")
            expect.append((f"viewB_hybrid_core{k}_total{total}", "viewB", "hybrid", total, k, members, core))
    check(len(sc_tsv) == len(expect) == m["n_scenarios"], "scenario count differs")
    for sid, view, strat, size, core_k, members, core in expect:
        cell = by.get(sid)
        if cell is None:
            check(False, f"missing scenario {sid}")
            continue
        ev = [em_rows[r["marker_id"]] for r in members]
        ins = [(r, e) for r, e in zip(members, ev) if e["scope"]]
        covs = [e["cov_tumor_informative_median"] for _, e in ins if e.get("cov_tumor_informative_median") is not None]
        miss = [e["missingness"] for _, e in ins]
        core_ids = {r["marker_id"] for r in core} if core else None
        ins_ids = {r["marker_id"] for r, _ in ins}
        want = {"view": view, "strategy": strat, "size": size, "core": core_k, "n_loci": len(members), "n_in_scope": len(ins),
                "fraction_in_scope": len(ins) / len(members), "n_out_of_scope_NOT_ASSESSED": len(members) - len(ins),
                "technical_ASSESSED": sum(1 for _, e in ins if e["technical_evidence"] == "ASSESSED"),
                "technical_INSUFFICIENT": sum(1 for _, e in ins if e["technical_evidence"] == "INSUFFICIENT"),
                "coverage_ASSESSED": sum(1 for _, e in ins if e["coverage_evidence"] == "ASSESSED"),
                "coverage_INSUFFICIENT": sum(1 for _, e in ins if e["coverage_evidence"] == "INSUFFICIENT"),
                "stability_ASSESSED": sum(1 for _, e in ins if e["stability_evidence"] == "ASSESSED"),
                "stability_INSUFFICIENT": sum(1 for _, e in ins if e["stability_evidence"] == "INSUFFICIENT"),
                "MSI_separation_NOT_ASSESSED": len(members),
                "reproducibility_INSUFFICIENT": sum(1 for _, e in ins if e["reproducibility_evidence"] == "INSUFFICIENT"),
                "reproducibility_NOT_ASSESSED": len(members) - len(ins),
                "missingness_median_in_scope": med(miss) if miss else None,
                "cov_tumor_informative_median_p10": quant(covs, 0.10), "cov_tumor_informative_median_p50": quant(covs, 0.50),
                "cov_tumor_informative_median_p90": quant(covs, 0.90),
                "in_scope_shared": sum(1 for r, _ in ins if r["_shared"]),
                "in_scope_viewB_or_A_only": sum(1 for r, _ in ins if not r["_shared"]),
                "core_in_scope": sum(1 for r in core if r["marker_id"] in ins_ids) if core_ids is not None else None,
                "extension_in_scope": sum(1 for r, _ in ins if r["marker_id"] not in core_ids) if core_ids is not None else None}
        for c, v in want.items():
            check(same(cell[c], v), f"scenario {sid}: {c} = {cell[c]!r}, recomputed {v!r}")
        check(cell["evidence_scope_caveat"].startswith("CHR17_ONLY"), f"scenario {sid}: caveat missing")

    # 9. sensitivity tables (descriptive) ----------------------------------------------
    for f in analysis:
        rec = m["parameter_sensitivity"][f["file_id"]]
        loci = {mid: v for (fid, mid), v in reads_store.items() if fid == f["file_id"]}

        def n_meet(a, b, d):
            return sum(1 for reads, _ in loci.values() if sum(1 for r in reads if informative(r, a, b)) >= d)
        check(rec["n_loci_with_at_least_depth_at_default_mapq_flank"] == {str(d): n_meet(mq, fl, d) for d in SENS_DEPTH},
              f"{f['file_id']}: depth sensitivity table differs")
        check(rec["n_loci_meeting_provisional_depth_by_mapq_flank"] ==
              {f"mapq{a}_flank{b}": n_meet(a, b, depth) for a in SENS_MAPQ for b in SENS_FLANK},
              f"{f['file_id']}: mapq/flank sensitivity table differs")
        base = {}
        for mid, (reads, L) in loci.items():
            h = hist_of(reads, mq, fl)
            if sum(h.values()) >= depth:
                base[mid] = allele(h, L)["modal"]
        exp = {}
        for a in SENS_MAPQ:
            for b in SENS_FLANK:
                both = eq = 0
                for mid, mod in base.items():
                    h = hist_of(loci[mid][0], a, b)
                    if sum(h.values()) >= depth:
                        both += 1
                        eq += allele(h, loci[mid][1])["modal"] == mod
                exp[f"mapq{a}_flank{b}"] = {"loci_assessed_in_both": both, "modal_length_unchanged": eq}
        check(rec["modal_length_agreement_with_default_by_mapq_flank"] == exp, f"{f['file_id']}: modal agreement table differs")
    return fails


def main(argv):
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    args = ap.parse_args(argv)
    fails = verify(args.manifest)
    if fails:
        print("VERIFICATION FAILED")
        for f in fails[:200]:
            print(" -", f)
        if len(fails) > 200:
            print(f" ... and {len(fails) - 200} more")
        return 1
    print("VERIFICATION PASSED: independent read extraction, evidence statuses, scenario counts, sensitivity tables and hashes agree")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
