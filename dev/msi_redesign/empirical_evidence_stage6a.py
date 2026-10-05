"""
GenRichi MSI Methodology Redesign -- DEVELOPMENT / RUO EMPIRICAL EVIDENCE, STAGE 6A.
NOT PRODUCTION CODE. NOT A CLINICALLY VALIDATED PANEL. No MSI classification.

Stage 6A = the LOCAL, READ-ONLY part of Stage 6:
  * technical evidence (coverage, callable fraction, missingness, mapping
    quality, strand, read-position) per locus x sample-file;
  * MSS-only allele-length descriptives for the single labelled biological
    unit available locally (HCC1395, MSS per an EXTERNAL annotation);
  * limited to the contigs present in the local BAMs (chr17), and every
    output says so.
There is NO MSI-vs-MSS discrimination, no panel selection, no panel size, no
composite candidate score, no threshold and no clinical claim. Stage 6B
(discrimination) is specified in the design document but is GATED and is NOT
executed here: it needs independent MSI-H and MSS labelled data that do not
exist locally.

Two distinct non-value states are never merged and never turned into 0 or
into "stable":
  NOT_ASSESSED : the analysis was not / could not be performed (locus outside
                 the data scope, or Stage 6B not run).
  INSUFFICIENT : the analysis was attempted but the data do not reach the
                 declared minimum.

PROVISIONAL DEVELOPMENT PARAMETERS. MIN_MAPQ, MIN_FLANK and MIN_LOCUS_DEPTH
are reused from the early prototype (bam_extract.py / msi_locus_model.py)
ONLY as explicitly provisional development parameters. They are not
validated, not calibrated and not tuned on HCC1395. A sensitivity table
(descriptive counts only) is written to the manifest.

CIRCULARITY GUARD. Stage-3 ranking attributes (|span-25|, Umap, RepeatMasker,
read-measurable fraction, uniqueness, Stage-3 ranks, tie hash) are NEVER used
to compute any empirical value. They appear only in separate `strat_*`
columns. The only Stage-3 columns read by the empirical code are the
coordinate/reference-allele representation (chrom, start, end, motif,
ref_repeat_length) and the scenario memberships (for the comparison table).

Immutable inputs: the Stage-5 manifest is pinned by hash, and everything
upstream of it is re-hashed before any work. BAMs are only read.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import re
import statistics
import subprocess
import sys
import tempfile
from bisect import bisect_left
from collections import Counter
from datetime import datetime, timezone

import tradeoff_pareto_stage5 as t5           # scenario definitions + hash-chain verification (unchanged)
from bam_extract import MIN_MAPQ, MIN_FLANK    # provisional prototype parameters (unchanged module)
from msi_locus_model import MIN_LOCUS_DEPTH    # provisional prototype parameter (unchanged module)

SCRIPT_VERSION = "empirical_evidence_stage6a.v1"
OUT_PREFIX = "empirical_evidence_stage6a.v1"
EXPECTED_STAGE5_MANIFEST_SHA256 = "9a5ecebf89e230cedf94e9bf25cce643dbb816829d0608190aa5ac4241442bef"

FLAG_UNMAPPED, FLAG_REVERSE, FLAG_SECONDARY, FLAG_DUPLICATE, FLAG_SUPPLEMENTARY = 0x4, 0x10, 0x100, 0x400, 0x800
FLAG_EXCLUDE = FLAG_SECONDARY | FLAG_DUPLICATE | FLAG_SUPPLEMENTARY   # unmapped (0x4) records are skipped outright

ST_ASSESSED = "ASSESSED"
ST_INSUFFICIENT = "INSUFFICIENT"
ST_NOT_ASSESSED = "NOT_ASSESSED"
FS_ASSESSED = "ASSESSED"
FS_INSUFFICIENT = "INSUFFICIENT_BELOW_PROVISIONAL_MIN_DEPTH"
FS_MISSING = "MISSING_NO_INFORMATIVE_READS"

SENS_MAPQ = [0, 20, 30, 40]
SENS_FLANK = [1, 5, 10, 20]
SENS_DEPTH = [1, 5, 10, 20, 50]

_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")

CLAIMS_NOT_MADE = [
    "no MSI-H / MSS classification and no clinical sensitivity or specificity",
    "no MSI-vs-MSS discrimination evidence (Stage 6B is gated and NOT executed)",
    "no clinical threshold; the provisional parameters are not validated thresholds",
    "no final panel, no panel membership, no panel size",
    "no statement that any strategy (quality-first / shared-first / hybrid) is best or clinically better",
    "no conversion of technical ranking into biological or clinical superiority",
    "HCC1395 is one biological unit (MSS, external label): it does not validate anything and nothing generalises from n=1",
    "results cover chr17 only; loci elsewhere are NOT_ASSESSED, not 'absent' or 'stable'",
    "INSUFFICIENT and NOT_ASSESSED are never zero and never evidence of stability",
]


# --------------------------------------------------------------------------- basics
def sha256_file(path):
    return t5.sha256_file(path)


def fmt(x):
    """Deterministic value formatting; None -> NA."""
    if x is None:
        return "NA"
    if isinstance(x, bool):
        return "True" if x else "False"
    if isinstance(x, float):
        return format(round(x, 6), ".6f")
    return str(x)


def median_or_none(values):
    values = list(values)
    return round(float(statistics.median(values)), 6) if values else None


def write_tsv(path, header, rows):
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for r in rows:
            w.writerow([fmt(v) for v in r])


def run(cmd, **kw):
    proc = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace", **kw)
    if proc.returncode != 0:
        raise SystemExit(f"command failed ({' '.join(cmd[:3])} ...): {proc.stderr.strip()[:400]}")
    return proc.stdout


# --------------------------------------------------------------------------- pure read analysis
def parse_cigar(cigar):
    if cigar in ("", "*"):
        return []
    return [(int(n), op) for n, op in _CIGAR_RE.findall(cigar)]


def analyse_alignment(flag, pos, mapq, cigar, start, end):
    """Pure. One SAM record vs one locus [start, end) (0-based half-open).

    Returns None if the record is unmapped (0x4), has no CIGAR, or its aligned
    reference span does not overlap the locus. Otherwise a dict:
      excluded_flag : secondary/duplicate/supplementary (0x100|0x400|0x800)
      mapq, reverse
      flank_before / flank_after : aligned (M,=,X) reference bases strictly before start / at-or-after end
      observed_len : read bases aligned inside [start,end) + inserted bases whose next reference
                     position lies in [start,end)
      offset       : query index (0-based, soft clips included) of the aligned base at the smallest
                     reference position >= start; read_len = sum of M,I,S,=,X
    """
    if flag & FLAG_UNMAPPED:
        return None
    ops = parse_cigar(cigar)
    if not ops:
        return None
    ref_pos = pos - 1
    ref_span = sum(n for n, op in ops if op in "MDN=X")
    if not (ref_pos < end and ref_pos + ref_span > start):
        return None
    read_len = sum(n for n, op in ops if op in "MIS=X")
    qpos = 0
    flank_before = flank_after = observed = 0
    offset = None
    for n, op in ops:
        if op in "M=X":
            s, e = ref_pos, ref_pos + n
            flank_before += max(0, min(e, start) - s)
            flank_after += max(0, e - max(s, end))
            observed += max(0, min(e, end) - max(s, start))
            if offset is None and e > start:
                offset = qpos + (max(s, start) - s)
            ref_pos += n
            qpos += n
        elif op == "I":
            if start <= ref_pos < end:
                observed += n
            qpos += n
        elif op == "S":
            qpos += n
        elif op in "DN":
            ref_pos += n
        # H, P: no reference or query contribution
    return {"excluded_flag": bool(flag & FLAG_EXCLUDE), "mapq": mapq, "reverse": bool(flag & FLAG_REVERSE),
            "flank_before": flank_before, "flank_after": flank_after, "observed_len": observed,
            "offset": offset, "read_len": read_len}


def is_informative(rec, min_mapq, min_flank):
    return (rec["mapq"] >= min_mapq and rec["flank_before"] >= min_flank and rec["flank_after"] >= min_flank)


def length_counts(reads, min_mapq, min_flank):
    return Counter(r["observed_len"] for r in reads if is_informative(r, min_mapq, min_flank))


def allele_stats(counts, ref_len):
    """Descriptive allele-length statistics for one locus x file."""
    n = sum(counts.values())
    top = max(counts.values())
    modal = min(length for length, c in counts.items() if c == top)      # tie -> shorter length
    mean = sum(length * c for length, c in counts.items()) / n
    var = sum(c * (length - mean) ** 2 for length, c in counts.items()) / n
    return {"modal_length": modal, "modal_fraction": counts[modal] / n,
            "ref_allele_fraction": counts.get(ref_len, 0) / n,
            "mean_delta_from_ref": mean - ref_len, "sd_length": var ** 0.5}


def tvd(counts_a, counts_b):
    na, nb = sum(counts_a.values()), sum(counts_b.values())
    keys = set(counts_a) | set(counts_b)
    return 0.5 * sum(abs(counts_a.get(k, 0) / na - counts_b.get(k, 0) / nb) for k in keys)


def hist_string(counts):
    return ";".join(f"{length}:{counts[length]}" for length in sorted(counts))


def summarise_locus_sample(reads, n_excl_flag, ref_len, min_mapq=MIN_MAPQ, min_flank=MIN_FLANK,
                           min_depth=MIN_LOCUS_DEPTH):
    """Pure. `reads` = analyse_alignment() dicts that were NOT flag-excluded."""
    n_mapq_low = sum(1 for r in reads if r["mapq"] < min_mapq)
    n_pass = len(reads) - n_mapq_low
    inf = [r for r in reads if is_informative(r, min_mapq, min_flank)]
    n_inf = len(inf)
    n_overlap = n_excl_flag + len(reads)
    counts = Counter(r["observed_len"] for r in inf)
    status = FS_MISSING if n_inf == 0 else (FS_INSUFFICIENT if n_inf < min_depth else FS_ASSESSED)
    fwd = sum(1 for r in inf if not r["reverse"])
    row = {
        "n_overlap_alignments": n_overlap, "n_excluded_flag": n_excl_flag, "n_excluded_mapq": n_mapq_low,
        "n_pass_filters": n_pass, "n_insufficient_flank": n_pass - n_inf, "n_informative": n_inf,
        "n_informative_fragments": len({r["qname"] for r in inf}),
        "callable_fraction": (n_inf / n_overlap) if n_overlap else None,
        "informative_fraction_of_pass": (n_inf / n_pass) if n_pass else None,
        "median_mapq": median_or_none(r["mapq"] for r in reads),
        "n_informative_fwd": fwd, "n_informative_rev": n_inf - fwd,
        "fwd_fraction": (fwd / n_inf) if n_inf else None,
        "median_rel_offset": median_or_none(r["offset"] / r["read_len"] for r in inf if r["offset"] is not None and r["read_len"]),
        "length_hist": hist_string(counts) if n_inf else None,
        "status": status,
        "modal_length": None, "modal_fraction": None, "ref_allele_fraction": None,
        "mean_delta_from_ref": None, "sd_length": None,
    }
    if status == FS_ASSESSED:
        row.update(allele_stats(counts, ref_len))
    return row


# --------------------------------------------------------------------------- registry / duplicates / reference checks
def load_registry(path):
    with open(path, encoding="utf-8") as fh:
        reg = json.load(fh)
    units = reg.get("biological_units", {})
    ids, seen_triplets = set(), set()
    if not reg.get("files"):
        raise SystemExit("registry lists no files")
    for f in reg["files"]:
        if f["file_id"] in ids:
            raise SystemExit(f"duplicate file_id in registry: {f['file_id']}")
        ids.add(f["file_id"])
        if f["biological_unit_id"] not in units:
            raise SystemExit(f"file {f['file_id']}: unknown biological unit {f['biological_unit_id']}")
        if f["role"] not in ("tumor", "normal"):
            raise SystemExit(f"file {f['file_id']}: role must be tumor or normal")
        trip = (f["biological_unit_id"], f["order_id"], f["role"])
        if trip in seen_triplets:
            raise SystemExit(f"registry lists (unit, order, role) twice: {trip}")
        seen_triplets.add(trip)
    for uid, u in units.items():
        if u.get("label") not in ("MSS", "MSI-H", "UNLABELLED"):
            raise SystemExit(f"unit {uid}: label must be MSS, MSI-H or UNLABELLED")
        if not u.get("label_source"):
            raise SystemExit(f"unit {uid}: label_source is required")
    if len({f["bam"] for f in reg["files"]}) != len(reg["files"]):
        raise SystemExit("registry lists the same BAM path twice")
    if reg.get("off_scope_policy") != "NOT_ASSESSED_residual_reads_ignored":
        raise SystemExit("registry must declare off_scope_policy = NOT_ASSESSED_residual_reads_ignored")
    return reg


def read_fai(path):
    out = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            out.append((p[0], int(p[1])))
    return out


def bam_header_contigs(samtools, bam):
    out = []
    for line in run([samtools, "view", "-H", bam]).splitlines():
        if line.startswith("@SQ"):
            f = dict(x.split(":", 1) for x in line.split("\t")[1:])
            out.append((f["SN"], int(f["LN"])))
    return out


def bam_pg_references(samtools, bam):
    text = run([samtools, "view", "-H", bam])
    return sorted(set(re.findall(r"(\S+\.fa(?:sta)?)\b", " ".join(l for l in text.splitlines() if l.startswith("@PG")))))


def check_reference_compatibility(samtools, registry, fai_contigs, files):
    """Fails loudly on any build / contig mismatch between FASTA, BAMs and scope."""
    fai = dict(fai_contigs)
    scope = registry["analysis_scope_contigs"]
    for c in scope:
        if c not in fai:
            raise SystemExit(f"scope contig {c} is not in the reference index")
    info = {}
    for f in files:
        contigs = bam_header_contigs(samtools, f["bam"])
        if contigs != fai_contigs:
            diff = sorted(set(contigs) ^ set(fai_contigs))[:5]
            raise SystemExit(f"{f['file_id']}: BAM @SQ list differs from the reference index (contig names/lengths/order); e.g. {diff}")
        fasta = registry["reference"]["fasta"]
        refs = bam_pg_references(samtools, f["bam"])
        if fasta not in refs:
            raise SystemExit(f"{f['file_id']}: BAM @PG does not name the registry reference {fasta} (found {refs})")
        stats = [l.split("\t") for l in run([samtools, "idxstats", f["bam"]]).splitlines() if l]
        mapped = {r[0]: int(r[2]) for r in stats if len(r) >= 3}
        primary = sorted(c for c in fai if re.fullmatch(r"chr([0-9]+|X|Y|M)", c))
        empty = [c for c in scope if mapped.get(c, 0) == 0]
        if empty:
            raise SystemExit(f"{f['file_id']}: declared scope contig(s) {empty} hold no reads")
        top = max(primary, key=lambda c: (mapped.get(c, 0), c))
        if top not in scope:
            raise SystemExit(f"{f['file_id']}: the contig holding most reads ({top}) is not in the declared scope {sorted(scope)}")
        in_scope = sum(mapped[c] for c in scope)
        off = sum(mapped.get(c, 0) for c in primary if c not in scope)
        info[f["file_id"]] = {"n_contigs_checked": len(contigs), "reference_in_PG": fasta,
                              "mapped_reads_total": sum(mapped.values()),
                              "mapped_reads_by_scope_contig": {c: mapped[c] for c in scope},
                              "mapped_reads_off_scope_primary_contigs": off,
                              "off_scope_read_fraction_of_primary": round(off / (off + in_scope), 6),
                              "off_scope_contigs_with_reads": sum(1 for c in primary if c not in scope and mapped.get(c, 0) > 0),
                              "off_scope_policy": registry["off_scope_policy"]}
    return info


def check_locus_coordinates(samtools, fasta, fai_contigs, rows, tmpdir):
    """Every union locus: end-start == ref_repeat_length, inside its contig, and the reference
    sequence is exactly the recorded motif repeated. Fails loudly (first offenders listed)."""
    lens = dict(fai_contigs)
    bad, regions = [], []
    for r in rows:
        s, e, L = int(r["start"]), int(r["end"]), int(r["ref_repeat_length"])
        if r["chrom"] not in lens:
            bad.append((r["marker_id"], "contig not in reference index"))
        elif not (0 <= s < e <= lens[r["chrom"]]) or e - s != L:
            bad.append((r["marker_id"], "coordinates inconsistent with ref_repeat_length or contig length"))
        else:
            regions.append((r["marker_id"], f"{r['chrom']}:{s + 1}-{e}", r["motif"]))
    reg_path = os.path.join(tmpdir, "locus_regions.txt")
    with open(reg_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(x[1] for x in regions) + "\n")
    out = run([samtools, "faidx", fasta, "-r", reg_path])
    seqs, cur = [], None
    for line in out.splitlines():
        if line.startswith(">"):
            seqs.append([line[1:], ""])
        else:
            seqs[-1][1] += line
    if len(seqs) != len(regions):
        raise SystemExit("faidx returned an unexpected number of sequences")
    for (mid, region, motif), (name, seq) in zip(regions, seqs):
        seq = seq.upper()
        m = motif.upper()
        if len(seq) % len(m) or seq != m * (len(seq) // len(m)):
            bad.append((mid, f"reference sequence at {region} is not a pure {m} repeat"))
    if bad:
        raise SystemExit(f"{len(bad)} locus coordinate/reference-allele mismatches, e.g. {bad[:5]}")
    return {"n_loci_checked": len(rows), "n_mismatch": 0}


def content_signature(bam, samtools="samtools"):
    """Hash of (QNAME, FLAG, POS, CIGAR) of every alignment in file order."""
    cmd = f"'{samtools}' view '{bam}' | cut -f1,2,4,6 | sha256sum"
    out = subprocess.run(["bash", "-c", "set -o pipefail; " + cmd], capture_output=True, text=True, encoding="utf-8", errors="replace")
    if out.returncode != 0:
        raise SystemExit(f"could not compute content signature of {bam}: {out.stderr.strip()[:200]}")
    return out.stdout.split()[0]


def qname_sample(samtools, bam, contig, n=2000):
    names = []
    proc = subprocess.Popen([samtools, "view", bam, contig], stdout=subprocess.PIPE, text=True, encoding="utf-8", errors="replace")
    for line in proc.stdout:
        names.append(line.split("\t", 1)[0])
        if len(names) >= n:
            break
    proc.stdout.close()
    proc.kill()
    proc.wait()
    return set(names)


def check_duplicates(files, samtools, scope_contig, precomputed=None):
    """Duplicate biological-unit checks. Returns per-file info incl. exact_duplicate_of.
    Fails if identical content / shared read names occur across DIFFERENT biological units."""
    info = {}
    for f in files:
        pre = (precomputed or {}).get(f["file_id"], {})
        info[f["file_id"]] = {
            "bam_sha256": pre.get("bam_sha256") or sha256_file(f["bam"]),
            "content_signature": pre.get("content_signature") or content_signature(f["bam"], samtools),
            "qnames": pre.get("qnames") or qname_sample(samtools, f["bam"], scope_contig),
        }
    for i, a in enumerate(files):
        for b in files[i + 1:]:
            ia, ib = info[a["file_id"]], info[b["file_id"]]
            same_unit = a["biological_unit_id"] == b["biological_unit_id"]
            if not same_unit:
                if ia["bam_sha256"] == ib["bam_sha256"] or ia["content_signature"] == ib["content_signature"]:
                    raise SystemExit(f"identical BAM content under different biological units: {a['file_id']} / {b['file_id']}")
                if ia["qnames"] & ib["qnames"]:
                    raise SystemExit(f"shared read names under different biological units (same sample counted twice?): "
                                     f"{a['file_id']} / {b['file_id']}")
    first_seen = {}
    for f in files:
        sig = (f["biological_unit_id"], f["role"], info[f["file_id"]]["content_signature"])
        info[f["file_id"]]["exact_duplicate_of"] = first_seen.get(sig)
        first_seen.setdefault(sig, f["file_id"])
    return info


# --------------------------------------------------------------------------- BAM extraction (one indexed pass per file)
def extract_file(samtools, bam, loci, tmpdir, tag):
    """loci: list of dict(start,end,chrom,marker_id) for ONE contig, sorted by start.
    One `samtools view -M -L bed` pass; each record is attributed to every overlapping locus.
    Returns {marker_id: (reads, n_excluded_flag)}."""
    bed = os.path.join(tmpdir, f"{tag}.bed")
    with open(bed, "w", encoding="utf-8", newline="\n") as fh:
        for l in loci:
            fh.write(f"{l['chrom']}\t{l['start']}\t{l['end']}\n")
    starts = [l["start"] for l in loci]
    prefix_max_end, cur = [], -1
    for l in loci:
        cur = max(cur, l["end"])
        prefix_max_end.append(cur)
    got = {l["marker_id"]: ([], [0]) for l in loci}
    proc = subprocess.Popen([samtools, "view", "-M", "-L", bed, bam], stdout=subprocess.PIPE, text=True, encoding="utf-8", errors="replace")
    for line in proc.stdout:
        f = line.split("\t", 6)
        flag, pos, mapq, cigar = int(f[1]), int(f[3]), int(f[4]), f[5]
        ops = parse_cigar(cigar)
        if not ops or flag & FLAG_UNMAPPED:
            continue
        r0 = pos - 1
        r1 = r0 + sum(n for n, op in ops if op in "MDN=X")
        i = bisect_left(starts, r1) - 1
        while i >= 0 and prefix_max_end[i] > r0:
            l = loci[i]
            if l["start"] < r1 and l["end"] > r0:
                rec = analyse_alignment(flag, pos, mapq, cigar, l["start"], l["end"])
                if rec is not None:
                    reads, excl = got[l["marker_id"]]
                    if rec["excluded_flag"]:
                        excl[0] += 1
                    else:
                        rec["qname"] = f[0]
                        reads.append(rec)
            i -= 1
    proc.stdout.close()
    if proc.wait() != 0:
        raise SystemExit(f"samtools view failed for {bam}")
    return {k: (v[0], v[1][0]) for k, v in got.items()}


# --------------------------------------------------------------------------- aggregation
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
            + [f"tumor_{c}" for c in ROLE_COLS] + [f"normal_{c}" for c in ROLE_COLS]
            + ["tn_pairs_assessed", "tn_tvd_min", "tn_tvd_median", "tn_tvd_max", "tumor_rerun_modal_length_concordant"]
            + STRAT_COLS)

EM_EVIDENCE_COLS = ["technical_evidence", "coverage_evidence", "stability_evidence", "MSI_separation_evidence",
                    "reproducibility_evidence", "missingness", "sample_count", "ground_truth_count",
                    "ground_truth_MSS_units", "ground_truth_MSI_H_units", "file_count"]
EM_SUPPORT_COLS = ["tech_median_mapq_median", "tech_fwd_fraction_median", "tech_rel_offset_median",
                   "cov_tumor_informative_median", "cov_normal_informative_median", "cov_callable_fraction_median",
                   "stab_tn_tvd_median", "stab_modal_fraction_median"]
EM_COLS = ["candidate", "chrom", "start", "end", "data_scope"] + EM_EVIDENCE_COLS + EM_SUPPORT_COLS + STRAT_COLS

SCOPE_IN = "IN_SCOPE_LOCAL_BAM"
SCOPE_OUT = "OUT_OF_SCOPE_NO_LOCAL_DATA"


def role_summary(rows, role):
    sel = [r for r in rows if r["role"] == role]
    ass = [r for r in sel if r["status"] == FS_ASSESSED]
    inf = [r["n_informative"] for r in sel]
    modal = sorted({r["modal_length"] for r in ass})
    return {
        "n_files": len(sel), "n_assessed": len(ass),
        "informative_min": min(inf) if inf else None, "informative_median": median_or_none(inf),
        "informative_max": max(inf) if inf else None,
        "callable_fraction_median": median_or_none(r["callable_fraction"] for r in sel if r["callable_fraction"] is not None),
        "fwd_fraction_median": median_or_none(r["fwd_fraction"] for r in ass),
        "median_mapq_median": median_or_none(r["median_mapq"] for r in sel if r["median_mapq"] is not None),
        "rel_offset_median": median_or_none(r["median_rel_offset"] for r in ass if r["median_rel_offset"] is not None),
        "modal_length_values": ";".join(str(x) for x in modal) if modal else None,
        "modal_fraction_median": median_or_none(r["modal_fraction"] for r in ass),
        "ref_allele_fraction_median": median_or_none(r["ref_allele_fraction"] for r in ass),
    }


def parse_hist(s):
    return {int(a): int(b) for a, b in (p.split(":") for p in s.split(";"))}


def locus_aggregate(ls_rows):
    """ls_rows: analysis-file rows for ONE in-scope locus (duplicates already removed)."""
    tumor = role_summary(ls_rows, "tumor")
    normal = role_summary(ls_rows, "normal")
    n_files = len(ls_rows)
    n_missing = sum(1 for r in ls_rows if r["status"] == FS_MISSING)
    by_order = {}
    for r in ls_rows:
        by_order.setdefault((r["biological_unit_id"], r["order_id"]), {})[r["role"]] = r
    tvds = []
    for pair in by_order.values():
        t, n = pair.get("tumor"), pair.get("normal")
        if t and n and t["status"] == FS_ASSESSED and n["status"] == FS_ASSESSED:
            tvds.append(round(tvd(parse_hist(t["length_hist"]), parse_hist(n["length_hist"])), 6))
    t_modals = [r["modal_length"] for r in ls_rows if r["role"] == "tumor" and r["status"] == FS_ASSESSED]
    return {
        "n_analysis_files": n_files, "n_missing_files": n_missing,
        "missing_fraction": (n_missing / n_files) if n_files else None,
        "tumor": tumor, "normal": normal,
        "tn_pairs_assessed": len(tvds), "tn_tvds": tvds,
        "tumor_rerun_modal_length_concordant": (len(set(t_modals)) == 1) if len(t_modals) >= 2 else None,
    }


def evidence_row(agg, ls_rows, units, min_units_repro=2):
    """Dimension statuses for one IN-SCOPE locus. Never merges INSUFFICIENT with NOT_ASSESSED."""
    any_assessed = any(r["status"] == FS_ASSESSED for r in ls_rows)
    cov_ok = agg["tumor"]["n_assessed"] >= 1 and agg["normal"]["n_assessed"] >= 1
    stab_ok = agg["tn_pairs_assessed"] >= 1
    assessed_units = {r["biological_unit_id"] for r in ls_rows if r["status"] == FS_ASSESSED}
    labels = Counter(units[u]["label"] for u in assessed_units)
    ass = [r for r in ls_rows if r["status"] == FS_ASSESSED]
    st = lambda ok: ST_ASSESSED if ok else ST_INSUFFICIENT
    return {
        "technical_evidence": st(any_assessed), "coverage_evidence": st(cov_ok), "stability_evidence": st(stab_ok),
        "MSI_separation_evidence": ST_NOT_ASSESSED,
        "reproducibility_evidence": (ST_INSUFFICIENT if len({r["biological_unit_id"] for r in ls_rows}) < min_units_repro
                                     else st(len(assessed_units) >= min_units_repro)),
        "missingness": agg["missing_fraction"], "sample_count": len(assessed_units),
        "ground_truth_count": sum(labels[k] for k in ("MSS", "MSI-H")),
        "ground_truth_MSS_units": labels.get("MSS", 0), "ground_truth_MSI_H_units": labels.get("MSI-H", 0),
        "file_count": len(ls_rows),
        "tech_median_mapq_median": median_or_none(r["median_mapq"] for r in ass) if any_assessed else None,
        "tech_fwd_fraction_median": median_or_none(r["fwd_fraction"] for r in ass) if any_assessed else None,
        "tech_rel_offset_median": median_or_none(r["median_rel_offset"] for r in ass if r["median_rel_offset"] is not None) if any_assessed else None,
        "cov_tumor_informative_median": agg["tumor"]["informative_median"] if cov_ok else None,
        "cov_normal_informative_median": agg["normal"]["informative_median"] if cov_ok else None,
        "cov_callable_fraction_median": median_or_none(r["callable_fraction"] for r in ass) if cov_ok else None,
        "stab_tn_tvd_median": median_or_none(agg["tn_tvds"]) if stab_ok else None,
        "stab_modal_fraction_median": median_or_none(r["modal_fraction"] for r in ass) if stab_ok else None,
    }


def out_of_scope_row():
    na = ST_NOT_ASSESSED
    return {"technical_evidence": na, "coverage_evidence": na, "stability_evidence": na, "MSI_separation_evidence": na,
            "reproducibility_evidence": na, "missingness": na, "sample_count": na, "ground_truth_count": na,
            "ground_truth_MSS_units": na, "ground_truth_MSI_H_units": na, "file_count": na,
            **{c: None for c in EM_SUPPORT_COLS}}


# --------------------------------------------------------------------------- scenarios
def scenario_definitions(rows, sizes, hybrid_cores, hybrid_totals):
    """Rebuilds the Stage-5 scenario memberships IN MEMORY ONLY (no membership file is written)."""
    a_q = t5.ordered(rows, "_inA", "rank_viewA_no_shared_key")
    a_s = t5.ordered(rows, "_inA", "rank_viewA")
    b_q = t5.ordered(rows, "_inB", "rank_viewB_no_shared_key")
    b_s = t5.ordered(rows, "_inB", "rank_viewB")
    defs = []
    for view, q, s in (("viewA", a_q, a_s), ("viewB", b_q, b_s)):
        for size in sizes:
            defs.append({"id": f"{view}_quality_first_{size}", "view": view, "strategy": "quality_first", "size": size,
                         "core": None, "members": q[:size], "core_members": None})
            defs.append({"id": f"{view}_shared_first_{size}", "view": view, "strategy": "shared_first", "size": size,
                         "core": None, "members": s[:size], "core_members": None})
    core_sorted = t5.core_order(rows)
    for k in hybrid_cores:
        for total in hybrid_totals:
            if total <= k:
                continue
            core, ext, members = t5.hybrid_members(core_sorted, b_q, k, total)
            defs.append({"id": f"viewB_hybrid_core{k}_total{total}", "view": "viewB", "strategy": "hybrid", "size": total,
                         "core": k, "members": members, "core_members": core})
    return defs


SC_COLS = ["scenario_id", "view", "strategy", "size", "core", "n_loci", "n_in_scope", "fraction_in_scope",
           "n_out_of_scope_NOT_ASSESSED", "technical_ASSESSED", "technical_INSUFFICIENT",
           "coverage_ASSESSED", "coverage_INSUFFICIENT", "stability_ASSESSED", "stability_INSUFFICIENT",
           "MSI_separation_NOT_ASSESSED", "reproducibility_INSUFFICIENT", "reproducibility_NOT_ASSESSED",
           "missingness_median_in_scope", "cov_tumor_informative_median_p10", "cov_tumor_informative_median_p50",
           "cov_tumor_informative_median_p90", "in_scope_shared", "in_scope_viewB_or_A_only",
           "core_in_scope", "extension_in_scope", "evidence_scope_caveat"]

CAVEAT = "CHR17_ONLY__ONE_BIOLOGICAL_UNIT_MSS_EXTERNAL_LABEL__NO_MSI_DISCRIMINATION__NOT_A_RANKING_OF_STRATEGIES"


def quantile(values, q):
    """Linear-interpolation quantile of a non-empty list (deterministic)."""
    v = sorted(values)
    if not v:
        return None
    pos = (len(v) - 1) * q
    lo = int(pos)
    hi = min(lo + 1, len(v) - 1)
    return round(v[lo] + (v[hi] - v[lo]) * (pos - lo), 6)


def scenario_row(sd, matrix_by_id):
    ms = sd["members"]
    ev = [matrix_by_id[r["marker_id"]] for r in ms]
    inscope = [e for e in ev if e["data_scope"] == SCOPE_IN]
    cnt = lambda key, val: sum(1 for e in inscope if e[key] == val)
    miss = [e["missingness"] for e in inscope]
    covs = [e["cov_tumor_informative_median"] for e in inscope if e["cov_tumor_informative_median"] is not None]
    core_ids = {r["marker_id"] for r in sd["core_members"]} if sd["core_members"] else None
    inscope_ids = {e["candidate"] for e in inscope}
    return [sd["id"], sd["view"], sd["strategy"], sd["size"], sd["core"], len(ms), len(inscope),
            (len(inscope) / len(ms)) if ms else None, len(ms) - len(inscope),
            cnt("technical_evidence", ST_ASSESSED), cnt("technical_evidence", ST_INSUFFICIENT),
            cnt("coverage_evidence", ST_ASSESSED), cnt("coverage_evidence", ST_INSUFFICIENT),
            cnt("stability_evidence", ST_ASSESSED), cnt("stability_evidence", ST_INSUFFICIENT),
            sum(1 for e in ev if e["MSI_separation_evidence"] == ST_NOT_ASSESSED),
            cnt("reproducibility_evidence", ST_INSUFFICIENT),
            sum(1 for e in ev if e["reproducibility_evidence"] == ST_NOT_ASSESSED),
            median_or_none(miss) if miss else None,
            quantile(covs, 0.10) if covs else None, quantile(covs, 0.50) if covs else None,
            quantile(covs, 0.90) if covs else None,
            sum(1 for r in ms if r["_shared"] and r["marker_id"] in inscope_ids),
            sum(1 for r in ms if not r["_shared"] and r["marker_id"] in inscope_ids),
            sum(1 for r in sd["core_members"] if r["marker_id"] in inscope_ids) if core_ids is not None else None,
            sum(1 for r in ms if r["marker_id"] not in core_ids and r["marker_id"] in inscope_ids) if core_ids is not None else None,
            CAVEAT]


# --------------------------------------------------------------------------- sensitivity (descriptive only)
def sensitivity(per_file_reads, ref_lens):
    """per_file_reads: {file_id: {marker_id: reads}}. Counts only; nothing is optimised or selected."""
    out = {}
    for fid, loci in per_file_reads.items():
        def n_loci_meeting(min_mapq, min_flank, depth):
            return sum(1 for reads in loci.values() if sum(1 for r in reads if is_informative(r, min_mapq, min_flank)) >= depth)
        depth_table = {str(d): n_loci_meeting(MIN_MAPQ, MIN_FLANK, d) for d in SENS_DEPTH}
        grid = {f"mapq{m}_flank{f}": n_loci_meeting(m, f, MIN_LOCUS_DEPTH) for m in SENS_MAPQ for f in SENS_FLANK}
        base_modal = {}
        for mid, reads in loci.items():
            c = length_counts(reads, MIN_MAPQ, MIN_FLANK)
            if sum(c.values()) >= MIN_LOCUS_DEPTH:
                base_modal[mid] = allele_stats(c, ref_lens[mid])["modal_length"]
        modal_same = {}
        for m in SENS_MAPQ:
            for f in SENS_FLANK:
                both = same = 0
                for mid, base in base_modal.items():
                    c = length_counts(loci[mid], m, f)
                    if sum(c.values()) >= MIN_LOCUS_DEPTH:
                        both += 1
                        same += allele_stats(c, ref_lens[mid])["modal_length"] == base
                modal_same[f"mapq{m}_flank{f}"] = {"loci_assessed_in_both": both, "modal_length_unchanged": same}
        out[fid] = {"n_loci_with_at_least_depth_at_default_mapq_flank": depth_table,
                    "n_loci_meeting_provisional_depth_by_mapq_flank": grid,
                    "modal_length_agreement_with_default_by_mapq_flank": modal_same}
    return out


# --------------------------------------------------------------------------- orchestration
def run_analysis(rows, registry, fai_contigs, samtools, out_dir, sizes, hybrid_cores, hybrid_totals,
                 precomputed_dup=None, tmpdir=None):
    """Everything after the upstream-hash verification. `rows` = load_union() rows.
    Returns (manifest_fragment, output_paths)."""
    own_tmp = tmpdir is None
    tmp = tempfile.TemporaryDirectory() if own_tmp else None
    tmpdir = tmp.name if own_tmp else tmpdir
    try:
        files = registry["files"]
        units = registry["biological_units"]
        scope = registry["analysis_scope_contigs"]
        ref_info = check_reference_compatibility(samtools, registry, fai_contigs, files)
        coord_info = check_locus_coordinates(samtools, registry["reference"]["fasta"], fai_contigs, rows, tmpdir)
        dup = check_duplicates(files, samtools, scope[0], precomputed_dup)
        analysis_files = [f for f in files if dup[f["file_id"]]["exact_duplicate_of"] is None]

        universe = {r["marker_id"]: r for r in rows}
        in_scope_loci = sorted((r for r in rows if r["chrom"] in scope), key=lambda r: (r["chrom"], int(r["start"])))
        ref_lens = {r["marker_id"]: int(r["ref_repeat_length"]) for r in in_scope_loci}

        ls_rows, per_file_reads = [], {}
        for f in analysis_files:
            unit = units[f["biological_unit_id"]]
            per_file_reads[f["file_id"]] = {}
            for contig in scope:
                loci = [{"marker_id": r["marker_id"], "chrom": r["chrom"], "start": int(r["start"]), "end": int(r["end"])}
                        for r in in_scope_loci if r["chrom"] == contig]
                got = extract_file(samtools, f["bam"], loci, tmpdir, f"{f['file_id']}.{contig}")
                for l in loci:
                    reads, n_excl = got[l["marker_id"]]
                    per_file_reads[f["file_id"]][l["marker_id"]] = reads
                    s = summarise_locus_sample(reads, n_excl, ref_lens[l["marker_id"]])
                    ls_rows.append({"file_id": f["file_id"], "biological_unit_id": f["biological_unit_id"], "role": f["role"],
                                    "order_id": f["order_id"], "dataset_origin": unit["dataset_origin"], "unit_label": unit["label"],
                                    "marker_id": l["marker_id"], "chrom": l["chrom"], "start": l["start"], "end": l["end"],
                                    "ref_repeat_length": ref_lens[l["marker_id"]], **s})
        by_locus = {}
        for r in ls_rows:
            by_locus.setdefault(r["marker_id"], []).append(r)

        lts_out, matrix, matrix_by_id = [], [], {}
        for r in rows:
            mid = r["marker_id"]
            strat = [r[c] for c in STRAT_SRC]
            if mid in by_locus:
                lr = by_locus[mid]
                agg = locus_aggregate(lr)
                ev = evidence_row(agg, lr, units)
                scope_tag = SCOPE_IN
                lts_row = [mid, r["chrom"], r["start"], r["end"], r["ref_repeat_length"], r["motif"], scope_tag,
                           agg["n_analysis_files"], agg["n_missing_files"], agg["missing_fraction"],
                           *[agg["tumor"][c] for c in ROLE_COLS], *[agg["normal"][c] for c in ROLE_COLS],
                           agg["tn_pairs_assessed"], min(agg["tn_tvds"]) if agg["tn_tvds"] else None,
                           median_or_none(agg["tn_tvds"]), max(agg["tn_tvds"]) if agg["tn_tvds"] else None,
                           agg["tumor_rerun_modal_length_concordant"], *strat]
                lts_out.append(lts_row)
            else:
                ev, scope_tag = out_of_scope_row(), SCOPE_OUT
            row = {"candidate": mid, "chrom": r["chrom"], "start": r["start"], "end": r["end"], "data_scope": scope_tag, **ev}
            matrix.append(row)
            matrix_by_id[mid] = row

        scen = [scenario_row(sd, matrix_by_id) for sd in scenario_definitions(rows, sizes, hybrid_cores, hybrid_totals)]

        os.makedirs(out_dir, exist_ok=True)
        paths = {}

        def emit(name, header, body):
            p = os.path.join(out_dir, f"{OUT_PREFIX}.{name}.tsv")
            write_tsv(p, header, body)
            paths[name] = p

        emit("dataset_inventory",
             ["file_id", "biological_unit_id", "role", "order_id", "bam_path", "bam_sha256", "bai_sha256",
              "mapped_reads_total", "mapped_reads_scope_contig", "mapped_reads_off_scope_primary",
              "off_scope_read_fraction", "content_signature", "exact_duplicate_of",
              "analysed", "unit_label", "label_source", "dataset_origin", "assay"],
             [[f["file_id"], f["biological_unit_id"], f["role"], f["order_id"], f["bam"], dup[f["file_id"]]["bam_sha256"],
               sha256_file(f["bam"] + ".bai"), ref_info[f["file_id"]]["mapped_reads_total"],
               ref_info[f["file_id"]]["mapped_reads_by_scope_contig"][scope[0]],
               ref_info[f["file_id"]]["mapped_reads_off_scope_primary_contigs"],
               ref_info[f["file_id"]]["off_scope_read_fraction_of_primary"], dup[f["file_id"]]["content_signature"],
               dup[f["file_id"]]["exact_duplicate_of"], dup[f["file_id"]]["exact_duplicate_of"] is None,
               units[f["biological_unit_id"]]["label"], units[f["biological_unit_id"]]["label_source"],
               units[f["biological_unit_id"]]["dataset_origin"], units[f["biological_unit_id"]]["assay"]] for f in files])
        emit("locus_sample_technical", LS_COLS, [[r[c] for c in LS_COLS] for r in ls_rows])
        emit("locus_technical_summary", LTS_COLS, lts_out)
        emit("evidence_matrix", EM_COLS, [[r["candidate"], r["chrom"], r["start"], r["end"], r["data_scope"],
                                            *[r[c] for c in EM_EVIDENCE_COLS], *[r[c] for c in EM_SUPPORT_COLS],
                                            *[universe[r["candidate"]][c] for c in STRAT_SRC]] for r in matrix])
        emit("scenario_comparison", SC_COLS, scen)

        frag = {
            "reference_checks": ref_info, "locus_coordinate_check": coord_info,
            "duplicate_checks": {fid: {"bam_sha256": v["bam_sha256"], "content_signature": v["content_signature"],
                                       "exact_duplicate_of": v["exact_duplicate_of"], "qname_sample_size": len(v["qnames"])}
                                 for fid, v in dup.items()},
            "analysis_file_ids": [f["file_id"] for f in analysis_files],
            "n_in_scope_loci": len(in_scope_loci), "n_union_loci": len(rows),
            "status_counts_by_file": {fid: dict(Counter(r["status"] for r in ls_rows if r["file_id"] == fid))
                                      for fid in per_file_reads},
            "parameter_sensitivity": sensitivity(per_file_reads, ref_lens),
            "dimension_status_counts": {dim: dict(Counter(r[dim] for r in matrix)) for dim in
                                        ("technical_evidence", "coverage_evidence", "stability_evidence",
                                         "MSI_separation_evidence", "reproducibility_evidence")},
            "n_scenarios": len(scen),
        }
        return frag, paths
    finally:
        if tmp is not None:
            tmp.cleanup()


def norm_path(path):
    """The Stage 1-5 manifests were written on Windows and record backslash separators. They are
    immutable, so separators are normalised at READ time only (the hashes still prove the identity
    of every file)."""
    return path.replace("\\", "/")


def _load_json(path):
    with open(norm_path(path), encoding="utf-8") as fh:
        return json.load(fh)


def _require_hash(path, expected, what):
    got = sha256_file(norm_path(path))
    if got != expected:
        raise SystemExit(f"{what} does not match its recorded hash")
    return got


def verify_stage5(stage5_manifest_path, expect_sha256):
    """Pins the Stage-5 manifest and re-hashes the WHOLE chain below it (stage 5 -> 4 -> 3 -> 2 -> 250 bp).
    Refuses (SystemExit) on any mismatch or any safety flag that is not False. Returns (m5, m4, m3, hashes)."""
    got = sha256_file(norm_path(stage5_manifest_path))
    if got != expect_sha256:
        raise SystemExit(f"stage-5 manifest hash {got} differs from the pinned value {expect_sha256}")
    m5 = _load_json(stage5_manifest_path)
    for flag in ("final_panel_size_fixed", "panel_size_selected", "membership_files_written", "hcc1395_results_used",
                 "classification_performed", "production_changes", "upstream_outputs_modified"):
        if m5[flag] is not False:
            raise SystemExit(f"stage-5 safety flag {flag} is not False")
    hashes = {"stage5_manifest": got}
    for name, info in m5["outputs"].items():
        hashes[f"stage5_{name}"] = _require_hash(info["path"], info["sha256"], f"stage-5 output {name}")

    h4 = _require_hash(m5["stage4_manifest_path"], m5["verified_upstream_hashes"]["stage4_manifest"],
                       "stage-4 manifest (changed since stage 5)")
    m4 = _load_json(m5["stage4_manifest_path"])
    for flag in ("final_panel_size_fixed", "panel_size_selected", "hcc1395_results_used", "classification_performed",
                 "production_changes", "ranking_algorithm_rewritten", "stage3_outputs_modified",
                 "chromosome_is_a_ranking_key", "chromosome_quotas", "chrY_removed"):
        if m4[flag] is not False:
            raise SystemExit(f"stage-4 safety flag {flag} is not False")
    hashes["stage4_manifest"] = h4
    for name, info in m4["outputs"].items():
        hashes[f"stage4_{name}"] = _require_hash(info["path"], info["sha256"], f"stage-4 output {name}")

    h3 = _require_hash(m4["stage3_manifest_path"], m4["verified_upstream_hashes"]["stage3_manifest"],
                       "stage-3 manifest (changed since stage 4)")
    m3 = _load_json(m4["stage3_manifest_path"])
    for flag in ("hcc1395_results_used", "classification_performed", "production_changes", "final_panel_size_fixed",
                 "panel_size_selected", "stage2_outputs_modified"):
        if m3[flag] is not False:
            raise SystemExit(f"stage-3 safety flag {flag} is not False")
    hashes["stage3_manifest"] = h3
    for name, info in m3["outputs"].items():
        hashes[f"stage3_{name}"] = _require_hash(info["path"], info["sha256"], f"stage-3 output {name}")

    hashes["stage2_manifest"] = _require_hash(m3["stage2_manifest_path"], m3["inputs_verified_hashes"]["stage2_manifest"],
                                              "stage-2 manifest (changed since stage 3)")
    m2 = _load_json(m3["stage2_manifest_path"])
    for name, info in m2["outputs"].items():
        hashes[f"stage2_{name}"] = _require_hash(info["path"], info["sha256"], f"stage-2 output {name}")

    hashes["spacing250_manifest"] = _require_hash(m2["prior_manifest"]["path"], m2["prior_manifest"]["sha256"],
                                                  "250 bp manifest")
    prior = _load_json(m2["prior_manifest"]["path"])
    for view, outs in prior["outputs"].items():
        for key, info in outs.items():
            hashes[f"spacing250_{view}_{key}"] = _require_hash(info["path"], info["sha256"], f"250 bp output {view}/{key}")
    return m5, m4, m3, hashes


def write_manifest(*, out_dir, samtools, registry_path, registry, fai_path, stage1_path, union_path, upstream,
                   stage5_manifest_path, expected_pin, sizes, hybrid_cores, hybrid_totals, frag, paths):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    src = {n: sha256_file(os.path.join(script_dir, n)) for n in
           ("empirical_evidence_stage6a.py", "bam_extract.py", "msi_locus_model.py", "tradeoff_pareto_stage5.py",
            "scenario_preview_stage4.py")}
    inputs = {"registry": {"path": registry_path, "sha256": sha256_file(registry_path)},
              "reference_fai": {"path": fai_path, "sha256": sha256_file(fai_path)},
              "stage1_scan_manifest": {"path": stage1_path, "sha256": sha256_file(stage1_path)},
              "stage3_ranked_union": {"path": union_path, "sha256": sha256_file(union_path)},
              "bams": {f["file_id"]: {"path": f["bam"], "bam_sha256": frag["duplicate_checks"][f["file_id"]]["bam_sha256"],
                                      "bai_sha256": sha256_file(f["bam"] + ".bai")} for f in registry["files"]}}
    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "stage": "Stage 6A empirical evidence: technical + MSS-only descriptives (development / RUO)",
        "status": "DEVELOPMENT_ANALYSIS_NO_MSI_DISCRIMINATION_NOT_A_SELECTED_PANEL",
        "stage_6B": {"status": "SPECIFIED_IN_DESIGN_GATED_NOT_EXECUTED",
                     "gate": "independent MSI-H and MSS labelled biological units, approved by the product owner",
                     "executed": False},
        "script": {"version": SCRIPT_VERSION, "path": os.path.abspath(__file__), "sha256": src["empirical_evidence_stage6a.py"],
                   "python": platform.python_version(), "samtools": run([samtools, "--version"]).splitlines()[0]},
        "source_files_sha256": src,
        "stage5_manifest_path": stage5_manifest_path, "expected_stage5_manifest_sha256": expected_pin,
        "verified_upstream_hashes": upstream,
        "inputs": inputs,
        "parameters": {
            "MIN_MAPQ": {"value": MIN_MAPQ, "provenance": "dev/msi_redesign/bam_extract.py (early prototype)", "validated": False, "status": "PROVISIONAL_DEVELOPMENT_PARAMETER"},
            "MIN_FLANK": {"value": MIN_FLANK, "provenance": "dev/msi_redesign/bam_extract.py (early prototype)", "validated": False, "status": "PROVISIONAL_DEVELOPMENT_PARAMETER"},
            "MIN_LOCUS_DEPTH": {"value": MIN_LOCUS_DEPTH, "provenance": "dev/msi_redesign/msi_locus_model.py (early prototype)", "validated": False, "status": "PROVISIONAL_DEVELOPMENT_PARAMETER"},
            "not_optimised_or_calibrated_on_HCC1395": True,
            "read_filter": "unmapped/secondary/duplicate/supplementary excluded by FLAG; MAPQ >= MIN_MAPQ; aligned flank >= MIN_FLANK on both sides",
            "sensitivity_grid": {"mapq": SENS_MAPQ, "flank": SENS_FLANK, "depth": SENS_DEPTH, "purpose": "descriptive counts only; no value is selected"},
            "scope_contigs": registry["analysis_scope_contigs"], "sizes": sizes, "hybrid_cores": hybrid_cores,
            "hybrid_totals": hybrid_totals},
        "definitions": {
            "NOT_ASSESSED": "analysis not performed or not possible (locus outside the declared local data scope, or Stage 6B not run). "
                            "Off-scope contigs hold only sparse residual reads (mates of in-scope reads); loci there are never queried and are not assessed",
            "INSUFFICIENT": "analysis attempted but data below the declared provisional minimum, or fewer than the required number of biological units",
            "technical_evidence": "ASSESSED if >=1 analysis file has >= MIN_LOCUS_DEPTH informative reads at the locus",
            "coverage_evidence": "ASSESSED if >=1 tumor AND >=1 normal analysis file each reach MIN_LOCUS_DEPTH informative reads",
            "stability_evidence": "MSS-only descriptive: ASSESSED if >=1 same-order tumor/normal pair are both ASSESSED (tumor-normal TVD computable); NOT a threshold and NOT MSI evidence",
            "MSI_separation_evidence": "always NOT_ASSESSED in Stage 6A",
            "reproducibility_evidence": "INSUFFICIENT in scope (needs >=2 independent biological units; one exists); re-processings of the same FASTQs are not independent",
            "missingness": "fraction of analysis files with zero informative reads at the locus; never imputed",
            "sample_count": "independent biological units with >=1 ASSESSED file at the locus (files of one unit are one sample)"},
        "evidence_columns": EM_EVIDENCE_COLS + EM_SUPPORT_COLS,
        "stratification_columns_descriptive_only": STRAT_COLS,
        "claims_not_made": CLAIMS_NOT_MADE,
        **frag,
        "final_panel_size_fixed": False, "panel_size_selected": False, "membership_files_written": False,
        "msi_classification_performed": False, "msi_discrimination_performed": False, "stage_6b_executed": False,
        "clinical_claims": False, "composite_score_computed": False, "production_changes": False,
        "upstream_outputs_modified": False, "downloads_performed": False,
        "outputs": {k: {"path": v, "sha256": sha256_file(v)} for k, v in paths.items()},
    }
    manifest_path = os.path.join(out_dir, f"{OUT_PREFIX}.manifest.json")
    with open(manifest_path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
        fh.write("\n")
    return manifest_path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--stage5-manifest", required=True)
    ap.add_argument("--registry", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--expect-stage5-manifest-sha256", default=EXPECTED_STAGE5_MANIFEST_SHA256)
    ap.add_argument("--samtools", default="samtools")
    args = ap.parse_args(argv)

    m5, m4, m3, hashes = verify_stage5(args.stage5_manifest, args.expect_stage5_manifest_sha256)
    registry = load_registry(args.registry)
    fasta = registry["reference"]["fasta"]
    fai_path = fasta + ".fai"
    stage1 = "dev/msi_redesign/output/msi_markers.v2.manifest.json"
    with open(stage1, encoding="utf-8") as fh:
        m_scan = json.load(fh)
    if m_scan["source_reference_path"] != fasta or m_scan["source_reference_fai_sha256"] != sha256_file(fai_path):
        raise SystemExit("registry reference differs from the reference used to scan the markers (path or .fai hash)")
    fai_contigs = read_fai(fai_path)

    union_path = norm_path(m3["outputs"]["ranked_union"]["path"])
    rows = t5.load_union(union_path)
    p = m5["parameters"]
    frag, paths = run_analysis(rows, registry, fai_contigs, args.samtools, args.out_dir, p["sizes"], p["hybrid_cores"],
                               p["hybrid_totals"])

    manifest_path = write_manifest(
        out_dir=args.out_dir, samtools=args.samtools, registry_path=args.registry, registry=registry, fai_path=fai_path,
        stage1_path=stage1, union_path=union_path, upstream=hashes,
        stage5_manifest_path=args.stage5_manifest, expected_pin=args.expect_stage5_manifest_sha256,
        sizes=p["sizes"], hybrid_cores=p["hybrid_cores"], hybrid_totals=p["hybrid_totals"], frag=frag, paths=paths)
    print(f"wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
