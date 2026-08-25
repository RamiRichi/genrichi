#!/usr/bin/env python3
"""
GenRichi Somatic Hotspot Panel — Web Application
Wraps the full somatic pipeline: fastp → BWA → MarkDup → BQSR →
Mutect2 → FilterMutectCalls → VEP → HTML Report
"""

import os, sys, uuid, json, sqlite3, threading, subprocess, gzip, re, base64
from datetime import datetime
from pathlib import Path
from flask import Flask, request, jsonify, render_template, send_file, abort

app = Flask(__name__)

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE        = Path("/home/rami/genrichi")
RESULTS     = BASE / "results"
LOGS_DIR    = BASE / "logs" / "app"
DB_FILE     = BASE / "genrichi_jobs.db"

REF = {
    "genome":      "/home/rami/genrichi/reference_db/ref/hg38.fa",
    "dbsnp":       "/home/rami/genrichi/reference_db/dbsnp/dbsnp_146.hg38.vcf.gz",
    "gnomad":      "/home/rami/genrichi/reference_db/gnomad/af-only-gnomad.hg38.vcf.gz",
    "clinvar":     "/home/rami/genrichi/reference_db/clinvar/clinvar.vcf.gz",
    "vep_cache":   "/home/rami/genrichi/reference_db/vep_cache",
    "default_bed": "/home/rami/genrichi/resources/panel/example_hotspots.bed",
}

STEPS = [
    {"id": "qc",       "name": "QC + Trimming (fastp)"},
    {"id": "align",    "name": "Alignment (BWA)"},
    {"id": "dedup",    "name": "Mark Duplicates (GATK)"},
    {"id": "bqsr",     "name": "BQSR (GATK)"},
    {"id": "mutect2",  "name": "Variant Calling (Mutect2)"},
    {"id": "filter",   "name": "Filter Variants"},
    {"id": "annotate", "name": "Annotation (VEP)"},
    {"id": "report",   "name": "HTML Report"},
]

RESULTS.mkdir(parents=True, exist_ok=True)
LOGS_DIR.mkdir(parents=True, exist_ok=True)

# ── Database ──────────────────────────────────────────────────────────────────
def db_connect():
    conn = sqlite3.connect(str(DB_FILE))
    conn.row_factory = sqlite3.Row
    return conn

def db_init():
    with db_connect() as conn:
        conn.execute("""
            CREATE TABLE IF NOT EXISTS jobs (
                id          TEXT PRIMARY KEY,
                sample_id   TEXT NOT NULL,
                r1          TEXT NOT NULL,
                r2          TEXT NOT NULL,
                bed         TEXT NOT NULL,
                status      TEXT DEFAULT 'queued',
                step        TEXT DEFAULT '',
                progress    INTEGER DEFAULT 0,
                created_at  TEXT,
                finished_at TEXT,
                error_msg   TEXT,
                report_path TEXT,
                vcf_path    TEXT
            )
        """)
        conn.execute("""
            CREATE TABLE IF NOT EXISTS job_logs (
                id      INTEGER PRIMARY KEY AUTOINCREMENT,
                job_id  TEXT NOT NULL,
                ts      TEXT,
                line    TEXT
            )
        """)

def job_get(job_id):
    with db_connect() as conn:
        row = conn.execute("SELECT * FROM jobs WHERE id=?", (job_id,)).fetchone()
    return dict(row) if row else None

def job_update(job_id, **kw):
    sets = ", ".join(f"{k}=?" for k in kw)
    vals = list(kw.values()) + [job_id]
    with db_connect() as conn:
        conn.execute(f"UPDATE jobs SET {sets} WHERE id=?", vals)

def log_append(job_id, line):
    ts = datetime.now().strftime("%H:%M:%S")
    with db_connect() as conn:
        conn.execute("INSERT INTO job_logs(job_id,ts,line) VALUES(?,?,?)",
                     (job_id, ts, line))

def logs_get(job_id, offset=0):
    with db_connect() as conn:
        rows = conn.execute(
            "SELECT ts, line FROM job_logs WHERE job_id=? ORDER BY id LIMIT 500 OFFSET ?",
            (job_id, offset)
        ).fetchall()
    return [{"ts": r["ts"], "line": r["line"]} for r in rows]

# ── Shell helper ──────────────────────────────────────────────────────────────
def run_cmd(job_id, cmd, step_name, env=None):
    log_append(job_id, f"▶ {step_name}")
    log_append(job_id, f"$ {' '.join(str(c) for c in cmd)}")
    e = os.environ.copy()
    if env:
        e.update(env)
    proc = subprocess.Popen(
        [str(c) for c in cmd],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, env=e
    )
    for line in proc.stdout:
        log_append(job_id, line.rstrip())
    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError(f"{step_name} failed (exit {proc.returncode})")
    log_append(job_id, f"✓ {step_name} done")

# ── VCF → TSV ─────────────────────────────────────────────────────────────────
def vcf_to_tsv(vcf_gz, tsv_out):
    rows = []
    opener = gzip.open if str(vcf_gz).endswith(".gz") else open
    with opener(str(vcf_gz), "rt") as fh:
        for raw in fh:
            if raw.startswith("#"):
                continue
            f = raw.rstrip().split("\t")
            if len(f) < 8:
                continue
            chrom, pos, _, ref, alt, _, filt, info = f[:8]
            if filt not in ("PASS", "."):
                continue
            fmt_keys = f[8].split(":") if len(f) > 8 else []
            fmt_vals = f[9].split(":") if len(f) > 9 else []
            fmt = dict(zip(fmt_keys, fmt_vals))
            # AF
            af = "."
            if "AF" in fmt:
                af = fmt["AF"]
            else:
                m = re.search(r"AF=([\d.eE+\-]+)", info)
                if m:
                    af = m.group(1)
            # DP
            dp = fmt.get("DP", ".")
            # Gene/consequence from CSQ (VEP)
            gene, cons, hgvsc, hgvsp, impact = ".", ".", ".", ".", "."
            m = re.search(r"CSQ=([^;]+)", info)
            if m:
                csq = m.group(1).split(",")[0].split("|")
                if len(csq) > 6:
                    cons   = csq[1]  if len(csq) > 1  else "."
                    impact = csq[2]  if len(csq) > 2  else "."
                    gene   = csq[3]  if len(csq) > 3  else "."
                    hgvsc  = csq[10] if len(csq) > 10 else "."
                    hgvsp  = csq[11] if len(csq) > 11 else "."
            rows.append([chrom, pos, ref, alt, filt, af, dp, gene, cons, impact, hgvsc, hgvsp])
    header = "CHROM\tPOS\tREF\tALT\tFILTER\tAF\tDP\tGENE\tCONSEQUENCE\tIMPACT\tHGVSc\tHGVSp\n"
    with open(str(tsv_out), "w") as out:
        out.write(header)
        for r in rows:
            out.write("\t".join(r) + "\n")
    return len(rows)

# ── HTML Report ───────────────────────────────────────────────────────────────
def make_html_report(job_id, sample_id, dirs):
    import json as _json

    # ── read QC stats ──────────────────────────────────────────────────────
    qc_stats = {}
    fastp_json = dirs["qc"] / f"{sample_id}_fastp.json"
    if fastp_json.exists():
        try:
            d = _json.loads(fastp_json.read_text())
            s = d.get("summary", {})
            qc_stats["total_reads"]     = s.get("before_filtering", {}).get("total_reads", "N/A")
            qc_stats["q30_rate"]        = round(s.get("before_filtering", {}).get("q30_rate", 0) * 100, 1)
            qc_stats["total_reads_out"] = s.get("after_filtering",  {}).get("total_reads", "N/A")
        except Exception:
            pass

    # ── read flagstat ──────────────────────────────────────────────────────
    flag_stats = {}
    flagstat_f = dirs["align"] / f"{sample_id}.flagstat"
    if flagstat_f.exists():
        for ln in flagstat_f.read_text().splitlines():
            if "mapped (" in ln:
                flag_stats["mapped"] = ln.split(" + ")[0].strip()
            if "in total" in ln:
                flag_stats["total"] = ln.split(" + ")[0].strip()
            if "duplicates" in ln:
                flag_stats["duplicates"] = ln.split(" + ")[0].strip()

    # ── mean coverage from mosdepth ────────────────────────────────────────
    mean_cov = "N/A"
    mossum = dirs["align"] / f"{sample_id}.mosdepth.summary.txt"
    if not mossum.exists():
        mossum = dirs["align"] / f"{sample_id}.mosdepth.global.dist.txt"
    # try regions summary
    mos2 = dirs["align"] / f"{sample_id}.regions.bed.gz"
    mossum_txt = dirs["align"] / f"{sample_id}.mosdepth.summary.txt"
    if mossum_txt.exists():
        for ln in mossum_txt.read_text().splitlines():
            if ln.startswith("total") or ln.startswith("Region"):
                continue
            parts = ln.split("\t")
            if len(parts) >= 4:
                try:
                    mean_cov = f"{float(parts[3]):.1f}x"
                    break
                except Exception:
                    pass

    # ── read variants TSV ──────────────────────────────────────────────────
    variants = []
    tsv_f = dirs["annotation"] / "variants.tsv"
    if tsv_f.exists():
        lines = tsv_f.read_text().splitlines()
        for ln in lines[1:]:
            parts = ln.split("\t")
            if len(parts) >= 12:
                variants.append(parts)

    # ── charts ─────────────────────────────────────────────────────────────
    charts = {}
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import io

        # VAF histogram
        afs = []
        for v in variants:
            try:
                afs.append(float(v[5]))
            except Exception:
                pass
        if afs:
            fig, ax = plt.subplots(figsize=(5, 3))
            ax.hist(afs, bins=20, color="#2563eb", edgecolor="white")
            ax.set_xlabel("Allele Frequency")
            ax.set_ylabel("Variant Count")
            ax.set_title("VAF Distribution")
            ax.spines[["top", "right"]].set_visible(False)
            buf = io.BytesIO()
            fig.tight_layout()
            fig.savefig(buf, format="png", dpi=100)
            plt.close(fig)
            charts["vaf"] = base64.b64encode(buf.getvalue()).decode()

        # Impact pie
        impacts = {}
        for v in variants:
            imp = v[9] if len(v) > 9 else "UNKNOWN"
            impacts[imp] = impacts.get(imp, 0) + 1
        if impacts:
            fig, ax = plt.subplots(figsize=(4, 3))
            cols = {"HIGH": "#dc2626", "MODERATE": "#f59e0b",
                    "LOW": "#22c55e", "MODIFIER": "#94a3b8", ".": "#cbd5e1"}
            clrs = [cols.get(k, "#94a3b8") for k in impacts]
            ax.pie(impacts.values(), labels=impacts.keys(),
                   colors=clrs, autopct="%1.0f%%", startangle=140)
            ax.set_title("Variant Impact")
            buf = io.BytesIO()
            fig.tight_layout()
            fig.savefig(buf, format="png", dpi=100)
            plt.close(fig)
            charts["impact"] = base64.b64encode(buf.getvalue()).decode()
    except Exception as e:
        pass  # charts are optional

    # ── build HTML ─────────────────────────────────────────────────────────
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    var_rows = ""
    for v in variants:
        chrom, pos, ref, alt, filt, af, dp, gene, cons, impact, hgvsc, hgvsp = (v + ["."]*12)[:12]
        badge_col = {"HIGH": "#dc2626", "MODERATE": "#f59e0b",
                     "LOW": "#22c55e", "MODIFIER": "#6b7280"}.get(impact, "#6b7280")
        try:
            af_pct = f"{float(af)*100:.1f}%"
        except Exception:
            af_pct = af
        var_rows += f"""
        <tr>
          <td>{gene}</td><td>{chrom}:{pos}</td><td>{ref}</td><td>{alt}</td>
          <td>{af_pct}</td><td>{dp}</td>
          <td><span style="background:{badge_col};color:#fff;padding:2px 8px;border-radius:9px;font-size:12px">{impact}</span></td>
          <td style="font-size:12px;color:#64748b">{cons}</td>
          <td style="font-size:11px">{hgvsc}</td>
          <td style="font-size:11px">{hgvsp}</td>
        </tr>"""

    vaf_img = f'<img src="data:image/png;base64,{charts["vaf"]}" style="max-width:100%">' if "vaf" in charts else ""
    imp_img = f'<img src="data:image/png;base64,{charts["impact"]}" style="max-width:100%">' if "impact" in charts else ""

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>GenRichi Report — {sample_id}</title>
<style>
  body {{font-family:'Segoe UI',Arial,sans-serif;margin:0;background:#f8fafc;color:#1e293b}}
  .header {{background:linear-gradient(135deg,#1e3a5f,#2563eb);color:#fff;padding:28px 40px;display:flex;align-items:center;gap:20px}}
  .header h1 {{margin:0;font-size:22px}}
  .header p  {{margin:4px 0 0;opacity:.75;font-size:13px}}
  .body  {{padding:32px 40px}}
  .cards {{display:grid;grid-template-columns:repeat(4,1fr);gap:16px;margin-bottom:28px}}
  .card  {{background:#fff;border-radius:10px;padding:20px;box-shadow:0 1px 4px rgba(0,0,0,.07)}}
  .card .val {{font-size:26px;font-weight:700;color:#2563eb}}
  .card .lbl {{font-size:12px;color:#64748b;margin-top:4px}}
  .section {{background:#fff;border-radius:10px;box-shadow:0 1px 4px rgba(0,0,0,.07);margin-bottom:24px}}
  .section-header {{padding:16px 24px;border-bottom:1px solid #e2e8f0;font-weight:600;color:#1e293b}}
  .charts {{display:grid;grid-template-columns:1fr 1fr;gap:16px;padding:20px}}
  table {{width:100%;border-collapse:collapse;font-size:13px}}
  th {{background:#f1f5f9;padding:10px 12px;text-align:left;font-weight:600;font-size:12px;color:#475569;position:sticky;top:0}}
  td {{padding:9px 12px;border-bottom:1px solid #f1f5f9}}
  tr:hover td {{background:#f8fafc}}
  .footer {{text-align:center;color:#94a3b8;font-size:12px;padding:20px}}
</style>
</head>
<body>
<div class="header">
  <div>
    <h1>GenRichi Diagnostics — Somatic Panel Report</h1>
    <p>Sample: <strong>{sample_id}</strong> &nbsp;|&nbsp; Generated: {now} &nbsp;|&nbsp; Pipeline: BWA / GATK Mutect2 / VEP</p>
  </div>
</div>
<div class="body">
  <div class="cards">
    <div class="card"><div class="val">{qc_stats.get('total_reads','N/A')}</div><div class="lbl">Total Reads (pre-QC)</div></div>
    <div class="card"><div class="val">{qc_stats.get('q30_rate','N/A')}%</div><div class="lbl">Q30 Rate</div></div>
    <div class="card"><div class="val">{mean_cov}</div><div class="lbl">Mean Panel Coverage</div></div>
    <div class="card"><div class="val">{len(variants)}</div><div class="lbl">PASS Variants</div></div>
  </div>

  <div class="section">
    <div class="section-header">QC &amp; Alignment Summary</div>
    <div style="padding:16px 24px;display:grid;grid-template-columns:repeat(3,1fr);gap:12px;font-size:13px">
      <div><span style="color:#64748b">Reads after trimming:</span> <strong>{qc_stats.get('total_reads_out','N/A')}</strong></div>
      <div><span style="color:#64748b">Mapped reads:</span> <strong>{flag_stats.get('mapped','N/A')}</strong></div>
      <div><span style="color:#64748b">Duplicate reads:</span> <strong>{flag_stats.get('duplicates','N/A')}</strong></div>
    </div>
  </div>

  <div class="section">
    <div class="section-header">Variant Charts</div>
    <div class="charts">
      <div>{vaf_img if vaf_img else '<p style="color:#94a3b8;padding:20px">No VAF data</p>'}</div>
      <div>{imp_img if imp_img else '<p style="color:#94a3b8;padding:20px">No impact data</p>'}</div>
    </div>
  </div>

  <div class="section">
    <div class="section-header">Detected Variants ({len(variants)} PASS)</div>
    <div style="overflow-x:auto">
      <table>
        <thead><tr>
          <th>Gene</th><th>Position</th><th>Ref</th><th>Alt</th>
          <th>AF</th><th>Depth</th><th>Impact</th><th>Consequence</th>
          <th>HGVSc</th><th>HGVSp</th>
        </tr></thead>
        <tbody>{var_rows if var_rows else '<tr><td colspan="10" style="text-align:center;color:#94a3b8;padding:24px">No PASS variants detected</td></tr>'}</tbody>
      </table>
    </div>
  </div>
</div>
<div class="footer">GenRichi Diagnostics &copy; {datetime.now().year} &nbsp;|&nbsp; For research use only</div>
</body></html>"""

    report_path = dirs["report"] / f"{sample_id}_report.html"
    dirs["report"].mkdir(parents=True, exist_ok=True)
    report_path.write_text(html)
    return report_path

# ── Pipeline worker ───────────────────────────────────────────────────────────
def run_pipeline(job_id, sample_id, r1, r2, bed):
    def step(sid, pct):
        job_update(job_id, step=sid, progress=pct, status="running")

    dirs = {
        "qc":         RESULTS / sample_id / "qc",
        "align":      RESULTS / sample_id / "align",
        "calling":    RESULTS / sample_id / "calling",
        "annotation": RESULTS / sample_id / "annotation",
        "report":     RESULTS / sample_id / "report",
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)

    try:
        # ── 1. QC ──────────────────────────────────────────────────────────
        step("qc", 5)
        r1_trim = dirs["qc"] / f"{sample_id}_R1_trimmed.fastq.gz"
        r2_trim = dirs["qc"] / f"{sample_id}_R2_trimmed.fastq.gz"
        run_cmd(job_id, [
            "fastp",
            "-i", r1, "-I", r2,
            "-o", r1_trim, "-O", r2_trim,
            "--json", dirs["qc"] / f"{sample_id}_fastp.json",
            "--html", dirs["qc"] / f"{sample_id}_fastp.html",
            "--thread", "4",
            "--detect_adapter_for_pe",
            "--qualified_quality_phred", "20",
            "--length_required", "50",
        ], "QC / fastp")

        # ── 2. Alignment ───────────────────────────────────────────────────
        step("align", 15)
        raw_bam = dirs["align"] / f"{sample_id}.raw.bam"
        rg = f"@RG\\tID:{sample_id}\\tSM:{sample_id}\\tPL:ILLUMINA\\tLB:{sample_id}_lib1\\tPU:{sample_id}"
        align_log = dirs["align"] / "bwa.log"
        log_append(job_id, "▶ Alignment (BWA mem)")
        bwa_proc = subprocess.Popen(
            ["bwa", "mem", "-t", "8", "-R", rg,
             REF["genome"], str(r1_trim), str(r2_trim)],
            stdout=subprocess.PIPE, stderr=open(str(align_log), "w")
        )
        sort_proc = subprocess.Popen(
            ["samtools", "sort", "-@", "8", "-o", str(raw_bam), "-"],
            stdin=bwa_proc.stdout
        )
        bwa_proc.stdout.close()
        sort_proc.wait()
        bwa_proc.wait()
        if bwa_proc.returncode not in (0, None) or sort_proc.returncode != 0:
            raise RuntimeError(f"BWA/sort failed: bwa={bwa_proc.returncode}, sort={sort_proc.returncode}")
        log_append(job_id, "✓ Alignment done")

        # ── 3. Mark Duplicates ─────────────────────────────────────────────
        step("dedup", 30)
        markdup_bam = dirs["align"] / f"{sample_id}.markdup.bam"
        markdup_bai = dirs["align"] / f"{sample_id}.markdup.bam.bai"
        markdup_met = dirs["align"] / f"{sample_id}.markdup_metrics.txt"
        run_cmd(job_id, [
            "gatk", "MarkDuplicates",
            "-I", raw_bam,
            "-O", markdup_bam,
            "-M", markdup_met,
            "--CREATE_INDEX", "true",
            "--VALIDATION_STRINGENCY", "SILENT",
        ], "Mark Duplicates")
        # ensure index
        bai_src = dirs["align"] / f"{sample_id}.markdup.bai"
        if bai_src.exists() and not markdup_bai.exists():
            bai_src.rename(markdup_bai)
        if not markdup_bai.exists():
            subprocess.run(["samtools", "index", str(markdup_bam), str(markdup_bai)], check=True)

        # ── 4. BQSR ────────────────────────────────────────────────────────
        step("bqsr", 43)
        recal_table = dirs["align"] / f"{sample_id}.recal.table"
        run_cmd(job_id, [
            "gatk", "BaseRecalibrator",
            "-I", markdup_bam,
            "-R", REF["genome"],
            "--known-sites", REF["dbsnp"],
            "-L", bed,
            "-O", recal_table,
        ], "BaseRecalibrator")

        final_bam = dirs["align"] / f"{sample_id}.final.bam"
        final_bai = dirs["align"] / f"{sample_id}.final.bam.bai"
        run_cmd(job_id, [
            "gatk", "ApplyBQSR",
            "-I", markdup_bam,
            "-R", REF["genome"],
            "--bqsr-recal-file", recal_table,
            "-O", final_bam,
        ], "ApplyBQSR")
        if not final_bai.exists():
            subprocess.run(["samtools", "index", str(final_bam), str(final_bai)], check=True)

        # ── Coverage (mosdepth) ────────────────────────────────────────────
        step("bqsr", 50)
        run_cmd(job_id, [
            "mosdepth",
            "--threads", "4",
            "--by", bed,
            str(dirs["align"] / sample_id),
            str(final_bam),
        ], "mosdepth coverage")

        # ── 5. Mutect2 ─────────────────────────────────────────────────────
        step("mutect2", 55)
        raw_vcf = dirs["calling"] / f"{sample_id}.mutect2.vcf.gz"
        mutect_stats = dirs["calling"] / f"{sample_id}.mutect2.vcf.gz.stats"
        cmd_mutect = [
            "gatk", "Mutect2",
            "-R", REF["genome"],
            "-I", final_bam,
            "-tumor", sample_id,
            "--germline-resource", REF["gnomad"],
            "-L", bed,
            "-O", raw_vcf,
        ]
        run_cmd(job_id, cmd_mutect, "Mutect2")

        # ── Read orientation model ─────────────────────────────────────────
        step("mutect2", 62)
        f1r2_tar = dirs["calling"] / f"{sample_id}.f1r2.tar.gz"
        rom_vcf   = dirs["calling"] / f"{sample_id}.mutect2.vcf.gz"
        # rerun Mutect2 with --f1r2-tar-gz if not present
        if not f1r2_tar.exists():
            cmd_mutect2 = cmd_mutect + ["--f1r2-tar-gz", str(f1r2_tar)]
            run_cmd(job_id, cmd_mutect2, "Mutect2 (with f1r2)")

        orient_model = dirs["calling"] / f"{sample_id}.orient_model.tar.gz"
        if f1r2_tar.exists():
            run_cmd(job_id, [
                "gatk", "LearnReadOrientationModel",
                "-I", f1r2_tar,
                "-O", orient_model,
            ], "LearnReadOrientationModel")

        # ── Contamination ──────────────────────────────────────────────────
        step("filter", 67)
        pileup_table = dirs["calling"] / f"{sample_id}.pileup.table"
        contam_table = dirs["calling"] / f"{sample_id}.contamination.table"
        run_cmd(job_id, [
            "gatk", "GetPileupSummaries",
            "-I", final_bam,
            "-V", REF["gnomad"],
            "-L", bed,
            "-O", pileup_table,
        ], "GetPileupSummaries")
        run_cmd(job_id, [
            "gatk", "CalculateContamination",
            "-I", pileup_table,
            "-O", contam_table,
        ], "CalculateContamination")

        # ── FilterMutectCalls ──────────────────────────────────────────────
        step("filter", 72)
        filtered_vcf = dirs["calling"] / f"{sample_id}.filtered.vcf.gz"
        filter_cmd = [
            "gatk", "FilterMutectCalls",
            "-R", REF["genome"],
            "-V", raw_vcf,
            "--stats", mutect_stats,
            "--contamination-table", contam_table,
            "-O", filtered_vcf,
        ]
        if orient_model.exists():
            filter_cmd += ["--ob-priors", orient_model]
        run_cmd(job_id, filter_cmd, "FilterMutectCalls")

        # PASS only
        pass_vcf = dirs["calling"] / f"{sample_id}.pass.vcf.gz"
        subprocess.run(
            f"bcftools view -f PASS {filtered_vcf} | bgzip -c > {pass_vcf}",
            shell=True, check=True
        )
        subprocess.run(["tabix", "-p", "vcf", str(pass_vcf)], check=True)
        log_append(job_id, f"✓ PASS VCF: {pass_vcf}")

        # ── 6. VEP Annotation ──────────────────────────────────────────────
        step("annotate", 78)
        ann_vcf = dirs["annotation"] / f"{sample_id}.vep.vcf.gz"
        vep_cmd = [
            "vep",
            "--input_file",  str(pass_vcf),
            "--output_file", str(ann_vcf),
            "--format", "vcf",
            "--vcf",
            "--compress_output", "bgzip",
            "--cache",
            "--dir_cache", REF["vep_cache"],
            "--assembly", "GRCh38",
            "--species", "homo_sapiens",
            "--offline",
            "--everything",
            "--fork", "4",
            "--force_overwrite",
        ]
        # ClinVar custom track
        clinvar = Path(REF["clinvar"])
        if clinvar.exists():
            vep_cmd += [
                "--custom", f"{clinvar},ClinVar,vcf,exact,0,CLNSIG,CLNDN"
            ]
        run_cmd(job_id, vep_cmd, "VEP annotation")
        subprocess.run(["tabix", "-p", "vcf", str(ann_vcf)], check=False)

        # VCF → TSV
        step("annotate", 85)
        tsv_path = dirs["annotation"] / "variants.tsv"
        n_vars = vcf_to_tsv(ann_vcf, tsv_path)
        log_append(job_id, f"✓ {n_vars} PASS variants exported to TSV")

        # ── 7. Report ──────────────────────────────────────────────────────
        step("report", 90)
        report_path = make_html_report(job_id, sample_id, dirs)
        log_append(job_id, f"✓ Report: {report_path}")

        job_update(
            job_id,
            status="done",
            step="report",
            progress=100,
            finished_at=datetime.now().isoformat(),
            report_path=str(report_path),
            vcf_path=str(ann_vcf),
        )

    except Exception as exc:
        import traceback
        log_append(job_id, f"✗ ERROR: {exc}")
        log_append(job_id, traceback.format_exc())
        job_update(
            job_id,
            status="error",
            finished_at=datetime.now().isoformat(),
            error_msg=str(exc),
        )

# ── Routes ────────────────────────────────────────────────────────────────────
@app.route("/")
def index():
    return render_template("dashboard.html", steps=STEPS)

@app.route("/api/submit", methods=["POST"])
def submit():
    data = request.json or {}
    sample_id = (data.get("sample_id") or "").strip()
    r1        = (data.get("r1")        or "").strip()
    r2        = (data.get("r2")        or "").strip()
    bed       = (data.get("bed")       or "").strip() or REF["default_bed"]

    errors = []
    if not sample_id:
        errors.append("Sample ID is required")
    if not r1 or not Path(r1).exists():
        errors.append(f"R1 FASTQ not found: {r1}")
    if not r2 or not Path(r2).exists():
        errors.append(f"R2 FASTQ not found: {r2}")
    if not Path(bed).exists():
        errors.append(f"BED file not found: {bed}")
    if errors:
        return jsonify({"ok": False, "errors": errors}), 400

    job_id = str(uuid.uuid4())
    with db_connect() as conn:
        conn.execute(
            "INSERT INTO jobs(id,sample_id,r1,r2,bed,status,created_at) VALUES(?,?,?,?,?,?,?)",
            (job_id, sample_id, r1, r2, bed, "queued", datetime.now().isoformat())
        )

    t = threading.Thread(target=run_pipeline,
                         args=(job_id, sample_id, r1, r2, bed), daemon=True)
    t.start()
    return jsonify({"ok": True, "job_id": job_id})

@app.route("/api/job/<job_id>")
def job_status(job_id):
    j = job_get(job_id)
    if not j:
        return jsonify({"error": "not found"}), 404
    offset = int(request.args.get("log_offset", 0))
    return jsonify({
        "id":          j["id"],
        "sample_id":   j["sample_id"],
        "status":      j["status"],
        "step":        j["step"],
        "progress":    j["progress"],
        "created_at":  j["created_at"],
        "finished_at": j["finished_at"],
        "error_msg":   j["error_msg"],
        "has_report":  bool(j["report_path"] and Path(j["report_path"]).exists()),
        "has_vcf":     bool(j["vcf_path"]    and Path(j["vcf_path"]).exists()),
        "logs":        logs_get(job_id, offset),
    })

@app.route("/api/jobs")
def list_jobs():
    with db_connect() as conn:
        rows = conn.execute(
            "SELECT id,sample_id,status,step,progress,created_at,finished_at,report_path,vcf_path "
            "FROM jobs ORDER BY created_at DESC LIMIT 50"
        ).fetchall()
    return jsonify([dict(r) for r in rows])

@app.route("/report/<job_id>")
def view_report(job_id):
    j = job_get(job_id)
    if not j or not j["report_path"]:
        abort(404)
    p = Path(j["report_path"])
    if not p.exists():
        abort(404)
    return send_file(str(p), mimetype="text/html")

@app.route("/api/download/<job_id>/<key>")
def download(job_id, key):
    j = job_get(job_id)
    if not j:
        abort(404)
    if key == "report":
        p = j.get("report_path")
    elif key == "vcf":
        p = j.get("vcf_path")
    else:
        abort(400)
    if not p or not Path(p).exists():
        abort(404)
    return send_file(str(p), as_attachment=True)

@app.route("/api/tools")
def tools_check():
    tools = ["fastp", "bwa", "samtools", "gatk", "bcftools",
             "mosdepth", "vep", "bgzip", "tabix"]
    result = {}
    for t in tools:
        ret = subprocess.run(["which", t], capture_output=True)
        result[t] = "ok" if ret.returncode == 0 else "missing"
    return jsonify(result)

# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    db_init()
    print("=" * 60)
    print("  GenRichi Somatic Panel — Web App")
    print("  http://localhost:5000")
    print("=" * 60)
    app.run(host="0.0.0.0", port=5000, debug=False, threaded=True)
