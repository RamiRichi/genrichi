"""
GenRichi Somatic Comprehensive Report
Combines: SNV/indel table + CNV calls + MSI score + TMB calculation
into a single self-contained HTML report.
"""

import base64
import csv
import io
import json
import math
import os
import re
from datetime import date

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd

sm = snakemake  # type: ignore[name-defined]

# ── Parameters ────────────────────────────────────────────────────────────────
SAMPLE_ID      = sm.params.sample_id
PATIENT_ID     = sm.params.patient_id
SEX            = sm.params.sex
TUMOR_TYPE     = sm.params.tumor_type
PANEL_NAME     = sm.params.panel_name
COMPANY        = sm.params.company
LOGO_PATH      = sm.params.logo
SHOW_SYN       = sm.params.show_synonymous
MSI_THRESHOLD  = float(sm.params.msi_threshold)
TMB_MB         = float(sm.params.tmb_coding_mb)
TMB_HIGH       = float(sm.params.tmb_high_threshold)
AMP_THR        = float(sm.params.cnv_amp_threshold)
DEL_THR        = float(sm.params.cnv_del_threshold)


# ── Helpers ───────────────────────────────────────────────────────────────────
def _safe_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _img_to_b64(path: str) -> str:
    try:
        with open(path, "rb") as f:
            return base64.b64encode(f.read()).decode()
    except Exception:
        return ""


def _fig_to_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=120)
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


# ── Load variants ─────────────────────────────────────────────────────────────
variants_df = pd.read_csv(sm.input.variants, sep="\t", dtype=str).fillna("")
if not SHOW_SYN:
    variants_df = variants_df[
        ~variants_df["consequence"].str.contains("synonymous", case=False, na=False)
    ]

# ── TMB calculation ───────────────────────────────────────────────────────────
coding_consequences = {
    "missense_variant", "stop_gained", "frameshift_variant",
    "splice_donor_variant", "splice_acceptor_variant",
    "start_lost", "stop_lost", "inframe_insertion", "inframe_deletion",
}
tmb_variants = variants_df[
    variants_df["consequence"].apply(
        lambda c: any(x in c.lower() for x in coding_consequences)
    )
]
tmb_value = round(len(tmb_variants) / TMB_MB, 2)
tmb_status = "TMB-High" if tmb_value >= TMB_HIGH else "TMB-Low"

# ── MSI score ─────────────────────────────────────────────────────────────────
msi_score = None
msi_status = "Unknown"
try:
    with open(sm.input.msi_score) as f:
        lines = f.readlines()
    # msisensor-pro output: last line has "Total_Number_of_Sites Somatic% ..."
    for line in lines:
        if not line.startswith("#") and line.strip():
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                try:
                    msi_score = float(parts[2])
                    msi_status = "MSI-High" if msi_score >= MSI_THRESHOLD else "MSS"
                except ValueError:
                    pass
except Exception:
    pass

# ── CNV calls ─────────────────────────────────────────────────────────────────
cnv_calls = []
try:
    cnv_df = pd.read_csv(sm.input.call_cns, sep="\t", comment="#")
    for _, row in cnv_df.iterrows():
        log2 = _safe_float(str(row.get("log2", "")))
        if log2 is None:
            continue
        if log2 >= AMP_THR:
            call_type = "Amplification"
        elif log2 <= DEL_THR:
            call_type = "Deletion"
        else:
            continue
        gene = str(row.get("gene", str(row.get("chromosome", "."))))
        cnv_calls.append({
            "gene": gene,
            "chrom": str(row.get("chromosome", ".")),
            "start": int(row.get("start", 0)),
            "end": int(row.get("end", 0)),
            "log2": round(log2, 3),
            "cn": str(row.get("cn", ".")),
            "type": call_type,
        })
except Exception:
    pass

# ── QC metrics ────────────────────────────────────────────────────────────────
def _parse_fastp(path):
    try:
        with open(path) as f:
            d = json.load(f)
        s = d.get("summary", {})
        br = s.get("before_filtering", {})
        q30 = round(br.get("q30_rate", 0) * 100, 1)
        return {
            "total_reads": br.get("total_reads", 0),
            "q30": q30,
        }
    except Exception:
        return {"total_reads": 0, "q30": 0}


def _parse_flagstat(path):
    try:
        with open(path) as f:
            txt = f.read()
        m = re.search(r"([\d,]+) \+ \d+ mapped \((.+?)%", txt)
        if m:
            return int(m.group(1).replace(",", "")), float(m.group(2))
        return 0, 0.0
    except Exception:
        return 0, 0.0


def _parse_mosdepth(path):
    try:
        df = pd.read_csv(path, sep="\t")
        row = df[df["chrom"] == "total_region"]
        if not row.empty:
            return round(float(row["mean"].values[0]), 1)
        row = df[df["chrom"] == "total"]
        if not row.empty:
            return round(float(row["mean"].values[0]), 1)
        return 0.0
    except Exception:
        return 0.0


def _parse_markdup(path):
    try:
        with open(path) as f:
            lines = f.readlines()
        for i, line in enumerate(lines):
            if line.startswith("ESTIMATED") or "PERCENT_DUPLICATION" in line:
                if i + 1 < len(lines):
                    vals = lines[i + 1].strip().split("\t")
                    headers = line.strip().split("\t")
                    d = dict(zip(headers, vals))
                    return round(float(d.get("PERCENT_DUPLICATION", 0)) * 100, 2)
        return 0.0
    except Exception:
        return 0.0


t_fastp  = _parse_fastp(sm.input.tumor_fastp)
n_fastp  = _parse_fastp(sm.input.normal_fastp)
t_mapped, t_map_pct = _parse_flagstat(sm.input.tumor_flagstat)
n_mapped, n_map_pct = _parse_flagstat(sm.input.normal_flagstat)
t_depth  = _parse_mosdepth(sm.input.tumor_mosdepth)
n_depth  = _parse_mosdepth(sm.input.normal_mosdepth)
t_dup    = _parse_markdup(sm.input.tumor_markdup)
n_dup    = _parse_markdup(sm.input.normal_markdup)

# ── Variant impact colours ────────────────────────────────────────────────────
IMPACT_COLOR = {
    "HIGH":     "#e74c3c",
    "MODERATE": "#e67e22",
    "LOW":      "#3498db",
    "MODIFIER": "#95a5a6",
}

# ── Charts ────────────────────────────────────────────────────────────────────
def _make_tmb_gauge(tmb_val, tmb_high):
    fig, ax = plt.subplots(figsize=(4, 2.5))
    ax.set_xlim(0, 50)
    ax.set_ylim(0, 1)
    ax.axvspan(0, tmb_high, alpha=0.15, color="#3498db")
    ax.axvspan(tmb_high, 50, alpha=0.15, color="#e74c3c")
    ax.axvline(tmb_val, color="#2c3e50", lw=3)
    ax.axvline(tmb_high, color="#e74c3c", lw=1.5, linestyle="--")
    ax.set_yticks([])
    ax.set_xlabel("TMB (mutations/Mb)", fontsize=10)
    ax.set_title(f"TMB = {tmb_val}  [{('TMB-High' if tmb_val >= tmb_high else 'TMB-Low')}]",
                 fontsize=11, fontweight="bold")
    ax.text(tmb_high + 0.5, 0.5, f"High ≥{tmb_high}", color="#e74c3c", fontsize=8, va="center")
    fig.tight_layout()
    return fig


def _make_msi_bar(score, threshold):
    fig, ax = plt.subplots(figsize=(4, 2))
    color = "#e74c3c" if (score or 0) >= threshold else "#2ecc71"
    ax.barh(["MSI"], [score or 0], color=color, height=0.4)
    ax.axvline(threshold, color="#e74c3c", lw=1.5, linestyle="--")
    ax.set_xlim(0, max(100, (score or 0) * 1.2))
    ax.set_xlabel("% Unstable sites", fontsize=10)
    status = "MSI-High" if (score or 0) >= threshold else "MSS"
    ax.set_title(f"MSI Score = {score:.1f}%  [{status}]" if score else "MSI Score = N/A",
                 fontsize=11, fontweight="bold")
    fig.tight_layout()
    return fig


def _make_qc_chart(t_depth, n_depth):
    fig, ax = plt.subplots(figsize=(5, 3))
    labels = ["Tumor", "Normal"]
    depths = [t_depth, n_depth]
    colors = ["#e74c3c", "#3498db"]
    bars = ax.bar(labels, depths, color=colors, width=0.5)
    ax.axhline(100, color="#95a5a6", lw=1, linestyle="--", label="100× target")
    for bar, val in zip(bars, depths):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 2,
                f"{val}×", ha="center", fontsize=10, fontweight="bold")
    ax.set_ylabel("Mean Panel Depth (×)", fontsize=10)
    ax.set_title("Panel Coverage", fontsize=11, fontweight="bold")
    ax.legend(fontsize=8)
    fig.tight_layout()
    return fig


tmb_gauge_b64  = _fig_to_b64(_make_tmb_gauge(tmb_value, TMB_HIGH))
msi_bar_b64    = _fig_to_b64(_make_msi_bar(msi_score, MSI_THRESHOLD))
qc_chart_b64   = _fig_to_b64(_make_qc_chart(t_depth, n_depth))
cnv_scatter_b64 = _img_to_b64(sm.input.cnv_scatter)
logo_b64       = _img_to_b64(LOGO_PATH)

# ── Variant rows HTML ─────────────────────────────────────────────────────────
def _variant_rows(df):
    if df.empty:
        return "<tr><td colspan='11' style='text-align:center;color:#999'>No variants in this category</td></tr>"
    rows = []
    for _, r in df.iterrows():
        impact = r.get("impact", "")
        color  = IMPACT_COLOR.get(impact, "#555")
        cosmic = r.get("cosmic_id", "")
        cosmic_link = (
            f'<a href="https://cancer.sanger.ac.uk/cosmic/mutation/overview?id={cosmic.replace("CDS","")}" '
            f'target="_blank">{cosmic}</a>'
            if cosmic else ""
        )
        vaf = r.get("vaf", "")
        try:
            vaf = f"{float(vaf)*100:.1f}%"
        except Exception:
            pass
        rows.append(
            f"<tr>"
            f"<td><b>{r.get('gene','')}</b></td>"
            f"<td>{r.get('hgvsc','')}</td>"
            f"<td>{r.get('hgvsp','')}</td>"
            f"<td><span style='color:{color}'>{r.get('consequence','')}</span></td>"
            f"<td><b style='color:{color}'>{impact}</b></td>"
            f"<td>{r.get('depth','')}</td>"
            f"<td>{vaf}</td>"
            f"<td>{r.get('gnomad_af','')}</td>"
            f"<td>{r.get('clinvar_sig','')}</td>"
            f"<td>{cosmic_link}</td>"
            f"</tr>"
        )
    return "\n".join(rows)


high_vars = variants_df[variants_df["impact"] == "HIGH"]
mod_vars  = variants_df[variants_df["impact"] == "MODERATE"]
other_vars = variants_df[~variants_df["impact"].isin(["HIGH", "MODERATE"])]

# ── CNV rows HTML ─────────────────────────────────────────────────────────────
def _cnv_rows(calls):
    if not calls:
        return "<tr><td colspan='6' style='text-align:center;color:#999'>No significant CNV calls</td></tr>"
    rows = []
    for c in calls:
        color = "#e74c3c" if c["type"] == "Amplification" else "#3498db"
        rows.append(
            f"<tr>"
            f"<td><b>{c['gene']}</b></td>"
            f"<td>{c['chrom']}:{c['start']:,}–{c['end']:,}</td>"
            f"<td>{c['log2']}</td>"
            f"<td>{c['cn']}</td>"
            f"<td><b style='color:{color}'>{c['type']}</b></td>"
            f"</tr>"
        )
    return "\n".join(rows)


# ── Render HTML ───────────────────────────────────────────────────────────────
logo_html = (
    f'<img src="data:image/png;base64,{logo_b64}" style="height:50px;vertical-align:middle;margin-right:12px">'
    if logo_b64 else ""
)

cnv_scatter_html = (
    f'<img src="data:image/png;base64,{cnv_scatter_b64}" style="max-width:100%;border-radius:8px">'
    if cnv_scatter_b64 else "<p style='color:#999'>CNV scatter plot not available</p>"
)

# Clinical summary banner
n_high = len(high_vars)
n_mod  = len(mod_vars)
summary_color = "#c0392b" if n_high > 0 else ("#e67e22" if n_mod > 0 else "#27ae60")
summary_text  = (
    f"{n_high} HIGH impact variant{'s' if n_high!=1 else ''} detected"
    if n_high > 0 else
    f"{n_mod} MODERATE impact variant{'s' if n_mod!=1 else ''} detected"
    if n_mod > 0 else
    "No HIGH or MODERATE impact somatic variants detected"
)

HTML = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{COMPANY} — Comprehensive Cancer Panel Report — {SAMPLE_ID}</title>
<style>
  body{{font-family:Arial,sans-serif;margin:0;background:#f5f5f5;color:#333}}
  .header{{background:#2c3e50;color:white;padding:20px 32px;display:flex;align-items:center}}
  .header h1{{margin:0;font-size:1.4em}}
  .header .sub{{font-size:0.85em;opacity:0.8;margin-top:4px}}
  .meta-grid{{display:grid;grid-template-columns:repeat(3,1fr);gap:0;background:white;border-bottom:2px solid #ecf0f1}}
  .meta-item{{padding:10px 20px;border-right:1px solid #ecf0f1}}
  .meta-item:last-child{{border-right:none}}
  .meta-item .label{{font-size:0.72em;color:#999;text-transform:uppercase;letter-spacing:.5px}}
  .meta-item .value{{font-size:0.95em;font-weight:bold;color:#2c3e50}}
  .banner{{padding:14px 32px;font-size:1em;font-weight:bold;color:white;background:{summary_color}}}
  .biomarker-row{{display:grid;grid-template-columns:repeat(3,1fr);gap:16px;padding:20px 32px;background:white;border-bottom:2px solid #ecf0f1}}
  .biomarker-card{{background:#f8f9fa;border-radius:8px;padding:16px;text-align:center}}
  .biomarker-card .bm-label{{font-size:0.8em;color:#777;text-transform:uppercase;letter-spacing:.5px;margin-bottom:6px}}
  .biomarker-card .bm-value{{font-size:1.6em;font-weight:bold}}
  .tmb-high{{color:#e74c3c}}.tmb-low{{color:#27ae60}}
  .msi-high{{color:#e74c3c}}.mss{{color:#27ae60}}
  .section{{background:white;margin:16px 32px;border-radius:8px;box-shadow:0 1px 3px rgba(0,0,0,.08)}}
  .section-header{{background:#34495e;color:white;padding:12px 20px;border-radius:8px 8px 0 0;font-weight:bold;font-size:0.95em}}
  .section-body{{padding:16px 20px}}
  table{{width:100%;border-collapse:collapse;font-size:0.85em}}
  th{{background:#ecf0f1;padding:8px 10px;text-align:left;font-size:0.8em;color:#555;text-transform:uppercase;letter-spacing:.3px}}
  td{{padding:7px 10px;border-bottom:1px solid #f0f0f0}}
  tr:last-child td{{border-bottom:none}}
  tr:hover td{{background:#fafafa}}
  .qc-grid{{display:grid;grid-template-columns:1fr 1fr;gap:16px}}
  .qc-table th,qc-table td{{font-size:0.82em}}
  .chart-row{{display:grid;grid-template-columns:repeat(3,1fr);gap:16px;padding:0 32px 20px}}
  .chart-card{{background:white;border-radius:8px;box-shadow:0 1px 3px rgba(0,0,0,.08);padding:12px;text-align:center}}
  .cnv-scatter{{background:white;margin:0 32px 20px;border-radius:8px;box-shadow:0 1px 3px rgba(0,0,0,.08);padding:16px}}
  .footer{{text-align:center;font-size:0.75em;color:#aaa;padding:24px;margin-top:8px}}
  .tag-high{{background:#fde8e6;color:#c0392b;padding:2px 8px;border-radius:4px;font-size:0.8em;font-weight:bold}}
  .tag-moderate{{background:#fef3e2;color:#d68910;padding:2px 8px;border-radius:4px;font-size:0.8em;font-weight:bold}}
</style>
</head>
<body>

<div class="header">
  {logo_html}
  <div>
    <h1>{COMPANY} — Comprehensive Cancer Panel Report</h1>
    <div class="sub">Clinical genomics reporting system &bull; Somatic variant analysis</div>
  </div>
</div>

<div class="meta-grid">
  <div class="meta-item"><div class="label">Sample ID</div><div class="value">{SAMPLE_ID}</div></div>
  <div class="meta-item"><div class="label">Patient ID</div><div class="value">{PATIENT_ID}</div></div>
  <div class="meta-item"><div class="label">Sex</div><div class="value">{SEX}</div></div>
  <div class="meta-item"><div class="label">Tumor Type</div><div class="value">{TUMOR_TYPE}</div></div>
  <div class="meta-item"><div class="label">Panel</div><div class="value">{PANEL_NAME}</div></div>
  <div class="meta-item"><div class="label">Report Date</div><div class="value">{date.today().isoformat()}</div></div>
</div>

<div class="banner">&#128203; {summary_text} &nbsp;|&nbsp; TMB: {tmb_value} mut/Mb [{tmb_status}] &nbsp;|&nbsp; MSI: {f"{msi_score:.1f}%" if msi_score is not None else "N/A"} [{msi_status}]</div>

<!-- Biomarker summary cards -->
<div class="biomarker-row">
  <div class="biomarker-card">
    <div class="bm-label">Tumor Mutational Burden</div>
    <div class="bm-value {'tmb-high' if tmb_value >= TMB_HIGH else 'tmb-low'}">{tmb_value} <span style="font-size:0.5em">mut/Mb</span></div>
    <div style="font-size:0.85em;margin-top:4px;color:{'#e74c3c' if tmb_value >= TMB_HIGH else '#27ae60'};font-weight:bold">{tmb_status}</div>
    <div style="font-size:0.75em;color:#999;margin-top:2px">Threshold: &ge;{TMB_HIGH} mut/Mb</div>
  </div>
  <div class="biomarker-card">
    <div class="bm-label">Microsatellite Instability</div>
    <div class="bm-value {'msi-high' if msi_status=='MSI-High' else 'mss'}">{f"{msi_score:.1f}%" if msi_score is not None else "N/A"}</div>
    <div style="font-size:0.85em;margin-top:4px;color:{'#e74c3c' if msi_status=='MSI-High' else '#27ae60'};font-weight:bold">{msi_status}</div>
    <div style="font-size:0.75em;color:#999;margin-top:2px">Threshold: &ge;{MSI_THRESHOLD}% unstable</div>
  </div>
  <div class="biomarker-card">
    <div class="bm-label">CNV Calls</div>
    <div class="bm-value" style="color:#8e44ad">{len(cnv_calls)}</div>
    <div style="font-size:0.85em;margin-top:4px;color:#555">
      {sum(1 for c in cnv_calls if c['type']=='Amplification')} amp &nbsp; {sum(1 for c in cnv_calls if c['type']=='Deletion')} del
    </div>
    <div style="font-size:0.75em;color:#999;margin-top:2px">Amp log2 &ge;{AMP_THR} &nbsp; Del log2 &le;{DEL_THR}</div>
  </div>
</div>

<!-- Charts row -->
<div class="chart-row">
  <div class="chart-card">
    <img src="data:image/png;base64,{tmb_gauge_b64}" style="max-width:100%">
  </div>
  <div class="chart-card">
    <img src="data:image/png;base64,{msi_bar_b64}" style="max-width:100%">
  </div>
  <div class="chart-card">
    <img src="data:image/png;base64,{qc_chart_b64}" style="max-width:100%">
  </div>
</div>

<!-- HIGH impact variants -->
<div class="section">
  <div class="section-header">&#128308; HIGH Impact Somatic Variants ({len(high_vars)})</div>
  <div class="section-body">
    <table>
      <tr><th>Gene</th><th>HGVSc</th><th>HGVSp</th><th>Consequence</th><th>Impact</th><th>Depth</th><th>VAF</th><th>gnomAD AF</th><th>ClinVar</th><th>COSMIC</th></tr>
      {_variant_rows(high_vars)}
    </table>
  </div>
</div>

<!-- MODERATE impact variants -->
<div class="section">
  <div class="section-header">&#128992; MODERATE Impact Somatic Variants ({len(mod_vars)})</div>
  <div class="section-body">
    <table>
      <tr><th>Gene</th><th>HGVSc</th><th>HGVSp</th><th>Consequence</th><th>Impact</th><th>Depth</th><th>VAF</th><th>gnomAD AF</th><th>ClinVar</th><th>COSMIC</th></tr>
      {_variant_rows(mod_vars)}
    </table>
  </div>
</div>

<!-- CNV calls -->
<div class="section">
  <div class="section-header">&#128200; Copy Number Variants ({len(cnv_calls)})</div>
  <div class="section-body">
    <table>
      <tr><th>Gene/Region</th><th>Location</th><th>log2 Ratio</th><th>Copy Number</th><th>Type</th></tr>
      {_cnv_rows(cnv_calls)}
    </table>
  </div>
</div>

<!-- CNV scatter -->
<div class="cnv-scatter">
  <div style="font-weight:bold;color:#34495e;margin-bottom:10px">&#128202; CNV Genome-Wide Scatter Plot</div>
  {cnv_scatter_html}
</div>

<!-- QC metrics -->
<div class="section">
  <div class="section-header">&#128202; Quality Control Metrics</div>
  <div class="section-body">
    <div class="qc-grid">
      <div>
        <b>Tumor</b>
        <table style="margin-top:8px">
          <tr><th>Metric</th><th>Value</th></tr>
          <tr><td>Total reads</td><td>{t_fastp['total_reads']:,}</td></tr>
          <tr><td>Q30 rate</td><td>{t_fastp['q30']} %</td></tr>
          <tr><td>Mapped reads</td><td>{t_mapped:,} ({t_map_pct:.1f} %)</td></tr>
          <tr><td>Duplication rate</td><td>{t_dup} %</td></tr>
          <tr><td>Mean panel depth</td><td>{t_depth} &times;</td></tr>
        </table>
      </div>
      <div>
        <b>Normal</b>
        <table style="margin-top:8px">
          <tr><th>Metric</th><th>Value</th></tr>
          <tr><td>Total reads</td><td>{n_fastp['total_reads']:,}</td></tr>
          <tr><td>Q30 rate</td><td>{n_fastp['q30']} %</td></tr>
          <tr><td>Mapped reads</td><td>{n_mapped:,} ({n_map_pct:.1f} %)</td></tr>
          <tr><td>Duplication rate</td><td>{n_dup} %</td></tr>
          <tr><td>Mean panel depth</td><td>{n_depth} &times;</td></tr>
        </table>
      </div>
    </div>
  </div>
</div>

<div class="footer">
  <b>{COMPANY}</b> &bull; {PANEL_NAME} &bull; Report generated {date.today().isoformat()}<br>
  <span style="color:#c0392b">FOR RESEARCH USE ONLY. Not validated for clinical diagnostic use. Requires expert review before clinical reporting.</span>
</div>

</body></html>"""

os.makedirs(os.path.dirname(sm.output.html), exist_ok=True)
with open(sm.output.html, "w", encoding="utf-8") as fh:
    fh.write(HTML)
