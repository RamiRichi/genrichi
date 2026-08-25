"""
GenRichi Multi-Panel Patient Dashboard

Reads results from all available GenRichi pipelines for one patient and
produces a single, self-contained HTML page.

Missing pipeline results are handled gracefully — sections are hidden rather
than raising errors, so the dashboard works even when only one or two
pipelines have run for a given patient.
"""

import base64
import csv
import io
import json
import math
import os
import re
from datetime import date
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

sm = snakemake  # type: ignore[name-defined]

# ── Parameters ────────────────────────────────────────────────────────────────
PATIENT_ID      = sm.params.patient_id
SEX             = sm.params.sex or "Unknown"
TUMOR_TYPE      = sm.params.tumor_type or "Unknown"
RESULTS_DIR     = sm.params.results_dir
COMPANY         = sm.params.company
LOGO_PATH       = sm.params.logo

HOTSPOT_SAMPLE  = str(sm.params.hotspot_sample  or "").strip()
GERMLINE_SAMPLE = str(sm.params.germline_sample or "").strip()
SOMATIC_SAMPLE  = str(sm.params.somatic_sample  or "").strip()

TMB_HIGH     = float(sm.params.tmb_high)
TMB_MB       = float(sm.params.tmb_coding_mb)
MSI_HIGH     = float(sm.params.msi_high)
AMP_THR      = float(sm.params.cnv_amp)
DEL_THR      = float(sm.params.cnv_del)
MIN_PROBES   = int(sm.params.min_probes)

OUT_HTML     = sm.output.html
REPORT_DATE  = date.today().strftime("%Y-%m-%d")

os.makedirs(os.path.dirname(OUT_HTML), exist_ok=True)


# ── File helpers ──────────────────────────────────────────────────────────────
def _p(sample: str, *parts) -> str:
    return str(Path(RESULTS_DIR) / sample / Path(*parts))


def _read_tsv(path: str) -> list[dict]:
    if not os.path.isfile(path):
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _embed_png(path: str) -> str:
    if not os.path.isfile(path):
        return ""
    with open(path, "rb") as fh:
        return "data:image/png;base64," + base64.b64encode(fh.read()).decode()


def _fig_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=110)
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def _safe_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


# ── Load Phase 1 — Hotspot ────────────────────────────────────────────────────
hotspot_rows = []
if HOTSPOT_SAMPLE:
    raw = _read_tsv(_p(HOTSPOT_SAMPLE, "annotation", f"{HOTSPOT_SAMPLE}.variants.tsv"))
    hotspot_rows = [
        r for r in raw
        if r.get("impact") in ("HIGH", "MODERATE") and r.get("filter") == "PASS"
    ]

# ── Load Phase 2 — Germline ACMG ─────────────────────────────────────────────
germline_rows = []
if GERMLINE_SAMPLE:
    raw = _read_tsv(_p(GERMLINE_SAMPLE, "annotation", f"{GERMLINE_SAMPLE}.acmg_classified.tsv"))
    # Keep Class 3 (VUS), 4 (LP), 5 (P) — drop benign
    germline_rows = [r for r in raw if r.get("acmg_code", "0") in ("3", "4", "5")]

n_pathogenic = sum(1 for r in germline_rows if r.get("acmg_code") == "5")
n_lp         = sum(1 for r in germline_rows if r.get("acmg_code") == "4")
n_vus        = sum(1 for r in germline_rows if r.get("acmg_code") == "3")

# ── Load Phase 3 — Somatic Comprehensive ─────────────────────────────────────
somatic_rows = []
tmb_value    = 0.0
tmb_status   = "N/A"
msi_score    = None
msi_status   = "N/A"
cnv_calls    = []
cnv_b64      = ""

CODING_CSQ = {
    "missense_variant", "stop_gained", "frameshift_variant",
    "splice_donor_variant", "splice_acceptor_variant",
    "start_lost", "stop_lost", "inframe_insertion", "inframe_deletion",
}

if SOMATIC_SAMPLE:
    # Variants
    all_somatic = _read_tsv(_p(SOMATIC_SAMPLE, "annotation", f"{SOMATIC_SAMPLE}.somatic_variants.tsv"))
    somatic_rows = [r for r in all_somatic if r.get("impact") in ("HIGH", "MODERATE")]

    # TMB — count coding somatic variants / Mb
    coding_vars = [
        r for r in all_somatic
        if any(c in r.get("consequence", "").lower() for c in CODING_CSQ)
    ]
    tmb_value  = round(len(coding_vars) / TMB_MB, 2)
    tmb_status = "TMB-High" if tmb_value >= TMB_HIGH else "TMB-Low"

    # MSI — read .msi file (tab-separated, header + data row)
    msi_file = _p(SOMATIC_SAMPLE, "msi", f"{SOMATIC_SAMPLE}.msi")
    if os.path.isfile(msi_file):
        with open(msi_file) as fh:
            lines = [l.strip() for l in fh if l.strip() and not l.startswith("#")]
        for ln in lines:
            parts = ln.split("\t")
            if len(parts) >= 3:
                try:
                    msi_score  = float(parts[2])
                    msi_status = "MSI-H" if msi_score >= MSI_HIGH else "MSS"
                except ValueError:
                    pass

    # CNV
    cnv_raw = _read_tsv(_p(SOMATIC_SAMPLE, "cnv", f"{SOMATIC_SAMPLE}.call.cns"))
    cnv_calls = [
        r for r in cnv_raw
        if r.get("type") in ("Amplification", "Deletion")
        and int(r.get("probes", 0) or 0) >= MIN_PROBES
    ]

    # CNV scatter PNG
    cnv_b64 = _embed_png(_p(SOMATIC_SAMPLE, "cnv", f"{SOMATIC_SAMPLE}-scatter.png"))

n_amp = sum(1 for c in cnv_calls if c.get("type") == "Amplification")
n_del = sum(1 for c in cnv_calls if c.get("type") == "Deletion")

msi_display = f"{msi_score:.1f}/Mb" if msi_score is not None else "N/A"


# ── Alert banner ──────────────────────────────────────────────────────────────
alerts = []
if n_pathogenic > 0:
    alerts.append(f"⚠ {n_pathogenic} Pathogenic germline variant(s)")
if n_lp > 0:
    alerts.append(f"⚠ {n_lp} Likely Pathogenic germline variant(s)")
if tmb_status == "TMB-High":
    alerts.append(f"⚠ TMB-High ({tmb_value} mut/Mb)")
if msi_status == "MSI-H":
    alerts.append(f"⚠ MSI-High ({msi_display})")
if n_amp > 0:
    alerts.append(f"⚠ {n_amp} amplification(s)")
if n_del > 0:
    alerts.append(f"⚠ {n_del} deletion(s)")
if not alerts:
    alerts = ["✓ No high-priority findings across all tested panels"]

alert_color = "#e74c3c" if any("⚠" in a for a in alerts) else "#27ae60"


# ── Charts ────────────────────────────────────────────────────────────────────
def _tmb_chart() -> str:
    fig, ax = plt.subplots(figsize=(4, 2.8))
    color = "#e74c3c" if tmb_value >= TMB_HIGH else "#3498db"
    ax.barh(["TMB"], [tmb_value], color=color, height=0.4)
    ax.axvline(TMB_HIGH, color="#e74c3c", lw=1.5, linestyle="--", alpha=0.7)
    ax.set_xlim(0, max(TMB_HIGH * 3, tmb_value * 1.4, 1))
    ax.text(TMB_HIGH + 0.3, 0, f"High ≥{TMB_HIGH:.0f}", va="center", fontsize=7,
            color="#e74c3c")
    ax.set_xlabel("Mutations per Mb", fontsize=8)
    ax.set_title(f"TMB = {tmb_value} [{tmb_status}]", fontsize=9, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    return _fig_b64(fig)


def _msi_chart() -> str:
    fig, ax = plt.subplots(figsize=(4, 2.8))
    score = msi_score if msi_score is not None else 0.0
    color = "#e74c3c" if msi_status == "MSI-H" else "#3498db"
    ax.barh(["MSI"], [score], color=color, height=0.4)
    ax.axvline(MSI_HIGH, color="#e74c3c", lw=1.5, linestyle="--", alpha=0.7)
    ax.set_xlim(0, max(MSI_HIGH * 3, score * 1.4, 1))
    ax.text(MSI_HIGH + 0.3, 0, f"High ≥{MSI_HIGH:.0f}", va="center", fontsize=7,
            color="#e74c3c")
    ax.set_xlabel("Indels per Mb", fontsize=8)
    title = f"MSI = {msi_display} [{msi_status}]" if msi_score is not None else "MSI = N/A"
    ax.set_title(title, fontsize=9, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    return _fig_b64(fig)


def _germline_donut() -> str:
    """Donut chart of ACMG classes."""
    data, labels, colors_ = [], [], []
    if n_pathogenic:
        data.append(n_pathogenic); labels.append(f"Pathogenic ({n_pathogenic})")
        colors_.append("#e74c3c")
    if n_lp:
        data.append(n_lp); labels.append(f"Likely Path. ({n_lp})")
        colors_.append("#e67e22")
    if n_vus:
        data.append(n_vus); labels.append(f"VUS ({n_vus})")
        colors_.append("#f1c40f")
    if not data:
        data = [1]; labels = ["No P/LP/VUS"]; colors_ = ["#2ecc71"]

    fig, ax = plt.subplots(figsize=(3.5, 3.0))
    wedges, _ = ax.pie(
        data, labels=None, colors=colors_,
        startangle=90,
        wedgeprops={"width": 0.5, "edgecolor": "white", "linewidth": 1.5},
    )
    ax.legend(wedges, labels, loc="lower center", fontsize=7,
              bbox_to_anchor=(0.5, -0.15), ncol=1)
    ax.set_title("Germline Findings", fontsize=9, fontweight="bold")
    fig.tight_layout()
    return _fig_b64(fig)


tmb_img      = _tmb_chart()      if SOMATIC_SAMPLE  else ""
msi_img      = _msi_chart()      if SOMATIC_SAMPLE  else ""
germline_img = _germline_donut() if GERMLINE_SAMPLE  else ""


# ── HTML helpers ──────────────────────────────────────────────────────────────
ACMG_CLASS_STYLE = {
    "5": "background:#e74c3c;color:#fff;padding:2px 7px;border-radius:4px;font-weight:bold",
    "4": "background:#e67e22;color:#fff;padding:2px 7px;border-radius:4px;font-weight:bold",
    "3": "background:#f1c40f;color:#333;padding:2px 7px;border-radius:4px;font-weight:bold",
    "2": "background:#95a5a6;color:#fff;padding:2px 7px;border-radius:4px",
    "1": "background:#2ecc71;color:#fff;padding:2px 7px;border-radius:4px",
}
ACMG_CLASS_LABEL = {
    "5": "Pathogenic",
    "4": "Likely Pathogenic",
    "3": "VUS",
    "2": "Likely Benign",
    "1": "Benign",
}
IMPACT_STYLE = {
    "HIGH":     "background:#e74c3c;color:#fff;padding:2px 7px;border-radius:4px;font-size:11px",
    "MODERATE": "background:#e67e22;color:#fff;padding:2px 7px;border-radius:4px;font-size:11px",
}
CNV_TYPE_STYLE = {
    "Amplification": "background:#e74c3c;color:#fff;padding:2px 8px;border-radius:4px",
    "Deletion":      "background:#3498db;color:#fff;padding:2px 8px;border-radius:4px",
}


def _badge(text: str, style: str) -> str:
    return f'<span style="{style}">{text}</span>'


def _img_tag(b64: str, alt: str = "", style: str = "max-width:100%;height:auto") -> str:
    if not b64:
        return f'<p style="color:#aaa;text-align:center;padding:20px">{alt} — not available</p>'
    return f'<img src="{b64}" alt="{alt}" style="{style}">'


def _section(title: str, icon: str, content: str, show: bool = True) -> str:
    if not show:
        return ""
    return f"""
    <div class="section">
      <div class="section-header">
        <span class="section-icon">{icon}</span> {title}
      </div>
      <div class="section-body">
        {content}
      </div>
    </div>"""


# ── Build tables ──────────────────────────────────────────────────────────────
def _germline_table(rows: list[dict]) -> str:
    if not rows:
        return '<p class="no-data">No P/LP/VUS variants detected in the germline panel.</p>'
    html = """<table><thead><tr>
        <th>Gene</th><th>HGVSc</th><th>HGVSp</th>
        <th>Consequence</th><th>Zygosity</th>
        <th>gnomAD AF</th><th>ACMG Class</th><th>Criteria</th>
    </tr></thead><tbody>"""
    for r in rows:
        code  = r.get("acmg_code", "3")
        label = ACMG_CLASS_LABEL.get(code, "VUS")
        style = ACMG_CLASS_STYLE.get(code, "")
        html += f"""<tr>
            <td><strong>{r.get('gene','')}</strong></td>
            <td class="mono">{r.get('hgvsc','')}</td>
            <td class="mono">{r.get('hgvsp','')}</td>
            <td>{r.get('consequence','').replace('_',' ')}</td>
            <td>{r.get('zygosity','')}</td>
            <td>{r.get('gnomad_af','')}</td>
            <td>{_badge(label, style)}</td>
            <td class="mono small">{r.get('acmg_criteria','')}</td>
        </tr>"""
    html += "</tbody></table>"
    return html


def _somatic_table(rows: list[dict], title_label: str = "somatic") -> str:
    if not rows:
        return f'<p class="no-data">No HIGH or MODERATE impact {title_label} variants detected.</p>'
    html = """<table><thead><tr>
        <th>Gene</th><th>HGVSc</th><th>HGVSp</th>
        <th>Consequence</th><th>Impact</th>
        <th>Depth</th><th>VAF</th>
        <th>gnomAD AF</th><th>ClinVar</th><th>COSMIC</th>
    </tr></thead><tbody>"""
    for r in rows:
        imp   = r.get("impact", "")
        style = IMPACT_STYLE.get(imp, "")
        html += f"""<tr>
            <td><strong>{r.get('gene','')}</strong></td>
            <td class="mono">{r.get('hgvsc','')}</td>
            <td class="mono">{r.get('hgvsp','')}</td>
            <td>{r.get('consequence','').replace('_',' ')}</td>
            <td>{_badge(imp, style)}</td>
            <td>{r.get('depth','')}</td>
            <td>{r.get('vaf','')}</td>
            <td>{r.get('gnomad_af','')}</td>
            <td>{r.get('clinvar_sig','')}</td>
            <td class="mono small">{r.get('cosmic_id','')}</td>
        </tr>"""
    html += "</tbody></table>"
    return html


def _cnv_table(rows: list[dict]) -> str:
    if not rows:
        return '<p class="no-data">No significant CNV calls (≥' + str(MIN_PROBES) + ' probes).</p>'
    html = """<table><thead><tr>
        <th>Gene / Region</th><th>Chromosome</th><th>Start</th><th>End</th>
        <th>log2 Ratio</th><th>Copy #</th><th>Probes</th><th>Type</th>
    </tr></thead><tbody>"""
    for r in rows:
        t     = r.get("type", "")
        style = CNV_TYPE_STYLE.get(t, "")
        start = f"{int(r.get('start', 0)):,}"
        end_  = f"{int(r.get('end', 0)):,}"
        html += f"""<tr>
            <td><strong>{r.get('gene','')}</strong></td>
            <td>{r.get('chromosome','')}</td>
            <td class="mono">{start}</td>
            <td class="mono">{end_}</td>
            <td class="mono">{r.get('log2','')}</td>
            <td>{r.get('cn','')}</td>
            <td>{r.get('probes','')}</td>
            <td>{_badge(t, style)}</td>
        </tr>"""
    html += "</tbody></table>"
    return html


# ── Pipeline availability pills ───────────────────────────────────────────────
def _pill(label: str, active: bool) -> str:
    bg = "#27ae60" if active else "#bdc3c7"
    return f'<span style="background:{bg};color:#fff;padding:3px 10px;border-radius:12px;font-size:11px;margin-right:6px">{label}</span>'


pipeline_pills = (
    _pill("Phase 1 · Hotspot",      bool(HOTSPOT_SAMPLE)) +
    _pill("Phase 2 · Germline",     bool(GERMLINE_SAMPLE)) +
    _pill("Phase 3 · Comprehensive", bool(SOMATIC_SAMPLE))
)


# ── Summary cards ─────────────────────────────────────────────────────────────
def _card(title: str, value: str, sub: str, color: str, available: bool = True) -> str:
    if not available:
        value = "—"
        sub   = "not run"
        color = "#bdc3c7"
    return f"""
    <div class="card">
      <div class="card-title">{title}</div>
      <div class="card-value" style="color:{color}">{value}</div>
      <div class="card-sub">{sub}</div>
    </div>"""


germline_card_color = (
    "#e74c3c" if n_pathogenic > 0 else
    "#e67e22" if n_lp > 0 else
    "#f1c40f" if n_vus > 0 else "#27ae60"
)
germline_card_value = f"{n_pathogenic + n_lp + n_vus}"
germline_card_sub   = f"{n_pathogenic} P · {n_lp} LP · {n_vus} VUS"

tmb_card_color = "#e74c3c" if tmb_status == "TMB-High" else "#27ae60"
msi_card_color = "#e74c3c" if msi_status == "MSI-H"    else "#27ae60"
cnv_card_color = "#e74c3c" if (n_amp + n_del) > 0      else "#27ae60"

cards_html = (
    _card("GERMLINE RISK",    germline_card_value,        germline_card_sub,
          germline_card_color, bool(GERMLINE_SAMPLE)) +
    _card("TUMOR MUTATIONAL BURDEN", f"{tmb_value} mut/Mb", tmb_status,
          tmb_card_color, bool(SOMATIC_SAMPLE)) +
    _card("MICROSATELLITE INSTABILITY", msi_display, msi_status,
          msi_card_color, bool(SOMATIC_SAMPLE)) +
    _card("COPY NUMBER VARIANTS", str(n_amp + n_del),
          f"{n_amp} amp · {n_del} del", cnv_card_color, bool(SOMATIC_SAMPLE))
)


# ── Alert banner HTML ─────────────────────────────────────────────────────────
alert_items = "".join(f"<li>{a}</li>" for a in alerts)
alert_html  = f"""
<div style="background:{alert_color};color:#fff;padding:12px 20px;
            border-radius:6px;margin-bottom:20px;font-size:13px">
  <ul style="margin:0;padding-left:18px">{alert_items}</ul>
</div>"""


# ── Charts section ────────────────────────────────────────────────────────────
charts_content = ""
if SOMATIC_SAMPLE or GERMLINE_SAMPLE:
    parts = []
    if germline_img:
        parts.append(f'<div class="chart-box">{_img_tag(germline_img,"Germline Summary")}</div>')
    if tmb_img:
        parts.append(f'<div class="chart-box">{_img_tag(tmb_img,"TMB")}</div>')
    if msi_img:
        parts.append(f'<div class="chart-box">{_img_tag(msi_img,"MSI")}</div>')
    charts_content = '<div class="chart-row">' + "".join(parts) + "</div>"

cnv_content = ""
if cnv_b64:
    cnv_content = _img_tag(cnv_b64, "CNV Genome-Wide Scatter",
                           "max-width:100%;height:auto;border-radius:6px")
elif SOMATIC_SAMPLE:
    cnv_content = '<p class="no-data">CNV scatter plot not available.</p>'


# ── Assemble HTML ─────────────────────────────────────────────────────────────
HTML = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{COMPANY} Multi-Panel Dashboard — {PATIENT_ID}</title>
<style>
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ font-family: "Segoe UI", Arial, sans-serif; background: #f0f2f5;
          color: #2c3e50; font-size: 13px; }}

  /* Header */
  .header {{ background: linear-gradient(135deg, #1a252f 0%, #2c3e50 100%);
             color: #fff; padding: 18px 28px; display: flex;
             justify-content: space-between; align-items: center; }}
  .header-left h1 {{ font-size: 20px; font-weight: 700; letter-spacing: .5px; }}
  .header-left p  {{ font-size: 11px; color: #95a5a6; margin-top: 2px; }}
  .header-right   {{ text-align: right; font-size: 12px; color: #bdc3c7; line-height: 1.8; }}
  .header-right strong {{ color: #ecf0f1; font-size: 14px; }}

  /* Patient meta bar */
  .meta-bar {{ background: #2c3e50; color: #ecf0f1; padding: 8px 28px;
               display: flex; gap: 32px; font-size: 12px; }}
  .meta-bar span {{ opacity: 0.75; }}
  .meta-bar strong {{ opacity: 1; }}

  /* Pipelines pill bar */
  .pills-bar {{ background: #34495e; padding: 6px 28px; }}

  /* Main */
  .main {{ max-width: 1300px; margin: 0 auto; padding: 20px 24px; }}

  /* Cards */
  .cards {{ display: grid; grid-template-columns: repeat(4, 1fr); gap: 14px;
            margin-bottom: 20px; }}
  @media (max-width: 900px) {{ .cards {{ grid-template-columns: repeat(2, 1fr); }} }}
  .card {{ background: #fff; border-radius: 8px; padding: 16px;
           box-shadow: 0 1px 4px rgba(0,0,0,.08); text-align: center; }}
  .card-title {{ font-size: 10px; text-transform: uppercase; letter-spacing: .8px;
                 color: #7f8c8d; margin-bottom: 6px; }}
  .card-value {{ font-size: 28px; font-weight: 700; margin-bottom: 2px; }}
  .card-sub   {{ font-size: 11px; color: #7f8c8d; }}

  /* Sections */
  .section {{ background: #fff; border-radius: 8px; margin-bottom: 16px;
              box-shadow: 0 1px 4px rgba(0,0,0,.08); overflow: hidden; }}
  .section-header {{ background: #2c3e50; color: #ecf0f1; padding: 10px 18px;
                     font-size: 13px; font-weight: 600; }}
  .section-icon {{ margin-right: 6px; }}
  .section-body {{ padding: 16px 18px; }}

  /* Tables */
  table {{ width: 100%; border-collapse: collapse; font-size: 12px; }}
  th {{ background: #ecf0f1; color: #2c3e50; font-weight: 600; padding: 7px 10px;
        text-align: left; font-size: 11px; text-transform: uppercase;
        letter-spacing: .5px; border-bottom: 2px solid #bdc3c7; }}
  td {{ padding: 7px 10px; border-bottom: 1px solid #f0f2f5; vertical-align: middle; }}
  tr:last-child td {{ border-bottom: none; }}
  tr:hover td {{ background: #f8f9fa; }}
  .mono  {{ font-family: "Consolas","Courier New",monospace; }}
  .small {{ font-size: 11px; }}

  /* Charts */
  .chart-row {{ display: flex; gap: 16px; flex-wrap: wrap; }}
  .chart-box {{ flex: 1; min-width: 200px; text-align: center; }}

  .no-data {{ color: #95a5a6; font-style: italic; padding: 12px 0; text-align: center; }}

  /* Footer */
  .footer {{ text-align: center; padding: 16px; color: #95a5a6; font-size: 11px; }}
</style>
</head>
<body>

<div class="header">
  <div class="header-left">
    <h1>🧬 {COMPANY} — Multi-Panel Patient Dashboard</h1>
    <p>Integrates Hotspot · Germline Hereditary · Somatic Comprehensive results</p>
  </div>
  <div class="header-right">
    Report Date: <strong>{REPORT_DATE}</strong><br>
    Pipeline: GenRichi v1 (Phase 4)
  </div>
</div>

<div class="meta-bar">
  <div><span>Patient ID</span><br><strong>{PATIENT_ID}</strong></div>
  <div><span>Sex</span><br><strong>{SEX}</strong></div>
  <div><span>Tumor Type</span><br><strong>{TUMOR_TYPE}</strong></div>
</div>

<div class="pills-bar">
  Panels run:&nbsp; {pipeline_pills}
</div>

<div class="main">

  {alert_html}

  <div class="cards">{cards_html}</div>

  {_section("Germline Cancer Risk", "🔬",
      _germline_table(germline_rows),
      show=bool(GERMLINE_SAMPLE)
  )}

  {_section("Somatic Driver Variants — Comprehensive Panel", "🧬",
      _somatic_table(somatic_rows, "somatic"),
      show=bool(SOMATIC_SAMPLE)
  )}

  {_section("Hotspot Panel — Driver Calls", "🎯",
      _somatic_table(hotspot_rows, "hotspot"),
      show=bool(HOTSPOT_SAMPLE)
  )}

  {_section("Copy Number Variants", "📊",
      _cnv_table(cnv_calls) + (
          f'<div style="margin-top:16px">{cnv_content}</div>'
          if cnv_b64 else ""
      ),
      show=bool(SOMATIC_SAMPLE)
  )}

  {_section("Biomarker Charts", "📈",
      charts_content,
      show=bool(charts_content)
  )}

</div>

<div class="footer">
  Generated by {COMPANY} GenRichi Pipeline · {REPORT_DATE} ·
  For research use only — not for clinical diagnosis without validation
</div>

</body>
</html>"""

with open(OUT_HTML, "w") as fh:
    fh.write(HTML)

with open(sm.log[0], "w") as fh:
    fh.write(f"Dashboard generated: {OUT_HTML}\n")
    fh.write(f"Germline: {GERMLINE_SAMPLE or 'not run'} — {len(germline_rows)} P/LP/VUS\n")
    fh.write(f"Somatic:  {SOMATIC_SAMPLE  or 'not run'} — TMB={tmb_value}, MSI={msi_display}\n")
    fh.write(f"Hotspot:  {HOTSPOT_SAMPLE  or 'not run'} — {len(hotspot_rows)} high/moderate\n")
    fh.write(f"CNV calls: {n_amp} amp, {n_del} del\n")
