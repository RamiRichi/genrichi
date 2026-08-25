"""
Generate a self-contained HTML clinical report for a hereditary cancer panel sample.

Inputs (via snakemake.input):
  variants          – ACMG-classified TSV (from acmg_classifier.py)
  fastp_json        – fastp QC JSON
  flagstat          – samtools flagstat output
  mosdepth_summary  – mosdepth summary file
  mosdepth_regions  – mosdepth per-region BED.gz
  markdup_metrics   – Picard MarkDuplicates metrics

Params (via snakemake.params):
  sample_id, panel_name, company, logo, show_benign,
  patient_id, sex, indication
"""

import base64
import csv
import gzip
import io
import json
import re
from datetime import date
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd

sm = snakemake  # type: ignore[name-defined]

SAMPLE_ID  = sm.params.sample_id
PANEL_NAME = sm.params.panel_name
COMPANY    = sm.params.company
LOGO_PATH  = sm.params.logo
SHOW_BENIGN = bool(sm.params.show_benign)
PATIENT_ID = sm.params.patient_id
SEX        = sm.params.sex
INDICATION = sm.params.indication
RUN_DATE   = date.today().isoformat()

ACMG_COLORS = {
    "Pathogenic":       "#c0392b",
    "Likely_Pathogenic":"#e67e22",
    "VUS":              "#f39c12",
    "Likely_Benign":    "#27ae60",
    "Benign":           "#2ecc71",
}
ACMG_LABELS = {
    "Pathogenic":        "Class 5 — Pathogenic",
    "Likely_Pathogenic": "Class 4 — Likely Pathogenic",
    "VUS":               "Class 3 — Variant of Uncertain Significance",
    "Likely_Benign":     "Class 2 — Likely Benign",
    "Benign":            "Class 1 — Benign",
}
DISPLAY_CLASSES = ["Pathogenic", "Likely_Pathogenic", "VUS"]
if SHOW_BENIGN:
    DISPLAY_CLASSES += ["Likely_Benign", "Benign"]


# ── Helpers ──────────────────────────────────────────────────────────────────

def fig_to_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def logo_b64(path: str | None) -> str | None:
    if not path or not Path(path).exists():
        return None
    ext  = Path(path).suffix.lstrip(".").lower()
    mime = {"png": "image/png", "jpg": "image/jpeg", "jpeg": "image/jpeg",
            "svg": "image/svg+xml"}.get(ext, "image/png")
    data = base64.b64encode(Path(path).read_bytes()).decode()
    return f"data:{mime};base64,{data}"


# ── Parse QC inputs ──────────────────────────────────────────────────────────

def parse_fastp(path: str) -> dict:
    with open(path) as f:
        j = json.load(f)
    filt = j.get("filtering_result", {})
    summ = j.get("summary", {})
    r1   = summ.get("before_filtering", {})
    return {
        "total_reads":     r1.get("total_reads", 0),
        "q30_rate":        round(r1.get("q30_rate", 0) * 100, 1),
        "passed_reads":    filt.get("passed_filter_reads", 0),
        "low_quality":     filt.get("low_quality_reads", 0),
    }


def parse_flagstat(path: str) -> dict:
    result = {"mapped_pct": "N/A", "mapped_reads": 0}
    with open(path) as f:
        for line in f:
            m = re.search(r"(\d+) \+ \d+ mapped \(([^%]+)%", line)
            if m:
                result["mapped_reads"] = int(m.group(1))
                result["mapped_pct"]   = f"{float(m.group(2)):.1f}"
    return result


def parse_mosdepth_summary(path: str) -> dict:
    result = {"mean_depth": 0.0, "pct_20x": "N/A"}
    with open(path) as f:
        for line in f:
            if line.startswith("total\t"):
                parts = line.split("\t")
                result["mean_depth"] = float(parts[3]) if len(parts) > 3 else 0.0
    return result


def parse_mosdepth_regions(path: str) -> pd.DataFrame:
    opener = gzip.open if path.endswith(".gz") else open
    rows = []
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.strip().split("\t")
            if len(p) >= 4:
                rows.append({"gene": p[3], "depth": float(p[4]) if len(p) > 4 else 0.0})
    if not rows:
        return pd.DataFrame(columns=["gene", "depth"])
    df = pd.DataFrame(rows).groupby("gene", as_index=False)["depth"].mean()
    df["depth"] = df["depth"].round(1)
    return df.sort_values("depth")


def parse_markdup(path: str) -> dict:
    with open(path) as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if line.startswith("ESTIMATED_LIBRARY_SIZE") or line.startswith("LIBRARY"):
            try:
                vals  = lines[i + 1].strip().split("\t")
                heads = line.strip().split("\t")
                d     = dict(zip(heads, vals))
                pct_dup = float(d.get("PERCENT_DUPLICATION", 0)) * 100
                return {"pct_duplication": round(pct_dup, 2)}
            except Exception:
                pass
    return {"pct_duplication": "N/A"}


# ── Charts ────────────────────────────────────────────────────────────────────

def make_coverage_chart(regions_df: pd.DataFrame) -> str:
    if regions_df.empty:
        return ""
    fig, ax = plt.subplots(figsize=(10, max(4, len(regions_df) * 0.3)))
    colors = ["#e74c3c" if d < 20 else "#f39c12" if d < 50 else "#27ae60"
              for d in regions_df["depth"]]
    ax.barh(regions_df["gene"], regions_df["depth"], color=colors)
    ax.axvline(20, color="red",    linestyle="--", linewidth=0.8, label="20x threshold")
    ax.axvline(50, color="orange", linestyle="--", linewidth=0.8, label="50x threshold")
    ax.set_xlabel("Mean Depth (×)", fontsize=10)
    ax.set_title("Gene Coverage", fontsize=12, fontweight="bold")
    ax.legend(fontsize=8)
    plt.tight_layout()
    return fig_to_b64(fig)


def make_acmg_pie(df: pd.DataFrame) -> str:
    counts = df["acmg_class"].value_counts()
    if counts.empty:
        return ""
    labels = [ACMG_LABELS.get(k, k) for k in counts.index]
    colors = [ACMG_COLORS.get(k, "#95a5a6") for k in counts.index]
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.pie(counts.values, labels=labels, colors=colors, autopct="%1.0f%%",
           startangle=140, textprops={"fontsize": 8})
    ax.set_title("ACMG Classification", fontweight="bold")
    plt.tight_layout()
    return fig_to_b64(fig)


# ── HTML helpers ─────────────────────────────────────────────────────────────

def _acmg_badge(cls: str) -> str:
    color = ACMG_COLORS.get(cls, "#95a5a6")
    label = cls.replace("_", " ")
    return (
        f'<span style="background:{color};color:#fff;padding:2px 8px;'
        f'border-radius:4px;font-size:0.8em;font-weight:bold">{label}</span>'
    )


def _gnomad_fmt(val: str) -> str:
    try:
        f = float(val)
        if f == 0:
            return "absent"
        if f < 0.0001:
            return f"{f:.2e}"
        return f"{f:.4f}"
    except (ValueError, TypeError):
        return val or "."


def _variant_table(rows: list[dict]) -> str:
    if not rows:
        return "<p><em>No variants in this category.</em></p>"
    cols = [
        ("gene", "Gene"), ("zygosity", "Zygosity"), ("hgvsc", "HGVSc"),
        ("hgvsp", "HGVSp"), ("consequence", "Consequence"),
        ("depth", "Depth"), ("gq", "GQ"), ("gnomad_af", "gnomAD AF"),
        ("clinvar_sig", "ClinVar"), ("acmg_criteria", "ACMG Criteria"),
    ]
    th = "".join(f'<th style="padding:6px 10px;text-align:left;border-bottom:2px solid #ddd">{h}</th>'
                 for _, h in cols)
    body = ""
    for i, r in enumerate(rows):
        bg = "#f9f9f9" if i % 2 == 0 else "#fff"
        td_parts = []
        for key, _ in cols:
            val = r.get(key, "")
            if key == "gnomad_af":
                val = _gnomad_fmt(str(val))
            elif key == "consequence":
                val = str(val).replace(",", ", ")
            td_parts.append(
                f'<td style="padding:5px 10px;font-size:0.85em">{val}</td>'
            )
        body += f'<tr style="background:{bg}">{"".join(td_parts)}</tr>\n'
    return (
        f'<div style="overflow-x:auto"><table style="width:100%;border-collapse:collapse">'
        f'<thead><tr>{th}</tr></thead><tbody>{body}</tbody></table></div>'
    )


def _section(title: str, cls: str, rows: list[dict]) -> str:
    color  = ACMG_COLORS.get(cls, "#95a5a6")
    label  = ACMG_LABELS.get(cls, cls)
    count  = len(rows)
    badge  = _acmg_badge(cls)
    table  = _variant_table(rows)
    return f"""
<div style="margin:20px 0;border-left:5px solid {color};padding-left:15px">
  <h3 style="color:{color};margin:0 0 6px">{label} ({count} variant{"s" if count != 1 else ""}) {badge}</h3>
  {table}
</div>
"""


# ── Build HTML ────────────────────────────────────────────────────────────────

def build_html(variants_df: pd.DataFrame, qc: dict, regions_df: pd.DataFrame) -> str:
    logo_tag = ""
    logo_src = logo_b64(LOGO_PATH)
    if logo_src:
        logo_tag = f'<img src="{logo_src}" style="height:50px;float:right">'

    # QC table
    qc_rows = [
        ("Total reads (before trim)",    f"{qc['fastp']['total_reads']:,}"),
        ("Reads passing QC",             f"{qc['fastp']['passed_reads']:,}"),
        ("Q30 rate",                     f"{qc['fastp']['q30_rate']} %"),
        ("Mapped reads",                 f"{qc['flagstat']['mapped_reads']:,} ({qc['flagstat']['mapped_pct']} %)"),
        ("Duplication rate",             f"{qc['markdup']['pct_duplication']} %"),
        ("Mean panel depth",             f"{qc['mosdepth']['mean_depth']} ×"),
    ]
    qc_table_html = "".join(
        f'<tr style="background:{"#f9f9f9" if i%2==0 else "#fff"}">'
        f'<td style="padding:6px 12px;font-weight:bold">{k}</td>'
        f'<td style="padding:6px 12px">{v}</td></tr>'
        for i, (k, v) in enumerate(qc_rows)
    )

    # Clinical sections
    clinical_html = ""
    for cls in DISPLAY_CLASSES:
        subset = variants_df[variants_df["acmg_class"] == cls].to_dict("records")
        clinical_html += _section(ACMG_LABELS.get(cls, cls), cls, subset)

    # Charts
    coverage_b64 = make_coverage_chart(regions_df)
    pie_b64      = make_acmg_pie(variants_df)

    coverage_img = (
        f'<img src="data:image/png;base64,{coverage_b64}" style="max-width:100%">'
        if coverage_b64 else "<p><em>No region coverage data.</em></p>"
    )
    pie_img = (
        f'<img src="data:image/png;base64,{pie_b64}" style="max-width:400px">'
        if pie_b64 else ""
    )

    # Variant count summary
    total_classified = len(variants_df)
    p_lp_count       = len(variants_df[variants_df["acmg_class"].isin(["Pathogenic", "Likely_Pathogenic"])])
    vus_count         = len(variants_df[variants_df["acmg_class"] == "VUS"])

    alert_style = (
        'style="background:#fdecea;border:1px solid #c0392b;border-radius:6px;'
        'padding:12px 18px;margin:10px 0"'
    )
    summary_banner = f"""
<div {alert_style}>
  <strong>Clinical Summary:</strong>
  &nbsp; {p_lp_count} Pathogenic/Likely Pathogenic variant{"s" if p_lp_count != 1 else ""}
  &nbsp;|&nbsp; {vus_count} VUS
  &nbsp;|&nbsp; {total_classified} total classified
</div>
""" if p_lp_count > 0 else f"""
<div style="background:#eafaf1;border:1px solid #27ae60;border-radius:6px;padding:12px 18px;margin:10px 0">
  <strong>Clinical Summary:</strong>
  No Pathogenic or Likely Pathogenic variants detected.
  &nbsp;|&nbsp; {vus_count} VUS &nbsp;|&nbsp; {total_classified} total classified
</div>
"""

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{COMPANY} Hereditary Report — {SAMPLE_ID}</title>
<style>
  body  {{ font-family: Arial, sans-serif; font-size: 14px; color: #333; margin: 0; padding: 0; }}
  .page {{ max-width: 1200px; margin: auto; padding: 30px; }}
  h1    {{ color: #2c3e50; border-bottom: 3px solid #2c3e50; padding-bottom: 8px; }}
  h2    {{ color: #2c3e50; margin-top: 30px; }}
  table {{ border-collapse: collapse; width: 100%; }}
  th    {{ background: #2c3e50; color: #fff; padding: 8px 12px; text-align: left; }}
  .footer {{ font-size: 0.78em; color: #777; margin-top: 40px; border-top: 1px solid #ddd; padding-top: 10px; }}
</style>
</head>
<body>
<div class="page">

<!-- Header -->
<div style="overflow:hidden">
  {logo_tag}
  <h1>{COMPANY} — Hereditary Cancer Panel Report</h1>
</div>

<table style="width:100%;margin-bottom:20px">
  <tr>
    <td><strong>Sample ID:</strong> {SAMPLE_ID}</td>
    <td><strong>Patient ID:</strong> {PATIENT_ID}</td>
    <td><strong>Sex:</strong> {SEX}</td>
  </tr>
  <tr>
    <td><strong>Panel:</strong> {PANEL_NAME}</td>
    <td><strong>Indication:</strong> {INDICATION}</td>
    <td><strong>Report date:</strong> {RUN_DATE}</td>
  </tr>
</table>

<!-- Summary Banner -->
<h2>Clinical Summary</h2>
{summary_banner}

<!-- QC Section -->
<h2>Quality Control Metrics</h2>
<table style="width:60%">
  <thead><tr><th>Metric</th><th>Value</th></tr></thead>
  <tbody>{qc_table_html}</tbody>
</table>

<!-- Clinical Findings -->
<h2>Clinical Findings</h2>
<p style="font-size:0.9em;color:#666">
  Variants classified per ACMG/AMP 2015 guidelines (Richards et al., PMID 25741868).
  Classification is based on available computational evidence (ClinVar, gnomAD, VEP consequence, SIFT/PolyPhen)
  and requires expert review before clinical reporting.
</p>
{clinical_html}

<!-- Charts -->
<h2>Classification Overview</h2>
<div style="display:flex;gap:30px;flex-wrap:wrap;align-items:flex-start">
  <div style="flex:1;min-width:300px">{pie_img}</div>
</div>

<h2>Gene Coverage</h2>
{coverage_img}

<!-- ACMG Legend -->
<h2>ACMG Criteria Legend</h2>
<table style="width:auto;font-size:0.85em">
  <tr><th>Code</th><th>Strength</th><th>Description</th></tr>
  <tr><td>PVS1</td><td>Very Strong</td><td>Predicted null in LOF-intolerant gene (stop/frameshift/splice)</td></tr>
  <tr style="background:#f9f9f9"><td>PS1</td><td>Strong</td><td>ClinVar: same amino-acid change classified Pathogenic</td></tr>
  <tr><td>PM1</td><td>Moderate</td><td>Missense in well-established functional domain</td></tr>
  <tr style="background:#f9f9f9"><td>PM2</td><td>Moderate</td><td>Absent or rare in gnomAD (AF &lt; 0.1 %)</td></tr>
  <tr><td>PP3</td><td>Supporting</td><td>SIFT deleterious + PolyPhen probably_damaging</td></tr>
  <tr style="background:#f9f9f9"><td>PP5</td><td>Supporting</td><td>ClinVar classified Pathogenic/LP in reputable source</td></tr>
  <tr><td>BA1</td><td>Benign Stand-Alone</td><td>gnomAD AF &gt; 5 %</td></tr>
  <tr style="background:#f9f9f9"><td>BS1/BS2</td><td>Benign Strong</td><td>Common in gnomAD / ClinVar Benign in reputable source</td></tr>
  <tr><td>BP4/BP6</td><td>Benign Supporting</td><td>SIFT tolerated + PolyPhen benign / ClinVar Benign/LB</td></tr>
</table>

<!-- Disclaimer -->
<div class="footer">
  <p><strong>Disclaimer:</strong>
  This report is generated by an automated bioinformatics pipeline and is intended for research and informational
  purposes only. ACMG classifications are based on computational evidence and must be reviewed by a qualified
  clinical geneticist or genetic counselor before clinical use. Variants not detected by this panel (including
  large structural rearrangements and deep intronic variants) are not reported.
  Pipeline: {COMPANY} Hereditary Cancer Panel v1 · Reference: GRCh38 · Generated: {RUN_DATE}
  </p>
</div>

</div>
</body>
</html>
"""


# ── Main ─────────────────────────────────────────────────────────────────────

variants_df = pd.read_csv(sm.input.variants, sep="\t", dtype=str)

qc = {
    "fastp":    parse_fastp(sm.input.fastp_json),
    "flagstat": parse_flagstat(sm.input.flagstat),
    "mosdepth": parse_mosdepth_summary(sm.input.mosdepth_summary),
    "markdup":  parse_markdup(sm.input.markdup_metrics),
}

regions_df = parse_mosdepth_regions(sm.input.mosdepth_regions)

html = build_html(variants_df, qc, regions_df)

Path(sm.output.html).write_text(html, encoding="utf-8")
