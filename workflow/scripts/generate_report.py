"""
Generate a self-contained HTML clinical report for one GenRichi sample.

Inputs (via snakemake.input):
  variants          – TSV from vcf_to_table.py
  fastp_json        – fastp QC JSON
  flagstat          – samtools flagstat output
  mosdepth_summary  – mosdepth summary file
  mosdepth_regions  – mosdepth per-region BED.gz
  contamination     – GATK CalculateContamination table
  markdup_metrics   – Picard MarkDuplicates metrics

Params (via snakemake.params):
  sample_id, panel_name, company, min_vaf, logo
"""

import base64
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

# ── Snakemake bindings ──────────────────────────────────────────────────────
sm = snakemake  # type: ignore[name-defined]
SAMPLE_ID   = sm.params.sample_id
PANEL_NAME  = sm.params.panel_name
COMPANY     = sm.params.company
MIN_VAF     = float(sm.params.min_vaf)
LOGO_PATH   = sm.params.logo
RUN_DATE    = date.today().isoformat()


# ── Helpers ─────────────────────────────────────────────────────────────────

def fig_to_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def logo_b64(path: str | None) -> str | None:
    if not path or not Path(path).exists():
        return None
    ext = Path(path).suffix.lstrip(".").lower()
    mime = {"png": "image/png", "jpg": "image/jpeg", "jpeg": "image/jpeg",
            "svg": "image/svg+xml"}.get(ext, "image/png")
    data = base64.b64encode(Path(path).read_bytes()).decode()
    return f"data:{mime};base64,{data}"


# ── Parse inputs ─────────────────────────────────────────────────────────────

def parse_fastp(path: str) -> dict:
    with open(path) as fh:
        d = json.load(fh)
    s = d.get("summary", {})
    before = s.get("before_filtering", {})
    after  = s.get("after_filtering", {})
    dup    = d.get("duplication", {})
    return {
        "total_reads_raw":     before.get("total_reads", 0),
        "total_bases_raw":     before.get("total_bases", 0),
        "q30_rate_raw":        round(before.get("q30_rate", 0) * 100, 2),
        "gc_content_raw":      round(before.get("gc_content", 0) * 100, 2),
        "total_reads_trimmed": after.get("total_reads", 0),
        "total_bases_trimmed": after.get("total_bases", 0),
        "q30_rate_trimmed":    round(after.get("q30_rate", 0) * 100, 2),
        "duplication_rate":    round(dup.get("rate", 0) * 100, 2),
    }


def parse_flagstat(path: str) -> dict:
    d: dict = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            m = re.match(r"^(\d+) \+ \d+ (.+)$", line)
            if m:
                d[m.group(2)] = int(m.group(1))
    total    = d.get("in total (QC-passed reads + QC-failed reads)", 0)
    mapped   = d.get("mapped", d.get("primary mapped", 0))
    properly = d.get("properly paired", 0)
    pct_map  = round(mapped / total * 100, 2) if total else 0
    pct_pp   = round(properly / total * 100, 2) if total else 0
    return {
        "total_reads":    total,
        "mapped_reads":   mapped,
        "pct_mapped":     pct_map,
        "properly_paired": properly,
        "pct_properly_paired": pct_pp,
    }


def parse_mosdepth_summary(path: str) -> dict:
    df = pd.read_csv(path, sep="\t")
    # "total_region" row has whole-panel stats
    row = df[df["chrom"] == "total_region"]
    if row.empty:
        row = df.tail(1)
    mean_depth = float(row["mean"].iloc[0]) if "mean" in row.columns else 0.0
    return {"mean_depth": round(mean_depth, 1)}


def parse_mosdepth_regions(path: str) -> pd.DataFrame:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        df = pd.read_csv(fh, sep="\t", header=None,
                         names=["chrom", "start", "end", "region", "mean_depth"])
    return df


def parse_contamination(path: str) -> float:
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            try:
                return round(float(parts[1]) * 100, 3)  # contamination fraction → %
            except (IndexError, ValueError):
                pass
    return 0.0


def parse_markdup(path: str) -> dict:
    with open(path) as fh:
        lines = fh.readlines()
    # Picard metrics block starts after two header lines
    for i, line in enumerate(lines):
        if line.startswith("ESTIMATED_LIBRARY_SIZE") or line.startswith("PERCENT_DUPLICATION"):
            header = lines[i].strip().split("\t")
            values = lines[i + 1].strip().split("\t") if i + 1 < len(lines) else []
            d = dict(zip(header, values))
            pct = float(d.get("PERCENT_DUPLICATION", 0)) * 100
            return {"pct_duplication": round(pct, 2),
                    "estimated_library_size": int(d.get("ESTIMATED_LIBRARY_SIZE", 0))}
    return {"pct_duplication": 0.0, "estimated_library_size": 0}


def load_variants(path: str, min_vaf: float) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    if df.empty:
        return df
    df["vaf"]       = pd.to_numeric(df["vaf"], errors="coerce").fillna(0)
    df["depth"]     = pd.to_numeric(df["depth"], errors="coerce").fillna(0).astype(int)
    df["alt_reads"] = pd.to_numeric(df["alt_reads"], errors="coerce").fillna(0).astype(int)
    return df[df["vaf"] >= min_vaf].reset_index(drop=True)


# ── Chart generators ─────────────────────────────────────────────────────────

IMPACT_COLORS = {
    "HIGH":     "#e74c3c",
    "MODERATE": "#f39c12",
    "LOW":      "#27ae60",
    "MODIFIER": "#95a5a6",
}


def chart_vaf_histogram(variants: pd.DataFrame) -> str:
    fig, ax = plt.subplots(figsize=(5, 3))
    if variants.empty or "vaf" not in variants.columns:
        ax.text(0.5, 0.5, "No variants", ha="center", va="center")
    else:
        ax.hist(variants["vaf"] * 100, bins=20, color="#2980b9", edgecolor="white",
                linewidth=0.5)
        ax.set_xlabel("Variant Allele Frequency (%)")
        ax.set_ylabel("Count")
        ax.set_title("VAF Distribution")
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig_to_b64(fig)


def chart_impact_pie(variants: pd.DataFrame) -> str:
    fig, ax = plt.subplots(figsize=(4, 4))
    if variants.empty or "impact" not in variants.columns:
        ax.text(0.5, 0.5, "No variants", ha="center", va="center")
    else:
        counts = variants["impact"].value_counts()
        colors = [IMPACT_COLORS.get(k, "#bdc3c7") for k in counts.index]
        wedges, texts, autotexts = ax.pie(
            counts, labels=counts.index, colors=colors,
            autopct="%1.0f%%", startangle=90,
            wedgeprops={"linewidth": 1, "edgecolor": "white"},
        )
        for t in autotexts:
            t.set_fontsize(9)
        ax.set_title("Variant Impact")
    fig.tight_layout()
    return fig_to_b64(fig)


def chart_coverage_bar(regions_df: pd.DataFrame) -> str:
    fig, ax = plt.subplots(figsize=(max(6, len(regions_df) * 0.3 + 2), 4))
    if regions_df.empty:
        ax.text(0.5, 0.5, "No coverage data", ha="center", va="center")
    else:
        df = regions_df.copy()
        df["label"] = df["region"].str.replace(r".*_", "", regex=True)
        colors = ["#e74c3c" if d < 30 else "#27ae60" for d in df["mean_depth"]]
        ax.bar(range(len(df)), df["mean_depth"], color=colors, width=0.7)
        ax.axhline(30, color="#e74c3c", linestyle="--", linewidth=1, label="30×")
        ax.axhline(100, color="#f39c12", linestyle="--", linewidth=1, label="100×")
        ax.set_xticks(range(len(df)))
        ax.set_xticklabels(df["label"], rotation=90, fontsize=7)
        ax.set_ylabel("Mean Depth (×)")
        ax.set_title("Coverage per Panel Region")
        ax.legend(fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig_to_b64(fig)


# ── HTML template ─────────────────────────────────────────────────────────────

HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width,initial-scale=1"/>
<title>{company} – Somatic Report – {sample_id}</title>
<style>
  :root {{
    --primary: #1a3a5c;
    --accent:  #2980b9;
    --warn:    #e74c3c;
    --ok:      #27ae60;
    --bg:      #f5f7fa;
  }}
  *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ font-family: 'Segoe UI', Arial, sans-serif; background: var(--bg);
          color: #2c3e50; font-size: 13px; }}
  header {{ background: var(--primary); color: white; padding: 16px 32px;
            display: flex; align-items: center; gap: 20px; }}
  header img {{ height: 48px; }}
  header .title-block h1 {{ font-size: 20px; font-weight: 700; }}
  header .title-block p  {{ font-size: 12px; opacity: 0.8; margin-top: 2px; }}
  .badge {{ background: rgba(255,255,255,0.15); border-radius: 4px;
            padding: 2px 8px; font-size: 11px; }}
  main {{ padding: 24px 32px; max-width: 1400px; margin: 0 auto; }}
  section {{ background: white; border-radius: 8px; box-shadow: 0 1px 4px rgba(0,0,0,.08);
             padding: 20px 24px; margin-bottom: 20px; }}
  section h2 {{ font-size: 15px; font-weight: 600; color: var(--primary);
                border-bottom: 2px solid var(--accent); padding-bottom: 6px;
                margin-bottom: 14px; }}
  .grid-2 {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }}
  .grid-3 {{ display: grid; grid-template-columns: repeat(3, 1fr); gap: 20px; }}
  .kv-table {{ width: 100%; border-collapse: collapse; }}
  .kv-table td {{ padding: 5px 10px; border-bottom: 1px solid #ecf0f1; }}
  .kv-table td:first-child {{ font-weight: 600; color: #555; width: 55%; }}
  .kv-table tr:last-child td {{ border-bottom: none; }}
  table.vt {{ width: 100%; border-collapse: collapse; font-size: 12px; }}
  table.vt thead th {{ background: var(--primary); color: white; padding: 7px 9px;
                       text-align: left; font-weight: 600; position: sticky; top: 0; }}
  table.vt tbody tr:nth-child(even) {{ background: #f8f9fa; }}
  table.vt tbody tr:hover {{ background: #ebf5fb; }}
  table.vt tbody td {{ padding: 6px 9px; border-bottom: 1px solid #ecf0f1;
                       vertical-align: top; }}
  .tag {{ display: inline-block; padding: 1px 6px; border-radius: 3px;
          font-size: 11px; font-weight: 600; }}
  .tag-HIGH     {{ background: #fde8e8; color: var(--warn); }}
  .tag-MODERATE {{ background: #fef6e4; color: #d35400; }}
  .tag-LOW      {{ background: #eafaf1; color: var(--ok); }}
  .tag-MODIFIER {{ background: #f0f0f0; color: #777; }}
  .metric {{ text-align: center; padding: 12px; }}
  .metric .val  {{ font-size: 26px; font-weight: 700; color: var(--accent); }}
  .metric .lbl  {{ font-size: 11px; color: #888; margin-top: 2px; }}
  .charts {{ display: flex; gap: 20px; flex-wrap: wrap; }}
  .charts img   {{ border-radius: 6px; border: 1px solid #eee; max-width: 100%; }}
  .warn-text {{ color: var(--warn); font-weight: 600; }}
  .ok-text   {{ color: var(--ok);   font-weight: 600; }}
  footer {{ text-align: center; font-size: 11px; color: #aaa; padding: 20px;
            border-top: 1px solid #e0e0e0; margin-top: 12px; }}
</style>
</head>
<body>

<header>
  {logo_tag}
  <div class="title-block">
    <h1>{company} — Somatic Hotspot Panel Report</h1>
    <p>{panel_name} &nbsp;|&nbsp; Sample: <strong>{sample_id}</strong>
       &nbsp;|&nbsp; Date: {run_date}</p>
  </div>
  <span class="badge">{variant_count} reportable variant(s)</span>
</header>

<main>

<!-- ── QC Summary ─────────────────────────────────────── -->
<section>
  <h2>Quality Control Summary</h2>
  <div class="grid-3">
    <table class="kv-table">
      <caption style="text-align:left;font-weight:700;margin-bottom:6px;color:#444">
        Sequencing QC (fastp)</caption>
      <tr><td>Raw reads</td><td>{total_reads_raw:,}</td></tr>
      <tr><td>Trimmed reads</td><td>{total_reads_trimmed:,}</td></tr>
      <tr><td>Q30 rate (trimmed)</td><td>{q30_rate_trimmed}%</td></tr>
      <tr><td>GC content</td><td>{gc_content_raw}%</td></tr>
      <tr><td>Duplication rate</td>
          <td class="{dup_cls}">{duplication_rate}%</td></tr>
    </table>
    <table class="kv-table">
      <caption style="text-align:left;font-weight:700;margin-bottom:6px;color:#444">
        Alignment (flagstat)</caption>
      <tr><td>Total reads</td><td>{flagstat_total:,}</td></tr>
      <tr><td>Mapped reads</td>
          <td class="{map_cls}">{mapped_reads:,} ({pct_mapped}%)</td></tr>
      <tr><td>Properly paired</td>
          <td>{properly_paired:,} ({pct_properly_paired}%)</td></tr>
      <tr><td>Estimated lib size</td><td>{estimated_library_size:,}</td></tr>
    </table>
    <table class="kv-table">
      <caption style="text-align:left;font-weight:700;margin-bottom:6px;color:#444">
        Coverage &amp; Contamination</caption>
      <tr><td>Mean panel depth</td>
          <td class="{depth_cls}">{mean_depth}×</td></tr>
      <tr><td>Contamination estimate</td>
          <td class="{cont_cls}">{contamination}%</td></tr>
      <tr><td>Min reportable VAF</td><td>{min_vaf_pct}%</td></tr>
      <tr><td>Panel</td><td>{panel_name_short}</td></tr>
    </table>
  </div>
</section>

<!-- ── Coverage Chart ─────────────────────────────────── -->
<section>
  <h2>Panel Coverage per Region</h2>
  <img src="data:image/png;base64,{chart_coverage}" style="width:100%;max-width:1100px"/>
</section>

<!-- ── Somatic Variants ───────────────────────────────── -->
<section>
  <h2>Reportable Somatic Variants (VAF &ge; {min_vaf_pct}%)</h2>
  {variant_table_html}
</section>

<!-- ── Charts ─────────────────────────────────────────── -->
<section>
  <h2>Variant Analytics</h2>
  <div class="charts">
    <img src="data:image/png;base64,{chart_vaf}" style="height:260px"/>
    <img src="data:image/png;base64,{chart_impact}" style="height:260px"/>
  </div>
</section>

</main>

<footer>
  Generated by {company} Somatic Pipeline &nbsp;|&nbsp; {run_date}
  &nbsp;|&nbsp; This report is for research use only unless validated for clinical use.
</footer>

</body>
</html>
"""


def build_variant_table(df: pd.DataFrame) -> str:
    if df.empty:
        return '<p style="color:#888;font-style:italic">No reportable variants detected.</p>'

    rows_html = []
    for _, row in df.iterrows():
        impact = row.get("impact", "")
        tag_cls = f"tag-{impact}" if impact in IMPACT_COLORS else "tag-MODIFIER"
        cosmic_cell = (
            f'<a href="https://cancer.sanger.ac.uk/cosmic/search?q={row["cosmic_id"]}"'
            f' target="_blank">{row["cosmic_id"]}</a>'
            if row.get("cosmic_id") else "–"
        )
        clnsig = row.get("clinvar_sig", "") or "–"
        gnomad = row.get("gnomad_af", "") or "–"
        rows_html.append(
            f"<tr>"
            f"<td><strong>{row.get('gene','')}</strong></td>"
            f"<td>{row.get('hgvsc','')}</td>"
            f"<td>{row.get('hgvsp','')}</td>"
            f"<td>{row.get('consequence','').replace(', ','<br/>')}</td>"
            f"<td><span class='tag {tag_cls}'>{impact}</span></td>"
            f"<td>{row.get('depth',0)}</td>"
            f"<td>{row.get('alt_reads',0)}</td>"
            f"<td><strong>{row['vaf']*100:.1f}%</strong></td>"
            f"<td>{cosmic_cell}</td>"
            f"<td>{clnsig}</td>"
            f"<td>{gnomad}</td>"
            f"</tr>"
        )

    return (
        "<div style='overflow-x:auto'>"
        "<table class='vt'>"
        "<thead><tr>"
        "<th>Gene</th><th>HGVSc</th><th>HGVSp</th><th>Consequence</th>"
        "<th>Impact</th><th>Depth</th><th>Alt reads</th><th>VAF</th>"
        "<th>COSMIC</th><th>ClinVar</th><th>gnomAD AF</th>"
        "</tr></thead>"
        "<tbody>"
        + "\n".join(rows_html)
        + "</tbody></table></div>"
    )


# ── Assemble report ───────────────────────────────────────────────────────────

fastp      = parse_fastp(sm.input.fastp_json)
flagstat   = parse_flagstat(sm.input.flagstat)
mos_sum    = parse_mosdepth_summary(sm.input.mosdepth_summary)
mos_reg    = parse_mosdepth_regions(sm.input.mosdepth_regions)
contam     = parse_contamination(sm.input.contamination)
markdup    = parse_markdup(sm.input.markdup_metrics)
variants   = load_variants(sm.input.variants, MIN_VAF)

chart_vaf      = chart_vaf_histogram(variants)
chart_impact   = chart_impact_pie(variants)
chart_coverage = chart_coverage_bar(mos_reg)

logo_tag = ""
logo_data = logo_b64(LOGO_PATH)
if logo_data:
    logo_tag = f'<img src="{logo_data}" alt="logo"/>'

dup_cls   = "warn-text" if fastp["duplication_rate"] > 30 else "ok-text"
map_cls   = "warn-text" if flagstat["pct_mapped"] < 90 else "ok-text"
depth_cls = "warn-text" if mos_sum["mean_depth"] < 100 else "ok-text"
cont_cls  = "warn-text" if contam > 3.0 else "ok-text"

html = HTML_TEMPLATE.format(
    company=COMPANY,
    sample_id=SAMPLE_ID,
    panel_name=PANEL_NAME,
    panel_name_short=PANEL_NAME[:40] + ("…" if len(PANEL_NAME) > 40 else ""),
    run_date=RUN_DATE,
    logo_tag=logo_tag,
    variant_count=len(variants),
    total_reads_raw=fastp["total_reads_raw"],
    total_reads_trimmed=fastp["total_reads_trimmed"],
    q30_rate_trimmed=fastp["q30_rate_trimmed"],
    gc_content_raw=fastp["gc_content_raw"],
    duplication_rate=fastp["duplication_rate"],
    dup_cls=dup_cls,
    flagstat_total=flagstat["total_reads"],
    mapped_reads=flagstat["mapped_reads"],
    pct_mapped=flagstat["pct_mapped"],
    properly_paired=flagstat["properly_paired"],
    pct_properly_paired=flagstat["pct_properly_paired"],
    estimated_library_size=markdup["estimated_library_size"],
    mean_depth=mos_sum["mean_depth"],
    depth_cls=depth_cls,
    contamination=contam,
    cont_cls=cont_cls,
    pct_duplication=markdup["pct_duplication"],
    min_vaf_pct=round(MIN_VAF * 100, 1),
    chart_coverage=chart_coverage,
    chart_vaf=chart_vaf,
    chart_impact=chart_impact,
    variant_table_html=build_variant_table(variants),
)

Path(sm.output.report).parent.mkdir(parents=True, exist_ok=True)
Path(sm.output.report).write_text(html, encoding="utf-8")
