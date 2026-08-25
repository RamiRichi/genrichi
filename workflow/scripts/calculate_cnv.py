"""
Panel-based CNV estimation from mosdepth per-region coverage.

Algorithm:
  1. Load per-region mean depth for tumor and normal (mosdepth .regions.bed.gz)
  2. Compute log2(tumor / normal) per region — the copy-number ratio
  3. Apply median smoothing across genomic windows
  4. Call segments: log2 >= amp_threshold → Amplification
                    log2 <= del_threshold → Deletion
  5. Write CNR (per-region), CNS (segments), call.cns, and scatter plot

No cnvkit or scipy required — uses only pandas + numpy + matplotlib,
all available in the existing report.yaml conda environment.
"""

import gzip
import io
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

sm = snakemake  # type: ignore[name-defined]

TUMOR_REGIONS  = sm.input.tumor_regions
NORMAL_REGIONS = sm.input.normal_regions
OUT_CNR        = sm.output.cnr
OUT_CNS        = sm.output.cns
OUT_CALL_CNS   = sm.output.call_cns
OUT_SCATTER    = sm.output.scatter

AMP_THR  = float(sm.params.amp_threshold)
DEL_THR  = float(sm.params.del_threshold)
MIN_PROB = int(sm.params.min_probes)   # minimum regions per segment
PSEUDO   = 1.0                          # pseudocount to avoid log(0)

CHROM_ORDER = [f"chr{i}" for i in list(range(1, 23)) + ["X", "Y"]]


# ── Load mosdepth regions files ───────────────────────────────────────────────
def _load_regions(path: str) -> pd.DataFrame:
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        df = pd.read_csv(
            fh, sep="\t",
            header=None,
            names=["chrom", "start", "end", "region", "mean_depth"],
            dtype={"chrom": str, "start": int, "end": int, "mean_depth": float},
            comment="#",
        )
    return df


tumor  = _load_regions(TUMOR_REGIONS)
normal = _load_regions(NORMAL_REGIONS)

# Merge on genomic coordinates
merged = pd.merge(
    tumor,
    normal[["chrom", "start", "end", "mean_depth"]].rename(
        columns={"mean_depth": "normal_depth"}
    ),
    on=["chrom", "start", "end"],
    how="inner",
).rename(columns={"mean_depth": "tumor_depth"})

# ── Log2 ratio ─────────────────────────────────────────────────────────────────
merged["log2"] = np.log2(
    (merged["tumor_depth"] + PSEUDO) / (merged["normal_depth"] + PSEUDO)
)
merged["gc"] = 0.5   # placeholder GC — no reference needed for ratio
merged["gene"] = merged.get("region", merged["chrom"])

# Chromosome order
merged["chrom_order"] = pd.Categorical(
    merged["chrom"], categories=CHROM_ORDER, ordered=True
)
merged = merged.sort_values(["chrom_order", "start"]).drop(columns="chrom_order")
merged = merged.reset_index(drop=True)

# ── Median smoothing (rolling 3-window) ────────────────────────────────────────
merged["log2_smooth"] = (
    merged.groupby("chrom")["log2"]
    .transform(lambda x: x.rolling(3, center=True, min_periods=1).median())
)

# ── Write CNR file (cnvkit-compatible) ────────────────────────────────────────
os.makedirs(os.path.dirname(OUT_CNR), exist_ok=True)
cnr_cols = ["chromosome", "start", "end", "gene", "log2", "depth", "weight"]
cnr_df = pd.DataFrame({
    "chromosome": merged["chrom"],
    "start":      merged["start"],
    "end":        merged["end"],
    "gene":       merged["gene"],
    "log2":       merged["log2_smooth"].round(4),
    "depth":      merged["tumor_depth"].round(2),
    "weight":     1.0,
})
cnr_df.to_csv(OUT_CNR, sep="\t", index=False)


# ── Simple segmentation: group consecutive same-direction regions ──────────────
def _segment(df: pd.DataFrame) -> pd.DataFrame:
    """
    Merge adjacent regions with log2_smooth in the same direction
    (both above amp_thr, both below del_thr, or both neutral).
    """
    segments = []
    for chrom, grp in df.groupby("chrom", sort=False):
        grp = grp.reset_index(drop=True)
        direction = np.where(
            grp["log2_smooth"] >= AMP_THR, 1,
            np.where(grp["log2_smooth"] <= DEL_THR, -1, 0)
        )
        i = 0
        while i < len(grp):
            j = i + 1
            while j < len(grp) and direction[j] == direction[i]:
                j += 1
            seg = grp.iloc[i:j]
            segments.append({
                "chromosome": chrom,
                "start":      int(seg["start"].iloc[0]),
                "end":        int(seg["end"].iloc[-1]),
                "gene":       seg["gene"].iloc[0],
                "log2":       round(float(seg["log2_smooth"].mean()), 4),
                "probes":     len(seg),
                "weight":     1.0,
            })
            i = j
    return pd.DataFrame(segments)


cns_df = _segment(merged)
cns_df.to_csv(OUT_CNS, sep="\t", index=False)

# ── Call segments (integer copy number) ───────────────────────────────────────
def _call_cn(log2: float) -> int:
    """Rough integer copy number from log2 ratio (assume diploid normal)."""
    ratio = 2 ** log2
    cn    = max(0, round(2 * ratio))
    return cn


call_df = cns_df.copy()
call_df["cn"]   = call_df["log2"].apply(_call_cn)
call_df["type"] = np.where(
    call_df["log2"] >= AMP_THR, "Amplification",
    np.where(call_df["log2"] <= DEL_THR, "Deletion", "Neutral")
)
call_df.to_csv(OUT_CALL_CNS, sep="\t", index=False)

# ── Scatter plot ───────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(14, 4))

# Assign x position per region
x_offset = 0
x_positions = []
chrom_ticks = []
chrom_labels = []
chrom_colors = {}
palette = ["#3498db", "#2ecc71"]

for ci, (chrom, grp) in enumerate(merged.groupby("chrom_order" if "chrom_order" in merged else "chrom", sort=False)):
    chrom = str(chrom)
    color = palette[ci % 2]
    chrom_colors[chrom] = color
    positions = list(range(x_offset, x_offset + len(grp)))
    x_positions.extend(positions)
    mid = x_offset + len(grp) // 2
    chrom_ticks.append(mid)
    chrom_labels.append(chrom.replace("chr", ""))
    x_offset += len(grp) + 5

merged["x_pos"] = x_positions
for chrom, grp in merged.groupby("chrom"):
    ax.scatter(
        grp["x_pos"], grp["log2_smooth"],
        c=chrom_colors.get(str(chrom), "#3498db"),
        s=8, alpha=0.6, linewidths=0,
    )

ax.axhline(0,       color="#555", lw=0.8, linestyle="--")
ax.axhline(AMP_THR, color="#e74c3c", lw=1, linestyle=":", label=f"Amp ({AMP_THR:+.2f})")
ax.axhline(DEL_THR, color="#3498db", lw=1, linestyle=":", label=f"Del ({DEL_THR:+.2f})")

ax.set_xticks(chrom_ticks)
ax.set_xticklabels(chrom_labels, fontsize=7, rotation=45)
ax.set_ylabel("log2 Ratio (T/N)", fontsize=10)
ax.set_title("Copy Number Profile — Panel Regions", fontsize=12, fontweight="bold")
ax.set_ylim(-4, 4)
ax.legend(fontsize=8, loc="upper right")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.tight_layout()
fig.savefig(OUT_SCATTER, dpi=120, bbox_inches="tight")
plt.close(fig)
