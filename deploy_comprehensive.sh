#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# GenRichi — Deploy Somatic Comprehensive Panel to /home/rami/genrichi/
#
# Run ONCE from WSL:
#   bash /mnt/c/GenRichi/deploy_comprehensive.sh
#
# What it does:
#   1. Copies all comprehensive workflow files (Snakefile, rules, scripts)
#   2. Creates workflow/envs → ../envs symlink so Snakemake reuses pre-built
#      conda environments — no new environment builds needed for this pipeline
#   3. Copies config, panel BED, and sample sheet
#   4. Prints the dry-run / execute commands
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

SRC="/mnt/c/GenRichi"
DST="/home/rami/genrichi"

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

ok()   { echo -e "${GREEN}  OK${NC}  $1"; }
warn() { echo -e "${YELLOW}WARN${NC}  $1"; }
err()  { echo -e "${RED} ERR${NC}  $1"; exit 1; }

echo ""
echo "=== GenRichi Comprehensive Panel Deploy ==="
echo "  Source : $SRC"
echo "  Target : $DST"
echo ""

# ── Guards ────────────────────────────────────────────────────────────────────
[ -d "$DST" ]                                  || err "$DST not found."
[ -f "$SRC/workflow/Snakefile_comprehensive" ]  || err "$SRC/workflow/Snakefile_comprehensive missing. Run from WSL."

# ── 1. Snakefile ──────────────────────────────────────────────────────────────
echo "[1] Snakefile"
cp "$SRC/workflow/Snakefile_comprehensive" "$DST/workflow/Snakefile_comprehensive"
ok "workflow/Snakefile_comprehensive"

# ── 2. Rule files ─────────────────────────────────────────────────────────────
echo ""
echo "[2] Rule files"
mkdir -p "$DST/workflow/rules"

RULES=(
    qc_paired.smk
    align_paired.smk
    somatic_calling_paired.smk
    cnv_calling.smk
    msi_scoring.smk
    annotation_comprehensive.smk
    report_comprehensive.smk
)

for f in "${RULES[@]}"; do
    cp "$SRC/workflow/rules/$f" "$DST/workflow/rules/$f"
    ok "rules/$f"
done

# ── 3. Python scripts ─────────────────────────────────────────────────────────
echo ""
echo "[3] Python scripts"
mkdir -p "$DST/workflow/scripts"

SCRIPTS=(
    calculate_cnv.py
    calculate_msi.py
    vcf_to_somatic_table.py
    generate_comprehensive_report.py
)

for f in "${SCRIPTS[@]}"; do
    cp "$SRC/workflow/scripts/$f" "$DST/workflow/scripts/$f"
    ok "scripts/$f"
done

# ── 4. Conda env symlink (reuse existing pre-built envs) ─────────────────────
echo ""
echo "[4] Conda env symlink"
WFENVS="$DST/workflow/envs"

if [ -L "$WFENVS" ]; then
    ok "workflow/envs symlink already exists → $(readlink "$WFENVS")"
elif [ -d "$WFENVS" ]; then
    warn "workflow/envs is a real directory (from somatic pipeline) — reusing"
else
    ln -s "$DST/envs" "$WFENVS"
    ok "Created: workflow/envs -> $DST/envs"
fi

echo "   Checking required env files:"
MISSING=0
for needed in align.yaml calling.yaml annotation.yaml report.yaml; do
    if [ -f "$DST/envs/$needed" ] || [ -f "$WFENVS/$needed" ]; then
        ok "  $needed"
    else
        warn "  $needed NOT FOUND"
        MISSING=$((MISSING + 1))
    fi
done

if [ $MISSING -gt 0 ]; then
    echo ""
    warn "Missing env YAML files — Snakemake will try to create them."
    warn "If builds fail, copy YAMLs manually:"
    warn "  ls $DST/envs/"
fi

# ── 5. Config files ───────────────────────────────────────────────────────────
echo ""
echo "[5] Config files"
mkdir -p "$DST/config"
cp "$SRC/config/comprehensive_config.yaml"  "$DST/config/comprehensive_config.yaml"
cp "$SRC/config/comprehensive_samples.tsv"  "$DST/config/comprehensive_samples.tsv"
ok "config/comprehensive_config.yaml"
ok "config/comprehensive_samples.tsv"

# ── 6. Panel BED ──────────────────────────────────────────────────────────────
echo ""
echo "[6] Panel BED"
mkdir -p "$DST/resources/panel"
cp "$SRC/resources/panel/comprehensive_genes.bed" \
   "$DST/resources/panel/comprehensive_genes.bed"
ok "resources/panel/comprehensive_genes.bed"

# ── 7. Output directories ─────────────────────────────────────────────────────
echo ""
echo "[7] Output directories"
mkdir -p "$DST/results"
mkdir -p "$DST/logs"
ok "results/ and logs/ ready"

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo -e "${GREEN}Deploy complete!${NC}"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Step 1 — update sample sheet with real paths:"
echo ""
echo "  nano $DST/config/comprehensive_samples.tsv"
echo ""
echo "  Columns: sample_id  tumor_r1  tumor_r2  normal_r1  normal_r2"
echo "           patient_id  sex  tumor_type"
echo ""
echo "Step 2 — activate environment and dry-run:"
echo ""
echo "  conda activate snakemake"
echo "  cd $DST"
echo ""
echo "  snakemake \\"
echo "    --snakefile workflow/Snakefile_comprehensive \\"
echo "    --configfile config/comprehensive_config.yaml \\"
echo "    --use-conda \\"
echo "    --conda-prefix $DST/.snakemake/conda \\"
echo "    --cores 8 --dry-run"
echo ""
echo "Step 3 — execute:"
echo ""
echo "  snakemake \\"
echo "    --snakefile workflow/Snakefile_comprehensive \\"
echo "    --configfile config/comprehensive_config.yaml \\"
echo "    --use-conda \\"
echo "    --conda-prefix $DST/.snakemake/conda \\"
echo "    --rerun-incomplete --printshellcmds \\"
echo "    --cores 8"
echo ""
echo "─── Key outputs per sample ─────────────────────────────────────"
echo ""
echo "  results/{sample}/report/{sample}_comprehensive_report.html"
echo "  results/{sample}/snv/{sample}.pass.vcf.gz       (somatic variants)"
echo "  results/{sample}/cnv/{sample}.call.cns           (CNV segments)"
echo "  results/{sample}/cnv/{sample}-scatter.png        (genome-wide CNV plot)"
echo "  results/{sample}/msi/{sample}.msi                (MSI score)"
echo ""
