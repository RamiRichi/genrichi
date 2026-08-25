#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# GenRichi — Deploy hereditary pipeline to /home/rami/genrichi/
#
# Run ONCE from WSL:
#   bash /mnt/c/GenRichi/deploy_hereditary.sh
#
# What it does:
#   1. Copies all hereditary workflow files into the genrichi home project
#   2. Creates workflow/envs → ../envs symlink so Snakemake reuses pre-built
#      conda environments (same YAML hashes as the somatic pipeline)
#   3. Prints the commands to dry-run and run the pipeline
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

SRC="/mnt/c/GenRichi"
DST="/home/rami/genrichi"
BACKUP_DIR="$DST/workflow/rules/_somatic_backup_$(date +%Y%m%d_%H%M%S)"

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

ok()   { echo -e "${GREEN}  OK${NC}  $1"; }
warn() { echo -e "${YELLOW}WARN${NC}  $1"; }
err()  { echo -e "${RED} ERR${NC}  $1"; exit 1; }

echo ""
echo "=== GenRichi Hereditary Deploy ==="
echo "  Source : $SRC"
echo "  Target : $DST"
echo ""

# ── Guards ────────────────────────────────────────────────────────────────────
[ -d "$DST" ]                              || err "$DST not found."
[ -f "$SRC/workflow/Snakefile_hereditary" ] || err "$SRC/workflow/Snakefile_hereditary missing. Run from a WSL shell."

# ── 1. Snakefile_hereditary ───────────────────────────────────────────────────
echo "[1] Snakefile"
cp "$SRC/workflow/Snakefile_hereditary" "$DST/workflow/Snakefile_hereditary"
ok "workflow/Snakefile_hereditary"

# ── 2. Rule files ─────────────────────────────────────────────────────────────
echo ""
echo "[2] Rule files"
mkdir -p "$DST/workflow/rules"

# Back up existing shared rules (align, qc) if they already exist from somatic
for shared in align.smk qc.smk; do
    target="$DST/workflow/rules/$shared"
    if [ -f "$target" ]; then
        mkdir -p "$BACKUP_DIR"
        cp "$target" "$BACKUP_DIR/$shared"
        warn "Backed up existing $shared → $BACKUP_DIR/$shared"
    fi
done

# Copy all rule files
for f in align.smk qc.smk germline_calling.smk annotation_germline.smk acmg_classify.smk report_hereditary.smk; do
    cp "$SRC/workflow/rules/$f" "$DST/workflow/rules/$f"
    ok "rules/$f"
done

# ── 3. Python scripts ─────────────────────────────────────────────────────────
echo ""
echo "[3] Python scripts"
mkdir -p "$DST/workflow/scripts"
for f in vcf_to_germline_table.py acmg_classifier.py generate_hereditary_report.py; do
    cp "$SRC/workflow/scripts/$f" "$DST/workflow/scripts/$f"
    ok "scripts/$f"
done

# ── 4. Conda env symlink ──────────────────────────────────────────────────────
echo ""
echo "[4] Conda env symlink"
WFENVS="$DST/workflow/envs"

if [ -L "$WFENVS" ]; then
    ok "workflow/envs symlink already exists → $(readlink "$WFENVS")"
elif [ -d "$WFENVS" ]; then
    warn "workflow/envs is a real directory — checking contents"
else
    # Create symlink: workflow/envs -> ../envs (i.e. /home/rami/genrichi/envs)
    ln -s "$DST/envs" "$WFENVS"
    ok "Created: workflow/envs -> $DST/envs"
fi

echo "   Checking required env files:"
MISSING_ENVS=0
for needed in calling.yaml annotation.yaml align.yaml qc.yaml report.yaml; do
    if [ -f "$WFENVS/$needed" ]; then
        ok "  $needed"
    else
        warn "  $needed NOT FOUND — conda will try to create a new environment"
        MISSING_ENVS=$((MISSING_ENVS + 1))
    fi
done

if [ $MISSING_ENVS -gt 0 ]; then
    echo ""
    warn "Some env YAML files are missing from $WFENVS"
    warn "Snakemake will attempt to create those environments from scratch."
    warn "If the build fails, copy the YAML files from $DST/envs/ manually:"
    warn "  ls $DST/envs/"
fi

# ── 5. Config ─────────────────────────────────────────────────────────────────
echo ""
echo "[5] Config files"
mkdir -p "$DST/config"
cp "$SRC/config/hereditary_config.yaml" "$DST/config/hereditary_config.yaml"
cp "$SRC/config/hereditary_samples.tsv" "$DST/config/hereditary_samples.tsv"
ok "config/hereditary_config.yaml"
ok "config/hereditary_samples.tsv"

# ── 6. Panel BED ──────────────────────────────────────────────────────────────
echo ""
echo "[6] Panel BED"
mkdir -p "$DST/resources/panel"
cp "$SRC/resources/panel/hereditary_genes.bed" "$DST/resources/panel/hereditary_genes.bed"
ok "resources/panel/hereditary_genes.bed"

# ── 7. Output directories ──────────────────────────────────────────────────────
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
echo "Step 1 — activate snakemake environment:"
echo ""
echo "  conda activate snakemake"
echo "  cd $DST"
echo ""
echo "Step 2 — dry-run (verify job graph, no execution):"
echo ""
echo "  snakemake \\"
echo "    --snakefile workflow/Snakefile_hereditary \\"
echo "    --configfile config/hereditary_config.yaml \\"
echo "    --use-conda \\"
echo "    --conda-prefix $DST/.snakemake/conda \\"
echo "    --cores 8 --dry-run --reason"
echo ""
echo "Step 3 — execute (when dry-run looks correct):"
echo ""
echo "  snakemake \\"
echo "    --snakefile workflow/Snakefile_hereditary \\"
echo "    --configfile config/hereditary_config.yaml \\"
echo "    --use-conda \\"
echo "    --conda-prefix $DST/.snakemake/conda \\"
echo "    --rerun-incomplete --printshellcmds \\"
echo "    --cores 8"
echo ""
