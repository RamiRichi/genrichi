#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# GenRichi — Deploy Multi-Panel Dashboard to /home/rami/genrichi/
#
# Run from WSL:
#   bash /mnt/c/GenRichi/deploy_dashboard.sh
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

SRC="/mnt/c/GenRichi"
DST="/home/rami/genrichi"

GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; NC='\033[0m'
ok()   { echo -e "${GREEN}  OK${NC}  $1"; }
warn() { echo -e "${YELLOW}WARN${NC}  $1"; }
err()  { echo -e "${RED} ERR${NC}  $1"; exit 1; }

echo ""
echo "=== GenRichi Dashboard Deploy ==="
echo "  Source : $SRC"
echo "  Target : $DST"
echo ""

[ -d "$DST" ] || err "$DST not found."
[ -f "$SRC/workflow/Snakefile_dashboard" ] || err "Snakefile_dashboard missing."

echo "[1] Snakefile"
cp "$SRC/workflow/Snakefile_dashboard" "$DST/workflow/Snakefile_dashboard"
ok "workflow/Snakefile_dashboard"

echo ""
echo "[2] Script"
cp "$SRC/workflow/scripts/generate_dashboard.py" "$DST/workflow/scripts/generate_dashboard.py"
ok "scripts/generate_dashboard.py"

echo ""
echo "[3] Config"
mkdir -p "$DST/config"
cp "$SRC/config/dashboard_config.yaml"   "$DST/config/dashboard_config.yaml"
cp "$SRC/config/dashboard_patients.tsv"  "$DST/config/dashboard_patients.tsv"
ok "config/dashboard_config.yaml"
ok "config/dashboard_patients.tsv"

echo ""
echo "[4] Output dir"
mkdir -p "$DST/results/dashboard"
ok "results/dashboard/"

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo -e "${GREEN}Deploy complete!${NC}"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Run the dashboard:"
echo ""
echo "  conda activate snakemake"
echo "  cd $DST"
echo ""
echo "  snakemake \\"
echo "    --snakefile workflow/Snakefile_dashboard \\"
echo "    --configfile config/dashboard_config.yaml \\"
echo "    --use-conda \\"
echo "    --conda-prefix $DST/.snakemake/conda \\"
echo "    --cores 4"
echo ""
echo "Output:"
echo "  results/dashboard/NA12878_dashboard.html"
echo ""
