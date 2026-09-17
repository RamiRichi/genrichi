#!/usr/bin/env bash
# GenRichi Hereditary Cancer Panel — pipeline launcher
#
# IMPORTANT: Run this script from the project root (the folder containing
#            workflow/, config/, resources/ etc.)
#
#   cd /path/to/GenRichi          ← must be here
#   ./run_hereditary.sh           # dry-run (preview only)
#   ./run_hereditary.sh --execute         # full run (local)
#   ./run_hereditary.sh --execute --cores 16
#   ./run_hereditary.sh --execute --slurm # HPC via SLURM profile

set -euo pipefail

# ── Guard: make sure we're in the right directory ────────────────────────────
if [[ ! -f "workflow/Snakefile_hereditary" ]]; then
  echo "ERROR: workflow/Snakefile_hereditary not found."
  echo "       Please run this script from the GenRichi project root, e.g.:"
  echo "         cd /path/to/GenRichi && ./run_hereditary.sh"
  exit 1
fi

CORES="${CORES:-8}"
SNAKEFILE="workflow/Snakefile_hereditary"
CONFIG="config/hereditary_config.yaml"

DRY=true
SLURM=false

for arg in "$@"; do
  case $arg in
    --execute) DRY=false ;;
    --slurm)   SLURM=true ;;
    --cores)   shift; CORES="$1" ;;
  esac
done

BASE_CMD=(
  snakemake
  --snakefile "$SNAKEFILE"
  --configfile "$CONFIG"
  --use-conda
  --conda-frontend mamba
  --rerun-incomplete
  --printshellcmds
  --cores "$CORES"
)

if $DRY; then
  echo "==> DRY RUN (hereditary panel) — no jobs will be executed"
  "${BASE_CMD[@]}" --dry-run --reason
  echo ""
  echo "Run with --execute to start the pipeline."
elif $SLURM; then
  echo "==> Submitting hereditary pipeline to SLURM cluster"
  "${BASE_CMD[@]}" \
    --cluster "sbatch --mem={resources.mem_mb}M --cpus-per-task={threads} --time=4:00:00" \
    --jobs 50 \
    --latency-wait 60
else
  echo "==> Running hereditary pipeline locally with $CORES cores"
  "${BASE_CMD[@]}"
fi
