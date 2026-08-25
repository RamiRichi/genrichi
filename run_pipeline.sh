#!/usr/bin/env bash
# GenRichi Somatic Hotspot Panel — pipeline launcher
# Usage:
#   ./run_pipeline.sh                   # dry-run (preview only)
#   ./run_pipeline.sh --execute         # full run (local)
#   ./run_pipeline.sh --execute --cores 16
#   ./run_pipeline.sh --execute --slurm  # HPC submission via SLURM profile

set -euo pipefail

CORES="${CORES:-8}"
SNAKEFILE="workflow/Snakefile"
CONFIG="config/config.yaml"

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
  echo "==> DRY RUN — no jobs will be executed"
  "${BASE_CMD[@]}" --dry-run --reason
  echo ""
  echo "Run with --execute to start the pipeline."
elif $SLURM; then
  echo "==> Submitting to SLURM cluster"
  "${BASE_CMD[@]}" \
    --cluster "sbatch --mem={resources.mem_mb}M --cpus-per-task={threads} --time=4:00:00" \
    --jobs 50 \
    --latency-wait 60
else
  echo "==> Running locally with $CORES cores"
  "${BASE_CMD[@]}"
fi
