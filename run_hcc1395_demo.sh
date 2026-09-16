#!/bin/bash
# GenRichi HCC1395 Demo Analysis
# Downloads chr17-only exome data (Griffith Lab, Washington University)
# HCC1395 breast cancer cell line; this chr17-only subset validates the TP53 p.Arg175His
# somatic hotspot mutation (chr17:7675088 C>T, VAF ~0.98) as its PASS variant.
# Run from: /home/rami/genrichi/

set -euo pipefail

cd /home/rami/genrichi

DEMO_DIR="/home/rami/genrichi/hcc1395_demo"

# ── 1. Download chr17 exome FASTQ (~680 MB total) ────────────────────────────
echo "=== Downloading HCC1395 chr17 exome FASTQ data ==="
mkdir -p "$DEMO_DIR"
cd "$DEMO_DIR"

[ -f Exome_Tumor.tar.gz ] || wget -c http://genomedata.org/pmbio-workshop/fastqs/chr17/Exome_Tumor.tar.gz
[ -f Exome_Norm.tar.gz  ] || wget -c http://genomedata.org/pmbio-workshop/fastqs/chr17/Exome_Norm.tar.gz

echo "=== Extracting archives ==="
tar -xzf Exome_Tumor.tar.gz
tar -xzf Exome_Norm.tar.gz

# ── 2. Locate FASTQ files ────────────────────────────────────────────────────
echo "=== Locating FASTQ files ==="
TUMOR_R1=$(find "$DEMO_DIR" -name "*.fastq.gz" | grep -i "tumor" | grep -iE "_R1|_1\." | sort | head -1)
TUMOR_R2=$(find "$DEMO_DIR" -name "*.fastq.gz" | grep -i "tumor" | grep -iE "_R2|_2\." | sort | head -1)
NORM_R1=$(find  "$DEMO_DIR" -name "*.fastq.gz" | grep -i "norm"  | grep -iE "_R1|_1\." | sort | head -1)
NORM_R2=$(find  "$DEMO_DIR" -name "*.fastq.gz" | grep -i "norm"  | grep -iE "_R2|_2\." | sort | head -1)

# Fallback: list all and let user check
if [[ -z "$TUMOR_R1" ]]; then
  echo "Could not auto-detect files. Files found in $DEMO_DIR:"
  find "$DEMO_DIR" -name "*.fastq.gz"
  echo ""
  echo "Please set TUMOR_R1, TUMOR_R2, NORM_R1, NORM_R2 and re-run."
  exit 1
fi

echo "  Tumor R1:  $TUMOR_R1"
echo "  Tumor R2:  $TUMOR_R2"
echo "  Normal R1: $NORM_R1"
echo "  Normal R2: $NORM_R2"

# ── 3. Write samples TSV ─────────────────────────────────────────────────────
cd /home/rami/genrichi

cat > config/hcc1395_samples.tsv << TSV
sample_id	tumor_r1	tumor_r2	normal_r1	normal_r2	patient_id	sex	tumor_type
HCC1395_demo	${TUMOR_R1}	${TUMOR_R2}	${NORM_R1}	${NORM_R2}	HCC1395	Female	Breast_Cancer
TSV

echo "=== Written: config/hcc1395_samples.tsv ==="

# ── 4. Run pipeline ──────────────────────────────────────────────────────────
echo ""
echo "=== Running somatic comprehensive pipeline ==="

snakemake \
  --snakefile workflow/Snakefile_comprehensive \
  --config samples=config/hcc1395_samples.tsv \
  --software-deployment-method conda \
  --conda-frontend mamba \
  --rerun-incomplete \
  --printshellcmds \
  --cores 4 \
  2>&1 | tee /tmp/hcc1395_run.log

echo ""
echo "=== Done! Report: results/HCC1395_demo/report/HCC1395_demo_comprehensive_report.html ==="
