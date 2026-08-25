#!/usr/bin/env bash
# Download and prepare reference resources for GenRichi pipeline (hg38)
# Run once before first pipeline execution.
# Requirements: wget, samtools, gatk, bwa-mem2, bcftools, tabix

set -euo pipefail

RESOURCES="resources/reference"
ANNOTATION="resources/annotation"
VEP_CACHE="resources/vep_cache"

mkdir -p "$RESOURCES" "$ANNOTATION" "$VEP_CACHE"

# ── 1. Reference genome (hg38) ──────────────────────────────────────────────
echo "[1/7] Downloading hg38 reference..."
if [ ! -f "$RESOURCES/hg38.fa.gz" ]; then
  wget -q -c \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ics/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz" \
    -O "$RESOURCES/hg38.fa.gz"
fi
if [ ! -f "$RESOURCES/hg38.fa" ]; then
  gunzip -k "$RESOURCES/hg38.fa.gz"
fi

# ── 2. Index reference ──────────────────────────────────────────────────────
echo "[2/7] Indexing reference (bwa-mem2 + samtools + gatk)..."
[ ! -f "$RESOURCES/hg38.fa.bwt.2bit.64" ] && bwa-mem2 index "$RESOURCES/hg38.fa"
[ ! -f "$RESOURCES/hg38.fa.fai" ]          && samtools faidx "$RESOURCES/hg38.fa"
[ ! -f "$RESOURCES/hg38.dict" ] && gatk CreateSequenceDictionary \
  -R "$RESOURCES/hg38.fa" -O "$RESOURCES/hg38.dict"

# ── 3. dbSNP ────────────────────────────────────────────────────────────────
echo "[3/7] Downloading dbSNP..."
if [ ! -f "$RESOURCES/dbsnp_146.hg38.vcf.gz" ]; then
  wget -q -c \
    "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz" \
    -O "$RESOURCES/dbsnp_146.hg38.vcf.gz"
  wget -q -c \
    "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi" \
    -O "$RESOURCES/dbsnp_146.hg38.vcf.gz.tbi"
fi

# ── 4. gnomAD AF-only VCF (for Mutect2) ────────────────────────────────────
echo "[4/7] Downloading gnomAD AF-only..."
if [ ! -f "$RESOURCES/af-only-gnomad.hg38.vcf.gz" ]; then
  wget -q -c \
    "https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz" \
    -O "$RESOURCES/af-only-gnomad.hg38.vcf.gz"
  wget -q -c \
    "https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi" \
    -O "$RESOURCES/af-only-gnomad.hg38.vcf.gz.tbi"
fi

# ── 5. ClinVar ──────────────────────────────────────────────────────────────
echo "[5/7] Downloading ClinVar..."
if [ ! -f "$ANNOTATION/clinvar_20240101.vcf.gz" ]; then
  wget -q -c \
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz" \
    -O "$ANNOTATION/clinvar_20240101.vcf.gz"
  tabix -p vcf "$ANNOTATION/clinvar_20240101.vcf.gz"
fi

# ── 6. VEP cache (GRCh38) ──────────────────────────────────────────────────
echo "[6/7] Downloading VEP cache (this is ~14GB, will take time)..."
if [ ! -d "$VEP_CACHE/homo_sapiens" ]; then
  vep_install \
    --AUTO cf \
    --SPECIES homo_sapiens \
    --ASSEMBLY GRCh38 \
    --CACHEMULTI \
    --NO_HTSLIB \
    --CACHEDIR "$VEP_CACHE"
fi

# ── 7. COSMIC (requires account at cancer.sanger.ac.uk) ────────────────────
echo "[7/7] COSMIC: Download manually from https://cancer.sanger.ac.uk/cosmic/download"
echo "      Place file at: $ANNOTATION/CosmicCodingMuts_v99_GRCh38.vcf.gz"
echo "      Then run: tabix -p vcf $ANNOTATION/CosmicCodingMuts_v99_GRCh38.vcf.gz"
echo ""
echo "==> Resource setup complete!"
echo "    Edit config/config.yaml to confirm all paths, then run:"
echo "    ./run_pipeline.sh --execute"
