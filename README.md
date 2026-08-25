# GenRichi NGS Diagnostic Pipeline

> **Clinical-grade Snakemake pipeline for somatic variant calling, hereditary panel analysis, and automated diagnostic reporting.**

Built by a clinical bioinformatician with 7+ years in ISO 15189-accredited molecular diagnostics. Designed for real-world NGS laboratory workflows — not just research.

---

## What GenRichi Does

GenRichi takes raw sequencing reads (FASTQ) and delivers a ready-to-review clinical HTML report — fully automated.

```
FASTQ → QC (fastp) → Alignment (BWA-MEM2) → Variant Calling (GATK Mutect2)
     → Annotation (VEP + COSMIC + ClinVar) → Clinical HTML Report
```

### Three Analysis Modes

| Mode | Use Case | Key Output |
|------|----------|-----------|
| **Somatic** | Cancer hotspot panel | Somatic variants, VAF, COSMIC annotation |
| **Hereditary** | BRCA, Lynch, cancer predisposition | ACMG-classified germline variants |
| **Comprehensive** | Combined somatic + TMB + MSI + CNV | Full oncology dashboard |

---

## Key Features

- **Somatic variant calling** — GATK Mutect2 with Panel-of-Normals support
- **Germline analysis** — ACMG variant classification built-in
- **VEP annotation** — with COSMIC v99, ClinVar, gnomAD allele frequencies
- **TMB & MSI scoring** — tumor mutational burden and microsatellite instability
- **CNV detection** — copy number variant calling
- **Automated HTML reports** — patient-ready, clinician-friendly output
- **Clinical portal** — web interface for order management and report delivery
- **Reproducible** — Conda environments pinned per rule
- **Configurable** — swap panels, references, thresholds via `config.yaml`

---

## Quick Start

```bash
# 1. Clone
git clone https://github.com/RamiRichi/genrichi.git
cd genrichi

# 2. Set up your samples
cp config/samples.tsv my_samples.tsv
# Edit my_samples.tsv with your sample IDs and FASTQ paths

# 3. Configure references
nano config/config.yaml
# Set paths to hg38 reference, dbSNP, gnomAD, VEP cache

# 4. Run (somatic panel)
snakemake --cores 8 --use-conda

# 5. View report
open results/{sample}/report/{sample}_report.html
```

---

## Repository Structure

```
genrichi/
├── workflow/
│   ├── Snakefile                    # Main pipeline (somatic)
│   ├── Snakefile_hereditary         # Hereditary panel
│   ├── Snakefile_comprehensive      # Full oncology mode
│   ├── rules/                       # Modular Snakemake rules
│   │   ├── qc.smk                   # fastp QC
│   │   ├── align.smk                # BWA-MEM2 alignment
│   │   ├── calling.smk              # GATK Mutect2
│   │   ├── annotation.smk           # VEP + COSMIC + ClinVar
│   │   ├── acmg_classify.smk        # ACMG germline classification
│   │   ├── cnv_calling.smk          # CNV detection
│   │   ├── msi_scoring.smk          # MSI calculation
│   │   └── report.smk               # HTML report generation
│   ├── scripts/                     # Python analysis scripts
│   └── envs/                        # Conda environment files
├── config/
│   ├── config.yaml                  # Main configuration
│   └── samples.tsv                  # Sample sheet template
├── resources/
│   └── panel/                       # BED files (hotspots, genes)
├── portal/                          # Clinical web portal (Flask)
└── test_data/                       # Test FASTQ files
```

---

## Configuration

Key parameters in `config/config.yaml`:

```yaml
calling:
  mutect2:
    extra: "--min-base-quality-score 20"
  filter:
    min_af: 0.02        # minimum allele frequency
    min_depth: 30       # minimum read depth
    min_alt_reads: 5    # minimum alt-supporting reads

annotation:
  vep:
    genome_build: GRCh38
    extra: "--everything --canonical --check_existing"
```

---

## Clinical Portal

GenRichi includes a web-based clinical portal for:
- Submitting sequencing orders
- Tracking pipeline status
- Viewing and downloading reports
- Multi-user access with role management

```bash
cd portal
pip install -r requirements.txt
python app.py
```

---

## Requirements

- Python ≥ 3.10
- Snakemake ≥ 7.0
- Conda / Mamba (for environment management)
- BWA-MEM2, GATK 4.x, VEP (auto-installed via Conda envs)

**Reference data required** (not included due to size):
- hg38 reference genome (BWA-indexed + GATK dict)
- dbSNP 146 (hg38)
- gnomAD AF-only VCF (hg38)
- VEP cache (GRCh38)
- COSMIC v99 (requires COSMIC license)

---

## Background

GenRichi was developed alongside clinical practice in an ISO 15189-accredited molecular diagnostics laboratory. The pipeline reflects real diagnostic requirements: reproducibility, auditability, and report quality suitable for clinical decision-making.

The name **GenRichi** reflects the mission: making genomic diagnostics richer — more complete, more accessible, and more actionable.

---

## Author

**Rami Richi** — Bioinformatician & Clinical Laboratory Specialist  
M.Sc. Bioinformatics & Molecular Biology (Karolinska Institute / Hochschule Mittweida)  
[LinkedIn](https://linkedin.com/in/rami-richi) · [GitHub](https://github.com/RamiRichi)

---

## License

MIT License — see [LICENSE](LICENSE) for details.

> For clinical deployment, collaboration, or licensing inquiries:  
> **rami_richi83@yahoo.com**
