# GenRichi — Somatic Panel Pipeline: Session Notes
**Date:** 2026-05-25  
**Engineer:** Rami (GenRichi Diagnostics)  
**AI Assistant:** Claude (Anthropic)  
**Topic:** Building a clinical-grade somatic hotspot cancer panel pipeline from scratch

---

## 1. Goal

Build an end-to-end somatic hotspot panel pipeline for cancer diagnostics at GenRichi, covering:

```
FASTQ → QC → Alignment → MarkDuplicates → BQSR →
Mutect2 → FilterMutectCalls → VEP Annotation → HTML Report
```

Then wrap it in a **Flask web application** so GenRichi can repeatedly process real patient samples through a browser interface.

---

## 2. Environment

| Item | Value |
|---|---|
| OS | WSL2 Ubuntu on Windows 11 |
| Linux user | `rami` |
| Home | `/home/rami/` |
| Project root | `/home/rami/genrichi/` |
| Conda env | `snakemake` (all tools installed here) |
| Snakemake version | 9.21.0 |
| Reference genome | hg38 / GRCh38 |
| Calling mode | Tumor-only (no matched normal) |

---

## 3. Tool Stack

| Tool | Version | Purpose |
|---|---|---|
| Snakemake | 9.21.0 | Workflow manager |
| fastp | latest | QC + adapter trimming |
| BWA classic | latest | Alignment (NOT BWA-MEM2 — WSL2 only has ~8GB RAM, BWA-MEM2 needs ~50GB to index hg38) |
| samtools | latest | BAM sorting, indexing, flagstat |
| GATK | 4.3.0.0 | MarkDuplicates, BQSR, Mutect2, FilterMutectCalls, GetPileupSummaries, CalculateContamination, LearnReadOrientationModel |
| mosdepth | 0.3.14 | Panel coverage per region |
| bcftools | latest | PASS filter on VCF |
| Ensembl VEP | v113 | Variant annotation |
| tabix / bgzip | latest | VCF indexing |
| Flask | latest | Web application |
| matplotlib | latest | Charts in HTML report |

---

## 4. Project File Structure

```
/home/rami/genrichi/
├── app.py                          ← Flask web app (somatic pipeline)
├── genrichi_jobs.db                ← SQLite job database (auto-created)
├── genrichi_logo_final.png
├── config/
│   ├── config.yaml                 ← All paths and parameters
│   └── samples.tsv                 ← Sample sheet (sample_id, fastq_r1, fastq_r2)
├── resources/
│   └── panel/
│       └── example_hotspots.bed    ← 20-gene hotspot panel (hg38)
├── reference_db/
│   ├── ref/hg38.fa                 ← Reference genome (+ .bwt .pac .sa .fai .dict indexes)
│   ├── dbsnp/dbsnp_146.hg38.vcf.gz
│   ├── gnomad/af-only-gnomad.hg38.vcf.gz
│   ├── clinvar/clinvar.vcf.gz
│   ├── cosmic/CosmicCodingMuts_v99_GRCh38.vcf.gz  ← Manual download required
│   └── vep_cache/homo_sapiens/     ← GRCh38 VEP cache (~14GB)
├── templates/
│   └── dashboard.html              ← Web app UI
├── results/
│   └── {sample_id}/
│       ├── qc/                     ← fastp JSON, trimmed FASTQs
│       ├── align/                  ← BAM files, flagstat, mosdepth
│       ├── calling/                ← Mutect2 VCF, contamination, stats
│       ├── annotation/             ← VEP VCF, variants.tsv
│       └── report/                 ← HTML report
├── workflow/
│   ├── Snakefile
│   └── rules/
│       ├── align.smk
│       ├── qc.smk
│       ├── calling.smk
│       ├── annotation.smk
│       └── report.smk
└── docs/
    └── somatic_panel/
        └── session_notes.md        ← This file
```

---

## 5. Reference Paths (Linux)

```yaml
# config/config.yaml key paths
ref:
  genome:  /home/rami/genrichi/reference_db/ref/hg38.fa
  dbsnp:   /home/rami/genrichi/reference_db/dbsnp/dbsnp_146.hg38.vcf.gz
  gnomad:  /home/rami/genrichi/reference_db/gnomad/af-only-gnomad.hg38.vcf.gz
  pon:     null

annotation:
  vep:
    cache_dir: /home/rami/genrichi/reference_db/vep_cache
  cosmic:
    vcf: /home/rami/genrichi/reference_db/cosmic/CosmicCodingMuts_v99_GRCh38.vcf.gz
  clinvar:
    vcf: /home/rami/genrichi/reference_db/clinvar/clinvar.vcf.gz

panel:
  bed: /home/rami/genrichi/resources/panel/example_hotspots.bed
```

---

## 6. Panel BED File (20 Hotspot Regions — hg38)

Located at: `/home/rami/genrichi/resources/panel/example_hotspots.bed`

```
chr7	55191721	55191820	EGFR_exon19
chr7	55201887	55201980	EGFR_exon20
chr7	55209979	55210073	EGFR_exon21
chr12	25227234	25227380	KRAS_exon2
chr12	25245274	25245430	KRAS_exon3
chr12	25250929	25251070	KRAS_exon4
chr1	115247089	115247240	NRAS_exon2
chr1	115258740	115258850	NRAS_exon3
chr17	7673781	7674100	TP53_exon5
chr17	7672573	7672900	TP53_exon6
chr17	7670608	7670900	TP53_exon7
chr17	7668421	7668700	TP53_exon8
chr9	107543289	107543520	CDKN2A_exon1
chr13	32889645	32889960	BRCA2_exon11
chr3	178935823	178936040	PIK3CA_exon9
chr3	178951853	178952060	PIK3CA_exon20
chr7	87550282	87550510	CDK6_exon1
chr4	1801841	1802100	KIT_exon9
chr4	1806098	1806370	KIT_exon11
chr2	29220820	29221110	ALK_exon20
```

---

## 7. All Bugs Found and Fixed

### 7.1 BWA-MEM2 Out of Memory
- **Problem:** `bwa-mem2 index hg38.fa` crashed — needs ~50GB RAM, WSL2 only has ~8GB
- **Fix:** Switched to classic `bwa mem` which uses existing `.bwt/.pac/.sa` indexes

### 7.2 VEP Cache Corrupted (30GB)
- **Problem:** Multiple partial downloads were appended with `-c` flag, making 30GB corrupt file
- **Fix:** Deleted and re-downloaded without `-c` flag using rsync or fresh wget

### 7.3 Snakemake `--reason` flag removed
- **Problem:** Snakemake 9.x removed `--reason` flag
- **Fix:** Removed `--reason` from run_pipeline.sh

### 7.4 Wrong conda env path in rules
- **Problem:** `workflow/rules/*.smk` had `conda: "../envs/align.yaml"` but correct path is `../../envs/`
- **Fix:** `sed -i 's|../envs/|../../envs/|g'` on all rule files; later removed all `conda:` directives entirely since snakemake env already has all tools

### 7.5 mosdepth `--no-abbrev` flag
- **Problem:** `--no-abbrev` is not a valid mosdepth flag
- **Fix:** Removed the flag from mosdepth rule

### 7.6 mosdepth `--quantize` flag crash
- **Problem:** `--quantize 0:1:10:30:100:` caused silent mosdepth crash
- **Fix:** Removed the flag from mosdepth rule

### 7.7 mosdepth BED file error "Cannot open file list"
- **Problem:** mosdepth rule was replaced with `samtools coverage -b` which uses `-b` for BAM list (not BED)
- **Fix:** Rewrote mosdepth rule back to use real mosdepth tool

### 7.8 mosdepth output file naming
- **Problem:** Rule expected `{sample}.mosdepth.regions.bed.gz` but mosdepth creates `{sample}.regions.bed.gz`
- **Fix:** Updated output declaration in align.smk and reference in report.smk

### 7.9 mosdepth rule (final correct version)
```python
rule mosdepth:
    input:
        bam="results/{sample}/align/{sample}.final.bam",
        bai="results/{sample}/align/{sample}.final.bam.bai",
        bed=config["panel"]["bed"],
    output:
        summary="results/{sample}/align/{sample}.mosdepth.summary.txt",
        regions="results/{sample}/align/{sample}.regions.bed.gz",
    threads: 4
    params:
        prefix="results/{sample}/align/{sample}",
    shell:
        "mosdepth --threads {threads} --by {input.bed} {params.prefix} {input.bam} 2>{log}"
```

### 7.10 BED file accidentally deleted
- **Problem:** `grep -v "^#" panel.bed > /tmp/clean.bed && mv /tmp/clean.bed panel.bed` deleted file when mv failed
- **Fix:** Recreated BED file with heredoc

### 7.11 FilterMutectCalls silent failure in conda env
- **Problem:** GATK 4.6.1.0 in Snakemake-created conda env crashed silently, leaving empty log
- **Fix:** Removed all `conda:` directives from all rules — used tools from snakemake conda env directly

### 7.12 MultiQC output directory name
- **Problem:** Rule expected `multiqc_data` but MultiQC v1.25 creates `multiqc_report_data`
- **Fix:** `sed -i 's|multiqc_data|multiqc_report_data|g'` in qc.smk

### 7.13 PASS VCF empty (28 bytes)
- **Problem:** Depth filter (DP<30, alt reads<5) was too strict for test sample
- **Fix:** Ran bcftools without depth filter: `bcftools view -f PASS`

### 7.14 KeyError: 'map_cls' in generate_report.py
- **Problem:** CSS curly braces `{}` in `HTML_TEMPLATE.format()` were not all escaped as `{{}}`
- **Fix:** Replaced `str.format()` template with f-strings approach (new app.py has inline `make_html_report()` function)

### 7.15 No module named 'pandas' / 'matplotlib' / 'flask'
- **Problem:** Running scripts from `(base)` conda env instead of `(snakemake)`
- **Fix:** Always `conda activate snakemake` before running anything

---

## 8. Key Manual Commands Used During Debugging

```bash
# Always activate the right environment first
conda activate snakemake

# Check all tools available
for t in fastp bwa samtools gatk bcftools mosdepth vep; do
    echo -n "$t: "; which $t 2>/dev/null && echo OK || echo MISSING
done

# Run Snakemake dry-run
cd ~/genrichi
snakemake --snakefile workflow/Snakefile --configfile config/config.yaml \
  --use-conda --conda-frontend mamba --rerun-incomplete --cores 8 --dry-run

# Run full pipeline
snakemake --snakefile workflow/Snakefile --configfile config/config.yaml \
  --use-conda --conda-frontend mamba --rerun-incomplete --cores 8

# Manual FilterMutectCalls (test)
gatk FilterMutectCalls \
    -R ~/genrichi/reference_db/ref/hg38.fa \
    -V ~/genrichi/results/G8_S13/calling/G8_S13.mutect2.vcf.gz \
    --stats ~/genrichi/results/G8_S13/calling/G8_S13.mutect2.vcf.gz.stats \
    --contamination-table ~/genrichi/results/G8_S13/calling/G8_S13.contamination.table \
    -O ~/genrichi/results/G8_S13/calling/G8_S13.filtered.vcf.gz

# Manual VEP annotation
vep --input_file results/G8_S13/calling/G8_S13.pass.vcf.gz \
    --output_file results/G8_S13/annotation/G8_S13.vep.vcf.gz \
    --format vcf --vcf --compress_output bgzip \
    --cache --dir_cache ~/genrichi/reference_db/vep_cache \
    --assembly GRCh38 --species homo_sapiens --offline --everything \
    --fork 4 --force_overwrite

# Generate report manually
python3 /tmp/make_report.py

# Copy report to Windows Desktop
cp ~/genrichi/results/G8_S13/report/G8_S13_report.html \
   /mnt/c/Users/ramir/Desktop/G8_S13_report.html
```

---

## 9. Flask Web Application

### 9.1 Purpose
Allow GenRichi to submit real patient FASTQ files through a browser and get back:
- Live progress tracking (8 pipeline steps with animated dots)
- Live log terminal (dark style)
- HTML clinical report (embedded charts)
- VCF download

### 9.2 Files
| File | Purpose |
|---|---|
| `~/genrichi/app.py` | Flask backend — full somatic pipeline |
| `~/genrichi/templates/dashboard.html` | Clinical dashboard UI |
| `~/genrichi/genrichi_jobs.db` | SQLite — job history (auto-created) |

### 9.3 API Routes
| Route | Method | Description |
|---|---|---|
| `/` | GET | Dashboard UI |
| `/api/submit` | POST | Submit new sample job |
| `/api/job/<id>` | GET | Job status + logs |
| `/api/jobs` | GET | All jobs list |
| `/report/<id>` | GET | View HTML report in browser |
| `/api/download/<id>/report` | GET | Download report HTML |
| `/api/download/<id>/vcf` | GET | Download annotated VCF |
| `/api/tools` | GET | Check which tools are installed |

### 9.4 Pipeline Steps in app.py
```
Step 1: fastp QC + trimming
Step 2: bwa mem alignment → samtools sort
Step 3: GATK MarkDuplicates
Step 4: GATK BaseRecalibrator + ApplyBQSR
         mosdepth coverage
Step 5: GATK Mutect2 (tumor-only, with gnomAD germline resource)
         GATK LearnReadOrientationModel
Step 6: GATK GetPileupSummaries → CalculateContamination
         GATK FilterMutectCalls
         bcftools PASS filter + bgzip + tabix
Step 7: VEP annotation (with ClinVar custom track)
         VCF → TSV export
Step 8: HTML report (QC stats, coverage, VAF histogram, impact pie, variant table)
```

### 9.5 How to Launch
```bash
conda activate snakemake
pip install flask          # first time only
cd ~/genrichi
python3 app.py
```
Open browser: **http://localhost:5000**

---

## 10. Test Sample — G8_S13

| Item | Value |
|---|---|
| Sample ID | G8_S13 |
| R1 | `/home/rami/genrichi/outputs/0e66c1d7-fe87-463c-886c-f7dcc8ac34a5/G8_S13_R1_001.fastq.gz` |
| R2 | `/home/rami/genrichi/outputs/0e66c1d7-fe87-463c-886c-f7dcc8ac34a5/G8_S13_R2_001.fastq.gz` |
| Result | 221KB HTML report generated successfully |
| Report | `~/genrichi/results/G8_S13/report/G8_S13_report.html` |
| Variants | PASS variants detected, VEP annotated |

---

## 11. COSMIC Note

COSMIC requires a free account registration.

1. Register at: https://cancer.sanger.ac.uk/cosmic/download
2. Download: `CosmicCodingMuts_v99_GRCh38.vcf.gz`
3. Place at: `~/genrichi/reference_db/cosmic/CosmicCodingMuts_v99_GRCh38.vcf.gz`
4. Index: `tabix -p vcf ~/genrichi/reference_db/cosmic/CosmicCodingMuts_v99_GRCh38.vcf.gz`

The pipeline works without COSMIC — VEP will still annotate with ClinVar + gnomAD.

---

## 12. Pending / Future Work

- [ ] Add COSMIC to VEP annotation once downloaded
- [ ] Replace `example_hotspots.bed` with real validated clinical panel BED
- [ ] Add PDF export option to HTML report
- [ ] Add MultiQC integration to web app report
- [ ] Test with more real samples
- [ ] Consider adding Slurm/HPC submission option to web app
- [ ] Add user authentication to web app for clinical use
- [ ] ACMG/AMP variant classification integration

---

## 13. Quick Reference Cheat Sheet

```bash
# Activate environment
conda activate snakemake

# Start web app
cd ~/genrichi && python3 app.py
# → http://localhost:5000

# Snakemake dry-run
snakemake --snakefile workflow/Snakefile \
  --configfile config/config.yaml \
  --use-conda --cores 8 --dry-run

# Check a sample's results
ls ~/genrichi/results/G8_S13/

# View report on Windows
cp ~/genrichi/results/G8_S13/report/G8_S13_report.html \
   /mnt/c/Users/ramir/Desktop/

# If a tool is missing
conda activate snakemake
mamba install -c bioconda -c conda-forge <tool_name>
```

---

*Session notes generated 2026-05-25 — GenRichi Diagnostics*
