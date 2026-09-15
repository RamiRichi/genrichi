# GenRichi Phase 1 — Software Baseline & Validation Record

**Status:** Computational/software migration validated end-to-end. **Not a clinical validation record** (see §6).
**Scope:** Phase 1 Solid Tumor SNV/Indel analysis only.
**Baseline frozen:** 2026-09-15, following the Miniforge3/Snakemake-managed conda environment migration.

This document records the exact software stack, environment identities, and validation evidence for the Phase 1 pipeline as of this freeze. It supersedes ad hoc `~/.local/bin` wrappers and system (`/usr/bin`) tool resolution as the basis for Phase 1 execution going forward.

---

## 1. Phase 1 Scope

**Analysis type:** Solid Tumor SNV/Indel (single-nucleotide variant and small insertion/deletion) sequence-level analysis, tumor/normal paired.

**Panel:** GenRichi Comprehensive Cancer Panel v1 (`panel.name` in `config/comprehensive_config.yaml`), authoritative Phase 1 solid-tumor target design.

| Property | Value |
|---|---|
| Gene count | 55 genes (verified — see `resources/panel/phase1_solid_tumor/phase1_audit.json`) |
| Transcript source | NCBI/Ensembl **MANE Select v1.4**, exactly one MANE_Select-tagged transcript per gene |
| Region definition | All CDS exons of the selected MANE Select transcript, each padded **±20 bp** (splice flank), merged within-gene only |
| Assembly | **GRCh38**, `chr`-prefixed contigs (chr1–22, chrX, chrY, chrM) |
| Total unique target bp | **200,785 bp** |
| Total BED intervals | **901** (BED4: chrom, start, end, gene_symbol) |
| BED file | `resources/panel/phase1_solid_tumor/solid_tumor_phase1_v1.bed` |
| Design provenance | `resources/panel/phase1_solid_tumor/README.md`, `phase1_audit.json` |
| Wired into config | Yes — `config/comprehensive_config.yaml` → `panel.bed` |

This BED is the panel actually used by `base_recalibrator_paired` (BQSR `-L`), `mosdepth_paired` (`--by`), and downstream coverage/CNV/TMB logic in the validated run.

**Entry point:** `workflow/Snakefile_comprehensive` (default `configfile: config/comprehensive_config.yaml`) — **not** `workflow/Snakefile` (which serves a separate, earlier germline/tumor-only pathway; see `docs/somatic_panel/session_notes.md` for that history).

---

## 2. Validated Software Stack

**Orchestrator:** Snakemake **9.21.0**, from the Miniforge3 distribution (`~/miniforge3/envs/snakemake`).
**Conda distribution:** Miniforge3 only (`~/miniforge3`). Every Phase 1 rule declares a `conda:` directive resolving to a Snakemake-managed, per-rule environment under `.snakemake/conda/`. No Phase 1 rule depends on `~/.local/bin` wrapper scripts or bare `/usr/bin` tool resolution.

| Env | Package | Version |
|---|---|---|
| **align** (`workflow/envs/align.yaml`) | BWA | 0.7.17 (classic — not bwa-mem2) |
| | samtools | 1.19.2 |
| | htslib | 1.21 (transitive; no explicit pin — see §7 rationale) |
| | GATK | 4.6.1.0 |
| | mosdepth | 0.3.10 |
| **calling** (`workflow/envs/calling.yaml`) | GATK | 4.6.1.0 |
| | bcftools | 1.19 |
| | htslib / tabix / bgzip | 1.19.1 |
| **qc** (`workflow/envs/qc.yaml`) | fastp | 0.23.4 |
| | MultiQC | 1.34 |
| **annotation** (`workflow/envs/annotation.yaml`) | Ensembl VEP | 113.x, cache version 113 (matched — see §7) |
| | bcftools | 1.19 |
| | htslib / tabix | 1.19.1 |
| | Python | 3.11.9 |
| | pandas | 2.2.3 |
| **report** (`workflow/envs/report.yaml`) | Python | 3.11.16 |
| | pandas | 2.2.3 |
| | numpy | 2.4.6 |
| | matplotlib | 3.9.4 |
| | jinja2 | 3.1.6 |

`cnv.yaml` and `msi.yaml` exist in `workflow/envs/` but are **not used** by any Phase 1 rule — `calculate_cnv` and `calculate_msi` are pure-Python implementations running under the `report` environment (no cnvkit/msisensor-pro dependency). Left untouched; out of scope for this freeze.

---

## 3. Snakemake Environment Hashes / Prefixes

Computed by Snakemake's own environment-hashing (`md5(realpath(--conda-prefix) + env yaml file content)`), rooted at `--conda-prefix .snakemake/conda`. Each hash is content-addressed to its exact `workflow/envs/*.yaml` file as frozen in this baseline.

| Env | Hash | Prefix |
|---|---|---|
| align | `b3e920885114bc7c28eb9e53d02501ec` | `.snakemake/conda/b3e920885114bc7c28eb9e53d02501ec_` |
| calling | `5a6a3d0b5d257e1cf1abb4ba0c8e5a66` | `.snakemake/conda/5a6a3d0b5d257e1cf1abb4ba0c8e5a66_` |
| qc | `5f7cf5f4c8def5afbc6ac17c69b5bd12` | `.snakemake/conda/5f7cf5f4c8def5afbc6ac17c69b5bd12_` |
| annotation | `c9a29be68648ef0994638fc9b534cd58` | `.snakemake/conda/c9a29be68648ef0994638fc9b534cd58_` |
| report | `07a78f4301110dfbbfa9784fb4680b9e` | `.snakemake/conda/07a78f4301110dfbbfa9784fb4680b9e_` |

If any `workflow/envs/*.yaml` file changes, its hash (and prefix) changes accordingly and Snakemake will create a new environment rather than silently reusing this one — this table is only valid for the exact YAML content frozen at this baseline.

---

## 4. Migration Validation Run

| Item | Value |
|---|---|
| Sample | `HCC1395_demo` (tumor/normal paired) only — no other sample in the DAG |
| Command basis | `snakemake --snakefile workflow/Snakefile_comprehensive --use-conda --conda-prefix .snakemake/conda --cores 4 --forceall results/HCC1395_demo/report/HCC1395_demo_comprehensive_report.html` |
| Cores | 4 |
| Jobs | **26 / 26 completed successfully** |
| Conda deployment | Every job logged `Activating conda environment: .snakemake/conda/<hash>_` against one of the 5 hashes in §3 |
| `/usr/bin` or `~/.local/bin` fallback | **None found** — full execution log searched, zero matches |
| Reference baseline preserved | `results_baseline_hcc1395_20260914/` (read-only copy of the original, pre-migration successful run + its Snakemake log) |
| New run output | `results/HCC1395_demo/` (current, conda-managed) |

`--forceall` was required because most upstream rule outputs (fastp, bwa_mem, markdup, BQSR, Mutect2, etc.) carried stale/absent provenance metadata from the original non-conda run and would not otherwise have been re-triggered by Snakemake's normal staleness detection.

---

## 5. Validation Checkpoints (new conda-managed run vs. `results_baseline_hcc1395_20260914/`)

| Checkpoint | Result |
|---|---|
| PASS somatic variant count | Match (1 PASS variant) |
| TP53 p.Arg175His — presence, VAF, depth, FILTER | Byte-identical (VAF 0.988, depth 83/82, FILTER=PASS) |
| BRCA1 intronic variant exclusion | Match (0 BRCA1 rows in both — correctly excluded by Phase 1 BED) |
| Tumor/normal pairing | Correct — verified via wildcard branches and matched flagstat read counts (no swap/cross-contamination) |
| Flagstat (tumor + normal) | Identical read counts, mapping %, duplicates, properly-paired % |
| mosdepth coverage (tumor + normal, total + on-target) | Identical |
| MSI output | Identical (1.5 sites, 0 somatic, 0.0%) |
| CNV output (`call.cns`) | Byte-identical |
| VEP annotation outcome | Identical (variants processed, filtered, overlapped genes/transcripts/regulatory features); VEP binary version now correctly 113.x (was 116.2 against a 113 cache pre-migration — see §7) |
| Final report generation | Generated successfully in both; only non-scientific differences (report-date timestamp, matplotlib-version metadata embedded in PNGs) |

---

## 6. Validation Scope — Explicit Distinction

**A. Computational / software validation: COMPLETED.**
This freeze certifies that the Phase 1 pipeline, running entirely on Miniforge3-managed, Snakemake-activated conda environments (no system or wrapper-script tool resolution), reproduces — scientifically, not just superficially — the previously-obtained HCC1395_demo result. This is a software/infrastructure migration validation.

**B. Clinical assay validation: NOT COMPLETED.**
This freeze is **not** a clinical validation of the Phase 1 assay. No claim is made here about analytical sensitivity/specificity, limit of detection, reproducibility across runs/operators/lots, or diagnostic performance. See §7 for what remains outstanding before any clinical claim can be made.

---

## 7. Known Limitations (carried forward, not resolved by this freeze)

- **No clinical Panel-of-Normals (PoN) validation yet** — `ref.pon` is `null` in `config/comprehensive_config.yaml`.
- **No orthogonal truth-set validation** — this run reproduces a prior pipeline result on HCC1395; it has not been benchmarked against an independent truth set (e.g. SEQC2/HCC1395 consensus calls, GIAB).
- **MSI is not clinically validated** — the MSI module is a panel-based indel-rate estimate (`msi.threshold` in config), explicitly documented in `workflow/rules/msi_scoring.smk` as a heuristic, not a clinically validated MSI assay.
- **CNV is not clinically validated** — `calculate_cnv` is a pure-Python mosdepth-ratio method, not a validated clinical CNV caller.
- **Fusion/structural rearrangement detection is not part of Phase 1.**
- **TERT promoter is not part of Phase 1** — the Phase 1 BED covers the TERT gene's CDS (for SNV/indel calling within coding sequence) but explicitly **excludes** the TERT promoter region; promoter hotspot mutations (e.g. C228T/C250T) are not detectable by this design.
- **EPCAM structural deletion detection is not part of Phase 1** — EPCAM CDS is covered for SNV/indel purposes only; EPCAM deletion (relevant to Lynch syndrome via MSH2 silencing) is a structural-variant capability not present here.
- **MANE_Select-only design** — 7 genes (BRAF, CDKN2A, FGFR2, HRAS, KRAS, MUTYH, SMARCA4) have an additional MANE_Plus_Clinical transcript not covered by this BED; see `resources/panel/phase1_solid_tumor/README.md` for the full rationale.
- Genes on this panel whose clinically relevant alterations may involve CNV, fusion, or other structural/regulatory regions (e.g. MET exon 14 skipping beyond standard splice flank, EGFR/ERBB2/MYC amplification, ALK/RET/ROS1-class fusions) **must not be represented as fully covered clinical capabilities** by this Phase 1 SNV/Indel design.
- **Phase 1 is sequence-level SNV/Indel analysis only.** Any broader claim (CNV, fusion, MSI-clinical, TERT promoter, structural variants) is out of scope until a separate, explicitly validated phase addresses it.

---

## References

- Panel design provenance: `resources/panel/phase1_solid_tumor/README.md`, `resources/panel/phase1_solid_tumor/phase1_audit.json`
- Pipeline configuration: `config/comprehensive_config.yaml`, `config/comprehensive_samples.tsv`
- Environment definitions: `workflow/envs/{align,calling,qc,annotation,report}.yaml`
- Preserved pre-migration baseline: `results_baseline_hcc1395_20260914/` (read-only)
- Earlier, separate-scope development log (20-gene tumor-only hotspot panel, different Snakefile/BED): `docs/somatic_panel/session_notes.md` — historical context only, not the Phase 1 baseline described here.
