# HCC1395 / SEQC2 Preliminary Concordance Check — v1

**Status:** Preliminary, limited technical concordance check. **Not a sensitivity/specificity validation, not a clinical validation, not a check of the full panel or of CNV/MSI/TMB.**
**Pipeline commit under test:** `9f0af05` (Phase 1: build TMB CDS from MANE Select and pair it with the 55-gene BED).
**Run date:** 2026-09-26/27. **Run location:** `~/seqc2_phase1_test/` on the development machine, entirely outside `~/genrichi` (the production tree) and outside this repository's `results/`. No production file, database, or deployment was touched.
**Raw sequencing data (FASTQ/BAM) are NOT in this repository** — see "Data used" below for public source links and checksums instead.

## What this checks, and what it does not

This is a first, small-scale comparison of GenRichi Phase 1 (55-gene) SNV/indel calls against a public somatic reference-sample truth set, restricted to the small fraction of the panel where that truth set actually applies. It shows that the pipeline can reproduce two known truth calls with one unmatched call, on public reference material — **not** that its sensitivity or specificity are known, and **not** that CNV, MSI, TMB, indel calling, or chrX regions have been checked at all.

## Data used

- **Sample:** HCC1395 (tumor) / HCC1395BL (normal), the SEQC2 reference cell-line pair.
- **Sequencing:** SEQC2 cross-site WES study, site pair `WES_FD_1` (`SRR7890879` tumor / `SRR7890880` normal), SureSelect V6+UTR, `PRJNA489865`. FASTQ from ENA; MD5 verified against ENA's published checksums (`docs/validation/HCC1395_SEQC2_v1_evidence/fastq_SOURCES.txt`, `fastq_expected_md5.txt`). **FASTQ files themselves are not stored here.**
- **Truth set:** SEQC2 Somatic Mutation Working Group release v1.2.1 — `high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz`, `high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz`, `High-Confidence_Regions_v1.2.bed` — downloaded from NCBI (`ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest/`), with source URLs, sizes and sha256/md5 recorded in `docs/validation/HCC1395_SEQC2_v1_evidence/truth_SOURCES_AND_CHECKSUMS.txt`.
- **Reference build check:** all contig names and lengths in both truth VCFs and the high-confidence BED were checked against `reference_db/ref/hg38.fa.fai` — full match, no renaming or coordinate shift needed (checked for every contig, not just the first).

## Evaluable region

- Phase 1 panel BED (55 genes, merged): 901 intervals, 200,785 bp.
- ∩ SEQC2 High-Confidence Regions v1.2: **873 intervals, 192,493 bp** — this intersection, not the full panel, is the only region this check speaks to.
- The panel's chrX regions fall entirely outside the SEQC2 high-confidence regions (autosomes only) and are **not evaluable** by this truth set at all.

## Truth variants inside the evaluable region

| Type | Count |
|---|---|
| SNV | 4 |
| Indel | **0** |

**No indel is evaluable by this check.** The four SNVs are TP53 (chr17, the classic HCC1395 R175H hotspot), and one variant each on chr2, chr8 and chr13.

## Primary comparison mask and result

Primary mask (fixed for this check, matching the pipeline's own `calling.filter.min_depth=20`): evaluable region ∩ {tumor depth ≥20 **and** normal depth ≥20}, depth from `samtools depth -Q20 -q20` on the pipeline's own final BAMs.

Under this mask, **2 of the 4 truth SNVs are evaluable**; the other two are excluded from this stratum as *not evaluable* (not counted as false negatives) because they do not meet the depth mask — see per-variant table below.

| Type | Evaluable truth variants | TP | FP | FN |
|---|---|---|---|---|
| SNV | **2** | **2** | **1** | 0 |
| Indel | 0 | — | 0 | — |

**Reading of this result:** at the primary mask, the two evaluable reference SNVs were both detected (2/2), and one additional PASS call did not match the truth set. With only two evaluable truth variants and one unmatched call among three total calls, **this is a limited, preliminary technical concordance, not an estimate of sensitivity or specificity** — the sample size is far too small to support such a claim, and it must not be described as "good accuracy" or similar. A formal illustrative Wilson 95% interval on 2/2 detected spans 0.34–1.0, which by itself shows how little this sample size constrains the true rate.

**The one unmatched (FP) call:** `chr16:23635069 G>A`, PASS by the pipeline's own filters, tumor VAF ≈8.3% (5/68 reads), not present in the SEQC2 high-confidence truth set. Because that truth set is itself a consensus superset (not an exhaustive list of every real variant), this call is recorded as *unmatched to the evaluable truth*, not asserted to be a genuine false call.

### Per-truth-variant detail

| Site | Type | Tumor depth | Normal depth | ≥10x/≥10x eligible | ≥20x/≥20x eligible (primary) | In final `pass.vcf.gz`? | Classification (primary, 20x/20x) |
|---|---|---|---|---|---|---|---|
| chr13:32339132 G>T | SNV | 45 | 37 | yes | yes | yes | **TP** |
| chr17:7675088 C>T (TP53) | SNV | 13 | 5 | no (normal <10) | no | no | not evaluable (insufficient coverage) |
| chr2:47806320 G>A | SNV | 37 | 54 | yes | yes | yes | **TP** |
| chr8:38428420 G>A | SNV | 18 | 46 | yes | no (tumor <20) | **no — excluded by `min_depth=20`** (Mutect2 itself marked this PASS, tumor AD 0,18, VAF 0.92) | not evaluable at this mask |

`chr8:38428420` is the clearest concrete illustration of the pipeline's `min_depth=20` filter acting at a borderline depth (tumor=18): Mutect2's own filter says PASS, but the pipeline's downstream depth threshold drops it before it reaches `pass.vcf.gz`. **This is being recorded as a known, deliberate consequence of the current threshold — not as a pipeline defect, and not as grounds to change the threshold from a single site.** `calling.filter.min_depth` stays at **20** until formal, pre-declared acceptance criteria for such a change are defined; see "Decisions and next steps" below.

At the more permissive ≥10x/≥10x additional stratum (176,476 bp evaluable), `chr17` (TP53) remains *not evaluable* (normal depth 5 <10), while `chr8` becomes evaluable and is the one FN in that stratum (TP=2, FP=1, FN=1). Full per-stratum numbers (including the ≥30x-tumor/≥20x-normal stratum and the no-coverage-filter transparency stratum) are in `HCC1395_SEQC2_v1_evidence/SEQC2_HCC1395_FD1.metrics.json`.

## Software and settings

Same environments as the frozen Phase 1 baseline (`docs/PHASE1_BASELINE.md` §2–3) — no new conda environment was created for this run. Versions observed at run time: BWA 0.7.17, GATK 4.6.1.0, samtools 1.19.2, bcftools 1.19, VEP 113.0, fastp 0.23.4. Panel BED sha256 `cb0ce560f78334ea1e57974c6bf4d52643d9f7a445bc893a6eea2974d3265b36` — the same committed 55-gene Phase 1 BED, unmodified.

The Snakemake run initially failed at the VEP-annotation step on a transient network error (could not reach `ensembldb.ensembl.org:3306`); it was resumed — not restarted — once connectivity was confirmed, re-using all already-completed upstream jobs. Full run log kept in the local test folder (not committed; see "Reproducibility" below).

## Evidence kept in this repository

`docs/validation/HCC1395_SEQC2_v1_evidence/`:
- `SEQC2_HCC1395_FD1.metrics.json` — full machine-readable result (all strata, all per-variant depths and classifications).
- `SEQC2_HCC1395_FD1.primary_FP_calls.tsv` — the one unmatched call's full record.
- `evaluate.py` + `evaluate.py.sha256_predeclared` — the comparison script and a dated log of its design and of two bugs found and fixed in it *before* any numeric result was produced (a `bcftools view -R`-needs-an-index issue on piped input, and the coverage-mask correction to ≥20x/≥20x primary requested on 2026-09-27). The design (regions, strata, matching rule) was fixed in writing before results were seen.
- `truth_SOURCES_AND_CHECKSUMS.txt`, `fastq_SOURCES.txt`, `fastq_expected_md5.txt` — provenance and checksums for every external file used. **No FASTQ, BAM, or VCF sequencing data file is included.**

## Reproducibility

The full pipeline run (all intermediate results, logs, and the `results/SEQC2_HCC1395_FD1/report/` HTML report) lives only in the local, non-committed test directory `~/seqc2_phase1_test/` on the machine this was run on. It is not part of this repository and can be regenerated from the public sources listed above plus this commit's `workflow/` and `config/comprehensive_config.yaml`.

## Decisions and next steps (open, not yet done)

1. **`calling.filter.min_depth` stays at 20** for now. Before any change, define in writing: the pre-declared acceptance criteria (what depth/VAF a change would need to satisfy, on what evidence) — not a threshold change made in response to one borderline site.
2. **Technical replicate:** run `WES_LL` (`SRR7890850`/`SRR7890851`, same HCC1395/HCC1395BL biological pair, different sequencing site) through the same check, to see whether the two evaluable SNVs and the one unmatched call reproduce across a run/site replicate. This is a **run/site reproducibility check, not a source of new reference-variant diversity** — it is the same biological pair, so it is not expected to add new truth variants or indels.
3. **Variant diversity, especially indels:** none of the above materials give any evaluable indel truth inside the panel. A different reference material (a different cell line or reference standard with known indels in panel genes) is needed before indel performance can be checked at all; this is a separate, still-open search.
4. No claim of sensitivity, specificity, or clinical performance should be made from this check alone, for SNVs, indels, CNV, MSI, or TMB.
