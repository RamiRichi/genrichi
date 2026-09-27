# HCC1395 / SEQC2 Preliminary Concordance Check — v1

**Status:** Preliminary, limited technical concordance check. **Not a sensitivity/specificity validation, not a clinical validation, not a check of the full panel or of CNV/MSI/TMB.**
**Pipeline commit under test:** `9f0af05` (Phase 1: build TMB CDS from MANE Select and pair it with the 55-gene BED).
**Run date:** 2026-09-26/27. **Run location:** `~/seqc2_phase1_test/` on the development machine, entirely outside `~/genrichi` (the production tree) and outside this repository's `results/`. No production file, database, or deployment was touched.
**Two runs are covered here:** `WES_FD` and `WES_LL`, both the same HCC1395/HCC1395BL biological pair sequenced at two different SEQC2 sites. **They are two technical replicates of one biological sample, not two independent samples** — see "WES_LL run and shared-mask comparison" below. Each run's own mask gives a different set of "evaluable" truth variants (each run has its own sequencing depth profile), so the two runs' individually-reported numbers are **not directly comparable position-by-position** without also computing a mask that both runs satisfy at once — which is what the shared-mask section does.
**Raw sequencing data (FASTQ/BAM) are NOT in this repository** — see "Data used" below for public source links and checksums instead.

## What this checks, and what it does not

This is a first, small-scale comparison of GenRichi Phase 1 (55-gene) SNV/indel calls against a public somatic reference-sample truth set, restricted to the small fraction of the panel where that truth set actually applies. It shows that the pipeline can reproduce two known truth calls with one unmatched call, on public reference material — **not** that its sensitivity or specificity are known, and **not** that CNV, MSI, TMB, indel calling, or chrX regions have been checked at all.

## Data used

- **Sample:** HCC1395 (tumor) / HCC1395BL (normal), the SEQC2 reference cell-line pair.
- **Sequencing:** SEQC2 cross-site WES study, `PRJNA489865`, SureSelect V6+UTR. Two site pairs, both HCC1395/HCC1395BL:
  - `WES_FD_1` (`SRR7890879` tumor / `SRR7890880` normal).
  - `WES_LL_1` (`SRR7890850` tumor / `SRR7890851` normal).
  FASTQ from ENA; MD5 verified against ENA's published checksums for both pairs (`docs/validation/HCC1395_SEQC2_v1_evidence/fastq_SOURCES.txt` + `fastq_expected_md5.txt` for FD, `fastq_SOURCES_LL.txt` + `fastq_expected_md5_LL.txt` for LL). **FASTQ files themselves are not stored here.**
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

## WES_LL run and shared-mask comparison

`WES_LL` (same HCC1395/HCC1395BL pair, different SEQC2 sequencing site) was run through the identical Phase 1 config and the identical `evaluate.py` (same script, same masks, unmodified) as a **technical/site replicate check**, not as an independent sample.

### WES_LL, its own primary mask (≥20x/≥20x)

`WES_LL` happens to have deeper sequencing at the two sites that were *not evaluable* under `WES_FD`'s own mask, so under **its own** mask all 4 truth SNVs are evaluable:

| Type | Evaluable truth variants | TP | FP | FN |
|---|---|---|---|---|
| SNV | 4 | 4 | 0 | 0 |
| Indel | 0 | — | 0 | — |

This is **not a better result than FD** — it reflects a deeper sequencing run at those two positions, not a difference in pipeline behaviour. Per-site depth comparison:

| Site | FD tumor/normal depth | LL tumor/normal depth |
|---|---|---|
| chr13:32339132 | 45 / 37 | 36 / 42 |
| chr17:7675088 (TP53) | 13 / 5 (excluded, primary) | 114 / 143 (evaluable) |
| chr2:47806320 | 37 / 54 | 36 / 57 |
| chr8:38428420 | 18 / 46 (excluded by `min_depth=20`) | 21 / 69 (just above 20x; evaluable) |

`WES_LL` also has one PASS call outside the SEQC2 high-confidence regions (`chr16:68823538 T>C`) — not evaluable, and **not comparable to `WES_FD`'s unmatched call** (different position, different run; see below).

### Shared-evaluable mask (both runs at once)

Because each run's own mask covers a different subset of positions, a **shared mask** — Phase 1 panel ∩ SEQC2 high-confidence regions ∩ {tumor ≥20x and normal ≥20x in **both** FD and LL} — was computed separately (`docs/validation/HCC1395_SEQC2_v1_evidence/compare_fd_ll.py`, standalone, does not modify `evaluate.py` or any pipeline/config file) to make the two runs directly comparable position-by-position.

| | Value |
|---|---|
| Shared mask | 1,322 intervals, **127,378 bp** |
| Truth SNVs inside shared mask | **2** (chr13:32339132, chr2:47806320 — TP53 and chr8 are excluded here because `WES_FD` alone does not meet ≥20x/≥20x at those two sites, regardless of `WES_LL`'s depth) |
| Truth indels inside shared mask | 0 |
| PASS calls inside shared mask, FD | 3 (chr13, chr2, chr16:23635069) |
| PASS calls inside shared mask, LL | 2 (chr13, chr2) |
| Called in both runs | chr13:32339132, chr2:47806320 — both match truth (**TP=2**) |
| Called in FD only | chr16:23635069 G>A |
| Called in LL only | none |

**Reading:** within the region both runs can actually be compared on, the two evaluable truth variants are detected identically in both runs. The one discrepancy is `WES_FD`'s unmatched call, which `WES_LL` does not reproduce.

### `chr16:23635069 G>A` (the FD-only unmatched call), investigated in the LL run

- FD depth at this exact position (via `samtools depth -Q20 -q20`): tumor 70, normal 43 — consistent with Mutect2's own reported tumor DP=68 (AD 63,5; VAF≈0.083) for this call.
- LL depth at the same exact position: tumor 106, normal 109 — deeper still, so lack of coverage cannot explain any difference.
- The position is **absent from all three of LL's VCF stages**: raw `mutect2.vcf.gz`, `filtered.vcf.gz`, and `pass.vcf.gz`. This means Mutect2 itself never proposed a candidate variant at this position in the LL run — it was not filtered out downstream, it was never called at all.
- **This cannot be attributed to a specific cause from the data examined here.** Possible explanations include a genuine low-level artifact specific to the FD sequencing/alignment (e.g. a strand- or batch-specific error), or a real but very low-level event not present (or below Mutect2's internal detection threshold) in the LL reads; distinguishing these would need a pileup-level look at the LL reads at this exact position, which was not done. It is **not** attributable to coverage.
- `WES_LL`'s own unrelated unmatched call, `chr16:68823538 T>C`, is a different position on the same chromosome and falls **outside** the SEQC2 high-confidence regions; it is unrelated to this call and is not itself evaluable.

## Software and settings

Same environments as the frozen Phase 1 baseline (`docs/PHASE1_BASELINE.md` §2–3) — no new conda environment was created for this run. Versions observed at run time: BWA 0.7.17, GATK 4.6.1.0, samtools 1.19.2, bcftools 1.19, VEP 113.0, fastp 0.23.4. Panel BED sha256 `cb0ce560f78334ea1e57974c6bf4d52643d9f7a445bc893a6eea2974d3265b36` — the same committed 55-gene Phase 1 BED, unmodified.

The Snakemake run initially failed at the VEP-annotation step on a transient network error (could not reach `ensembldb.ensembl.org:3306`); it was resumed — not restarted — once connectivity was confirmed, re-using all already-completed upstream jobs. Full run log kept in the local test folder (not committed; see "Reproducibility" below).

## Evidence kept in this repository

`docs/validation/HCC1395_SEQC2_v1_evidence/`:
- `SEQC2_HCC1395_FD1.metrics.json`, `SEQC2_HCC1395_LL1.metrics.json` — full machine-readable per-run results (all strata, all per-variant depths and classifications).
- `SEQC2_HCC1395_FD1.primary_FP_calls.tsv`, `SEQC2_HCC1395_LL1.primary_FP_calls.tsv` — each run's unmatched-call records at the primary mask (LL's is empty: 0 unmatched calls there).
- `fd_vs_ll_shared_mask_comparison.json` — the full shared-mask comparison result (regions, per-run call sets, agreement, and the chr16:23635069 investigation).
- `evaluate.py` + `evaluate.py.sha256_predeclared` — the per-run comparison script and a dated log of its design and of two bugs found and fixed in it *before* any numeric result was produced (a `bcftools view -R`-needs-an-index issue on piped input, and the coverage-mask correction to ≥20x/≥20x primary requested on 2026-09-27). The design (regions, strata, matching rule) was fixed in writing before results were seen.
- `compare_fd_ll.py` + `compare_fd_ll.py.sha256.txt` — the separate, standalone shared-mask comparison script (does not modify `evaluate.py`, its masks, or any pipeline/config file) and its dated hash.
- `truth_SOURCES_AND_CHECKSUMS.txt`, `fastq_SOURCES.txt`/`fastq_expected_md5.txt` (FD), `fastq_SOURCES_LL.txt`/`fastq_expected_md5_LL.txt` (LL) — provenance and checksums for every external file used. **No FASTQ, BAM, or VCF sequencing data file is included.**

## Reproducibility

Both full pipeline runs (all intermediate results, logs, and the `results/SEQC2_HCC1395_FD1/report/` and `results/SEQC2_HCC1395_LL1/report/` HTML reports) live only in the local, non-committed test directory `~/seqc2_phase1_test/` on the machine this was run on. They are not part of this repository and can be regenerated from the public sources listed above plus this commit's `workflow/` and `config/comprehensive_config.yaml`.

## Decisions and next steps (open, not yet done)

1. **`calling.filter.min_depth` stays at 20** for now. Before any change, define in writing: the pre-declared acceptance criteria (what depth/VAF a change would need to satisfy, on what evidence) — not a threshold change made in response to one borderline site. `chr8:38428420` remains documented as a site excluded by this threshold in both runs' own masks (and, in FD's case, in the shared mask too), not as a defect.
2. **Technical replicate — done.** `WES_LL` was run through the identical pipeline and check. Result: on the shared-evaluable mask, both runs agree on both evaluable truth SNVs (TP=2/2 in each); the one discrepancy is `WES_FD`'s unmatched `chr16:23635069` call, which `WES_LL` does not reproduce and which is not explained by coverage (see "WES_LL run and shared-mask comparison" above). This confirms `WES_LL` was correctly used as a **run/site reproducibility check**, not as a source of new reference-variant diversity — it added zero new truth variants and zero indels, as expected for the same biological pair.
3. **Variant diversity, especially indels:** still open. None of the material used so far (`WES_FD`, `WES_LL`, or the SEQC2 v1.2.1 high-confidence truth set) gives any evaluable indel truth inside the panel, and the two runs together still cover only **4 reference SNVs and 0 indels** in total. A different reference material (a different cell line or reference standard with known indels in panel genes) is needed before indel performance can be checked at all; this is the next, still-open step.
4. No claim of sensitivity, specificity, or clinical performance should be made from this check alone, for SNVs, indels, CNV, MSI, or TMB — the total evidence base remains 4 reference SNVs and 0 indels across two technical replicates of one biological sample.
