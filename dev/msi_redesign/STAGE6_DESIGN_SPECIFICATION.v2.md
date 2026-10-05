# Stage 6 Design Specification v2 — Empirical MSI Validation (development / RUO)

Status: DESIGN ONLY. Amends the approved Stage 6 design (v1, 19 sections) and the approved Stage 6A implementation.
No code, no downloads, no data requests, no panel membership, no panel size, no thresholds, no production change.

Retitle: **Stage 6 = Empirical MSI Validation.** The goal is defensible evidence for a future diagnostic product, treated as
clinical development. The word "validation" names the *goal*. Until 6B has run on independent, pre-registered data, no output, file
or document may say anything is "validated"; everything stays development / RUO.

## 0. Principles (product-owner requirements, binding)
P1 No locus is chosen because it looks good in one dataset.
P2 No threshold is set on the samples later used to establish performance.
P3 No label produced by GenRichi is used to prove GenRichi.
P4 Biological replicates are never mixed with independent samples.
P5 Every sample has explicit provenance and an explicit label source.
P6 Every later panel decision is traceable to evidence.
P7 Selection data are separated from an independent validation set.
P8 Insufficient results stay `NOT_ASSESSED` / `INSUFFICIENT`; never a guess.

Where the current implementation already satisfies them: P1, P2 (provisional parameters are labelled and untuned), P3
(`label_source` mandatory), P4 (one unit for all HCC1395 orders; content-duplicate detection), P8 (verified). Gaps closed by the
amendments below: P5 (fuller provenance), P6 (decision trace), P7 (locked validation set + pre-registration).

## A. Pre-registered analysis plan (closes P1, P2, P7)
1. Before any MSI-H data are opened, a **Stage 6B Analysis Plan** document is written and frozen (SHA-256 recorded in a
   `stage6b_plan.manifest.json`). It fixes: metrics and their definitions (Section 8 of v1), missing-data rules, minimum
   evidence per dimension, partition rules, multiplicity control, and the exact comparison tables to be produced.
2. Numeric minimums (units per class per partition; informative reads per locus) are **proposed by the product owner** and
   recorded in the plan; none is derived from HCC1395 or from the data being analysed.
3. Any change after the plan hash is recorded is a **deviation**: appended to a deviation log with date, reason and whether
   any outcome had been seen. A deviation made after seeing validation results invalidates the validation run.
4. The provisional prototype parameters (MIN_MAPQ 20, MIN_FLANK 5, MIN_LOCUS_DEPTH 20) remain provisional. They may be replaced
   only inside the plan, using **development-partition data only**, never validation data.

## B. Partitions and the locked validation set (closes P2, P4, P7)
1. Unit of partitioning = biological patient/sample (v1 Section 4). All files, orders, replicates and tumor/normal of one unit
   share one partition.
2. Partitions: **development** (metric exploration, locus evidence, parameter choice), **internal validation** (optional,
   used once to check the frozen development decisions), **independent validation** (ideally a different cohort / source /
   sequencing batch than development). If unit counts do not allow all three, the run states which are missing; a split is
   never invented.
3. Assignment is deterministic (seeded hash of `biological_unit_id`) and written to `split_assignment.tsv` **before**
   any metric is computed on labelled data.
4. The independent-validation partition is **locked**: opened once, on the frozen plan and frozen development decisions.
   Every access is written to an access log inside the manifest (who/what script/hash/time). A second look at validation data
   after a decision change turns that partition into development data; a new independent set is then required.
5. A verifier check refuses any selection, ranking, parameter or threshold step that reads validation-partition labels.
6. Technical replicates and re-processings (like the four HCC1395 orders) never count as independent samples; cohort
   overlap between sources (same patient in TCGA and CPTAC, cell lines shared between databases) is checked before assignment.

## C. Provenance schema for every sample (closes P3, P5)
Required registry fields per biological unit (missing required field = refuse to register):
`biological_unit_id`, `dataset_origin` (cohort / vendor / cell line), `source_url_or_accession`, `access_tier` (open /
controlled / commercial), `data_use_terms_checked` (date + who), `consent_deidentification_status` (for patient data),
`tumor_type`, `specimen_type` (fresh-frozen / FFPE / cell line / contrived), `assay` (WES / WGS / panel + capture kit),
`reference_build`, `aligner_and_version`, `matched_normal` (yes/no + id), `tumor_purity` (value + method, or
`NOT_ASSESSED`), `sequencing_depth`, `label` (`MSS` | `MSI-H` | `UNLABELLED`), `label_source`, `label_method`, `label_date`,
`label_independent_of_genrichi` (must be true for any unit used for performance).
Label-source hierarchy, recorded not assumed: (1) clinical reference methods (PCR/CE panel, MMR IHC) → (2) published expert
annotation → (3) published *computational* call from another tool (e.g. MANTIS, MSIsensor) → (4) GenRichi output (**forbidden
as ground truth**). Tier-3 labels are usable only if the publication and tool are recorded and the limitation "label is itself
a computational call" is carried into every result table.

## D. Decision trace (closes P6)
1. Every evidence row already carries candidate id, dataset, unit, metric definition and `manifest` hash. Add an
   **evidence id** (hash of candidate + metric + partition + input-file hashes + plan hash).
2. Any future decision record (a locus kept, dropped, or ranked differently; a panel composition proposal) must cite the
   evidence ids it relies on, the plan hash and the partition it used. A decision citing validation-partition evidence for
   *selection* is rejected by the verifier.
3. Decision records are append-only files with hashes; nothing is edited in place.
4. `NOT_ASSESSED` / `INSUFFICIENT` evidence can be cited only as "no evidence"; a decision may not treat it as support.

## E. Data-governance gate (before any patient-derived data enters the environment)
Nothing is requested or downloaded until, per source: licence / data-use terms read and recorded; access tier and any
data-access-committee requirement understood; consent and de-identification status known; storage location and retention
decided (outside Git); the product owner has approved that source by name. Controlled-access patient data (BAM/VCF) must not
be copied into the shared or production areas of the environment.

## F. Candidate independent data sources (NOT requested; each needs the gate in E)
Verified from public pages during this design step (not yet from the terms of use themselves):

| Source | What it offers | Access / caveat | Label-source class |
|---|---|---|---|
| TCGA (colorectal, endometrial, gastric, etc.) via GDC | Tumor + matched-normal WES with BAM/VCF; MSI status widely analysed with MANTIS and other tools (Bonneville 2017: 458 pairs, six cancer types; Cortes-Ciriano 2017) | BAM/VCF are controlled access (dbGaP authorisation); label is a published *computational* call, and some cohorts have PCR/IHC annotations — must be checked per sample | Tier 3 (some Tier 1/2, to verify) |
| CPTAC colon and endometrial (via GDC / Proteomic Data Commons) | WES and WGS of tumor and blood normal, harmonised to GRCh38; separate proteogenomic cohort from TCGA | GDC controlled VCF/BAM; MSI label source per case must be located and recorded | to be determined per case |
| Seraseq MSI reference materials (SeraCare / LGC) | Contrived MSI-High gDNA (cell-line-based) and MSI reference panel mixes at 5% and 20% allele fraction (BAT-25/26, NR-21/24, MONO-27) | Commercial purchase; contrived material, not patient tumor; suited to assay-level checks and limit-of-detection, not to cohort discrimination | Tier 1-like (manufacturer-declared) |
| MSI-H / MSS cancer cell lines with public sequencing (e.g. HCT116-type MSI-H, MSS lines) | Cheap positive/negative controls | Usually tumor-only (no matched normal), single-line n, labels from database annotation; useful for development, never independent validation | Tier 2 |
| Existing GenRichi orders (NA12878, G8_S13) | None for MSI | Unlabelled; not usable as ground truth | — |

Realistic reading: one source alone will not give a defensible validation. A defensible design most likely combines a patient
cohort with clinical-method labels (development + independent validation drawn from different cohorts) and contrived
reference material for analytical checks. Which sources to pursue, and the numeric sample-size minimums, are product-owner
decisions.

## G. What changes in the v1 deliverables
- v1 Sections 4, 10, 11, 12 are amended by A, B and D above.
- v1 Section 3 (ground truth) is extended by C.
- v1 Section 17 (manifest) gains: `plan_hash`, `partition_access_log`, `deviation_log_ref`, `evidence_id` scheme,
  provenance completeness per unit.
- v1 Section 19 (will-not-claim) gains: no statement of "validated", "clinical performance" or "regulatory readiness" before
  6B runs on locked independent data under a frozen plan, and none even then without the regulatory pathway.
- Stage 6A outputs are unchanged and remain technical characterisation only.

## H. Regulatory note
Treating this as clinical development implies design-control style records and an IVDR-compatible evidence trail. The IVDR
pathway has not been started (per project records). The evidence structure above is chosen so it can later be reused rather
than rebuilt; it is not itself a regulatory submission.

## I. Open decisions for the product owner (nothing proceeds without them)
1. Which data sources (by name) may go through the governance gate in E.
2. The numeric minimums: independent MSI-H and MSS units per partition, informative reads per locus, number of partitions.
3. Whether tier-3 (computational) labels are acceptable for development only, for validation, or not at all.
4. Who signs the frozen Stage 6B Analysis Plan.

Sources used for Section F: Bonneville et al. / MANTIS performance evaluation (PMC5352334); Cortes-Ciriano et al. pan-cancer MSI
portrait (PMC5467167); GDC release of CPTAC genomic data (proteomics.cancer.gov); CPTAC Pan-Cancer data (Proteomic Data
Commons); SeraCare Seraseq MSI reference materials (seracare.com).
