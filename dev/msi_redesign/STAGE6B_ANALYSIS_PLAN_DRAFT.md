# Stage 6B Analysis Plan — Conditional Draft

**Status: DRAFT FOR PRODUCT-OWNER REVIEW — NOT APPROVED, NOT FROZEN, NOT REGISTERED.**  
Date: 2026-10-02.  
Scope: planning only. No source is selected, no cohort is constructed, no patient or controlled-access data are accessed, no labels are inspected, and no Stage 6B analysis is run by this document. This draft is conditional on governance closure and a later role-specific Source Selection Review.

## 1. Purpose and permitted conclusion

Stage 6B is intended to estimate how well the MSI redesign's candidate locus-based method distinguishes independently labelled MSS from MSI-H specimens. Until the development pipeline is completely verified, an eligible authorized source, and a locked cohort exist, the analysis is **not executable**. Any eventual result describes only the evaluated specimens, assay, and analysis conditions; it does not by itself establish clinical validity, clinical utility, diagnostic readiness, or regulatory acceptance.

Stage 6A data, including HCC1395 and any data used to choose loci, metrics, parameters, or thresholds, are development evidence and cannot serve as locked validation. The development directory contains marker scanners (`scan_markers.py`, `scan_markers_v2.py`), a CIGAR-aware BAM reader (`bam_extract.py`), a runner (`run_msi_dev.py`), and separate synthetic unit tests. The runner currently expects the v1 marker TSV schema (including `gene`); v2 is genome-wide and omits that field, so version, scope, and schema compatibility require an explicit decision. Synthetic integration tests now pass a marker TSV through the runner, SAM parsing, read filters, CIGAR repeat-length extraction, classifier, and report writer while simulating the external samtools response. The samtools-backed synthetic BAM tests remain skipped where samtools is unavailable, so real samtools I/O is not verified in this environment.

## 2. Preconditions — all must pass before data access

For each proposed source and role, record evidence and obtain the required owner decisions:

1. Close applicable data-use, requester-eligibility, access-route, consent, and intended-use questions in the Stage 6 governance tracker. Public availability or commercial-product wording alone is not authorization.
2. Complete a role-specific Source Selection Review. Document a lawful route, current terms, permitted intended use, responsible access account, storage/retention controls, and sample provenance.
3. Establish independent reference labels and their per-specimen provenance. Labels produced by GenRichi are never ground truth. Computational MSI calls alone are supporting evidence, not reference labels, under the current source-selection rules.
4. Resolve cohort-to-accession mapping, patient-level duplicate structure, and overlap with every development or candidate validation source. Unknown overlap remains `NOT_ASSESSED`; it cannot support an independence claim.
5. Identify an independent owner/reviewer and obtain written approval of source, intended role, and this plan before opening any validation data.
6. Complete the technical verification gate: select and document the intended marker-file version, scope, reference, panel, and schema; verify the scanner-to-reader-to-classifier handoff with fabricated inputs; and run the synthetic BAM I/O test when a samtools-enabled development environment is available. Do not treat these checks as evidence of biological accuracy.

Failure of any precondition means **PAUSE**. No analysis proceeds by substituting a more accessible source or by assuming missing facts.

## 3. Analysis question, units, and labels

- **Primary question:** among eligible, independently labelled specimens with sufficient assay data, how does the frozen GenRichi MSI method classify MSS and MSI-H specimens?
- **Biological unit:** patient. All specimens, tumour/normal pairs, aliquots, tumour pieces, technical replicates, sequencing runs, and reprocessings from one patient stay together and count as one independent unit.
- **Classes:** preserve source labels as `MSS`, `MSI-H`, `INTERMEDIATE/OTHER`, `INDETERMINATE`, or `UNLABELLED`. Do not silently merge intermediate or indeterminate cases into either primary class. Only MSS and MSI-H enter the primary binary comparison; other states are reported separately with counts and reasons.
- **Reference label:** retain all available evidence and provenance per patient. Use the predeclared source adjudication rule; keep discordant evidence visible. The rule and labels must be frozen before locked validation is opened.
- **Proposed discordance handling (owner confirmation required):** the target construct is MSI status by a qualified molecular MSI reference assay. A technically valid, appropriately documented PCR/fragment-analysis result is not overridden by MMR-IHC, because IHC measures protein expression/MMR status and is orthogonal evidence. Record IHC disagreement as discordance and report it separately. If two qualified molecular MSI reference assays disagree, do not resolve by majority vote or by GenRichi output: hold the case out of the primary reference-labelled analysis as `INDETERMINATE` pending independent, blinded adjudication under the reference laboratory's documented SOP; retain it in the flow and report a sensitivity analysis if adjudicated labels become available before lock. POLE-associated hypermutation or other biological context may be recorded as a covariate, but must not silently redefine the PCR-based MSI label.
- **Analysis record:** retain only the minimum authorized pseudonymous patient key needed for grouping, split assignment, and duplicate/overlap checks. Identifiers and controlled files must remain in their approved storage location, outside Git and shared production areas.

## 4. Data partitions and leakage controls

Use the design's three conceptual roles: `DEVELOPMENT`, optional `INTERNAL_VALIDATION`, and `LOCKED_INDEPENDENT_VALIDATION`. Do not force three partitions if the eligible cohort cannot support them; report the absent role and stop short of a validation claim.

- Stage 6A and any data used to develop the method remain `DEVELOPMENT`.
- A locked validation source must be independent at patient level from development data. If cross-source overlap cannot be resolved, the sources may not be assigned opposite development/validation roles.
- Assign whole patients deterministically using a predeclared, label-blind rule. Record the rule and assignments before label-dependent calculations. Keep each patient's tumour and matched normal and all replicates in one partition.
- Freeze membership, exclusions, labels and provenance, method version, parameters, thresholds, missingness rules, and plan before opening locked validation. Log each authorized access and its purpose.
- Any tuning or decision informed by locked-validation outcomes converts the affected data to development evidence; a new independent validation set is then required.

**Partition proportions, randomization seed/hash rule, and minimum patient counts per class are intentionally unresolved.** The owner must set and justify them from feasibility and the intended precision before the cohort is opened; they must not be chosen after inspecting validation outcomes. Do not claim adequate power or precision until a justified sample-size rationale is recorded.

## 5. Frozen method and parameter handling

Use the development classification core in `dev/msi_redesign/msi_locus_model.py` as the candidate classifier, identified by repository revision and file hashes at plan freeze. At actual freeze, record the Git commit hash **and SHA-256 hashes** for at least `dev/msi_redesign/msi_locus_model.py`, `dev/msi_redesign/bam_extract.py`, the selected marker scanner, and `dev/msi_redesign/run_msi_dev.py`; also hash the selected marker TSV, reference FASTA, panel BED (for v1), and relevant configuration. Values without pinned code and input versions do not define a reproducible classifier. This is not the Stage 6A empirical-evidence program and is not the production MSI implementation. The supporting scanner, reader, and runner exist; v1 is the currently compatible runner input, while v2 needs schema/scope resolution before it can be used. Synthetic unit and runner integration checks are in place; samtools-backed synthetic I/O remains to be run in an environment with samtools. Do not alter production MSI code as part of Stage 6B.

The classifier's existing provisional values are `MIN_LOCUS_DEPTH=20`, `INSTABILITY_DISTANCE_CUTOFF=0.20`, `MIN_SHIFT_UNITS=1`, `MIN_INFORMATIVE_LOCI=5`, and `MSI_H_FRACTION_THRESHOLD=0.40`. The read-extraction filters `MIN_MAPQ=20` and `MIN_FLANK=5` belong to the separate Stage 6A evidence code and are not automatically Stage 6B parameters. All relevant extraction and classifier values remain provisional until the owner explicitly freezes them. Any permissible selection or tuning must use development data only and be completed before locked validation access. Record exact values and method versions in the frozen manifest; this draft does not ratify or calibrate them.

For a tied read-length mode, the current development implementation now deterministically selects the smaller repeat length so locus calls do not depend on read order. This tie-break is a provisional algorithm choice and must be reviewed and frozen with the other classifier rules before any labelled cohort is opened.

The current per-locus rule requires **both** a modal repeat-length shift (`shift >= MIN_SHIFT_UNITS`) **and** a tumor/normal total-variation distance (`distance >= INSTABILITY_DISTANCE_CUTOFF`) to call a locus `UNSTABLE`. Consequently, a synthetic 50:50 tumor distribution of repeat lengths 10 and 11 against a normal distribution entirely at length 10 has `mode(tumor)=10` under the smaller-length tie-break, `shift=0`, and `TVD=0.50`; the locus is therefore `STABLE` under the current conjunction, despite half of tumor reads showing the longer repeat. This is a conservative behavior that can miss bimodal distributions. The conjunction is provisional and must be explicitly reviewed and frozen before labelled-data access; replacing `AND` with `OR` would change the classifier and its sensitivity/specificity and is not implied or authorized by this note.

## 6. Outcomes and reporting

Report results by partition and reference-label tier, without pooling patient replicates as independent observations:

1. A patient-level confusion matrix using the classifier's actual calls: `VALID_MSS`, `VALID_MSI_H`, `INDETERMINATE`, and `INSUFFICIENT_DATA`. Preserve `NOT_ASSESSED` separately for records where analysis could not be performed; do not force any state into MSS or MSI-H.
2. Primary sensitivity = reference MSI-H patients called `VALID_MSI_H` divided by **all eligible reference MSI-H patients**; primary specificity = reference MSS patients called `VALID_MSS` divided by **all eligible reference MSS patients**. Thus, an `INDETERMINATE` or `INSUFFICIENT_DATA` result remains in the denominator and is a non-success for its corresponding primary endpoint. Report both with two-sided 95% confidence intervals and numerator/denominator. Also report valid-call-only sensitivity/specificity as secondary, explicitly conditional on a valid call, so no-call cases cannot be hidden by a reduced denominator.
3. Overall correct classification and balanced accuracy with confidence intervals, clearly secondary to sensitivity and specificity.
4. Number and fraction of `INDETERMINATE`, `INSUFFICIENT_DATA`, and `NOT_ASSESSED` cases, plus reasons and the full eligible-cohort denominators. Any exclusion before the eligible-cohort denominator must follow a frozen, outcome-blind rule and be shown in the flow table.
5. Per-patient and per-locus coverage/informativeness summaries; report locus-level missingness and class-specific evaluability.
6. Stratified descriptive summaries by label evidence tier, tumour type, assay, preservation, purity/cellularity, and sequencing quality where metadata permit. These are exploratory unless separately powered and predeclared.

The statistical implementation, exact confidence-interval method, handling of multiple secondary comparisons, and clinically meaningful acceptance targets remain **owner decisions** and must be fixed before lock. No post-hoc threshold, subgroup, or exclusion selection is allowed. An inconclusive interval or insufficient evaluable count is reported as inconclusive, not as success.

## 7. Missingness, exclusions, and deviations

- Register all candidate units before exclusions with explicit reasons; preserve an auditable flow count from candidate through analysis.
- Exclude a unit from the primary evaluable set only under predeclared, label-blind technical/data-quality rules. Record the rule, evidence, and counts by class after the rule is frozen.
- Missing required provenance, unresolved label conflict, inadequate repeat-locus evidence, or unresolved patient duplication yields `NOT_ASSESSED`/`INSUFFICIENT` or a predeclared exclusion; never impute a class or infer a missing value.
- Do not revise labels, membership, exclusions, thresholds, or metrics after locked outcomes are seen. Record every deviation with date, rationale, approver, affected records, and whether outcomes were visible. A material post-outcome change invalidates the original locked-validation claim.

## 8. Multiplicity and interpretation

The primary comparison is one patient-level MSS-versus-MSI-H analysis using the frozen method and primary metrics. All secondary locus, subgroup, or alternative-metric analyses are labelled exploratory; apply a predeclared multiplicity procedure if formal inferential claims are planned. The procedure and family of comparisons must be specified before lock. Do not select a “best” locus or cutoff using locked-validation labels.

## 9. Deliverables after approval and eligibility

Only after all preconditions pass and a separately approved execution request exists, produce: (a) frozen plan and SHA-256 manifest containing the Git commit hash plus SHA-256 hashes for the classifier, BAM reader, selected scanner, runner, marker TSV, reference FASTA, and (where applicable) panel BED/configuration; (b) source/role decision and provenance manifest; (c) patient-level split and access logs in the authorized location; (d) reproducible analysis outputs and evidence IDs; (e) a results report that states limitations and any deviations. Do not place controlled data or identifying keys in Git.

This draft does not create these execution artifacts, register a plan hash, or authorize data access.

## 10. Owner decisions needed to finalize

### Proposed defaults for owner review (not yet approved)

1. **Reference labels:** require specimen-level PCR/fragment-analysis or another accepted clinical reference assay with documented lab, method, class rule, and provenance for locked validation. MMR-IHC is orthogonal supporting evidence and should be reported separately. Published computational labels may support development/sensitivity analyses only; they do not qualify as validation ground truth under the current source-selection specification.
2. **Source posture:** do not nominate TCGA, CPTAC, Sanger/Broad cell lines, or MSK as eligible now. Their applicable gates remain open or the evidence is insufficient for the validation role. Seraseq-like contrived materials may be considered later for analytical checks only, not patient-level validation.
3. **Partition posture:** require a genuinely independent, patient-level locked cohort for any validation-performance statement. Keep all Stage 6A and tuning-exposed samples in development. Do not split one cohort into development and validation and call the latter independent.
4. **Primary outcomes:** sensitivity and specificity with two-sided 95% score confidence intervals, shown as counts and fractions; report the 2×2 table and insufficient-call rate. Report balanced accuracy only as secondary. FDA's diagnostic-test reporting guidance recommends sensitivity/specificity with two-sided 95% confidence intervals and cautions that estimates depend on the study population and reference standard. This is reporting guidance, not a sample-size or regulatory approval rule. [FDA guidance](https://www.fda.gov/regulatory-information/search-fda-guidance-documents/statistical-guidance-reporting-results-studies-evaluating-diagnostic-tests-guidance-industry-and-fda).
5. **Sample-size rule:** do not set a fixed minimum until the owner states the intended precision or a performance goal. Calculate the required number of *evaluable independent patients separately for MSI-H and MSS* using a predeclared binomial confidence-interval method and a justified expected performance; then inflate for expected unevaluable specimens. MSI-H prevalence varies by tumour type and affects recruitment feasibility: for a consecutive-cohort design, estimate the total recruitment requirement from the expected MSI-H prevalence and the required MSI-H count; a case-control/enriched design may estimate sensitivity and specificity but does not directly estimate prevalence-dependent PPV/NPV. As an illustration only, estimating a proportion near 90% to roughly ±5 percentage points by the simple normal approximation requires about 139 evaluable patients in that class (`1.96² × 0.90 × 0.10 / 0.05²`); this is not a recommended or sufficient clinical-validation sample size, and the final calculation should use the selected score/exact method and intended acceptance criterion. FDA's reporting guidance itself does not prescribe design sample sizes.
6. **Analysis method:** predeclare Wilson score intervals (or exact binomial intervals if the owner prefers a conservative small-sample method); use patient as the independent unit; no threshold tuning or subgroup selection on locked data. Define any acceptance targets only after intended use and clinical risk are specified.
7. **Classifier baseline:** freeze the current development behavior as a versioned experimental baseline, if approved: `MIN_LOCUS_DEPTH=20`, `INSTABILITY_DISTANCE_CUTOFF=0.20`, `MIN_SHIFT_UNITS=1`, `MIN_INFORMATIVE_LOCI=5`, `MSI_H_FRACTION_THRESHOLD=0.40`, read filters `MIN_MAPQ=20` and `MIN_FLANK=5`, smaller-repeat-length mode tie-break, and the conjunctive `shift AND TVD` locus rule. This freezes a reproducible starting point only; it does not calibrate or validate any value. Any later tuning must use development data only and be completed before the locked validation set is opened. Do not change the classifier to `shift OR TVD` based on locked-validation outcomes.

### Decisions still requiring owner confirmation

1. Name and role of the plan approver and independent reviewer.
2. Whether to accept the proposed reference-label hierarchy and source posture.
3. Accept or revise the proposed PCR/IHC and molecular-reference discordance rule; identify the reference-laboratory SOP and independent adjudicator required for unresolved cases.
4. Intended precision or performance target that will drive separate MSS/MSI-H minimums; expected MSI-H prevalence for the intended tumour type(s); partition feasibility and whether a distinct internal-validation cohort is practical.
5. Accept or revise the proposed classifier baseline in item 7 above; separately select the confidence-interval method, multiplicity approach, and any acceptance targets.
6. Confirmation that marker-version/schema selection, synthetic scanner-reader-classifier integration verification, and (when available) samtools-backed synthetic BAM I/O are required gates before labelled-data access.
7. Authorized storage location, access logging, retention/destruction, and execution operator.

**Current disposition: PLAN DRAFTED; APPROVAL AND SOURCE ELIGIBILITY PENDING; STAGE 6B NOT AUTHORIZED TO RUN.**
