# Stage 6 — Source / Label Evidence Review v3: Evidence Gap Classification (development / RUO)

Scope: classification of the `UNVERIFIED` items already documented in `STAGE6_SOURCE_LABEL_EVIDENCE_REVIEW.v2.md`. No new literature
search, no data downloaded/requested/opened, no source or cohort selected or recommended, no development/validation cohort defined,
no thresholds, minimum counts, partitions or Stage 6B Analysis Plan. The scientific design is unchanged. Evidence IDs are those of v2;
items in v2 that had no ID (Section 8) receive a provisional ID here (`CL-01`, `CL-02`). Where a v2 line bundled several questions it is
split into separate items (`UV-nn`) so that each carries exactly one classification. Date: 2026-09-25.

## 1. Classification definitions (as used here)
* **BLOCKER** — must be resolved before the affected source/cohort can be selected: it determines whether that source's labels can be
  assigned a provenance tier, whether it may lawfully/practically be used or reached, or whether its independence from another source can be judged.
* **REQUIRED BEFORE DATA ACCESS** — may stay open while sources are compared, but must be resolved before any data are requested or accessed.
* **NON-BLOCKING** — may remain a documented limitation without preventing the next stage.
* **OUT OF SCOPE** — not needed for the current Stage 6 decision; deferred.

Stages named below: **SSR** = Source Selection Review (next gate); **DAG** = Data Access Gate (pre-access verification, before any request);
**6B-Plan** = Stage 6B Analysis Plan (not started; not defined here).

## 2. Classification table

| Item | Evidence ID | Exact unresolved question | Classification | Rationale (one sentence) | Resolve at |
|---|---|---|---|---|---|
| UV-01 | T-01, T-03 | Do the primary TCGA documents (supplements) confirm the assay panel and the MSI-H/MSI-L/MSS definition per cohort (COAD, READ, UCEC, STAD), currently known only via a secondary paper? | REQUIRED BEFORE DATA ACCESS | The label tier is already assignable from the secondary description, so comparison can proceed, but no sample should be requested on a secondary-only label definition. | DAG |
| UV-02 | T-02 | What are the names of the seven repeat loci used for UCEC MSI testing? | NON-BLOCKING | The MSI-H definition is known at cohort level; locus identity refines interpretation and marker overlap but does not decide source eligibility. | 6B-Plan (only if TCGA UCEC is used) |
| UV-03 | T-04, T-05 | Can each sample-level label for COAD/UCEC be traced to the original assay record (not only the Broad GDAC copy), and what do the STAD "Public/Level 1 Microsatellite Instability Data" files contain? | REQUIRED BEFORE DATA ACCESS | Sample-level traceability is a design principle (every sample with explicit provenance) but is a per-sample property checked when a request list exists, not a cohort-level selection criterion. | DAG |
| UV-04 | T-06 | Is any orthogonal pathology evidence (MMR-IHC or per-marker PCR result) documented for TCGA cases? | NON-BLOCKING | The TCGA label is already an assay-derived label, so orthogonal IHC would corroborate it and help with discordant cases but is not required to assign its tier. | Documented limitation; revisit at 6B-Plan |
| UV-05 | C-04 | Does any orthogonal clinical evidence (PCR/IHC) exist for the CPTAC colon and CPTAC endometrial labels beyond the 85 colon PCR comparisons (colon Table S2; endometrial Table S3)? | BLOCKER | The CPTAC labels are computational (T3), so without orthogonal clinical evidence they cannot qualify as ground truth under the standing ruling, and CPTAC's eligibility cannot be judged. | SSR |
| UV-06 | C-03 | Which PCR assay, loci and laboratory produced the 85 colon PCR results? | REQUIRED BEFORE DATA ACCESS | It matters only if UV-05 shows the PCR subset is usable, and then before those samples are requested. | DAG |
| UV-07 | T-07 | Are the sample lists/barcodes (supplementary tables) available so patient-level leakage, including one patient sequenced at several centres, can be checked? | REQUIRED BEFORE DATA ACCESS | Leakage checks need sample identifiers, which are needed once candidate samples exist, not to compare sources. | DAG |
| UV-08 | S-02 | Why do 257 of 1,087 samples in the five MSI-prone TCGA types have no reference label (830 labelled exomes)? | REQUIRED BEFORE DATA ACCESS | The number and identity of labelled samples must be known before requesting, but the cohort-level label method is unaffected. | DAG |
| UV-09 | T-11 | What label details does Hause 2016 (paywalled) give? | OUT OF SCOPE | Tool-leakage on TCGA labels is already established from two read papers, so this adds nothing to the decision. | Deferred |
| UV-10 | T-10 | What was the specimen preservation (frozen/FFPE) for TCGA COAD and UCEC samples? | REQUIRED BEFORE DATA ACCESS | Pre-analytical comparability must be known before samples are pooled or requested, while STAD (frozen) and CPTAC (frozen/OCT) are verified. | DAG |
| UV-11 | T-10, C-06 | Are per-sample tumour purity/cellularity values available for TCGA (Aran supplement, ABSOLUTE for 11 types) and CPTAC (ESTIMATE colon; ABSOLUTE and methylation-based UCEC)? | REQUIRED BEFORE DATA ACCESS | Purity affects interpretation of any empirical result, so its per-sample availability must be confirmed before samples are chosen. | DAG |
| UV-12 | C-02 | Why does the CPTAC endometrial count differ between the primary paper (25 MSI) and a 2026 secondary paper (21 MSI-H, 73 MSS, 1 missing)? | REQUIRED BEFORE DATA ACCESS | The CPTAC label is T3 either way, so the discrepancy does not change eligibility, but the authoritative per-sample label set must be fixed before requesting. | DAG |
| UV-13 | C-02 | Does the second CPTAC endometrial set (N=108, secondary source only) exist, and how were its labels made? | NON-BLOCKING | It is not a documented candidate in v2 and would need its own evidence if it were ever added. | SSR (only if proposed as a source) |
| UV-14 | C-07 | What is the access tier of the colon raw genomics at SRA BioProject PRJNA514017? | REQUIRED BEFORE DATA ACCESS | The GDC/dbGaP controlled route is verified, and this SRA tier only affects which route to request. | DAG |
| UV-15 | C-07 | Do CPTAC-specific data-use limitations permit the intended development use? | BLOCKER | A source whose terms exclude the intended use cannot be selected. | SSR |
| UV-16 | T-09 | Which dbGaP consent groups and data-use limitations apply to TCGA (phs000178) and do they permit the intended development use? | BLOCKER | The tier and process are verified, but permitted use is unverified and decides whether TCGA can be selected. | SSR |
| UV-17 | T-09, C-07 | What are the operational obligations of controlled-data users (security, retention, sharing, publication, term)? | REQUIRED BEFORE DATA ACCESS | These govern how access is requested and handled and do not affect which source is preferable. | DAG |
| UV-18 | X-02 | Is there patient-level overlap between TCGA and the CPTAC prospective colon and endometrial cohorts? | BLOCKER | Independence between candidate sources cannot be judged without it, and a source cannot be assigned an independent role while overlap is unknown. | SSR |
| UV-19 | X-03 | Is there overlap between MSK-IMPACT and TCGA, CPTAC or other cohorts? | BLOCKER | The same independence requirement applies, and no published statement was found. | SSR |
| UV-20 | X-01 | What does the full text of Zhang 2014 say about the 95 TCGA tumours (90 patients) in the retrospective CPTAC colorectal study? | OUT OF SCOPE | That study is not a candidate source, and the live overlap question is carried by UV-18. | Deferred |
| UV-21 | M-04 | Is there any sequence-level access route for MSK-IMPACT data (the published statements cover processed results and clinical data only)? | BLOCKER | A source with no route to its raw data cannot be a candidate. | SSR |
| UV-22 | M-04 | How does the 12,288-patient MSI study map to the 10,000-patient public release, and are MSIsensor scores in that release? | REQUIRED BEFORE DATA ACCESS | It matters only if UV-21 yields an access route, and then before anything is requested. | DAG (conditional on UV-21) |
| UV-23 | CL-01 (v2 §8) | Which cell line(s) underlie the Seraseq "cell line-based" MSI-High material? | NON-BLOCKING | Seraseq is contrived analytical material and not a cohort, so the question only affects overlap with cell-line sources. | SSR (only if Seraseq is compared) |
| UV-24 | CL-02 (v2 §8) | What are the data-use terms of the Sanger Cell Model Passports/DepMap and CCLE resources? | BLOCKER | Permitted use decides whether the cell-line resources can be selected, as for UV-15 and UV-16. | SSR |

## 3. Counts
| Classification | Items | Count |
|---|---|---|
| BLOCKER | UV-05, UV-15, UV-16, UV-18, UV-19, UV-21, UV-24 | **7** |
| REQUIRED BEFORE DATA ACCESS | UV-01, UV-03, UV-06, UV-07, UV-08, UV-10, UV-11, UV-12, UV-14, UV-17, UV-22 | **11** |
| NON-BLOCKING | UV-02, UV-04, UV-13, UV-23 | **4** |
| OUT OF SCOPE | UV-09, UV-20 | **2** |
| Total | | **24** |

Every one of the 24 items appears exactly once. The 12 numbered "remaining UNVERIFIED" entries in v2 Section 9 map onto these items
(v2 no.1→UV-01, UV-02; no.2→UV-03; no.3→UV-04, UV-05; no.4→UV-07, UV-08; no.5→UV-09; no.6→UV-10, UV-11; no.7→UV-12, UV-13; no.8→UV-06;
no.9→UV-14, UV-15, UV-16, UV-17; no.10→UV-18, UV-19, UV-20; no.11→UV-21, UV-22; no.12→UV-23, UV-24).

## 4. Answers to the focus questions (classification only)
* **Is label provenance sufficient?** For comparing sources: yes at cohort level — TCGA labels are assay-derived (partially verified), CPTAC labels are
  computational (verified), MSK labels are MSIsensor calls checked against PCR/IHC in a subset, cell-line labels are computational. For
  using any sample: no — primary-source confirmation and per-sample traceability remain open (UV-01, UV-03), so these are gated at DAG.
* **Is orthogonal clinical evidence necessary?** It is a BLOCKER only where the label is otherwise purely computational (CPTAC, UV-05); where the label is
  already assay-derived (TCGA, UV-04) its absence is a documented limitation.
* **Is sample-level traceability required?** Yes, but as a pre-access condition (UV-03, UV-07, UV-08), not as a source-selection criterion.
* **Patient/cohort overlap and leakage:** unresolved overlap is a BLOCKER because independence cannot be judged without it (UV-18, UV-19); the
  within-TCGA duplicate-sequencing leakage check waits for sample identifiers (UV-07).
* **Specimen/purity metadata:** required before access (UV-10, UV-11); not a selection criterion.
* **Controlled-access and data-use obligations:** what use is permitted is a BLOCKER (UV-15, UV-16, UV-24); how access is operated is pre-access (UV-17).
* **Count discrepancies:** required before access (UV-08, UV-12); a hypothetical extra cohort is non-blocking (UV-13).
* **Raw-data accessibility:** a BLOCKER only for MSK-IMPACT, where no route is documented (UV-21); TCGA and CPTAC access routes are verified.

## 5. Resolution-route caveat (noted, not acted on)
Several pre-access items (UV-01, UV-03, UV-05, UV-06, UV-07, UV-08, UV-11, UV-12) can probably be resolved only by reading per-sample
supplementary tables or sample lists. Those tables contain patient-derived identifiers and, in some cases, MSI labels. Opening them
would conflict with the current instruction not to open MSI-H/MSS or patient data, so I have not opened any; an explicit owner ruling on whether
such documentation-level tables may be read is needed before those items can be worked.

## Next Gate: Source Selection Review
Purpose: use the resolved blockers to decide, source by source, whether each candidate is eligible to be compared (label tier, permitted use,
access route, independence). It needs the seven BLOCKER items (UV-05, UV-15, UV-16, UV-18, UV-19, UV-21, UV-24) resolved or explicitly recorded as
unresolved. It does not select cohorts, define partitions or thresholds, or start the Stage 6B Analysis Plan. The 11 pre-access items wait for the
Data Access Gate. Nothing proceeds without the owner's instruction.

Decision Status: NO SOURCE SELECTED
