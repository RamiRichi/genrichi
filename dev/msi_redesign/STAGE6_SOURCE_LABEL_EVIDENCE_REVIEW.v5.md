# Stage 6 — Source / Label Evidence Review v5: Final Gap Disposition (development / RUO)

Scope: disposition of every remaining `UNVERIFIED` or `PARTIALLY VERIFIED` item, using **only** the evidence already documented in
`STAGE6_SOURCE_LABEL_EVIDENCE_REVIEW.v4.md`. No further searching, downloading, data requests or attempts to bypass 403/anti-automation
restrictions were made; no controlled, raw-sequence, patient or authorisation-only data were touched. No source or sample is selected, no cohort constructed,
no threshold, minimum count, partition or Stage 6B Analysis Plan defined. No other file was modified. Date: 2026-09-25.

## 1. Rules applied
* Statuses are carried over from v4 unchanged. **No item was upgraded** to VERIFIED; small remaining uncertainty is not treated as resolved.
* Nothing is inferred: independence between cohorts, permission for commercial or diagnostic development under any data-use agreement, raw-data accessibility,
  and clinical validity from computational concordance are all **not** inferred anywhere below.
* Exactly one disposition per item; no additional category:
  1. **MUST RESOLVE BEFORE SOURCE SELECTION** — the point decides whether a source is eligible at all.
  2. **MUST RESOLVE BEFORE DATA ACCESS** — eligibility can be assessed provisionally, but the point must be closed before anything is requested or accessed.
  3. **DOCUMENTED LIMITATION — ACCEPTABLE** — the limitation does not prevent consideration, provided it is explicitly carried forward.
  4. **NOT USABLE FOR THIS STAGE** — the evidence limitation makes the source or evidence unsuitable for the specific role named, with no statement about its general scientific value.
* Item numbering (UV-01…UV-24) and evidence IDs are those of v3/v4.

## 2. Final dispositions

| Item | Evidence ID | v4 status | Final disposition | One-sentence justification | Exact consequence for the next gate |
|---|---|---|---|---|---|
| UV-01 | T-01, T-03 | PARTIALLY VERIFIED | **2. MUST RESOLVE BEFORE DATA ACCESS** | The assay type is documented for STAD and the marker panels/≥40% rule for colorectal and endometrial only through a secondary paper, which is enough to compare sources but not to request samples. | Source Selection Review may proceed treating the TCGA MSI-H definition as secondary-sourced; no TCGA sample may be requested until the primary assay description and threshold per cohort are documented. |
| UV-02 | T-02 | UNVERIFIED | **3. DOCUMENTED LIMITATION — ACCEPTABLE** | The names of the seven UCEC loci refine interpretation but do not decide eligibility, since the cohort-level definition is known. | Carry forward "UCEC marker panel not identified"; revisit only if TCGA UCEC labels are ever used. |
| UV-03 | T-04, T-05 | PARTIALLY VERIFIED | **2. MUST RESOLVE BEFORE DATA ACCESS** | Traceability to the assay files is documented structurally for STAD only, and sample-level provenance is a pre-access condition. | Before any request, sample-level label-to-assay traceability must be documented for the cohort concerned; COAD/UCEC currently have none. |
| UV-04 | T-06 | PARTIALLY VERIFIED | **3. DOCUMENTED LIMITATION — ACCEPTABLE** | Site-recorded MSI/MMR-IHC fields exist for a minority of colorectal cases and none for endometrial or gastric, and their agreement with the consortium label was never tested. | Carry forward "orthogonal pathology evidence for TCGA is incomplete and untested against the consortium label"; no clinical-validity claim may rest on it. |
| UV-05 | C-04 | PARTIALLY VERIFIED | **4. NOT USABLE FOR THIS STAGE** | CPTAC labels are computational, the orthogonal evidence covers only 36 of 95 endometrial tumours (17 of the 25 MSI-H not interpretable) and a colon PCR subset whose per-sample table was not read, and computational concordance is not clinical validity. | For the role "ground-truth MSI labels", CPTAC colon and endometrial labels are not usable at this stage; this says nothing about CPTAC's other scientific value, and the item reopens only if the owner supplies the gated per-sample PCR/IHC tables. |
| UV-06 | C-03 | UNVERIFIED | **3. DOCUMENTED LIMITATION — ACCEPTABLE** | The colon PCR assay details matter only if CPTAC colon labels were reconsidered, which UV-05 already excludes for the ground-truth role. | Carry forward "PCR assay, loci and laboratory behind 85 colon comparisons unknown". |
| UV-07 | T-07 | PARTIALLY VERIFIED | **2. MUST RESOLVE BEFORE DATA ACCESS** | Patient linkage within TCGA is visible (multiple samples per patient) but cross-cohort linkage and supplementary sample lists are not, and leakage checks need identifiers. | Before any request, a patient-level leakage check design using documented identifiers must exist; one patient sequenced as several exome pairs must be treated as one unit. |
| UV-08 | S-02 | UNVERIFIED | **2. MUST RESOLVE BEFORE DATA ACCESS** | Which of the 1,087 samples carry a consortium reference label is unknown and must be known before choosing anything to request. | Before any request, the labelled/unlabelled status per sample must be documented (830 labelled exomes vs 1,087 samples is currently unexplained). |
| UV-09 | T-11 | UNVERIFIED | **3. DOCUMENTED LIMITATION — ACCEPTABLE** | Hause 2016 label details are unread but tool-leakage on TCGA labels is already documented from two read papers. | Carry forward "Hause 2016 not read". |
| UV-10 | T-10 | PARTIALLY VERIFIED | **2. MUST RESOLVE BEFORE DATA ACCESS** | FFPE-flagged tumours exist in TCGA colorectal (13 of 476) and endometrial (4 of 552) data and preservation is not otherwise documented, so pre-analytical status must be settled per sample before requesting. | Before any request, per-sample specimen preservation must be documented; "not flagged FFPE" is not to be read as "frozen". |
| UV-11 | T-10, C-06 | PARTIALLY VERIFIED | **2. MUST RESOLVE BEFORE DATA ACCESS** | Pathology tumour-nuclei percentages are public but sometimes below the stated policy thresholds and genomic purity estimates were not obtainable. | Before any request, per-sample purity/cellularity availability (and its source) must be documented; policy thresholds are not assumed to hold per sample. |
| UV-12 | C-02 | PARTIALLY VERIFIED | **3. DOCUMENTED LIMITATION — ACCEPTABLE** | The authoritative endometrial count is established and the remaining differences do not decide eligibility. | Carry forward, undiscarded: CPTAC endometrial **25 MSI-H of 95** (primary paper and public table agree) versus the 2026 secondary source's **21 MSI-H / 73 MSS / 1 missing**, cause unknown; CPTAC colon primary text 24 MSI-H / 82 MSS (106 exomes) versus public table 24 / 81 with 5 not available (110 patients), cause unknown. |
| UV-13 | C-02 | UNVERIFIED | **4. NOT USABLE FOR THIS STAGE** | The second endometrial set (N=108) rests on a single secondary source with no primary paper, so it is unusable as evidence. | The 108-tumour set is not to be cited or counted as a source; it may return only if a primary publication is supplied. |
| UV-14 | C-07 | PARTIALLY VERIFIED | **2. MUST RESOLVE BEFORE DATA ACCESS** | The BioProject page states no public data is linked and no tier, which is not proof of a route. | Before any request, the actual access route and tier for CPTAC colon raw genomics must be documented; accessibility is not inferred. |
| UV-15 | C-07 | UNVERIFIED | **1. MUST RESOLVE BEFORE SOURCE SELECTION** | CPTAC-specific data-use limitations are unread, so it is unknown whether any use consistent with the project is permitted. | CPTAC cannot be assessed for eligibility at Source Selection Review until its data-use terms are documented; nothing about commercial or diagnostic use is assumed. |
| UV-16 | T-09 | PARTIALLY VERIFIED | **1. MUST RESOLVE BEFORE SOURCE SELECTION** | The 2014 TCGA DUC ties use to an approved research project, is silent on commercial use, and the consent-group limitations and current version were unreadable. | TCGA cannot be assessed for eligibility until the current use limitations are documented and it is established what the approved-use route covers; permission is not inferred from silence. |
| UV-17 | T-09, C-07 | PARTIALLY VERIFIED | **2. MUST RESOLVE BEFORE DATA ACCESS** | Operational obligations are documented for TCGA (2014 DUC) but not for CPTAC and concern how access is handled, not which source is eligible. | Before any request, the current obligations (reporting, security, destruction, publication) for the source concerned must be documented. |
| UV-18 | X-02 | UNVERIFIED | **1. MUST RESOLVE BEFORE SOURCE SELECTION** | No published statement and no usable identifier comparison establishes whether TCGA and the CPTAC prospective cohorts share patients, and independence is not to be inferred. | No statement of independence between TCGA and CPTAC may be made, and they cannot be assessed together in an independent-validation role, until overlap is documented. |
| UV-19 | X-03 | UNVERIFIED | **3. DOCUMENTED LIMITATION — ACCEPTABLE** | MSK-IMPACT is already not usable as a data source (UV-21), so its unresolved overlap with other cohorts matters only if that changes. | Carry forward "MSK-IMPACT overlap with any cohort unknown; independence not established"; it reopens as a pre-selection item only if a documented access route later reinstates MSK as a candidate. |
| UV-20 | X-01 | UNVERIFIED | **3. DOCUMENTED LIMITATION — ACCEPTABLE** | The retrospective CPTAC–TCGA colorectal overlap rests on search-result text only and that study is not a candidate. | Carry forward "documented TCGA–CPTAC overlap exists for the retrospective study (95 TCGA tumours); full text not read" as a warning against assuming independence. |
| UV-21 | M-04 | UNVERIFIED | **4. NOT USABLE FOR THIS STAGE** | No sequence-level access route for MSK-IMPACT is documented, raw-data accessibility is not inferred, and only processed data are public. | For the role "data source for this project", MSK-IMPACT is not usable at this stage; its published clinical evidence (12,288 patients; MSIsensor ≥10; PCR/IHC comparison) remains citable as literature only, and the item reopens only if a documented route is supplied. |
| UV-22 | M-04 | PARTIALLY VERIFIED | **3. DOCUMENTED LIMITATION — ACCEPTABLE** | The public MSK release has no MSI field and its relation to the 12,288-patient cohort is undocumented, which is moot while MSK is not a usable source. | Carry forward "public MSK release: 10,945 samples/10,336 patients, no MSI score, about 82% FFPE; relation to the Middha cohort unknown". |
| UV-23 | CL-01 | UNVERIFIED | **3. DOCUMENTED LIMITATION — ACCEPTABLE** | The cell line(s) behind the Seraseq material are unknown, which limits only overlap assessment with cell-line sources, not consideration of the material. | Carry forward "Seraseq underlying cell line unknown; overlap with cell-line resources not assessable". |
| UV-24 | CL-02 | PARTIALLY VERIFIED | **1. MUST RESOLVE BEFORE SOURCE SELECTION** | The Sanger policy allows internal proprietary research but states that commercial use of the data is not permitted without prior consent, and the Broad DepMap/CCLE terms were unreadable. | Cell-line resources cannot be assessed for eligibility until it is documented whether the project's use is permitted (for Sanger, whether prior consent exists; for Broad/CCLE, what the terms say); permission is not assumed. |

## 3. Final disposition counts
| Disposition | Items | Count |
|---|---|---|
| 1. MUST RESOLVE BEFORE SOURCE SELECTION | UV-15, UV-16, UV-18, UV-24 | **4** |
| 2. MUST RESOLVE BEFORE DATA ACCESS | UV-01, UV-03, UV-07, UV-08, UV-10, UV-11, UV-14, UV-17 | **8** |
| 3. DOCUMENTED LIMITATION — ACCEPTABLE | UV-02, UV-04, UV-06, UV-09, UV-12, UV-19, UV-20, UV-22, UV-23 | **9** |
| 4. NOT USABLE FOR THIS STAGE | UV-05, UV-13, UV-21 | **3** |
| Total | | **24** |

## 4. The seven v3 blockers
| Item | Disposition | Effect |
|---|---|---|
| UV-05 CPTAC orthogonal evidence | 4. NOT USABLE FOR THIS STAGE | CPTAC labels are not usable as ground truth at this stage. |
| UV-15 CPTAC data-use terms | 1. MUST RESOLVE BEFORE SOURCE SELECTION | Open; needs CPTAC terms documented. |
| UV-16 TCGA data-use terms | 1. MUST RESOLVE BEFORE SOURCE SELECTION | Open; needs the current use limitations and what the approved-use route covers. |
| UV-18 TCGA–CPTAC overlap | 1. MUST RESOLVE BEFORE SOURCE SELECTION | Open; no independence claim allowed. |
| UV-19 MSK overlap | 3. DOCUMENTED LIMITATION — ACCEPTABLE | Moot while MSK is not usable; reopens only if a route is documented. |
| UV-21 MSK sequence-level access | 4. NOT USABLE FOR THIS STAGE | MSK-IMPACT is not usable as a data source at this stage. |
| UV-24 cell-line data-use terms | 1. MUST RESOLVE BEFORE SOURCE SELECTION | Open; Sanger consent and Broad/CCLE terms needed. |

## 5. Carry-forward register (limitations that must appear in any later document)
Orthogonal evidence incomplete (TCGA, CPTAC); CPTAC labels computational; CPTAC endometrial count 25/95 versus secondary 21 (kept, unexplained); CPTAC colon 24/82 versus 24/81 + 5 unavailable;
UCEC loci unnamed; Hause 2016 unread; colon PCR details unknown; TCGA–CPTAC retrospective overlap documented, prospective overlap unknown; MSK overlap unknown and MSK raw data not accessible from documentation;
MSK public release has no MSI field and is mostly FFPE; Seraseq underlying line unknown. None of these is a claim about clinical validity or about any source's general scientific value.

## 6. What the four "before source selection" items need (no action taken)
They can only be closed by documents that I could not read without further searching or bypassing restrictions: the current dbGaP data-use text and consent-group limitations for phs000178 and phs001287
(UV-15, UV-16), a published statement or identifier-level documentation of TCGA–CPTAC overlap (UV-18), and Sanger consent status plus the Broad DepMap/CCLE terms (UV-24). Obtaining them, and deciding
how, is the owner's call.

Decision Status: NO SOURCE SELECTED

Next Gate: Source Selection Review — only after all items classified as MUST RESOLVE BEFORE SOURCE SELECTION are resolved.
