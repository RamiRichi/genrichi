# Stage 6 — TCGA ↔ CPTAC Independence Resolution (UV-18) (documentation only)

Question: **Can we establish, with sufficient evidence, whether the relevant candidate CPTAC cohorts contain patients/samples that overlap TCGA?**
Only public documentation and open metadata were read. No controlled-access data (dbGaP, EGA, controlled SRA/GDC files) was accessed; no patient-level linkage was built; no cohort or source selected or preferred; no provider contacted;
no code, threshold, classifier or production file touched. Aggregate counts of open GDC case metadata and identifier *shapes* were inspected; no identifier is reproduced. The TCGA, CPTAC and cell-line governance documents were not reopened or modified.
Independence is not inferred from different accessions, papers or institutions; overlap is not inferred from shared TCGA-derived comparison. Date: 2026-09-25.

## 1. Scope
The overlap question is evaluated for: (a) **CPTAC-2 Retrospective** versus TCGA colorectal; (b) the **CPTAC colon** cohort of Vasaikar 2019 (candidate for Stage 6); (c) the **CPTAC endometrial** cohort of Dou 2020 (candidate for Stage 6); each versus TCGA COAD/READ/UCEC.
One further CPTAC cohort surfaced by open GDC metadata (gastric, CPTAC-3) is noted in §4 but not reviewed. The question is *not* which source is preferable.

## 2. Evidence table
Evidence-quality labels: **PRIMARY** (paper or official database text read directly), **OPEN METADATA** (public GDC/dbGaP metadata), **SECONDARY** (search-result text, descriptive pages). Status words: ESTABLISHED / NOT ESTABLISHED / PARTIAL / NOT APPLICABLE.

| Cohort/material | TCGA relationship | Evidence | Evidence quality | Independence status |
|---|---|---|---|---|
| **CPTAC-2 Retrospective** | **ESTABLISHED overlap**: proteomics on TCGA colorectal specimens. Sci Data descriptor: "We analyzed 95 samples representing 90 TCGA colon and rectal (CRC) tumors, which comprise a subset of the 224 sample collection subjected to multiplatform genomic analyses by the TCGA"; "Of the 276 samples described in the TCGA study, 95 specimens from 90 patients were available … (5 patients' tumors were provided as duplicate pieces)"; "Specimens for proteomic study were sectioned sequentially from the tumor blocks used for genomic studies". Zhang 2014: "the Clinical Proteomic Tumor Analysis Consortium (CPTAC) is performing proteomic analyses of TCGA tumor specimens". CGC: this CPTAC-2 Retrospective set "consists of four TCGA cancer types (Ovarian …, Breast …, Colon adenocarcinoma, Rectum adenocarcinoma)". | PRIMARY (Zhang 2014, PMC4249766; Scientific Data 2015, PMC4477697); SECONDARY (CGC page) | **NOT APPLICABLE — overlap ESTABLISHED** (same tumors) |
| **CPTAC colon** (prospective; Vasaikar 2019; 110 patients) | **NOT ESTABLISHED** either way. Provenance documented: "collected by several tissue source sites in strict accordance to the CPTAC-2 colon procurement protocol … with an informed consent from the patients"; inclusion "newly diagnosed, untreated patients undergoing primary surgery for colon adenocarcinoma"; "prospectively collected". dbGaP phs000892 describes confirmatory samples "designed to confirm CPTAC findings from the TCGA samples … collected via a protocol optimized for proteomics" (ischaemic time under 30 minutes). No source states that these patients are, or are not, TCGA participants. TCGA itself also describes a prospective collection ("most of the tumours … were derived from a prospective collection", COAD 2012), so "prospective" does not discriminate. | PRIMARY (Vasaikar 2019; dbGaP study page phs000892; TCGA COAD 2012) | **NOT ESTABLISHED** (PARTIAL supporting provenance) |
| **CPTAC endometrial** (prospective; Dou 2020; 95 tumours) | **NOT ESTABLISHED** either way. Provenance documented: "Biospecimens were collected from newly diagnosed patients with endometrial cancer … undergoing surgical resection … no prior treatment"; tissue source sites in the US, Europe and Asia; cold ischaemic time within 30 min. dbGaP phs001287 (CPTAC-3, launched September 2016): inclusion "Newly diagnosed, untreated patients undergoing primary cytoreductive surgery … Informed consent provided for the provision of biospecimens and data". No source states whether any patient is a TCGA participant. | PRIMARY (Dou 2020; dbGaP study page phs001287) | **NOT ESTABLISHED** (PARTIAL supporting provenance) |

## 3. Known overlap (documented; not extrapolated)
* **What is established:** the CPTAC-2 Retrospective colorectal proteomics study used **95 specimens from 90 TCGA patients** (64 COAD + 31 READ, from search-result text of the paper's record, so SECONDARY; 5 patients contributed duplicate pieces), cut sequentially from the same tumor blocks TCGA used for genomics; they are a subset of TCGA's 224 (Sci Data) / 276 (TCGA text) sample collections. Their genomic data are TCGA's own.
* **Which CPTAC material this is:** a *retrospective* proteomic study of TCGA specimens. It is **not** the GDC project CPTAC-2 (342 cases with "CPTAC" IDs, none of whose submitter IDs contains "TCGA": open GDC metadata) and **not** the prospective colon cohort of Vasaikar 2019.
* **Relevance to the candidate prospective cohorts:** irrelevant as evidence of overlap for them, relevant as a **warning**: any data labelled "CPTAC colorectal" must be traced to the retrospective set (TCGA patients) or the prospective set before use.
* No statement extends this overlap to any other cohort.

## 4. Independence determination (answered separately)
**CPTAC colon — independence NOT ESTABLISHED.**
* Provenance is documented (new procurement protocol, informed consent, newly diagnosed patients, tissue source sites) but nothing states that no patient was also a TCGA participant, and the CPTAC-2 phase (launched September 2011 per dbGaP) is not shown to be temporally disjoint from TCGA colorectal collection (collection dates for TCGA were not established here), so timing cannot be used.
* **Cohort-membership discrepancies preserved, not resolved:** the paper reports 110 patients and 106 exomes (24 MSI-H, 82 MSS); the public cBioPortal table lists 110 patients of whom 106 are flagged sequenced; open GDC metadata for project CPTAC-2 lists 105 colon-site cases plus 1 rectum-site case (106 in total), among 342 cases of five tumour sites, with dbGaP phs000892 covering 348 consented subjects across four tumour types. The numerical coincidence with 106 is **not** treated as a mapping. These discrepancies do not favour overlap or independence; they mean the exact patient set of "the colon cohort" cannot yet be defined, which any patient-level check would require.

**CPTAC endometrial — independence NOT ESTABLISHED.**
* Provenance documented (newly diagnosed, consented for CPTAC biospecimens; phase started September 2016, after the 2013 TCGA endometrial publication), which is supporting circumstance but not a statement of non-overlap; "prospectively collected" is not discriminating (see above).
* **Discrepancy preserved:** the paper's cohort is 95 tumours, while open GDC metadata lists **241 "Uterus, NOS" cases** in CPTAC-3 — the relationship between the 95 and the 241 is not documented in what I read.

**Other CPTAC cohort identified by authoritative evidence:** open GDC metadata lists **165 stomach cases** in CPTAC-3 (relevant to TCGA-STAD) and 1 colon case in CPTAC-3. Their provenance was **not reviewed**; independence from TCGA is **NOT ESTABLISHED** by default, not by finding.

**Overall for the prospective cohorts:** `NO OVERLAP DOCUMENTED — INDEPENDENCE NOT PROVEN` (the classification used in the gate closure), meaning neither overlap nor independence is established.

## 5. Final gate
TCGA-CPTAC INDEPENDENCE STATUS: BLOCKED

Evidence missing (nothing else is claimed):
1. An **authoritative provenance statement** from NCI/CPTAC/the Biospecimen Core Resource stating, per cohort and accession, whether any CPTAC patient was also a TCGA participant.
2. Or an **authorized patient-level identity cross-check** (see §7) — not performed and not permitted at this stage.
3. A defined **cohort ↔ accession mapping** for the colon cohort (110/106/105+1; phs000892 vs phs001287 vs SRA PRJNA514017, from the CPTAC governance document) and for the endometrial cohort (95 vs 241 GDC cases), so that "the cohort" is a defined set of patients.

## 6. Consequence for Stage 6 (methodological only)
* **No independent-validation claim** may currently rest on TCGA and the prospective CPTAC cohorts being independent.
* The prospective CPTAC colon and endometrial cohorts can only be treated as **potentially overlapping with TCGA**; any CPTAC-2 Retrospective colorectal material must be treated as **the same patients** as TCGA COAD/READ.
* An independent cohort or a provenance confirmation is still required before any independence claim. No final validation source is chosen here.

## 7. What evidence would suffice (Objective D; none of it performed)
| Route | Sufficient if… | Status now |
|---|---|---|
| Authoritative provenance statement | NCI/CPTAC/BCR states explicitly, for a named cohort and accession, that none (or which) of its patients are TCGA participants | Not found in public documentation |
| Authorized patient-level linkage | Performed by an Approved User under DAC-approved access, using subject-level identity evidence (documented subject identifiers or germline genotype concordance between the two programs' normal samples), with the method and result recorded | Requires controlled access; not performed |
| Explicit provider statement | Written response from NCI/CPTAC to a specific enquiry | Not requested; no provider contacted |
| Public unique-sample crosswalk | A legally accessible table linking CPTAC case IDs to TCGA barcodes, if one exists | None found |
Sufficiency is judged against the intended role of each source, not in the abstract.

## 8. Next action (not performed)
The minimum action is a **written provenance confirmation from NCI/CPTAC (with the Biospecimen Core Resource)** stating, for the CPTAC colon cohort (with its defining accession) and the CPTAC endometrial cohort (with its defining accession and the 95-versus-241 relationship), whether any patient is also a TCGA participant. No provider is contacted by this document.
