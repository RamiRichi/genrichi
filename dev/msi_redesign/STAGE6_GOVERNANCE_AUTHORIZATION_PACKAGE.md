# Stage 6 — Governance Authorization Package (preparation only; development / RUO)

Nothing in this document requests access, submits an application, contacts any provider or authority, or opens controlled or patient data. No source or cohort is selected; no partition,
threshold, minimum sample count or Stage 6B Analysis Plan is defined. It restates only findings already established in `STAGE6_SOURCE_LABEL_EVIDENCE_REVIEW.v5.md` and
`STAGE6_PRESELECTION_GATE_CLOSURE.md` (no further searching). It does not interpret any agreement beyond its text and is **not legal advice**. Date: 2026-09-25.
The words "approved", "permitted" and "validated" are used only where the cited documentation says so.

## 1. TCGA / dbGaP
| Item | Established (source: dbGaP Data Use Certification Agreement Ver. 25 Jan 2025 for phs000178.v11.p8 and its Addendum; GDC pages) |
|---|---|
| Accession / study | dbGaP **phs000178** — The Cancer Genome Atlas (TCGA); data access committee named in the addendum: **NCI DAC** |
| Current access tier | **Controlled access** for primary sequence data (BAM) and similar files, via dbGaP Authorized Access with eRA Commons authentication; clinical data and somatic mutation calls are open-tier (dbGaP study page, GDC). Consent group 1 "General Research Use (GRU)". |
| Data-use restrictions (as written) | "Use of the data is limited only by the terms of the model Data Use Certification." Research use "solely in connection with the approved research project described in the DAR" (Research Use Statement); new uses require a new DAR; access granted for one year, renewable, with annual Progress Update; collaborators at other institutions submit their own DAR; datasets and Data Derivatives may not be redistributed without NIH approval and "may not be sold to any individual at any point in time for any purpose"; publication encouraged with acknowledgement and accession number; "Public Posting of Genomic Summary Results — Not Allowed"; security best practices; destruction at close-out. Requester = "home institution or organization" through an Institutional Signing Official; PI = permanent employee at a level equivalent to a tenure-track professor or senior scientist. |
| Commercial-product language (verified, quoted) | NIH "considers these data as pre-competitive and urges Approved Users to avoid making IP claims derived directly from the dataset(s)"; "these NIH-provided data, and conclusions derived therefrom, will remain freely available, without requirement for licensing"; "there is no restriction on development of commercial products resulting from the knowledge gained from the research project". |
| Unresolved question | Whether GenRichi's intended use — using TCGA controlled data to develop and validate a proprietary MSI panel for a diagnostic product, including any regulatory use — is an "approved research project" and is compatible with the IP language above. The agreement contains no language on diagnostic or clinical development, IVD use, regulatory submission, or for-profit requester eligibility. Silence is not treated as permission. |
| Determination required | **(a) NCI DAC:** a decision on a Data Access Request whose Research Use Statement describes the intended use (the DAC, not GenRichi, determines whether that use is approved). **(b) Legal/compliance:** a reading of DUC terms 1 (Research Use), 6 (no sale/distribution), and 9 (Intellectual Property) against the intended product and any regulatory submission. **(c) Legal/compliance:** confirmation that the company can act as Requester (Signing Official, eRA Commons) and that an eligible PI exists. |

## 2. CPTAC / dbGaP
| Item | Established (sources: DUC Ver. 25 Jan 2025 for phs001287.v23.p7 with Addendum; GDC, CGC pages; Dou 2020; Vasaikar 2019) |
|---|---|
| Accession / study | dbGaP **phs001287** — CPTAC Proteogenomic Study (cited by the endometrial paper for raw genomic data); **NCI DAC**. The colon paper cites **SRA BioProject PRJNA514017** for raw genomic data, which states no public data linked. |
| Access route as established | dbGaP controlled access and GDC (controlled BAM/VCF; some open-tier data such as CPTAC-3 RNA-seq quantification). Consent group 1 "General Research Use (GRU)". Same DUC terms as §1 (research use within the approved DAR, one-year renewable access, no sale of Data Derivatives, security, close-out). Public posting of genomic summary results "Allowed"; publication rules: freedom to publish if a global-analysis publication exists, or 22 months after all samples were received at the data production site, or CPTAC Steering Committee approval; required acknowledgement and citation of phs001287 with retrieval date. |
| Unresolved CPTAC-specific terms | The addendum's "IC Specific Access Term" requires users to recognise restrictions "stated on the CPTAC Data Use Agreement page"; that page now redirects to an unrelated NCI overview, so those restrictions are **unread**. The retired CPTAC portal FAQ states a 15-month embargo variant and no commercial or diagnostic-development terms. The dataset-specific publication status (study page "Available Data" tab) was not read. No diagnostic-development language found. |
| Unresolved raw-data access route | For **CPTAC colon** raw sequence, the applicable repository, accession and access tier are not established (SRA BioProject PRJNA514017 vs dbGaP phs001287). For **CPTAC endometrial**, the paper cites GDC or dbGaP phs001287 on request. |
| Determination required | **(a) Data provider (NCI/CPTAC):** the current CPTAC Data Use Agreement text and the applicable access route/tier for the colon raw sequence data. **(b) NCI DAC:** a decision on a DAR describing the intended use. **(c) Legal/compliance:** review of the DUC and the CPTAC-specific term against the intended product. |

Note (from v5, not a new finding): CPTAC labels are computational and were classified NOT USABLE FOR THIS STAGE as ground truth (UV-05); this package does not change that.

## 3. TCGA–CPTAC independence
| Item | Established |
|---|---|
| Documented overlap | **CPTAC-2 Retrospective consists of TCGA tumours**: the CGC page describes its MS data as "four TCGA cancer types (Ovarian …, Breast …, Colon adenocarcinoma, Rectum adenocarcinoma)"; Zhang et al. 2014 analysed 95 TCGA tumour samples from 90 patients (64 COAD, 31 READ) (search-result text of the paper and data record; the full text was not read). |
| Retrospective vs prospective | The CPTAC colon cohort of Vasaikar 2019 and the endometrial cohort of Dou 2020 are described as **prospective** cohorts (110 colon patients; 95 endometrial patients) and are compared with TCGA cohorts as external comparators. The TCGA-derived CPTAC-2 Retrospective material is a different data set from these prospective cohorts. |
| Why independence of the prospective cohorts is unproven | Neither paper, the GDC page nor the CGC page states that the prospective patients are absent from TCGA; public metadata tables use different identifier formats, which neither shows nor excludes overlap; a secondary 2026 paper's statement that identifiers were cross-checked is unpublished and secondary. Absence of an overlap statement is not evidence of independence. |
| Evidence that would be sufficient (no analysis performed here) | Either (i) a written provenance statement from the data-holding authority (NCI / CPTAC) stating whether any CPTAC prospective-cohort patient is a TCGA participant, or (ii) a patient-level identity cross-check (subject identifiers or germline genotype concordance) carried out under authorised access by an Approved User. Sufficiency must be judged against the role intended for each source. |

## 4. Cell-line evidence (kept separate; neither source chosen)
| | **Sanger / Cell Model Passports** | **Broad / DepMap / CCLE** |
|---|---|---|
| Documented licence / terms | Sanger Cancer Dependency Map Data Usage Policy: "non-exclusive, non-transferable right to use data files for internal proprietary research and educational purposes, including target, biomarker and drug discovery"; excludes "resale … or … provision of commercial services"; "Commercial use … is not permitted without prior consent." The Cell Model Passports page itself states no terms; which datasets the policy covers is not enumerated. | Figshare records of DepMap Public 23Q2, 23Q4 and 24Q2 (published by Broad DepMap) carry **CC BY 4.0** (later releases not checked). The DepMap/CCLE portal terms are behind a Cloudflare challenge and **unread**; a forum post says most CCLE files use CC BY 4.0 with some unclear; raw CCLE sequence terms are not established; a search snippet asserting a non-commercial restriction is unverified and conflicts with the Figshare licence (neither adopted nor dismissed). |
| Unresolved commercial/diagnostic question | Consent for commercial use is required by the text; diagnostic development is not addressed. | Whether portal or raw-data terms add restrictions to the CC BY 4.0 release files; diagnostic development is not addressed by the licence. |
| Suitability | The labels are computational (MSIsensor-pro score ≥7 = MSI), so the evidence is suitable **only as development/control material, and only pending authorization** for the intended use. | CCLE labels are computational (deletion-count classification); suitable **only as development/control material, pending resolution of the terms**. |
| Determination required | Data provider terms: Sanger's written position on the intended use (contact named in its policy; not contacted). | Data provider terms: the Broad portal and CCLE terms as they apply to the files actually contemplated; legal/compliance review of any conflict between portal terms and the CC BY 4.0 licence. |

## 5. Required external decisions (checklist)
| ID | Decision or document needed | From | Status |
|---|---|---|---|
| G-01 | NCI DAC decision on a TCGA (phs000178) Data Access Request describing the intended use | NCI DAC | **OPEN — DAC DECISION REQUIRED** |
| G-02 | Reading of DUC terms 1, 6 and 9 against a proprietary diagnostic-panel product and any regulatory submission (TCGA) | GenRichi legal/compliance | **OPEN — LEGAL/COMPLIANCE REVIEW REQUIRED** |
| G-03 | Confirmation that GenRichi can act as Requester (Signing Official, eRA Commons) with an eligible PI | GenRichi legal/compliance | **OPEN — LEGAL/COMPLIANCE REVIEW REQUIRED** |
| G-04 | Current CPTAC Data Use Agreement text and dataset-specific publication status | NCI / CPTAC | **OPEN — DATA PROVIDER TERMS REQUIRED** |
| G-05 | Applicable repository, accession and access tier for CPTAC colon raw sequence data | NCI / CPTAC | **OPEN — DATA PROVIDER TERMS REQUIRED** |
| G-06 | NCI DAC decision on a CPTAC (phs001287) Data Access Request describing the intended use | NCI DAC | **OPEN — DAC DECISION REQUIRED** |
| G-07 | Review of the DUC and the CPTAC-specific term against the intended product | GenRichi legal/compliance | **OPEN — LEGAL/COMPLIANCE REVIEW REQUIRED** |
| G-08 | Independence evidence between the TCGA and the CPTAC prospective cohorts (statement of provenance or authorised patient-level cross-check), sufficient for each source's intended role | NCI / CPTAC or authorised Approved User | **OPEN — INDEPENDENCE EVIDENCE REQUIRED** |
| G-09 | Sanger's written position on commercial/diagnostic use of Cell Model Passports / Dependency Map data | Sanger | **OPEN — DATA PROVIDER TERMS REQUIRED** |
| G-10 | Broad DepMap portal terms and CCLE terms (release files and raw sequence) as they apply to the contemplated files | Broad Institute | **OPEN — DATA PROVIDER TERMS REQUIRED** |
| G-11 | Review of any conflict between DepMap portal terms and the CC BY 4.0 licence | GenRichi legal/compliance | **OPEN — LEGAL/COMPLIANCE REVIEW REQUIRED** |

No item above is approved, permitted or validated.

## 6. Gate logic
SOURCE SELECTION STATUS: BLOCKED
NO DATA ACCESS AUTHORIZED
NO SOURCE SELECTED

**Condition for reopening Source Selection:** all applicable data-use questions must be resolved by the competent authority or provider, and independence must be established sufficiently to support the intended role of each source.
The v5 items classified "must resolve before data access" and the v5 NOT USABLE FOR THIS STAGE dispositions (UV-05, UV-13, UV-21) are unchanged by this package.

Next Gate: Governance Resolution → Source Selection Review
