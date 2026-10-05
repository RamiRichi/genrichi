# Stage 6 — CPTAC / dbGaP Governance Resolution (governance documentation only)

Question: **Is the governance basis for using CPTAC controlled-access data for the intended GenRichi Stage 6 activity sufficiently established?**
Only public documentation and public metadata were read. No patient-level or controlled-access CPTAC data was downloaded, requested, opened, copied or processed; nobody (NCI, dbGaP, CPTAC, SRA, any provider) was contacted; no data-access request or DAC request was submitted;
no cohort was selected; CPTAC is not compared with any other source; no code, classifier, threshold or production file was touched. Nothing here is legal advice. Public files or a published paper are not treated as permission.
`STAGE6_TCGA_GOVERNANCE_RESOLUTION.md` was not reopened or modified. Date: 2026-09-25.

## 1. Scope
**Use being evaluated:** access to CPTAC controlled-access, sequence-level data (exome/genome reads and related files) by GenRichi to generate Stage 6 empirical MSI evidence for developing a proprietary, locus-based MSI panel intended for a future diagnostic product. Categories, kept separate:

| Category | Applies? |
|---|---|
| Research/development validation (empirical evidence on candidate MSI loci, development/RUO, no clinical claim) | **Yes — the activity evaluated** |
| Diagnostic product development (the programme the evidence feeds) | **Yes, as the purpose behind the activity; a distinct category** |
| Clinical validation (clinical sensitivity/specificity, clinical thresholds) | **No — not covered, no claim made** |

Research-use permission is not treated as clinical-validation permission, and commercial-product language is not treated as diagnostic-development permission.

## 2. CPTAC access routes (each treated separately; no equivalence assumed)
Which accession holds which cohort matters, so this review checked GDC project metadata. Findings that widen the picture of the requested routes are flagged.

### A. dbGaP `phs001287` — "CPTAC Proteogenomic Study" (GDC project CPTAC-3)
| Item | Finding | Source |
|---|---|---|
| Accession / project | dbGaP phs001287.v23.p7; GDC project `CPTAC-3` ("CPTAC-Brain, Head and Neck, Kidney, Lung, Pancreas, Uterus"), dbGaP accession recorded as phs001287 | dbGaP study page `…/study.cgi?study_id=phs001287`; GDC API `https://api.gdc.cancer.gov/projects/CPTAC-3` |
| Access tier | Controlled access via dbGaP "Authorized Access" with a Data Use Certification; GDC lists CPTAC-3 files as open or controlled (110,077 files: 85,906 controlled, 24,171 open at the time of query); controlled BAM/VCF need dbGaP authorization | GDC API file counts by `access`; GDC `https://docs.gdc.cancer.gov/Encyclopedia/pages/Controlled_Access` |
| Current agreement | dbGaP DUC Agreement "Ver. January 25, 2025" with the "Study specific DUC addendum phs001287" | `https://dbgap.ncbi.nlm.nih.gov/aa/wga.cgi?view_pdf=&stacc=phs001287.v23.p7` (text read) |
| Authorization requirements | PI with eRA Commons account; Data Access Request signed with an Institutional Signing Official; approval by the **NCI DAC**; collaborators at other institutions submit their own DAR | DUC; GDC access page |
| Approved-use restriction | Consent group 1 "General Research Use (GRU)": "Use of the data is limited only by the terms of the model Data Use Certification." "Research use will occur solely in connection with the approved research project described in the DAR" (Research Use Statement); new uses need a new DAR; one-year renewable access | DUC Term 1; Addendum |
| Publication / IP | Public posting of genomic summary results "Allowed"; publication encouraged with acknowledgement and accession phs001287 plus retrieval date; **CPTAC freedom-to-publish rule**: no limitation once (1) a global-analysis publication exists, or (2) 22 months after all samples were received at the data production site, or (3) CPTAC Steering Committee approval (dataset status shown on the study page "Available Data" tab, not read). DUC Term 9: NIH "considers these data as pre-competitive and urges Approved Users to avoid making IP claims derived directly from the dataset(s)"; data and conclusions "will remain freely available, without requirement for licensing". Datasets and Data Derivatives "may not be sold to any individual at any point in time for any purpose". | DUC Terms 6, 9, 10; Addendum |
| Commercial-use language | DUC Term 9: "there is no restriction on development of commercial products resulting from the knowledge gained from the research project." No language on diagnostic or clinical development, IVD or regulatory use, or for-profit requester eligibility. | DUC Term 9 |
| CPTAC-specific term | The addendum's "IC Specific Access Term" requires users to "recognize any restrictions on data use established by the submitting institution through the Institutional Certification and stated on the CPTAC Data Use Agreement page" (`https://proteomics.cancer.gov/data-portal/about/data-use-agreement`). **That URL now redirects (HTTP 301) to `https://dctd.cancer.gov/programs/occpr`, an overview page containing no data-use terms.** | Addendum; page fetch |
| Requester eligibility | Requester = "home institution or organization" via an Institutional Signing Official; PI = "a permanent employee of their institution at a level equivalent to a tenure-track professor or senior scientist"; for-profit status not addressed | DUC Definitions |
| Terms actually available and authoritative? | DUC and addendum: **available, read, authoritative**. The CPTAC Data Use Agreement referenced by the addendum: **not available at its cited URL**. | as above |

### B. SRA BioProject `PRJNA514017` (cited by the CPTAC colon paper)
| Item | Finding | Source |
|---|---|---|
| Accession | BioProject PRJNA514017 "Proteogenomic characterization of human colon cancer reveals new therapeutic opportunities"; data type "Raw sequence reads"; registered 2019/01/09; submitter Clinical Proteomic Tumor Analysis Consortium / Baylor College of Medicine; `sequencing_status: Unknown` | NCBI E-utilities BioProject summary (public metadata) |
| Access tier | **Not stated.** Public NCBI indexes return **0 SRA records** for `PRJNA514017[BioProject]` and for `SRP178677`, and 0 BioSample records (the BioProject page mentions 220 linked BioSamples and "No public data is linked to this project"). This shows no publicly indexed reads; it does **not** establish whether the data are unreleased, controlled through dbGaP, or otherwise restricted. | E-utilities counts, queried 2026-09-25; BioProject page |
| Current applicable terms | **None identified.** No agreement, data-use limitation or consent group could be linked to this accession. NCBI documentation states that sequence data for dbGaP studies must be submitted through dbGaP and that human data without explicit consent for public release are kept behind a controlled-access firewall — general statements that do not say what applies to this project. | NCBI SRA docs `https://www.ncbi.nlm.nih.gov/sra/docs/submitdbgap/` (search-result text) |
| Authorization / approved use / publication-IP / commercial language / requester eligibility | **Not established** for this route. | — |
| Where the route is cited | The colon paper's data-availability statement: raw genomics "available at the Sequence Read Archive, BioProject ID PRJNA514017", with a review-only FTP link. A publication is not treated as proof of current terms or access. | Vasaikar 2019, Data Availability (PMC6768830) |
| Terms available and authoritative? | **No.** | — |

### C. dbGaP `phs000892` — identified during this review (GDC project CPTAC-2)
Not in the requested list, but it appears to be the controlled-access study for colon exome data, so it is documented separately instead of being assumed equivalent to A or B.
| Item | Finding | Source |
|---|---|---|
| Accession / project | dbGaP phs000892.v7.p1 "CPTAC Proteogenomic Confirmatory Study of Breast, Colon, Lung, and Ovarian Tumors" (348 consented subjects); GDC project `CPTAC-2` ("CPTAC-Breast, Colon, Ovary"), dbGaP accession recorded as phs000892. The study page says that since November 2019 primary data are available at the GDC, with genotyping files still hosted at dbGaP, and describes controlled ischaemic time under 30 minutes, whole-exome from ACGT Inc., germline DNA from blood, adjacent colon tissue as normal for colon cases. | dbGaP `…/study.cgi?study_id=phs000892`; GDC API `https://api.gdc.cancer.gov/projects/CPTAC-2` |
| Access tier | Controlled via dbGaP Authorized Access; GDC CPTAC-2: 9,244 files, 7,956 controlled and 1,288 open at the time of query. | GDC API |
| Agreement | DUC "Ver. January 25, 2025" with "Study specific DUC addendum phs000892" (`…stacc=phs000892.v7.p1`): consent group 1 GRU ("limited only by the terms of the model Data Use Certification"); public posting of genomic summary results "Allowed"; acknowledgement statement with grant numbers; **no CPTAC "IC Specific Access Term"** in this addendum. Terms 1, 6 and 9 as in A. | DUC PDF (text read) |
| Whether the 110-patient colon cohort of the CPTAC colon paper is (part of) phs000892 | **Not stated in any authoritative source read.** The design description resembles the paper's cohort, but the paper cites SRA, not this accession. Mapping is unconfirmed. | — |

## 3. Evidence quality
| Conclusion | Class |
|---|---|
| phs001287 and phs000892 sequence-level data are controlled-access with NCI DAC approval; GRU consent group; research use limited to the approved DAR | **DIRECT AUTHORITATIVE EVIDENCE** (DUC + addenda, dbGaP/GDC pages) |
| DUC Term 9 wording on commercial products, IP, free availability; no sale of Data Derivatives | **DIRECT AUTHORITATIVE EVIDENCE** |
| DUC contains no diagnostic/clinical/regulatory or for-profit-eligibility language | **DIRECT AUTHORITATIVE EVIDENCE** (absence in the text read; not treated as permission) |
| What GenRichi's intended use means under Terms 1, 6 and 9, and whether Stage 6 evidence tables are Data Derivatives or "genomic summary results" | **REQUIRES LEGAL REVIEW** |
| Whether the intended use is an "approved research project" | **REQUIRES PROVIDER/DAC CONFIRMATION** (NCI DAC) |
| Content of the CPTAC Data Use Agreement cited by the phs001287 addendum | **UNVERIFIED** (cited page unavailable) — **REQUIRES PROVIDER/DAC CONFIRMATION** |
| Which accession/route holds the colon cohort's exome data | **UNVERIFIED** — **REQUIRES PROVIDER/DAC CONFIRMATION** |
| Access tier and terms of SRA `PRJNA514017` | **UNVERIFIED** — **REQUIRES PROVIDER/DAC CONFIRMATION** |
| The colon paper's statement of where raw data are | **SECONDARY/PUBLISHED EVIDENCE** (not proof of current terms) |
| 22-month freedom-to-publish rule (addendum) vs 15-month rule (retired CPTAC portal FAQ) | **Discrepancy preserved**; the addendum text is the current authoritative statement and the FAQ is a retired-portal page (**SECONDARY**) |
| Requester (for-profit organization) eligibility | **UNVERIFIED** — **REQUIRES ADDITIONAL AUTHORITATIVE EVIDENCE** |

## 4. Decision matrix
| Question | Evidence | Interpretation | Status |
|---|---|---|---|
| Is controlled access required? | GDC/dbGaP: controlled BAM/VCF; DUC for phs001287 and phs000892 | Yes for the dbGaP/GDC routes (A, C). For route B the tier is not established (see §4b). | **RESOLVED** (routes A, C) |
| Is the intended Stage 6 use covered by the applicable research purpose? | DUC Term 1; consent group GRU; approval by NCI DAC | No approved research purpose exists yet; coverage cannot be determined from the text | **REQUIRES DAC/PROVIDER DECISION** |
| Are commercial/for-profit uses addressed? | Term 9 addresses commercial products resulting from the research; nothing on for-profit requesters; the CPTAC-specific term (phs001287) is unread | Commercial products are addressed only in the Term 9 sense; requester status and the CPTAC-specific term remain unaddressed | **REQUIRES ADDITIONAL AUTHORITATIVE EVIDENCE** |
| Does any commercial language resolve diagnostic-development use? | Term 9 only | No; it does not address diagnostic or clinical development, nor how "avoid making IP claims" and "freely available" apply to a proprietary panel | **REQUIRES LEGAL REVIEW** |
| Is DAC/provider approval required? | DUC effective date is the DAC approval date; GDC: DAC approves or disapproves | Yes for routes A and C; not established for route B | **RESOLVED** (requirement established; no decision obtained) |
| Is legal review required? | Terms 1, 6, 9; diagnostic-product purpose; Data Derivatives and summary-results posting questions | Interpretations beyond the text are a matter for counsel/compliance | **REQUIRES LEGAL REVIEW** |
| Is requester eligibility established? | DUC Definitions; no for-profit statement; GenRichi's own facts not documented in the reviewed material | Not established | **REQUIRES ADDITIONAL AUTHORITATIVE EVIDENCE** |
| Are the applicable current terms actually available? | DUC and addenda for A and C available; CPTAC Data Use Agreement (A) unavailable at its URL; no terms found for B | Partly: the CPTAC-specific term for A and all terms for B are not available | **REQUIRES ADDITIONAL AUTHORITATIVE EVIDENCE** |
| *(added)* Is the access route for the colon cohort's exome data established? | Paper cites SRA (B); GDC/dbGaP hold a colon project under phs000892 (C); phs001287 (A) holds mainly other tumour types | Not established which route governs the colon cohort | **REQUIRES ADDITIONAL AUTHORITATIVE EVIDENCE** |

### 4b. Route B (SRA `PRJNA514017`) separately
| Question | Status |
|---|---|
| Is controlled access required? | **REQUIRES ADDITIONAL AUTHORITATIVE EVIDENCE** |
| Applicable terms available? | **REQUIRES ADDITIONAL AUTHORITATIVE EVIDENCE** |
| DAC/provider approval required? | **REQUIRES ADDITIONAL AUTHORITATIVE EVIDENCE** |

## 5. Final gate
CPTAC GOVERNANCE STATUS: BLOCKED

Nothing prohibits the use and nothing states that it is permitted; CPTAC is **not** approved, permitted or validated for GenRichi. Exact remaining blockers:
1. **C1 — CPTAC-specific term (phs001287):** the CPTAC Data Use Agreement cited by the addendum is not available at its URL (UNVERIFIED).
2. **C2 — Route identification:** which of `phs001287`, `phs000892` and SRA `PRJNA514017` governs the colon cohort's exome data, and the access tier and terms of `PRJNA514017` (UNVERIFIED).
3. **C3 — DAC decision:** whether the intended use is an approved research project (NCI DAC).
4. **C4 — Legal review:** DUC Terms 1, 6 and 9 for a diagnostic-product context, including Data Derivatives and genomic-summary-results questions.
5. **C5 — Requester eligibility evidence:** the provider's position on for-profit requesters and GenRichi's own eligibility facts.

## 6. Next action (not performed)
The minimum concrete action is a **written provider confirmation from NCI/CPTAC** (Office of Cancer Clinical Proteomics Research) covering two facts: (a) the current text of the CPTAC Data Use Agreement referred to in the `phs001287` addendum, and (b) which accession/route (`phs001287`, `phs000892` or SRA `PRJNA514017`) holds the CPTAC colon cohort's exome data and its access tier and terms. The internal legal/compliance determination on use description, Research Use Statement and requester eligibility (already named as the next action in the TCGA resolution) is a separate, independent step; no request is made by this document.
