# Stage 6 — TCGA / dbGaP Governance Resolution (phs000178) (governance documentation only)

Question: **Is the governance basis for using TCGA controlled-access data for the intended GenRichi Stage 6 activity sufficiently established?**
This document only records what authoritative public documents say. No patient-level or controlled-access data was downloaded, requested, opened, copied or processed; nobody was contacted; no DAC
request was submitted; no cohort selected; no classifier, threshold, validation code or production code touched; TCGA is not compared with CPTAC or MSK; nothing here is legal advice. Where an interpretation needs a DAC
decision or legal review, it is marked so. Date: 2026-09-25.

## 1. Scope
**Use being evaluated:** access to TCGA **controlled-access, sequence-level data** (e.g. BAM files of tumour/normal exomes; dbGaP `phs000178`) by GenRichi, to generate the Stage 6 empirical MSI evidence for the development of a
proprietary, locus-based MSI panel intended for a future diagnostic product. It does **not** cover open-tier TCGA clinical or somatic-mutation tables (not needed to answer this question), other cohorts, or which cohort would be chosen.

## 2. Authoritative access terms
| Statement | Authoritative source |
|---|---|
| Study/accession: **phs000178**, The Cancer Genome Atlas; the study page identifies **phs000178.v11.p8** as the current version and states that a Data Use Certification Agreement is available for it | dbGaP study page `https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000178` (read through a fetch-tool summary) |
| Agreement read: dbGaP **Data Use Certification (DUC) Agreement, "Ver. January 25, 2025"**, with the "Study specific DUC addendum phs000178" | DUC PDF for phs000178.v11.p8: `https://dbgap.ncbi.nlm.nih.gov/aa/wga.cgi?view_pdf=&stacc=phs000178.v11.p8` (text read directly). Whether a still newer version exists was not verified. The 20 Aug 2014 TCGA DUC read earlier is superseded and is not relied on. |
| Access tier: controlled-access data include "Primary sequence data (.bam files)"; access requires "user certification through dbGaP Authorized Access"; DAC = **NCI DAC** | dbGaP study page (fetch summary); DUC Addendum; GDC `https://docs.gdc.cancer.gov/Encyclopedia/pages/Controlled_Access` ("raw sequencing data (BAM, FASTQ files), VCF files, and protected MAF files"; dbGaP DAC approval and an eRA Commons ID) |
| Access steps: NIH eRA Commons account for the PI; dbGaP authorization; a research project created in dbGaP; then GDC login once approved | GDC `https://gdc.cancer.gov/access-data/obtaining-access-controlled-data` (fetch summary) |
| Consent group 1, "General Research Use (GRU)": "Use of the data is limited only by the terms of the model Data Use Certification." | DUC Addendum phs000178 |
| Approved-use restriction: "Research use will occur solely in connection with the approved research project described in the DAR, which includes a 1-2 paragraph description of the proposed research (i.e., a Research Use Statement)"; new uses "will require submission of a new DAR" | DUC Term 1 (Research Use) |
| Commercial-product language (verbatim): NIH "considers these data as pre-competitive and urges Approved Users to avoid making IP claims derived directly from the dataset(s)"; "these NIH-provided data, and conclusions derived therefrom, will remain freely available, without requirement for licensing"; "there is no restriction on development of commercial products resulting from the knowledge gained from the research project"; ownership of IP "will be governed by applicable patent law" | DUC Term 9 (Intellectual Property) |
| Sale/distribution: datasets and Data Derivatives "may not be sold to any individual at any point in time for any purpose"; no distribution to entities not identified in the approved request without NIH written approval. "Data Derivative": "Data derived from controlled-access datasets … Examples … imputed datasets and single nucleotide polymorphisms, or any data explicitly designated as Data Derivatives by NIH." | DUC Terms 6/7 and Definitions |
| Publication and sharing: publication encouraged; acknowledgement including the dbGaP accession; addendum: "Public Posting of Genomic Summary Results — Not Allowed" (the definition of "genomic summary results" was not confirmed in the text read) | DUC Term 10; Addendum |
| Duration and reporting: access for one year, renewable; annual Progress Update reporting "any downstream intellectual property generated from the data"; destruction at Project Close-out; security best practices; non-identification | DUC Terms 1, 3, 7, 11, 13 |
| Requester and PI requirements: Requester = "home institution or organization" acting through an Institutional Signing Official; PI = "a permanent employee of their institution at a level equivalent to a tenure-track professor or senior scientist"; collaborators at other institutions submit their own DAR; contractors must be described in the Research Use Statement | DUC Definitions; Term 1 |
| Institutional Certification: certifies "the appropriate research uses of the data and the uses that are specifically excluded by the relevant informed consent documents" | DUC Definitions |
| The DUC contains no language on diagnostic or clinical development, IVD use, regulatory submission, or for-profit eligibility (searches for these terms in the text returned nothing relevant) | DUC text as read |

**Not readable by my tools (recorded, not bypassed):** the NIH Genomic Data Sharing Policy notice (NOT-OD-14-124) and `sharing.nih.gov` (HTTP 403); the dbGaP FAQ and "Applying for Controlled Access Data" pages (bot challenge); the older dbGaP request-procedures PDF read shows no commercial-entity language.

## 3. Intended GenRichi use (categories kept separate)
| Category | Applies to this evaluation? |
|---|---|
| **Research/development validation** — generating empirical evidence on candidate MSI loci from existing exome data, development/RUO, no clinical performance claim | **Yes — this is the activity being evaluated (Stage 6).** |
| **Diagnostic product development** — the programme this evidence feeds (a future proprietary MSI panel for a diagnostic product) | **Yes, as the purpose behind the activity.** It is a distinct category from the line above, and how it is described in any request is a decision for GenRichi and its legal/compliance function. |
| **Clinical validation** — establishing clinical sensitivity/specificity or clinical thresholds | **No.** Stage 6 makes no such claim and this evaluation does not cover it. |
| Regulatory-evidence use of the data (e.g. inclusion in a regulatory submission) | **Not evaluated** here; a separate, later question. |

## 4. Decision matrix
| Question | Evidence | Interpretation | Status |
|---|---|---|---|
| Is controlled access required? | GDC and dbGaP list BAM/FASTQ/VCF files as controlled; DUC Terms of Access | Yes, for sequence-level data of the kind Stage 6 needs; this is what the agreement text states | **RESOLVED** |
| Is GenRichi's intended use covered by the approved research purpose? | DUC Term 1: use "solely in connection with the approved research project described in the DAR"; consent group GRU; approval is by the NCI DAC | No approved research purpose exists yet; whether the described use is covered is a determination the DAC makes on a submitted Research Use Statement, and cannot be decided from the text | **REQUIRES DAC DECISION** |
| Does commercial-product language resolve diagnostic-development use? | DUC Term 9 (quoted in §2) | No. The language concerns commercial products "resulting from the knowledge gained from the research project"; it says nothing on whether controlled data may be used inside diagnostic-product development or validation, on clinical or regulatory use, or on how "urges … to avoid making IP claims derived directly from the dataset(s)" and "conclusions derived therefrom will remain freely available" apply to a proprietary panel. Also open: whether Stage 6 evidence tables would be "Data Derivatives" or "genomic summary results" | **REQUIRES LEGAL REVIEW** |
| Is DAC approval required? | GDC ("Each project has a Data Access Committee (DAC) that will approve or disapprove data access requests"); DUC (effective date = DAC approval date) | Yes, approval by the NCI DAC of a Data Access Request is required before access; it has not been sought or obtained | **RESOLVED** (requirement established; the decision itself is not obtained — see row 2) |
| Is legal review required? | DUC binds the Requester through its Signing Official; Terms 1, 6, 9 as above; the intended use is diagnostic-product development | The agreement does not answer the diagnostic-development question, and any reading beyond its text is a legal interpretation for GenRichi's counsel/compliance, not an assumption | **REQUIRES LEGAL REVIEW** |
| Is requester eligibility established? | DUC definitions (Requester = home institution or organization via Signing Official; PI credential definition; eRA Commons); no statement on for-profit organizations; GenRichi's own facts (Signing Official, eRA Commons registration, a PI meeting the definition) are not documented in the material reviewed | Not established either way; the missing statement is the provider's position on commercial/for-profit requesters, plus GenRichi's own facts | **REQUIRES ADDITIONAL AUTHORITATIVE EVIDENCE** |

No row is "NOT APPLICABLE".

## 5. Final gate
TCGA GOVERNANCE STATUS: BLOCKED

The basis for using TCGA controlled-access data for the intended Stage 6 activity is **not** sufficiently established. Nothing in the documents prohibits it, and nothing states that it is permitted. Exact remaining blockers:
1. **B1 — DAC decision** on whether the described use is an approved research project (row 2).
2. **B2 — Legal review** of DUC Terms 1, 6 and 9 for a diagnostic-product development context, including whether Stage 6 evidence tables are Data Derivatives or "genomic summary results" (rows 3 and 5).
3. **B3 — Requester eligibility evidence:** the provider's statement on commercial/for-profit requesters and GenRichi's own eligibility facts (row 6).

Not declared: TCGA is **not** approved, permitted or validated for GenRichi.

## 6. Next action (not performed)
The minimum concrete action is an internal **legal/compliance determination for GenRichi** that (i) states in one agreed sentence what TCGA controlled data would be used for, (ii) records whether that description can be truthfully submitted as a dbGaP Research Use Statement, and (iii) confirms or denies the eligibility facts (Signing Official, eRA Commons, qualifying PI). Only after that can the owner decide whether to submit a DAR, which is the only route to the DAC decision (B1). No request to NCI/dbGaP is made by this document.
