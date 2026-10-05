# Stage 6 — Governance Gate Status (consolidation; documentation/control only)

This document consolidates four completed governance reviews. It performs no new research and no legal analysis beyond accurately summarizing them. No patient-level, controlled-access, raw-sequencing or cell-line data were
accessed; no access was requested; no provider (NCI, CPTAC, TCGA, Sanger, Broad, dbGaP, SRA) was contacted; no DAC request was submitted; no code, classifier, threshold, validation or production file was touched; no validation cohort was selected;
**no dataset is claimed to be authorized for GenRichi.** It is a status consolidation, **not legal advice**. Statuses below are reproduced exactly as recorded in the source documents.

## 1. Current status
**STAGE 6 GOVERNANCE GATE: NOT CLOSED**

Meaning: the **internal** governance review is complete — the evidence has been gathered, the terms read where they were readable, and the gaps and discrepancies documented — but **external authorization, provider confirmation and legal/compliance determinations remain outstanding**.
Nothing in the four reviews states that any use is approved, permitted or validated for GenRichi; commercial-product language in the data-use agreements has **not** been read as permission for diagnostic-product development.

## 2. Gate summary
| Gate | Current status | Remaining condition (taken from the corresponding document) |
|---|---|---|
| TCGA data-use (`phs000178`) | **BLOCKED** | **B1** NCI DAC decision on whether the described use is an approved research project. **B2** legal review of DUC Terms 1, 6 and 9 for a diagnostic-product context, including whether Stage 6 evidence tables are Data Derivatives or "genomic summary results". **B3** requester-eligibility evidence: the provider's position on commercial/for-profit requesters and GenRichi's own eligibility facts. |
| CPTAC data-use (`phs001287`; also `phs000892` and SRA `PRJNA514017` as separate routes) | **BLOCKED** | **C1** the CPTAC Data Use Agreement cited by the phs001287 addendum is unavailable at its URL. **C2** which of phs001287, phs000892 and SRA PRJNA514017 governs the colon cohort's exome data, and the tier and terms of PRJNA514017. **C3** NCI DAC decision on whether the intended use is an approved research project. **C4** legal review of DUC Terms 1, 6 and 9 for a diagnostic-product context, including Data Derivatives and genomic-summary-results questions. **C5** requester-eligibility evidence. |
| Cell-line data-use (Sanger and Broad/DepMap, kept separate) | **CONDITIONAL** | All six conditions must hold: (1) cell-line data used only as a development/control resource, never as patient-level ground truth; (2) DepMap portal Terms, CCLE terms and the release README (unread) checked and shown not to add restrictions to CC BY 4.0 for the files contemplated; (3) Sanger data used only within "internal proprietary research and educational purposes", with a legal/compliance determination that the intended use falls within it, and prior Sanger consent for any product incorporation, commercial service or diagnostic use; (4) no raw sequencing data used without its specific terms (EGA DAC/DAA for EGAD00001001039; Broad/submitter terms for PRJNA523380); (5) no diagnostic-development use assumed permitted; (6) third-party-rights language reviewed before any use beyond internal research. |
| TCGA-CPTAC independence | **BLOCKED** | Missing: (1) an authoritative provenance statement from NCI/CPTAC/the Biospecimen Core Resource stating, per cohort and accession, whether any CPTAC patient is also a TCGA participant; (2) or an authorized patient-level identity cross-check (not permitted now); (3) a defined cohort↔accession mapping for the colon cohort (110 patients / 106 exomes / 105 colon-site + 1 rectum-site GDC cases; phs000892 vs phs001287 vs SRA PRJNA514017) and for the endometrial cohort (95 tumours vs 241 GDC "Uterus, NOS" cases). |

## 3. External decisions required (consolidated, no duplication)
**A. DAC / provider decisions**
| ID | Decision or document | From |
|---|---|---|
| A1 | Decision on a TCGA (`phs000178`) Data Access Request describing the intended use | NCI DAC |
| A2 | Decision on a CPTAC Data Access Request describing the intended use, for the applicable accession(s) | NCI DAC |
| A3 | Current text of the CPTAC Data Use Agreement referenced in the phs001287 addendum | NCI / CPTAC |
| A4 | Which accession/route (phs001287, phs000892, SRA PRJNA514017) holds the CPTAC colon cohort's exome data, with its access tier and terms (this also supplies the colon half of the cohort↔accession mapping in D2) | NCI / CPTAC |

**B. Legal / compliance decisions**
| ID | Decision | Scope |
|---|---|---|
| B1 | Reading of DUC Terms 1 (research use within the approved DAR), 6 (no sale of Data Derivatives) and 9 (IP; "avoid making IP claims derived directly from the dataset(s)"; data and conclusions "freely available") against a proprietary diagnostic-panel product and any regulatory submission; whether Stage 6 evidence tables are Data Derivatives or "genomic summary results" (posting is "Not Allowed" in the TCGA addendum and "Allowed" in the CPTAC addenda) | TCGA and CPTAC dbGaP agreements |
| B2 | Review of the CPTAC-specific term (item A3) against the intended product once its text is available | CPTAC |

**C. Requester eligibility**
| ID | Decision or evidence | From |
|---|---|---|
| C1 | The provider's position on commercial/for-profit requesters (the agreements are silent) | dbGaP / NIH |
| C2 | GenRichi's own eligibility facts: an Institutional Signing Official, eRA Commons registration, and a PI meeting the agreement's definition | GenRichi legal/compliance |

**D. Provenance / independence confirmation**
| ID | Decision or evidence | From |
|---|---|---|
| D1 | Written provenance statement stating, per cohort and accession, whether any CPTAC patient is also a TCGA participant — or an authorized patient-level identity cross-check performed by an Approved User under approved access | NCI / CPTAC / Biospecimen Core Resource; or authorized Approved User |
| D2 | Cohort↔accession mapping: colon (see A4) and endometrial (relationship between the paper's 95 tumours and the 241 GDC cases) | NCI / CPTAC |

**E. Cell-line authorization conditions**
| ID | Condition | From |
|---|---|---|
| E1 | Role limited to development/control; never patient-level ground truth | GenRichi (internal control) |
| E2 | DepMap portal Terms, CCLE terms and release README read and shown not to restrict the CC BY 4.0 release files contemplated | Broad / DepMap terms (read by a person with browser access) |
| E3 | Legal/compliance determination that the intended use is "internal proprietary research" under the Sanger policy; Sanger prior consent for any product incorporation, commercial service or diagnostic use | GenRichi legal/compliance; Sanger |
| E4 | Raw-data terms: EGA DAC/Data Access Agreement for EGAD00001001039; Broad/submitter terms for SRA PRJNA523380 — needed only if raw data are ever contemplated | Sanger DAC; Broad |
| E5 | Diagnostic-development use of cell-line resources not assumed permitted; kept as an open legal question per resource | GenRichi legal/compliance |
| E6 | Third-party-rights review (CC BY 4.0 rights not covered; cell-bank/CCLE rights) before any use beyond internal research | GenRichi legal/compliance |

## 4. What is already closed internally
* The governance **evidence review has been performed** for TCGA, CPTAC, the cell-line resources and TCGA–CPTAC independence.
* **Known limitations and discrepancies are documented and preserved**, including: the CPTAC endometrial count (25 MSI-H of 95 in the primary paper and public table versus 21 in a 2026 secondary source); the CPTAC colon 110/106/105+1 differences; 95 versus 241 CPTAC-3 uterus cases in open GDC metadata; the 22-month (addendum) versus 15-month (retired portal FAQ) CPTAC publication rule; the unreadable CPTAC Data Use Agreement page; unreadable Broad DepMap/CCLE portal terms.
* **No patient-level controlled data were accessed.**
* **No validation cohort was selected.**
* **No data were used as ground truth.**
* **No code or classifier changes** were made as part of this governance work.
* Internal contradictions between earlier and later documents were reconciled by supersession: the gate closure's "UNRESOLVED" for UV-24 (cell lines) is superseded by the cell-line document's CONDITIONAL, and its "UNRESOLVED" for UV-18 by the independence document's BLOCKED, each on the basis of additional evidence documented there.

## 5. What must NOT happen yet
* **Cohort selection** (development or validation).
* **Patient-level data access** of any kind (including requesting or opening controlled data).
* **MSI empirical validation** (Stage 6B) on any patient-derived or cell-line data.
* **Threshold tuning** (the provisional prototype parameters remain provisional and untuned).
* **Classifier training** or modification.
* **Treating any public label as validation ground truth** (all public MSI labels reviewed are computational or of unverified provenance and none is authorized for this use).

## 6. Next gate
**GOVERNANCE RESOLUTION → SOURCE SELECTION REVIEW**

Source Selection Review **may begin only after the applicable external authorization and provenance conditions above are resolved** by the competent authority or provider (and, where noted, by GenRichi's legal/compliance function), and independence is established sufficiently to support the intended role of each source.
Source Selection is **not** performed here.

## 7. Decision log
| Document | Recorded status |
|---|---|
| `STAGE6_TCGA_GOVERNANCE_RESOLUTION.md` | TCGA GOVERNANCE STATUS: **BLOCKED** |
| `STAGE6_CPTAC_GOVERNANCE_RESOLUTION.md` | CPTAC GOVERNANCE STATUS: **BLOCKED** |
| `STAGE6_CELLLINE_GOVERNANCE_RESOLUTION.md` | CELLLINE GOVERNANCE STATUS: **CONDITIONAL** |
| `STAGE6_TCGA_CPTAC_INDEPENDENCE_RESOLUTION.md` | TCGA-CPTAC INDEPENDENCE STATUS: **BLOCKED** |
| `STAGE6_GOVERNANCE_GATE_STATUS.md` (this document) | STAGE 6 GOVERNANCE GATE: **NOT CLOSED** |

Date of consolidation: **2026-09-25**. This document is a status consolidation, **not legal advice**, and does not claim that any dataset is authorized for GenRichi.

STAGE 6 GOVERNANCE GATE: NOT CLOSED — EXTERNAL RESOLUTION REQUIRED
