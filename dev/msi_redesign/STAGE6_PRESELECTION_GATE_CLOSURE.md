# Stage 6 — Pre-selection Gate Closure (UV-15, UV-16, UV-18, UV-24) (development / RUO)

Scope: only the four items classified in v5 as MUST RESOLVE BEFORE SOURCE SELECTION. Only public documentation was read. No patient-derived, controlled-access or sequence
data were accessed or downloaded; no sample selected; no cohort constructed; no threshold, minimum, partition or Stage 6B Analysis Plan defined. Nothing was inferred from silence in
an agreement, and independence was not inferred from the absence of an overlap statement. No 403, anti-automation or authentication control was bypassed. This is a documentation
finding, **not legal advice**. Date: 2026-09-25.

**Access limits encountered (recorded, not circumvented):** the DepMap portal terms pages (`depmap.org/portal/terms/`, `/portal/ccle/terms_and_conditions`) are behind a Cloudflare challenge;
the NIH policy notice (NOT-OD-14-124), `sharing.nih.gov` and the MSK data catalogue returned 403; PMC supplementary downloads remain gated; the "CPTAC Data Use Agreement" page cited in the
CPTAC addendum redirects to an unrelated NCI programme page.

**Statuses used:** `RESOLVED — PERMITTED`, `RESOLVED — PROHIBITED`, `RESOLVED — NOT REQUIRED FOR SOURCE SELECTION`, `UNRESOLVED — REQUIRES FURTHER AUTHORIZATION/REVIEW`.
For UV-18 the overlap classification requested (`OVERLAP DOCUMENTED` / `NO OVERLAP DOCUMENTED — INDEPENDENCE NOT PROVEN` / `OVERLAP STATUS UNRESOLVED`) is given as well.

## Authoritative documents read (new in this step)
* **dbGaP Data Use Certification Agreement, Ver. January 25, 2025, for phs000178.v11.p8 (TCGA)** and **for phs001287.v23.p7 (CPTAC)**, each with its study-specific Addendum — public PDFs served by dbGaP
  (`dbgap.ncbi.nlm.nih.gov/aa/wga.cgi?view_pdf=&stacc=<accession>`). Text extracted locally; quoted terms below are from these documents. (These supersede the 2014 TCGA DUC read in v4.)
* Sanger Cancer Dependency Map Data Usage Policy (v4); Figshare API license metadata for the Broad DepMap Public releases; CGC Knowledge Center CPTAC page; GDC CPTAC page; retired CPTAC portal FAQ; the CPTAC papers already read (Dou 2020, Vasaikar 2019).

---

## UV-16 — TCGA data-use terms (evidence ID T-09)
**Status: UNRESOLVED — REQUIRES FURTHER AUTHORIZATION/REVIEW.**

| Aspect | What the authoritative text says (DUC Ver. 25 Jan 2025, phs000178.v11.p8) |
|---|---|
| Access tier | Controlled-access via dbGaP Authorized Access; DAC = NCI DAC. (Open tier for clinical/somatic MAF per v4; BAM/VCF controlled per GDC.) |
| Consent group / use limitation | Consent group 1 "General Research Use (GRU)": "Use of the data is limited only by the terms of the model Data Use Certification." |
| Approved-project limitation | "Research use will occur solely in connection with the approved research project described in the DAR, which includes a 1-2 paragraph description of the proposed research (i.e., a Research Use Statement)"; new uses need a new DAR; access "granted for a period of one (1) year", renewable, with annual Progress Update; scientific collaborators at another institution must submit their own DAR; contractors must be described. |
| Requester / PI eligibility | Requester = "home institution or organization" acting through an Institutional Signing Official; PI = "an investigator who is a permanent employee of their institution at a level equivalent to a tenure-track professor or senior scientist". For-profit status is **not addressed**. |
| Publication / data-sharing | Publication encouraged; acknowledgement with the dbGaP accession required; "Public Posting of Genomic Summary Results — Not Allowed" for TCGA; data and Data Derivatives may not be redistributed without NIH approval and "may not be sold to any individual at any point in time for any purpose"; destruction/close-out at end; security best practices. |
| Commercial-use language | §9 Intellectual Property: NIH "considers these data as pre-competitive and urges Approved Users to avoid making IP claims derived directly from the dataset(s)"; "these NIH-provided data, and conclusions derived therefrom, will remain freely available, without requirement for licensing"; "there is no restriction on development of commercial products resulting from the knowledge gained from the research project". |
| Diagnostic-development language | **None.** The words diagnostic, clinical use, IVD or regulatory submission do not appear as permitted or prohibited uses. |

**Finding.** Commercial product development is not prohibited by the DUC and the consent group carries no additional limitation, but permission for GenRichi's intended use (developing and validating a proprietary MSI panel for a diagnostic product) is **not stated**: it depends on
a DAC-approved Research Use Statement and on how "approved research project" is read against a company's product-validation work, and the NIH expectation that data and conclusions stay freely available and that IP claims directly derived from the data be avoided may bear on a proprietary panel. Permission is not inferred.
**What is needed to resolve:** (1) an NCI DAC decision on a submitted DAR whose Research Use Statement describes the intended use; (2) legal review of DUC terms 1, 6 and 9 against the intended product and any regulatory submission; (3) confirmation that the company can act as Requester (Signing Official / eRA Commons) with an eligible PI.
**Consequence for Source Selection Review:** TCGA cannot be declared eligible; it stays a candidate with the use-permission gate open.

## UV-15 — CPTAC data-use terms (evidence ID C-07)
**Status: UNRESOLVED — REQUIRES FURTHER AUTHORIZATION/REVIEW.**

| Aspect | What the authoritative text says (DUC Ver. 25 Jan 2025, phs001287.v23.p7, "CPTAC Proteogenomic Study") |
|---|---|
| Access tier | Controlled-access via dbGaP; DAC = NCI DAC; the GDC lists CPTAC genomic data as open or controlled (per CGC and GDC pages). CPTAC-3 open RNA-seq quantification is open on the GDC/AWS. **The applicable route for CPTAC colon raw sequence is unclear:** the colon paper cites SRA BioProject PRJNA514017 (which states no public data linked), whereas the endometrial paper cites dbGaP phs001287. |
| Consent group / use limitation | Consent group 1 "General Research Use (GRU)": "Use of the data is limited only by the terms of the model Data Use Certification." |
| Approved-project limitation | Same DUC terms as TCGA (research use solely within the approved DAR; one-year renewable access; collaborator DARs; destruction at close-out). |
| CPTAC-specific term | "IC Specific Access Term": approved users must "recognize any restrictions on data use established by the submitting institution through the Institutional Certification and stated on the CPTAC Data Use Agreement page" (`proteomics.cancer.gov/data-portal/about/data-use-agreement`) — **that page now redirects to an unrelated NCI overview, so the referenced restrictions could not be read.** |
| Publication | Public posting of genomic summary results "Allowed"; freedom-to-publish criteria: a global-analysis publication exists, or 22 months after all samples were received at the data production site, or CPTAC Steering Committee approval; required acknowledgement text and citation of accession phs001287 with retrieval date. (The retired CPTAC portal FAQ states a 15-month variant, so the rule has changed over time; the applicable status per dataset is on the study page's "Available Data" tab, which I did not read.) |
| Commercial-use language | The DUC's general §9 text (as above) applies; the retired CPTAC FAQ does not state commercial or diagnostic-development terms. |
| Diagnostic-development language | **None found.** |

**Finding.** Same position as TCGA plus an unread CPTAC-specific restriction document and an unclear colon access route. Permission is not inferred.
**What is needed to resolve:** the current CPTAC Data Use Agreement text from NCI (the cited URL no longer serves it); the actual access tier/route for the CPTAC colon raw sequences; a DAC decision on a DAR describing the intended use; legal review of the DUC and the CPTAC term against the intended product.
**Consequence for Source Selection Review:** CPTAC cannot be declared eligible; note also v5 UV-05 (CPTAC labels not usable as ground truth at this stage), which is unaffected.

## UV-18 — TCGA–CPTAC overlap (evidence ID X-02)
**Overlap classification: `NO OVERLAP DOCUMENTED — INDEPENDENCE NOT PROVEN`** for the cohorts under consideration (CPTAC prospective colon, Vasaikar 2019; CPTAC endometrial, Dou 2020, against TCGA COAD/READ/UCEC).
**Gate status: UNRESOLVED — REQUIRES FURTHER AUTHORIZATION/REVIEW.**

| Finding | Evidence | Status |
|---|---|---|
| **Overlap is documented for a different CPTAC set.** "CPTAC-2 Retrospective" consists of TCGA tumours: the CGC page states its MS data "consists of four TCGA cancer types (… Colon adenocarcinoma, Rectum adenocarcinoma)"; Zhang et al. 2014 analysed 95 TCGA tumour samples from 90 patients (64 COAD, 31 READ) (search-result text of the paper and data record). | CGC Knowledge Center CPTAC page; Zhang 2014 (Nature 513:382) record | `OVERLAP DOCUMENTED` for CPTAC-2 Retrospective vs TCGA COAD/READ |
| **Prospective cohorts:** the papers describe "a prospective cohort of 95 EC patients" and "prospective colon tumor samples" and use TCGA cohorts only as external comparators (correlation of average profiles; comparison of mutation spectra). Neither states that patients do not overlap. | Dou 2020 (Introduction, Results); Vasaikar 2019 (Fig. S1, Results), texts read in v2/v4 | `NO OVERLAP DOCUMENTED — INDEPENDENCE NOT PROVEN` |
| Public metadata: CPTAC prospective tables use identifier formats different from TCGA barcodes and contain no TCGA identifier (v4); that neither proves nor excludes overlap. A 2026 secondary paper says it cross-checked identifiers, without publishing the check. | v4; EMMA-STRAT | not sufficient |
| GDC CPTAC page and CGC page: no statement of overlap or relationship. | GDC, CGC pages | — |

**Finding.** Independence between TCGA and the CPTAC prospective cohorts is not proven; a real overlap exists between TCGA colorectal and the CPTAC-2 Retrospective set. Different names or repositories are not treated as independence.
**What is needed to resolve:** a written provenance statement from NCI/CPTAC on patient overlap between CPTAC-2 (prospective) / CPTAC-3 cohorts and TCGA, or a patient-level cross-check that requires authorised access (e.g. comparing subject identifiers or germline genotypes), which is not permitted at this stage.
**Consequence for future cohort partitioning (no partition is constructed now):** (a) any CPTAC-2 Retrospective colorectal data must be treated as the same patients as TCGA COAD/READ and kept in one partition; (b) TCGA and the CPTAC prospective cohorts cannot be assigned opposite roles (development vs independent validation) on the strength of independence, because independence is unproven; (c) any assignment must first pass a patient-level check once authorised access exists.
**Consequence for Source Selection Review:** no TCGA–CPTAC independence claim may be made.

## UV-24 — Cell-line data-use terms (evidence ID CL-02)
**Status: UNRESOLVED — REQUIRES FURTHER AUTHORIZATION/REVIEW.**

| Resource | What the authoritative source says | Internal research | Commercial research | Diagnostic development |
|---|---|---|---|---|
| **Sanger Cancer Dependency Map / Cell Model Passports data** (Sanger DepMap Data Usage Policy) | "non-exclusive, non-transferable right to use data files for internal proprietary research and educational purposes, including target, biomarker and drug discovery"; excluded: "resale … or for provision of commercial services"; "Commercial use of the Cancer Dependency Map at Sanger API and/or data is not permitted without prior consent." The Cell Model Passports "about" page states no terms and the policy does not enumerate covered datasets. | Explicitly permitted | **Restricted — prior consent required** | **Not addressed**; product incorporation and commercial services are excluded without consent |
| **Broad DepMap Public release files** (Figshare records for 23Q2, 23Q4 and 24Q2, published by Broad DepMap) | Repository licence metadata: **CC BY 4.0**. Later releases were not checked. | Permitted by the licence | Permitted by the licence (attribution) | Not addressed in the licence |
| **Broad DepMap / CCLE portal terms** | Not readable (Cloudflare challenge). A moderator forum post (v4) says most CCLE files in the public release use CC BY 4.0 but some are unclear, and quotes language that users are responsible for third-party rights. A search snippet asserting a non-commercial restriction is **unverified and conflicts** with the Figshare licence; it is not adopted and not dismissed. | — | **Unresolved** | **Unresolved** |
| **CCLE raw sequence data** (e.g. AWS registry entry "by downloading this data you agree to our Terms and Conditions") | Terms not established from readable text | — | Unresolved | Unresolved |

**Finding.** The only explicit permission found is the CC BY 4.0 licence on the Figshare DepMap Public releases, but the portal's own terms, the raw-sequence terms and the CC BY/portal conflict are unresolved, and the Sanger policy explicitly restricts commercial use without consent. "Internal proprietary research" is not read as permission for diagnostic development.
**What is needed to resolve:** (1) Sanger's written response to a consent request (`depmap@sanger.ac.uk`, per the policy) covering the intended use; (2) the Broad DepMap portal Terms and CCLE Terms as shown in a browser (which I cannot read) and confirmation of which apply to the files actually contemplated (labels/annotations versus raw sequence); (3) legal review of any conflict between the portal terms and the CC BY 4.0 licence.
**Consequence for Source Selection Review:** cell-line resources cannot be declared eligible for any commercial or diagnostic-development purpose; if the contemplated evidence were limited to DepMap Public release annotations, the CC BY 4.0 licence is the only explicit basis and it remains subject to the unread portal terms.

---

## Gate summary
| Gate | Evidence ID | Status | Key reason |
|---|---|---|---|
| UV-15 CPTAC data-use terms | C-07 | **UNRESOLVED — REQUIRES FURTHER AUTHORIZATION/REVIEW** | DUC has no diagnostic-development language; CPTAC-specific restriction document unreadable; colon access route unclear |
| UV-16 TCGA data-use terms | T-09 | **UNRESOLVED — REQUIRES FURTHER AUTHORIZATION/REVIEW** | Commercial products not restricted, but intended use is not stated as permitted; depends on DAC-approved use and legal review |
| UV-18 TCGA–CPTAC overlap | X-02 | **UNRESOLVED — REQUIRES FURTHER AUTHORIZATION/REVIEW** (classification: NO OVERLAP DOCUMENTED — INDEPENDENCE NOT PROVEN; OVERLAP DOCUMENTED for CPTAC-2 Retrospective) | Independence cannot be shown without authorised sample-level access or an NCI/CPTAC statement |
| UV-24 cell-line data-use terms | CL-02 | **UNRESOLVED — REQUIRES FURTHER AUTHORIZATION/REVIEW** | Sanger requires consent for commercial use; Broad portal terms unreadable; only explicit basis is the CC BY 4.0 Figshare licence |

None of the four gates is RESOLVED — PERMITTED, RESOLVED — PROHIBITED or RESOLVED — NOT REQUIRED FOR SOURCE SELECTION.

Source Selection Status: BLOCKED
