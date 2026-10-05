# Stage 6 — Governance Resolution Tracker (tracker construction only)

**STAGE 6 GOVERNANCE GATE: NOT CLOSED**
**SOURCE SELECTION STATUS: NOT STARTED — NO VALIDATION SOURCE SELECTED**

Canonical source of blockers: `STAGE6_GOVERNANCE_GATE_STATUS.md` (16 canonical blockers: A1–A4, B1–B2, C1–C2, D1–D2, E1–E6 in that document's own lettering). Details and evidence are taken only from
`STAGE6_TCGA_GOVERNANCE_RESOLUTION.md`, `STAGE6_CPTAC_GOVERNANCE_RESOLUTION.md`, `STAGE6_CELLLINE_GOVERNANCE_RESOLUTION.md`, `STAGE6_TCGA_CPTAC_INDEPENDENCE_RESOLUTION.md` and `STAGE6_GOVERNANCE_AUTHORIZATION_PACKAGE.md`.
No existing file was modified. No web search, provider contact, email, DAC request, data request, login, or access to controlled data occurred. No source or cohort selected; no legal conclusion is written here; nothing is marked approved.
Not legal advice. Created 2026-09-25.

## 0. How to read this tracker
* Tracker rows are `TR-<domain><n>`. Each maps to one or more **canonical IDs** of the gate-status document (`GS-…`, written here as `GS-A1` etc. to avoid confusion with the tracker's domain letters). Two canonical blockers (GS-E3, GS-E4) span two providers or a provider plus a legal question, so each is split into two rows; there are therefore **18 rows for 16 canonical blockers**, plus **1 row explicitly marked non-canonical** (TR-E6, added to cover the data-protection topic you listed; its content comes from obligations already documented in the DUC terms and the Source Selection Review specification).
* Tracker domains: **A** TCGA/NCI/dbGaP · **B** CPTAC/NCI/dbGaP/SRA · **C** TCGA–CPTAC independence · **D** Cell lines (Sanger and Broad/DepMap kept apart) · **E** Legal/compliance (questions only; no conclusions).
* Status values used: `OPEN`, `WAITING FOR PROVIDER`, `WAITING FOR DAC`, `WAITING FOR LEGAL/COMPLIANCE`, `WAITING FOR REQUESTER ELIGIBILITY`. `RESOLVED`, `REJECTED` and `NOT APPLICABLE` are not used: no authoritative closing evidence exists and no blocker has been rejected or found inapplicable.
* "Internal owner" is `UNASSIGNED` for every row: the governance documents name roles (product owner; GenRichi legal/compliance function) but no individuals, and none is invented here.
* "Date opened" is the date the blocker was documented (2026-09-25); "Date resolved" is empty until an authoritative closing record exists.
* No DAR has been drafted, submitted or requested for any row that mentions a DAC; "WAITING FOR DAC" means the decision belongs to the DAC, not that a request is pending.

## 1. Summary
| Row | Canonical ID | Domain | Topic | Status |
|---|---|---|---|---|
| TR-A1 | GS-A1 | A TCGA | DAC decision on intended use (phs000178) | WAITING FOR DAC |
| TR-A2 | GS-C1 | A TCGA (also applies to B) | Provider position on commercial/for-profit requesters | WAITING FOR PROVIDER |
| TR-A3 | GS-C2 | A TCGA (also applies to B) | GenRichi requester-eligibility facts | WAITING FOR REQUESTER ELIGIBILITY |
| TR-B1 | GS-A3 | B CPTAC | Current CPTAC Data Use Agreement text | WAITING FOR PROVIDER |
| TR-B2 | GS-A4 | B CPTAC | Colon accession, access route, tier and terms | WAITING FOR PROVIDER |
| TR-B3 | GS-A2 | B CPTAC | DAC decision on intended use (applicable accession(s)) | WAITING FOR DAC |
| TR-C1 | GS-D1 | C independence | Provenance statement or authorized identity cross-check | WAITING FOR PROVIDER |
| TR-C2 | GS-D2 | C independence | Cohort↔accession mapping (endometrial 95 vs 241; colon via TR-B2) | WAITING FOR PROVIDER |
| TR-D1 | GS-E3 (Sanger consent part) | D Sanger | Prior consent: commercial use, product incorporation, diagnostic use | WAITING FOR PROVIDER |
| TR-D2 | GS-E4 (Sanger part) | D Sanger | Raw exome access terms (EGA DAC / DAA) | WAITING FOR PROVIDER |
| TR-D3 | GS-E2 | D Broad/DepMap | Portal terms, CCLE terms, release README, CC BY 4.0 scope | OPEN |
| TR-D4 | GS-E4 (Broad part) | D Broad/DepMap | CCLE raw-read terms (SRA PRJNA523380) | OPEN |
| TR-D5 | GS-E1 | D cell lines (internal control) | Role limited to development/control | OPEN |
| TR-E1 | GS-B1 | E legal | DUC Terms 1, 6, 9; derived outputs; diagnostic and commercial context (TCGA and CPTAC dbGaP agreements) | WAITING FOR LEGAL/COMPLIANCE |
| TR-E2 | GS-B2 | E legal | CPTAC-specific term review (once its text exists) | OPEN |
| TR-E3 | GS-E3 (legal part) | E legal | Sanger "internal proprietary research" determination | WAITING FOR LEGAL/COMPLIANCE |
| TR-E4 | GS-E5 | E legal | Diagnostic-development use of cell-line resources (per resource) | WAITING FOR LEGAL/COMPLIANCE |
| TR-E5 | GS-E6 | E legal | Third-party rights before any use beyond internal research | WAITING FOR LEGAL/COMPLIANCE |
| TR-E6 | *(non-canonical)* | E legal | Data-protection obligations after any authorization | OPEN |

**Counts:** 16 canonical blockers → 18 tracker rows (+1 non-canonical). By tracker domain (each canonical blocker counted once by its primary row): **A 3 · B 3 · C 2 · D 4 · E 4** (= 16).
Status counts over all 19 rows: `OPEN` 5 · `WAITING FOR PROVIDER` 7 · `WAITING FOR DAC` 2 · `WAITING FOR LEGAL/COMPLIANCE` 4 · `WAITING FOR REQUESTER ELIGIBILITY` 1 · `RESOLVED` 0 · `REJECTED` 0 · `NOT APPLICABLE` 0.

## 2. Topic coverage check (every topic you listed → row)
| Domain | Topic | Row(s) |
|---|---|---|
| A | intended-use compatibility; DAC requirement/decision | TR-A1 |
| A | diagnostic-product development; commercial/for-profit interpretation; derived data/genomic summary restrictions | TR-E1 (legal question); requester side TR-A2, TR-A3 |
| A | requester eligibility | TR-A2, TR-A3 |
| B | current CPTAC-specific terms | TR-B1 (text), TR-E2 (review) |
| B | colon accession/access route; access tier | TR-B2 |
| B | intended-use compatibility; DAC requirement/decision | TR-B3 |
| B | requester eligibility | TR-A2, TR-A3 (shared) |
| C | prospective colon provenance; prospective endometrial provenance; patient-level overlap; independence determination | TR-C1 |
| C | cohort-to-accession mapping | TR-C2 (endometrial), TR-B2 (colon) |
| D Sanger | commercial use; diagnostic development; product incorporation; prior consent | TR-D1 (consent), TR-E3, TR-E4 (legal) |
| D Sanger | raw exome access | TR-D2 |
| D Broad | exact applicable licence; portal terms; release README; CC BY 4.0 scope; commercial/diagnostic use | TR-D3, TR-E4 |
| D Broad | CCLE raw-read terms | TR-D4 |
| D Broad | third-party restrictions | TR-E5 |
| E | DUA compatibility; diagnostic development; commercial context; derived outputs | TR-E1 |
| E | third-party rights | TR-E5 |
| E | data protection after authorization | TR-E6 (non-canonical) |
Gap noted: "data protection after authorization" has no canonical blocker in the gate-status document; it is carried only as the clearly marked non-canonical TR-E6.

## 3. Blocker cards
Common fields for every row: **Internal owner:** UNASSIGNED · **Date opened:** 2026-09-25 · **Date resolved:** — .

### Domain A — TCGA / NCI / dbGaP
**TR-A1 — canonical GS-A1**
* **Domain / provider:** A · NCI DAC (dbGaP `phs000178`)
* **Question:** Is the described GenRichi use — TCGA controlled-access sequence data for Stage 6 empirical MSI evidence for a proprietary diagnostic-product programme — an "approved research project"?
* **Current status:** WAITING FOR DAC
* **Why internal review cannot close it:** the DUC limits use "solely in connection with the approved research project described in the DAR"; approval of that project is the DAC's determination, and no DAR exists.
* **Required authoritative evidence:** a DAC decision on a Data Access Request whose Research Use Statement describes the intended use (dbGaP authorization record).
* **Decision-maker:** NCI DAC
* **Resolution condition:** an authoritative DAC decision on the described use is recorded.
* **Evidence location:** `STAGE6_TCGA_GOVERNANCE_RESOLUTION.md` §2, §4 (rows 2 and 4), §5 (B1); `STAGE6_GOVERNANCE_GATE_STATUS.md` §3 A1
* **Notes:** prerequisite: an agreed internal one-sentence use description (see TR-E1). No DAR has been prepared or submitted.

**TR-A2 — canonical GS-C1**
* **Domain / provider:** A (also applies to B) · NIH / dbGaP
* **Question:** May a commercial/for-profit organization act as a Requester under the DUC? The agreements are silent.
* **Current status:** WAITING FOR PROVIDER
* **Why internal review cannot close it:** eligibility rules are the provider's; the DUC text read defines the Requester as a "home institution or organization" and does not address for-profit status; the NIH pages that might state it were unreadable.
* **Required authoritative evidence:** a written NIH/dbGaP statement (or readable authoritative policy text) on commercial/for-profit requester eligibility.
* **Decision-maker:** NIH / dbGaP (provider); DAC in individual cases
* **Resolution condition:** the provider's position is recorded from an authoritative source and applies to the relevant accessions.
* **Evidence location:** `STAGE6_TCGA_GOVERNANCE_RESOLUTION.md` §2, §4 row 6, §5 (B3); `STAGE6_CPTAC_GOVERNANCE_RESOLUTION.md` §4
* **Notes:** shared with the CPTAC route; silence is not treated as permission.

**TR-A3 — canonical GS-C2**
* **Domain / provider:** A (also applies to B) · GenRichi (facts) checked against the DUC definitions
* **Question:** Does GenRichi have an Institutional Signing Official, an eRA Commons registration and a PI meeting the DUC definition ("a permanent employee … at a level equivalent to a tenure-track professor or senior scientist")?
* **Current status:** WAITING FOR REQUESTER ELIGIBILITY
* **Why internal review cannot close it:** these are facts about GenRichi not documented in the governance material, and the provider's rule (TR-A2) is not yet known.
* **Required authoritative evidence:** documented GenRichi eligibility facts (registrations and roles) assessed against the provider's stated criteria.
* **Decision-maker:** GenRichi legal/compliance for the facts; provider for the criteria
* **Resolution condition:** the facts are documented and shown to satisfy the criteria confirmed under TR-A2.
* **Evidence location:** `STAGE6_TCGA_GOVERNANCE_RESOLUTION.md` §2, §4 row 6; `STAGE6_GOVERNANCE_GATE_STATUS.md` §3 C2
* **Notes:** shared with the CPTAC route.

### Domain B — CPTAC / NCI / dbGaP / SRA
**TR-B1 — canonical GS-A3**
* **Domain / provider:** B · NCI / CPTAC
* **Question:** What is the current text of the "CPTAC Data Use Agreement" cited by the `phs001287` addendum, and the dataset-specific publication status?
* **Current status:** WAITING FOR PROVIDER
* **Why internal review cannot close it:** the cited URL now redirects to an unrelated NCI page, so the term is unread.
* **Required authoritative evidence:** the current agreement text from NCI/CPTAC and the study page's dataset status.
* **Decision-maker / provider:** NCI / CPTAC (Office of Cancer Clinical Proteomics Research)
* **Resolution condition:** the current authoritative text is obtained and recorded.
* **Evidence location:** `STAGE6_CPTAC_GOVERNANCE_RESOLUTION.md` §2A, §3, §5 (C1)
* **Notes:** also records the 22-month (addendum) versus 15-month (retired FAQ) publication-rule discrepancy.

**TR-B2 — canonical GS-A4**
* **Domain / provider:** B · NCI / CPTAC (dbGaP, GDC, SRA)
* **Question:** Which accession/route governs the CPTAC colon cohort's exome data — `phs001287`, `phs000892` or SRA `PRJNA514017` — and what are its access tier and terms?
* **Current status:** WAITING FOR PROVIDER
* **Why internal review cannot close it:** the colon paper cites SRA; GDC/dbGaP hold a colon project under `phs000892`; no authoritative source maps the 110-patient cohort to an accession; the SRA record has no public data and no stated tier.
* **Required authoritative evidence:** provider confirmation of accession, repository, tier and applicable agreement for the colon cohort.
* **Decision-maker / provider:** NCI / CPTAC
* **Resolution condition:** cohort ↔ accession ↔ tier ↔ agreement recorded from an authoritative source.
* **Evidence location:** `STAGE6_CPTAC_GOVERNANCE_RESOLUTION.md` §2A–C, §4 (added row), §5 (C2)
* **Notes:** supplies the colon half of TR-C2; the 110/106/105+1 differences are preserved, not resolved.

**TR-B3 — canonical GS-A2**
* **Domain / provider:** B · NCI DAC (applicable accession(s) after TR-B2)
* **Question:** Is the described use an approved research project for the applicable CPTAC accession(s)?
* **Current status:** WAITING FOR DAC
* **Why internal review cannot close it:** same DUC Term 1 mechanism as TR-A1; the applicable accession is not yet established.
* **Required authoritative evidence:** a DAC decision on a DAR for the applicable accession(s).
* **Decision-maker:** NCI DAC
* **Resolution condition:** authoritative DAC decision recorded for each applicable accession.
* **Evidence location:** `STAGE6_CPTAC_GOVERNANCE_RESOLUTION.md` §4, §5 (C3)
* **Notes:** depends on TR-B2; no DAR prepared or submitted.

### Domain C — TCGA–CPTAC independence
**TR-C1 — canonical GS-D1**
* **Domain / provider:** C · NCI / CPTAC / Biospecimen Core Resource (or an authorized Approved User)
* **Question:** Is any patient in the CPTAC prospective colon cohort or the CPTAC endometrial cohort also a TCGA participant? (Provenance, patient-level overlap and the independence determination all follow from this.)
* **Current status:** WAITING FOR PROVIDER
* **Why internal review cannot close it:** public documents establish overlap only for CPTAC-2 Retrospective; for the prospective cohorts neither overlap nor independence is stated; identifier formats differ, which proves neither.
* **Required authoritative evidence:** a written provenance statement per cohort and accession, or an authorized patient-level identity cross-check under approved access.
* **Decision-maker / provider:** NCI / CPTAC / Biospecimen Core Resource; or an Approved User under DAC-approved access
* **Resolution condition:** overlap status for each cohort is documented by one of those routes and sufficient for the intended role.
* **Evidence location:** `STAGE6_TCGA_CPTAC_INDEPENDENCE_RESOLUTION.md` §2–§5, §7
* **Notes:** independence must not be claimed until then; CPTAC-2 Retrospective colorectal = TCGA patients (established).

**TR-C2 — canonical GS-D2**
* **Domain / provider:** C · NCI / CPTAC
* **Question:** How do the cohorts map to accessions and cases — endometrial: the paper's 95 tumours versus 241 GDC "Uterus, NOS" cases; colon: see TR-B2?
* **Current status:** WAITING FOR PROVIDER
* **Why internal review cannot close it:** the relationship is undocumented in public sources read.
* **Required authoritative evidence:** provider documentation defining each cohort as a set of cases within accessions.
* **Decision-maker / provider:** NCI / CPTAC
* **Resolution condition:** each cohort is defined as a documented set of cases and accessions.
* **Evidence location:** `STAGE6_TCGA_CPTAC_INDEPENDENCE_RESOLUTION.md` §4, §5 (item 3)
* **Notes:** the CPTAC endometrial count discrepancy (25 MSI-H of 95 primary/public table vs 21 in a 2026 secondary source) is preserved, not resolved.

### Domain D — Cell lines
**Sanger**

**TR-D1 — canonical GS-E3 (consent part)**
* **Domain / provider:** D Sanger · Wellcome Sanger Institute (contact named in its policy)
* **Question:** Does Sanger consent to commercial use, product incorporation, or diagnostic-product use of its Cancer Dependency Map/Cell Model Passports data? The policy states commercial use "is not permitted without prior consent"; diagnostic development is not addressed.
* **Current status:** WAITING FOR PROVIDER
* **Why internal review cannot close it:** consent can only be granted by Sanger; the policy does not enumerate which Cell Model Passports files it covers.
* **Required authoritative evidence:** Sanger's written position or consent covering the intended use.
* **Decision-maker / provider:** Wellcome Sanger Institute
* **Resolution condition:** written consent or a written refusal is recorded.
* **Evidence location:** `STAGE6_CELLLINE_GOVERNANCE_RESOLUTION.md` §2A, §5 (condition 3)
* **Notes:** not contacted; no consent assumed.

**TR-D2 — canonical GS-E4 (Sanger part)**
* **Domain / provider:** D Sanger · Sanger DAC (WTSI CGP, EGAC00001000000) via EGA
* **Question:** What are the Data Access Agreement and any commercial-use terms for the Sanger exome dataset EGAD00001001039? (Needed only if raw data are ever contemplated.)
* **Current status:** WAITING FOR PROVIDER
* **Why internal review cannot close it:** the dataset is DAC-controlled with general-research-use, publication-required, approved-user and approved-institution modifiers; the agreement text and the "Data Sharing Policy" were not read.
* **Required authoritative evidence:** the DAA and policy text from the DAC/EGA.
* **Decision-maker / provider:** Sanger DAC
* **Resolution condition:** terms recorded (or the raw route formally excluded from Stage 6).
* **Evidence location:** `STAGE6_CELLLINE_GOVERNANCE_RESOLUTION.md` §2B, §5 (condition 4)
* **Notes:** conditional relevance.

**Broad / DepMap**

**TR-D3 — canonical GS-E2**
* **Domain / provider:** D Broad/DepMap · Broad DepMap portal, CCLE terms, Figshare release records
* **Question:** Do the DepMap portal Terms, the CCLE terms and the release `README.txt` add restrictions to the CC BY 4.0 licence recorded for DepMap 23Q2, 23Q4 and 24Q2, for the specific files contemplated (exact applicable licence, portal vs downloadable terms, CC BY 4.0 scope, commercial/diagnostic use)?
* **Current status:** OPEN
* **Why internal review cannot close it:** the portal terms are behind a Cloudflare challenge my tools could not pass (not bypassed) and the README was not downloaded; a search snippet claiming a non-commercial restriction conflicts with the licence and is unverified.
* **Required authoritative evidence:** the terms and README as shown to a person with browser access, and confirmation of which apply to the contemplated files; a provider clarification if ambiguous.
* **Decision-maker / provider:** Broad Institute / DepMap (terms author); GenRichi reads and records
* **Resolution condition:** the applicable licence and terms are recorded for the specific files, with any conflict identified.
* **Evidence location:** `STAGE6_CELLLINE_GOVERNANCE_RESOLUTION.md` §3A–B, §5 (condition 2)
* **Notes:** later releases than 24Q2 were not checked.

**TR-D4 — canonical GS-E4 (Broad part)**
* **Domain / provider:** D Broad/DepMap · Broad (submitter) / NCBI SRA `PRJNA523380`
* **Question:** Which terms govern the CCLE raw reads (publicly indexed, 4,550 SRA records)? NCBI places no restrictions itself, but submitters may claim rights and Broad's own terms are not stated. (Needed only if raw data are ever contemplated.)
* **Current status:** OPEN
* **Why internal review cannot close it:** the submitter's terms were not found in readable authoritative text.
* **Required authoritative evidence:** Broad/submitter terms for the raw reads.
* **Decision-maker / provider:** Broad Institute
* **Resolution condition:** terms recorded (or the raw route formally excluded from Stage 6).
* **Evidence location:** `STAGE6_CELLLINE_GOVERNANCE_RESOLUTION.md` §3C, §5 (condition 4)
* **Notes:** conditional relevance.

**Internal control (applies to both resources)**

**TR-D5 — canonical GS-E1**
* **Domain / provider:** D cell lines · GenRichi (internal control)
* **Question:** Is cell-line material confined to a development/control role and never used as patient-level ground truth?
* **Current status:** OPEN
* **Why internal review cannot close it:** it is an ongoing condition on future use, verified at Source Selection Review and at every later use, not a one-off fact; it is recorded in the Source Selection Review specification but nothing external closes it.
* **Required authoritative evidence:** the recorded role assignment in the future source-selection record.
* **Decision-maker / provider:** GenRichi product owner
* **Resolution condition:** the role restriction is recorded in the source-selection record for any cell-line use.
* **Evidence location:** `STAGE6_CELLLINE_GOVERNANCE_RESOLUTION.md` §1, §5 (condition 1); `STAGE6_SOURCE_SELECTION_REVIEW.md` §1
* **Notes:** cell-line MSI labels are computational (Sanger MSIsensor-pro score ≥7; CCLE deletion-count classification).

### Domain E — Legal / compliance (questions only; no conclusions)
**TR-E1 — canonical GS-B1**
* **Domain / provider:** E · GenRichi legal/compliance (qualified reviewer)
* **Question:** How do DUC Term 1 (research use within the approved DAR), Term 6 (no sale of Data Derivatives) and Term 9 (IP; "avoid making IP claims derived directly from the dataset(s)"; data and conclusions "freely available"; "no restriction on development of commercial products resulting from the knowledge gained from the research project") apply to a proprietary diagnostic-panel product and any regulatory submission, for the TCGA and CPTAC dbGaP agreements? Are Stage 6 evidence tables "Data Derivatives" or "genomic summary results" (posting "Not Allowed" in the TCGA addendum, "Allowed" in the CPTAC addenda)?
* **Current status:** WAITING FOR LEGAL/COMPLIANCE
* **Why internal review cannot close it:** the agreements contain no diagnostic/clinical/regulatory language; reading them beyond their text is a legal interpretation that documentation review cannot supply, and the commercial-product sentence is not treated as diagnostic-development permission.
* **Required authoritative evidence:** a written determination by qualified legal/compliance reviewers on these questions.
* **Decision-maker / provider:** qualified legal/compliance reviewer
* **Resolution condition:** a written determination is recorded for each question.
* **Evidence location:** `STAGE6_TCGA_GOVERNANCE_RESOLUTION.md` §2, §4 row 3, §5 (B2); `STAGE6_CPTAC_GOVERNANCE_RESOLUTION.md` §2A, §4; `STAGE6_GOVERNANCE_GATE_STATUS.md` §3 B1
* **Notes:** includes the internal use-description sentence needed for TR-A1 and TR-B3.

**TR-E2 — canonical GS-B2**
* **Domain / provider:** E · GenRichi legal/compliance
* **Question:** How does the CPTAC-specific term relate to the intended product once its text is available?
* **Current status:** OPEN
* **Why internal review cannot close it:** the term's text is unavailable (TR-B1); nothing can be reviewed.
* **Required authoritative evidence:** the text from TR-B1 and a legal/compliance review of it.
* **Decision-maker / provider:** qualified legal/compliance reviewer
* **Resolution condition:** review completed and recorded after TR-B1 closes.
* **Evidence location:** `STAGE6_CPTAC_GOVERNANCE_RESOLUTION.md` §2A, §5 (C4)
* **Notes:** blocked by TR-B1.

**TR-E3 — canonical GS-E3 (legal part)**
* **Domain / provider:** E · GenRichi legal/compliance
* **Question:** Does the intended Stage 6 use fall within the Sanger policy's "internal proprietary research and educational purposes", and not "commercial use"?
* **Current status:** WAITING FOR LEGAL/COMPLIANCE
* **Why internal review cannot close it:** the policy does not define the boundary between internal proprietary research and commercial use.
* **Required authoritative evidence:** a written legal/compliance determination.
* **Decision-maker / provider:** qualified legal/compliance reviewer
* **Resolution condition:** determination recorded; if the use is not within the internal-research language, TR-D1 consent is needed.
* **Evidence location:** `STAGE6_CELLLINE_GOVERNANCE_RESOLUTION.md` §2A, §5 (condition 3)
* **Notes:** pairs with TR-D1.

**TR-E4 — canonical GS-E5**
* **Domain / provider:** E · GenRichi legal/compliance
* **Question:** For each cell-line resource (Sanger, Broad/DepMap), does any term address or permit use for diagnostic-product development, and how do the commercial-use terms apply?
* **Current status:** WAITING FOR LEGAL/COMPLIANCE
* **Why internal review cannot close it:** neither the Sanger policy nor the CC BY 4.0 licence addresses diagnostic development; permission is not assumed.
* **Required authoritative evidence:** a written legal/compliance determination per resource.
* **Decision-maker / provider:** qualified legal/compliance reviewer
* **Resolution condition:** determination recorded per resource.
* **Evidence location:** `STAGE6_CELLLINE_GOVERNANCE_RESOLUTION.md` §4, §5 (condition 5)
* **Notes:** covers "Sanger diagnostic development" and "Broad commercial/diagnostic use".

**TR-E5 — canonical GS-E6**
* **Domain / provider:** E · GenRichi legal/compliance
* **Question:** What third-party rights apply (CC BY 4.0 does not cover them; cell-bank/CCLE rights) before any use beyond internal research?
* **Current status:** WAITING FOR LEGAL/COMPLIANCE
* **Why internal review cannot close it:** the licence text states that third-party rights may limit use; assessing them is a legal question.
* **Required authoritative evidence:** a written legal/compliance review.
* **Decision-maker / provider:** qualified legal/compliance reviewer
* **Resolution condition:** review recorded.
* **Evidence location:** `STAGE6_CELLLINE_GOVERNANCE_RESOLUTION.md` §3A, §5 (condition 6)
* **Notes:** applies to Broad/DepMap release files and the CCLE material.

**TR-E6 — NON-CANONICAL (added for the "data protection after authorization" topic)**
* **Domain / provider:** E · GenRichi legal/compliance and IT/security (roles only)
* **Question:** What storage/processing location and data-protection obligations (security best practices, retention and destruction at close-out, redistribution limits, incident reporting) would apply after any authorization, and can GenRichi meet them?
* **Current status:** OPEN
* **Why internal review cannot close it:** the obligations are conditional on authorization that has not been granted; the DUC lists the obligations but which agreement and version apply is not yet known (TR-B2, TR-A1).
* **Required authoritative evidence:** the applicable agreement text plus a GenRichi compliance plan reviewed by qualified staff.
* **Decision-maker / provider:** qualified legal/compliance and security reviewers; NCI DAC for the agreement terms
* **Resolution condition:** obligations identified from the applicable agreement and a compliant storage/processing arrangement documented.
* **Evidence location:** `STAGE6_TCGA_GOVERNANCE_RESOLUTION.md` §2 (duration, security, destruction); `STAGE6_SOURCE_SELECTION_REVIEW.md` §2D
* **Notes:** **not a canonical blocker of the gate-status document**; listed only because it was requested. It is not counted in the 16.

## 4. Gate logic (unchanged)
STAGE 6 GOVERNANCE GATE: NOT CLOSED — EXTERNAL RESOLUTION REQUIRED
SOURCE SELECTION STATUS: NOT STARTED
NO VALIDATION SOURCE SELECTED
NEXT GATE: GOVERNANCE RESOLUTION → FORMAL SOURCE SELECTION REVIEW
