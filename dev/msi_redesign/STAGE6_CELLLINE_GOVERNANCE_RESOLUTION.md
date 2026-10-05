# Stage 6 — Cell-line Data-use Governance Resolution (UV-24) (governance documentation only)

Question: **Can cell-line resources be used legitimately as a development/control resource for GenRichi Stage 6, and under what conditions?**
Only public documentation and public metadata were read. **No cell-line data was downloaded** (no release file, no README file, no sequence data), no raw sequence data was opened, no cell line was selected, no cohort constructed,
no code, threshold, classifier or production file was touched, no provider (Sanger, Broad, DepMap, EGA, NCBI) was contacted and no access was requested. TCGA and CPTAC governance were not reopened. Nothing here is legal advice.
Sanger and Broad/DepMap are kept completely separate. Date: 2026-09-25.

## 1. Scope
**Use evaluated:** use of cancer cell-line **annotations and, if ever needed, sequencing data** (for example MSI status annotations and exome/genome reads of established cell lines) as a **development/control resource** for Stage 6 — i.e. internal research and development on candidate MSI loci.
Categories kept apart: **research/development use** (evaluated), **commercial use** (evaluated as a term category), **diagnostic-product development** (the purpose behind the programme; evaluated only as far as terms address it), **clinical validation** (not evaluated as a use; excluded).
**Fixed role:** cell-line material is a development/control resource and **not** a substitute for independent patient-level MSI validation. Cell-line MSI labels are computational (Sanger: MSIsensor-pro score ≥7; CCLE: deletion-count classification, per v2), never patient-level clinical ground truth.

## 2. Sanger (Wellcome Sanger Institute: Cancer Dependency Map at Sanger / Cell Model Passports / GDSC / EGA)

**A. Processed data and annotations — "Data Usage Policy" (Cancer Dependency Map at Sanger, version 2.0.0 page; read directly at `https://depmap.sanger.ac.uk/documentation/data-usage-policy/`).**
| Point | Text as written |
|---|---|
| Research-use permission | "Users have a non-exclusive, non-transferable right to use data files for internal proprietary research and educational purposes, including target, biomarker and drug discovery." |
| Excluded uses | "Excluded from this licence are use of the data (in whole or any significant part) for resale either alone or in combination with additional data/product offerings, or for provision of commercial services." |
| Commercial use | "Commercial use of the Cancer Dependency Map at Sanger API and/or data is not permitted without prior consent. Please contact depmap@sanger.ac.uk for enquiries." |
| Product incorporation | "If you are interested in incorporating results or software into a product, or have questions, please contact depmap@sanger.ac.uk." |
| Regulatory / diagnostic | "The data files are experimental and academic in nature and are not licensed or certified by any regulatory body." Diagnostic-product development is not otherwise addressed. |
| Warranty | Provided "on an 'as is' basis and excludes all warranties of any kind." |
| Redistribution | The right is "non-transferable"; direct use of the API in third-party websites "without prior permission is not permitted"; API use is free "for non-commercial use". |
| Which datasets are covered | The policy applies to "data files" in the Cancer Dependency Map at Sanger; it does not enumerate Cell Model Passports datasets. The Sanger DepMap "Websites" page describes Cell Model Passports as part of the Sanger Dependency Map and states that "Links to raw data are provided" and that MSI, ploidy and mutational burden are reported. Whether this policy governs every Cell Model Passports file is **not stated**. The Cell Model Passports web application was not readable by my tools (JavaScript page). |

**B. Raw sequencing data (Sanger).**
| Point | Evidence |
|---|---|
| Repository and tier | The CCLE paper cites Sanger GDSC exome BAMs at EGA accession **EGAD00001001039**. The public EGA dataset page: "Cancer Cell Line Exome Sequencing", 1,072 samples, DAC **EGAC00001000000 (WTSI CGP Data Access Committee)**, "Request Access" required — i.e. controlled access through a DAC (`https://ega-archive.org/datasets/EGAD00001001039`, `https://ega-archive.org/dacs/EGAC00001000000`). |
| Data-use modifiers on the dataset | GRU (general research use: "use is allowed for general research use for any research purpose"), PUB (publication required), US (use limited to approved users), IS (use limited "within an approved institution"). No non-commercial or non-profit modifier is listed. The dataset cites the "Wellcome Trust Sanger Institute Cancer Genome Group Data Sharing Policy" (text not read). |
| Request procedure | EGA states that a request "may require additional approvals or agreements depending on the datasets and data providers involved" and that the requester should also contact the DAC (`datasharing@sanger.ac.uk` per the DAC page). |
| Not established | The Data Access Agreement text, any commercial-use or diagnostic-development term for this dataset, and whether the Sanger Data Usage Policy above also governs these raw files. |

## 3. Broad / DepMap / CCLE

**A. Release files (processed data and annotations).**
| Point | Evidence |
|---|---|
| Licence of specific releases | Figshare records **DepMap 23Q2, DepMap 23Q4 and DepMap 24Q2 Public**, published by "Broad DepMap": licence **CC BY 4.0** (Figshare API metadata; e.g. `https://api.figshare.com/v2/articles/25880521`). Later releases were not found or checked. |
| What CC BY 4.0 states (`https://creativecommons.org/licenses/by/4.0/`) | Reusers may "copy and redistribute the material in any medium or format for any purpose, even commercially" and "remix, transform, and build upon the material for any purpose, even commercially", with attribution; no added legal or technological restrictions on what the licence permits; the licence "does not cover" publicity, privacy or moral rights, gives "No warranties", and patent, trademark and third-party rights may limit use. |
| Release description | The Figshare description refers to `README.txt` for more information; the README (a file inside the dataset) was **not** read, to respect the no-download rule. |
| Diagnostic / product development | Not addressed by the licence or the description. |

**B. DepMap portal and CCLE terms.** The portal terms pages (`depmap.org/portal/terms/`, `/portal/ccle/terms_and_conditions`) are served behind a Cloudflare challenge and are **unavailable to my tools; not bypassed**. The AWS Registry entry for the DepMap/CCLE omics data states "By downloading this data you agree to our Terms and Conditions" (that link is the unreadable portal page). A community-forum thread (secondary) says most CCLE files in the public DepMap release use CC BY 4.0, some are unclear, and quotes third-party-rights language. A search snippet asserting that portal data are "not intended for clinical or commercial uses" was **not verified** from the primary page and conflicts with the Figshare licence; it is neither adopted nor dismissed. Whether portal downloads and the Figshare files are governed by different terms is **not established**.

**C. Raw sequencing data (CCLE).**
| Point | Evidence |
|---|---|
| Repository | The CCLE paper's data-availability statement: "Raw sequencing data are available at Sequence Read Archive (SRA) under accession number **PRJNA523380**" (Ghandi 2019, PMC6697103). |
| Tier | Public SRA index returns **4,550 SRA records** for `PRJNA523380[BioProject]` (public metadata count, queried 2026-09-25), so publicly indexed reads exist. |
| Terms | NCBI policy: "NCBI itself places no restrictions on the use or distribution of the data contained therein"; "Some submitters of the original data … may claim patent, copyright, or other intellectual property rights", which NCBI cannot assess (`https://www.ncbi.nlm.nih.gov/home/about/policies/`). The submitter's (Broad's) own terms for these reads are **not stated** in what I read. |

## 4. Comparison matrix
Status values: **RESOLVED** (explicit text supports the answer), **PARTIALLY RESOLVED**, **UNRESOLVED**. No status is inferred from the existence of a public release.
| Question | Sanger | Broad/DepMap | Status |
|---|---|---|---|
| Research/development use | Processed data: explicitly permitted for "internal proprietary research … including target, biomarker and drug discovery". Raw EGA exomes: GRU/PUB/US/IS via DAC request. | Release files: CC BY 4.0 permits use "for any purpose"; portal terms unread. Raw: public SRA reads; submitter terms not stated. | Sanger: **RESOLVED** (processed) / **UNRESOLVED** (raw). Broad: **PARTIALLY RESOLVED** |
| Commercial use | **Restricted**: "not permitted without prior consent"; resale and commercial services excluded. Raw EGA: no commercial term stated. | CC BY 4.0 explicitly allows commercial use of the release files; portal/CCLE terms unread, and a conflicting unverified snippet exists. | Sanger: **RESOLVED (consent required)**. Broad: **PARTIALLY RESOLVED** |
| Diagnostic-product development | Not addressed; product incorporation requires contact ("incorporating results … into a product"); data "not licensed or certified by any regulatory body". | Not addressed by the licence, README unread, portal terms unread. | **UNRESOLVED** (both) |
| Clinical validation | Not a permitted role for this resource by project rule; Sanger states the data are experimental/academic and not regulator-certified. | Not a permitted role by project rule; no Broad statement read. | **RESOLVED** as excluded by rule (no terms needed) |
| Raw sequencing-data access | EGAD00001001039: controlled via WTSI CGP DAC, request required; DAA unread. | SRA PRJNA523380: publicly indexed (4,550 records); NCBI places no restrictions; submitter claims and Broad terms not stated. | Sanger: **PARTIALLY RESOLVED**. Broad: **PARTIALLY RESOLVED** |
| Redistribution | "Non-transferable" right; API use in third-party sites needs permission; EGA data limited to approved users/institution (US/IS). | CC BY 4.0 permits redistribution with attribution for the release files; portal terms unread. | Sanger: **RESOLVED (restricted)**. Broad: **PARTIALLY RESOLVED** |
| Additional authorization required | **Yes** for commercial use/product incorporation (Sanger consent) and for raw exomes (DAC request). | Not established: none stated by the licence for release files; portal terms unread. | Sanger: **RESOLVED (yes)**. Broad: **UNRESOLVED** |

## 5. Final gate
CELLLINE GOVERNANCE STATUS: CONDITIONAL

The resource may be considered for the fixed role only (development/control, never patient-level ground truth), and only if **all** of the following conditions are met. Nothing is approved, permitted or validated beyond what is quoted above.
1. **Role condition:** cell-line data are used only as a development/control resource; their computational MSI labels are never treated as clinical ground truth.
2. **Broad/DepMap condition (processed release files):** the DepMap portal Terms and Conditions, the CCLE terms and the release `README.txt` — currently unread — are checked and shown not to add restrictions to the CC BY 4.0 licence for the specific files contemplated; the attribution condition and the third-party-rights limits are recorded.
3. **Sanger condition (processed data):** Sanger data are used only within "internal proprietary research and educational purposes"; a legal/compliance determination confirms that the intended Stage 6 use falls within that language and is not "commercial use"; **any** incorporation of Sanger-derived results into a product or commercial service, or any diagnostic-product use, is treated as requiring prior consent from Sanger (`depmap@sanger.ac.uk`, named in the policy) and is **not** covered by this gate.
4. **Raw-data condition:** no raw sequencing data (Sanger EGA or Broad SRA) is used unless its specific terms are obtained and reviewed (EGA DAC/DAA for EGAD00001001039; Broad/submitter terms for PRJNA523380); the raw-data route is otherwise out of scope of this gate.
5. **Diagnostic-development condition:** no diagnostic-product-development use of cell-line resources is assumed to be permitted by any text above; it stays an open legal/compliance question for each resource.
6. **Third-party rights condition:** legal review of the third-party-rights language (CC BY 4.0 rights not covered; CCLE/cell-bank rights) before any use that goes beyond internal research.

## 6. Next action (not performed)
The minimum concrete action for the Broad/DepMap route is a **browser reading of the DepMap portal Terms and Conditions, the CCLE terms and the release README.txt** (which my tools could not read) with a recorded conclusion on whether they add restrictions to CC BY 4.0. For the Sanger route the minimum action is the legal/compliance determination in condition 3; if Sanger-derived results are ever to be incorporated into a product, Sanger's prior consent is the required step. No provider is contacted by this document.
