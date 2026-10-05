# Stage 6 — Source / Label Evidence Review v4: public-documentation verification of v3 gap items (development / RUO)

**Decision Status: NO SOURCE SELECTED.** No source or cohort is recommended, excluded or ranked; no development/validation cohort, threshold,
minimum count, partition or Stage 6B Analysis Plan is defined. v1–v3 are unchanged. The v3 classifications (BLOCKER / REQUIRED BEFORE DATA ACCESS /
NON-BLOCKING / OUT OF SCOPE) are **not changed here**; this document only updates the verification status of each item. Date: 2026-09-25.

## 0. Boundary followed
Authorised and used: publicly accessible supplementary/metadata tables, published labels and assay annotations in them, public identifiers only to
judge provenance/traceability/overlap, and official public documentation. **Not done:** no controlled-access or dbGaP-controlled record was accessed, no
FASTQ/BAM/CRAM/VCF was downloaded or processed, nothing requiring authorisation was opened, no individual sample was selected, no cohort constructed.
**Privacy:** per-sample tables were parsed in memory for column structure, counts and cross-tabulations only; no sample or patient identifier is reproduced here.
The downloaded public tables were kept only in the session scratch area and deleted after this review.

Public tables actually read (aggregate use only): cBioPortal datahub clinical tables for CPTAC colon (`coad_cptac_2019`), CPTAC endometrial
(`ucec_cptac_2020`), TCGA legacy studies (`coadread_tcga`, `ucec_tcga`, `stad_tcga`) and `msk_impact_2017`; open-tier TCGA BCR "biotab" clinical/biospecimen
tables for COAD, READ, UCEC, STAD served by the public GDC API (all reported `access: open`); the GDC-hosted "Public Microsatellite Instability Data" documentation
archive of the TCGA STAD paper (MAGE-TAB IDF/SDRF).

**Blocked or unreadable (recorded, not circumvented):** PMC supplementary-file downloads (Vasaikar Table S2, Dou Table S3, Aran Supplementary Data 1, the TCGA
supplements) are behind an anti-automation proof-of-work challenge; I did not bypass it. NIH pages `sharing.nih.gov` and the dbGaP request page returned HTTP 403; the
MSK data catalogue returned 403; the Proteomic Data Commons guidelines page rendered no text; the two dbGaP DUC PDFs for phs001287/phs000178 were downloaded but their fonts
could not be decoded with the tools available, so **their content was not read**.

## 1. Verification status by v3 item
Status: **VERIFIED** / **PARTIALLY VERIFIED** / **UNVERIFIED**. "v3 class" is shown for reference only.

| Item | Evidence ID | v3 class | Status now | What was established (sources below) |
|---|---|---|---|---|
| UV-01 | T-01, T-03 | REQUIRED BEFORE DATA ACCESS | **PARTIALLY VERIFIED** | STAD: the GDC-hosted TCGA MAGE-TAB documentation names the investigator site (Nationwide Children's Hospital) and the protocol types `fragment_analysis` ("Fragment analysis to produce .fsa files") and `fsa2txt` ("Conversion of .fsa files to txt report") — i.e. a capillary fragment-analysis assay, as stated by Cortes-Ciriano. Marker panel and the MSI-H threshold are **not** in that document. COAD/READ/UCEC: no equivalent primary document found (GDC legacy archive returns 410 Gone; paper supplements gated). The ≥40%-of-markers definition and the two marker panels remain known only from Cortes-Ciriano (secondary). [P1, P2] |
| UV-02 | T-02 | NON-BLOCKING | **UNVERIFIED** | The seven UCEC loci are not named in any document I could read. |
| UV-03 | T-04, T-05 | REQUIRED BEFORE DATA ACCESS | **PARTIALLY VERIFIED** | STAD: the SDRF maps each extract/TCGA barcode to derived "Level 1" data files with an "Include for Analysis" flag (734 rows), so sample-level traceability to the assay files is documented structurally; the Level-1 files themselves were not accessed. COAD/UCEC: no such map found. [P2] |
| UV-04 | T-06 | NON-BLOCKING | **PARTIALLY VERIFIED** | Open-tier TCGA clinical tables carry site-recorded MSI and MMR-IHC fields for COAD and READ only: COAD (459 cases) `microsatellite_instability` YES 11 / NO 81 / unknown 17 / not available 350, with loci-tested/abnormal counts for a subset, and MMR-IHC "tested" YES 56 / NO 287; READ (about 170 cases): MSI field for 26, IHC tested YES 8. The UCEC and STAD patient tables have **no** MSI or IHC fields. These are clinical-record fields (site-reported); their agreement with the TCGA consortium label was **not** assessed (that would need per-sample consortium labels). [P3, P4] |
| UV-05 | C-04 | BLOCKER | **PARTIALLY VERIFIED** | CPTAC endometrial: the public table carries MMR-IHC for MLH1, MSH6, PMS2 and a column labelled "MLH2" (presumably MSH2; not assumed); interpretable for 36 of 95 tumours (57 "cannot be determined", 2 unknown). Aggregate cross-tab against the label: 25 MSI-H → 8 with any loss, 17 not interpretable, 0 with all intact; 70 MSS → 27 all intact, 1 with loss, 42 not interpretable (35 of 36 interpretable tumours agree; descriptive only, not a performance claim). CPTAC colon: no IHC field; PCR for 85 of 106 tumours is stated in the paper, but the per-sample PCR table (Table S2) is gated and was not read. So orthogonal clinical evidence exists for a subset in both cohorts; labels remain computational. [P5, P6, P7] |
| UV-06 | C-03 | REQUIRED BEFORE DATA ACCESS | **UNVERIFIED** | The PCR assay, loci and laboratory behind the 85 colon comparisons are in Table S2, which could not be downloaded. |
| UV-07 | T-07 | REQUIRED BEFORE DATA ACCESS | **PARTIALLY VERIFIED** | Public TCGA tables carry barcodes; the GDC barcode documentation states the participant segment links all samples of a patient. Within-cohort patient linkage is visible: COAD/READ 640 samples for 636 patients, UCEC 549/548, STAD 478/478 (cBioPortal legacy tables). Cross-cohort checks are not possible from these tables (see UV-18). The Cortes-Ciriano supplementary sample lists were not obtained. [P8, P9] |
| UV-08 | S-02 | REQUIRED BEFORE DATA ACCESS | **UNVERIFIED** | Which of the 1,087 samples have a consortium reference label cannot be determined without the consortium label table / supplementary lists, which were not obtained. |
| UV-09 | T-11 | OUT OF SCOPE | **UNVERIFIED** | Not pursued (paywalled; deferred as classified). |
| UV-10 | T-10 | REQUIRED BEFORE DATA ACCESS | **PARTIALLY VERIFIED** | Open TCGA biospecimen sample tables carry an `is_ffpe` flag. Primary tumours flagged FFPE: COAD 13 of 476, UCEC 4 of 552, READ 0 of 172, STAD 0 of 478. `preservation_method` is not populated (994 of 995 COAD rows "not available"); an `oct_embedded` flag is partial. So FFPE material is present in TCGA colorectal and endometrial data, and "not FFPE" is not the same as a documented frozen specimen. [P10] |
| UV-11 | T-10, C-06 | REQUIRED BEFORE DATA ACCESS | **PARTIALLY VERIFIED** | TCGA: pathology-recorded `tumor_nuclei_percent` is public for essentially every tumour (COAD 462, READ 168, UCEC 544, STAD 443–444 rows); observed range COAD 40–100 (median 70), READ 60–100 (70), UCEC 40–100 (80), STAD 20–100 (75) — the minima are below the 60%/80% policy thresholds quoted by Aran et al., so the policy is not a per-sample guarantee. The genomic purity estimates (ABSOLUTE etc.) are in an Aran supplement that was not obtainable. CPTAC endometrial public table: `PURITY_CANCER` for 91 of 95 and an ESTIMATE score for 95 of 95. CPTAC colon public table: no purity column (ESTIMATE purity is in Supplementary Table 1, not obtainable). MSK public release: `TUMOR_PURITY` non-missing for 10,475 of 10,945 samples. [P4, P6, P10, P11, P12] |
| UV-12 | C-02 | REQUIRED BEFORE DATA ACCESS | **PARTIALLY VERIFIED** | Endometrial: the primary paper (25 MSI, of 95) and the public cBioPortal table (25 MSI-H, 70 MSS, 95 with a label) agree; the 2026 secondary paper's 21 MSI-H / 73 MSS / 1 missing disagrees with both, so the authoritative count is 25 — the cause of the secondary paper's figure is **UNVERIFIED** (not needed). Colon: the primary text gives 24 MSI-H / 82 MSS (106 exomes) while the public table gives 24 MSI-H / 81 MSS / 5 not available (110 patients, 106 sequenced) — a one-tumour difference whose cause is unknown. [P5, P6, P7] |
| UV-13 | C-02 | NON-BLOCKING | **UNVERIFIED** | Existence and label method of the second endometrial set (N=108, secondary source) not examined. |
| UV-14 | C-07 | REQUIRED BEFORE DATA ACCESS | **PARTIALLY VERIFIED** | BioProject PRJNA514017: raw sequence reads, 220 BioSample records, "No public data is linked to this project"; no access tier stated. Consistent with restricted access but not proof. [P13] |
| UV-15 | C-07 | BLOCKER | **UNVERIFIED** | CPTAC-specific data-use limitations were not established: the phs001287 DUC PDF is undecodable with my tools and the PDC guidelines page rendered no text. Only the framework is documented: the CPTAC-3 open RNA-seq data on the AWS registry are "governed by the NIH Genomic Data Sharing Policy" with no commercial-use statement; GDC controlled access needs dbGaP authorisation. [P14, P15] |
| UV-16 | T-09 | BLOCKER | **PARTIALLY VERIFIED** | The TCGA Data Use Certification Agreement (version dated 20 Aug 2014, hosted by NCI) was read: access is granted to a named PI and institution; "Research use will occur solely in connection with the approved research project described in the DAR"; new uses need a new DAR; the text contains **no** statement about commercial or for-profit use (word not present). Whether developing a diagnostic product falls within an approved research use is therefore **not answered by the document**; the dbGaP consent groups/data-use limitations for phs000178 were not readable. The current DUC version (dbGaP PDF) could not be decoded. [P16] |
| UV-17 | T-09, C-07 | REQUIRED BEFORE DATA ACCESS | **PARTIALLY VERIFIED** | Operational obligations in the 2014 TCGA DUC: annual progress update or final report, data destruction (confirmed in writing by the Signing Official) on expiry/non-renewal/DAC request/inconsistent use, non-identification, security best practices (approved user information posted publicly on dbGaP), collaborators at other institutions need their own DAR, citations may be posted. GDC open-data users must acknowledge datasets and not attempt identification. CPTAC-specific obligations: not established. [P16, P17] |
| UV-18 | X-02 | BLOCKER | **UNVERIFIED** | The public CPTAC tables use identifier formats different from TCGA barcodes and contain no TCGA identifier; that neither shows nor excludes overlap. No published overlap statement found. |
| UV-19 | X-03 | BLOCKER | **UNVERIFIED** | No published statement of overlap between MSK-IMPACT and TCGA/CPTAC found. |
| UV-20 | X-01 | OUT OF SCOPE | **UNVERIFIED** | Not pursued (deferred as classified). |
| UV-21 | M-04 | BLOCKER | **UNVERIFIED** | No public sequence-level access route for MSK-IMPACT was found; the MSK data catalogue returned 403. What is public is processed data on cBioPortal (below). |
| UV-22 | M-04 | REQUIRED BEFORE DATA ACCESS | **PARTIALLY VERIFIED** | The public `msk_impact_2017` clinical sample table has 10,945 samples from 10,336 patients and 18 columns, **none an MSI score/status**; specimens: FFPE 8,982, DNA 1,941, cell pellet 20, other 1, FNA 1 (about 82% FFPE); sample type primary 6,213 / metastasis 4,732. The Middha cohort (12,288 patients, 2014–2016) and this release (10,336 patients) have different sizes; their relationship is not documented. [P12, P18] |
| UV-23 | CL-01 | NON-BLOCKING | **UNVERIFIED** | Not examined. |
| UV-24 | CL-02 | BLOCKER | **PARTIALLY VERIFIED** | Sanger Cancer Dependency Map Data Usage Policy (official page): a "non-exclusive, non-transferable right to use data files for internal proprietary research and educational purposes, including target, biomarker and drug discovery"; excluded: "resale … or for provision of commercial services"; "Commercial use of the Cancer Dependency Map at Sanger API and/or data is not permitted without prior consent." Which datasets (e.g. Cell Model Passports) it covers is not enumerated. Broad DepMap/CCLE: the terms page rendered no text; a moderator's forum post says most CCLE files in the public DepMap release use CC-BY 4.0 with some unclear, and quotes third-party-rights language; commercial-use terms **UNVERIFIED** (a search snippet claims a non-commercial restriction; not verified from the primary page). [P19, P20] |

## 2. Counts
| Status | Items | Count |
|---|---|---|
| VERIFIED | (none) | **0** |
| PARTIALLY VERIFIED | UV-01, 03, 04, 05, 07, 10, 11, 12, 14, 16, 17, 22, 24 | **13** |
| UNVERIFIED | UV-02, 06, 08, 09, 13, 15, 18, 19, 20, 21, 23 | **11** |
| Total | | **24** |

Among the seven v3 BLOCKERs: UV-05 and UV-16 and UV-24 are PARTIALLY VERIFIED; UV-15, UV-18, UV-19 and UV-21 remain UNVERIFIED. None is fully closed.

## 3. Observations that bear on later gates (facts only; no recommendation)
* The Sanger policy is restrictive for commercial use without consent, and the TCGA DUC ties use to an approved research project without addressing commercial purposes.
* FFPE material exists in TCGA colorectal/endometrial data and dominates the public MSK release, so pre-analytical variables differ across sources.
* Labels: CPTAC endometrial labels are computational (verified in v2) but the public table shows partial MMR-IHC for 36 of 95 tumours; CPTAC colon labels are computational with PCR for 85 of 106 per the paper.
* The TCGA consortium's MSI assay is documented as capillary fragment analysis at Nationwide Children's Hospital for STAD; per-cohort marker panels and thresholds remain secondary-only.
* Site-recorded MSI/IHC fields in TCGA are limited to COAD/READ and cover a minority of cases; their relation to the consortium label is untested.

## 4. Remaining UNVERIFIED (exact)
UV-02 seven UCEC loci · UV-06 colon PCR assay details · UV-08 which of 1,087 samples have labels · UV-09 Hause 2016 · UV-13 the 108-tumour endometrial set · UV-15 CPTAC data-use limitations
(DUC/PDC content unread) · UV-18 TCGA–CPTAC prospective overlap · UV-19 MSK overlap · UV-20 Zhang 2014 full text · UV-21 MSK sequence-level access route · UV-23 Seraseq underlying cell line.
Partially verified items still carry residual gaps named in Section 1 (chiefly: TCGA COAD/READ/UCEC assay descriptions and thresholds in primary supplements, the consortium label tables,
CPTAC supplementary Tables S2/S3, current dbGaP DUC and consent-group text, Broad DepMap/CCLE terms).

## 5. Sources (public documents actually read; P-numbers used above)
P1 GDC publication page, TCGA STAD 2014 (links the public MSI documentation archive) · P2 TCGA STAD MSI MAGE-TAB IDF/SDRF (GDC-hosted public archive, "Public Microsatellite Instability Data") ·
P3 GDC API file listing and open-tier BCR biotab clinical tables, TCGA-COAD/READ/UCEC/STAD · P4 Aran 2015 PMC4671203 (v2) and BCR biospecimen/SSF tables ·
P5 Dou 2020 PMC7233456 (v2) · P6 cBioPortal datahub `ucec_cptac_2020` · P7 Vasaikar 2019 PMC6768830 (v2) and cBioPortal datahub `coad_cptac_2019` ·
P8 GDC TCGA barcode documentation · P9 cBioPortal datahub legacy TCGA studies · P10 open TCGA BCR biospecimen sample/slide tables ·
P11 Aran 2015 (v2) · P12 cBioPortal datahub `msk_impact_2017` · P13 NCBI BioProject PRJNA514017 · P14 AWS Registry of Open Data, CPTAC-3 · P15 GDC data access pages ·
P16 TCGA Data Use Certification Agreement, NCI-hosted PDF dated 20 Aug 2014 · P17 GDC Data Access Policy / GDC Policies · P18 Middha 2017 PMC6130812 and Zehir 2017 PMC5461196 (v2) ·
P19 Sanger Cancer Dependency Map Data Usage Policy · P20 DepMap community forum thread on CCLE terms; AWS registry entry for DepMap/CCLE omics.

Decision Status: NO SOURCE SELECTED
