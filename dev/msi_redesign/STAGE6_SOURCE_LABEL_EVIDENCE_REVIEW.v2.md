# Stage 6 — Source / Label Evidence Review v2 (primary-source verification; development / RUO)

**Decision Status: NO SOURCE SELECTED.**
Nothing in this document selects, recommends or excludes a cohort. No thresholds, minimum sample counts, partitions or Stage 6B
Analysis Plan are defined. No data (patient, MSI-H/MSS, controlled-access) was downloaded, requested, opened or processed. Only
publications, official documentation pages and vendor documents were read. Code, tests, pipelines, registry schemas, manifests
and datasets were not touched. Supersedes the "UNVERIFIED" items of `STAGE6_SOURCE_LABEL_EVIDENCE_REVIEW.md` (v1, kept unchanged).
Date: 2026-09-25.

## 0. Method and legend
Open-access papers were read as full text (Europe PMC XML / PMC author-manuscript pages) and searched directly for the quoted
statements. Paywalled or bot-blocked pages (Cell, ScienceDirect, PubMed, Nature: HTTP 403 or challenge) were not read; where only a
search-result snippet or a fetch-tool summary was available, the item is at most PARTIALLY VERIFIED. Fetch-tool summaries were
found to be unreliable at least once (item S-01), so any claim resting only on a summary is labelled as such.

* **VERIFIED** — stated in a primary/official text that I read directly.
* **PARTIALLY VERIFIED** — supported by a secondary text, a snippet, or only part of the claim is confirmed.
* **UNVERIFIED** — not established from anything I could read. Not assumed.

Evidence ID = `<source>-<nn>`. Label-tier vocabulary (T1 clinical reference method … T4 GenRichi output, forbidden) is from
`STAGE6_DESIGN_SPECIFICATION.v2.md`, Section C.

## 1. Corrections to v1 (errors found during this verification)
| # | v1 statement | Correction | Evidence |
|---|---|---|---|
| K1 | Cortes-Ciriano: "281 samples had reference labels" (fetch summary) | 281 is the Table-1 total of cases **predicted** MSI-H by their model (confidence 0.75) across all 23 types, not a count of labelled samples | S-01 |
| K2 | MSK-IMPACT validation was "138 colorectal and endometrial cases" | 138 colorectal (24 MSI-H, 114 MSS) **plus** 40 uterine endometrioid (15 MSI-H, 25 MSS) | M-03 |
| K3 | 12,288-case MSK study cited as JCO PO 17.00189 | The 12,288-patient study is Middha et al., JCO Precision Oncology (DOI 10.1200/PO.17.00084, PMC6130812). DOI 10.1200/PO.17.00189 is a different article (Haraldsdottir; not open access) | Europe PMC DOI lookups |
| K4 | CPTAC endometrial MSI method "UNVERIFIED" | Now VERIFIED: computational consensus (C-01) | C-01 |

## 2. Cortes-Ciriano discrepancy (item 5) — RESOLVED
* **S-01 (VERIFIED, primary text: Cortes-Ciriano et al. 2017, Nat Commun, [PMC5467167](https://pmc.ncbi.nlm.nih.gov/articles/PMC5467167/)).**
  Table 1 ("Tumour samples utilized to profile MSI") lists 7,919 tumour/normal pairs over 23 types; its last column total, **281**,
  is the number of cases "predicted as MSI-H at a confidence level of 0.75" (asterisked values). Figure 4 legend: "190 MSI-H,
  118 MSI-L and 522 MSS exomes" and "25 MSI-H, 19 MSI-L and 105 MSS whole genomes". The 190 MSI-H equals the sum of the
  non-asterisked TCGA-labelled MSI-H values in Table 1 for the five MSI-prone types (75 UCEC + 64 STAD + 45 COAD + 3 READ + 3 ESCA).
  So 281 and 190/118/522 are different quantities and do not conflict.
* **S-02 (UNVERIFIED).** The five MSI-prone types have 1,087 samples in Table 1 (265+292+271+76+183) versus 830 exomes with a
  TCGA label (190+118+522). Why 257 samples lack a label is not stated in the text I read; it is not assumed.

## 3. TCGA (colorectal COAD/READ, endometrial UCEC, gastric STAD)
| ID | Question | Finding | Status | Source |
|---|---|---|---|---|
| T-01 | Original MSI assay | "MSI status was evaluated by the TCGA consortium for COAD, READ, ESCA, STAD and UCEC tumours using a panel of four mononucleotide repeats (BAT25, BAT26, BAT40 and TGFBRII) and three dinucleotide repeats (D2S123, D5S346 and D17S250), except for a subset of COAD/READ genomes evaluated by five mononucleotide markers (BAT25, BAT26, NR21, NR24 and MONO27)"; described as a "capillary sequencing-based fragment length assay" | **PARTIALLY VERIFIED** (stated by a secondary analysis paper; the TCGA supplements that hold the assay descriptions were not read) | Cortes-Ciriano [PMC5467167](https://pmc.ncbi.nlm.nih.gov/articles/PMC5467167/), Methods |
| T-02 | Primary-paper statements | COAD 2012: main text gives no MSI assay; reports 23 of 30 hypermutated tumours MSI-H. UCEC 2013: "MSI testing performed on all samples using seven repeat loci" (loci not named), MSI in 40% of endometrioid and 2% of serous tumours. STAD 2014: "Microsatellite instability (MSI) testing was performed on all tumour DNA"; MSI-high = 22% of 295 tumours; no assay detail | **VERIFIED** for these statements; assay detail in main text **UNVERIFIED** | [COAD PMC3401966](https://pmc.ncbi.nlm.nih.gov/articles/PMC3401966/), [UCEC PMC3704730](https://pmc.ncbi.nlm.nih.gov/articles/PMC3704730/), [STAD PMC4170219](https://pmc.ncbi.nlm.nih.gov/articles/PMC4170219/) |
| T-03 | MSI-H definition | "MSI-H (≥40% of markers altered), MSI-L (<40% of markers altered) and MSS (no marker altered)" (Cortes-Ciriano). Bonneville: MSI-PCR classes MSS 0/5, MSI-L 1/5, MSI-H ≥2/5 loci (general Bethesda description, not a TCGA-specific statement). Per-cohort definitions in the TCGA papers: **UNVERIFIED** | **PARTIALLY VERIFIED** | Cortes-Ciriano Methods; Bonneville [PMC5352334](https://pmc.ncbi.nlm.nih.gov/articles/PMC5352334/), Introduction |
| T-04 | Are labels the original clinical assay or a later computational call? | The TCGA labels are assay-derived (above). Cortes-Ciriano obtained them by download: "The MSI status … were downloaded from the GDAC" (Broad Genome Data Analysis Center) — a repackaged copy, one step removed from the assay. Bonneville uses "MSI-PCR status as a gold standard" for TCGA pairs. Computational tools (MANTIS, MSIsensor, the Cortes-Ciriano model) are **evaluated against** these labels, not the source of them | **VERIFIED** that they are assay-derived; **PARTIALLY VERIFIED** for traceability to the original assay file per sample | Cortes-Ciriano Methods; Bonneville Methods, "Tool performance evaluation" |
| T-05 | Where sample-level assay data are documented | GDC publication page for STAD lists "Public Microsatellite Instability Data [tar]" and "Level 1 Microsatellite Instability Data [tar]" (contents not opened). The COAD and UCEC pages list clinical data and biospecimen XML but no MSI file | STAD: **PARTIALLY VERIFIED** (files exist; contents UNVERIFIED). COAD/UCEC: **UNVERIFIED** | GDC [STAD](https://gdc.cancer.gov/about-data/publications/stad_2014), [COAD/READ](https://gdc.cancer.gov/about-data/publications/coadread_2012), [UCEC](https://gdc.cancer.gov/about-data/publications/ucec_2013) |
| T-06 | Orthogonal pathology / MMR-IHC | No MMR immunohistochemistry is described in the main text of the COAD, UCEC or STAD papers. Whether TCGA clinical files carry IHC or per-marker PCR fields: not found | **UNVERIFIED** | as above; my searches did not locate field documentation |
| T-07 | Identifiers and patient-level leakage checks | Barcode: "A parent barcode prefixes any of its descendent barcodes"; the participant segment links all samples/analytes of one patient. Cortes-Ciriano states the full sample list is in its supplementary tables; Bonneville names cases (e.g. TCGA-V5-A7RE). Bonneville also documents that 4 COAD/READ MSI-H pairs were sequenced at both BI and BCM and 20 UCEC MSI-H pairs at both BI and WUGSC — i.e. **one patient can appear as several exome pairs** | **VERIFIED** that patient-level linkage is possible from the barcode and that duplicate sequencing of one patient exists; the supplementary sample tables were **not opened** | GDC [TCGA barcode](https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/); Bonneville Methods |
| T-08 | Sample counts (as published; not natural prevalence) | COAD 2012: 224 tumour-normal exome pairs; UCEC 2013: 373 carcinomas; STAD 2014: 295 patients, fresh-frozen. Bonneville TCGA subsets (curated): COAD/READ 76 pairs (38 MSI-H/38 MSS), UCEC 99 (49/50), STAD 100 (50/50), ESCA 71 (2/69), UCS 53 (2/51), PRAD 59 (1/58); COAD/READ sequenced at Baylor; PRAD from dbGaP, others from CGHub | **VERIFIED** | papers above; Bonneville "Sample data" |
| T-09 | Controlled-access requirements | GDC: "Data in the GDC is considered either open or controlled access"; controlled = raw sequencing (BAM, FASTQ), VCF, protected MAF; two steps: dbGaP project application through a Data Access Committee, and an eRA Commons ID; governed by the NIH Genomic Data Sharing Policy. dbGaP TCGA study phs000178 page: controlled tier requires "user certification through dbGaP Authorized Access" and a Data Use Certification Agreement; controlled types listed include "Primary sequence data (.bam files)"; "somatic mutations or clinical data are open access" (these page statements come from a fetch-tool summary of the dbGaP page, not from my own reading of it). Consent-group and data-use-limitation details, and user obligations (security, retention, sharing, publication): not stated on the pages I could read (the NIH dbGaP request page returned 403) | **VERIFIED** for tiers and process; obligations **UNVERIFIED** | GDC [Controlled Access](https://docs.gdc.cancer.gov/Encyclopedia/pages/Controlled_Access), [Data Access Policy](https://docs.gdc.cancer.gov/Encyclopedia/pages/Data_Access_Policy); dbGaP [phs000178.v10.p8](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000178.v10.p8) |
| T-10 | Specimen type and tumour purity | STAD: "primary tumour tissue (fresh frozen)". TCGA quality threshold "at least 80% tumour nuclei … later reduced to 60%" (visual). Purity estimates (ABSOLUTE, ESTIMATE, LUMP, IHC image analysis) exist for 9,364 TCGA tumours in Aran et al. Supplementary Data 1 (ABSOLUTE for only 11 cancer types); agreement between genomic methods and IHC is low but positive. Specimen preservation for COAD and UCEC is not stated in the main texts I read; FFPE not mentioned in any | STAD frozen and the threshold history: **VERIFIED**. Per-sample purity: **PARTIALLY VERIFIED** (published, table not opened). COAD/UCEC preservation: **UNVERIFIED** | STAD PMC4170219; Aran 2015 [PMC4671203](https://pmc.ncbi.nlm.nih.gov/articles/PMC4671203/), Results "Purity of TCGA tumour samples" |
| T-11 | Tool-leakage risk | The same labelled TCGA samples were used to develop or evaluate MANTIS (Bonneville), the Cortes-Ciriano model (10-fold CV on TCGA labels: sensitivity 92%, specificity 99%), and others; Hause 2016 (Nat Med) is paywalled and was not read | Bonneville, Cortes-Ciriano: **VERIFIED**; Hause: **UNVERIFIED** | PMC5352334; PMC5467167 |

## 4. CPTAC (colon, endometrial)
| ID | Question | Finding | Status | Source |
|---|---|---|---|---|
| C-01 | UCEC MSI method (Dou 2020, Cell) | **Computational consensus**, not a clinical assay. STAR Methods: five criteria — mutation load, MMR-gene mutation, MSIsensor v0.2 score, MSMuTect v1.0 score, MLH1 methylation; two-cluster K-means on four of them; "A sample was officially called MSI-H if it was predicted to be MSI-H by no fewer than 3 of 5 methods (Table S3)". This is a **T3** label | **VERIFIED** | Dou 2020 author manuscript [PMC7233456](https://pmc.ncbi.nlm.nih.gov/articles/PMC7233456/), STAR Methods "Microsatellite Instability Prediction" |
| C-02 | UCEC counts | 95 prospectively collected tumours; "Our cohort included 7 POLE, 25 MSI, 43 CNV-low, and 20 CNV-high tumors". A 2026 secondary paper (EMMA-STRAT, [PMC13474893](https://pmc.ncbi.nlm.nih.gov/articles/PMC13474893/)) lists, from LinkedOmics metadata, 21 MSI-H / 73 MSS / 1 missing for the 95-tumour cohort and a second CPTAC endometrial set of 108 (36 MSI-H / 72 MSS). The 25 vs 21 difference is **unreconciled**; the 108-tumour set was not verified from a primary paper | 95 and 25: **VERIFIED**. 21/73: **PARTIALLY VERIFIED (secondary)**, discrepancy **UNVERIFIED**. 108-set: **UNVERIFIED** | Dou 2020 Results; EMMA-STRAT Table 1 |
| C-03 | Colon MSI method (Vasaikar 2019, Cell) | 110 patients enrolled; whole-exome sequencing of 106 tumours. "The number of MS INDELs showed a clear bimodal distribution, which allowed us to separate the samples into a MSI-H group (n=24) and a MSS group (n=82)" — **computational (T3)** from WXS. "For the 85 samples with PCR-based MSI testing results, WXS-based assignment agreed completely with PCR assignment (Table S2)" | **VERIFIED** (colon text). PCR assay details (loci, who ran it) and Table S2 contents **UNVERIFIED** | Vasaikar author manuscript [PMC6768830](https://pmc.ncbi.nlm.nih.gov/articles/PMC6768830/), Results "Somatic mutations and their proteomic consequences" |
| C-04 | Orthogonal IHC / PCR | Colon: PCR for 85 of 106 tumours (C-03); no MMR-IHC (the only IHC in the text is Human Protein Atlas protein staining). Endometrial: no PCR/IHC statement in the text | Colon PCR: **PARTIALLY VERIFIED**. Colon IHC and endometrial PCR/IHC: **UNVERIFIED** (not found; Table S3/S2 not read) | PMC6768830; PMC7233456 |
| C-05 | Experimental vs computational labels | UCEC: computational (C-01). Colon: computational label with PCR concordance reported for a subset (C-03). Neither cohort's assigned label is a clinical assay result | **VERIFIED** | as above |
| C-06 | Specimens / purity | Colon: OCT-embedded, ischemic time <30 min, ">300 mg … at least 60% tumor cell nuclei and less than 20% necrosis"; tumour purity by ESTIMATE (expression-based, Supplementary Table 1); WXS with the Nextera Rapid Capture Exome kit. UCEC: acceptable segments ">80% viable tumor nuclei, total cellularity >50%, necrosis <20%", cryopulverized from liquid-nitrogen storage; purity by methylation deconvolution and ABSOLUTE; three tumours with purity <10% were excluded. FFPE not used in either | **VERIFIED** (primary) | PMC6768830; PMC7233456 |
| C-07 | Access and data use | UCEC: raw genomic data "from the Genomic Data Commons or upon request from dbGaP (phs001287)"; colon: raw genomics at SRA BioProject PRJNA514017 (access tier not stated), proteomics at the CPTAC Data Portal. GDC CPTAC page: open vs controlled tiers, controlled requires token authentication. dbGaP phs001287: "Authorized Access" portal and a Data Use Certification Agreement; consent groups and detailed limitations "NOT STATED" on the page. TCIA CPTAC-UCEC images (not genomic data) are CC BY 3.0 | Tiers: **VERIFIED**. SRA tier for PRJNA514017 and CPTAC-specific obligations: **UNVERIFIED** | PMC7233456 Data availability; PMC6768830 Data availability; GDC [CPTAC-3_2019_2](https://gdc.cancer.gov/about-data/publications/CPTAC-3_2019_2); dbGaP [phs001287](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001287); [TCIA](https://wiki.cancerimagingarchive.net/pages/viewpage.action?pageId=33948263) |

## 5. Cross-cohort overlap (item 3)
| ID | Finding | Status |
|---|---|---|
| X-01 | **TCGA–CPTAC (retrospective study): documented overlap.** Zhang et al. 2014 (Nature 513:382) profiled proteomes of TCGA colorectal tumours: 95 TCGA tumour samples from 90 patients (64 COAD, 31 READ) — per search-result text of the paper and its data record; the Nature full text was not read | **PARTIALLY VERIFIED** ([Nature](https://www.nature.com/articles/nature13438)) |
| X-02 | **TCGA–CPTAC prospective cohorts (Vasaikar 2019, Dou 2020): no explicit patient-overlap statement found.** The colon paper treats the prospective cohort (n=106) and TCGA colorectal as separate comparators (Supplementary Fig. S1); the endometrial paper has no such statement. A 2026 secondary paper says it cross-checked sample identifiers "to ensure that there was no overlap" among TCGA and two CPTAC sets, without publishing the check | Independence **UNVERIFIED** (secondary statement only: PARTIALLY) |
| X-03 | **MSK-IMPACT vs TCGA/CPTAC and other cohorts:** no published statement found | **UNVERIFIED** |
| X-04 | **Same patient, several files inside TCGA:** documented (T-07) | **VERIFIED** |
| X-05 | **Same tumours used to build tools** (T-11) | **VERIFIED** |

Nothing was inferred; without accessing sample tables, patient-level overlap between cohorts remains open.

## 6. MSK-IMPACT (item 6)
| ID | Finding | Status | Source |
|---|---|---|---|
| M-01 | Published cohort: "data from 12,288 patients who underwent molecular testing with MSK-IMPACT between January 1, 2014, and December 31, 2016"; 13,091 cancer samples, 66 principal cancer types; 11,553 samples from 10,900 patients had "sufficient coverage of 200× and TP of ≥ 25%" for the MSI analysis ("TP" is not expanded in the passage I read; not assumed to mean purity); MSK-IMPACT uses tumour and matched normal DNA | **VERIFIED** | Middha 2017 [PMC6130812](https://pmc.ncbi.nlm.nih.gov/articles/PMC6130812/), Methods/Results |
| M-02 | MSI-H definition: **MSIsensor score ≥ 10** (percentage of unstable microsatellites). The cutoff was established by cross-validation against MSI PCR and MMR IHC, and sensitivity/specificity were estimated on the same validation set (so the concordance figure is not an out-of-sample estimate) | **VERIFIED** | PMC6130812 Methods, Abstract |
| M-03 | Orthogonal comparison: MSI PCR (five loci; "Instability at two or more of the five microsatellite loci defines MSI-H") and/or MMR IHC in 138 CRC (24 MSI-H, 114 MSS) and 40 UEC (15 MSI-H, 25 MSS): concordance 99.4%. Non-CRC/UEC: 68 MSIsensor-high; 49 had material; 46 concordant (39 PCR, 7 IHC), 3 MSI-L by PCR but MMR-deficient by IHC. Two cancers were PCR/IHC MSI-H but <10 by MSIsensor | **VERIFIED** | PMC6130812 |
| M-04 | Raw-data accessibility: the Middha text contains no data-availability statement that I could find (IRB approval only). The related 10,000-patient MSK-IMPACT paper states that "All genomic results and associated clinical data for all patients in this study are publically available in the cBioPortal" — processed results and clinical data, not stated to include BAM files or the MSIsensor scores; >85% of patients consented under an IRB protocol (NCT01775072). Whether the Middha 12,288 cohort is a superset of the 10,000-patient release, and any access route for sequence files, is not established | Middha access: **UNVERIFIED**. cBioPortal release of processed data: **VERIFIED** (Zehir [PMC5461196](https://pmc.ncbi.nlm.nih.gov/articles/PMC5461196/)). Sequence-level access, cohort mapping: **UNVERIFIED** | PMC6130812; PMC5461196 |
| M-05 | Status for this project | Published clinical evidence for a **targeted-panel** MSIsensor method; **not** an available validation cohort. Independence from other cohorts is not established (X-03) and raw-data access is not established (M-04). Not treated as usable | as above |

## 7. Tumour purity and specimen information (item 4) — what the public metadata establishes
* TCGA: fresh-frozen primary tumour tissue is stated for STAD (VERIFIED); nuclei threshold 80% then 60% (VERIFIED, visual); per-sample
  purity estimates from four methods exist in a published supplement, ABSOLUTE for only 11 cancer types (PARTIALLY VERIFIED, table not opened);
  purity/preservation for every individual sample: **UNVERIFIED**.
* CPTAC: frozen/OCT specimens with pathology thresholds stated per cohort (VERIFIED, C-06); sample-level purity estimates are published
  as supplementary tables (ESTIMATE colon; ABSOLUTE and methylation-based UCEC) — **not opened**, so availability per sample is PARTIALLY VERIFIED.
* No source's documentation, as read, provides FFPE status, cold-ischaemia time or fixation for every sample; nothing is assumed.

## 8. Items carried from v1 (not re-verified in this pass)
Seraseq MSI reference materials: VERIFIED from the vendor product sheet (research use only; MSI-High call by TSO500 with score >20%; ddPCR/qPCR-CE for AF products;
marker coordinates in hg19; contrived cell-line/plasmid material) — [product sheet](https://www.seracare.com/globalassets/seracare-resources/ps-mkt-00534_0710-1670.0710-2236.0710-1675.0710-1676_seraseq_msi-reference_materials.pdf).
Cell-line resources: PARTIALLY VERIFIED (Sanger MSIsensor-pro score ≥7 = MSI; CCLE classification from deletion counts — both **T3**; data-use terms **UNVERIFIED**).

## 9. What remains UNVERIFIED (complete list)
1. Original TCGA assay descriptions and MSI-H definitions in the TCGA supplements (COAD, UCEC, STAD); the seven UCEC loci names (T-01/T-02/T-03).
2. Per-sample traceability from label to original assay file for COAD and UCEC; contents of STAD "Public/Level 1 Microsatellite Instability Data" (T-04/T-05).
3. Any orthogonal pathology/MMR-IHC evidence for TCGA (T-06); and for CPTAC (C-04: colon IHC; endometrial PCR/IHC; Table S2/S3 contents).
4. Contents of the TCGA/Cortes-Ciriano supplementary sample tables, hence actual patient-level leakage checks (T-07); why 257 of 1,087 samples in the five MSI-prone types lack a label (S-02).
5. Hause 2016 label details (paywalled) (T-11).
6. TCGA COAD and UCEC specimen preservation; per-sample purity/cellularity availability for TCGA and CPTAC (T-10, C-06 tables not opened).
7. CPTAC UCEC count discrepancy 25 (primary) vs 21 (secondary) and the existence and label method of the 108-tumour endometrial set (C-02).
8. PCR assay details (loci, laboratory) behind the 85 colon comparisons (C-03).
9. Access tier of SRA BioProject PRJNA514017; CPTAC-specific data-use limits; dbGaP consent groups, data-use limitations and user obligations for TCGA and CPTAC (T-09, C-07); the NIH dbGaP request page and Genomic Data Sharing Policy were not readable (HTTP 403).
10. Patient-level overlap between TCGA and the CPTAC prospective cohorts, and between MSK-IMPACT and any other cohort (X-02, X-03); the Zhang 2014 full text (X-01) was not read.
11. MSK-IMPACT: any data-availability route for sequence files; mapping between the 12,288-patient and the 10,000-patient releases; MSIsensor scores in the public release (M-04).
12. Vendor/cell-line items in Section 8 that were not re-examined (which cell line underlies Seraseq material; Sanger/CCLE data-use terms).

## 9a. Decision Status
**NO SOURCE SELECTED.** No source is recommended, excluded, or ranked. The Stage 6B Analysis Plan, minimum sample counts, partitions and
thresholds are not defined. Stage 6 remains `DESIGN FROZEN FOR REVIEW`; nothing is validated.

## 10. Reference list (documents actually read)
Cortes-Ciriano 2017 PMC5467167 · Bonneville 2017 PMC5352334 · TCGA COAD 2012 PMC3401966 · TCGA UCEC 2013 PMC3704730 · TCGA STAD 2014 PMC4170219 ·
Aran 2015 PMC4671203 · Vasaikar 2019 PMC6768830 · Dou 2020 PMC7233456 · Middha 2017 PMC6130812 · Zehir 2017 PMC5461196 · EMMA-STRAT PMC13474893 (secondary) ·
GDC Data Access Policy / Controlled Access / TCGA Barcode pages · GDC publication pages (COAD/READ 2012, UCEC 2013, STAD 2014, CPTAC-3 2019) ·
dbGaP study pages phs000178.v10.p8 and phs001287 · TCIA CPTAC-UCEC wiki · SeraCare product sheet MKT-00534-03 ·
search-result text only (not read in full): Zhang 2014 Nature 513:382; Hause 2016 Nat Med; Sanger Cell Model Passports / DepMap documentation.
