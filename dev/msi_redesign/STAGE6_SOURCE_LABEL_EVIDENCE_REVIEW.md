# Stage 6 — Source / Label Evidence Review (development / RUO)

Status: REVIEW ONLY. Literature and vendor/documentation review. No data downloaded, requested or opened; no patient,
MSI-H or MSS data accessed; no thresholds, partitions, minimums or Stage 6B Analysis Plan. Nothing here selects a source.
Date: 2026-09-25. Every claim carries a citation; anything not established from a source I could actually read is
marked **UNVERIFIED**. Several publisher pages (PubMed, Nature, PMC HTML) were blocked by bot checks; text was read via Europe PMC
XML where possible. Summaries produced by the fetch tool are noted as such and should be re-read against the primary text before
any decision.

## 0. Label-tier reminder (from v2, Section C)
T1 clinical reference method (PCR/CE panel, MMR IHC) · T2 published expert annotation · T3 computational call by another tool ·
T4 GenRichi output (forbidden). Presence of a label does not make a sample a clinical reference; the label's method is part of
the evidence chain.

## 1. TCGA (colorectal COAD/READ, endometrial UCEC, gastric STAD, others)

| Item | Finding | Source / status |
|---|---|---|
| Sample type | Primary tumor tissue with matched normal, whole-exome sequencing. COAD 2012: 224 tumor-normal pairs in exome analysis; UCEC 2013: 373 carcinomas. Specimen preservation (fresh-frozen vs FFPE) | COAD [PMC3401966](https://pmc.ncbi.nlm.nih.gov/articles/PMC3401966/), UCEC [PMC3704730](https://pmc.ncbi.nlm.nih.gov/articles/PMC3704730/). Preservation: **UNVERIFIED** |
| MSI label and how produced | Labels for five tumor types (COAD, READ, ESCA, STAD, UCEC) are described as clinically determined by capillary fragment-length assay: four mononucleotide (BAT25, BAT26, BAT40, TGFBRII) and three dinucleotide repeats; MSI-H if ≥40% of markers altered, MSI-L below that, MSS if none. UCEC paper: "MSI testing performed on all samples using seven repeat loci" (loci not named in the text I read). Bonneville et al. use "MSI-PCR status as a gold standard" for TCGA pairs | Cortes-Ciriano 2017 [PMC5467167](https://pmc.ncbi.nlm.nih.gov/articles/PMC5467167/) (via Europe PMC, tool summary); UCEC PMC3704730; Bonneville [PMC5352334](https://pmc.ncbi.nlm.nih.gov/articles/PMC5352334/). Exact marker panel per cancer type and testing centre: **UNVERIFIED** (COAD 2012 text: the MSI assay itself was not identified in what I read; it reports 23 of 30 hypermutated tumors MSI-H) |
| Clinical / pathology evidence | MMR IHC availability per case, pathology review of the tested sample, Lynch/germline status: **UNVERIFIED**. UCEC paper links MSI to MLH1 promoter methylation | PMC3704730 for methylation only |
| Sample-level availability | Bonneville's TCGA subsets: COAD/READ 76 pairs (38 MSI-H, 38 MSS), UCEC 99 (49/50), STAD 100 (50/50), ESCA 71 (2/69), UCS 53 (2/51), PRAD 59 (1 MSI-positive/58) — a curated, partly balanced selection, not natural prevalence. Another paper reports 78/85/156 MSI-H and 245/265/274 MSS for COAD/STAD/UCEC with hg38 BAMs from the GDC (secondary source, not checked) | PMC5352334; Briefings in Bioinformatics 2024 [bbae390](https://academic.oup.com/bib/article/25/5/bbae390/7731494) (search snippet only). Full per-case label table location: **UNVERIFIED** |
| Cortes-Ciriano counts | Fetch summary gave 281 labelled samples with 190/118/522 by class, which do not add up | **UNVERIFIED** — re-read the primary text |
| Terms of use | BAM/FASTQ/VCF/protected MAF are controlled access: dbGaP authorisation via a Data Access Committee under the NIH Genomic Data Sharing Policy; open-access users must not attempt to identify participants. Detailed controlled-data obligations (storage, retention, sharing, publication): **UNVERIFIED** | GDC [Data Access Policy](https://docs.gdc.cancer.gov/Encyclopedia/pages/Data_Access_Policy), [Controlled Access](https://docs.gdc.cancer.gov/Encyclopedia/pages/Controlled_Access) |
| Development vs independent validation | Development: plausible only after DAC approval. Independent validation: **not by default** — the same TCGA labels and samples have been used to build/evaluate MANTIS, MSIsensor, the Hause 2016 classifier and the Cortes-Ciriano method, so any locus or threshold chosen with TCGA labels cannot be validated on TCGA. Regulatory acceptability of research-lab PCR labels as a clinical reference: **UNVERIFIED** | PMC5352334; PMC5467167; Hause 2016 [Nature Medicine](https://www.nature.com/articles/nm.4191) (paywalled; abstract only) |
| Leakage / splitting | Split by patient. Cohort composition differs by cancer type (few MSI-H in ESCA/UCS/PRAD). MSI-L is a third class that must be handled explicitly. Hypermutated MSS (POLE) tumors are expected in endometrial data. Tumor purity per sample: **UNVERIFIED**. Overlap of patients with CPTAC: **UNVERIFIED** (CPTAC is described as a separate prospective cohort) | as above |

## 2. CPTAC (colon, endometrial and others; via GDC / Proteomic Data Commons)

| Item | Finding | Source / status |
|---|---|---|
| Sample type | Prospectively collected tumors with adjacent tissue and blood normal; DNA from whole-genome and whole-exome sequencing, harmonised to GRCh38. Endometrial study: 95 tumors (83 endometrioid, 12 serous). Colon study: cohort size **UNVERIFIED** (110 was my assumption, not confirmed) | GDC/CPTAC [release note](https://proteomics.cancer.gov/news_and_announcements/gdc-releases-harmonized-genomic-data-cptac); Dou 2020 [Cell](https://www.cell.com/cell/fulltext/S0092-8674(20)30107-0) (search summary); Vasaikar 2019 [GDC page](https://gdc.cancer.gov/about-data/publications/CPTAC-3_2019_2) |
| MSI label and how produced | Dou 2020 has a supplementary table titled "MSI Status Determination" (Table S3), so the method is documented there. The method itself (PCR, IHC, or computational) is **UNVERIFIED**. Colon: **UNVERIFIED** | Dou 2020 supplementary listing (search result); full text non-open-access, fetch failed |
| Clinical / pathology evidence | **UNVERIFIED** | — |
| Availability | Controlled VCF/TSV/BAM at the GDC; number of MSI-H vs MSS per cohort **UNVERIFIED** | GDC release note |
| Terms | Controlled data: dbGaP/DAC. CPTAC-specific terms **UNVERIFIED** | GDC policy pages |
| Suitability | Potentially a different cohort from TCGA (helpful for independence), but not assessable until the label method is read. **UNVERIFIED** | — |
| Leakage | Patient-level overlap with TCGA **UNVERIFIED** | — |

## 3. Seraseq MSI reference materials (LGC SeraCare)

| Item | Finding | Source |
|---|---|---|
| Sample type | Manufactured, "cell line or plasmid-based" reference material. Products: gDNA MSI-High Mix (tumor-only), FFPE MSI-High RM (one 10 µm curl, >200 ng), and MSI Reference Panel Mix AF5%/AF20% (tumor + normal gDNA; normal background GM24385, "known to be microsatellite stable") | [Product sheet MKT-00534-03](https://www.seracare.com/globalassets/seracare-resources/ps-mkt-00534_0710-1670.0710-2236.0710-1675.0710-1676_seraseq_msi-reference_materials.pdf), [product page](https://www.seracare.com/Seraseq-gDNA-MSIHigh-Mix-0710-1670/) |
| MSI label and method | MSI-High call comes from targeted NGS (TSO500): average MSI score 77.1 (gDNA) and 75.6/71.4 (FFPE); "the value must be >20%" for High. AF5%/AF20% products: variants quantitated by ddPCR and qPCR/CE fragment-length analysis. Markers BAT-25, BAT-26, NR-21, NR-24, MONO-27, coordinates given in **hg19** (build conversion needed) | product sheet |
| Clinical / pathology evidence | None: contrived material. Which cell line(s) underlie the "human diseased cell line-based" material: **UNVERIFIED** | product sheet |
| Availability | Commercial, lot-based, material numbers 0710-1670/-2236/-1675/-1676 | product sheet |
| Terms | Product sheet states "For research use only. Not for use in diagnostic procedures"; cGMP / ISO 13485 manufacturing. Data-sharing/publication terms **UNVERIFIED** | product sheet |
| Suitability | Analytical use (positive control, limit of detection, tumor-only behaviour). Not patient-representative, not a cohort, single material. The MSI-High label is the manufacturer's NGS call (T3-like for the mix, ddPCR/CE-quantitated for the AF products). Not an independent-validation set | product sheet |
| Leakage | Lots of one material are replicates of one biological source; the same normal background DNA appears across products | product sheet |

## 4. Cancer cell lines with public sequencing (Cell Model Passports / Sanger; DepMap/CCLE)

| Item | Finding | Source |
|---|---|---|
| Sample type | Immortalised cancer cell lines (e.g. HCC1395 for MSS in this project; MSI-H colorectal/gastric lines). Matched normal generally absent except paired lines such as HCC1395/HCC1395BL; for others **UNVERIFIED** | Cell Model Passports [about](https://cellmodelpassports.sanger.ac.uk/documentation/cellmodelpassports/about) |
| MSI label and method | Sanger DepMap: MSI estimated with MSIsensor-pro, score ≥7 = MSI, otherwise MSS (**T3, computational**). CCLE/DepMap (Ghandi 2019): MSI/MSS/indeterminate classified from the number of short deletions and the fraction of deletions in microsatellite regions across WES/WGS/hybrid-capture data (**T3**) | Sanger [MSI, Ploidy & Mutational Burden](https://depmap.sanger.ac.uk/documentation/cell-models/msi-ploidy-mutational-burden/); [DepMap forum thread](https://forum.depmap.org/t/msi-annotations/644); [CCLE on AWS](https://registry.opendata.aws/depmap-omics-ccle/) |
| Clinical / pathology evidence | Cell-line-level, not patient-level; MMR status can be checked from mismatch-repair gene alterations (not done) | — |
| Availability | Roughly 1,000 cell lines with WGS/WES/RNA-seq (CCLE dataset); per-line raw-data access tier **UNVERIFIED** | AWS registry entry |
| Terms | A Sanger DepMap "Data Usage Policy" exists; its content and CCLE terms: **UNVERIFIED** | Sanger DepMap documentation index |
| Suitability | Development, positive/negative controls, method checks. Not independent validation (T3 labels, in vitro drift, no clinical pathology, tumor-only) | as above |
| Leakage | The same lines appear in several databases and are used in many MSI tool papers; line identity and passage must be tracked as one biological unit | as above |

## 5. Institutional clinical NGS datasets (MSK-IMPACT and similar)

| Item | Finding | Source |
|---|---|---|
| Sample type | Clinical tumor samples on a large targeted panel (468 genes, ~1,000 microsatellites); 12,288 advanced solid cancers scored with MSIsensor | Middha 2017 [JCO PO](https://ascopubs.org/doi/10.1200/PO.17.00189) (search summary) |
| Label | MSIsensor call (T3); validated in 138 colorectal and endometrial cases against MMR IHC and PCR (99.4% concordance, NGS slightly more sensitive) | same |
| Availability / terms | Raw sequence access for external users: **UNVERIFIED** (likely not public) | — |
| Suitability | Not assessable without access; panel data would not match WES/WGS loci anyway. **UNVERIFIED** | — |

## 6. Local data (for completeness)
HCC1395: one MSS unit, external label (T2/T3, database annotation), chr17 subset, four orders = one unit. NA12878 / G8_S13 /
NA12878_somatic: no MSI label (see Stage 6A report). No local MSI-H sample exists.

## 7. What this review does and does not establish
- Established (with citations): how TCGA labels were produced in principle (capillary PCR assays), that BAMs are controlled-access,
  that Seraseq is RUO contrived material with an NGS/ddPCR-derived label in hg19, and that cell-line MSI labels are computational.
- Not established (UNVERIFIED): per-sample TCGA/CPTAC label tables and testing methods per cohort, CPTAC MSI method, pathology/IHC
  availability, patient overlap between cohorts, purity, specimen preservation, and the detailed controlled-data obligations.
  Each needs a primary-source reading (supplementary tables, DAC terms), not a data download.
- No source is recommended or excluded here. Choice, partitions, minimums and the 6B plan remain product-owner decisions after
  the UNVERIFIED items are resolved.
