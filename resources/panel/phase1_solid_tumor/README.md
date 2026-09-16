# GenRichi Phase-1 Solid Tumor SNV/Indel Target Design — v1

**Status:** Draft authoritative Phase-1 resource. NOT yet wired into the pipeline or config. Does NOT replace `resources/panel/comprehensive_genes.bed` (left untouched).

## Provenance
- **Gene list:** the verified 55-gene GenRichi Solid Tumor Comprehensive Panel (see project design-audit conversation history for the removal/addition reconciliation from the original 62-gene placeholder list).
- **Transcript source:** NCBI/Ensembl MANE Select v1.4 GTF (`resources/reference/mane/MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz`, see its README for download provenance). Exactly one MANE_Select-tagged transcript per gene was used.
- **Assembly:** GRCh38 primary assembly, `chr`-prefixed contigs only (`chr1-22, chrX, chrY, chrM`).
- **Region definition:** all CDS exons of the selected MANE Select transcript, each independently padded ±20 bp, then merged (overlapping/adjacent) **within each gene only** — intervals are never merged across different genes, so every row stays attributable to exactly one gene/transcript.
- **Generation method:** one-off script (`build_phase1_bed.py`, run outside the repo's pipeline code, not committed) parsing the MANE GTF in two passes: (1) resolve the single MANE_Select transcript per gene and verify all 55 map cleanly, (2) extract, pad, merge, sort, and dedupe CDS intervals.
- **Generated:** 2026-09-14.

## Files
- `solid_tumor_phase1_v1.bed` — BED4 (`chrom`, `start`, `end`, `gene_symbol`), 0-based half-open, karyotype-sorted (chr1–22, X, Y, M; then position).
- `phase1_audit.json` — full per-gene audit (transcript ID, exon count, raw CDS bp, final padded/merged bp, final interval count) plus totals and the "unusual transcript situation" flags.

## Known limitations of this v1 design (carried forward from the design audit)
- MANE_Select only — 7 genes (BRAF, CDKN2A, FGFR2, HRAS, KRAS, MUTYH, SMARCA4) have an additional MANE_Plus_Clinical transcript not covered by this design; deliberate per this step's instructions, revisit if a clinically relevant alternate-transcript region falls outside MANE_Select CDS.
- This is a **SNV/Indel-only** design. It does not include: TERT promoter, MET intron13/14 splice-flank beyond standard ±20bp, EPCAM deletion-detection region, or any CNV/fusion capability. Those are separate later phases per the GenRichi phased architecture and must not be assumed present.
- Not yet validated against any truth set (HCC1395, GIAB NA12878, reference standards) — validation is a Phase-1 task still pending, not yet performed.
- Sort order used here is karyotype numeric order (chr1, chr2, …, chr22, chrX, chrY, chrM). This must be checked against the actual `hg38.fa.fai`/sequence-dictionary contig order before being used as a GATK `-L` interval file — not yet verified in this step (pipeline not run per instruction).
