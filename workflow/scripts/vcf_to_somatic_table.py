"""Parse VEP-annotated somatic VCF (Mutect2) and emit a flat TSV for reporting."""

import csv
import gzip
import re

vcf_path: str = snakemake.input.vcf   # type: ignore[name-defined]
out_path: str = snakemake.output.tsv  # type: ignore[name-defined]

FIELDNAMES = [
    "chrom", "pos", "ref", "alt",
    "gene", "transcript", "consequence", "impact",
    "hgvsc", "hgvsp", "exon", "biotype",
    "gt", "depth", "alt_reads", "vaf",
    "clinvar_sig", "clinvar_disease",
    "cosmic_id", "cosmic_count",
    "gnomad_af",
    "sift", "polyphen",
    "filter",
]


def _parse_csq_header(header_lines: list[str]) -> list[str]:
    for line in header_lines:
        if line.startswith("##INFO=<ID=CSQ"):
            m = re.search(r"Format: ([^\"]+)", line)
            if m:
                return m.group(1).rstrip('"').split("|")
    return []


def _canonical(csq_entries: list[dict]) -> dict:
    canonical = [e for e in csq_entries if e.get("CANONICAL") == "YES"]
    ranked = canonical or csq_entries
    coding = [e for e in ranked if e.get("BIOTYPE") == "protein_coding"]
    return (coding or ranked)[0] if ranked else {}


def _parse_vcf(path: str) -> list[dict]:
    opener = gzip.open if path.endswith(".gz") else open
    header_lines: list[str] = []
    variants: list[dict] = []

    with opener(path, "rt") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("##"):
                header_lines.append(line)
                continue
            if line.startswith("#CHROM"):
                col_names = line.lstrip("#").split("\t")
                continue

            cols = line.split("\t")
            record = dict(zip(col_names, cols))

            chrom = record["CHROM"]
            pos   = int(record["POS"])
            ref   = record["REF"]
            alt   = record["ALT"]
            filt  = record["FILTER"]

            # FORMAT / tumor sample (first sample column after NORMAL in Mutect2 output)
            fmt_keys    = record["FORMAT"].split(":")
            sample_cols = col_names[9:]                 # all sample columns
            tumor_col   = sample_cols[0]                # tumor is first in Mutect2 paired output
            fmt_vals    = record[tumor_col].split(":")
            fmt         = dict(zip(fmt_keys, fmt_vals))

            gt  = fmt.get("GT", "./.")
            dp  = int(fmt.get("DP", 0) or 0)
            ad_raw    = fmt.get("AF", "0").split(",")   # Mutect2 uses AF, not AD
            # Prefer AD (allele depth) if available
            if "AD" in fmt:
                ad_split  = fmt["AD"].split(",")
                alt_reads = int(ad_split[1]) if len(ad_split) > 1 else 0
                total     = sum(int(x) for x in ad_split) or dp
                vaf       = round(alt_reads / total, 4) if total else 0.0
            else:
                vaf       = float(ad_raw[0]) if ad_raw[0] not in (".", "") else 0.0
                alt_reads = int(round(vaf * dp))
                total     = dp

            # VEP CSQ
            info_str = record.get("INFO", "")
            info: dict[str, str] = {}
            for token in info_str.split(";"):
                if "=" in token:
                    k, v = token.split("=", 1)
                    info[k] = v
                else:
                    info[token] = "true"

            csq_fields = _parse_csq_header(header_lines)
            csq_raw    = info.get("CSQ", "")
            anno: dict[str, str] = {}
            if csq_raw and csq_fields:
                entries = [
                    dict(zip(csq_fields, entry.split("|")))
                    for entry in csq_raw.split(",")
                ]
                anno = _canonical(entries)

            variants.append({
                "chrom":            chrom,
                "pos":              pos,
                "ref":              ref,
                "alt":              alt,
                "gene":             anno.get("SYMBOL", ""),
                "transcript":       anno.get("Feature", ""),
                "consequence":      anno.get("Consequence", "").replace("&", ", "),
                "impact":           anno.get("IMPACT", ""),
                "hgvsc":            anno.get("HGVSc", "").split(":", 1)[-1],
                "hgvsp":            anno.get("HGVSp", "").split(":", 1)[-1],
                "exon":             anno.get("EXON", ""),
                "biotype":          anno.get("BIOTYPE", ""),
                "gt":               gt,
                "depth":            total,
                "alt_reads":        alt_reads,
                "vaf":              vaf,
                "clinvar_sig":      anno.get("ClinVar_CLNSIG", "").replace("&", "/"),
                "clinvar_disease":  anno.get("ClinVar_CLNDN", "").replace("&", "/"),
                "cosmic_id":        anno.get("COSMIC_CDS", ""),
                "cosmic_count":     anno.get("COSMIC_CNT", ""),
                "gnomad_af":        anno.get("gnomADg_AF", ""),
                "sift":             anno.get("SIFT", ""),
                "polyphen":         anno.get("PolyPhen", ""),
                "filter":           filt,
            })

    return variants


rows = _parse_vcf(vcf_path)

with open(out_path, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=FIELDNAMES, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
