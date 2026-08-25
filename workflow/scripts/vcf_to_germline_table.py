"""Parse VEP-annotated germline VCF and emit a flat TSV for ACMG classification."""

import csv
import gzip
import re

vcf_path: str = snakemake.input.vcf   # type: ignore[name-defined]
out_path: str = snakemake.output.tsv  # type: ignore[name-defined]

FIELDNAMES = [
    "chrom", "pos", "ref", "alt",
    "gene", "transcript", "consequence", "impact",
    "hgvsc", "hgvsp", "exon", "biotype",
    "gt", "zygosity", "depth", "alt_reads", "vaf", "gq",
    "clinvar_sig", "clinvar_disease",
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


def _zygosity(gt: str) -> str:
    alleles = [a for a in re.split(r"[/|]", gt) if a not in (".", "")]
    if not alleles:
        return "Unknown"
    unique = set(alleles)
    if "0" not in unique:
        return "Homozygous"
    if len(unique) > 1:
        return "Heterozygous"
    return "Ref/Ref"


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

            # Core fields
            chrom = record["CHROM"]
            pos   = int(record["POS"])
            ref   = record["REF"]
            alt   = record["ALT"]
            filt  = record["FILTER"]

            # FORMAT / sample fields
            fmt_keys   = record["FORMAT"].split(":")
            sample_col = col_names[-1]
            fmt_vals   = record[sample_col].split(":")
            fmt        = dict(zip(fmt_keys, fmt_vals))

            gt  = fmt.get("GT", "./.")
            gq  = fmt.get("GQ", ".")
            dp  = int(fmt.get("DP", 0) or 0)
            ad_raw    = fmt.get("AD", "0,0").split(",")
            alt_reads = int(ad_raw[1]) if len(ad_raw) > 1 else 0
            total     = int(ad_raw[0]) + alt_reads if ad_raw else dp
            vaf       = round(alt_reads / total, 4) if total else 0.0

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
                "chrom":           chrom,
                "pos":             pos,
                "ref":             ref,
                "alt":             alt,
                "gene":            anno.get("SYMBOL", ""),
                "transcript":      anno.get("Feature", ""),
                "consequence":     anno.get("Consequence", "").replace("&", ", "),
                "impact":          anno.get("IMPACT", ""),
                "hgvsc":           anno.get("HGVSc", "").split(":", 1)[-1],
                "hgvsp":           anno.get("HGVSp", "").split(":", 1)[-1],
                "exon":            anno.get("EXON", ""),
                "biotype":         anno.get("BIOTYPE", ""),
                "gt":              gt,
                "zygosity":        _zygosity(gt),
                "depth":           total,
                "alt_reads":       alt_reads,
                "vaf":             vaf,
                "gq":              gq,
                "clinvar_sig":     anno.get("ClinVar_CLNSIG", "").replace("&", "/"),
                "clinvar_disease": anno.get("ClinVar_CLNDN", "").replace("&", "/"),
                "gnomad_af":       anno.get("gnomADg_AF", ""),
                "sift":            anno.get("SIFT", ""),
                "polyphen":        anno.get("PolyPhen", ""),
                "filter":          filt,
            })

    return variants


rows = _parse_vcf(vcf_path)

with open(out_path, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=FIELDNAMES, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
