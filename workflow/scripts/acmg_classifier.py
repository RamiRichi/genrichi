"""
ACMG/AMP 2015 variant classification for hereditary cancer panel variants.

Criteria implemented from available VEP fields:
  PVS1  — predicted null in LOF-intolerant gene (stop/frameshift/splice)
  PS1   — ClinVar: same amino-acid change classified Pathogenic
  PM1   — missense in well-established functional domain (gene-level proxy)
  PM2   — absent / rare in gnomAD (AF < pm2_af threshold)
  PP3   — SIFT deleterious + PolyPhen probably_damaging
  PP5   — ClinVar Pathogenic/LP in a reputable database
  BA1   — gnomAD AF > ba1_af (5 %)  → Benign Stand-Alone
  BS1   — gnomAD AF > bs1_af (1 %)  → Benign Strong
  BS2   — ClinVar Benign in a reputable database (strong)
  BP4   — SIFT tolerated + PolyPhen benign
  BP6   — ClinVar Benign/LB in a reputable database (supporting)

Classification logic follows Table 5 of Richards et al. 2015 (PMID 25741868).

Additional rules applied here:
  - Variants outside the panel gene list (lof_intolerant_genes) are dropped.
  - ClinVar B/LB overrides weak computational evidence (PP/PM only) when no
    strong pathogenic criterion (PVS1 or PS-level) is present.  Expert-curated
    database evidence outweighs in-silico prediction at supporting/moderate level.
  - ClinVar compound entries ("Benign/Likely_benign") are split on "/" so that
    the most confident interpretation (Benign) drives the criterion assigned.
"""

import csv

sm = snakemake  # type: ignore[name-defined]

BA1_AF = float(sm.params.ba1_af)
BS1_AF = float(sm.params.bs1_af)
PM2_AF = float(sm.params.pm2_af)
LOF_INTOLERANT: set[str] = set(sm.params.lof_intolerant_genes)
PANEL_GENES: set[str] = LOF_INTOLERANT   # same list for this 25-gene panel

LOF_CONSEQUENCES = {
    "stop_gained", "frameshift_variant",
    "splice_donor_variant", "splice_acceptor_variant",
    "start_lost", "stop_lost", "transcript_ablation",
}


def _safe_float(val: str):
    try:
        return float(val) if val and val not in (".", "") else None
    except ValueError:
        return None


def _parse_clnsig(raw: str) -> tuple[bool, bool, bool, bool]:
    """
    Parse ClinVar significance string into (is_path, is_lp, is_ben, is_lb).

    ClinVar uses "/" to join multiple submitter classifications, e.g.:
        "Benign/Likely_benign"
        "Pathogenic/Likely_pathogenic"
    We split on "/" and evaluate each part, then take the most confident
    interpretation: Benign wins over Likely_benign, Pathogenic over LP.
    """
    if not raw:
        return False, False, False, False

    norm = raw.lower().replace(" ", "_")

    # Bail out on conflicting / uncertain entries — treat as no evidence
    if "conflicting" in norm or "uncertain" in norm:
        return False, False, False, False

    parts = {p.strip() for p in norm.split("/")}

    is_path = "pathogenic" in parts          # standalone "pathogenic"
    is_lp   = "likely_pathogenic" in parts and not is_path
    is_ben  = "benign" in parts              # "benign" (covers "benign/likely_benign")
    is_lb   = "likely_benign" in parts and not is_ben

    return is_path, is_lp, is_ben, is_lb


def _collect_criteria(row: dict) -> list[str]:
    criteria: list[str] = []
    gene        = row.get("gene", "")
    consequence = row.get("consequence", "").lower()
    impact      = row.get("impact", "")
    clnsig_raw  = row.get("clinvar_sig", "")
    sift        = row.get("sift", "").lower()
    polyphen    = row.get("polyphen", "").lower()
    af          = _safe_float(row.get("gnomad_af", ""))

    is_path, is_lp, is_ben, is_lb = _parse_clnsig(clnsig_raw)

    # ── Population frequency ────────────────────────────────────────────
    if af is not None and af > BA1_AF:
        criteria.append("BA1")
    elif af is not None and af > BS1_AF:
        criteria.append("BS1")
    elif af is None or af == 0.0 or (af is not None and af < PM2_AF):
        criteria.append("PM2")

    # ── ClinVar evidence ────────────────────────────────────────────────
    if is_path:
        criteria.append("PS1")
        criteria.append("PP5")
    elif is_lp:
        criteria.append("PP5")
    elif is_ben:
        # "Benign" or "Benign/Likely_benign" → strong benign
        criteria.append("BS2")
        criteria.append("BP6")
    elif is_lb:
        # Pure "Likely_benign" → supporting benign
        criteria.append("BP6")

    # ── Loss-of-function (PVS1) ─────────────────────────────────────────
    is_lof = any(c in consequence for c in LOF_CONSEQUENCES)
    if is_lof and gene in LOF_INTOLERANT:
        criteria.append("PVS1")

    # ── Computational evidence ──────────────────────────────────────────
    if "deleterious" in sift and "probably_damaging" in polyphen:
        criteria.append("PP3")
    elif "tolerated" in sift and "benign" in polyphen:
        criteria.append("BP4")

    # ── Missense in critical domain (PM1, gene-level proxy) ─────────────
    if "missense_variant" in consequence and gene in LOF_INTOLERANT and impact == "MODERATE":
        criteria.append("PM1")

    return criteria


def _apply_rules(criteria: list[str]) -> tuple[str, int]:
    pvs1 = "PVS1" in criteria
    ps   = sum(1 for c in criteria if c.startswith("PS"))
    pm   = sum(1 for c in criteria if c.startswith("PM"))
    pp   = sum(1 for c in criteria if c.startswith("PP"))
    ba1  = "BA1" in criteria
    bs   = sum(1 for c in criteria if c.startswith("BS"))
    bp   = sum(1 for c in criteria if c.startswith("BP"))

    # ── Benign stand-alone ──────────────────────────────────────────────
    if ba1:
        return "Benign", 1

    # ── ClinVar override ────────────────────────────────────────────────
    # When no strong pathogenic evidence exists (PVS1 or any PS criterion),
    # expert-curated ClinVar B/LB classification takes precedence over
    # moderate/supporting computational signals (PM, PP only).
    has_strong_path = pvs1 or ps >= 1
    if not has_strong_path:
        if "BS2" in criteria:   # ClinVar Benign or Benign/Likely_benign
            return "Benign", 1
        if "BP6" in criteria:   # ClinVar Likely_benign
            return "Likely_Benign", 2

    # ── Pathogenic ──────────────────────────────────────────────────────
    if pvs1 and ps >= 1:
        return "Pathogenic", 5
    if ps >= 2:
        return "Pathogenic", 5
    if ps >= 1 and pm >= 3:
        return "Pathogenic", 5
    if ps >= 1 and pm >= 2 and pp >= 2:
        return "Pathogenic", 5
    if ps >= 1 and pm >= 1 and pp >= 4:
        return "Pathogenic", 5

    # ── Likely Pathogenic ───────────────────────────────────────────────
    if pvs1 and pm >= 1:
        return "Likely_Pathogenic", 4
    if pvs1 and pp >= 2:
        return "Likely_Pathogenic", 4
    if ps >= 1 and 1 <= pm <= 2:
        return "Likely_Pathogenic", 4
    if ps >= 1 and pp >= 2:
        return "Likely_Pathogenic", 4
    if pm >= 3:
        return "Likely_Pathogenic", 4
    if pm >= 2 and pp >= 2:
        return "Likely_Pathogenic", 4
    if pm >= 1 and pp >= 4:
        return "Likely_Pathogenic", 4

    # ── Benign (two strong) ─────────────────────────────────────────────
    if bs >= 2:
        return "Benign", 1

    # ── Likely Benign ───────────────────────────────────────────────────
    if bs >= 1 and bp >= 1:
        return "Likely_Benign", 2
    if bp >= 2:
        return "Likely_Benign", 2

    return "VUS", 3


def classify(row: dict) -> tuple[str, int, str]:
    criteria = _collect_criteria(row)
    acmg_class, acmg_code = _apply_rules(criteria)
    return acmg_class, acmg_code, ",".join(criteria) if criteria else "."


# ── Main ─────────────────────────────────────────────────────────────────────
with open(sm.input.tsv, newline="") as fin, open(sm.output.tsv, "w", newline="") as fout:
    reader = csv.DictReader(fin, delimiter="\t")
    out_fields = (reader.fieldnames or []) + ["acmg_class", "acmg_code", "acmg_criteria"]
    writer = csv.DictWriter(fout, fieldnames=out_fields, delimiter="\t")
    writer.writeheader()
    for row in reader:
        # Drop variants outside the 25-gene hereditary panel
        if row.get("gene", "") not in PANEL_GENES:
            continue
        acmg_class, acmg_code, acmg_criteria = classify(row)
        row["acmg_class"]    = acmg_class
        row["acmg_code"]     = acmg_code
        row["acmg_criteria"] = acmg_criteria
        writer.writerow(row)
