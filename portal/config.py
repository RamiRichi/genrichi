"""GenRichi Portal — Configuration"""

import os

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR       = os.path.dirname(os.path.abspath(__file__))
GENRICHI_DIR   = os.path.dirname(BASE_DIR)          # /home/rami/genrichi
UPLOADS_DIR    = os.path.join(BASE_DIR, "uploads")  # FASTQs uploaded via browser
RESULTS_DIR    = os.path.join(GENRICHI_DIR, "results")
WORKFLOW_DIR   = os.path.join(GENRICHI_DIR, "workflow")
CONFIG_DIR     = os.path.join(GENRICHI_DIR, "config")
DB_PATH        = os.path.join(BASE_DIR, "genrichi_orders.db")
LOG_DIR        = os.path.join(BASE_DIR, "portal_logs")

# ── Flask ─────────────────────────────────────────────────────────────────────
SECRET_KEY     = "genrichi-change-this-in-production-2026"
MAX_CONTENT_LENGTH = 20 * 1024 * 1024 * 1024  # 20 GB max upload

# ── Auth (simple single-user, upgrade to LDAP/OAuth later) ───────────────────
PORTAL_USER    = "admin"
PORTAL_PASS    = "GenRichi2026!"   # Change before going live

# ── Snakemake ─────────────────────────────────────────────────────────────────
SNAKEMAKE_CMD     = "/home/rami/miniforge3/envs/snakemake/bin/snakemake"
CONDA_PREFIX      = os.path.join(GENRICHI_DIR, ".snakemake", "conda")
CONDA_FRONTEND    = "mamba"
DEFAULT_CORES     = 8

# ── Pipelines ─────────────────────────────────────────────────────────────────
# "hotspot" and "hereditary" are intentionally NOT exposed here.
#
# hotspot (workflow/Snakefile): a Phase 5.5 audit (Step 12-C/D) found a
# deterministic, unconditional DAG failure -- align.smk's mosdepth rule
# produces {sample}.regions.bed.gz, but report.smk's generate_report rule
# expects {sample}.mosdepth.regions.bed.gz as input, so no rule satisfies
# it. align.smk's mosdepth call also still carries the
# `--quantize 0:1:10:30:100:` flag already documented elsewhere in this
# repo (docs/somatic_panel/session_notes.md) as causing a silent mosdepth
# crash. Removed from customer-facing exposure until fixed and validated.
#
# hereditary (workflow/Snakefile_hereditary): produces ACMG 2015
# pathogenicity classifications for a 25-gene hereditary panel but has no
# regression baseline and no test coverage of the classifier logic
# (workflow/scripts/acmg_classifier.py). Not modified or deleted -- only
# withheld from customer-facing execution through the portal until a
# validation baseline exists.
#
# Neither pipeline's own code was changed by this restriction.
PIPELINE_MAP = {
    "comprehensive": {
        "label":      "Somatic Comprehensive Panel",
        "snakefile":  "Snakefile_comprehensive",
        "configfile": "config/comprehensive_config.yaml",
        "description": "Full somatic workup — SNV/indel + CNV + MSI + TMB (55 cancer genes)",
    },
}

# Panel types that need a matched normal
PAIRED_PANELS = {"comprehensive"}

# ── Pricing (EUR) ─────────────────────────────────────────────────────────────
PANEL_PRICES = {
    "comprehensive": 750.00,
}
COMPANY_NAME    = "GenRichi GmbH"
COMPANY_ADDRESS = "Musterstraße 1, 10115 Berlin, Deutschland"
COMPANY_EMAIL   = "info@genrichi.de"
COMPANY_WEB     = "www.genrichi.de"
COMPANY_TAX_ID  = "DE123456789"
BANK_IBAN       = "DE89 3704 0044 0532 0130 00"
BANK_BIC        = "COBADEFFXXX"
PAYMENT_DAYS    = 30

# ── Email notifications (optional) ───────────────────────────────────────────
# Set SMTP_ENABLED = True and fill in your Gmail credentials to send emails.
# For Gmail: enable 2FA → create App Password at myaccount.google.com/apppasswords
SMTP_ENABLED  = True
SMTP_HOST     = "smtp.ionos.de"
SMTP_PORT     = 587
SMTP_USER     = "info@genrichi.de"    # your Gmail address
SMTP_PASS     = "Rami83.com"     # 16-char Gmail App Password
SMTP_FROM     = "GenRichi Portal <info@genrichi.de>"
PORTAL_URL    = "https://portal.genrichi.de"
