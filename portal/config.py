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
PIPELINE_MAP = {
    "hotspot": {
        "label":      "Somatic Hotspot Panel",
        "snakefile":  "Snakefile",
        "configfile": "config/config.yaml",
        "description": "Quick driver mutation screen — KRAS, BRAF, TP53, PIK3CA, EGFR …",
    },
    "hereditary": {
        "label":      "Hereditary Germline Panel",
        "snakefile":  "Snakefile_hereditary",
        "configfile": "config/hereditary_config.yaml",
        "description": "ACMG-classified germline variants — BRCA1/2, ATM, PALB2, MLH1 … (25 genes)",
    },
    "comprehensive": {
        "label":      "Somatic Comprehensive Panel",
        "snakefile":  "Snakefile_comprehensive",
        "configfile": "config/comprehensive_config.yaml",
        "description": "Full somatic workup — SNV/indel + CNV + MSI + TMB (65 cancer genes)",
    },
}

# Panel types that need a matched normal
PAIRED_PANELS = {"comprehensive"}

# ── Pricing (EUR) ─────────────────────────────────────────────────────────────
PANEL_PRICES = {
    "hotspot":       350.00,
    "hereditary":    490.00,
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
