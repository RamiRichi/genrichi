"""GenRichi Portal — Configuration"""

import json
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
INSTANCE_DIR   = os.path.join(BASE_DIR, "instance")             # gitignored
INSTANCE_SETTINGS_PATH = os.path.join(INSTANCE_DIR, "settings.json")


# ── Environment / secrets ──────────────────────────────────────────────────────
# Real credentials never live in this tracked file. They are read from the
# process environment, optionally populated from portal/.env (gitignored --
# see portal/.env.example for the expected variable names; that file must
# never contain a real value). A real shell `export` or a systemd
# `Environment=` directive works the same way and always takes precedence
# over portal/.env, since dotenv values only fill in names that aren't
# already set.
def _load_dotenv(path):
    """Minimal, dependency-free .env loader (KEY=VALUE per line, '#'
    comments and blank lines ignored). Never overrides a variable already
    present in the real process environment."""
    if not os.path.isfile(path):
        return
    with open(path, "r", encoding="utf-8") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, _, value = line.partition("=")
            key = key.strip()
            value = value.strip()
            if len(value) >= 2 and value[0] == value[-1] and value[0] in "\"'":
                value = value[1:-1]
            os.environ.setdefault(key, value)


_load_dotenv(os.path.join(BASE_DIR, ".env"))


def _load_instance_settings(path):
    """Load non-secret runtime settings the portal's own Settings page has
    persisted (currently: SMTP_ENABLED, SMTP_USER). JSON, gitignored
    (portal/instance/ -- Flask's own "instance folder" convention, already
    reserved for exactly this in .gitignore before Stage 2). This file
    NEVER holds a secret -- SMTP_PASS is never read from or written to it;
    it remains governed solely by the environment / portal/.env mechanism
    from Stage 1. Returns {} if the file is absent or unreadable, so a
    fresh checkout with no instance/ directory behaves identically to
    before this mechanism existed."""
    if not os.path.isfile(path):
        return {}
    try:
        with open(path, "r", encoding="utf-8") as fh:
            data = json.load(fh)
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


_INSTANCE_SETTINGS = _load_instance_settings(INSTANCE_SETTINGS_PATH)


def _require_env(name):
    """Read a required secret from the environment. Raises RuntimeError
    (never a silent placeholder/default) if it is unset or empty. The
    error message names the missing variable only -- never its value."""
    value = os.environ.get(name)
    if not value:
        raise RuntimeError(
            f"{name} is not set. The portal will not start without it. "
            f"Set it in the process environment or in portal/.env "
            f"(gitignored, never committed) -- see portal/.env.example "
            f"for the expected variable names."
        )
    return value


# ── Flask ─────────────────────────────────────────────────────────────────────
SECRET_KEY     = _require_env("SECRET_KEY")
MAX_CONTENT_LENGTH = 20 * 1024 * 1024 * 1024  # 20 GB max upload

# ── Auth (simple single-user, upgrade to LDAP/OAuth later) ───────────────────
PORTAL_USER    = os.environ.get("PORTAL_USER", "admin")
PORTAL_PASS    = _require_env("PORTAL_PASS")

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
# SMTP_ENABLED and SMTP_USER are plain settings (not secrets). Precedence,
# highest first: a real environment variable, then a value the portal's own
# Settings page has persisted to the gitignored portal/instance/settings.json
# (Stage 2 -- see _load_instance_settings above), then the hardcoded default
# that used to be this file's only value (preserving prior behavior for a
# fresh checkout with neither of the other two present).
#
# SMTP_PASS remains exclusively environment/.env-sourced (Stage 1) and is
# NEVER read from or written to settings.json -- the Settings page cannot
# change it; only required, and only validated at startup, when SMTP is
# actually enabled, so the portal still starts cleanly with email unconfigured
# (set SMTP_ENABLED=false in portal/.env to opt out).
SMTP_ENABLED  = os.environ.get(
    "SMTP_ENABLED", str(_INSTANCE_SETTINGS.get("SMTP_ENABLED", "true"))
).strip().lower() in ("1", "true", "yes", "on")
SMTP_HOST     = "smtp.ionos.de"
SMTP_PORT     = 587
SMTP_USER     = os.environ.get("SMTP_USER", _INSTANCE_SETTINGS.get("SMTP_USER") or "info@genrichi.de")
SMTP_PASS     = _require_env("SMTP_PASS") if SMTP_ENABLED else os.environ.get("SMTP_PASS", "")
SMTP_FROM     = f"GenRichi Portal <{SMTP_USER}>"
PORTAL_URL    = "https://portal.genrichi.de"
