"""
GenRichi Portal — Background Pipeline Runner

Picks up Queued orders from the DB, writes a per-order sample sheet,
launches Snakemake in a subprocess, and updates the DB with status / logs.

Usage: runs as a daemon thread inside the Flask process.
"""

import os
import subprocess
import threading
import time
import csv
import logging
from pathlib import Path

import models
import mailer
from config import (
    GENRICHI_DIR, WORKFLOW_DIR, RESULTS_DIR, LOG_DIR,
    SNAKEMAKE_CMD, CONDA_PREFIX, CONDA_FRONTEND, DEFAULT_CORES,
    PIPELINE_MAP, PAIRED_PANELS, PORTAL_URL
)

logger = logging.getLogger("runner")
_lock  = threading.Lock()   # only one pipeline at a time per server


def _write_sample_sheet(order: dict, panel_type: str) -> str:
    """Write a per-order samples TSV and return its path."""
    sample_id = order["order_id"].replace("-", "_")
    cfg_dir   = os.path.join(GENRICHI_DIR, "config")
    os.makedirs(cfg_dir, exist_ok=True)
    tsv_path  = os.path.join(cfg_dir, f"{sample_id}_samples.tsv")

    if panel_type == "hotspot":
        headers = ["sample_id", "fastq_r1", "fastq_r2"]
        row     = [sample_id, order["fastq_r1"], order["fastq_r2"]]

    elif panel_type == "hereditary":
        headers = ["sample_id", "fastq_r1", "fastq_r2",
                   "patient_id", "sex", "indication"]
        row     = [sample_id,
                   order["fastq_r1"], order["fastq_r2"],
                   order["patient_id"], order["sex"] or "Unknown",
                   order["tumor_type"] or "Hereditary Cancer Panel"]

    elif panel_type == "comprehensive":
        headers = ["sample_id", "tumor_r1", "tumor_r2",
                   "normal_r1", "normal_r2",
                   "patient_id", "sex", "tumor_type"]
        row     = [sample_id,
                   order["fastq_r1"],        order["fastq_r2"],
                   order["fastq_normal_r1"], order["fastq_normal_r2"],
                   order["patient_id"],      order["sex"] or "Unknown",
                   order["tumor_type"] or "Unknown"]
    else:
        raise ValueError(f"Unknown panel_type: {panel_type}")

    with open(tsv_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(headers)
        writer.writerow(row)

    return tsv_path


def _run_snakemake(order_id: str):
    order     = dict(models.get_order(order_id))
    panel     = order["panel_type"]
    pipeline  = PIPELINE_MAP[panel]
    sample_id = order_id.replace("-", "_")

    os.makedirs(LOG_DIR, exist_ok=True)
    log_path  = os.path.join(LOG_DIR, f"{order_id}.log")

    models.update_status(order_id, "Running", log_path=log_path)

    try:
        # Write sample sheet for this order
        tsv_path  = _write_sample_sheet(order, panel)

        # Patch the config to point at this order's sample sheet
        # We pass it via --config samples=... override
        snakefile = os.path.join(WORKFLOW_DIR, pipeline["snakefile"])
        configfile = os.path.join(GENRICHI_DIR, pipeline["configfile"])

        # Tools (bwa, samtools, bcftools, mosdepth, gatk4) are installed
        # system-wide (apt + snakemake env) — no --use-conda needed.
        snakemake_args = " ".join([
            "--snakefile",  snakefile,
            "--configfile", configfile,
            "--config",     f"samples={tsv_path}",
            "--rerun-incomplete",
            "--cores",      str(DEFAULT_CORES),
            "--printshellcmds",
        ])

        CONDA_BASE    = "/home/rami/miniforge3"
        SNAKEMAKE_BIN = f"{CONDA_BASE}/envs/snakemake/bin"
        VEP_DIR       = "/home/rami/ensembl-vep"
        VEP_PERL5LIB  = (
            f"{VEP_DIR}/ensembl/modules:"
            f"{VEP_DIR}/ensembl-variation/modules:"
            f"{VEP_DIR}/ensembl-funcgen/modules:"
            f"{VEP_DIR}/ensembl-io/modules:"
            "/home/rami/perl5/lib/perl5:"
            "/home/rami/perl5/lib/perl5/x86_64-linux-gnu-thread-multi"
        )
        conda_init    = (
            f"source {CONDA_BASE}/etc/profile.d/conda.sh && "
            f"conda activate snakemake && "
            f"export PATH={VEP_DIR}:{SNAKEMAKE_BIN}:/usr/local/bin:/usr/bin:/bin:$PATH && "
            f"export PERL5LIB={VEP_PERL5LIB}:$PERL5LIB"
        )
        shell_cmd = f"{conda_init} && {SNAKEMAKE_CMD} {snakemake_args}"

        logger.info("Starting (via bash): %s", shell_cmd)

        with open(log_path, "w") as lf:
            process = subprocess.Popen(
                ["bash", "-c", shell_cmd],
                cwd=GENRICHI_DIR,
                stdout=lf,
                stderr=subprocess.STDOUT,
                text=True,
            )
            # Store PID so the user can cancel the order
            models.update_status(order_id, "Running", pid=process.pid)
            proc_result = process.wait()

        class _R:
            returncode = proc_result
        proc = _R()

        if proc.returncode == 0:
            report = _find_report(sample_id, panel)
            models.update_status(order_id, "Done", report_path=report or "")
            logger.info("Order %s completed. Report: %s", order_id, report)
            # Send completion email
            updated_order = dict(models.get_order(order_id))
            report_url = f"{PORTAL_URL}/order/{order_id}/report" if report else ""
            mailer.send_completion_email(updated_order, report_url)
        else:
            models.update_status(order_id, "Failed",
                                 error_msg=f"Snakemake exit code {proc.returncode}")
            logger.error("Order %s FAILED (exit %s)", order_id, proc.returncode)
            # Send failure email
            updated_order = dict(models.get_order(order_id))
            mailer.send_completion_email(updated_order)

    except Exception as exc:
        models.update_status(order_id, "Failed", error_msg=str(exc))
        logger.exception("Runner exception for %s", order_id)


def _find_report(sample_id: str, panel: str) -> str:
    """Return path to the generated HTML report, or empty string."""
    candidates = {
        "hotspot":      f"results/{sample_id}/report/{sample_id}_report.html",
        "hereditary":   f"results/{sample_id}/report/{sample_id}_hereditary_report.html",
        "comprehensive":f"results/{sample_id}/report/{sample_id}_comprehensive_report.html",
    }
    rel = candidates.get(panel, "")
    full = os.path.join(GENRICHI_DIR, rel)
    return full if os.path.isfile(full) else ""


def _worker():
    """Daemon loop: pick up one Queued order at a time and run it."""
    logger.info("Runner daemon started.")
    while True:
        try:
            orders = models.list_orders(50)
            queued = [o for o in orders if o["status"] == "Queued"]
            if queued and _lock.acquire(blocking=False):
                try:
                    _run_snakemake(queued[0]["order_id"])
                finally:
                    _lock.release()
        except Exception:
            logger.exception("Worker loop error")
        time.sleep(10)


def _reset_stuck_orders():
    """On startup: any order stuck in 'Running' was interrupted — reset to Queued."""
    try:
        orders = models.list_orders(200)
        stuck  = [o for o in orders if o["status"] == "Running"]
        for o in stuck:
            models.update_status(o["order_id"], "Queued",
                                 error_msg="Reset after unexpected restart")
            logger.warning("Reset stuck order %s → Queued", o["order_id"])
        # Also unlock Snakemake directory
        os.system(f"rm -rf {GENRICHI_DIR}/.snakemake/locks/")
        if stuck:
            logger.info("Unlocked Snakemake directory.")
    except Exception:
        logger.exception("Error resetting stuck orders")


def start():
    """Start the background runner thread."""
    _reset_stuck_orders()
    t = threading.Thread(target=_worker, daemon=True, name="pipeline-runner")
    t.start()
    logger.info("Runner thread started.")
