"""
GenRichi Portal — Background Pipeline Runner

Picks up Queued orders from the DB, writes a per-order sample sheet,
launches Snakemake in a subprocess, and updates the DB with status / logs.

Also runs the durable input-pinning worker (see "Input pinning" below), which
turns "Pinning" orders into "Queued" ones by copying and verifying each
selected SFTP file before Snakemake is ever allowed to touch it.

Usage: runs as a daemon thread inside the Flask process.
"""

import hashlib
import os
import subprocess
import threading
import time
import csv
import logging
from pathlib import Path

import models
import mailer
from sftp_paths import resolve_sftp_token, path_is_within_sftp_root
from config import (
    GENRICHI_DIR, WORKFLOW_DIR, RESULTS_DIR, LOG_DIR, PINNED_INPUTS_DIR,
    SNAKEMAKE_CMD, CONDA_PREFIX, CONDA_FRONTEND, DEFAULT_CORES,
    PIPELINE_MAP, PAIRED_PANELS,
)

logger = logging.getLogger("runner")
_lock  = threading.Lock()   # only one Snakemake run at a time per server


# ═══════════════════════════════════════════════════════════════════════════
# Input pinning
#
# An order never runs against a live path in the client-writable SFTP tree.
# Selecting a file (app.py's /order/new) only records WHICH file was chosen
# (an SFTP client + filename, re-verified against sftp_paths.resolve_sftp_token
# independently of that first check) and puts the order in status "Pinning".
# This worker is what actually copies each selected file into
# PINNED_INPUTS_DIR/<order_id>/, under the portal service account's own
# ownership -- the SFTP account can never write there -- verifies the source
# had genuinely finished uploading and did not change while being copied,
# hashes the copy, and only then lets the order become "Queued". _run_snakemake
# re-hashes every pinned file again immediately before launching Snakemake, and
# Snakemake itself is given only pinned paths, never anything under
# /srv/genrichi-sftp.
#
# Durability: all state (which order is in "Pinning", which of its inputs are
# done/pending/copying) lives in the database, not in this thread's memory. A
# portal restart loses no progress -- _reset_stuck_pinning() (called from
# start(), like _reset_stuck_orders() already was) treats any input left
# 'copying' as an interrupted, untrustworthy partial file, deletes it, and
# resets it to 'pending' so the worker below picks it back up from scratch.
# ═══════════════════════════════════════════════════════════════════════════

_pin_lock = threading.Lock()               # one order's inputs pinned at a time
_PIN_POLL_INTERVAL_SEC = 5
_STABILITY_CHECK_INTERVAL_SEC = 2
_STABILITY_MAX_ATTEMPTS = 60               # ≈2 minutes waiting for an upload to stop changing
_COPY_CHUNK_SIZE = 4 * 1024 * 1024         # 4 MiB

# Disk-space estimate for the pinning pre-flight check. NOT a flat "3x" rule --
# derived from real Phase 1 "comprehensive" runs (see docs/validation/
# HCC1395_SEQC2_v1_preliminary_concordance.md): on-disk results/<sample>/ size
# divided by that run's total tumor+normal input FASTQ size.
#   SEQC2_HCC1395_FD1: 30,941,738,307 B results / 10,574,786,978 B input ≈ 2.93x
#   SEQC2_HCC1395_LL1: 26,518,777,841 B results /  9,159,815,570 B input ≈ 2.90x
# Rounded up from the ~2.9x measured average for sample-to-sample variance
# (coverage, duplication rate) the two references may not capture.
_MEASURED_OUTPUT_INPUT_RATIO = {"comprehensive": 3.2}
# No equivalent measurement exists yet for these panels. Reusing the
# comprehensive ratio is a deliberate, labelled over-estimate (the safe
# direction for a free-space check), not a real measurement -- replace it as
# soon as either panel has completed real runs to measure from.
_UNMEASURED_PANEL_FALLBACK_RATIO = _MEASURED_OUTPUT_INPUT_RATIO["comprehensive"]
_SPACE_SAFETY_MARGIN_BYTES = 5 * 1024 ** 3  # flat 5 GB buffer, not a multiplier


def _estimated_output_bytes(panel_type: str, total_input_bytes: int) -> int:
    ratio = _MEASURED_OUTPUT_INPUT_RATIO.get(panel_type, _UNMEASURED_PANEL_FALLBACK_RATIO)
    return int(total_input_bytes * ratio)


def _check_disk_space(total_input_bytes: int, panel_type: str):
    """required = the pinned copy itself (== input size) + the estimated
    pipeline output for this panel + a small flat safety margin. Checked
    against the filesystem that holds RESULTS_DIR, where both
    PINNED_INPUTS_DIR and the eventual Snakemake results/ actually live.
    Returns (ok, needed_bytes, free_bytes)."""
    import shutil
    needed = total_input_bytes + _estimated_output_bytes(panel_type, total_input_bytes) + _SPACE_SAFETY_MARGIN_BYTES
    free = shutil.disk_usage(RESULTS_DIR).free
    return free >= needed, needed, free


def _dest_suffix(filename: str) -> str:
    """Preserve the source's compression suffix (.fastq.gz / .fq.gz / .fastq /
    .fq / a bare .gz) so fastp/bwa still auto-detect gzip correctly -- but
    nothing else about the original name is kept."""
    suf = Path(filename).suffixes[-2:]
    return "".join(suf) if suf else ".fastq.gz"


def _safe_unlink(path: str) -> None:
    try:
        os.remove(path)
    except OSError:
        pass


def _fingerprint(st: os.stat_result) -> tuple:
    """(size, mtime_ns, ctime_ns, dev, ino) -- as strong an "is this still the
    very same, unchanged file" signal as a stat()/fstat() can give without
    reading its content.

    size + nanosecond mtime/ctime alone can still coincide by chance, or on a
    filesystem with coarse timestamp resolution; (dev, ino) adds file
    *identity* -- a same-size, same-second replacement (the common SFTP
    "upload as a temp name, then rename over the final name" pattern) almost
    always gets a fresh inode even when size and both timestamps happen to
    collide, because the old and new files are two different inodes that
    merely share a name for an instant.

    Known, deliberately-not-hidden residual limitation: this is still a
    stat()-level signal, not a content checksum. A storage layer that does
    not update mtime/ctime on a content-changing write (unusual, and would
    violate normal POSIX semantics, but some non-local/virtual filesystems
    behave this way), or an in-place overwrite that happens to reuse the same
    inode with byte-identical resulting size, cannot be distinguished from "no
    change" by any stat()-based check, including this one -- only re-hashing
    the source content itself would catch that, which pinning does not do
    (see _copy_and_hash()'s own docstring for where identity is actually
    enforced against this fingerprint, and _run_snakemake()'s pre-launch
    re-hash of the *pinned copy*, which is a separate, later check)."""
    return (st.st_size, st.st_mtime_ns, st.st_ctime_ns, st.st_dev, st.st_ino)


def _wait_for_stable_source(client: str, filename: str):
    """Re-resolve the token (independently of whatever selection-time or a
    previous pinning attempt found) and poll its fingerprint (see
    _fingerprint()) until two consecutive reads agree, i.e. the upload has
    genuinely finished -- not just that a file with this name currently
    exists, and not fooled by a same-size, same-mtime replacement the way a
    plain (size, mtime) comparison could be. Returns
    (True, (path, fingerprint)) or (False, error_message)."""
    token = f"{client}/{filename}"
    path = resolve_sftp_token(token, allowed_client=client)
    if path is None:
        return False, "source file not found, or no longer authorized for its recorded SFTP account"
    try:
        prev = _fingerprint(os.stat(path))
    except OSError:
        return False, "source file disappeared before pinning could start"
    for _ in range(_STABILITY_MAX_ATTEMPTS):
        time.sleep(_STABILITY_CHECK_INTERVAL_SEC)
        path2 = resolve_sftp_token(token, allowed_client=client)  # catches a symlink swap mid-wait too
        if path2 is None:
            return False, "source file disappeared, or became unauthorized, while waiting for the upload to finish"
        try:
            cur = _fingerprint(os.stat(path2))
        except OSError:
            return False, "source file disappeared while waiting for the upload to finish"
        if cur == prev:
            return True, (path2, cur)
        prev = cur
    return False, "source file did not stop changing (the upload may still be in progress); please retry once it finishes"


def _copy_and_hash(order_id: str, field_name: str, source_path: str, source_filename: str,
                   expected_fingerprint: tuple):
    """Stream source_path -> a fresh file under PINNED_INPUTS_DIR/<order_id>/,
    hashing exactly the bytes written to the destination.

    Opened with O_NOFOLLOW: if the final path component was swapped for a
    symlink in the (necessarily nonzero) time between the stability check
    that produced expected_fingerprint and this call, the open itself fails
    (ELOOP) instead of silently following it -- a TOCTOU window that a
    path-based os.stat() re-check could not have closed, because re-stat-ing
    the same path is exactly as swappable as the original stat was.

    Identity/change checks are all done via fstat() on the one open file
    descriptor, not by re-stat()-ing the path again -- an fd, once open,
    keeps referring to the same inode no matter what later happens to the
    path (a rename or unlink of the name does not redirect or invalidate an
    already-open fd). That makes both of the following genuinely race-free
    relative to anything happening to the *path* after open() succeeds:
      1. immediately after opening, fstat() must match expected_fingerprint
         exactly -- catches a source swapped between the stability check and
         the copy starting, including a same-size/same-mtime replacement
         (different inode) that a size/mtime-only check would have missed;
      2. immediately after the read loop reaches EOF, fstat() on the SAME fd
         must *still* match expected_fingerprint -- catches genuine in-place
         modification of the same inode's content during the copy (which
         updates its size and/or mtime/ctime).
    Byte count read must also equal the recorded size, as a basic sanity
    check independent of both stat calls.

    Residual, deliberately-not-hidden limitation: an open fd cannot detect a
    rename that points the *name* elsewhere mid-copy while leaving the
    original inode (and this fd's view of it) untouched -- but that is not a
    integrity problem here, since what gets pinned and hashed is still
    exactly the stable, fingerprint-verified bytes this function was asked to
    copy, not whatever now happens to sit at that name. What cannot be fully
    eliminated, on any POSIX filesystem, is the combination of (a) a storage
    layer that fails to update mtime/ctime on a content-changing write, or
    (b) the astronomically unlikely case of a deleted-and-recreated file
    reusing the exact same inode number *and* ending up the exact same size,
    both within the single fstat-to-fstat window above -- see
    _fingerprint()'s docstring for the same point in more detail. Returns
    (True, {"path","sha256","size"}) or (False, error_message)."""
    dest_dir = os.path.join(PINNED_INPUTS_DIR, order_id)
    os.makedirs(dest_dir, exist_ok=True)
    dest_path = os.path.join(dest_dir, f"{field_name}{_dest_suffix(source_filename)}")
    _safe_unlink(dest_path)  # a previous failed attempt may have left a partial file

    o_nofollow = getattr(os, "O_NOFOLLOW", 0)  # POSIX-only; production runs on Linux (see docstring)
    try:
        fd = os.open(source_path, os.O_RDONLY | o_nofollow)
    except OSError as exc:
        return False, f"source file could not be opened safely: {exc}"

    try:
        try:
            pre = _fingerprint(os.fstat(fd))
        except OSError as exc:
            return False, f"source file could not be inspected after opening: {exc}"
        if pre != expected_fingerprint:
            return False, "source file changed between the stability check and the start of copying; please retry"

        h = hashlib.sha256()
        written = 0
        try:
            with os.fdopen(fd, "rb", closefd=False) as src, open(dest_path, "xb") as dst:
                while True:
                    chunk = src.read(_COPY_CHUNK_SIZE)
                    if not chunk:
                        break
                    dst.write(chunk)
                    h.update(chunk)
                    written += len(chunk)
        except OSError as exc:
            _safe_unlink(dest_path)
            return False, f"copy failed: {exc}"

        try:
            post = _fingerprint(os.fstat(fd))
        except OSError as exc:
            _safe_unlink(dest_path)
            return False, f"source file could not be re-inspected after copying: {exc}"
        if post != expected_fingerprint or written != expected_fingerprint[0]:
            _safe_unlink(dest_path)
            return False, "source file changed while it was being pinned; please retry"
    finally:
        os.close(fd)

    return True, {"path": dest_path, "sha256": h.hexdigest(), "size": written}


def _fail_order_pinning(order_id: str) -> None:
    """One field failed, or the pre-copy space check did -- the order can
    never run partially pinned, so undo every sibling field that DID already
    succeed too (delete its file, forget its recorded pin) before marking the
    whole order Failed."""
    rows = models.get_order_inputs(order_id)
    for row in rows:
        if row["status"] == "done" and row["pinned_path"]:
            _safe_unlink(row["pinned_path"])
            models.clear_order_input_pin(order_id, row["field_name"])
    reasons = "; ".join(f"{r['field_name']}: {r['error_msg']}" for r in models.get_order_inputs(order_id)
                        if r["error_msg"])
    models.update_status(order_id, "Failed", error_msg=f"Input pinning failed — {reasons}")
    logger.error("Order %s: pinning failed (%s)", order_id, reasons)


def _pin_order_inputs(order_id: str) -> None:
    order = models.get_order(order_id)
    if order is None:
        return
    rows = [r for r in models.get_order_inputs(order_id) if r["status"] != "done"]
    if not rows:
        # Nothing left to pin (e.g. a retry where everything was already
        # done) -- go straight to Queued.
        models.update_status(order_id, "Queued")
        return

    # Phase 1: every pending source must be a genuinely finished upload before
    # anything is copied. Collected sizes drive the space check in phase 2.
    stable = {}
    for row in rows:
        models.mark_input_copying(order_id, row["field_name"])
        ok, info = _wait_for_stable_source(row["source_client"], row["source_filename"])
        if not ok:
            models.mark_input_failed(order_id, row["field_name"], info)
            _fail_order_pinning(order_id)
            return
        stable[row["field_name"]] = info  # (path, fingerprint) -- see _fingerprint()

    # Phase 2: one combined free-space check for the whole order, before any
    # copy starts -- never begin a multi-GB copy only to run out mid-way.
    total_input_bytes = sum(fingerprint[0] for _, fingerprint in stable.values())  # fingerprint[0] == size
    ok, needed, free = _check_disk_space(total_input_bytes, order["panel_type"])
    if not ok:
        msg = (f"insufficient disk space to pin and process this order: "
              f"need ~{needed / 1e9:.1f} GB, have ~{free / 1e9:.1f} GB free")
        for row in rows:
            models.mark_input_failed(order_id, row["field_name"], msg)
        _fail_order_pinning(order_id)
        return

    # Phase 3: copy + hash each field, re-verified against its own phase-1 fingerprint.
    for row in rows:
        field = row["field_name"]
        source_path, fingerprint = stable[field]
        ok, result = _copy_and_hash(order_id, field, source_path, row["source_filename"], fingerprint)
        if not ok:
            models.mark_input_failed(order_id, field, result)
            _fail_order_pinning(order_id)
            return
        models.mark_input_done(order_id, field, result["path"], result["sha256"], result["size"])

    models.update_status(order_id, "Queued")
    logger.info("Order %s: all inputs pinned, now Queued.", order_id)


def _pinning_worker() -> None:
    """Daemon loop: pick up one Pinning order at a time and pin its inputs."""
    logger.info("Pinning worker started.")
    while True:
        try:
            pending = models.list_pinning_orders(50)
            if pending and _pin_lock.acquire(blocking=False):
                try:
                    _pin_order_inputs(pending[0]["order_id"])
                finally:
                    _pin_lock.release()
        except Exception:
            logger.exception("Pinning worker loop error")
        time.sleep(_PIN_POLL_INTERVAL_SEC)


def _reset_stuck_pinning() -> None:
    """On startup: any order_inputs row left 'copying' means the portal
    process died mid-copy. Never trust a partial file left on disk from
    before the crash -- delete anything for that field (whatever extension it
    was using) and reset the row to 'pending' so the worker retries it
    cleanly from scratch."""
    try:
        stuck = models.stuck_copying_inputs()
        for row in stuck:
            dest_dir = os.path.join(PINNED_INPUTS_DIR, row["order_id"])
            if os.path.isdir(dest_dir):
                for f in os.listdir(dest_dir):
                    if f.startswith(row["field_name"] + "."):
                        _safe_unlink(os.path.join(dest_dir, f))
            models.reset_input_to_pending(row["order_id"], row["field_name"])
            logger.warning("Reset stuck pinning input %s/%s -> pending",
                           row["order_id"], row["field_name"])
        if stuck:
            logger.info("Reset %d stuck pinning input(s) after unexpected restart.", len(stuck))
    except Exception:
        logger.exception("Error resetting stuck pinning inputs")


def _verify_pinned_inputs(order_id: str):
    """Re-hash every input this order has recorded as pinned, right before
    Snakemake is allowed to start, independent of the hash computed when it
    was originally pinned (which could in principle have been altered
    afterward by anything else with access to PINNED_INPUTS_DIR). Returns
    (True, "") or (False, reason)."""
    for row in models.get_order_inputs(order_id):
        if row["status"] != "done":
            continue
        path = row["pinned_path"]
        try:
            h = hashlib.sha256()
            with open(path, "rb") as f:
                while True:
                    chunk = f.read(_COPY_CHUNK_SIZE)
                    if not chunk:
                        break
                    h.update(chunk)
        except OSError as exc:
            return False, f"{row['field_name']}: pinned copy missing or unreadable ({exc})"
        if h.hexdigest() != row["sha256"]:
            return False, f"{row['field_name']}: pinned copy no longer matches its recorded SHA-256"
    return True, ""


def _check_no_untracked_live_sftp_paths(order: dict) -> tuple:
    """Refuse to run if a non-empty orders.fastq_* value has no matching
    'done' order_inputs row for that field AND resolves inside the live,
    client-writable SFTP tree.

    This is the gap a pre-pinning-era order can fall into: it was created
    before this whole mechanism existed, so it was never routed through
    'Pinning' and has zero order_inputs rows -- _verify_pinned_inputs() above
    has nothing to re-hash for it, because there is nothing recorded to
    re-hash. Silently trusting such a value just because the column happens
    to be non-empty would let the pipeline read directly from a path the
    SFTP client can still write to at any moment, exactly what pinning exists
    to prevent -- an untracked field is not automatically safe, it is
    automatically unverified.

    A genuinely input-free legacy order, or one whose fastq_* points
    somewhere else entirely (e.g. config.UPLOADS_DIR, the portal's own local
    upload area from before SFTP selection existed), is unaffected: this only
    rejects a value that is actually inside SFTP_ROOT right now. Nothing here
    ever rewrites or deletes the order's own data -- an order this rejects
    can still be inspected, and a fresh order can be created against the
    same/updated source through the normal pinned path.

    Returns (True, "") or (False, reason)."""
    tracked_done = {r["field_name"] for r in models.get_order_inputs(order["order_id"])
                    if r["status"] == "done"}
    for field in models.ORDER_FASTQ_FIELDS:
        value = order.get(field) or ""
        if not value or field in tracked_done:
            continue
        if path_is_within_sftp_root(value):
            return False, (f"{field}: this order predates input pinning and its recorded path is "
                          "still inside the live SFTP upload tree, which is never trusted directly "
                          "-- create a new order and reselect the file so it is pinned and verified")
    return True, ""


# ═══════════════════════════════════════════════════════════════════════════
# Snakemake execution
# ═══════════════════════════════════════════════════════════════════════════

def _write_sample_sheet(order: dict, panel_type: str) -> str:
    """Write a per-order samples TSV and return its path.

    order["fastq_r1"] etc. are, by the time an order ever reaches this
    function, either empty (no file was selected for that slot) or a pinned
    path under PINNED_INPUTS_DIR -- never a live path under /srv/genrichi-sftp.
    Snakemake reads only from here.
    """
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

    # Last-mile integrity gate: re-hash every pinned input before Snakemake
    # ever runs. This is independent of, and in addition to, the hash taken
    # when the file was pinned -- see _verify_pinned_inputs()'s docstring.
    ok, reason = _verify_pinned_inputs(order_id)
    # Second, separate gate: catch a pre-pinning-era order whose fastq_* still
    # points directly at the live SFTP tree with no order_inputs row backing
    # it at all -- see _check_no_untracked_live_sftp_paths()'s docstring.
    if ok:
        ok, reason = _check_no_untracked_live_sftp_paths(order)
    if not ok:
        models.update_status(order_id, "Failed", error_msg=f"Input integrity check failed: {reason}")
        logger.error("Order %s: input integrity check failed: %s", order_id, reason)
        return

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
            # Send a generic status-change notice -- no report link or order
            # data goes into the email itself (see mailer.send_completion_email).
            updated_order = dict(models.get_order(order_id))
            mailer.send_completion_email(updated_order)
        else:
            models.update_status(order_id, "Failed",
                                 error_msg=f"Snakemake exit code {proc.returncode}")
            logger.error("Order %s FAILED (exit %s)", order_id, proc.returncode)
            # Send a generic status-change notice -- no error text goes into
            # the email itself (see mailer.send_completion_email).
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
    """Start the background runner threads."""
    _reset_stuck_orders()
    _reset_stuck_pinning()
    t = threading.Thread(target=_worker, daemon=True, name="pipeline-runner")
    t.start()
    p = threading.Thread(target=_pinning_worker, daemon=True, name="input-pinning")
    p.start()
    logger.info("Runner and pinning threads started.")
