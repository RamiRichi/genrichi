"""
GenRichi Portal — Flask Web Application
Multi-user edition (admin + lab_staff roles)
"""

import json
import os
import sys
import logging
from functools import wraps
from pathlib import Path

from flask import (
    Flask, render_template, request, redirect, url_for,
    session, flash, send_file, abort, jsonify
)
from werkzeug.middleware.proxy_fix import ProxyFix

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import config as cfg
import models
import runner
from models import SFTP_ACCOUNT_RE

# ── Logging ───────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
)
logger = logging.getLogger("portal")

# ── Flask app ─────────────────────────────────────────────────────────────────
app = Flask(__name__)
app.wsgi_app = ProxyFix(app.wsgi_app, x_for=1, x_proto=1, x_host=1, x_prefix=1)
app.secret_key = cfg.SECRET_KEY
app.config["MAX_CONTENT_LENGTH"] = cfg.MAX_CONTENT_LENGTH

os.makedirs(cfg.UPLOADS_DIR, exist_ok=True)
os.makedirs(cfg.LOG_DIR,     exist_ok=True)
os.makedirs(cfg.PINNED_INPUTS_DIR, exist_ok=True)


# ── Auth decorators ───────────────────────────────────────────────────────────
def login_required(f):
    @wraps(f)
    def decorated(*args, **kwargs):
        if not session.get("logged_in"):
            return redirect(url_for("login", next=request.url))
        # Old session missing role — force re-login
        if "role" not in session:
            session.clear()
            return redirect(url_for("login"))
        return f(*args, **kwargs)
    return decorated


def admin_required(f):
    @wraps(f)
    def decorated(*args, **kwargs):
        if not session.get("logged_in"):
            return redirect(url_for("login"))
        if session.get("role") != "admin":
            flash("Access denied — admin only.", "error")
            return redirect(url_for("dashboard"))
        return f(*args, **kwargs)
    return decorated


def _session_lab_id():
    """The logged-in user's lab, read fresh from the DB (never from the request)."""
    user = models.get_user(session.get("username", ""))
    return user["lab_id"] if user else None


def _may_access_lab_data(lab_id):
    """Admins see every lab; everyone else only their own single lab.
    Used for billing data, which is deliberately lab_id-only (no legacy fallback):
    an order with no lab has no billing lab to check against."""
    if session.get("role") == "admin":
        return True
    mine = _session_lab_id()
    return mine is not None and lab_id is not None and mine == lab_id


def _may_access_order(order):
    """The single access rule for every order/report/download/status/log/action
    route. Admins see every order. A lab user sees an order if it is linked to
    their own lab (shared with every colleague at that lab) -- never another
    lab's. An order with no lab_id (legacy, created before lab-linking existed)
    falls back to creator-only access, exactly as before this fix: this keeps a
    user's own historical orders visible without ever widening a no-lab order's
    visibility to a whole lab it was never actually linked to."""
    if session.get("role") == "admin":
        return True
    if order["lab_id"] is not None:
        mine = _session_lab_id()
        return mine is not None and mine == order["lab_id"]
    return order["created_by"] == session.get("username")


# ── SFTP upload isolation ─────────────────────────────────────────────────────
# The trusted server-side link from a lab to its SFTP directory is
# labs.sftp_account (models.set_lab_sftp_account), set by an admin -- never
# inferred here from a lab name, a username, or anything a request sends. A
# lab with no sftp_account has no authorized SFTP files at all: the browser
# API returns an empty, explicitly-unlinked list for it rather than falling
# back to showing every client's files. The browser is only ever given
# 'client/filename' tokens, never a real filesystem path; every token is
# re-resolved and re-validated here again at order-creation time, independent
# of what the listing endpoint returned, so a hand-crafted form POST gets the
# same enforcement as the picker UI. Resolution itself now lives in
# sftp_paths.py (not here) so runner.py can re-run the identical check again,
# a third time, right before it copies the file into the pinned-inputs area --
# see "Input pinning" below and sftp_paths.py's own docstring.
from sftp_paths import (
    sftp_upload_root as _sftp_upload_root,
    list_sftp_clients as _list_sftp_clients,
    list_sftp_files as _list_sftp_files,
    resolve_sftp_token as _resolve_sftp_token,
)


def _resolve_order_fastq_field(raw_token, lab_row):
    """Validate one fastq_* form field for /order/new. Empty stays empty (an
    order may be created before a file is staged) -- returns (True, None).
    A non-empty value must resolve, via resolve_sftp_token, to a real file
    inside the authorized SFTP root: the caller's own lab's account for a lab
    user, or any real client for an admin.

    This is the SELECTION-time check only -- it decides whether the order may
    be created at all, and what (client, filename) gets recorded to pin. It is
    NOT the last check on these bytes: runner.py's pinning worker re-runs
    resolve_sftp_token on the exact same (client, filename) again, completely
    independently, right before it copies the file (see "Input pinning"
    below) -- this function's result is never trusted as still valid by then.

    Returns (True, None) for an empty field, (True, {"client", "name"}) for a
    resolved one, or (False, flash_message) on any failure.
    """
    raw_token = (raw_token or "").strip()
    if not raw_token:
        return True, None
    if session.get("role") == "admin":
        allowed_client = None
    else:
        allowed_client = lab_row["sftp_account"] if lab_row else None
        if not allowed_client:
            return False, ("Your lab is not yet linked to an SFTP upload account, so no "
                           "uploaded file can be selected. Ask an administrator to link it.")
    path = _resolve_sftp_token(raw_token, allowed_client=allowed_client)
    if path is None:
        return False, "The selected file was not found in your authorized SFTP uploads."
    client, _, name = raw_token.partition("/")
    return True, {"client": client, "name": name}


# ── Login / Logout ────────────────────────────────────────────────────────────
@app.route("/login", methods=["GET", "POST"])
def login():
    if request.method == "POST":
        u = request.form.get("username", "").strip().lower()
        p = request.form.get("password", "")
        user = models.verify_user(u, p)
        if user:
            session["logged_in"] = True
            session["username"]  = user["username"]
            session["role"]      = user["role"]
            session["full_name"] = user["full_name"] or user["username"]
            logger.info("Login: %s (%s)", u, user["role"])
            return redirect(request.args.get("next") or url_for("dashboard"))
        flash("Invalid username or password.", "error")
    return render_template("login.html")


@app.route("/logout")
def logout():
    session.clear()
    return redirect(url_for("login"))


# ── Dashboard ─────────────────────────────────────────────────────────────────
@app.route("/")
@login_required
def dashboard():
    uname = session["username"]
    role  = session["role"]
    orders = models.list_orders(200, username=uname, role=role, lab_id=_session_lab_id())
    stats = {
        "total":   len(orders),
        "pinning": sum(1 for o in orders if o["status"] == "Pinning"),
        "queued":  sum(1 for o in orders if o["status"] == "Queued"),
        "running": sum(1 for o in orders if o["status"] == "Running"),
        "done":    sum(1 for o in orders if o["status"] == "Done"),
        "failed":  sum(1 for o in orders if o["status"] == "Failed"),
    }
    return render_template("dashboard.html", orders=orders, stats=stats,
                           pipelines=cfg.PIPELINE_MAP)


# ── Statistics ────────────────────────────────────────────────────────────────
@app.route("/stats")
@login_required
def stats():
    from collections import Counter, defaultdict
    from datetime import date, timedelta

    uname  = session["username"]
    role   = session["role"]
    orders = models.list_orders(500, username=uname, role=role, lab_id=_session_lab_id())

    stats = {
        "total":   len(orders),
        "done":    sum(1 for o in orders if o["status"] == "Done"),
        "failed":  sum(1 for o in orders if o["status"] == "Failed"),
        "running": sum(1 for o in orders if o["status"] == "Running"),
        "queued":  sum(1 for o in orders if o["status"] == "Queued"),
        "pinning": sum(1 for o in orders if o["status"] == "Pinning"),
    }

    status_data = {"Done": stats["done"], "Running": stats["running"],
                   "Failed": stats["failed"], "Queued": stats["queued"]}

    panel_counts = Counter(o["panel_type"] for o in orders)
    panel_data   = {cfg.PIPELINE_MAP.get(k, {}).get("label", k): v
                    for k, v in panel_counts.items()}

    day_counts = defaultdict(int)
    for o in orders:
        day = o["created_at"][:10]
        day_counts[day] += 1
    today = date.today()
    labels, values = [], []
    for i in range(29, -1, -1):
        d = str(today - timedelta(days=i))
        labels.append(d[5:])
        values.append(day_counts.get(d, 0))
    timeline_data = {"labels": labels, "values": values}

    return render_template("stats.html", stats=stats, orders=orders,
                           status_data=status_data, panel_data=panel_data,
                           timeline_data=timeline_data, pipelines=cfg.PIPELINE_MAP)


# ── New order ─────────────────────────────────────────────────────────────────
@app.route("/order/new", methods=["GET", "POST"])
@login_required
def new_order():
    if request.method == "GET":
        return render_template("new_order.html", pipelines=cfg.PIPELINE_MAP,
                               paired_panels=list(cfg.PAIRED_PANELS),
                               labs=models.list_labs() if session["role"] == "admin" else [])

    panel_type   = request.form.get("panel_type", "comprehensive")
    patient_id   = request.form.get("patient_id", "").strip()
    patient_name = request.form.get("patient_name", "").strip()
    sex          = request.form.get("sex", "")
    dob          = request.form.get("dob", "")
    tumor_type   = request.form.get("tumor_type", "").strip()
    notes        = request.form.get("notes", "").strip()
    notify_email = request.form.get("notify_email", "").strip()

    if not patient_id:
        flash("Patient ID is required.", "error")
        return redirect(url_for("new_order"))

    r1        = request.form.get("fastq_r1",        "").strip()
    r2        = request.form.get("fastq_r2",        "").strip()
    normal_r1 = request.form.get("fastq_normal_r1", "").strip()
    normal_r2 = request.form.get("fastq_normal_r2", "").strip()

    # The billing lab is decided server-side: a lab user's order always carries
    # that user's own lab; only an admin may choose one (never taken from a
    # lab user's form).
    if session["role"] == "admin":
        raw_lab = request.form.get("lab_id", "").strip()
        lab_id = int(raw_lab) if raw_lab.isdigit() else None
        if lab_id is not None and models.get_lab(lab_id) is None:
            flash("Unknown lab.", "error")
            return redirect(url_for("new_order"))
    else:
        lab_id = _session_lab_id()
    lab_row = models.get_lab(lab_id) if lab_id is not None else None

    # Every non-empty FASTQ field must resolve, server-side, to a real file inside
    # the authorized SFTP uploads root -- for a lab user, their own lab's linked
    # account only; for an admin, any real client. This runs regardless of
    # whether the value came from the file-browser picker or was typed by hand,
    # and independent of whatever /api/sftp-files previously returned. It only
    # decides whether the order may be created and which (client, filename) get
    # recorded to pin -- it is NOT what the pipeline ends up reading from (see
    # "Input pinning" below): resolved[*] here is either None (no file for that
    # slot) or a {"client","name"} pair, never a filesystem path.
    resolved = {}
    for field, raw in (("fastq_r1", r1), ("fastq_r2", r2),
                       ("fastq_normal_r1", normal_r1), ("fastq_normal_r2", normal_r2)):
        ok, value = _resolve_order_fastq_field(raw, lab_row)
        if not ok:
            flash(value, "error")
            return redirect(url_for("new_order"))
        resolved[field] = value

    # Input pinning: the order NEVER stores a live SFTP path. It is created with
    # empty fastq_* columns and, if any field was selected, status "Pinning" --
    # runner.py's pinning worker (durable, DB-driven, survives a portal restart;
    # see runner.py's "Input pinning" section) copies each selected file into a
    # server-account-owned directory the SFTP client can never write to,
    # verifying the source was not still uploading and did not change during the
    # copy, and hashes the copy. Only once every selected field is pinned does
    # the order become "Queued"; if any fails, none of the order's pinned copies
    # are kept and the order is marked "Failed" without ever running. An order
    # with nothing selected has nothing to pin and goes straight to "Queued", as
    # before this change.
    # Order creation and (if any file was selected) every one of its
    # order_inputs rows must land together, in one transaction -- see
    # models.new_order_with_inputs()'s docstring for why a separate insert per
    # field, after the order row's own insert, would risk leaving a 'Pinning'
    # order with only some of its selected fields recorded.
    inputs = [(field, info["client"], info["name"])
             for field, info in resolved.items() if info is not None]
    needs_pinning = bool(inputs)
    if needs_pinning:
        order_id = models.new_order_with_inputs(
            patient_id=patient_id, patient_name=patient_name,
            sex=sex, dob=dob, tumor_type=tumor_type,
            panel_type=panel_type, inputs=inputs,
            notes=notes, notify_email=notify_email,
            created_by=session["username"], lab_id=lab_id,
        )
    else:
        order_id = models.new_order(
            patient_id=patient_id, patient_name=patient_name,
            sex=sex, dob=dob, tumor_type=tumor_type,
            panel_type=panel_type,
            fastq_r1="", fastq_r2="", fastq_normal_r1="", fastq_normal_r2="",
            notes=notes, notify_email=notify_email,
            created_by=session["username"], lab_id=lab_id,
            status="Queued",
        )

    if needs_pinning:
        flash(f"Order {order_id} created. Verifying and copying the selected file(s)...", "success")
    else:
        flash(f"Order {order_id} created and queued.", "success")
    return redirect(url_for("order_detail", order_id=order_id))


# ── Order detail ──────────────────────────────────────────────────────────────
@app.route("/order/<order_id>")
@login_required
def order_detail(order_id):
    order = models.get_order(order_id)
    if not order:
        abort(404)

    # Lab staff can only see their own orders
    if not _may_access_order(order):
        abort(403)

    log_tail = ""
    if order["log_path"] and os.path.isfile(order["log_path"]):
        with open(order["log_path"]) as f:
            lines    = f.readlines()
            log_tail = "".join(lines[-60:])

    return render_template("order.html", order=order, log_tail=log_tail,
                           pipeline=cfg.PIPELINE_MAP.get(order["panel_type"], {}),
                           invoice_issued=models.get_invoice(order_id) is not None,
                           labs=models.list_labs() if session["role"] == "admin" else [])


# ── Report viewer ─────────────────────────────────────────────────────────────
@app.route("/order/<order_id>/report/embed")
@login_required
def report_embed(order_id):
    order = models.get_order(order_id)
    if not order or not order["report_path"]:
        abort(404)
    if not _may_access_order(order):
        abort(403)
    if not os.path.isfile(order["report_path"]):
        abort(404)
    return render_template("report_view.html", order=order,
                           pipeline=cfg.PIPELINE_MAP.get(order["panel_type"], {}))


@app.route("/order/<order_id>/report")
@login_required
def view_report(order_id):
    order = models.get_order(order_id)
    if not order or not order["report_path"]:
        abort(404)
    if not _may_access_order(order):
        abort(403)
    if not os.path.isfile(order["report_path"]):
        abort(404)
    return send_file(order["report_path"], mimetype="text/html")


@app.route("/order/<order_id>/report/download")
@login_required
def download_report(order_id):
    order = models.get_order(order_id)
    if not order or not order["report_path"]:
        abort(404)
    if not _may_access_order(order):
        abort(403)
    if not os.path.isfile(order["report_path"]):
        abort(404)
    fname = f"{order_id}_{order['panel_type']}_report.html"
    return send_file(order["report_path"], mimetype="text/html",
                     as_attachment=True, download_name=fname)


# ── Invoice ───────────────────────────────────────────────────────────────────
# Recipient = the ordering lab (never the patient). GET never writes: without an
# issued invoice it renders a number-less admin-only draft. An invoice exists
# (with a fixed, gapless number) only after an explicit admin POST to
# /order/<id>/invoice/issue.
def _issuer():
    return {
        "name": cfg.COMPANY_NAME, "owner": cfg.COMPANY_OWNER,
        "address": cfg.COMPANY_ADDRESS, "email": cfg.COMPANY_EMAIL,
        "web": cfg.COMPANY_WEB, "tax_no": cfg.COMPANY_TAX_NO,
        "tax_note": cfg.INVOICE_TAX_NOTE,
    }


@app.route("/order/<order_id>/invoice")
@login_required
def invoice(order_id):
    from datetime import date, timedelta
    order = models.get_order(order_id)
    if not order:
        abort(404)
    if not _may_access_lab_data(order["lab_id"]):
        abort(403)

    inv = models.get_invoice(order_id)
    if inv is None and session["role"] != "admin":
        abort(404)          # lab users only ever see issued invoices

    net_price = cfg.PANEL_PRICES.get(order["panel_type"], 0.0)
    pipeline_label = cfg.PIPELINE_MAP.get(order["panel_type"], {}).get("label", order["panel_type"])
    if inv is not None:
        recipient = json.loads(inv["recipient_json"])
        issuer = json.loads(inv["issuer_json"])
        net_price = inv["net_cents"] / 100
        idate = date.fromisoformat(inv["invoice_date"])
        ctx = dict(draft=False, invoice_no=inv["invoice_no"],
                   invoice_date=idate.strftime("%d.%m.%Y"),
                   due_date=(idate + timedelta(days=cfg.PAYMENT_DAYS)).strftime("%d.%m.%Y"),
                   service_date=date.fromisoformat(inv["service_date"]).strftime("%d.%m.%Y"),
                   description=inv["description"])
    else:
        lab = models.get_lab(order["lab_id"])
        recipient = {k: lab[k] for k in models.LAB_FIELDS} if lab else None
        issuer = _issuer()
        finished = (order["finished_at"] or "")[:10]
        ctx = dict(draft=True, invoice_no=None, invoice_date=None, due_date=None,
                   service_date=None, description=pipeline_label,
                   default_service_date=finished,
                   can_issue=bool(cfg.INVOICE_ISSUING_ENABLED and recipient
                                  and order["status"] == "Done" and net_price > 0),
                   issuing_enabled=cfg.INVOICE_ISSUING_ENABLED,
                   order_done=order["status"] == "Done")

    return render_template("invoice.html", order=order, cfg=cfg, issuer=issuer,
                           recipient=recipient, net_price=net_price, gross_price=net_price,
                           **ctx)


@app.route("/order/<order_id>/invoice/issue", methods=["POST"])
@admin_required
def issue_invoice(order_id):
    from datetime import date
    order = models.get_order(order_id)
    if not order:
        abort(404)
    if not cfg.INVOICE_ISSUING_ENABLED:
        flash("Invoice issuing is disabled (set INVOICE_ISSUING_ENABLED=true "
              "once the tax treatment has been confirmed).", "error")
        return redirect(url_for("invoice", order_id=order_id))
    if order["status"] != "Done":
        flash("Only completed orders can be invoiced.", "error")
        return redirect(url_for("invoice", order_id=order_id))
    net = cfg.PANEL_PRICES.get(order["panel_type"], 0.0)
    if net <= 0:
        flash("No price is configured for this panel.", "error")
        return redirect(url_for("invoice", order_id=order_id))
    try:
        service_date = date.fromisoformat(request.form.get("service_date", "").strip())
    except ValueError:
        flash("A valid service date (Leistungsdatum) is required.", "error")
        return redirect(url_for("invoice", order_id=order_id))
    today = date.today()
    if service_date > today or service_date < date.fromisoformat(order["created_at"][:10]):
        flash("The service date must be between the order date and today.", "error")
        return redirect(url_for("invoice", order_id=order_id))
    try:
        inv = models.issue_invoice(
            order_id, issued_by=session["username"], invoice_date=today.isoformat(),
            service_date=service_date.isoformat(),
            description=cfg.PIPELINE_MAP.get(order["panel_type"], {}).get("label", order["panel_type"]),
            net_cents=round(net * 100), issuer=_issuer())
    except models.LabError as exc:
        flash(str(exc), "error")
        return redirect(url_for("invoice", order_id=order_id))
    logger.info("Invoice %s issued for %s by %s", inv["invoice_no"], order_id, session["username"])
    flash(f"Invoice {inv['invoice_no']} issued.", "success")
    return redirect(url_for("invoice", order_id=order_id))


# ── SFTP file browser API ────────────────────────────────────────────────────
# Response shape: {"files": [...], "linked": bool}. "linked" tells the picker UI
# whether the browser is actually usable ("false" means: no admin-configured
# SFTP link exists yet for this account, not "you have zero files"). Every file
# entry carries a "token" ('client/filename'), never a real filesystem path --
# see the "SFTP upload isolation" helpers above for how a token is resolved and
# re-validated server-side at order-creation time.
@app.route("/api/sftp-files")
@login_required
def api_sftp_files():
    if session["role"] == "admin":
        files = []
        for client in _list_sftp_clients():
            files.extend(_list_sftp_files(client))
        return jsonify({"files": files, "linked": True})

    lab = models.get_lab(_session_lab_id())
    account = lab["sftp_account"] if lab else None
    if not account:
        return jsonify({"files": [], "linked": False})
    return jsonify({"files": _list_sftp_files(account), "linked": True})


# ── Status API ────────────────────────────────────────────────────────────────
@app.route("/api/order/<order_id>/status")
@login_required
def api_status(order_id):
    order = models.get_order(order_id)
    if not order:
        abort(404)
    if not _may_access_order(order):
        abort(403)
    return jsonify({"status": order["status"],
                    "report_ready": bool(order["report_path"])})


# ── Log download ──────────────────────────────────────────────────────────────
@app.route("/order/<order_id>/log")
@login_required
def download_log(order_id):
    order = models.get_order(order_id)
    if not order or not order["log_path"] or not os.path.isfile(order["log_path"]):
        abort(404)
    if not _may_access_order(order):
        abort(403)
    return send_file(order["log_path"], as_attachment=True,
                     download_name=f"{order_id}.log")


# ── Retry ─────────────────────────────────────────────────────────────────────
@app.route("/order/<order_id>/retry", methods=["POST"])
@login_required
def retry_order(order_id):
    order = models.get_order(order_id)
    if not order:
        abort(404)
    if not _may_access_order(order):
        abort(403)
    if order["status"] != "Failed":
        flash("Only failed orders can be retried.", "error")
        return redirect(url_for("order_detail", order_id=order_id))
    # A failure can happen either during pinning (inputs never got copied) or
    # during Snakemake itself (inputs were already pinned fine). Retrying must
    # not skip straight to "Queued" in the first case -- the runner would then
    # try to read fastq_* columns that are still empty. If every recorded input
    # for this order is already pinned, go straight back to Queued (re-running
    # the pipeline is enough); otherwise reset the unfinished ones and route
    # back through "Pinning" so the pinning worker retries them -- against
    # their original SFTP source, which may itself no longer exist or may have
    # changed since order creation, in which case this will (correctly) fail
    # again with a clear reason rather than silently running on stale bytes.
    if models.order_inputs_all_done(order_id):
        models.update_status(order_id, "Queued", error_msg="")
        flash(f"Order {order_id} re-queued.", "success")
    else:
        models.reset_order_inputs_for_retry(order_id)
        models.update_status(order_id, "Pinning", error_msg="")
        flash(f"Order {order_id} — re-verifying and copying the selected file(s)...", "success")
    return redirect(url_for("order_detail", order_id=order_id))


# ── Cancel ────────────────────────────────────────────────────────────────────
@app.route("/order/<order_id>/cancel", methods=["POST"])
@login_required
def cancel_order(order_id):
    import signal
    order = models.get_order(order_id)
    if not order:
        abort(404)
    if not _may_access_order(order):
        abort(403)
    if order["status"] != "Running":
        flash("Only running orders can be cancelled.", "error")
        return redirect(url_for("order_detail", order_id=order_id))

    pid = order["pid"]
    if pid:
        try:
            os.kill(pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
        except Exception as e:
            logger.warning("Could not kill PID %s: %s", pid, e)

    import subprocess as sp
    sp.run(["rm", "-rf", "/home/rami/genrichi/.snakemake/locks/"], check=False)

    models.update_status(order_id, "Failed", error_msg="Cancelled by user", pid=None)
    flash(f"Order {order_id} cancelled.", "warning")
    return redirect(url_for("order_detail", order_id=order_id))


# ── User Management (admin only) ──────────────────────────────────────────────
@app.route("/admin/users")
@admin_required
def user_management():
    users = models.list_users()
    return render_template("users.html", users=users, labs=models.list_labs())


@app.route("/admin/users/create", methods=["POST"])
@admin_required
def create_user():
    username  = request.form.get("username", "").strip().lower()
    password  = request.form.get("password", "").strip()
    role      = request.form.get("role", "lab_staff")
    full_name = request.form.get("full_name", "").strip()
    email     = request.form.get("email", "").strip()

    if not username or not password:
        flash("Username and password are required.", "error")
        return redirect(url_for("user_management"))
    if len(password) < 6:
        flash("Password must be at least 6 characters.", "error")
        return redirect(url_for("user_management"))

    raw_lab = request.form.get("lab_id", "").strip()
    lab_id = int(raw_lab) if raw_lab.isdigit() and models.get_lab(int(raw_lab)) else None
    ok = models.create_user(username, password, role, full_name, email, lab_id)
    if ok:
        flash(f"User '{username}' created successfully.", "success")
        logger.info("Admin created user: %s (%s)", username, role)
    else:
        flash(f"Username '{username}' already exists.", "error")
    return redirect(url_for("user_management"))


@app.route("/admin/users/<username>/toggle", methods=["POST"])
@admin_required
def toggle_user(username):
    if username == session["username"]:
        flash("You cannot deactivate your own account.", "error")
        return redirect(url_for("user_management"))
    user = models.get_user(username)
    if not user:
        abort(404)
    new_active = 0 if user["active"] else 1
    models.update_user(username, active=new_active)
    state = "activated" if new_active else "deactivated"
    flash(f"User '{username}' {state}.", "success")
    return redirect(url_for("user_management"))


@app.route("/admin/users/<username>/reset_password", methods=["POST"])
@admin_required
def reset_user_password(username):
    if username == session["username"]:
        flash("Use Settings to change your own password.", "error")
        return redirect(url_for("user_management"))
    new_pass = request.form.get("new_pass", "").strip()
    if len(new_pass) < 6:
        flash("Password must be at least 6 characters.", "error")
        return redirect(url_for("user_management"))
    models.update_user(username, password=new_pass)
    flash(f"Password for '{username}' reset successfully.", "success")
    return redirect(url_for("user_management"))


@app.route("/admin/users/<username>/delete", methods=["POST"])
@admin_required
def delete_user(username):
    if username == session["username"]:
        flash("You cannot delete your own account.", "error")
        return redirect(url_for("user_management"))
    models.delete_user(username)
    flash(f"User '{username}' deleted.", "warning")
    return redirect(url_for("user_management"))


# ── Labs / billing profiles ───────────────────────────────────────────────────
@app.route("/admin/labs")
@admin_required
def lab_management():
    labs = models.list_labs()
    members = {l["id"]: models.lab_members(l["id"]) for l in labs}
    return render_template("labs.html", labs=labs, members=members)


def _lab_form():
    return {k: request.form.get(k, "") for k in models.LAB_FIELDS}


@app.route("/admin/labs/create", methods=["POST"])
@admin_required
def create_lab():
    try:
        models.create_lab(**_lab_form())
        flash("Lab created.", "success")
    except models.LabError as exc:
        flash(str(exc), "error")
    return redirect(url_for("lab_management"))


@app.route("/admin/labs/<int:lab_id>/edit", methods=["POST"])
@admin_required
def edit_lab(lab_id):
    if models.get_lab(lab_id) is None:
        abort(404)
    try:
        models.update_lab(lab_id, **_lab_form())
        flash("Lab updated. Already issued invoices keep their original recipient.", "success")
    except models.LabError as exc:
        flash(str(exc), "error")
    return redirect(url_for("lab_management"))


@app.route("/admin/labs/<int:lab_id>/sftp", methods=["POST"])
@admin_required
def set_lab_sftp(lab_id):
    """Link (or clear) a lab's trusted SFTP upload account. This directory name
    must be confirmed out-of-band (e.g. against the actual SFTP server config) --
    it is never guessed from the lab's legal name or any user's username."""
    if models.get_lab(lab_id) is None:
        abort(404)
    account = request.form.get("sftp_account", "").strip()
    try:
        models.set_lab_sftp_account(lab_id, account)
        flash("Lab's SFTP account updated." if account else "Lab's SFTP link cleared.", "success")
    except models.LabError as exc:
        flash(str(exc), "error")
    return redirect(url_for("lab_management"))


@app.route("/admin/users/<username>/lab", methods=["POST"])
@admin_required
def set_user_lab(username):
    raw = request.form.get("lab_id", "").strip()
    lab_id = int(raw) if raw.isdigit() else None
    try:
        models.set_user_lab(username, lab_id)
        flash(f"Lab of '{username}' updated.", "success")
    except models.LabError as exc:
        flash(str(exc), "error")
    return redirect(url_for("user_management"))


@app.route("/admin/orders/<order_id>/lab", methods=["POST"])
@admin_required
def set_order_lab(order_id):
    raw = request.form.get("lab_id", "").strip()
    try:
        models.set_order_lab(order_id, int(raw) if raw.isdigit() else None)
        flash("Order billing lab set.", "success")
    except models.LabError as exc:
        flash(str(exc), "error")
    return redirect(url_for("order_detail", order_id=order_id))


@app.route("/lab")
@login_required
def my_lab():
    """A lab user's read-only view of their own lab's billing profile.
    There is deliberately no lab id in this URL: nobody can address another lab."""
    if session["role"] == "admin":
        return redirect(url_for("lab_management"))
    lab = models.get_lab(_session_lab_id())
    members = models.lab_members(lab["id"]) if lab else []
    return render_template("lab.html", lab=lab, members=members)


# ── Settings ──────────────────────────────────────────────────────────────────
@app.route("/settings")
@login_required
def settings():
    import shutil
    try:
        total, used, free = shutil.disk_usage(cfg.RESULTS_DIR)
        disk_usage = f"{used//1024//1024//1024:.1f} GB used / {total//1024//1024//1024:.1f} GB total"
    except Exception:
        disk_usage = "N/A"

    # SMTP_PASS is never read here and never passed into the template
    # context -- not even masked. It must stay server-side only (Phase
    # 5.5 Stage 2 review). The template shows a static, non-secret status
    # line instead of any value derived from the real password.
    return render_template("settings.html",
        smtp_enabled  = cfg.SMTP_ENABLED,
        smtp_user     = cfg.SMTP_USER,
        portal_url    = cfg.PORTAL_URL,
        admin_user    = cfg.PORTAL_USER,
        db_path       = cfg.DB_PATH,
        workflow_dir  = cfg.WORKFLOW_DIR,
        total_orders  = len(models.list_orders(1000, username=session["username"],
                                               role=session["role"], lab_id=_session_lab_id())),
        disk_usage    = disk_usage,
        pipelines     = cfg.PIPELINE_MAP,
    )


@app.route("/settings/save", methods=["POST"])
@login_required
def settings_save():
    action = request.form.get("action")

    if action == "change_password":
        current = request.form.get("current_pass", "")
        new_p   = request.form.get("new_pass", "")
        confirm = request.form.get("confirm_pass", "")

        # Verify current password against DB
        user = models.verify_user(session["username"], current)
        if not user:
            flash("Current password is incorrect.", "error")
        elif new_p != confirm:
            flash("New passwords do not match.", "error")
        elif len(new_p) < 8:
            flash("Password must be at least 8 characters.", "error")
        else:
            models.update_user(session["username"], password=new_p)
            flash("Password updated successfully.", "success")

    elif action == "email_settings" and session["role"] == "admin":
        enabled   = "smtp_enabled" in request.form
        smtp_user = request.form.get("smtp_user", "").strip()
        smtp_pass = request.form.get("smtp_pass", "").strip()

        # SMTP_ENABLED / SMTP_USER are non-secret runtime settings: persist
        # them to the gitignored portal/instance/settings.json, never to
        # any tracked Python source file (see config.py's
        # _load_instance_settings / Stage 2 of the Phase 5.5 security
        # hardening). SMTP_PASS is a secret and is never written by this
        # route at all -- it stays exclusively environment/.env-sourced
        # (Stage 1); this form's password field is display-only.
        new_settings = dict(cfg._INSTANCE_SETTINGS)
        new_settings["SMTP_ENABLED"] = enabled
        new_settings["SMTP_USER"] = smtp_user or cfg.SMTP_USER

        os.makedirs(cfg.INSTANCE_DIR, exist_ok=True)
        with open(cfg.INSTANCE_SETTINGS_PATH, "w", encoding="utf-8") as f:
            json.dump(new_settings, f, indent=2, sort_keys=True)
            f.write("\n")

        cfg._INSTANCE_SETTINGS = new_settings
        cfg.SMTP_ENABLED = enabled
        cfg.SMTP_USER = new_settings["SMTP_USER"]
        cfg.SMTP_FROM = f"GenRichi Portal <{cfg.SMTP_USER}>"

        if smtp_pass:
            flash(
                "Email address and enabled/disabled state saved. The SMTP "
                "password itself is not editable here -- set SMTP_PASS in "
                "portal/.env (or your process environment) and restart the "
                "portal to change it.",
                "warning",
            )
        else:
            flash("Email settings saved.", "success")

    return redirect(url_for("settings"))


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    models.init_db()
    runner.start()
    logger.info("GenRichi Portal starting on http://0.0.0.0:5000")
    app.run(host="0.0.0.0", port=5000, debug=False, threaded=True)
