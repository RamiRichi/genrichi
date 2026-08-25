"""
GenRichi Portal — Flask Web Application
Multi-user edition (admin + lab_staff roles)
"""

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
    orders = models.list_orders(200, username=uname, role=role)
    stats = {
        "total":   len(orders),
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
    orders = models.list_orders(500, username=uname, role=role)

    stats = {
        "total":   len(orders),
        "done":    sum(1 for o in orders if o["status"] == "Done"),
        "failed":  sum(1 for o in orders if o["status"] == "Failed"),
        "running": sum(1 for o in orders if o["status"] == "Running"),
        "queued":  sum(1 for o in orders if o["status"] == "Queued"),
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
                               paired_panels=list(cfg.PAIRED_PANELS))

    panel_type   = request.form.get("panel_type", "hotspot")
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

    order_id = models.new_order(
        patient_id=patient_id, patient_name=patient_name,
        sex=sex, dob=dob, tumor_type=tumor_type,
        panel_type=panel_type,
        fastq_r1=r1, fastq_r2=r2,
        fastq_normal_r1=normal_r1, fastq_normal_r2=normal_r2,
        notes=notes, notify_email=notify_email,
        created_by=session["username"],
    )

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
    if session["role"] != "admin" and order["created_by"] != session["username"]:
        abort(403)

    log_tail = ""
    if order["log_path"] and os.path.isfile(order["log_path"]):
        with open(order["log_path"]) as f:
            lines    = f.readlines()
            log_tail = "".join(lines[-60:])

    return render_template("order.html", order=order, log_tail=log_tail,
                           pipeline=cfg.PIPELINE_MAP.get(order["panel_type"], {}))


# ── Report viewer ─────────────────────────────────────────────────────────────
@app.route("/order/<order_id>/report/embed")
@login_required
def report_embed(order_id):
    order = models.get_order(order_id)
    if not order or not order["report_path"]:
        abort(404)
    if session["role"] != "admin" and order["created_by"] != session["username"]:
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
    if session["role"] != "admin" and order["created_by"] != session["username"]:
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
    if session["role"] != "admin" and order["created_by"] != session["username"]:
        abort(403)
    if not os.path.isfile(order["report_path"]):
        abort(404)
    fname = f"{order_id}_{order['panel_type']}_report.html"
    return send_file(order["report_path"], mimetype="text/html",
                     as_attachment=True, download_name=fname)


# ── Invoice ───────────────────────────────────────────────────────────────────
@app.route("/order/<order_id>/invoice")
@login_required
def invoice(order_id):
    from datetime import date, timedelta
    order = models.get_order(order_id)
    if not order:
        abort(404)
    if session["role"] != "admin" and order["created_by"] != session["username"]:
        abort(403)

    net_price   = cfg.PANEL_PRICES.get(order["panel_type"], 0.0)
    vat         = round(net_price * 0.19, 2)
    gross_price = round(net_price + vat, 2)
    today       = date.today()
    invoice_no  = f"GR-{today.strftime('%Y%m')}-{order_id[-6:]}"
    pipeline_label = cfg.PIPELINE_MAP.get(order["panel_type"], {}).get("label", order["panel_type"])

    return render_template("invoice.html",
        order=order,
        cfg=cfg,
        invoice_no=invoice_no,
        invoice_date=today.strftime("%d.%m.%Y"),
        due_date=(today + timedelta(days=cfg.PAYMENT_DAYS)).strftime("%d.%m.%Y"),
        net_price=net_price,
        vat=vat,
        gross_price=gross_price,
        pipeline_label=pipeline_label,
    )


# ── SFTP file browser API ────────────────────────────────────────────────────
@app.route("/api/sftp-files")
@login_required
def api_sftp_files():
    sftp_root = "/srv/genrichi-sftp"
    files = []
    try:
        for user_dir in sorted(Path(sftp_root).iterdir()):
            uploads = user_dir / "uploads"
            if not uploads.is_dir():
                continue
            for f in sorted(uploads.iterdir()):
                if f.suffix.lower() in (".gz", ".fastq", ".fq") or \
                   f.name.endswith(".fastq.gz") or f.name.endswith(".fq.gz"):
                    files.append({
                        "client": user_dir.name,
                        "name":   f.name,
                        "path":   str(f),
                        "size_mb": round(f.stat().st_size / 1024 / 1024, 1),
                    })
    except Exception:
        pass
    return jsonify(files)


# ── Status API ────────────────────────────────────────────────────────────────
@app.route("/api/order/<order_id>/status")
@login_required
def api_status(order_id):
    order = models.get_order(order_id)
    if not order:
        abort(404)
    if session["role"] != "admin" and order["created_by"] != session["username"]:
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
    if session["role"] != "admin" and order["created_by"] != session["username"]:
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
    if session["role"] != "admin" and order["created_by"] != session["username"]:
        abort(403)
    if order["status"] != "Failed":
        flash("Only failed orders can be retried.", "error")
        return redirect(url_for("order_detail", order_id=order_id))
    models.update_status(order_id, "Queued", error_msg="")
    flash(f"Order {order_id} re-queued.", "success")
    return redirect(url_for("order_detail", order_id=order_id))


# ── Cancel ────────────────────────────────────────────────────────────────────
@app.route("/order/<order_id>/cancel", methods=["POST"])
@login_required
def cancel_order(order_id):
    import signal
    order = models.get_order(order_id)
    if not order:
        abort(404)
    if session["role"] != "admin" and order["created_by"] != session["username"]:
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
    return render_template("users.html", users=users)


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

    ok = models.create_user(username, password, role, full_name, email)
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

    return render_template("settings.html",
        smtp_enabled  = cfg.SMTP_ENABLED,
        smtp_user     = cfg.SMTP_USER,
        smtp_pass     = "••••••••" if cfg.SMTP_PASS else "",
        portal_url    = cfg.PORTAL_URL,
        admin_user    = cfg.PORTAL_USER,
        db_path       = cfg.DB_PATH,
        workflow_dir  = cfg.WORKFLOW_DIR,
        total_orders  = len(models.list_orders(1000)),
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

        config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "config.py")
        with open(config_path, "r") as f:
            content = f.read()

        content = content.replace(
            f'SMTP_ENABLED  = {cfg.SMTP_ENABLED}',
            f'SMTP_ENABLED  = {enabled}'
        )
        if smtp_user:
            content = content.replace(
                f'SMTP_USER     = "{cfg.SMTP_USER}"',
                f'SMTP_USER     = "{smtp_user}"'
            )
            content = content.replace(
                f'SMTP_FROM     = "GenRichi Portal <{cfg.SMTP_USER}>"',
                f'SMTP_FROM     = "GenRichi Portal <{smtp_user}>"'
            )
        if smtp_pass and smtp_pass != "••••••••":
            content = content.replace(
                f'SMTP_PASS     = "{cfg.SMTP_PASS}"',
                f'SMTP_PASS     = "{smtp_pass}"'
            )

        with open(config_path, "w") as f:
            f.write(content)

        cfg.SMTP_ENABLED = enabled
        if smtp_user: cfg.SMTP_USER = smtp_user
        if smtp_pass and smtp_pass != "••••••••": cfg.SMTP_PASS = smtp_pass

        flash("Email settings saved.", "success")

    return redirect(url_for("settings"))


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    models.init_db()
    runner.start()
    logger.info("GenRichi Portal starting on http://0.0.0.0:5000")
    app.run(host="0.0.0.0", port=5000, debug=False, threaded=True)
