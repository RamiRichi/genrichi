"""GenRichi Portal — Email notifications"""

import smtplib
import logging
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

import config as cfg

logger = logging.getLogger("mailer")


def send_completion_email(order: dict, report_url: str = "") -> bool:
    """Send an email when an order completes (Done or Failed)."""
    if not cfg.SMTP_ENABLED:
        return False

    recipient = order.get("notify_email", "")
    if not recipient:
        return False

    status  = order.get("status", "Unknown")
    oid     = order.get("order_id", "")
    pid     = order.get("patient_id", "")
    panel   = order.get("panel_type", "")
    panel_label = cfg.PIPELINE_MAP.get(panel, {}).get("label", panel)

    subject = f"GenRichi — Order {oid} {status}"

    if status == "Done":
        color  = "#28a745"
        icon   = "✅"
        body   = f"""
        <p>The analysis for <strong>{pid}</strong> has completed successfully.</p>
        <p><strong>Panel:</strong> {panel_label}</p>
        {"<p><a href='" + report_url + "' style='background:#28a745;color:#fff;padding:10px 20px;text-decoration:none;border-radius:4px;'>View Report</a></p>" if report_url else ""}
        """
    else:
        color  = "#dc3545"
        icon   = "❌"
        error  = order.get("error_msg", "Unknown error")
        body   = f"""
        <p>The analysis for <strong>{pid}</strong> has <strong>failed</strong>.</p>
        <p><strong>Panel:</strong> {panel_label}</p>
        <p><strong>Error:</strong> {error}</p>
        <p><a href='{cfg.PORTAL_URL}/order/{oid}'>View details and retry</a></p>
        """

    html = f"""
    <!DOCTYPE html>
    <html>
    <body style="font-family:Arial,sans-serif;background:#f4f7f9;padding:30px;">
      <div style="max-width:600px;margin:auto;background:#fff;border-radius:8px;overflow:hidden;box-shadow:0 2px 8px rgba(0,0,0,.1)">
        <div style="background:#1a2332;padding:20px 30px;">
          <h2 style="color:#fff;margin:0">🧬 GenRichi — Bioinformatics Partner Portal</h2>
        </div>
        <div style="padding:30px;">
          <h3 style="color:{color}">{icon} Order {oid} — {status}</h3>
          {body}
          <hr style="margin:20px 0;border:none;border-top:1px solid #eee">
          <p style="color:#888;font-size:12px">GenRichi Bioinformatics Partner Portal &bull; {cfg.PORTAL_URL}</p>
        </div>
      </div>
    </body>
    </html>
    """

    try:
        msg = MIMEMultipart("alternative")
        msg["Subject"] = subject
        msg["From"]    = cfg.SMTP_FROM
        msg["To"]      = recipient
        msg.attach(MIMEText(html, "html"))

        with smtplib.SMTP(cfg.SMTP_HOST, cfg.SMTP_PORT, timeout=15) as server:
            server.starttls()
            server.login(cfg.SMTP_USER, cfg.SMTP_PASS)
            server.sendmail(cfg.SMTP_USER, [recipient], msg.as_string())

        logger.info("Email sent to %s for order %s", recipient, oid)
        return True

    except Exception as exc:
        logger.warning("Email failed for %s: %s", oid, exc)
        return False
