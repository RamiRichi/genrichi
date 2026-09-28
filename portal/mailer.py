"""GenRichi Portal — Email notifications"""

import smtplib
import logging
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

import config as cfg

logger = logging.getLogger("mailer")


def send_completion_email(order: dict) -> bool:
    """Send a generic notice that an order's status changed.

    Deliberately carries no order data: no patient identifier, no panel/
    analysis type, no file path, no error text, and no direct link to the
    order or its report. It names only the opaque order_id and points to the
    portal's login page -- the recipient sees the order's real content only
    if, after signing in, the same lab-membership check every order/report
    route already applies (_may_access_order / the invoice lab check) lets
    them. The email itself is never a source of authorization: the address
    typed into notify_email at order creation is not assumed to be entitled
    to see the order's data, so nothing here can leak to it, correct
    recipient or not.
    """
    if not cfg.SMTP_ENABLED:
        return False

    recipient = order.get("notify_email", "")
    if not recipient:
        return False

    oid = order.get("order_id", "")
    subject = "GenRichi Portal — order update"
    login_url = f"{cfg.PORTAL_URL}/login"

    html = f"""
    <!DOCTYPE html>
    <html>
    <body style="font-family:Arial,sans-serif;background:#f4f7f9;padding:30px;">
      <div style="max-width:600px;margin:auto;background:#fff;border-radius:8px;overflow:hidden;box-shadow:0 2px 8px rgba(0,0,0,.1)">
        <div style="background:#1a2332;padding:20px 30px;">
          <h2 style="color:#fff;margin:0">🧬 GenRichi — Bioinformatics Partner Portal</h2>
        </div>
        <div style="padding:30px;">
          <h3 style="color:#1a2332">Order update</h3>
          <p>There is an update on an order in the GenRichi portal (reference {oid}).</p>
          <p>Sign in to the portal to view it, if you are entitled to:</p>
          <p><a href="{login_url}" style="background:#1a2332;color:#fff;padding:10px 20px;text-decoration:none;border-radius:4px;">Sign in</a></p>
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

        logger.info("Notification email sent to %s for order %s", recipient, oid)
        return True

    except Exception as exc:
        logger.warning("Email failed for %s: %s", oid, exc)
        return False
