"""portal/mailer.py: completion/failure notifications must be fully generic.

Verifies send_completion_email() never puts patient data, analysis/panel type,
file paths, error text, or a direct order/report link into the email it sends,
and instead points only to the portal login page. Mocks smtplib.SMTP entirely
-- no real network connection, no real SMTP credentials, no real send.
"""
import email
import os
import sys
import unittest
from pathlib import Path
from unittest import mock

PORTAL_DIR = Path(__file__).resolve().parents[1] / "portal"
sys.path.insert(0, str(PORTAL_DIR))
os.environ.setdefault("SECRET_KEY", "unit-test-placeholder-not-a-real-secret")
os.environ.setdefault("PORTAL_PASS", "unit-test-placeholder-admin-pass")
os.environ.setdefault("SMTP_ENABLED", "true")
os.environ.setdefault("SMTP_PASS", "unit-test-placeholder-smtp-pass")

import config as cfg  # noqa: E402
import mailer  # noqa: E402

SENSITIVE = {
    "patient_id": "PT-SECRET-99",
    "patient_name": "Jane Q. Testpatient",
    "panel_type": "comprehensive",
    "fastq_r1": "/srv/genrichi-sftp/acme_lab/uploads/tumor_R1.fastq.gz",
    "report_path": "/home/rami/genrichi/results/GR-20260927-ABCDEF/report/x.html",
    "error_msg": "Snakemake exit code 1: MuTect2 failed on chr17 -- see log",
    "notify_email": "labstaff@yourlab.de",
}


def make_order(status):
    return {"order_id": "GR-20260927-ABCDEF", "status": status, **SENSITIVE}


class TestMailerPrivacy(unittest.TestCase):
    def send_and_capture(self, order):
        sent = {}

        class FakeSMTP:
            def __init__(self, *a, **kw):
                pass

            def __enter__(self):
                return self

            def __exit__(self, *a):
                return False

            def starttls(self):
                pass

            def login(self, *a):
                pass

            def sendmail(self, from_addr, to_addrs, msg_str):
                sent["from"] = from_addr
                sent["to"] = to_addrs
                sent["raw"] = msg_str

        with mock.patch("smtplib.SMTP", FakeSMTP), \
             mock.patch.object(cfg, "SMTP_ENABLED", True), \
             mock.patch.object(cfg, "PORTAL_URL", "https://portal.genrichi.de"):
            ok = mailer.send_completion_email(order)
        self.assertTrue(ok)
        # Decode the actual MIME payload (base64 HTML body) rather than matching
        # the raw wire format -- otherwise every assertion below is vacuous.
        parsed = email.message_from_string(sent["raw"])
        html_part = next(p for p in parsed.walk() if p.get_content_type() == "text/html")
        body = html_part.get_payload(decode=True).decode("utf-8")
        subject = str(email.header.make_header(email.header.decode_header(parsed["Subject"])))
        return subject + "\n" + body

    def assert_generic_and_private(self, raw):
        self.assertIn("GR-20260927-ABCDEF", raw)  # the opaque order_id is fine
        self.assertIn("https://portal.genrichi.de/login", raw)  # points at login, not the order
        for secret in (SENSITIVE["patient_id"], SENSITIVE["patient_name"], SENSITIVE["panel_type"],
                      SENSITIVE["fastq_r1"], SENSITIVE["report_path"], SENSITIVE["error_msg"],
                      "comprehensive", "acme_lab", "MuTect2", "chr17"):
            self.assertNotIn(secret, raw, f"leaked: {secret!r}")
        self.assertNotIn("/order/GR-20260927-ABCDEF", raw)  # no direct order/report link at all
        self.assertNotIn("/report", raw)

    def test_done_order_email_is_generic(self):
        self.assert_generic_and_private(self.send_and_capture(make_order("Done")))

    def test_failed_order_email_is_generic_and_has_no_error_text(self):
        self.assert_generic_and_private(self.send_and_capture(make_order("Failed")))

    def test_no_recipient_sends_nothing(self):
        order = make_order("Done")
        order["notify_email"] = ""
        with mock.patch("smtplib.SMTP") as smtp:
            self.assertFalse(mailer.send_completion_email(order))
            smtp.assert_not_called()

    def test_smtp_disabled_sends_nothing(self):
        with mock.patch.object(cfg, "SMTP_ENABLED", False), mock.patch("smtplib.SMTP") as smtp:
            self.assertFalse(mailer.send_completion_email(make_order("Done")))
            smtp.assert_not_called()

    def test_signature_no_longer_accepts_a_report_url(self):
        # A stray report_url argument (the old call shape) must not silently
        # be accepted and re-exposed; the function takes only the order.
        with self.assertRaises(TypeError):
            mailer.send_completion_email(make_order("Done"), "https://portal.genrichi.de/order/x/report")


if __name__ == "__main__":
    unittest.main()
