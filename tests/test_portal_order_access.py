"""Portal order access: the unified lab-sharing rule for order/report/status/log/action routes.

Reuses the BillingBase fixture from test_portal_billing.py (same throw-away SQLite DB,
same fake labs/users; never touches the real portal DB, uploads, logs, runner or a
real pipeline). Covers order_detail, report views, download, the status API, the log
download, retry and cancel, plus dashboard/stats list-leak checks -- the routes the
billing route (/order/<id>/invoice, tested separately) was deliberately NOT extended
to change.
"""
import os
import re
import sys
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_portal_billing import BillingBase, HAVE_FLASK, LAB_A, LAB_B, USER_PW  # noqa: E402

if HAVE_FLASK:
    import models


@unittest.skipUnless(HAVE_FLASK, "Flask is not installed in this Python")
class TestOrderAccess(BillingBase):
    def body(self, resp):
        html = resp.get_data(as_text=True)
        html = re.sub(r"(?is)<(script|style)\b.*?</\1>", " ", html)
        return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html))

    def assert_ok(self, user, path, msg=""):
        self.assertEqual(self.client(user).get(path).status_code, 200, f"{user} {path} {msg}")

    def assert_blocked(self, user, path, msg=""):
        self.assertEqual(self.client(user).get(path).status_code, 403, f"{user} {path} {msg}")

    def test_lab_colleague_shares_order_detail_and_status(self):
        # a2 did not create order_a (a1 did) but is in the same lab -- this is the bug being fixed.
        self.assert_ok("a1", f"/order/{self.order_a}")
        self.assert_ok("a2", f"/order/{self.order_a}")
        self.assert_ok("a2", f"/api/order/{self.order_a}/status")
        self.assertEqual(self.client("a2").get(f"/api/order/{self.order_a}/status").get_json()["status"], "Done")

    def test_other_lab_and_no_lab_user_are_blocked_from_order_detail(self):
        self.assert_blocked("b1", f"/order/{self.order_a}")
        self.assert_blocked("c1", f"/order/{self.order_a}")
        self.assert_blocked("a1", f"/order/{self.order_b}")

    def test_admin_sees_every_order(self):
        adm = self.admin()
        self.assertEqual(adm.get(f"/order/{self.order_a}").status_code, 200)
        self.assertEqual(adm.get(f"/order/{self.order_b}").status_code, 200)

    def test_report_status_log_retry_cancel_all_share_the_same_rule(self):
        # Give order_a a log file and mark it Running so retry/cancel have something to act on,
        # and drop a fake report file so the report/download routes have output to serve.
        order = models.get_order(self.order_a)
        log_dir = Path(self._tmp.name)
        log_path = log_dir / "order_a.log"
        log_path.write_text("line1\nline2\n")
        report_path = log_dir / "order_a_report.html"
        report_path.write_text("<html>report</html>")
        models.update_status(self.order_a, "Running", log_path=str(log_path), report_path=str(report_path), pid=None)

        for user, code in (("a2", 200), ("b1", 403), ("c1", 403)):
            for path in (f"/order/{self.order_a}/report", f"/order/{self.order_a}/report/embed",
                        f"/order/{self.order_a}/report/download", f"/order/{self.order_a}/log"):
                self.assertEqual(self.client(user).get(path).status_code, code, f"{user} {path}")

        # retry/cancel are POST-only actions guarded by the same rule (via order status
        # preconditions they may redirect with a flash rather than 403 on success, but a
        # blocked user must never reach the point of changing status).
        b1 = self.client("b1")
        r = b1.post(f"/order/{self.order_a}/cancel")
        self.assertEqual(r.status_code, 403)
        self.assertEqual(models.get_order(self.order_a)["status"], "Running")

        a2 = self.client("a2")
        r = a2.post(f"/order/{self.order_a}/cancel")
        self.assertEqual(r.status_code, 302)
        self.assertEqual(models.get_order(self.order_a)["status"], "Failed")

        c1 = self.client("c1")
        r = c1.post(f"/order/{self.order_a}/retry")
        self.assertEqual(r.status_code, 403)
        self.assertEqual(models.get_order(self.order_a)["status"], "Failed")

        r = a2.post(f"/order/{self.order_a}/retry")
        self.assertEqual(r.status_code, 302)
        self.assertEqual(models.get_order(self.order_a)["status"], "Queued")

    def test_legacy_no_lab_order_is_creator_only_not_lab_wide(self):
        """An order predating lab-linking (lab_id IS NULL) must not be widened to the
        creator's whole lab -- it stays exactly as narrow as it always was."""
        c1_order = self.make_order("c1")  # c1 has no lab -> lab_id is None
        self.assertIsNone(models.get_order(c1_order)["lab_id"])
        self.assert_ok("c1", f"/order/{c1_order}")
        # a1/a2 are in a real lab but that is irrelevant here: c1_order has no lab at all.
        self.assert_blocked("a1", f"/order/{c1_order}")
        self.assert_blocked("b1", f"/order/{c1_order}")

        # Same check for a *lab* user's legacy order: give a1 a no-lab order directly in
        # the DB (simulating data that predates lab_id existing) and confirm a2 (same lab)
        # still cannot see it -- legacy access is never retroactively widened to the lab.
        legacy = self.make_order("a1")
        with models._conn() as conn:
            conn.execute("UPDATE orders SET lab_id=NULL WHERE order_id=?", (legacy,))
        self.assertIsNone(models.get_order(legacy)["lab_id"])
        self.assert_ok("a1", f"/order/{legacy}")
        self.assert_blocked("a2", f"/order/{legacy}")

    def test_dashboard_and_stats_never_list_another_labs_orders(self):
        for path in ("/", "/stats"):
            text = self.body(self.client("a2").get(path))
            self.assertIn(self.order_a, text)
            self.assertNotIn(self.order_b, text)
            text_b = self.body(self.client("b1").get(path))
            self.assertIn(self.order_b, text_b)
            self.assertNotIn(self.order_a, text_b)

    def test_settings_total_orders_count_does_not_leak_other_labs(self):
        # Before this fix, /settings called models.list_orders(1000) with no
        # username/role, which defaulted to the admin-wide total for every
        # logged-in user regardless of lab. Both pages must still render.
        self.assertEqual(self.client("a2").get("/settings").status_code, 200)
        self.assertEqual(self.client("b1").get("/settings").status_code, 200)
        # Assert via the exact model call the route now makes, per session identity.
        import config as cfg
        self.assertEqual(len(models.list_orders(1000, username="a2", role="lab_staff", lab_id=self.lab_a)), 1)
        self.assertEqual(len(models.list_orders(1000, username="b1", role="lab_staff", lab_id=self.lab_b)), 1)
        self.assertEqual(len(models.list_orders(1000, username=cfg.PORTAL_USER, role="admin", lab_id=None)), 2)


if __name__ == "__main__":
    unittest.main()
