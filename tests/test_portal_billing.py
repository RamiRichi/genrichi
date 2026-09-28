"""Portal billing: lab profiles, lab-scoped access, and the /order/<id>/invoice route.

Uses the Flask test client against a throw-away SQLite file with fake users and
orders. It never touches the real portal DB, uploads or logs, never starts the
runner or a pipeline, and issues no real invoice. Skipped if Flask is missing.
"""
import os
import re
import sqlite3
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

PORTAL_DIR = Path(__file__).resolve().parents[1] / "portal"
sys.path.insert(0, str(PORTAL_DIR))
os.environ.setdefault("SECRET_KEY", "unit-test-placeholder-not-a-real-secret")
os.environ.setdefault("PORTAL_PASS", "unit-test-placeholder-admin-pass")
os.environ.setdefault("SMTP_ENABLED", "false")

try:
    import flask  # noqa: F401
    HAVE_FLASK = True
except ImportError:  # pragma: no cover
    HAVE_FLASK = False

if HAVE_FLASK:
    # app.py creates upload/log directories at import time; keep that off the disk.
    with mock.patch("os.makedirs"):
        import app as portal_app
    import models
    import config as cfg

ADMIN_PW = os.environ["PORTAL_PASS"]
USER_PW = "test-pass-123"
PATIENT = {"patient_id": "PT-SECRET-77", "patient_name": "Zed Testpatient", "dob": "1970-01-01",
           "tumor_type": "Testtumor-Xyz"}
LAB_A = dict(legal_name="Alpha Labor GmbH", address_line1="Alphastr. 1", address_line2="Haus 2",
             postal_code="10115", city="Berlin", country="Deutschland")
LAB_B = dict(legal_name="Beta Diagnostik AG", address_line1="Betaweg 9", address_line2="",
             postal_code="80331", city="München", country="Deutschland")


def counts(db):
    con = sqlite3.connect(db)
    try:
        return {t: con.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
                for t in ("invoices", "invoice_counters", "audit_log", "orders", "labs")}
    finally:
        con.close()


@unittest.skipUnless(HAVE_FLASK, "Flask is not installed in this Python")
class BillingBase(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.db = os.path.join(self._tmp.name, "test.db")
        self._db_patch = mock.patch.object(models, "DB_PATH", self.db)
        self._db_patch.start()
        self._issue_patch = mock.patch.object(cfg, "INVOICE_ISSUING_ENABLED", False)
        self._issue_patch.start()
        models.init_db()
        portal_app.app.config["TESTING"] = True
        self.lab_a = models.create_lab(**LAB_A)
        self.lab_b = models.create_lab(**LAB_B)
        models.create_user("a1", USER_PW, "lab_staff", "Anna A", "", self.lab_a)
        models.create_user("a2", USER_PW, "lab_staff", "Arne A", "", self.lab_a)
        models.create_user("b1", USER_PW, "lab_staff", "Bea B", "", self.lab_b)
        models.create_user("c1", USER_PW, "lab_staff", "Cem C", "", None)
        self.order_a = self.make_order("a1")
        self.order_b = self.make_order("b1")

    def tearDown(self):
        self._issue_patch.stop()
        self._db_patch.stop()
        self._tmp.cleanup()

    def client(self, user, pw=USER_PW):
        c = portal_app.app.test_client()
        r = c.post("/login", data={"username": user, "password": pw})
        self.assertEqual(r.status_code, 302, f"login failed for {user}")
        return c

    def admin(self):
        return self.client(cfg.PORTAL_USER, ADMIN_PW)

    def make_order(self, user, done=True, **extra):
        c = self.client(user)
        # Fastq fields are left empty: since the SFTP-token server-side resolution
        # fix (test_portal_sftp_isolation.py), a non-empty value must resolve to a
        # real file inside the caller's authorized SFTP uploads root, which these
        # billing/access tests have no need to set up. Empty stays empty and is
        # always accepted (an order may be created before a file is staged).
        data = {"panel_type": "comprehensive", **PATIENT, **extra}
        r = c.post("/order/new", data=data)
        self.assertEqual(r.status_code, 302)
        oid = r.headers["Location"].rsplit("/", 1)[-1]
        if done:
            models.update_status(oid, "Done")
        return oid

    def issue(self, oid, service_date=None):
        from datetime import date
        with mock.patch.object(cfg, "INVOICE_ISSUING_ENABLED", True):
            return self.admin().post(f"/order/{oid}/invoice/issue",
                                     data={"service_date": date.today().isoformat() if service_date is None else service_date})

    def body(self, resp):
        """Visible text only: scripts, styles and tags removed, whitespace collapsed."""
        html = resp.get_data(as_text=True)
        html = re.sub(r"(?is)<(script|style)\b.*?</\1>", " ", html)
        return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html))


class TestMigration(unittest.TestCase):
    @unittest.skipUnless(HAVE_FLASK, "Flask is not installed in this Python")
    def test_existing_db_is_migrated_without_data_loss_and_is_idempotent(self):
        with tempfile.TemporaryDirectory() as d:
            db = os.path.join(d, "old.db")
            con = sqlite3.connect(db)
            con.executescript("""
                CREATE TABLE users (id INTEGER PRIMARY KEY AUTOINCREMENT, username TEXT UNIQUE NOT NULL,
                    password_hash TEXT NOT NULL, role TEXT DEFAULT 'lab_staff', created_at TEXT NOT NULL);
                CREATE TABLE orders (id INTEGER PRIMARY KEY AUTOINCREMENT, order_id TEXT UNIQUE NOT NULL,
                    patient_id TEXT NOT NULL, panel_type TEXT NOT NULL, status TEXT DEFAULT 'Queued',
                    created_at TEXT NOT NULL);
                INSERT INTO users (username, password_hash, role, created_at) VALUES ('old', 'x', 'lab_staff', '2026-01-01');
                INSERT INTO orders (order_id, patient_id, panel_type, created_at) VALUES ('GR-OLD-1', 'P1', 'comprehensive', '2026-01-01');
            """)
            con.commit()
            con.close()
            with mock.patch.object(models, "DB_PATH", db):
                models.init_db()
                models.init_db()          # second run must be a no-op
            con = sqlite3.connect(db)
            self.assertEqual(con.execute("SELECT username FROM users WHERE username='old'").fetchone()[0], "old")
            self.assertEqual(con.execute("SELECT patient_id, lab_id FROM orders").fetchone(), ("P1", None))
            self.assertIn("lab_id", {r[1] for r in con.execute("PRAGMA table_info(users)")})
            self.assertIn("lab_id", {r[1] for r in con.execute("PRAGMA table_info(orders)")})
            tables = {r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type='table'")}
            self.assertTrue({"labs", "invoices", "invoice_counters"} <= tables)
            con.close()


class TestInvoiceAccess(BillingBase):
    def test_get_never_writes_and_draft_is_admin_only(self):
        before = counts(self.db)
        r = self.admin().get(f"/order/{self.order_a}/invoice")
        self.assertEqual(r.status_code, 200)
        text = self.body(r)
        self.assertIn("Rechnungsentwurf", text)
        self.assertIn("wird bei Ausstellung vergeben", text)
        self.assertNotRegex(text, r"GR-\d{4}-\d{4}")
        self.assertEqual(counts(self.db), before, "GET must not write anything")
        self.assertEqual(self.client("a1").get(f"/order/{self.order_a}/invoice").status_code, 404)
        self.assertEqual(counts(self.db), before)

    def test_issuing_is_off_by_default_and_admin_only(self):
        oid = self.order_a
        r = self.admin().post(f"/order/{oid}/invoice/issue", data={"service_date": "2026-01-01"})
        self.assertEqual(r.status_code, 302)
        self.assertEqual(counts(self.db)["invoices"], 0)
        with mock.patch.object(cfg, "INVOICE_ISSUING_ENABLED", True):
            r = self.client("a1").post(f"/order/{oid}/invoice/issue", data={"service_date": "2026-01-01"})
        self.assertEqual(r.status_code, 302)
        self.assertEqual(counts(self.db)["invoices"], 0)

    def test_issued_invoice_content_recipient_and_no_patient_data(self):
        self.issue(self.order_a)
        self.assertEqual(counts(self.db)["invoices"], 1)
        r = self.client("a1").get(f"/order/{self.order_a}/invoice")
        self.assertEqual(r.status_code, 200)
        text = self.body(r)
        for needle in ("Alpha Labor GmbH", "Alphastr. 1", "Haus 2", "10115 Berlin", self.order_a,
                       "Referenz (Auftragsnr.)", "Leistungsdatum", "§ 19 UStG", "Steuernummer: 114/262/00939",
                       "GenRichi Diagnostics", "[IBAN nicht konfiguriert]"):
            self.assertIn(needle, text)
        for secret in (PATIENT["patient_id"], PATIENT["patient_name"], PATIENT["tumor_type"], "1970", "Anna"):
            self.assertNotIn(secret, text)
        for gone in ("MwSt", "USt-IdNr", "Beta Diagnostik"):
            self.assertNotIn(gone, text)
        self.assertRegex(text, r"Gesamtbetrag: 750\.00 €")
        self.assertNotIn("Nettobetrag", text)
        draft = self.body(self.admin().get(f"/order/{self.order_a}/invoice"))
        for secret in (PATIENT["patient_id"], PATIENT["patient_name"], PATIENT["tumor_type"]):
            self.assertNotIn(secret, draft)

    def test_lab_members_share_and_other_labs_are_blocked(self):
        self.issue(self.order_a)
        self.issue(self.order_b)
        for user in ("a1", "a2"):                       # same lab, incl. non-creator
            self.assertEqual(self.client(user).get(f"/order/{self.order_a}/invoice").status_code, 200, user)
            self.assertEqual(self.client(user).get(f"/order/{self.order_b}/invoice").status_code, 403, user)
        self.assertEqual(self.client("b1").get(f"/order/{self.order_a}/invoice").status_code, 403)
        self.assertEqual(self.client("b1").get(f"/order/{self.order_b}/invoice").status_code, 200)
        self.assertEqual(self.client("c1").get(f"/order/{self.order_a}/invoice").status_code, 403)
        anon = portal_app.app.test_client().get(f"/order/{self.order_a}/invoice")
        self.assertEqual(anon.status_code, 302)
        self.assertIn("/login", anon.headers["Location"])
        self.assertEqual(self.admin().get("/order/GR-NOPE/invoice").status_code, 404)

    def test_lab_is_taken_from_the_user_not_the_form(self):
        oid = self.make_order("a1", lab_id=str(self.lab_b))
        self.assertEqual(models.get_order(oid)["lab_id"], self.lab_a)
        self.assertEqual(models.get_order(self.order_b)["lab_id"], self.lab_b)
        oid_c = self.make_order("c1")
        self.assertIsNone(models.get_order(oid_c)["lab_id"])

    def test_moving_a_user_to_another_lab_moves_access_immediately(self):
        self.issue(self.order_a)
        a2 = self.client("a2")
        self.assertEqual(a2.get(f"/order/{self.order_a}/invoice").status_code, 200)
        self.admin().post("/admin/users/a2/lab", data={"lab_id": str(self.lab_b)})
        self.assertEqual(a2.get(f"/order/{self.order_a}/invoice").status_code, 403)

    def test_order_without_lab_is_unbillable_for_staff_and_flagged_for_admin(self):
        oid = self.make_order("c1")
        self.assertEqual(self.client("c1").get(f"/order/{oid}/invoice").status_code, 403)
        text = self.body(self.admin().get(f"/order/{oid}/invoice"))
        self.assertIn("kein Labor zugeordnet", text)
        self.assertEqual(self.issue(oid).status_code, 302)
        self.assertEqual(counts(self.db)["invoices"], 0)
        self.admin().post(f"/admin/orders/{oid}/lab", data={"lab_id": str(self.lab_a)})
        self.assertEqual(models.get_order(oid)["lab_id"], self.lab_a)


class TestInvoiceNumbering(BillingBase):
    def test_number_is_fixed_unique_gapless_and_snapshot_is_immutable(self):
        from datetime import date
        self.issue(self.order_a)
        n1 = models.get_invoice(self.order_a)["invoice_no"]
        self.assertEqual(n1, f"GR-{date.today().year}-0001")
        self.issue(self.order_a)                        # re-issue: idempotent
        self.assertEqual(models.get_invoice(self.order_a)["invoice_no"], n1)
        self.assertEqual(counts(self.db)["invoices"], 1)
        self.issue(self.order_b)
        self.assertEqual(models.get_invoice(self.order_b)["invoice_no"], f"GR-{date.today().year}-0002")
        shown = [re.search(r"GR-\d{4}-\d{4}", self.body(self.client("a1").get(f"/order/{self.order_a}/invoice"))).group(0)
                 for _ in range(2)]
        self.assertEqual(shown, [n1, n1])
        # editing the lab profile must not change an issued invoice
        r = self.admin().post(f"/admin/labs/{self.lab_a}/edit", data={**LAB_A, "legal_name": "Alpha Renamed GmbH"})
        self.assertEqual(r.status_code, 302)
        text = self.body(self.client("a1").get(f"/order/{self.order_a}/invoice"))
        self.assertIn("Alpha Labor GmbH", text)
        self.assertNotIn("Alpha Renamed", text)

    def test_service_date_and_status_are_validated(self):
        from datetime import date, timedelta
        bad = ("", "not-a-date", (date.today() + timedelta(days=1)).isoformat(), "2000-01-01")
        for value in bad:
            self.assertEqual(self.issue(self.order_a, service_date=value).status_code, 302)
        self.assertEqual(counts(self.db)["invoices"], 0)
        running = self.make_order("a1", done=False)
        self.issue(running)
        self.assertEqual(counts(self.db)["invoices"], 0)


class TestLabAdminPermissions(BillingBase):
    def test_lab_staff_cannot_manage_labs_or_relink_themselves(self):
        a1 = self.client("a1")
        before = counts(self.db)
        self.assertEqual(a1.get("/admin/labs").status_code, 302)
        r = a1.post("/admin/labs/create", data={**LAB_A, "legal_name": "Evil Lab"})
        self.assertEqual(r.status_code, 302)
        r = a1.post(f"/admin/labs/{self.lab_b}/edit", data={**LAB_B, "legal_name": "Hijacked"})
        self.assertEqual(r.status_code, 302)
        a1.post("/admin/users/a1/lab", data={"lab_id": str(self.lab_b)})
        a1.post(f"/admin/orders/{self.order_a}/lab", data={"lab_id": str(self.lab_b)})
        self.assertEqual(counts(self.db), before)
        self.assertEqual(models.get_user("a1")["lab_id"], self.lab_a)
        self.assertEqual(models.get_lab(self.lab_b)["legal_name"], LAB_B["legal_name"])
        self.assertEqual(models.get_order(self.order_a)["lab_id"], self.lab_a)

    def test_my_lab_shows_only_own_lab(self):
        text = self.body(self.client("a1").get("/lab"))
        self.assertIn("Alpha Labor GmbH", text)
        self.assertIn("a2", text)
        self.assertNotIn("Beta Diagnostik", text)
        self.assertNotIn("Bea B", text)
        self.assertNotRegex(text, r"b1")
        self.assertIn("not linked to a lab", self.body(self.client("c1").get("/lab")))
        self.assertEqual(self.client("a1").get("/lab?id=2").status_code, 200)   # parameter is ignored
        self.assertNotIn("Beta Diagnostik", self.body(self.client("a1").get("/lab?id=2")))

    def test_admin_lab_crud_and_validation(self):
        adm = self.admin()
        r = adm.get("/admin/labs")
        self.assertEqual(r.status_code, 200)
        text = self.body(r)
        self.assertIn("Alpha Labor GmbH", text)
        self.assertIn("Beta Diagnostik AG", text)
        n = counts(self.db)["labs"]
        adm.post("/admin/labs/create", data={**LAB_A})                              # duplicate name
        adm.post("/admin/labs/create", data={**LAB_A, "legal_name": "", })          # missing name
        adm.post("/admin/labs/create", data={**LAB_A, "legal_name": "X", "city": ""})
        self.assertEqual(counts(self.db)["labs"], n)
        adm.post("/admin/labs/create", data={**LAB_A, "legal_name": "Gamma  Labor\nGmbH"})
        self.assertEqual(counts(self.db)["labs"], n + 1)
        self.assertIn("Gamma Labor GmbH", [l["legal_name"] for l in models.list_labs()])
        self.assertEqual(adm.post("/admin/labs/999/edit", data=LAB_A).status_code, 404)

    def test_admin_can_link_user_and_html_is_escaped(self):
        adm = self.admin()
        adm.post("/admin/labs/create", data={**LAB_A, "legal_name": "<script>alert(1)</script> Lab"})
        html = adm.get("/admin/labs").get_data(as_text=True)
        self.assertNotIn("<script>alert(1)</script>", html)
        lab = [l for l in models.list_labs() if "script" in l["legal_name"]][0]
        adm.post("/admin/users/c1/lab", data={"lab_id": str(lab["id"])})
        self.assertEqual(models.get_user("c1")["lab_id"], lab["id"])
        adm.post("/admin/users/c1/lab", data={"lab_id": ""})
        self.assertIsNone(models.get_user("c1")["lab_id"])
        adm.post("/admin/users/c1/lab", data={"lab_id": "9999"})
        self.assertIsNone(models.get_user("c1")["lab_id"])


if __name__ == "__main__":
    unittest.main()
