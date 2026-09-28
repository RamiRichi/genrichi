"""Portal SFTP upload isolation: /api/sftp-files and the server-side FASTQ-field
resolution in /order/new.

Reuses the BillingBase fixture (throw-away SQLite DB, fake labs/users; never
touches the real portal DB, uploads, logs, runner or a real pipeline). Builds a
fake SFTP_ROOT under a temp directory instead of /srv/genrichi-sftp.
"""
import os
import sys
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_portal_billing import BillingBase, HAVE_FLASK, LAB_A, LAB_B  # noqa: E402

if HAVE_FLASK:
    import app as portal_app
    import models
    import sftp_paths


@unittest.skipUnless(HAVE_FLASK, "Flask is not installed in this Python")
class TestSftpIsolation(BillingBase):
    def setUp(self):
        super().setUp()
        self._sftp_tmpdir = _TempSftpRoot()
        self.sftp_root = self._sftp_tmpdir.__enter__()
        # SFTP_ROOT and the path-safety functions live in sftp_paths.py (app.py
        # and runner.py both just import them from there), so that is the
        # module to patch -- patching portal_app.SFTP_ROOT would no longer have
        # any effect on sftp_paths.sftp_upload_root()'s own module global.
        self._root_patch = mock.patch.object(sftp_paths, "SFTP_ROOT", self.sftp_root)
        self._root_patch.start()

        def mk(client, *names):
            d = Path(self.sftp_root) / client / "uploads"
            d.mkdir(parents=True, exist_ok=True)
            for n in names:
                (d / n).write_bytes(b"@read\nACGT\n+\nIIII\n")
            return d

        mk("acme_lab", "tumor_R1.fastq.gz", "tumor_R2.fastq.gz")
        mk("other_lab", "secret_patient_R1.fastq.gz")
        models.set_lab_sftp_account(self.lab_a, "acme_lab")
        # lab_b deliberately left unlinked (no sftp_account) for the "not linked" tests.

    def tearDown(self):
        self._root_patch.stop()
        self._sftp_tmpdir.__exit__(None, None, None)
        super().tearDown()

    # ── /api/sftp-files ──────────────────────────────────────────────────────
    def test_lab_user_sees_only_their_own_lab_files_no_absolute_path(self):
        data = self.client("a1").get("/api/sftp-files").get_json()
        self.assertTrue(data["linked"])
        names = {f["name"] for f in data["files"]}
        self.assertEqual(names, {"tumor_R1.fastq.gz", "tumor_R2.fastq.gz"})
        for f in data["files"]:
            self.assertEqual(f["client"], "acme_lab")
            self.assertEqual(f["token"], f"acme_lab/{f['name']}")
            self.assertNotIn("path", f)
            self.assertNotIn(self.sftp_root, str(f))  # never the real filesystem path

    def test_unlinked_lab_gets_explicit_not_linked_not_a_silent_empty_list(self):
        data = self.client("b1").get("/api/sftp-files").get_json()
        self.assertFalse(data["linked"])
        self.assertEqual(data["files"], [])

    def test_no_lab_user_gets_not_linked_too(self):
        data = self.client("c1").get("/api/sftp-files").get_json()
        self.assertFalse(data["linked"])
        self.assertEqual(data["files"], [])

    def test_admin_sees_every_clients_files(self):
        data = self.admin().get("/api/sftp-files").get_json()
        self.assertTrue(data["linked"])
        clients = {f["client"] for f in data["files"]}
        self.assertEqual(clients, {"acme_lab", "other_lab"})

    # ── /order/new server-side resolution ───────────────────────────────────
    def order_form(self, **extra):
        return {"panel_type": "comprehensive", "patient_id": "PT-1", **extra}

    def test_lab_user_can_create_order_with_own_labs_token(self):
        # The order is created immediately with the (client, filename) recorded
        # for pinning -- it does NOT store a live SFTP path (that never happens
        # at all any more; see test_portal_input_pinning.py for the full
        # select -> pin -> Queued -> pinned-path flow this feeds into).
        r = self.client("a1").post("/order/new", data=self.order_form(fastq_r1="acme_lab/tumor_R1.fastq.gz"))
        self.assertEqual(r.status_code, 302)
        oid = r.headers["Location"].rsplit("/", 1)[-1]
        order = models.get_order(oid)
        self.assertEqual(order["status"], "Pinning")
        self.assertEqual(order["fastq_r1"], "")
        row = models.get_order_inputs(oid)[0]
        self.assertEqual((row["field_name"], row["source_client"], row["source_filename"], row["status"]),
                         ("fastq_r1", "acme_lab", "tumor_R1.fastq.gz", "pending"))

    def test_lab_user_cannot_reference_another_labs_token_even_hand_crafted(self):
        c = self.client("a1")  # one client instance, reused, so the flash message survives the redirect
        r = c.post("/order/new", data=self.order_form(fastq_r1="other_lab/secret_patient_R1.fastq.gz"))
        self.assertEqual(r.status_code, 302)
        text = c.get(r.headers["Location"]).get_data(as_text=True)
        self.assertIn("not found in your authorized", text)
        self.assertEqual({o["order_id"] for o in models.list_orders(10, username="a1", role="lab_staff", lab_id=self.lab_a)}, {self.order_a})  # only the pre-existing order_a

    def test_unlinked_lab_user_cannot_submit_any_sftp_token(self):
        r = self.client("b1").post("/order/new", data=self.order_form(fastq_r1="other_lab/secret_patient_R1.fastq.gz"))
        self.assertEqual(r.status_code, 302)
        self.assertEqual({o["order_id"] for o in models.list_orders(10, username="b1", role="lab_staff", lab_id=self.lab_b)}, {self.order_b})  # only the pre-existing order_b

    def test_unlinked_lab_user_can_still_submit_with_no_fastq_at_all(self):
        r = self.client("b1").post("/order/new", data=self.order_form())
        self.assertEqual(r.status_code, 302)
        self.assertIn("/order/GR-", r.headers["Location"])

    def test_path_traversal_token_is_rejected(self):
        for bad in ("acme_lab/../other_lab/secret_patient_R1.fastq.gz",
                    "acme_lab/../../etc/passwd", "../acme_lab/tumor_R1.fastq.gz",
                    "acme_lab/tumor_R1.fastq.gz/../../other_lab/secret_patient_R1.fastq.gz"):
            r = self.client("a1").post("/order/new", data=self.order_form(fastq_r1=bad))
            self.assertEqual(r.status_code, 302)
        self.assertEqual({o["order_id"] for o in models.list_orders(10, username="a1", role="lab_staff", lab_id=self.lab_a)}, {self.order_a})  # only the pre-existing order_a

    def test_admin_can_reference_any_real_client(self):
        r = self.admin().post("/order/new", data=self.order_form(
            fastq_r1="other_lab/secret_patient_R1.fastq.gz", lab_id=str(self.lab_b)))
        self.assertEqual(r.status_code, 302)
        oid = r.headers["Location"].rsplit("/", 1)[-1]
        self.assertEqual(models.get_order(oid)["status"], "Pinning")
        row = models.get_order_inputs(oid)[0]
        self.assertEqual((row["source_client"], row["source_filename"]), ("other_lab", "secret_patient_R1.fastq.gz"))

    def test_nonexistent_file_token_is_rejected(self):
        r = self.client("a1").post("/order/new", data=self.order_form(fastq_r1="acme_lab/does_not_exist.fastq.gz"))
        self.assertEqual(r.status_code, 302)
        self.assertEqual({o["order_id"] for o in models.list_orders(10, username="a1", role="lab_staff", lab_id=self.lab_a)}, {self.order_a})  # only the pre-existing order_a

    # ── uniqueness of labs.sftp_account (review finding 2) ──────────────────
    def test_two_labs_cannot_share_the_same_sftp_account(self):
        with self.assertRaises(models.LabError):
            models.set_lab_sftp_account(self.lab_b, "acme_lab")  # already lab_a's
        self.assertIsNone(models.get_lab(self.lab_b)["sftp_account"])  # write never took effect

    def test_uniqueness_enforced_via_the_admin_route_too(self):
        adm = self.admin()
        r = adm.post(f"/admin/labs/{self.lab_b}/sftp", data={"sftp_account": "acme_lab"})
        self.assertEqual(r.status_code, 302)
        self.assertIsNone(models.get_lab(self.lab_b)["sftp_account"])
        self.assertEqual(models.get_lab(self.lab_a)["sftp_account"], "acme_lab")  # lab_a's link untouched

    def test_changing_a_labs_link_frees_the_old_name_without_leaving_a_duplicate(self):
        models.set_lab_sftp_account(self.lab_a, "acme_lab_v2")
        self.assertEqual(models.get_lab(self.lab_a)["sftp_account"], "acme_lab_v2")
        models.set_lab_sftp_account(self.lab_b, "acme_lab")  # now free -- must succeed
        self.assertEqual(models.get_lab(self.lab_b)["sftp_account"], "acme_lab")

    def test_multiple_labs_may_stay_unlinked_at_once(self):
        # NULL must never collide with NULL: two labs with no link is the normal
        # starting state for every lab and must not trip the unique index.
        models.set_lab_sftp_account(self.lab_a, "")
        models.set_lab_sftp_account(self.lab_b, "")
        self.assertIsNone(models.get_lab(self.lab_a)["sftp_account"])
        self.assertIsNone(models.get_lab(self.lab_b)["sftp_account"])

    # ── symlink hardening (review finding 3) ─────────────────────────────────
    def test_account_directory_symlinked_to_another_account_is_rejected(self):
        # evil_lab/ itself is a symlink straight to acme_lab's real account dir.
        os.symlink(Path(self.sftp_root) / "acme_lab", Path(self.sftp_root) / "evil_lab",
                  target_is_directory=True)
        models.set_lab_sftp_account(self.lab_b, "evil_lab")
        self.assertEqual(self.client("b1").get("/api/sftp-files").get_json()["files"], [])
        r = self.client("b1").post("/order/new", data=self.order_form(fastq_r1="evil_lab/tumor_R1.fastq.gz"))
        self.assertEqual(r.status_code, 302)
        self.assertEqual({o["order_id"] for o in models.list_orders(10, username="b1", role="lab_staff", lab_id=self.lab_b)}, {self.order_b})

    def test_uploads_directory_symlinked_to_another_accounts_uploads_is_rejected(self):
        # real_client/ is a genuine directory, but its uploads/ is a symlink into other_lab.
        real_client = Path(self.sftp_root) / "real_client"
        real_client.mkdir()
        os.symlink(Path(self.sftp_root) / "other_lab" / "uploads", real_client / "uploads",
                  target_is_directory=True)
        models.set_lab_sftp_account(self.lab_b, "real_client")
        self.assertEqual(self.client("b1").get("/api/sftp-files").get_json()["files"], [])
        r = self.client("b1").post("/order/new",
                                   data=self.order_form(fastq_r1="real_client/secret_patient_R1.fastq.gz"))
        self.assertEqual(r.status_code, 302)
        self.assertEqual({o["order_id"] for o in models.list_orders(10, username="b1", role="lab_staff", lab_id=self.lab_b)}, {self.order_b})

    def test_symlinked_file_inside_a_legitimate_uploads_dir_is_never_listed_or_selectable(self):
        link = Path(self.sftp_root) / "acme_lab" / "uploads" / "linked.fastq.gz"
        os.symlink(Path(self.sftp_root) / "other_lab" / "uploads" / "secret_patient_R1.fastq.gz", link)
        data = self.client("a1").get("/api/sftp-files").get_json()
        self.assertEqual({f["name"] for f in data["files"]}, {"tumor_R1.fastq.gz", "tumor_R2.fastq.gz"})
        r = self.client("a1").post("/order/new", data=self.order_form(fastq_r1="acme_lab/linked.fastq.gz"))
        self.assertEqual(r.status_code, 302)
        self.assertEqual({o["order_id"] for o in models.list_orders(10, username="a1", role="lab_staff", lab_id=self.lab_a)}, {self.order_a})

    def test_dot_and_dotdot_account_and_file_names_are_rejected_explicitly(self):
        for bad_client in (".", ".."):
            self.assertIsNone(portal_app._sftp_upload_root(bad_client))
        for bad_name in (".", ".."):
            self.assertIsNone(portal_app._resolve_sftp_token(f"acme_lab/{bad_name}", allowed_client="acme_lab"))
        with self.assertRaises(models.LabError):
            models.set_lab_sftp_account(self.lab_b, "..")

    # ── XSS hardening of the file picker (review finding 1) ──────────────────
    def test_new_order_template_never_interpolates_file_data_into_innerhtml(self):
        src = (Path(__file__).resolve().parents[1] / "portal" / "templates" / "new_order.html").read_text(encoding="utf-8")
        self.assertNotIn("tr.innerHTML", src)
        self.assertNotIn("${f.name}", src)
        self.assertNotIn("${f.client}", src)
        self.assertIn("textContent = f.name", src)
        self.assertIn("textContent = f.client", src)


class _TempSftpRoot:
    def __enter__(self):
        import tempfile
        self._d = tempfile.TemporaryDirectory()
        return self._d.name

    def __exit__(self, *a):
        self._d.cleanup()


if __name__ == "__main__":
    unittest.main()
