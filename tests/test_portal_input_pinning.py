"""Input pinning: /order/new -> "Pinning" -> a durable copy-and-hash worker ->
"Queued", with Snakemake reading only the pinned copy and re-verifying its
SHA-256 immediately before launch.

Reuses the BillingBase fixture (throw-away SQLite DB, fake labs/users) and the
fake-SFTP-root pattern from test_portal_sftp_isolation.py. Runs the pinning
worker function directly (never the real daemon thread/its sleep loop) so
tests are fast and deterministic; stability-check timing constants are
patched to milliseconds. Never touches the real portal DB, uploads, SFTP
root, or runs a real Snakemake subprocess (subprocess.Popen is mocked in the
one test that reaches _run_snakemake).
"""
import hashlib
import os
import sqlite3
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_portal_billing import BillingBase, HAVE_FLASK, LAB_A, LAB_B  # noqa: E402

if HAVE_FLASK:
    import app as portal_app
    import models
    import runner
    import sftp_paths


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


@unittest.skipUnless(HAVE_FLASK, "Flask is not installed in this Python")
class PinningBase(BillingBase):
    def setUp(self):
        super().setUp()
        self._sftp_tmpdir = tempfile.TemporaryDirectory()
        self.sftp_root = self._sftp_tmpdir.name
        self._pin_tmpdir = tempfile.TemporaryDirectory()
        self.pin_root = self._pin_tmpdir.name

        self._patches = [
            mock.patch.object(sftp_paths, "SFTP_ROOT", self.sftp_root),
            mock.patch.object(runner, "PINNED_INPUTS_DIR", self.pin_root),
            mock.patch.object(runner, "_STABILITY_CHECK_INTERVAL_SEC", 0.02),
            mock.patch.object(runner, "_STABILITY_MAX_ATTEMPTS", 5),
        ]
        for p in self._patches:
            p.start()

        d = Path(self.sftp_root) / "acme_lab" / "uploads"
        d.mkdir(parents=True)
        self.tumor_r1 = d / "tumor_R1.fastq.gz"
        self.tumor_r2 = d / "tumor_R2.fastq.gz"
        self.tumor_r1.write_bytes(b"@read1\nACGTACGT\n+\nIIIIIIII\n" * 100)
        self.tumor_r2.write_bytes(b"@read1\nTGCATGCA\n+\nIIIIIIII\n" * 120)
        models.set_lab_sftp_account(self.lab_a, "acme_lab")

    def tearDown(self):
        for p in reversed(self._patches):
            p.stop()
        self._pin_tmpdir.cleanup()
        self._sftp_tmpdir.cleanup()
        super().tearDown()

    def order_form(self, **extra):
        return {"panel_type": "comprehensive", "patient_id": "PT-1", **extra}

    def create_pinning_order(self, user="a1", **fastq):
        r = self.client(user).post("/order/new", data=self.order_form(**fastq))
        self.assertEqual(r.status_code, 302)
        return r.headers["Location"].rsplit("/", 1)[-1]


class TestOrderCreationGoesToPinning(PinningBase):
    def test_order_with_a_selected_file_starts_in_pinning_not_queued(self):
        oid = self.create_pinning_order(fastq_r1="acme_lab/tumor_R1.fastq.gz")
        order = models.get_order(oid)
        self.assertEqual(order["status"], "Pinning")
        self.assertEqual(order["fastq_r1"], "")  # no live SFTP path is ever stored
        rows = models.get_order_inputs(oid)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["field_name"], "fastq_r1")
        self.assertEqual(rows[0]["status"], "pending")
        self.assertEqual(rows[0]["source_client"], "acme_lab")
        self.assertEqual(rows[0]["source_filename"], "tumor_R1.fastq.gz")

    def test_order_with_no_files_selected_goes_straight_to_queued_unchanged(self):
        oid = self.create_pinning_order()
        order = models.get_order(oid)
        self.assertEqual(order["status"], "Queued")
        self.assertEqual(models.get_order_inputs(oid), [])


class TestPinningWorkerSuccess(PinningBase):
    def test_pinning_copies_hashes_and_moves_order_to_queued(self):
        oid = self.create_pinning_order(
            fastq_r1="acme_lab/tumor_R1.fastq.gz", fastq_r2="acme_lab/tumor_R2.fastq.gz")
        runner._pin_order_inputs(oid)

        order = models.get_order(oid)
        self.assertEqual(order["status"], "Queued")

        expected_r1 = sha256_of(self.tumor_r1)
        expected_r2 = sha256_of(self.tumor_r2)
        r1_row = next(r for r in models.get_order_inputs(oid) if r["field_name"] == "fastq_r1")
        r2_row = next(r for r in models.get_order_inputs(oid) if r["field_name"] == "fastq_r2")

        self.assertEqual(r1_row["status"], "done")
        self.assertEqual(r1_row["sha256"], expected_r1)
        self.assertEqual(r1_row["size_bytes"], self.tumor_r1.stat().st_size)
        self.assertTrue(r1_row["pinned_path"].startswith(self.pin_root))
        self.assertIn(oid, r1_row["pinned_path"])
        self.assertNotIn("tumor_R1", os.path.basename(r1_row["pinned_path"]))  # neutral filename

        self.assertEqual(order["fastq_r1"], r1_row["pinned_path"])
        self.assertEqual(order["fastq_r2"], r2_row["pinned_path"])
        self.assertEqual(sha256_of(order["fastq_r1"]), expected_r1)
        self.assertEqual(sha256_of(order["fastq_r2"]), expected_r2)

        # The pinned copy is a real independent file, not a hard link: mutating
        # the original SFTP source afterward must not touch the pinned bytes.
        self.tumor_r1.write_bytes(b"TAMPERED-AFTER-PIN")
        self.assertEqual(sha256_of(order["fastq_r1"]), expected_r1)

    def test_pinned_filename_preserves_gzip_suffix(self):
        oid = self.create_pinning_order(fastq_r1="acme_lab/tumor_R1.fastq.gz")
        runner._pin_order_inputs(oid)
        row = models.get_order_inputs(oid)[0]
        self.assertTrue(row["pinned_path"].endswith("fastq_r1.fastq.gz"))


class TestPinningFailureCases(PinningBase):
    def test_source_deleted_before_pinning_fails_the_whole_order_no_run(self):
        # Deleted AFTER the order is created (selection succeeded) but BEFORE
        # the pinning worker gets to it -- the exact TOCTOU window this
        # feature exists to close. (Deleting it before creation would instead
        # be caught by the pre-existing selection-time check in app.py, which
        # is a different, already-tested code path.)
        oid = self.create_pinning_order(
            fastq_r1="acme_lab/tumor_R1.fastq.gz", fastq_r2="acme_lab/tumor_R2.fastq.gz")
        os.remove(self.tumor_r1)
        runner._pin_order_inputs(oid)
        order = models.get_order(oid)
        self.assertEqual(order["status"], "Failed")
        self.assertIn("fastq_r1", order["error_msg"])
        # The sibling field (r2) that WAS still valid must be cleaned up too --
        # an order never runs partially pinned.
        r2_row = next(r for r in models.get_order_inputs(oid) if r["field_name"] == "fastq_r2")
        self.assertEqual(r2_row["status"], "pending")
        self.assertIsNone(r2_row["pinned_path"])
        self.assertEqual(list(Path(self.pin_root).rglob("*.fastq.gz")), [])  # nothing left on disk

    def test_source_that_keeps_changing_never_stabilizes_and_fails_cleanly(self):
        oid = self.create_pinning_order(fastq_r1="acme_lab/tumor_R1.fastq.gz")

        # Keep touching the file so its mtime never settles across two checks.
        import threading as th
        keep_going = th.Event()
        keep_going.set()

        def churn():
            while keep_going.is_set():
                self.tumor_r1.write_bytes(self.tumor_r1.read_bytes() + b"x")
                time.sleep(0.01)

        t = th.Thread(target=churn, daemon=True)
        t.start()
        try:
            runner._pin_order_inputs(oid)
        finally:
            keep_going.clear()
            t.join(timeout=2)

        order = models.get_order(oid)
        self.assertEqual(order["status"], "Failed")
        self.assertIn("did not stop changing", order["error_msg"])

    def test_source_changed_between_stability_check_and_copy_is_rejected(self):
        oid = self.create_pinning_order(fastq_r1="acme_lab/tumor_R1.fastq.gz")
        # Directly exercise _copy_and_hash with a fingerprint that no longer
        # matches reality -- deterministic stand-in for "it changed between
        # the stability check and copy start" without relying on real-time
        # races. Only the size field is wrong; everything else is correct, so
        # this specifically exercises the fingerprint mismatch path, not a
        # malformed tuple.
        real = runner._fingerprint(self.tumor_r1.stat())
        bogus = (real[0] + 999,) + real[1:]
        ok, result = runner._copy_and_hash(
            oid, "fastq_r1", str(self.tumor_r1), "tumor_R1.fastq.gz", bogus)
        self.assertFalse(ok)
        self.assertIn("changed between the stability check", result)
        # No partial file left behind.
        dest_dir = Path(self.pin_root) / oid
        self.assertEqual(list(dest_dir.glob("fastq_r1.*")) if dest_dir.exists() else [], [])

    def test_insufficient_disk_space_fails_before_any_copy(self):
        oid = self.create_pinning_order(fastq_r1="acme_lab/tumor_R1.fastq.gz")
        with mock.patch("shutil.disk_usage") as du:
            du.return_value = mock.Mock(total=10**12, used=0, free=1)  # ~1 byte free
            runner._pin_order_inputs(oid)
        order = models.get_order(oid)
        self.assertEqual(order["status"], "Failed")
        self.assertIn("insufficient disk space", order["error_msg"])
        self.assertEqual(list(Path(self.pin_root).rglob("*.fastq.gz")), [])

    def test_symlink_swapped_in_after_selection_is_still_rejected_at_pin_time(self):
        # Selection-time resolution happened against a real file; simulate the
        # account's uploads dir being replaced by a symlink to another
        # account's directory before the pinning worker actually runs.
        oid = self.create_pinning_order(fastq_r1="acme_lab/tumor_R1.fastq.gz")
        other = Path(self.sftp_root) / "other_lab" / "uploads"
        other.mkdir(parents=True)
        acme_uploads = Path(self.sftp_root) / "acme_lab" / "uploads"
        import shutil as sh
        sh.rmtree(acme_uploads)
        os.symlink(other, acme_uploads, target_is_directory=True)

        runner._pin_order_inputs(oid)
        order = models.get_order(oid)
        self.assertEqual(order["status"], "Failed")
        self.assertIn("not found", order["error_msg"])


class TestRestartResumption(PinningBase):
    def test_stuck_copying_row_is_reset_and_its_partial_file_removed_on_startup(self):
        oid = self.create_pinning_order(fastq_r1="acme_lab/tumor_R1.fastq.gz")
        models.mark_input_copying(oid, "fastq_r1")
        dest_dir = Path(self.pin_root) / oid
        dest_dir.mkdir(parents=True)
        partial = dest_dir / "fastq_r1.fastq.gz.part-garbage"
        # Match the cleanup's own prefix rule ("field_name" + ".").
        partial = dest_dir / "fastq_r1.fastq.gz"
        partial.write_bytes(b"PARTIAL-INTERRUPTED-COPY")

        runner._reset_stuck_pinning()

        row = models.get_order_inputs(oid)[0]
        self.assertEqual(row["status"], "pending")
        self.assertFalse(partial.exists())

    def test_worker_then_resumes_and_completes_the_reset_order(self):
        oid = self.create_pinning_order(fastq_r1="acme_lab/tumor_R1.fastq.gz")
        models.mark_input_copying(oid, "fastq_r1")
        runner._reset_stuck_pinning()
        runner._pin_order_inputs(oid)
        self.assertEqual(models.get_order(oid)["status"], "Queued")


class TestPreflightHashRecheckBeforeSnakemake(PinningBase):
    def test_tampered_pinned_file_blocks_the_run_and_snakemake_is_never_launched(self):
        oid = self.create_pinning_order(fastq_r1="acme_lab/tumor_R1.fastq.gz")
        runner._pin_order_inputs(oid)
        self.assertEqual(models.get_order(oid)["status"], "Queued")

        pinned_path = models.get_order(oid)["fastq_r1"]
        with open(pinned_path, "ab") as f:
            f.write(b"TAMPERED-AFTER-PINNING-BEFORE-RUN")

        with mock.patch("subprocess.Popen") as popen:
            runner._run_snakemake(oid)
            popen.assert_not_called()

        order = models.get_order(oid)
        self.assertEqual(order["status"], "Failed")
        self.assertIn("Input integrity check failed", order["error_msg"])
        self.assertIn("fastq_r1", order["error_msg"])

    def test_untampered_pinned_file_passes_the_recheck_and_snakemake_is_launched(self):
        oid = self.create_pinning_order(fastq_r1="acme_lab/tumor_R1.fastq.gz")
        runner._pin_order_inputs(oid)

        class FakeProcess:
            pid = 12345
            def wait(self):
                return 0

        with mock.patch("subprocess.Popen", return_value=FakeProcess()) as popen, \
             mock.patch.object(runner, "_find_report", return_value=""):
            runner._run_snakemake(oid)
            popen.assert_called_once()

        self.assertEqual(models.get_order(oid)["status"], "Done")


class TestRetryRouting(PinningBase):
    def test_retry_of_a_pipeline_failure_with_inputs_already_pinned_goes_straight_to_queued(self):
        oid = self.create_pinning_order(fastq_r1="acme_lab/tumor_R1.fastq.gz")
        runner._pin_order_inputs(oid)
        models.update_status(oid, "Failed", error_msg="Snakemake exit code 1")

        r = self.client("a1").post(f"/order/{oid}/retry")
        self.assertEqual(r.status_code, 302)
        self.assertEqual(models.get_order(oid)["status"], "Queued")
        # The already-pinned file/hash must survive a pipeline-only retry untouched.
        self.assertEqual(len(models.get_order_inputs(oid)), 1)
        self.assertEqual(models.get_order_inputs(oid)[0]["status"], "done")

    def test_retry_of_a_pinning_failure_routes_back_through_pinning_not_queued(self):
        oid = self.create_pinning_order(fastq_r1="acme_lab/tumor_R1.fastq.gz")
        os.remove(self.tumor_r1)  # after selection succeeded, before pinning runs
        runner._pin_order_inputs(oid)
        self.assertEqual(models.get_order(oid)["status"], "Failed")

        r = self.client("a1").post(f"/order/{oid}/retry")
        self.assertEqual(r.status_code, 302)
        self.assertEqual(models.get_order(oid)["status"], "Pinning")
        self.assertEqual(models.get_order_inputs(oid)[0]["status"], "pending")

    def test_other_labs_user_cannot_retry_someone_elses_order(self):
        oid = self.create_pinning_order(fastq_r1="acme_lab/tumor_R1.fastq.gz")
        runner._pin_order_inputs(oid)
        models.update_status(oid, "Failed", error_msg="x")
        self.assertEqual(self.client("b1").post(f"/order/{oid}/retry").status_code, 403)


class _FlakyConn:
    """Wraps one real sqlite3 connection so a chosen INSERT can be made to
    raise partway through a transaction -- used to prove new_order_with_inputs()
    really is atomic, not just documented as such. __enter__/__exit__ forward
    to the real connection's own (commit on clean exit, rollback on
    exception), which is genuine sqlite3 behaviour, not simulated."""
    def __init__(self, real, fail_on_sql_substring):
        self._real = real
        self._fail_on = fail_on_sql_substring

    def __enter__(self):
        self._real.__enter__()
        return self

    def __exit__(self, exc_type, exc, tb):
        return self._real.__exit__(exc_type, exc, tb)

    def execute(self, sql, *args, **kwargs):
        if self._fail_on in sql:
            raise sqlite3.OperationalError("simulated mid-transaction failure")
        return self._real.execute(sql, *args, **kwargs)


class TestAtomicOrderCreation(PinningBase):
    """Finding 2: order creation and every selected order_inputs row must
    commit together in one transaction, never some-but-not-all."""

    def _order_and_input_counts(self):
        con = sqlite3.connect(self.db)
        try:
            n_orders = con.execute("SELECT COUNT(*) FROM orders").fetchone()[0]
            n_inputs = con.execute("SELECT COUNT(*) FROM order_inputs").fetchone()[0]
            return n_orders, n_inputs
        finally:
            con.close()

    def test_failure_on_the_first_order_inputs_insert_leaves_neither_order_nor_any_input(self):
        before = self._order_and_input_counts()
        real_conn = models._conn

        def flaky():
            return _FlakyConn(real_conn(), fail_on_sql_substring="INSERT INTO order_inputs")

        with mock.patch.object(models, "_conn", side_effect=flaky):
            with self.assertRaises(sqlite3.OperationalError):
                models.new_order_with_inputs(
                    patient_id="PT-ATOMIC", patient_name="", sex="", dob="", tumor_type="",
                    panel_type="comprehensive",
                    inputs=[("fastq_r1", "acme_lab", "tumor_R1.fastq.gz"),
                            ("fastq_r2", "acme_lab", "tumor_R2.fastq.gz")],
                    created_by="a1", lab_id=self.lab_a)

        # Not "the order exists but with 0 or 1 inputs" -- NEITHER the order
        # row NOR any order_inputs row survives; the whole attempt vanished.
        self.assertEqual(self._order_and_input_counts(), before)

    def test_failure_on_the_second_order_inputs_insert_still_rolls_back_the_first(self):
        # The failure this time happens AFTER one order_inputs row already
        # inserted successfully inside the same transaction -- proving the
        # rollback undoes an already-succeeded sibling insert too, not just
        # the order row.
        before = self._order_and_input_counts()
        real_conn = models._conn
        calls = {"order_inputs": 0}

        def flaky():
            real = real_conn()
            wrapper = _FlakyConn(real, fail_on_sql_substring="__never__")
            real_execute = real.execute

            def execute(sql, *args, **kwargs):
                if "INSERT INTO order_inputs" in sql:
                    calls["order_inputs"] += 1
                    if calls["order_inputs"] == 2:
                        raise sqlite3.OperationalError("simulated mid-transaction failure")
                return real_execute(sql, *args, **kwargs)
            wrapper.execute = execute
            return wrapper

        with mock.patch.object(models, "_conn", side_effect=flaky):
            with self.assertRaises(sqlite3.OperationalError):
                models.new_order_with_inputs(
                    patient_id="PT-ATOMIC2", patient_name="", sex="", dob="", tumor_type="",
                    panel_type="comprehensive",
                    inputs=[("fastq_r1", "acme_lab", "tumor_R1.fastq.gz"),
                            ("fastq_r2", "acme_lab", "tumor_R2.fastq.gz")],
                    created_by="a1", lab_id=self.lab_a)

        self.assertEqual(self._order_and_input_counts(), before)

    def test_all_selected_fields_committed_together_on_success(self):
        oid = models.new_order_with_inputs(
            patient_id="PT-OK", patient_name="", sex="", dob="", tumor_type="",
            panel_type="comprehensive",
            inputs=[("fastq_r1", "acme_lab", "tumor_R1.fastq.gz"),
                    ("fastq_r2", "acme_lab", "tumor_R2.fastq.gz"),
                    ("fastq_normal_r1", "acme_lab", "normal_R1.fastq.gz"),
                    ("fastq_normal_r2", "acme_lab", "normal_R2.fastq.gz")],
            created_by="a1", lab_id=self.lab_a)
        order = models.get_order(oid)
        self.assertEqual(order["status"], "Pinning")
        rows = {r["field_name"]: r for r in models.get_order_inputs(oid)}
        self.assertEqual(set(rows), {"fastq_r1", "fastq_r2", "fastq_normal_r1", "fastq_normal_r2"})
        self.assertTrue(all(r["status"] == "pending" for r in rows.values()))

    def test_rejects_being_called_with_no_inputs(self):
        with self.assertRaises(ValueError):
            models.new_order_with_inputs(
                patient_id="PT-EMPTY", patient_name="", sex="", dob="", tumor_type="",
                panel_type="comprehensive", inputs=[], created_by="a1", lab_id=self.lab_a)

    def test_http_route_still_creates_the_order_and_all_its_inputs_together(self):
        oid = self.create_pinning_order(
            fastq_r1="acme_lab/tumor_R1.fastq.gz", fastq_r2="acme_lab/tumor_R2.fastq.gz")
        self.assertEqual(models.get_order(oid)["status"], "Pinning")
        self.assertEqual({r["field_name"] for r in models.get_order_inputs(oid)},
                         {"fastq_r1", "fastq_r2"})


class TestPathIsWithinSftpRoot(PinningBase):
    """Direct unit tests of sftp_paths.path_is_within_sftp_root(), the guard
    that decides whether a legacy fastq_* value is untrusted. Covers both
    directions: a path that lexically sits inside SFTP_ROOT but resolves
    outside it (via a symlink placed inside the tree -- the case that a
    resolved-target-only check would previously have missed), and a path
    that lexically sits outside SFTP_ROOT but resolves into it (via a
    symlink placed outside the tree -- the pre-existing, still-required
    direction)."""

    def setUp(self):
        super().setUp()
        self.outside_dir = Path(self._pin_tmpdir.name).parent / "outside_sftp_root"
        self.outside_dir.mkdir(exist_ok=True)

    def test_a_real_file_directly_inside_the_tree_is_within(self):
        self.assertTrue(sftp_paths.path_is_within_sftp_root(str(self.tumor_r1)))

    def test_a_genuinely_unrelated_outside_path_is_not_within(self):
        # "Preserve existing behavior for genuinely outside paths."
        outside_file = self.outside_dir / "old_upload_R1.fastq.gz"
        outside_file.write_bytes(b"@r\nACGT\n+\nIIII\n")
        self.assertFalse(sftp_paths.path_is_within_sftp_root(str(outside_file)))

    def test_empty_or_none_path_is_not_within(self):
        self.assertFalse(sftp_paths.path_is_within_sftp_root(""))
        self.assertFalse(sftp_paths.path_is_within_sftp_root(None))

    def test_file_symlink_inside_the_tree_pointing_outside_is_still_within(self):
        # The exact scenario this fix addresses: the symlink's own location
        # is inside SFTP_ROOT (so the SFTP client can repoint or replace it
        # at any time) even though it currently resolves to a file that is
        # genuinely outside the tree.
        outside_target = self.outside_dir / "secret_elsewhere.fastq.gz"
        outside_target.write_bytes(b"SHOULD-BE-TREATED-AS-LIVE-SFTP-CONTROLLED")
        uploads = Path(self.sftp_root) / "acme_lab" / "uploads"
        link = uploads / "looks_like_a_normal_upload.fastq.gz"
        os.symlink(outside_target, link)

        self.assertTrue(sftp_paths.path_is_within_sftp_root(str(link)))
        # And, for contrast, the resolved target alone is correctly NOT
        # within the tree -- confirming this is genuinely the lexical check
        # doing the work, not a coincidence.
        self.assertFalse(sftp_paths.path_is_within_sftp_root(str(outside_target)))

    def test_directory_symlink_inside_the_tree_pointing_outside_is_still_within(self):
        # A whole subdirectory under SFTP_ROOT is itself a symlink to a
        # directory outside the tree; a file reached *through* it is still
        # lexically inside SFTP_ROOT as a path string, even though every
        # component past the symlink physically resolves elsewhere.
        outside_subdir = self.outside_dir / "real_directory_elsewhere"
        outside_subdir.mkdir()
        (outside_subdir / "tumor_R1.fastq.gz").write_bytes(b"SHOULD-BE-TREATED-AS-LIVE-SFTP-CONTROLLED")
        uploads = Path(self.sftp_root) / "acme_lab" / "uploads"
        linked_subdir = uploads / "linked_elsewhere"
        os.symlink(outside_subdir, linked_subdir, target_is_directory=True)

        path_through_link = str(linked_subdir / "tumor_R1.fastq.gz")
        self.assertTrue(sftp_paths.path_is_within_sftp_root(path_through_link))
        self.assertFalse(sftp_paths.path_is_within_sftp_root(
            str(outside_subdir / "tumor_R1.fastq.gz")))

    def test_symlink_outside_the_tree_pointing_into_it_is_still_within(self):
        # The complementary, pre-existing direction: the path is lexically
        # OUTSIDE SFTP_ROOT, but resolves INTO it -- still unsafe, and must
        # stay caught by the resolved-target check.
        link = self.outside_dir / "alias_of_a_real_upload.fastq.gz"
        os.symlink(self.tumor_r1, link)
        self.assertTrue(sftp_paths.path_is_within_sftp_root(str(link)))

    def test_legacy_order_via_a_symlink_pointing_outside_is_blocked_not_run(self):
        # End-to-end: the same scenario as the two symlink unit tests above,
        # exercised through the actual order-running guard, not just the
        # bare path function.
        outside_target = self.outside_dir / "secret_elsewhere.fastq.gz"
        outside_target.write_bytes(b"SHOULD-NEVER-BE-READ-BY-THE-PIPELINE")
        uploads = Path(self.sftp_root) / "acme_lab" / "uploads"
        link = uploads / "sneaky.fastq.gz"
        os.symlink(outside_target, link)

        oid = models.new_order(
            patient_id="PT-SYMLINK", patient_name="", sex="", dob="", tumor_type="",
            panel_type="comprehensive", fastq_r1=str(link), fastq_r2="",
            fastq_normal_r1="", fastq_normal_r2="",
            created_by="a1", lab_id=self.lab_a, status="Queued")

        with mock.patch("subprocess.Popen") as popen:
            runner._run_snakemake(oid)
            popen.assert_not_called()

        order = models.get_order(oid)
        self.assertEqual(order["status"], "Failed")
        self.assertIn("fastq_r1", order["error_msg"])
        self.assertIn("SFTP", order["error_msg"])


class TestLegacyOrdersWithoutInputRows(PinningBase):
    """Finding 3: an order created before input pinning existed has no
    order_inputs rows at all. Its fastq_* value must still be refused if it
    resolves into the live SFTP tree (nothing re-verifies it), and must still
    be preserved/run normally if it points somewhere else entirely."""

    def _legacy_order(self, fastq_r1, status="Queued"):
        return models.new_order(
            patient_id="PT-LEGACY", patient_name="", sex="", dob="", tumor_type="",
            panel_type="comprehensive", fastq_r1=fastq_r1, fastq_r2="",
            fastq_normal_r1="", fastq_normal_r2="",
            created_by="a1", lab_id=self.lab_a, status=status)

    def test_legacy_fastq_pointing_into_the_live_sftp_tree_is_refused_not_run(self):
        oid = self._legacy_order(str(self.tumor_r1))  # a real file, genuinely inside sftp_root
        self.assertEqual(models.get_order_inputs(oid), [])  # the untracked case this finding is about

        with mock.patch("subprocess.Popen") as popen:
            runner._run_snakemake(oid)
            popen.assert_not_called()

        order = models.get_order(oid)
        self.assertEqual(order["status"], "Failed")
        self.assertIn("fastq_r1", order["error_msg"])
        self.assertIn("SFTP", order["error_msg"])

    def test_genuinely_input_free_legacy_order_with_a_non_sftp_path_still_runs(self):
        # A path that is NOT inside sftp_root at all (e.g. the portal's own
        # local upload area from before SFTP selection existed) must be
        # completely unaffected by this finding's new check.
        other_dir = Path(self._pin_tmpdir.name).parent / "legacy_local_uploads"
        other_dir.mkdir(exist_ok=True)
        legacy_fastq = other_dir / "old_upload_R1.fastq.gz"
        legacy_fastq.write_bytes(b"@r\nACGT\n+\nIIII\n")
        self.assertFalse(sftp_paths.path_is_within_sftp_root(str(legacy_fastq)))
        oid = self._legacy_order(str(legacy_fastq))

        class FakeProcess:
            pid = 1
            def wait(self):
                return 0

        with mock.patch("subprocess.Popen", return_value=FakeProcess()) as popen, \
             mock.patch.object(runner, "_find_report", return_value=""):
            runner._run_snakemake(oid)
            popen.assert_called_once()

        self.assertEqual(models.get_order(oid)["status"], "Done")

    def test_order_with_no_fastq_selected_at_all_is_unaffected(self):
        oid = self._legacy_order("")  # no input of any kind -- nothing to check
        ok, reason = runner._check_no_untracked_live_sftp_paths(dict(models.get_order(oid)))
        self.assertTrue(ok)
        self.assertEqual(reason, "")

    def test_a_field_already_tracked_and_done_is_never_re_flagged_even_if_its_path_is_inside_sftp_root(self):
        # Defence in depth / no false positive: if a field DOES have a 'done'
        # order_inputs row, this check must not re-examine orders.<field> at
        # all -- that path is exactly the pinned copy under PINNED_INPUTS_DIR,
        # which legitimately is not inside SFTP_ROOT, but this test makes the
        # "already tracked" branch explicit rather than relying on that.
        oid = self.create_pinning_order(fastq_r1="acme_lab/tumor_R1.fastq.gz")
        runner._pin_order_inputs(oid)
        self.assertEqual(models.get_order(oid)["status"], "Queued")
        ok, reason = runner._check_no_untracked_live_sftp_paths(dict(models.get_order(oid)))
        self.assertTrue(ok)
        self.assertEqual(reason, "")

    def test_retry_of_a_blocked_legacy_order_goes_to_queued_and_fails_closed_again_not_silently(self):
        oid = self._legacy_order(str(self.tumor_r1))
        with mock.patch("subprocess.Popen") as popen:
            runner._run_snakemake(oid)
            popen.assert_not_called()
        self.assertEqual(models.get_order(oid)["status"], "Failed")

        r = self.client("a1").post(f"/order/{oid}/retry")
        self.assertEqual(r.status_code, 302)
        # Zero order_inputs rows -> order_inputs_all_done() is vacuously true
        # -> retry_order() sends it straight to Queued, not through Pinning
        # (there is nothing recorded to pin). Documented, expected behaviour
        # for this case -- not a data rewrite, and not a silent bypass,
        # because running it again hits the exact same guard again below.
        self.assertEqual(models.get_order(oid)["status"], "Queued")

        with mock.patch("subprocess.Popen") as popen:
            runner._run_snakemake(oid)
            popen.assert_not_called()
        order = models.get_order(oid)
        self.assertEqual(order["status"], "Failed")
        self.assertIn("fastq_r1", order["error_msg"])
        self.assertIn("SFTP", order["error_msg"])
        # No order data was rewritten or deleted by any of this.
        self.assertEqual(order["patient_id"], "PT-LEGACY")
        self.assertEqual(order["fastq_r1"], str(self.tumor_r1))


class TestStrengthenedSourceChangeDetection(PinningBase):
    """Finding 4: the copy check must catch more than a plain size+mtime
    comparison would -- same-size replacement (via file identity) and
    in-place modification during the copy itself (via fstat on the open fd,
    both before and after the read) -- and must refuse a symlink swapped in
    for the source file at the exact moment of opening it."""

    def test_fingerprint_distinguishes_two_same_size_same_mtime_files_by_identity(self):
        # Two genuinely different files, forced to the exact same size and
        # the exact same nanosecond mtime -- the one case a plain (size,
        # mtime) comparison cannot tell apart. Only file identity (dev, ino)
        # can, and this proves _fingerprint() actually uses it.
        d = Path(self.sftp_root) / "acme_lab" / "uploads"
        a, b = d / "twin_a.fastq.gz", d / "twin_b.fastq.gz"
        a.write_bytes(b"SAME-SIZE-CONTENT-A")
        b.write_bytes(b"SAME-SIZE-CONTENT-B")
        self.assertEqual(a.stat().st_size, b.stat().st_size)
        same_ns = 1_700_000_000_123456789
        os.utime(a, ns=(same_ns, same_ns))
        os.utime(b, ns=(same_ns, same_ns))
        fa, fb = runner._fingerprint(a.stat()), runner._fingerprint(b.stat())
        self.assertEqual(fa[0], fb[0])                    # size: equal
        self.assertEqual(fa[1], fb[1])                    # mtime_ns: forced equal
        self.assertNotEqual((fa[3], fa[4]), (fb[3], fb[4]))  # (dev, ino): still distinguishes them
        self.assertNotEqual(fa, fb)

    def test_same_size_source_replacement_is_rejected_by_copy_and_hash(self):
        oid = self.create_pinning_order(fastq_r1="acme_lab/tumor_R1.fastq.gz")
        original_fp = runner._fingerprint(self.tumor_r1.stat())
        original_bytes = self.tumor_r1.read_bytes()
        replacement = bytes(byte ^ 0xFF for byte in original_bytes)  # same length, different content
        self.assertEqual(len(replacement), len(original_bytes))
        self.tumor_r1.unlink()
        self.tumor_r1.write_bytes(replacement)  # a fresh file -- new inode, same size

        new_fp = runner._fingerprint(self.tumor_r1.stat())
        self.assertEqual(new_fp[0], original_fp[0])  # same size (the scenario under test)
        self.assertNotEqual(new_fp, original_fp)     # yet still detected as a different file

        ok, result = runner._copy_and_hash(
            oid, "fastq_r1", str(self.tumor_r1), "tumor_R1.fastq.gz", original_fp)
        self.assertFalse(ok)
        self.assertIn("changed between the stability check", result)
        dest_dir = Path(self.pin_root) / oid
        self.assertEqual(list(dest_dir.glob("fastq_r1.*")) if dest_dir.exists() else [], [])

    def test_in_place_content_edit_during_copy_is_detected_by_post_copy_fingerprint(self):
        oid = self.create_pinning_order(fastq_r1="acme_lab/tumor_R1.fastq.gz")
        fp = runner._fingerprint(self.tumor_r1.stat())
        real_fdopen = os.fdopen
        mutated = {"done": False}
        source_path = self.tumor_r1

        class _MutatingReader:
            """Delegates to the real wrapped file, but on its first read()
            (i.e. right after the copy has consumed the first chunk) edits
            the source's content in place through a SEPARATE, path-based
            file object -- the exact "content changed under an fd that is
            already open" race this finding's post-copy fstat exists to
            catch."""
            def __init__(self, real_file):
                self._real = real_file

            def read(self, n=-1):
                chunk = self._real.read(n)
                if not mutated["done"]:
                    mutated["done"] = True
                    with open(source_path, "r+b") as f:
                        f.seek(0)
                        f.write(b"EDITED-DURING-COPY"[:max(1, len(chunk))])
                return chunk

            def __enter__(self):
                return self

            def __exit__(self, *exc):
                self._real.close()
                return False

        def fake_fdopen(fd, mode, *a, **kw):
            if fd_holder.get("fd") == fd:
                return _MutatingReader(real_fdopen(fd, mode, *a, **kw))
            return real_fdopen(fd, mode, *a, **kw)

        fd_holder = {}
        real_open = os.open

        def spying_open(path, flags, *a, **kw):
            f = real_open(path, flags, *a, **kw)
            if path == str(source_path):
                fd_holder["fd"] = f
            return f

        with mock.patch("os.open", side_effect=spying_open), \
             mock.patch("os.fdopen", side_effect=fake_fdopen):
            ok, result = runner._copy_and_hash(
                oid, "fastq_r1", str(self.tumor_r1), "tumor_R1.fastq.gz", fp)

        self.assertTrue(mutated["done"], "test setup did not actually trigger the in-place edit")
        self.assertFalse(ok)
        self.assertIn("changed while it was being pinned", result)
        dest_dir = Path(self.pin_root) / oid
        self.assertEqual(list(dest_dir.glob("fastq_r1.*")) if dest_dir.exists() else [], [])

    @unittest.skipUnless(hasattr(os, "O_NOFOLLOW"), "O_NOFOLLOW is POSIX-only")
    def test_file_swapped_for_a_symlink_after_the_stability_check_is_rejected_at_open_time(self):
        oid = self.create_pinning_order(fastq_r1="acme_lab/tumor_R1.fastq.gz")
        fp = runner._fingerprint(self.tumor_r1.stat())
        # Simulate the source being replaced by a symlink AFTER the
        # stability check already captured fp, but before the copy opens it
        # -- narrower and later than the whole-directory symlink swap the
        # pre-existing test above covers, and the specific TOCTOU window
        # O_NOFOLLOW exists to close.
        other = Path(self.sftp_root) / "other_lab" / "uploads"
        other.mkdir(parents=True, exist_ok=True)
        secret = other / "secret.fastq.gz"
        secret.write_bytes(b"SHOULD-NEVER-BE-READ-OR-PINNED")
        self.tumor_r1.unlink()
        os.symlink(secret, self.tumor_r1)

        ok, result = runner._copy_and_hash(
            oid, "fastq_r1", str(self.tumor_r1), "tumor_R1.fastq.gz", fp)
        self.assertFalse(ok)
        self.assertIn("could not be opened safely", result)
        dest_dir = Path(self.pin_root) / oid
        self.assertEqual(list(dest_dir.glob("fastq_r1.*")) if dest_dir.exists() else [], [])


if __name__ == "__main__":
    unittest.main()
