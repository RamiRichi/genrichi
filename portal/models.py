"""GenRichi Portal — SQLite database models"""

import re
import sqlite3
import uuid
from datetime import datetime
from werkzeug.security import generate_password_hash, check_password_hash
from config import DB_PATH


def _conn():
    conn = sqlite3.connect(DB_PATH, check_same_thread=False, timeout=30)
    conn.row_factory = sqlite3.Row
    return conn


def init_db():
    with _conn() as conn:
        conn.executescript("""
        CREATE TABLE IF NOT EXISTS users (
            id            INTEGER PRIMARY KEY AUTOINCREMENT,
            username      TEXT    UNIQUE NOT NULL,
            password_hash TEXT    NOT NULL,
            role          TEXT    DEFAULT 'lab_staff',
            full_name     TEXT    DEFAULT '',
            email         TEXT    DEFAULT '',
            active        INTEGER DEFAULT 1,
            created_at    TEXT    NOT NULL
        );

        CREATE TABLE IF NOT EXISTS orders (
            id              INTEGER PRIMARY KEY AUTOINCREMENT,
            order_id        TEXT    UNIQUE NOT NULL,
            patient_id      TEXT    NOT NULL,
            patient_name    TEXT    DEFAULT '',
            sex             TEXT    DEFAULT '',
            dob             TEXT    DEFAULT '',
            tumor_type      TEXT    DEFAULT '',
            panel_type      TEXT    NOT NULL,
            status          TEXT    DEFAULT 'Queued',
            created_at      TEXT    NOT NULL,
            started_at      TEXT,
            finished_at     TEXT,
            fastq_r1        TEXT    DEFAULT '',
            fastq_r2        TEXT    DEFAULT '',
            fastq_normal_r1 TEXT    DEFAULT '',
            fastq_normal_r2 TEXT    DEFAULT '',
            report_path     TEXT    DEFAULT '',
            log_path        TEXT    DEFAULT '',
            notes           TEXT    DEFAULT '',
            error_msg       TEXT    DEFAULT '',
            notify_email    TEXT    DEFAULT '',
            pid             INTEGER DEFAULT NULL,
            created_by      TEXT    DEFAULT 'admin'
        );

        CREATE TABLE IF NOT EXISTS labs (
            id            INTEGER PRIMARY KEY AUTOINCREMENT,
            legal_name    TEXT    NOT NULL UNIQUE COLLATE NOCASE,
            address_line1 TEXT    NOT NULL,
            address_line2 TEXT    NOT NULL DEFAULT '',
            postal_code   TEXT    NOT NULL,
            city          TEXT    NOT NULL,
            country       TEXT    NOT NULL DEFAULT 'Deutschland',
            created_at    TEXT    NOT NULL,
            updated_at    TEXT    NOT NULL
        );

        CREATE TABLE IF NOT EXISTS invoice_counters (
            year    INTEGER PRIMARY KEY,
            last_no INTEGER NOT NULL
        );

        -- An invoice row exists only once an admin has explicitly issued it.
        -- Recipient and issuer are stored as a snapshot so a later edit of
        -- the lab profile never rewrites an issued invoice.
        CREATE TABLE IF NOT EXISTS invoices (
            id             INTEGER PRIMARY KEY AUTOINCREMENT,
            invoice_no     TEXT    NOT NULL UNIQUE,
            order_id       TEXT    NOT NULL UNIQUE,
            lab_id         INTEGER NOT NULL,
            invoice_date   TEXT    NOT NULL,
            service_date   TEXT    NOT NULL,
            description    TEXT    NOT NULL,
            net_cents      INTEGER NOT NULL,
            recipient_json TEXT    NOT NULL,
            issuer_json    TEXT    NOT NULL,
            issued_by      TEXT    NOT NULL,
            created_at     TEXT    NOT NULL
        );

        CREATE TABLE IF NOT EXISTS audit_log (
            id         INTEGER PRIMARY KEY AUTOINCREMENT,
            order_id   TEXT,
            event      TEXT,
            detail     TEXT,
            ts         TEXT
        );

        -- One row per FASTQ field an order actually selected (0-4 rows; an
        -- order with none needing no pinning has none). source_client/
        -- source_filename are the exact SFTP token components recorded at
        -- selection time, re-resolved independently by the pinning worker --
        -- never trusted as still valid just because they were valid once.
        -- pinned_path/sha256/size_bytes are NULL until status='done'.
        CREATE TABLE IF NOT EXISTS order_inputs (
            id               INTEGER PRIMARY KEY AUTOINCREMENT,
            order_id         TEXT    NOT NULL,
            field_name       TEXT    NOT NULL,
            source_client    TEXT    NOT NULL,
            source_filename  TEXT    NOT NULL,
            pinned_path      TEXT,
            sha256           TEXT,
            size_bytes       INTEGER,
            status           TEXT    NOT NULL DEFAULT 'pending',
            error_msg        TEXT    DEFAULT '',
            created_at       TEXT    NOT NULL,
            pinned_at        TEXT,
            UNIQUE(order_id, field_name)
        );
        CREATE INDEX IF NOT EXISTS idx_order_inputs_order  ON order_inputs(order_id);
        CREATE INDEX IF NOT EXISTS idx_order_inputs_status ON order_inputs(status);
        """)
    # Migrate existing DB: add missing columns without breaking data
    _migrate()
    # Ensure admin user exists
    _ensure_admin()


def _migrate():
    """Safe migrations — add new columns to existing tables if absent."""
    with _conn() as conn:
        existing_orders = {row[1] for row in conn.execute("PRAGMA table_info(orders)")}
        for col, ddl in [
            ("notify_email", "TEXT DEFAULT ''"),
            ("pid",          "INTEGER DEFAULT NULL"),
            ("created_by",   "TEXT DEFAULT 'admin'"),
            ("lab_id",       "INTEGER DEFAULT NULL REFERENCES labs(id)"),
        ]:
            if col not in existing_orders:
                conn.execute(f"ALTER TABLE orders ADD COLUMN {col} {ddl}")

        existing_users = {row[1] for row in conn.execute("PRAGMA table_info(users)")}
        for col, ddl in [
            ("full_name", "TEXT DEFAULT ''"),
            ("email",     "TEXT DEFAULT ''"),
            ("active",    "INTEGER DEFAULT 1"),
            ("lab_id",    "INTEGER DEFAULT NULL REFERENCES labs(id)"),
        ]:
            if col not in existing_users:
                conn.execute(f"ALTER TABLE users ADD COLUMN {col} {ddl}")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_users_lab  ON users(lab_id)")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_orders_lab ON orders(lab_id)")

        existing_labs = {row[1] for row in conn.execute("PRAGMA table_info(labs)")}
        for col, ddl in [
            # The trusted, server-side link from a lab to its SFTP upload account
            # (the directory name under SFTP_ROOT/<account>/uploads on the portal
            # host). NULL until an admin sets it explicitly -- never inferred from
            # a lab name, a username, or anything the browser sends. A lab with no
            # sftp_account has no authorized SFTP uploads at all (see app.py).
            ("sftp_account", "TEXT DEFAULT NULL"),
        ]:
            if col not in existing_labs:
                conn.execute(f"ALTER TABLE labs ADD COLUMN {col} {ddl}")
        # SQLite UNIQUE indexes treat every NULL as distinct from every other NULL,
        # so any number of unlinked labs (sftp_account IS NULL) coexist fine --
        # only two labs sharing the same non-null account name is rejected, at
        # the database level, so it can never happen even under a race.
        conn.execute("CREATE UNIQUE INDEX IF NOT EXISTS idx_labs_sftp_account "
                    "ON labs(sftp_account)")


def _ensure_admin():
    from config import PORTAL_USER, PORTAL_PASS
    with _conn() as conn:
        exists = conn.execute(
            "SELECT 1 FROM users WHERE username=?", (PORTAL_USER,)
        ).fetchone()
        if not exists:
            conn.execute("""
                INSERT INTO users (username, password_hash, role, full_name, created_at)
                VALUES (?, ?, 'admin', 'Administrator', ?)
            """, (PORTAL_USER, generate_password_hash(PORTAL_PASS),
                  datetime.now().isoformat()))


# ── User management ────────────────────────────────────────────────────────────
def get_user(username: str):
    with _conn() as conn:
        return conn.execute(
            "SELECT * FROM users WHERE username=?", (username,)
        ).fetchone()


def verify_user(username: str, password: str):
    """Return user row if credentials valid and active, else None."""
    user = get_user(username)
    if user and user["active"] and check_password_hash(user["password_hash"], password):
        return user
    return None


def list_users() -> list:
    with _conn() as conn:
        return conn.execute(
            "SELECT * FROM users ORDER BY created_at"
        ).fetchall()


def create_user(username, password, role="lab_staff", full_name="", email="",
                lab_id=None) -> bool:
    try:
        with _conn() as conn:
            conn.execute("""
                INSERT INTO users (username, password_hash, role, full_name, email,
                                   created_at, lab_id)
                VALUES (?, ?, ?, ?, ?, ?, ?)
            """, (username, generate_password_hash(password), role,
                  full_name, email, datetime.now().isoformat(), lab_id))
        return True
    except sqlite3.IntegrityError:
        return False


def update_user(username, **kwargs):
    if "password" in kwargs:
        kwargs["password_hash"] = generate_password_hash(kwargs.pop("password"))
    if not kwargs:
        return
    set_clause = ", ".join(f"{k}=?" for k in kwargs)
    vals = list(kwargs.values()) + [username]
    with _conn() as conn:
        conn.execute(f"UPDATE users SET {set_clause} WHERE username=?", vals)


def delete_user(username: str):
    with _conn() as conn:
        conn.execute("DELETE FROM users WHERE username=?", (username,))


# ── Orders ─────────────────────────────────────────────────────────────────────
def new_order(patient_id, patient_name, sex, dob, tumor_type,
              panel_type, fastq_r1, fastq_r2,
              fastq_normal_r1="", fastq_normal_r2="",
              notes="", notify_email="", created_by="admin", lab_id=None,
              status="Queued") -> str:
    """status defaults to "Queued" (unchanged behaviour for every existing
    caller/test). app.py passes "Pinning" when the order has at least one
    SFTP-selected input still to be copied and verified -- see
    create_order_input() and runner.py's pinning worker."""
    order_id = "GR-" + datetime.now().strftime("%Y%m%d") + "-" + uuid.uuid4().hex[:6].upper()
    with _conn() as conn:
        conn.execute("""
            INSERT INTO orders
              (order_id, patient_id, patient_name, sex, dob, tumor_type,
               panel_type, status, created_at,
               fastq_r1, fastq_r2, fastq_normal_r1, fastq_normal_r2,
               notes, notify_email, created_by, lab_id)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """, (order_id, patient_id, patient_name, sex, dob, tumor_type,
              panel_type, status, datetime.now().isoformat(),
              fastq_r1, fastq_r2, fastq_normal_r1, fastq_normal_r2,
              notes, notify_email, created_by, lab_id))
        _audit(conn, order_id, "ORDER_CREATED", f"Panel: {panel_type}")
    return order_id


def new_order_with_inputs(patient_id, patient_name, sex, dob, tumor_type,
                          panel_type, inputs, notes="", notify_email="",
                          created_by="admin", lab_id=None) -> str:
    """Atomic variant of new_order() for the SFTP-pinning path: the order row
    (status='Pinning', fastq_* columns empty) and every selected order_inputs
    row are inserted in ONE transaction (one `with _conn() as conn:` block, so
    sqlite3 commits all of it together or rolls all of it back on any
    exception -- see _conn()/new_order() for the same pattern used elsewhere).

    This exists because order_inputs_all_done() treats an order with ZERO
    order_inputs rows as vacuously fully pinned (see its own docstring) --
    correct for an order that never selected any file, but wrong for one that
    did and would otherwise be left, after a crash or DB error between
    separate inserts, in status 'Pinning' with only SOME of its selected
    fields recorded. The pinning worker would then pin only those rows, see
    "all done", and route the order to 'Queued' to run with the missing
    field(s) still empty. Committing the order and all of its inputs in one
    transaction makes that partial state impossible: either every selected
    field is recorded alongside the order, or (on any failure) neither the
    order nor any of its inputs exist at all, and the caller/HTTP request can
    be retried cleanly from scratch.

    inputs: a non-empty list of (field_name, source_client, source_filename)
    triples -- use new_order() instead for an order with nothing to pin.
    """
    if not inputs:
        raise ValueError("new_order_with_inputs() requires at least one input; "
                         "use new_order() for an order with nothing to pin")
    order_id = "GR-" + datetime.now().strftime("%Y%m%d") + "-" + uuid.uuid4().hex[:6].upper()
    now = datetime.now().isoformat()
    with _conn() as conn:
        conn.execute("""
            INSERT INTO orders
              (order_id, patient_id, patient_name, sex, dob, tumor_type,
               panel_type, status, created_at,
               fastq_r1, fastq_r2, fastq_normal_r1, fastq_normal_r2,
               notes, notify_email, created_by, lab_id)
            VALUES (?,?,?,?,?,?,?,'Pinning',?,'','','','',?,?,?,?)
        """, (order_id, patient_id, patient_name, sex, dob, tumor_type,
              panel_type, now, notes, notify_email, created_by, lab_id))
        for field_name, source_client, source_filename in inputs:
            conn.execute("""
                INSERT INTO order_inputs (order_id, field_name, source_client, source_filename,
                                          status, created_at)
                VALUES (?,?,?,?, 'pending', ?)
            """, (order_id, field_name, source_client, source_filename, now))
        _audit(conn, order_id, "ORDER_CREATED", f"Panel: {panel_type}")
    return order_id


def get_order(order_id: str):
    with _conn() as conn:
        return conn.execute(
            "SELECT * FROM orders WHERE order_id=?", (order_id,)
        ).fetchone()


def list_orders(limit: int = 100, username: str = None, role: str = "admin", lab_id=None) -> list:
    """Admin sees all. Lab staff see their lab's orders (shared with every colleague
    linked to the same lab_id) plus, for backward compatibility, their own orders
    that predate lab-linking and so carry no lab_id at all -- never another lab's,
    and legacy no-lab access is never widened beyond the original creator."""
    with _conn() as conn:
        if role == "admin" or username is None:
            return conn.execute(
                "SELECT * FROM orders ORDER BY created_at DESC LIMIT ?", (limit,)
            ).fetchall()
        else:
            return conn.execute(
                "SELECT * FROM orders WHERE "
                "(lab_id IS NOT NULL AND lab_id=?) OR (lab_id IS NULL AND created_by=?) "
                "ORDER BY created_at DESC LIMIT ?",
                (lab_id, username, limit)
            ).fetchall()


def update_status(order_id: str, status: str, **kwargs):
    fields = {"status": status}
    if status == "Running":
        fields["started_at"] = datetime.now().isoformat()
    elif status in ("Done", "Failed"):
        fields["finished_at"] = datetime.now().isoformat()
    fields.update(kwargs)
    set_clause = ", ".join(f"{k}=?" for k in fields)
    vals = list(fields.values()) + [order_id]
    with _conn() as conn:
        conn.execute(f"UPDATE orders SET {set_clause} WHERE order_id=?", vals)
        _audit(conn, order_id, f"STATUS_{status.upper()}", str(kwargs))


# ── Order input pinning ──────────────────────────────────────────────────────
# See runner.py's "Input pinning" section for the durable worker that consumes
# these rows. Every write here that also needs to touch orders.<field_name> is
# done in a single connection/transaction so a reader never observes "input
# marked done" without the order's fastq_* column already pointing at it.
#
# field_name always comes from this fixed set in practice (app.py only ever
# passes one of these four literals; runner.py only ever reads field_name back
# out of rows this module itself wrote). The two functions below that splice
# it into an UPDATE ... SET {field_name}=? column name still assert against
# this whitelist explicitly rather than relying on that being true by
# construction -- an f-string column name is exactly the kind of thing a later
# edit could accidentally widen.
ORDER_FASTQ_FIELDS = ("fastq_r1", "fastq_r2", "fastq_normal_r1", "fastq_normal_r2")


def create_order_input(order_id: str, field_name: str, source_client: str, source_filename: str) -> None:
    with _conn() as conn:
        conn.execute("""
            INSERT INTO order_inputs (order_id, field_name, source_client, source_filename,
                                      status, created_at)
            VALUES (?,?,?,?, 'pending', ?)
        """, (order_id, field_name, source_client, source_filename, datetime.now().isoformat()))


def get_order_inputs(order_id: str) -> list:
    with _conn() as conn:
        return conn.execute(
            "SELECT * FROM order_inputs WHERE order_id=? ORDER BY field_name", (order_id,)
        ).fetchall()


def list_pinning_orders(limit: int = 50) -> list:
    """Orders currently in 'Pinning', oldest first -- what the pinning worker polls."""
    with _conn() as conn:
        return conn.execute(
            "SELECT * FROM orders WHERE status='Pinning' ORDER BY created_at ASC LIMIT ?", (limit,)
        ).fetchall()


def mark_input_copying(order_id: str, field_name: str) -> None:
    with _conn() as conn:
        conn.execute(
            "UPDATE order_inputs SET status='copying', error_msg='' WHERE order_id=? AND field_name=?",
            (order_id, field_name))


def mark_input_done(order_id: str, field_name: str, pinned_path: str, sha256: str, size_bytes: int) -> None:
    assert field_name in ORDER_FASTQ_FIELDS, field_name
    now = datetime.now().isoformat()
    with _conn() as conn:
        conn.execute("""
            UPDATE order_inputs SET status='done', pinned_path=?, sha256=?, size_bytes=?,
                   pinned_at=?, error_msg='' WHERE order_id=? AND field_name=?
        """, (pinned_path, sha256, size_bytes, now, order_id, field_name))
        conn.execute(f"UPDATE orders SET {field_name}=? WHERE order_id=?", (pinned_path, order_id))
        _audit(conn, order_id, "INPUT_PINNED", f"{field_name}: {size_bytes} bytes, sha256={sha256}")


def mark_input_failed(order_id: str, field_name: str, error_msg: str) -> None:
    with _conn() as conn:
        conn.execute(
            "UPDATE order_inputs SET status='failed', error_msg=? WHERE order_id=? AND field_name=?",
            (error_msg, order_id, field_name))
        _audit(conn, order_id, "INPUT_FAILED", f"{field_name}: {error_msg}")


def order_inputs_all_done(order_id: str) -> bool:
    """True if every recorded input for this order is pinned. Also true,
    vacuously, for an order with no order_inputs rows at all (it never had
    anything to pin -- e.g. no FASTQ was selected -- so there is nothing
    outstanding); retry_order() relies on that to send such an order straight
    back to Queued instead of routing it through Pinning for no reason."""
    return all(r["status"] == "done" for r in get_order_inputs(order_id))


def order_inputs_any_failed(order_id: str) -> bool:
    return any(r["status"] == "failed" for r in get_order_inputs(order_id))


def clear_order_input_pin(order_id: str, field_name: str) -> None:
    """Forget a successful pin's recorded path/hash (the caller deletes the
    actual file) -- used when a sibling field fails, since an order can only
    run with ALL of its selected inputs pinned, never some of them."""
    assert field_name in ORDER_FASTQ_FIELDS, field_name
    with _conn() as conn:
        conn.execute("""
            UPDATE order_inputs SET status='pending', pinned_path=NULL, sha256=NULL,
                   size_bytes=NULL, pinned_at=NULL WHERE order_id=? AND field_name=?
        """, (order_id, field_name))
        conn.execute(f"UPDATE orders SET {field_name}='' WHERE order_id=?", (order_id,))


def reset_order_inputs_for_retry(order_id: str) -> None:
    """Before routing a Failed order back through Pinning: any input not
    already 'done' is reset to 'pending' so the worker retries it (against its
    original source, which may no longer exist -- that failure, if it recurs,
    is correct and expected, not a bug)."""
    with _conn() as conn:
        conn.execute("""
            UPDATE order_inputs SET status='pending', error_msg=''
            WHERE order_id=? AND status!='done'
        """, (order_id,))


def stuck_copying_inputs() -> list:
    """order_inputs rows left in 'copying' -- i.e. the portal process died
    mid-copy. Used only at startup; see runner.py's _reset_stuck_pinning()."""
    with _conn() as conn:
        return conn.execute("SELECT * FROM order_inputs WHERE status='copying'").fetchall()


def reset_input_to_pending(order_id: str, field_name: str) -> None:
    with _conn() as conn:
        conn.execute(
            "UPDATE order_inputs SET status='pending' WHERE order_id=? AND field_name=?",
            (order_id, field_name))


def _audit(conn, order_id, event, detail=""):
    conn.execute(
        "INSERT INTO audit_log (order_id, event, detail, ts) VALUES (?,?,?,?)",
        (order_id, event, detail, datetime.now().isoformat())
    )


# ── Labs (shared billing profile) ─────────────────────────────────────────────
LAB_FIELDS = ("legal_name", "address_line1", "address_line2",
              "postal_code", "city", "country")
_LAB_REQUIRED = ("legal_name", "address_line1", "postal_code", "city", "country")
_LAB_MAXLEN = 200

# The SFTP upload account name is deliberately NOT in LAB_FIELDS: LAB_FIELDS is
# also what issue_invoice() snapshots onto a printed invoice, and this value must
# never appear there. It is a bare directory-name-equivalent token, never a path.
SFTP_ACCOUNT_RE = re.compile(r"^[A-Za-z0-9_.-]+$")


class LabError(ValueError):
    """Invalid lab data; the message is safe to show to the admin."""


def _clean_lab(fields: dict) -> dict:
    out = {}
    for k in LAB_FIELDS:
        v = " ".join(str(fields.get(k, "") or "").split())   # trim + collapse whitespace/newlines
        if len(v) > _LAB_MAXLEN:
            raise LabError(f"{k} is longer than {_LAB_MAXLEN} characters.")
        out[k] = v
    missing = [k for k in _LAB_REQUIRED if not out[k]]
    if missing:
        raise LabError("Required: " + ", ".join(missing) + ".")
    return out


def create_lab(**fields) -> int:
    f = _clean_lab(fields)
    now = datetime.now().isoformat()
    try:
        with _conn() as conn:
            cur = conn.execute("""
                INSERT INTO labs (legal_name, address_line1, address_line2, postal_code,
                                  city, country, created_at, updated_at)
                VALUES (?,?,?,?,?,?,?,?)
            """, (f["legal_name"], f["address_line1"], f["address_line2"],
                  f["postal_code"], f["city"], f["country"], now, now))
            return cur.lastrowid
    except sqlite3.IntegrityError:
        raise LabError("A lab with this legal name already exists.")


def update_lab(lab_id: int, **fields) -> None:
    f = _clean_lab(fields)
    try:
        with _conn() as conn:
            cur = conn.execute("""
                UPDATE labs SET legal_name=?, address_line1=?, address_line2=?,
                       postal_code=?, city=?, country=?, updated_at=?
                WHERE id=?
            """, (f["legal_name"], f["address_line1"], f["address_line2"],
                  f["postal_code"], f["city"], f["country"],
                  datetime.now().isoformat(), lab_id))
            if cur.rowcount == 0:
                raise LabError("Lab not found.")
    except sqlite3.IntegrityError:
        raise LabError("A lab with this legal name already exists.")


def set_lab_sftp_account(lab_id: int, account) -> None:
    """Set (or clear, with None/'') the lab's trusted SFTP upload account name.

    This is the ONLY server-side source of truth for 'which SFTP client
    directory belongs to which lab' -- portal/app.py never infers it from a
    lab name, a username, or anything a request supplies. Admin-only at the
    route level. Rejects anything but a bare token (letters/digits/._-): no
    slashes, no '..', so it can never be used to escape SFTP_ROOT.
    """
    account = (account or "").strip()
    if account and (account in (".", "..") or not SFTP_ACCOUNT_RE.match(account)):
        raise LabError("SFTP account must contain only letters, digits, '.', '_' or '-' "
                       "(and cannot be '.' or '..').")
    try:
        with _conn() as conn:
            cur = conn.execute("UPDATE labs SET sftp_account=?, updated_at=? WHERE id=?",
                               (account or None, datetime.now().isoformat(), lab_id))
            if cur.rowcount == 0:
                raise LabError("Lab not found.")
    except sqlite3.IntegrityError:
        # idx_labs_sftp_account (UNIQUE) rejected this write atomically -- the
        # update never took effect, so no lab is left with a duplicate link.
        raise LabError(f"SFTP account '{account}' is already linked to another lab.")


def get_lab(lab_id):
    if lab_id is None:
        return None
    with _conn() as conn:
        return conn.execute("SELECT * FROM labs WHERE id=?", (lab_id,)).fetchone()


def list_labs() -> list:
    with _conn() as conn:
        return conn.execute("""
            SELECT l.*, (SELECT COUNT(*) FROM users u WHERE u.lab_id = l.id) AS n_users
            FROM labs l ORDER BY l.legal_name COLLATE NOCASE
        """).fetchall()


def lab_members(lab_id: int) -> list:
    with _conn() as conn:
        return conn.execute(
            "SELECT username, full_name, role, active FROM users WHERE lab_id=? ORDER BY username",
            (lab_id,)).fetchall()


def set_user_lab(username: str, lab_id) -> None:
    """Link a user to a lab (or unlink with None). Admin-only at the route level."""
    if lab_id is not None and get_lab(lab_id) is None:
        raise LabError("Lab not found.")
    with _conn() as conn:
        cur = conn.execute("UPDATE users SET lab_id=? WHERE username=?", (lab_id, username))
        if cur.rowcount == 0:
            raise LabError("User not found.")


def set_order_lab(order_id: str, lab_id) -> None:
    """Attach a lab to an order. An order that already has an issued invoice is frozen."""
    if get_lab(lab_id) is None:
        raise LabError("Lab not found.")
    with _conn() as conn:
        if conn.execute("SELECT 1 FROM invoices WHERE order_id=?", (order_id,)).fetchone():
            raise LabError("This order already has an issued invoice.")
        cur = conn.execute("UPDATE orders SET lab_id=? WHERE order_id=?", (lab_id, order_id))
        if cur.rowcount == 0:
            raise LabError("Order not found.")
        _audit(conn, order_id, "ORDER_LAB_SET", f"lab_id={lab_id}")


# ── Invoices ──────────────────────────────────────────────────────────────────
def get_invoice(order_id: str):
    with _conn() as conn:
        return conn.execute("SELECT * FROM invoices WHERE order_id=?", (order_id,)).fetchone()


def issue_invoice(order_id: str, *, issued_by: str, invoice_date: str, service_date: str,
                  description: str, net_cents: int, issuer: dict):
    """Issue exactly one invoice for an order with a gapless per-year number.

    Called only from an explicit POST. The number, recipient snapshot and amounts
    are written in one immediate transaction; a second call for the same order
    returns the existing invoice unchanged.
    """
    import json
    conn = _conn()
    conn.isolation_level = None            # manual transaction control
    try:
        conn.execute("BEGIN IMMEDIATE")
        existing = conn.execute("SELECT * FROM invoices WHERE order_id=?", (order_id,)).fetchone()
        if existing:
            conn.execute("COMMIT")
            return existing
        order = conn.execute("SELECT lab_id FROM orders WHERE order_id=?", (order_id,)).fetchone()
        if not order:
            raise LabError("Order not found.")
        lab = None
        if order["lab_id"] is not None:
            lab = conn.execute("SELECT * FROM labs WHERE id=?", (order["lab_id"],)).fetchone()
        if not lab:
            raise LabError("The order has no billing lab assigned.")
        year = int(invoice_date[:4])
        row = conn.execute("SELECT last_no FROM invoice_counters WHERE year=?", (year,)).fetchone()
        n = (row["last_no"] if row else 0) + 1
        if row:
            conn.execute("UPDATE invoice_counters SET last_no=? WHERE year=?", (n, year))
        else:
            conn.execute("INSERT INTO invoice_counters (year, last_no) VALUES (?,?)", (year, n))
        invoice_no = f"GR-{year}-{n:04d}"
        recipient = {k: lab[k] for k in LAB_FIELDS}
        now = datetime.now().isoformat()
        conn.execute("""
            INSERT INTO invoices (invoice_no, order_id, lab_id, invoice_date, service_date,
                                  description, net_cents, recipient_json, issuer_json,
                                  issued_by, created_at)
            VALUES (?,?,?,?,?,?,?,?,?,?,?)
        """, (invoice_no, order_id, lab["id"], invoice_date, service_date, description,
              net_cents, json.dumps(recipient, ensure_ascii=False, sort_keys=True),
              json.dumps(issuer, ensure_ascii=False, sort_keys=True), issued_by, now))
        _audit(conn, order_id, "INVOICE_ISSUED", invoice_no)
        inv = conn.execute("SELECT * FROM invoices WHERE order_id=?", (order_id,)).fetchone()
        conn.execute("COMMIT")
        return inv
    except BaseException:
        if conn.in_transaction:
            conn.execute("ROLLBACK")
        raise
    finally:
        conn.close()
