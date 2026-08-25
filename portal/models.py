"""GenRichi Portal — SQLite database models"""

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

        CREATE TABLE IF NOT EXISTS audit_log (
            id         INTEGER PRIMARY KEY AUTOINCREMENT,
            order_id   TEXT,
            event      TEXT,
            detail     TEXT,
            ts         TEXT
        );
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
        ]:
            if col not in existing_orders:
                conn.execute(f"ALTER TABLE orders ADD COLUMN {col} {ddl}")

        existing_users = {row[1] for row in conn.execute("PRAGMA table_info(users)")}
        for col, ddl in [
            ("full_name", "TEXT DEFAULT ''"),
            ("email",     "TEXT DEFAULT ''"),
            ("active",    "INTEGER DEFAULT 1"),
        ]:
            if col not in existing_users:
                conn.execute(f"ALTER TABLE users ADD COLUMN {col} {ddl}")


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


def create_user(username, password, role="lab_staff", full_name="", email="") -> bool:
    try:
        with _conn() as conn:
            conn.execute("""
                INSERT INTO users (username, password_hash, role, full_name, email, created_at)
                VALUES (?, ?, ?, ?, ?, ?)
            """, (username, generate_password_hash(password), role,
                  full_name, email, datetime.now().isoformat()))
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
              notes="", notify_email="", created_by="admin") -> str:
    order_id = "GR-" + datetime.now().strftime("%Y%m%d") + "-" + uuid.uuid4().hex[:6].upper()
    with _conn() as conn:
        conn.execute("""
            INSERT INTO orders
              (order_id, patient_id, patient_name, sex, dob, tumor_type,
               panel_type, status, created_at,
               fastq_r1, fastq_r2, fastq_normal_r1, fastq_normal_r2,
               notes, notify_email, created_by)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """, (order_id, patient_id, patient_name, sex, dob, tumor_type,
              panel_type, "Queued", datetime.now().isoformat(),
              fastq_r1, fastq_r2, fastq_normal_r1, fastq_normal_r2,
              notes, notify_email, created_by))
        _audit(conn, order_id, "ORDER_CREATED", f"Panel: {panel_type}")
    return order_id


def get_order(order_id: str):
    with _conn() as conn:
        return conn.execute(
            "SELECT * FROM orders WHERE order_id=?", (order_id,)
        ).fetchone()


def list_orders(limit: int = 100, username: str = None, role: str = "admin") -> list:
    """Admin sees all; lab staff see only their own."""
    with _conn() as conn:
        if role == "admin" or username is None:
            return conn.execute(
                "SELECT * FROM orders ORDER BY created_at DESC LIMIT ?", (limit,)
            ).fetchall()
        else:
            return conn.execute(
                "SELECT * FROM orders WHERE created_by=? ORDER BY created_at DESC LIMIT ?",
                (username, limit)
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


def _audit(conn, order_id, event, detail=""):
    conn.execute(
        "INSERT INTO audit_log (order_id, event, detail, ts) VALUES (?,?,?,?)",
        (order_id, event, detail, datetime.now().isoformat())
    )
