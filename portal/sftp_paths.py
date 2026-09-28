"""GenRichi Portal — SFTP upload path safety.

Shared by app.py (selection-time: browsing and resolving a form token when an
order is created) and runner.py (pin-time: re-verifying the exact same source
right before it is copied). Kept in its own module, not in app.py, so runner.py
can import it without a circular import (app.py imports runner.py at startup).

The trusted server-side link from a lab to its SFTP directory is
labs.sftp_account (models.set_lab_sftp_account), set by an admin -- never
inferred here from a lab name, a username, or anything a request sends.
"""
import os
from pathlib import Path

from models import SFTP_ACCOUNT_RE

SFTP_ROOT = "/srv/genrichi-sftp"


def sftp_upload_root(client):
    """The real, existing uploads directory for exactly one SFTP client account,
    or None. Refuses anything that isn't a bare account-name token (no '/', and
    '.'/'..' rejected explicitly even though SFTP_ACCOUNT_RE would otherwise
    accept them as "all allowed characters"). Refuses if EITHER the account
    directory or its 'uploads' subdirectory is a symlink -- even one that still
    resolves to somewhere under SFTP_ROOT -- which would otherwise let one
    account's directory silently alias another account's real files. islink()
    is an lstat, so it never follows the link it is inspecting."""
    if not client or client in (".", "..") or not SFTP_ACCOUNT_RE.match(client):
        return None
    base = os.path.realpath(SFTP_ROOT)
    account_dir = os.path.join(base, client)
    uploads_dir = os.path.join(account_dir, "uploads")
    try:
        if os.path.islink(account_dir) or os.path.islink(uploads_dir):
            return None
        if not os.path.isdir(uploads_dir):
            return None
    except OSError:
        return None
    return uploads_dir


def list_sftp_clients():
    """Every real client-account directory name directly under SFTP_ROOT
    (admin-only listing). A plain wrapper kept in this module -- not
    reimplemented against a bare SFTP_ROOT reference elsewhere -- purely so
    that patching SFTP_ROOT here (e.g. in tests) affects this too."""
    try:
        return sorted(d.name for d in Path(SFTP_ROOT).iterdir() if d.is_dir())
    except OSError:
        return []


def list_sftp_files(client):
    root = sftp_upload_root(client)
    if root is None:
        return []
    out = []
    try:
        for f in sorted(Path(root).iterdir()):
            if f.is_symlink():
                continue  # never list a link, even one pointing to a real fastq file
            if f.is_file() and (f.suffix.lower() in (".gz", ".fastq", ".fq") or
                                f.name.endswith(".fastq.gz") or f.name.endswith(".fq.gz")):
                out.append({"client": client, "name": f.name, "token": f"{client}/{f.name}",
                            "size_mb": round(f.stat().st_size / 1024 / 1024, 1)})
    except OSError:
        pass
    return out


def resolve_sftp_token(token, allowed_client):
    """Resolve a 'client/filename' token to a real file path.

    allowed_client=None means "any real client directory" (admin only, at
    selection time). A non-None allowed_client restricts the token's client
    segment to exactly that value. Returns None on any mismatch, missing link,
    traversal attempt, symlink (account, uploads dir, or the file itself), or
    missing file; never raises, never partially trusts the input.

    Used twice per file, independently: once when a lab user selects it (app.py,
    allowed_client = their lab's own account), and again right before it is
    copied into the pinned-inputs area (runner.py, allowed_client = the exact
    client recorded for that input at selection time) -- the second call is a
    fresh re-derivation, not a reuse of any value cached from the first.
    """
    if not token or "/" not in token:
        return None
    client, _, name = token.partition("/")
    if not name or name in (".", "..") or "/" in name or "\\" in name:
        return None
    if allowed_client is not None and client != allowed_client:
        return None
    root = sftp_upload_root(client)  # already confirmed free of symlinks at both levels
    if root is None:
        return None
    candidate = os.path.join(root, name)
    try:
        if os.path.islink(candidate):
            return None  # never dereference a file-level symlink either
        real = os.path.realpath(candidate)
    except OSError:
        return None
    if not real.startswith(root + os.sep) or not os.path.isfile(real):
        return None
    return real


def path_is_within_sftp_root(path) -> bool:
    """True if `path` is inside SFTP_ROOT -- the client-writable upload tree
    -- either LEXICALLY (its own location, as a plain string, before
    following any symlink anywhere in it) or by where it actually RESOLVES
    to. Either one is enough to call it "inside".

    Checking only the resolved target (as an earlier version of this
    function did) is not sufficient: a path that lexically lives inside
    SFTP_ROOT is, by definition, something the SFTP client can write to --
    including replacing it with, or repointing, a symlink -- even if today
    it happens to resolve to somewhere entirely outside the root. The
    client can change that target at any later moment (classic TOCTOU), so a
    location the client controls is untrusted regardless of what it
    currently points at. Conversely, a path that lexically sits *outside*
    SFTP_ROOT but resolves (via a symlink somewhere in it) *into* SFTP_ROOT
    is exactly as unsafe as one recorded directly -- that is what the
    resolved-target check still catches.

    Used by runner.py to catch a pre-pinning-era order whose orders.fastq_*
    column still points at (or through) a live SFTP path (no order_inputs
    row was ever created for it, so nothing re-verifies it): such a value
    must never be treated as safe just because it happens to be non-empty.
    Only compares locations; it does not require the path to currently exist
    -- a path that no longer exists is handled by the caller, not here."""
    if not path:
        return False
    try:
        lexical_root = os.path.abspath(SFTP_ROOT)
        lexical_path = os.path.abspath(path)
    except OSError:
        return False
    if lexical_path == lexical_root or lexical_path.startswith(lexical_root + os.sep):
        return True
    try:
        real_root = os.path.realpath(SFTP_ROOT)
        real_path = os.path.realpath(path)
    except OSError:
        return False
    return real_path == real_root or real_path.startswith(real_root + os.sep)
