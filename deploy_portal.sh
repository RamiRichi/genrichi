#!/usr/bin/env bash
# GenRichi portal deployment -- DRY-RUN by default.
#
# Copies ONLY the explicit allowlist of portal application files (code,
# templates, static image) from a Git commit into /home/rami/genrichi/portal.
# It never restarts the running portal and never changes any systemd unit.
#
# Usage:  bash deploy_portal.sh [options]
#   (no option)         dry-run: run every check, print the plan, write nothing
#   --execute           perform the deployment (backup first, atomic writes)
#   --from-worktree     deploy files from the working tree instead of HEAD
#   --allow-dirty       with --from-worktree: allow uncommitted changes under portal/
#   --expect-commit SHA refuse unless HEAD is exactly this full 40-hex commit SHA
#                       (in --from-worktree mode it pins HEAD only, not the file bytes)
#   --help              show this text
#
# Exit code: 0 ok, 1 refused / a check failed, 2 usage error.
# See scripts/portal_status.sh for read-only verification after a deployment.

set -euo pipefail
umask 022

# ── fixed production locations ───────────────────────────────────────────────
readonly PROD_HOME="/home/rami/genrichi"
readonly PROD_REPO="/mnt/c/GenRichi"
readonly PROD_PYTHON="/home/rami/miniforge3/envs/snakemake/bin/python"
readonly PORTAL_UNIT="genrichi-portal.service"
readonly DUPLICATE_UNIT="genrichi.service"
readonly EXPECTED_REMOTE_MATCH="RamiRichi/genrichi"
readonly KNOWN_GOOD_FLASK="3.1.3"
readonly KNOWN_GOOD_WERKZEUG="3.1.8"

# The ONLY files this script may write (relative to portal/).
readonly ALLOWLIST=(
    app.py
    config.py
    mailer.py
    models.py
    runner.py
    templates/base.html
    templates/dashboard.html
    templates/invoice.html
    templates/login.html
    templates/new_order.html
    templates/order.html
    templates/report_view.html
    templates/settings.html
    templates/stats.html
    templates/users.html
    static/img/logo.png
)
# Tracked under portal/ but deliberately not deployed.
readonly NOT_DEPLOYED_TRACKED=(.env.example)

# ── test mode (unit tests only; every path must live under /tmp) ─────────────
TEST_MODE=0
HOME_DIR="$PROD_HOME"
REPO="$PROD_REPO"
PY="$PROD_PYTHON"
CGROUP_ROOT="/sys/fs/cgroup"
REQUIRED_MODULES="flask werkzeug"
if [ -n "${GENRICHI_DEPLOY_TEST_HOME:-}" ]; then
    TEST_MODE=1
    HOME_DIR="$GENRICHI_DEPLOY_TEST_HOME"
    REPO="${GENRICHI_DEPLOY_TEST_REPO:?test mode needs GENRICHI_DEPLOY_TEST_REPO}"
    PY="${GENRICHI_DEPLOY_TEST_PYTHON:?test mode needs GENRICHI_DEPLOY_TEST_PYTHON}"
    CGROUP_ROOT="${GENRICHI_DEPLOY_TEST_CGROUP_ROOT:?test mode needs GENRICHI_DEPLOY_TEST_CGROUP_ROOT}"
    REQUIRED_MODULES="${GENRICHI_DEPLOY_TEST_MODULES:-flask werkzeug}"
    for p in "$HOME_DIR" "$REPO" "$CGROUP_ROOT"; do
        case "$p" in /tmp/*) ;; *) echo "test mode only accepts paths under /tmp: refusing" >&2; exit 2 ;; esac
    done
fi
DST="$HOME_DIR/portal"
BACKUP_PARENT="$HOME_DIR/deploy_backups/portal"
ENV_FILE="$DST/.env"
LOCK_FILE="${TMPDIR:-/tmp}/genrichi-portal-deploy.lock"

# ── output helpers ───────────────────────────────────────────────────────────
FAILURES=()
WARNINGS=()
section() { printf '\n%s\n' "$*"; }
info()    { printf '  %s\n' "$*"; }
pass()    { printf '  ok    %s\n' "$*"; }
fail()    { FAILURES+=("$*"); printf '  FAIL  %s\n' "$*"; }
warn()    { WARNINGS+=("$*"); printf '  warn  %s\n' "$*"; }
die()     { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

usage() {
    cat <<'USAGE'
GenRichi portal deployment (dry-run unless --execute is given)

  bash deploy_portal.sh [--execute] [--from-worktree [--allow-dirty]]
                        [--expect-commit FULL_40_HEX_SHA] [--help]

Deploys only the fixed allowlist of portal files. It never writes .env, the
database, uploads, results or Snakemake data, never installs packages, never
restarts the portal and never changes a systemd unit.
USAGE
}

# ── arguments ────────────────────────────────────────────────────────────────
EXECUTE=0
FROM_WORKTREE=0
ALLOW_DIRTY=0
EXPECT_COMMIT=""
while [ $# -gt 0 ]; do
    case "$1" in
        --execute)       EXECUTE=1 ;;
        --from-worktree) FROM_WORKTREE=1 ;;
        --allow-dirty)   ALLOW_DIRTY=1 ;;
        --expect-commit)
            [ $# -ge 2 ] && [ -n "$2" ] || { echo "--expect-commit needs a value" >&2; exit 2; }
            EXPECT_COMMIT="$2"; shift ;;
        -h|--help)       usage; exit 0 ;;
        *)               echo "unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
    shift
done
if [ -n "$EXPECT_COMMIT" ]; then
    case "$EXPECT_COMMIT" in *[!0-9a-fA-F]*) echo "--expect-commit needs a hex SHA" >&2; exit 2 ;; esac
    if [ "${#EXPECT_COMMIT}" -ne 40 ]; then
        echo "--expect-commit needs the full 40-character Git SHA (got ${#EXPECT_COMMIT}); short prefixes are not accepted" >&2
        exit 2
    fi
    EXPECT_COMMIT="${EXPECT_COMMIT,,}"
fi

# ── small utilities ──────────────────────────────────────────────────────────
sha_stdin() { sha256sum | cut -d' ' -f1; }

path_has_symlink() {
    local p="$1" cur="" comp
    local IFS=/
    for comp in $p; do
        [ -z "$comp" ] && continue
        cur="$cur/$comp"
        if [ -L "$cur" ]; then return 0; fi
    done
    return 1
}

in_list() {  # in_list VALUE ITEM...
    local needle="$1" item
    shift
    for item in "$@"; do [ "$item" = "$needle" ] && return 0; done
    return 1
}

# --no-optional-locks: read-only queries must never refresh (rewrite) .git/index.
repo_git() { git --no-optional-locks -c core.fileMode=false -C "$REPO" "$@"; }

src_cat() {  # relative path under portal/ -> file bytes on stdout
    if [ "$FROM_WORKTREE" -eq 1 ]; then
        cat -- "$REPO/portal/$1"
    else
        repo_git cat-file blob "$COMMIT:portal/$1"
    fi
}

# Lines (numbers only) that assign a quoted literal to a secret-named variable.
CRED_PATTERN='(pass|passwd|password|secret|token|api_?key)[a-z0-9_]*[[:space:]]*=[[:space:]]*["'"'"'][^"'"'"']{6,}["'"'"']'
CRED_ALLOWED='os[.]environ|getenv|_require_env|request[.]form|placeholder=|<input|type='
cred_lines() {
    src_cat "$1" | CRED_PATTERN="$CRED_PATTERN" CRED_ALLOWED="$CRED_ALLOWED" \
        awk '{ l = tolower($0) } l ~ ENVIRON["CRED_PATTERN"] && l !~ ENVIRON["CRED_ALLOWED"] { print NR }'
}

# ── banner ───────────────────────────────────────────────────────────────────
if [ "$EXECUTE" -eq 1 ]; then MODE_LABEL="EXECUTE"; else MODE_LABEL="DRY-RUN (nothing will be written)"; fi
echo "GenRichi portal deployment -- $MODE_LABEL"
if [ "$TEST_MODE" -eq 1 ]; then echo "  *** TEST MODE: paths are redirected under /tmp ***"; fi
info "destination : $DST"
info "source mode : $([ "$FROM_WORKTREE" -eq 1 ] && echo 'working tree (opt-in)' || echo 'Git HEAD commit (default)')"
if [ "$ALLOW_DIRTY" -eq 1 ] && [ "$FROM_WORKTREE" -eq 0 ]; then warn "--allow-dirty has no effect without --from-worktree"; fi

# ── 1. source ────────────────────────────────────────────────────────────────
section "1. Source"
COMMIT=""
if ! command -v git >/dev/null 2>&1; then
    fail "git is not available"
elif [ ! -d "$REPO" ]; then
    fail "source repository not found: $REPO"
else
    TOP="$(repo_git rev-parse --show-toplevel 2>/dev/null || true)"
    if [ "$TEST_MODE" -eq 0 ] && [ "$TOP" != "$PROD_REPO" ]; then
        fail "repository toplevel is '${TOP:-unknown}', expected $PROD_REPO"
    elif [ -z "$TOP" ]; then
        fail "$REPO is not a Git repository"
    else
        COMMIT="$(repo_git rev-parse --verify 'HEAD^{commit}' 2>/dev/null || true)"
        if [ -z "$COMMIT" ]; then
            fail "cannot resolve HEAD in $REPO"
        else
            if [ "$FROM_WORKTREE" -eq 1 ]; then
                pass "HEAD is $COMMIT ($(repo_git log -1 --format='%cs' "$COMMIT" 2>/dev/null)); file bytes come from the WORKING TREE, not from this commit"
            else
                pass "commit $COMMIT ($(repo_git log -1 --format='%cs' "$COMMIT" 2>/dev/null)); file bytes are read from this commit"
            fi
            info "subject: $(repo_git log -1 --format='%s' "$COMMIT" 2>/dev/null | cut -c1-90)"
            if [ -n "$EXPECT_COMMIT" ]; then
                if [ "$COMMIT" = "$EXPECT_COMMIT" ]; then
                    pass "HEAD matches --expect-commit $EXPECT_COMMIT"
                else
                    fail "HEAD $COMMIT does not match --expect-commit $EXPECT_COMMIT"
                fi
                if [ "$FROM_WORKTREE" -eq 1 ]; then
                    warn "--expect-commit pins HEAD only: the deployed bytes come from the working tree and are not pinned to any commit"
                fi
            else
                warn "no --expect-commit given: pass the reviewed 40-hex SHA to pin the source"
            fi
            if [ "$TEST_MODE" -eq 0 ]; then
                REMOTE_URL="$(repo_git config --get remote.origin.url 2>/dev/null || true)"
                case "$REMOTE_URL" in
                    *"$EXPECTED_REMOTE_MATCH"*) pass "origin remote is the expected repository" ;;
                    *)                          fail "origin remote does not look like $EXPECTED_REMOTE_MATCH" ;;
                esac
            fi
            if [ -z "$(repo_git branch -r --contains "$COMMIT" 2>/dev/null)" ]; then
                warn "this commit is not on any remote branch (unpushed)"
            fi
        fi
    fi
fi

DIRTY_SCOPE=""
if [ -n "$COMMIT" ]; then
    DIRTY_SCOPE="$(repo_git status --porcelain --untracked-files=normal -- portal 2>/dev/null || true)"
    if [ "$FROM_WORKTREE" -eq 1 ]; then
        if [ -n "$DIRTY_SCOPE" ] && [ "$ALLOW_DIRTY" -eq 0 ]; then
            fail "uncommitted changes under portal/ and --allow-dirty was not given ($(printf '%s\n' "$DIRTY_SCOPE" | wc -l) entries)"
        elif [ -n "$DIRTY_SCOPE" ]; then
            warn "deploying uncommitted working-tree changes (--allow-dirty)"
        else
            pass "portal/ working tree is clean"
        fi
    elif [ -n "$DIRTY_SCOPE" ]; then
        info "portal/ has uncommitted changes; they are IGNORED because HEAD blobs are deployed"
    else
        pass "portal/ working tree is clean"
    fi
fi

# tracked files under portal/ must all be classified (fail closed on unknown files)
SOURCE_OK=0
if [ -n "$COMMIT" ]; then
    if [ "$FROM_WORKTREE" -eq 1 ]; then
        TRACKED="$(repo_git ls-files -- portal 2>/dev/null || true)"
    else
        TRACKED="$(repo_git ls-tree -r --name-only "$COMMIT" -- portal 2>/dev/null || true)"
    fi
    unknown=0
    while IFS= read -r f; do
        [ -n "$f" ] || continue
        rel="${f#portal/}"
        if in_list "$rel" "${ALLOWLIST[@]}" || in_list "$rel" "${NOT_DEPLOYED_TRACKED[@]}"; then continue; fi
        fail "unknown tracked file under portal/: $rel (add it to the allowlist deliberately or remove it)"
        unknown=$((unknown + 1))
    done <<< "$TRACKED"
    missing=0
    for rel in "${ALLOWLIST[@]}"; do
        if [ "$FROM_WORKTREE" -eq 1 ]; then
            [ -f "$REPO/portal/$rel" ] && [ ! -L "$REPO/portal/$rel" ] || { fail "allowlisted file missing in the working tree: $rel"; missing=$((missing + 1)); }
        else
            repo_git cat-file -e "$COMMIT:portal/$rel" 2>/dev/null || { fail "allowlisted file missing in commit: $rel"; missing=$((missing + 1)); }
        fi
    done
    if [ "$unknown" -eq 0 ] && [ "$missing" -eq 0 ]; then
        pass "all ${#ALLOWLIST[@]} allowlisted files present; no unknown tracked files under portal/"
        SOURCE_OK=1
    fi
fi

# content checks on the exact bytes that would be deployed
declare -A NEW_SHA
if [ "$SOURCE_OK" -eq 1 ]; then
    bad_content=0
    for rel in "${ALLOWLIST[@]}"; do
        NEW_SHA["$rel"]="$(src_cat "$rel" | sha_stdin)"
        case "$rel" in
            *.py|*.html)
                crs="$(src_cat "$rel" | tr -cd '\r' | wc -c)"
                if [ "$crs" -ne 0 ]; then fail "carriage returns in $rel (CRLF would corrupt the deployed file)"; bad_content=1; fi
                hits="$(cred_lines "$rel" | tr '\n' ' ')"
                if [ -n "$hits" ]; then fail "credential-like literal in $rel at line(s): $hits(value not shown)"; bad_content=1; fi
                ;;
        esac
        case "$rel" in
            *.py)
                if ! src_cat "$rel" | "$PY" -B -c 'import ast, sys; ast.parse(sys.stdin.buffer.read())' 2>/dev/null; then
                    fail "Python syntax check failed for $rel"; bad_content=1
                fi ;;
        esac
    done
    if [ "$bad_content" -eq 0 ]; then pass "source bytes: LF only, no credential-like literals, Python files parse"; fi
fi

# ── 2. destination ───────────────────────────────────────────────────────────
section "2. Destination and production guards"
if [ "$TEST_MODE" -eq 0 ] && [ "$DST" != "/home/rami/genrichi/portal" ]; then
    fail "destination is not exactly /home/rami/genrichi/portal"
elif [ ! -d "$DST" ]; then
    fail "destination directory does not exist: $DST"
elif path_has_symlink "$DST"; then
    fail "a component of the destination path is a symlink"
else
    pass "destination is the expected directory and contains no symlink component"
fi
if [ "$(id -u)" -eq 0 ]; then
    fail "refusing to run as root (deployed files must stay owned by the service user)"
elif [ -d "$DST" ] && [ "$(stat -c '%U' "$DST")" != "$(id -un)" ]; then
    fail "destination is owned by $(stat -c '%U' "$DST"), not by $(id -un)"
else
    pass "running as $(id -un), the owner of the destination"
fi
for sub in templates static static/img; do
    if [ -L "$DST/$sub" ]; then fail "destination sub-directory is a symlink: $sub"
    elif [ -e "$DST/$sub" ] && [ ! -d "$DST/$sub" ]; then fail "destination path is not a directory: $sub"; fi
done
[ -d "$DST" ] && [ ! -w "$DST" ] && fail "destination is not writable"

if [ -L "$ENV_FILE" ] || [ ! -f "$ENV_FILE" ]; then
    fail "production .env is missing (or not a regular file); it is never created or copied by this script"
else
    env_mode="$(stat -c '%a' "$ENV_FILE")"
    if [ $(( 8#$env_mode & 8#7137 )) -ne 0 ]; then
        fail "production .env mode is $env_mode; it must not be wider than 640"
    else
        pass "production .env exists with mode $env_mode (contents never read)"
    fi
fi

# ── 3. services and running work (read-only queries only) ────────────────────
section "3. Services and running analyses"
if ! command -v systemctl >/dev/null 2>&1; then
    fail "systemctl is not available: cannot verify the service state"
else
    dup_state="$(systemctl is-active "$DUPLICATE_UNIT" 2>/dev/null || true)"
    case "$dup_state" in
        active|activating|reloading|deactivating) fail "$DUPLICATE_UNIT is $dup_state; it must not run" ;;
        inactive|failed)                          pass "$DUPLICATE_UNIT is $dup_state" ;;
        *)                                        fail "cannot determine the state of $DUPLICATE_UNIT (got '${dup_state:-nothing}')" ;;
    esac
    case "$(systemctl is-enabled "$DUPLICATE_UNIT" 2>/dev/null || true)" in
        enabled|enabled-runtime|static|alias) warn "$DUPLICATE_UNIT is enabled; it should stay disabled" ;;
    esac

    portal_state="$(systemctl is-active "$PORTAL_UNIT" 2>/dev/null || true)"
    info "$PORTAL_UNIT is ${portal_state:-unknown} (never restarted by this script)"
    if [ "$portal_state" = "active" ]; then
        cg="$(systemctl show "$PORTAL_UNIT" -p ControlGroup --value 2>/dev/null || true)"
        procs_file="$CGROUP_ROOT$cg/cgroup.procs"
        if [ -z "$cg" ] || [ ! -r "$procs_file" ]; then
            fail "cannot read the cgroup of $PORTAL_UNIT: unable to prove that no analysis is running"
        else
            nprocs="$(wc -l < "$procs_file")"
            if [ "$nprocs" -gt 1 ]; then fail "an analysis appears to be running ($nprocs processes in the portal cgroup)"
            else pass "the portal cgroup holds only the portal itself (no pipeline child)"; fi
        fi
    elif [ "$portal_state" != "inactive" ] && [ "$portal_state" != "failed" ]; then
        fail "cannot determine the state of $PORTAL_UNIT (got '${portal_state:-nothing}')"
    fi
fi
# pgrep exits 0 (match) or 1 (no match, prints 0 with -c); anything else is an error.
# An error or unusable output must never be read as "nothing is running".
if ! command -v pgrep >/dev/null 2>&1; then
    fail "pgrep is not available: unable to prove that no Snakemake process is running"
else
    pgrep_rc=0
    snake_count="$(pgrep -f -c '[s]nakemake --snakefile' 2>/dev/null)" || pgrep_rc=$?
    if [ "$pgrep_rc" -gt 1 ]; then
        fail "pgrep failed (exit $pgrep_rc): unable to prove that no Snakemake process is running"
    elif ! [[ "$snake_count" =~ ^[0-9]+$ ]] || { [ "$pgrep_rc" -eq 1 ] && [ "$snake_count" -ne 0 ]; } || { [ "$pgrep_rc" -eq 0 ] && [ "$snake_count" -eq 0 ]; }; then
        fail "pgrep gave unusable output: unable to prove that no Snakemake process is running"
    elif [ "$snake_count" -ne 0 ]; then
        fail "$snake_count Snakemake process(es) are running"
    else
        pass "no Snakemake process is running"
    fi
fi
if [ -e "$HOME_DIR/.snakemake/locks" ]; then fail "$HOME_DIR/.snakemake/locks exists (a run may hold the directory)"; else pass "no .snakemake/locks directory"; fi

# ── 4. dependencies (verified, never installed) ──────────────────────────────
section "4. Python dependencies (verification only)"
if [ ! -x "$PY" ]; then
    fail "production Python is not executable: $PY"
else
    dep_out="$("$PY" -B - $REQUIRED_MODULES <<'PYDEPS'
import importlib.util, sys
from importlib import metadata
for name in sys.argv[1:]:
    if importlib.util.find_spec(name) is None:
        print("MISSING", name)
        continue
    version = "unknown"
    for cand in (name, name.capitalize(), name.title()):
        try:
            version = metadata.version(cand)
            break
        except metadata.PackageNotFoundError:
            pass
    print("OK", name, version)
PYDEPS
)" || dep_out=""
    if [ -z "$dep_out" ]; then
        fail "could not run the dependency check with $PY"
    else
        while read -r status name version; do
            case "$status" in
                MISSING) fail "required Python module is missing: $name (nothing is installed by this script)" ;;
                OK)
                    pass "$name $version"
                    if [ "$name" = "flask" ] && [ "$version" != "$KNOWN_GOOD_FLASK" ]; then
                        warn "flask $version differs from the validated $KNOWN_GOOD_FLASK"
                    fi
                    if [ "$name" = "werkzeug" ] && [ "$version" != "$KNOWN_GOOD_WERKZEUG" ]; then
                        warn "werkzeug $version differs from the validated $KNOWN_GOOD_WERKZEUG"
                    fi
                    ;;
            esac
        done <<< "$dep_out"
    fi
fi

# ── 5. systemd unit drift (report only) ──────────────────────────────────────
section "5. systemd unit (report only, never modified)"
if command -v systemctl >/dev/null 2>&1; then
    unit_text="$(systemctl cat "$PORTAL_UNIT" 2>/dev/null || true)"
    if [ -z "$unit_text" ]; then
        info "unit text not available; skipped"
    else
        drift=0
        expect_line() {  # KEY VALUE
            if [ "$(printf '%s\n' "$unit_text" | grep -c -x -F -- "$1=$2")" -eq 0 ]; then
                info "UNIT DRIFT: expected '$1=$2'"; drift=1
            fi
        }
        expect_line User rami
        expect_line WorkingDirectory "$DST"
        expect_line ExecStart "$PY app.py"
        expect_line Restart on-failure
        expect_line WantedBy multi-user.target
        if [ "$drift" -eq 0 ]; then
            pass "installed unit matches the expected settings"
        else
            warn "the installed unit differs from the expected settings; change it separately if intended"
        fi
    fi
fi

# ── 6. plan ──────────────────────────────────────────────────────────────────
section "6. Plan"
declare -A OLD_SHA ACTION
TO_WRITE=()
n_new=0; n_changed=0; n_same=0
if [ "$SOURCE_OK" -eq 1 ] && [ -d "$DST" ]; then
    printf '  %-9s %-30s %-12s %s\n' STATUS FILE "OLD SHA-256" "NEW SHA-256"
    for rel in "${ALLOWLIST[@]}"; do
        dst="$DST/$rel"
        if [ ! -e "$dst" ]; then
            ACTION["$rel"]="NEW"; OLD_SHA["$rel"]="-"; n_new=$((n_new + 1)); TO_WRITE+=("$rel")
        elif [ -L "$dst" ] || [ ! -f "$dst" ]; then
            fail "destination is not a regular file: $rel"; continue
        else
            OLD_SHA["$rel"]="$(sha256sum "$dst" | cut -d' ' -f1)"
            if [ "${OLD_SHA[$rel]}" = "${NEW_SHA[$rel]}" ]; then ACTION["$rel"]="same"; n_same=$((n_same + 1))
            else ACTION["$rel"]="CHANGED"; n_changed=$((n_changed + 1)); TO_WRITE+=("$rel"); fi
        fi
        printf '  %-9s %-30s %-12s %s\n' "${ACTION[$rel]}" "$rel" "$(printf '%.12s' "${OLD_SHA[$rel]}")" "$(printf '%.12s' "${NEW_SHA[$rel]}")"
    done
    info "$n_changed changed, $n_new new, $n_same unchanged, ${#ALLOWLIST[@]} allowlisted"
    info "never written: .env, .env.example, *.db, uploads/, portal_logs/, results/, .snakemake/, workflow/, config/, resources/, scripts/, tests/"
fi

# ── decision ─────────────────────────────────────────────────────────────────
section "Result"
if [ "${#WARNINGS[@]}" -gt 0 ]; then info "${#WARNINGS[@]} warning(s)"; fi
if [ "${#FAILURES[@]}" -gt 0 ]; then
    if [ "$EXECUTE" -eq 1 ]; then
        echo "  REFUSED: ${#FAILURES[@]} check(s) failed; nothing was written."
    else
        echo "  DRY-RUN: ${#FAILURES[@]} check(s) failed; --execute would be refused. Nothing was written."
    fi
    exit 1
fi
if [ "$EXECUTE" -eq 0 ]; then
    echo "  DRY-RUN complete: all checks passed; nothing was written."
    echo "  To deploy the plan above run again with --execute (pin the source with --expect-commit)."
    exit 0
fi
if [ "${#TO_WRITE[@]}" -eq 0 ]; then
    echo "  Nothing to deploy: production already matches the source. No backup created."
    exit 0
fi

# ── execute ──────────────────────────────────────────────────────────────────
exec 9>"$LOCK_FILE" || die "cannot open the deployment lock $LOCK_FILE"
command -v flock >/dev/null 2>&1 || die "flock is not available"
flock -n 9 || die "another deployment is running (lock $LOCK_FILE)"

TMP_FILES=()
cleanup() {
    local f
    for f in "${TMP_FILES[@]:-}"; do
        if [ -n "$f" ]; then rm -f -- "$f"; fi
    done
}
trap cleanup EXIT

TS="$(date -u +%Y%m%dT%H%M%SZ)"
BACKUP_DIR="$BACKUP_PARENT/$TS"
section "Backup"
[ ! -e "$BACKUP_DIR" ] || die "backup directory already exists: $BACKUP_DIR"
( umask 077; mkdir -p -- "$BACKUP_DIR/files" ) || die "cannot create the backup directory"
chmod 700 "$HOME_DIR/deploy_backups" "$BACKUP_PARENT" "$BACKUP_DIR" "$BACKUP_DIR/files" 2>/dev/null || die "cannot restrict backup permissions"
{
    echo "GenRichi portal deployment manifest"
    echo "utc_time=$TS"
    if [ "$FROM_WORKTREE" -eq 1 ]; then
        echo "source_mode=worktree"
        echo "source_bytes=WORKING-TREE files (not guaranteed to equal any commit)"
        echo "source_commit=none (working-tree bytes are not a commit)"
        echo "worktree_base_head=$COMMIT"
        echo "worktree_uncommitted_entries=$(printf '%s' "$DIRTY_SCOPE" | grep -c . || true)"
    else
        echo "source_mode=head"
        echo "source_bytes=exact blobs of commit $COMMIT"
        echo "source_commit=$COMMIT"
    fi
    echo "source_repo=$REPO"
    echo "operator=$(id -un)"
    echo "destination=$DST"
    echo "# action relpath old_sha256 new_sha256 (backups exist only for CHANGED files)"
} > "$BACKUP_DIR/MANIFEST.txt"
for rel in "${TO_WRITE[@]}"; do
    echo "${ACTION[$rel]} $rel ${OLD_SHA[$rel]} ${NEW_SHA[$rel]}" >> "$BACKUP_DIR/MANIFEST.txt"
    if [ "${ACTION[$rel]}" = "CHANGED" ]; then
        ( umask 077; mkdir -p -- "$(dirname "$BACKUP_DIR/files/$rel")" )
        cp -p -- "$DST/$rel" "$BACKUP_DIR/files/$rel" || die "backup copy failed for $rel; nothing was deployed"
        [ "$(sha256sum "$BACKUP_DIR/files/$rel" | cut -d' ' -f1)" = "${OLD_SHA[$rel]}" ] || die "backup verification failed for $rel; nothing was deployed"
    fi
done
chmod 600 "$BACKUP_DIR/MANIFEST.txt"
pass "backup verified byte-for-byte: $BACKUP_DIR"

section "Deploy"
for rel in "${TO_WRITE[@]}"; do
    dst="$DST/$rel"
    ddir="$(dirname "$dst")"
    if [ ! -d "$ddir" ]; then
        path_has_symlink "$ddir" && die "refusing to create a directory through a symlink: $ddir"
        mkdir -p -- "$ddir"
    fi
    tmp="$(mktemp "$ddir/.deploy.XXXXXX")"
    TMP_FILES+=("$tmp")
    src_cat "$rel" > "$tmp" || die "could not read $rel from the source; stopped"
    chmod 0644 "$tmp"
    [ "$(sha256sum "$tmp" | cut -d' ' -f1)" = "${NEW_SHA[$rel]}" ] || die "source changed while deploying $rel; stopped"
    mv -f -- "$tmp" "$dst" || die "could not replace $rel; stopped"
    [ "$(sha256sum "$dst" | cut -d' ' -f1)" = "${NEW_SHA[$rel]}" ] || die "verification failed after writing $rel; stopped"
    pass "${ACTION[$rel]} $rel"
done

section "Done"
if [ "$FROM_WORKTREE" -eq 1 ]; then
    echo "  Deployed ${#TO_WRITE[@]} file(s) from the WORKING TREE (HEAD was $COMMIT; the bytes are not guaranteed to equal any commit)."
else
    echo "  Deployed ${#TO_WRITE[@]} file(s) from commit $COMMIT."
fi
echo "  The running portal was NOT restarted and no service or unit was changed."
echo "  Python modules stay old in memory, and templates not yet rendered by the running"
echo "  process load the new file on first use, until a controlled restart is performed"
echo "  separately."
echo "  Verify with: bash $PROD_REPO/scripts/portal_status.sh"
echo "  Backup and manifest: $BACKUP_DIR"
echo "  To roll back the CHANGED files:"
echo "    (cd \"$BACKUP_DIR/files\" && find . -type f -exec cp -p -- {} \"$DST/{}\" \\;)"
echo "  (files reported as NEW have no backup; remove them by hand if needed)"
