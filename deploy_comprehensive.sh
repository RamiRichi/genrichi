#!/usr/bin/env bash
# GenRichi Comprehensive pipeline deployment -- DRY-RUN by default.
#
# Replaces the old, unhardened deploy_comprehensive.sh. Copies ONLY the
# explicit allowlist of Comprehensive pipeline code, config and panel
# resource files from a Git commit into /home/rami/genrichi. It never
# touches reference_db/, envs/, results/, .snakemake/, the portal, .env,
# the database, or any other pipeline's files, and it never restarts
# anything.
#
# Usage:  bash deploy_comprehensive.sh [options]
#   (no option)         dry-run: run every check, print the plan, write nothing
#   --execute           perform the deployment (backup first, atomic writes,
#                        then write a Comprehensive-scoped build stamp)
#   --expect-commit SHA refuse unless HEAD is exactly this full 40-hex commit SHA
#   --help               show this text
#
# Exit code: 0 ok, 1 refused / a check failed, 2 usage error.

set -euo pipefail
umask 022

# ── fixed production locations ───────────────────────────────────────────────
readonly PROD_HOME="/home/rami/genrichi"
readonly PROD_REPO="/mnt/c/GenRichi"
readonly EXPECTED_REMOTE_MATCH="RamiRichi/genrichi"
readonly PORTAL_UNIT="genrichi-portal.service"
readonly PIPELINE_ID="comprehensive"
readonly BUILD_STAMP_RELPATH="workflow/build_info.json"

# A. Comprehensive source code: the Snakefile, every rule file it includes,
#    and every script those rules (or provenance.py) import. This list is
#    cross-checked against the actual `include:` lines in the deployed
#    Snakefile at check time -- it is not trusted blindly forever.
readonly ALLOWLIST_CODE=(
    workflow/Snakefile_comprehensive
    workflow/rules/runtime_versions.smk
    workflow/rules/qc_paired.smk
    workflow/rules/align_paired.smk
    workflow/rules/somatic_calling_paired.smk
    workflow/rules/cnv_calling.smk
    workflow/rules/msi_scoring.smk
    workflow/rules/annotation_comprehensive.smk
    workflow/rules/report_comprehensive.smk
    workflow/scripts/calculate_cnv.py
    workflow/scripts/calculate_msi.py
    workflow/scripts/vcf_to_somatic_table.py
    workflow/scripts/generate_comprehensive_report.py
    workflow/scripts/provenance.py
    workflow/scripts/resource_inspect.py
    workflow/scripts/qc_status.py
    workflow/scripts/sample_columns.py
    workflow/scripts/runtime_probe.py
)
# B. Configuration + panel content the code above depends on. The two BED
#    paths are the exact values comprehensive_config.yaml must reference --
#    checked below -- specifically to catch the class of drift found in the
#    provenance investigation (config pointing at one panel BED, the old
#    unhardened script deploying a different, retired one).
readonly ALLOWLIST_CONFIG=(
    config/comprehensive_config.yaml
    config/comprehensive_samples.tsv
    resources/panel/phase1_solid_tumor/solid_tumor_phase1_v1.bed
    resources/panel/comprehensive_genes_cds.bed
)
readonly ALLOWLIST=("${ALLOWLIST_CODE[@]}" "${ALLOWLIST_CONFIG[@]}")
readonly EXPECTED_PANEL_BED="resources/panel/phase1_solid_tumor/solid_tumor_phase1_v1.bed"
readonly EXPECTED_TMB_BED="resources/panel/comprehensive_genes_cds.bed"

# ── test mode (unit tests only; every path must live under /tmp) ─────────────
TEST_MODE=0
HOME_DIR="$PROD_HOME"
REPO="$PROD_REPO"
CGROUP_ROOT="/sys/fs/cgroup"
if [ -n "${GENRICHI_DEPLOY_TEST_HOME:-}" ]; then
    TEST_MODE=1
    HOME_DIR="$GENRICHI_DEPLOY_TEST_HOME"
    REPO="${GENRICHI_DEPLOY_TEST_REPO:?test mode needs GENRICHI_DEPLOY_TEST_REPO}"
    CGROUP_ROOT="${GENRICHI_DEPLOY_TEST_CGROUP_ROOT:?test mode needs GENRICHI_DEPLOY_TEST_CGROUP_ROOT}"
    for p in "$HOME_DIR" "$REPO" "$CGROUP_ROOT"; do
        case "$p" in /tmp/*) ;; *) echo "test mode only accepts paths under /tmp: refusing" >&2; exit 2 ;; esac
    done
fi
BACKUP_PARENT="$HOME_DIR/deploy_backups/$PIPELINE_ID"
ENV_FILE="$HOME_DIR/portal/.env"
LOCK_FILE="${TMPDIR:-/tmp}/genrichi-comprehensive-deploy.lock"
STAMP_PATH="$HOME_DIR/$BUILD_STAMP_RELPATH"

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
GenRichi Comprehensive pipeline deployment (dry-run unless --execute is given)

  bash deploy_comprehensive.sh [--execute] [--expect-commit FULL_40_HEX_SHA] [--help]

Deploys only the fixed allowlist of Comprehensive pipeline code, config and
panel resource files, from Git HEAD blobs. It never writes .env, the
database, reference_db/, envs/, results/, uploads, logs or another
pipeline's files, never installs packages, and never restarts anything.
On --execute it also writes a Comprehensive-scoped build stamp
(workflow/build_info.json) recording the exact source commit and the
hashes of every file actually deployed.
USAGE
}

# ── arguments ────────────────────────────────────────────────────────────────
EXECUTE=0
EXPECT_COMMIT=""
while [ $# -gt 0 ]; do
    case "$1" in
        --execute) EXECUTE=1 ;;
        --expect-commit)
            [ $# -ge 2 ] && [ -n "$2" ] || { echo "--expect-commit needs a value" >&2; exit 2; }
            EXPECT_COMMIT="$2"; shift ;;
        -h|--help) usage; exit 0 ;;
        *)         echo "unknown option: $1" >&2; usage >&2; exit 2 ;;
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

src_cat() { repo_git cat-file blob "$COMMIT:$1"; }  # repo-root-relative path -> bytes on stdout

# Lines (numbers only) that assign a quoted literal to a secret-named variable.
CRED_PATTERN='(pass|passwd|password|secret|token|api_?key)[a-z0-9_]*[[:space:]]*=[[:space:]]*["'"'"'][^"'"'"']{6,}["'"'"']'
CRED_ALLOWED='os[.]environ|getenv|_require_env|placeholder='
cred_lines() {
    src_cat "$1" | CRED_PATTERN="$CRED_PATTERN" CRED_ALLOWED="$CRED_ALLOWED" \
        awk '{ l = tolower($0) } l ~ ENVIRON["CRED_PATTERN"] && l !~ ENVIRON["CRED_ALLOWED"] { print NR }'
}

# ── banner ───────────────────────────────────────────────────────────────────
if [ "$EXECUTE" -eq 1 ]; then MODE_LABEL="EXECUTE"; else MODE_LABEL="DRY-RUN (nothing will be written)"; fi
echo "GenRichi Comprehensive pipeline deployment -- $MODE_LABEL"
if [ "$TEST_MODE" -eq 1 ]; then echo "  *** TEST MODE: paths are redirected under /tmp ***"; fi
info "destination : $HOME_DIR"
info "pipeline    : $PIPELINE_ID"
info "source mode : Git HEAD commit"

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
            pass "commit $COMMIT ($(repo_git log -1 --format='%cs' "$COMMIT" 2>/dev/null)); file bytes are read from this commit"
            info "subject: $(repo_git log -1 --format='%s' "$COMMIT" 2>/dev/null | cut -c1-90)"
            if [ -n "$EXPECT_COMMIT" ]; then
                if [ "$COMMIT" = "$EXPECT_COMMIT" ]; then
                    pass "HEAD matches --expect-commit $EXPECT_COMMIT"
                else
                    fail "HEAD $COMMIT does not match --expect-commit $EXPECT_COMMIT"
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

# dirty-scope limited to exactly the paths this script can deploy
DIRTY_SCOPE=""
if [ -n "$COMMIT" ]; then
    DIRTY_SCOPE="$(repo_git status --porcelain --untracked-files=normal -- \
        workflow/Snakefile_comprehensive workflow/rules workflow/scripts \
        config/comprehensive_config.yaml config/comprehensive_samples.tsv \
        resources/panel/phase1_solid_tumor resources/panel/comprehensive_genes_cds.bed \
        2>/dev/null || true)"
    if [ -n "$DIRTY_SCOPE" ]; then
        info "uncommitted changes in the deployable scope; they are IGNORED because HEAD blobs are deployed"
    else
        pass "deployable scope's working tree is clean"
    fi
fi

# every allowlisted file must exist in the commit (fail closed on anything missing)
SOURCE_OK=0
if [ -n "$COMMIT" ]; then
    missing=0
    for rel in "${ALLOWLIST[@]}"; do
        repo_git cat-file -e "$COMMIT:$rel" 2>/dev/null || { fail "allowlisted file missing in commit: $rel"; missing=$((missing + 1)); }
    done

    # the allowlist must exactly match the Snakefile's own declared rule-file
    # dependency graph -- if the Snakefile now includes a rule file this
    # script doesn't know about (or no longer includes one it still lists),
    # fail closed rather than silently deploying a partial or stale pipeline.
    include_mismatch=0
    if repo_git cat-file -e "$COMMIT:workflow/Snakefile_comprehensive" 2>/dev/null; then
        DECLARED_RULES="$(src_cat workflow/Snakefile_comprehensive \
            | grep -oE '^include:[[:space:]]*"rules/[^"]+\.smk"' \
            | sed -E 's/^include:[[:space:]]*"rules\/([^"]+)\.smk"/workflow\/rules\/\1.smk/')"
        ALLOWLISTED_RULES="$(printf '%s\n' "${ALLOWLIST_CODE[@]}" | grep '^workflow/rules/' | sort)"
        DECLARED_SORTED="$(printf '%s\n' "$DECLARED_RULES" | sort)"
        if [ "$DECLARED_SORTED" != "$ALLOWLISTED_RULES" ]; then
            fail "Snakefile_comprehensive's include: list no longer matches ALLOWLIST_CODE -- unexpected source dependency, refusing to guess"
            diff <(printf '%s\n' "$ALLOWLISTED_RULES") <(printf '%s\n' "$DECLARED_SORTED") | sed 's/^/    /' || true
            include_mismatch=1
        fi
    fi

    if [ "$missing" -eq 0 ] && [ "$include_mismatch" -eq 0 ]; then
        pass "all ${#ALLOWLIST[@]} allowlisted files present; Snakefile's rule includes match ALLOWLIST_CODE exactly"
        SOURCE_OK=1
    fi
fi

# config content must reference exactly the allowlisted panel BEDs --
# the exact drift class found in the provenance investigation.
if [ "$SOURCE_OK" -eq 1 ]; then
    cfg_panel_bed="$(src_cat config/comprehensive_config.yaml | awk '/^panel:/{p=1;next} p && /^[a-zA-Z]/{p=0} p && /^[[:space:]]+bed:/{sub(/^[[:space:]]+bed:[[:space:]]*/,""); sub(/[[:space:]]*#.*/,""); print; exit}')"
    cfg_tmb_bed="$(src_cat config/comprehensive_config.yaml | awk '/^tmb:/{p=1;next} p && /^[a-zA-Z]/{p=0} p && /^[[:space:]]+coding_bed:/{sub(/^[[:space:]]+coding_bed:[[:space:]]*/,""); sub(/[[:space:]]*#.*/,""); print; exit}')"
    if [ "$cfg_panel_bed" != "$EXPECTED_PANEL_BED" ]; then
        fail "comprehensive_config.yaml panel.bed is '$cfg_panel_bed', expected the allowlisted '$EXPECTED_PANEL_BED'"
    elif [ "$cfg_tmb_bed" != "$EXPECTED_TMB_BED" ]; then
        fail "comprehensive_config.yaml tmb.coding_bed is '$cfg_tmb_bed', expected the allowlisted '$EXPECTED_TMB_BED'"
    else
        pass "comprehensive_config.yaml panel/TMB BED references match the allowlisted resource files exactly"
    fi
fi

# content checks on the exact bytes that would be deployed
declare -A NEW_SHA
if [ "$SOURCE_OK" -eq 1 ]; then
    bad_content=0
    for rel in "${ALLOWLIST[@]}"; do
        NEW_SHA["$rel"]="$(src_cat "$rel" | sha_stdin)"
        crs="$(src_cat "$rel" | tr -cd '\r' | wc -c)"
        if [ "$crs" -ne 0 ]; then fail "carriage returns in $rel (CRLF would corrupt the deployed file)"; bad_content=1; fi
        case "$rel" in
            *.py)
                hits="$(cred_lines "$rel" | tr '\n' ' ')"
                if [ -n "$hits" ]; then fail "credential-like literal in $rel at line(s): $hits(value not shown)"; bad_content=1; fi
                if ! src_cat "$rel" | python3 -B -c 'import ast, sys; ast.parse(sys.stdin.buffer.read())' 2>/dev/null; then
                    fail "Python syntax check failed for $rel"; bad_content=1
                fi ;;
            *.yaml)
                hits="$(cred_lines "$rel" | tr '\n' ' ')"
                if [ -n "$hits" ]; then fail "credential-like literal in $rel at line(s): $hits(value not shown)"; bad_content=1; fi ;;
        esac
    done
    if [ "$bad_content" -eq 0 ]; then pass "source bytes: LF only, no credential-like literals, Python files parse"; fi
fi

# ── 2. destination ───────────────────────────────────────────────────────────
section "2. Destination and production guards"
if [ "$TEST_MODE" -eq 0 ] && [ "$HOME_DIR" != "/home/rami/genrichi" ]; then
    fail "destination is not exactly /home/rami/genrichi"
elif [ ! -d "$HOME_DIR" ]; then
    fail "destination directory does not exist: $HOME_DIR"
elif path_has_symlink "$HOME_DIR"; then
    fail "a component of the destination path is a symlink"
else
    pass "destination is the expected directory and contains no symlink component"
fi
if [ "$(id -u)" -eq 0 ]; then
    fail "refusing to run as root (deployed files must stay owned by the service user)"
elif [ -d "$HOME_DIR" ] && [ "$(stat -c '%U' "$HOME_DIR")" != "$(id -un)" ]; then
    fail "destination is owned by $(stat -c '%U' "$HOME_DIR"), not by $(id -un)"
else
    pass "running as $(id -un), the owner of the destination"
fi
for rel in "${ALLOWLIST[@]}"; do
    dst="$HOME_DIR/$rel"
    if path_has_symlink "$(dirname "$dst")"; then fail "destination directory for $rel contains a symlink component"; fi
    if [ -e "$dst" ] && ! [ -f "$dst" ]; then fail "destination path is not a regular file: $rel"; fi
done
[ -d "$HOME_DIR" ] && [ ! -w "$HOME_DIR" ] && fail "destination is not writable"

if [ -L "$ENV_FILE" ] || [ ! -f "$ENV_FILE" ]; then
    fail "production .env is missing (or not a regular file); it is never touched by this script"
else
    pass "production .env exists (contents never read, never touched by this script)"
fi

# ── 3. running analyses (read-only queries only; never kill, never touch locks) ──
section "3. Active pipeline execution"
if ! command -v systemctl >/dev/null 2>&1; then
    fail "systemctl is not available: cannot verify the portal/pipeline state"
else
    portal_state="$(systemctl is-active "$PORTAL_UNIT" 2>/dev/null || true)"
    info "$PORTAL_UNIT is ${portal_state:-unknown} (never restarted by this script)"
    if [ "$portal_state" = "active" ]; then
        cg="$(systemctl show "$PORTAL_UNIT" -p ControlGroup --value 2>/dev/null || true)"
        procs_file="$CGROUP_ROOT$cg/cgroup.procs"
        if [ -z "$cg" ] || [ ! -r "$procs_file" ]; then
            fail "cannot read the cgroup of $PORTAL_UNIT: unable to prove that no Comprehensive run is active"
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
if [ -e "$HOME_DIR/.snakemake/locks" ]; then
    fail "$HOME_DIR/.snakemake/locks exists (a run may hold the directory) -- refusing; this script never removes or touches locks"
else
    pass "no .snakemake/locks directory"
fi

# ── 4. plan ──────────────────────────────────────────────────────────────────
section "4. Plan"
declare -A OLD_SHA ACTION
TO_WRITE=()
n_new=0; n_changed=0; n_same=0
if [ "$SOURCE_OK" -eq 1 ] && [ -d "$HOME_DIR" ]; then
    printf '  %-9s %-55s %-12s %s\n' STATUS FILE "OLD SHA-256" "NEW SHA-256"
    for rel in "${ALLOWLIST[@]}"; do
        dst="$HOME_DIR/$rel"
        if [ ! -e "$dst" ]; then
            ACTION["$rel"]="NEW"; OLD_SHA["$rel"]="-"; n_new=$((n_new + 1)); TO_WRITE+=("$rel")
        elif [ -L "$dst" ] || [ ! -f "$dst" ]; then
            fail "destination is not a regular file: $rel"; continue
        else
            OLD_SHA["$rel"]="$(sha256sum "$dst" | cut -d' ' -f1)"
            if [ "${OLD_SHA[$rel]}" = "${NEW_SHA[$rel]}" ]; then ACTION["$rel"]="same"; n_same=$((n_same + 1))
            else ACTION["$rel"]="CHANGED"; n_changed=$((n_changed + 1)); TO_WRITE+=("$rel"); fi
        fi
        printf '  %-9s %-55s %-12s %s\n' "${ACTION[$rel]}" "$rel" "$(printf '%.12s' "${OLD_SHA[$rel]}")" "$(printf '%.12s' "${NEW_SHA[$rel]}")"
    done
    info "$n_changed changed, $n_new new, $n_same unchanged, ${#ALLOWLIST[@]} allowlisted"
    info "never written: .env, database, reference_db/, envs/, results/, uploads/, portal_logs/, logs/, .snakemake/, retired resources/panel/comprehensive_genes.bed"
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
    echo "  Nothing to deploy: production already matches the source. No backup created, no build stamp updated."
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
    echo "GenRichi Comprehensive pipeline deployment manifest"
    echo "pipeline=$PIPELINE_ID"
    echo "utc_time=$TS"
    echo "source_mode=head"
    echo "source_bytes=exact blobs of commit $COMMIT"
    echo "source_commit=$COMMIT"
    echo "source_repo=$REPO"
    echo "operator=$(id -un)"
    echo "destination=$HOME_DIR"
    echo "# action relpath old_sha256 new_sha256 (backups exist only for CHANGED files)"
} > "$BACKUP_DIR/MANIFEST.txt"
for rel in "${TO_WRITE[@]}"; do
    echo "${ACTION[$rel]} $rel ${OLD_SHA[$rel]} ${NEW_SHA[$rel]}" >> "$BACKUP_DIR/MANIFEST.txt"
    if [ "${ACTION[$rel]}" = "CHANGED" ]; then
        ( umask 077; mkdir -p -- "$(dirname "$BACKUP_DIR/files/$rel")" )
        cp -p -- "$HOME_DIR/$rel" "$BACKUP_DIR/files/$rel" || die "backup copy failed for $rel; nothing was deployed"
        [ "$(sha256sum "$BACKUP_DIR/files/$rel" | cut -d' ' -f1)" = "${OLD_SHA[$rel]}" ] || die "backup verification failed for $rel; nothing was deployed"
    fi
done
chmod 600 "$BACKUP_DIR/MANIFEST.txt"
pass "backup verified byte-for-byte: $BACKUP_DIR"

section "Deploy"
for rel in "${TO_WRITE[@]}"; do
    dst="$HOME_DIR/$rel"
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

section "Build stamp"
# Scoped strictly to the Comprehensive pipeline -- never a repo-wide claim.
# Compatible with provenance.py's read_build_stamp()/resolve_pipeline_identity():
# git_commit_sha / git_describe / git_dirty / stamped_at_utc are the exact
# keys those functions read; the rest are additive and ignored by them.
GIT_DESCRIBE="$(repo_git describe --tags --always --dirty 2>/dev/null || true)"
REPO_DIRTY="false"; [ -n "$DIRTY_SCOPE" ] && REPO_DIRTY="true"
CODE_FP_INPUT=""
for rel in "${ALLOWLIST_CODE[@]}"; do
    CODE_FP_INPUT="${CODE_FP_INPUT}${rel}"$'\t'"${NEW_SHA[$rel]}"$'\n'
done
COMPREHENSIVE_CODE_FINGERPRINT="$(printf '%s' "$CODE_FP_INPUT" | sort | sha256sum | cut -d' ' -f1)"
{
    printf '{\n'
    printf '  "pipeline": "%s",\n' "$PIPELINE_ID"
    printf '  "deployment_id": "%s",\n' "$TS"
    printf '  "deployed_at_utc": "%s",\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf '  "stamped_at_utc": "%s",\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf '  "source_repo": "%s",\n' "$REPO"
    printf '  "source_commit": "%s",\n' "$COMMIT"
    printf '  "git_commit_sha": "%s",\n' "$COMMIT"
    printf '  "git_describe": %s,\n' "$([ -n "$GIT_DESCRIBE" ] && printf '"%s"' "$GIT_DESCRIBE" || echo null)"
    printf '  "git_dirty": %s,\n' "$REPO_DIRTY"
    printf '  "comprehensive_allowlist_count": %d,\n' "${#ALLOWLIST[@]}"
    printf '  "comprehensive_code_fingerprint_sha256": "%s",\n' "$COMPREHENSIVE_CODE_FINGERPRINT"
    printf '  "deployment_manifest": "%s",\n' "$BACKUP_DIR/MANIFEST.txt"
    printf '  "deployed_files": {\n'
    n=${#ALLOWLIST[@]}; i=0
    for rel in "${ALLOWLIST[@]}"; do
        i=$((i + 1))
        sep=","; [ "$i" -eq "$n" ] && sep=""
        printf '    "%s": "%s"%s\n' "$rel" "${NEW_SHA[$rel]}" "$sep"
    done
    printf '  }\n'
    printf '}\n'
} > "${STAMP_PATH}.tmp"
mv -f -- "${STAMP_PATH}.tmp" "$STAMP_PATH"
chmod 0644 "$STAMP_PATH"
pass "wrote $STAMP_PATH (pipeline=$PIPELINE_ID, commit ${COMMIT:0:12}, code fingerprint ${COMPREHENSIVE_CODE_FINGERPRINT:0:12})"

section "Done"
echo "  Deployed ${#TO_WRITE[@]} file(s) from commit $COMMIT for the Comprehensive pipeline."
echo "  No service was restarted and no unit was changed."
echo "  Snakemake picks up the new files on its next invocation (each run reads files fresh from disk)."
echo "  Build stamp: $STAMP_PATH"
echo "  Backup and manifest: $BACKUP_DIR"
echo "  To roll back the CHANGED files:"
echo "    (cd \"$BACKUP_DIR/files\" && find . -type f -exec cp -p -- {} \"$HOME_DIR/{}\" \\;)"
echo "  (files reported as NEW have no backup; remove them by hand if needed. The build stamp is"
echo "   not rolled back automatically -- rerun this script against the prior commit to restamp.)"
