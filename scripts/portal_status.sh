#!/usr/bin/env bash
# GenRichi portal status -- STRICTLY READ-ONLY.
#
# Only queries state: systemd properties, listening sockets and HTTP GET
# requests. It changes nothing (no service, process, file, or database).
# systemd is the only supervisor of the production portal and tunnel.
#
# Usage:  bash scripts/portal_status.sh [--wait SECONDS] [--no-public]
#   --wait N     poll for up to N seconds until the local /login answers 200
#   --no-public  skip the request to the public URL
#
# Exit code: 0 = healthy, 1 = not healthy. The exit code is informational;
# nothing is repaired by this script.

set -u

PORTAL_UNIT="genrichi-portal.service"
DUPLICATE_UNIT="genrichi.service"
TUNNEL_UNIT="cloudflared.service"
PORT=5000
LOCAL_URL="http://127.0.0.1:${PORT}/login"
PUBLIC_URL="https://portal.genrichi.de/login"

WAIT_SECONDS=0
CHECK_PUBLIC=1
while [ $# -gt 0 ]; do
    case "$1" in
        --wait)
            if [ $# -lt 2 ]; then echo "usage error: --wait needs a whole number of seconds, e.g. --wait 90" >&2; exit 2; fi
            WAIT_SECONDS="$2"; shift 2 ;;
        --no-public) CHECK_PUBLIC=0; shift ;;
        -h|--help)   sed -n '2,15p' "$0"; exit 0 ;;
        *)           echo "unknown option: $1" >&2; exit 2 ;;
    esac
done
case "$WAIT_SECONDS" in ''|*[!0-9]*) echo "usage error: --wait needs a whole number of seconds, e.g. --wait 90" >&2; exit 2 ;; esac

problems=0
note() { echo "  $*"; }
bad()  { echo "  PROBLEM: $*"; problems=$((problems + 1)); }

unit_state() {
    local unit="$1" active enabled pid restarts
    active=$(systemctl is-active "$unit" 2>/dev/null || true)
    enabled=$(systemctl is-enabled "$unit" 2>/dev/null || true)
    pid=$(systemctl show "$unit" -p MainPID --value 2>/dev/null || true)
    restarts=$(systemctl show "$unit" -p NRestarts --value 2>/dev/null || true)
    printf '  %-26s active=%-11s enabled=%-9s MainPID=%-7s NRestarts=%s\n' \
        "$unit" "${active:-unknown}" "${enabled:-unknown}" "${pid:-?}" "${restarts:-?}"
    UNIT_ACTIVE="$active"
    UNIT_ENABLED="$enabled"
    UNIT_PID="$pid"
}

http_code() { curl -s -o /dev/null -m 10 -w '%{http_code}' "$1" 2>/dev/null || echo 000; }

echo "GenRichi portal status (read-only)  $(date '+%Y-%m-%d %H:%M:%S')"
echo

echo "systemd units"
if ! command -v systemctl >/dev/null 2>&1; then
    bad "systemctl not available in this shell"
else
    unit_state "$PORTAL_UNIT";    portal_active="$UNIT_ACTIVE"; portal_pid="$UNIT_PID"
    unit_state "$TUNNEL_UNIT";    tunnel_active="$UNIT_ACTIVE"
    unit_state "$DUPLICATE_UNIT"; dup_active="$UNIT_ACTIVE";    dup_enabled="$UNIT_ENABLED"
    [ "$portal_active" = "active" ] || bad "$PORTAL_UNIT is not active"
    [ "$tunnel_active" = "active" ] || bad "$TUNNEL_UNIT is not active"
    case "$dup_active" in active|activating) bad "$DUPLICATE_UNIT is running; it must stay disabled" ;; esac
    case "$dup_enabled" in disabled|masked|"not-found"|"") ;; *) bad "$DUPLICATE_UNIT is '$dup_enabled'; it must stay disabled" ;; esac
fi
echo

echo "port ${PORT}"
owners=$(ss -ltnp 2>/dev/null | awk -v p=":${PORT}\$" '$4 ~ p' | grep -o 'pid=[0-9]*' | cut -d= -f2 | sort -u)
count=$(printf '%s\n' "$owners" | grep -c .)
if [ "$count" -eq 0 ]; then
    bad "nothing is listening on port ${PORT}"
else
    for pid in $owners; do
        unit=$(grep -o -E '[A-Za-z0-9@_.-]+\.service' "/proc/${pid}/cgroup" 2>/dev/null | head -1)
        note "pid=${pid}  unit=${unit:-unknown}"
        if [ "${unit:-}" != "$PORTAL_UNIT" ]; then bad "port ${PORT} is held by ${unit:-an unknown unit}, not ${PORTAL_UNIT}"; fi
    done
    [ "$count" -eq 1 ] || bad "port ${PORT} has ${count} owners; expected exactly one"
fi
echo

echo "HTTP"
local_code=$(http_code "$LOCAL_URL")
waited=0
while [ "$local_code" != "200" ] && [ "$waited" -lt "$WAIT_SECONDS" ]; do
    sleep 2
    waited=$((waited + 2))
    local_code=$(http_code "$LOCAL_URL")
done
if [ "$waited" -gt 0 ]; then
    note "GET ${LOCAL_URL} -> ${local_code} (waited ${waited}s)"
else
    note "GET ${LOCAL_URL} -> ${local_code}"
fi
[ "$local_code" = "200" ] || bad "local /login did not answer 200"
if [ "$CHECK_PUBLIC" -eq 1 ]; then
    public_code=$(http_code "$PUBLIC_URL")
    note "GET ${PUBLIC_URL} -> ${public_code}"
    [ "$public_code" = "200" ] || bad "public /login did not answer 200"
fi
echo

if [ "$problems" -eq 0 ]; then
    echo "RESULT: healthy"
    exit 0
fi
echo "RESULT: ${problems} problem(s) found (report only; nothing was changed)"
exit 1
