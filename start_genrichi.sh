#!/usr/bin/env bash
# GenRichi — Start portal + tunnel
set -euo pipefail

PORTAL="/home/rami/genrichi/portal"
PYTHON="/home/rami/miniforge3/envs/snakemake/bin/python"

echo "=============================="
echo "  GenRichi Clinical Portal"
echo "=============================="

# Kill any old instances
pkill -f "python app.py"    2>/dev/null || true
pkill -f "cloudflared"      2>/dev/null || true
sleep 1

# Start Flask
cd "$PORTAL"
nohup "$PYTHON" app.py > "$PORTAL/portal.log" 2>&1 &
echo "Flask started (PID $!)"

# Start Tunnel
nohup cloudflared tunnel run genrichi-portal > "$PORTAL/tunnel.log" 2>&1 &
echo "Tunnel started (PID $!)"

sleep 3

# Verify
if pgrep -f "python app.py" > /dev/null; then
    echo "OK  Flask running"
else
    echo "ERR Flask failed — check portal.log"
fi

if pgrep -f "cloudflared" > /dev/null; then
    echo "OK  Tunnel running"
else
    echo "ERR Tunnel failed — check tunnel.log"
fi

echo ""
echo "Portal: https://portal.genrichi.de"
echo "Login:  admin / GenRichi2026!"
echo ""
