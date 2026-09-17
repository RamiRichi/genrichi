#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# GenRichi Portal — Deploy to /home/rami/genrichi/portal/
#
# Run from WSL:
#   bash /mnt/c/GenRichi/deploy_portal.sh
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

SRC="/mnt/c/GenRichi/portal"
DST="/home/rami/genrichi/portal"

GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; NC='\033[0m'
ok()   { echo -e "${GREEN}  OK${NC}  $1"; }
warn() { echo -e "${YELLOW}WARN${NC}  $1"; }
err()  { echo -e "${RED} ERR${NC}  $1"; exit 1; }

echo ""
echo "=== GenRichi Portal Deploy ==="
echo "  Source : $SRC"
echo "  Target : $DST"
echo ""

[ -d "/home/rami/genrichi" ] || err "/home/rami/genrichi not found."

echo "[1] Portal files"
mkdir -p "$DST/templates" "$DST/static" "$DST/uploads" "$DST/portal_logs"
cp "$SRC/app.py"     "$DST/app.py"
cp "$SRC/models.py"  "$DST/models.py"
cp "$SRC/runner.py"  "$DST/runner.py"
cp "$SRC/config.py"  "$DST/config.py"
ok "Python files"

cp "$SRC/templates/"*.html "$DST/templates/"
ok "HTML templates"

echo ""
echo "[2] Install Flask (into snakemake conda env)"
/home/rami/miniforge3/envs/snakemake/bin/pip install flask werkzeug --quiet
ok "Flask installed"

echo ""
echo "[3] Systemd service (optional — manual start also works)"
cat > /tmp/genrichi-portal.service << 'SVCEOF'
[Unit]
Description=GenRichi Clinical Genomics Portal
After=network.target

[Service]
User=rami
WorkingDirectory=/home/rami/genrichi/portal
ExecStart=/home/rami/miniforge3/envs/snakemake/bin/python app.py
Restart=on-failure
RestartSec=5

[Install]
WantedBy=multi-user.target
SVCEOF

if command -v systemctl &>/dev/null && systemctl --version &>/dev/null 2>&1; then
    sudo cp /tmp/genrichi-portal.service /etc/systemd/system/
    sudo systemctl daemon-reload
    sudo systemctl enable genrichi-portal --quiet
    ok "Systemd service installed (genrichi-portal)"
else
    warn "systemd not available (WSL) — use manual start below"
fi

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo -e "${GREEN}Deploy complete!${NC}"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Start the portal:"
echo ""
echo "  conda activate snakemake"
echo "  cd /home/rami/genrichi/portal"
echo "  python app.py"
echo ""
echo "Then open in your browser:"
echo "  http://localhost:5000"
echo ""
echo "Login:"
echo "  Username : admin"
echo "  Password : GenRichi2026!"
echo ""
echo "⚠  Change PORTAL_PASS in portal/config.py before sharing with others."
echo ""
