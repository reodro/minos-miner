#!/usr/bin/env bash
# Gets/renews the dashboard's Let's Encrypt certificate for the bare IP (short-lived, ~6 days; renews at half-life)
# and installs it into dashboard/.certs/. Run twice a day by the systemd timer minos-dashboard-cert.timer.
# The http-01 challenge is served by dashboard/redirect.py on port 80.
set -euo pipefail
IP=${DASHBOARD_IP:-5.196.64.61}
HERE=$(cd "$(dirname "$0")" && pwd)
LEGO=${LEGO:-/opt/lego/lego}
DATA=/opt/lego/data
mkdir -p "$HERE/.certs/webroot"
"$LEGO" run --accept-tos --path "$DATA" -d "$IP" --http --http.webroot "$HERE/.certs/webroot" --profile shortlived
if ! cmp -s "$DATA/certificates/$IP.crt" "$HERE/.certs/cert.pem"; then
    install -m 600 "$DATA/certificates/$IP.crt" "$HERE/.certs/cert.pem"
    install -m 600 "$DATA/certificates/$IP.key" "$HERE/.certs/key.pem"
    systemctl try-restart minos-dashboard 2>/dev/null || true   # the dashboard, when one is installed
    echo "installed new certificate"
fi
