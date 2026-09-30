#!/usr/bin/env bash
# Self-signed TLS certificate for main's fleet API. Run on main; copy fleet.crt (NOT fleet.key) to every node.
#   bash scripts/make_fleet_cert.sh <main public IP or DNS name> [out_dir]
set -euo pipefail
HOST="${1:?usage: make_fleet_cert.sh <main IP or DNS name> [out_dir]}"
OUT="${2:-$HOME/.minos/fleet}"
mkdir -p "$OUT" && chmod 700 "$OUT"
if [[ "$HOST" =~ ^[0-9.]+$ ]]; then SAN="IP:$HOST"; else SAN="DNS:$HOST"; fi
openssl req -x509 -newkey rsa:3072 -nodes -days 365 -subj "/CN=$HOST" -addext "subjectAltName=$SAN" \
  -keyout "$OUT/fleet.key" -out "$OUT/fleet.crt" 2>/dev/null
chmod 600 "$OUT/fleet.key"
echo "wrote $OUT/fleet.crt and $OUT/fleet.key - copy fleet.crt to the nodes, keep fleet.key on main"
