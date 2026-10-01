# dashboard

Read-only web page showing what the miner is doing: the live round (answer set, evaluations, best score, nodes,
submissions, warnings) and the history of finished rounds (with true score and rank when known).

- **Address**: `https://<main's IP>` (port 443). The certificate is self-signed: the browser warns once.
- **Password**: `dashboard/.streamlit/secrets.toml` → `password = "..."` (git-ignored). Change it, then
  `systemctl restart minos-dashboard`.
- **Data**: main writes `<DATA_DIR>/status/live.json` and `rounds.jsonl` (`main/status.py`); the dashboard only reads
  them, so it can never slow down or break mining. Node names are anonymized (public page).
- **Service**: `systemctl status|restart minos-dashboard` (unit in `/etc/systemd/system/minos-dashboard.service`,
  runs `.venv/bin/streamlit run app.py` from this folder; config in `.streamlit/config.toml`).
- **Install** (main only, not the fleet nodes): `.venv/bin/pip install -r dashboard/requirements.txt`;
  certificate: `openssl req -x509 -newkey rsa:2048 -nodes -days 825 -subj "/CN=<IP>"
  -addext "subjectAltName=IP:<IP>" -keyout .certs/key.pem -out .certs/cert.pem` (inside `dashboard/`).

## http → https redirect

`dashboard/redirect.py` (systemd service `minos-dashboard-redirect`) answers on port 80 and redirects to `https://`,
so typing the bare address works. It is plain python3, no extra packages.
