# dashboard

Read-only web page showing what the miner is doing:

- **Live round**: round and time left; evaluations (evaluated = configs that ran through GATK, i.e. used core time
  and budget; useful, wasted = scored 0, failed, running, evaluated per minute, configs per GATK run); one chart of
  every score over time with the best-so-far line and the submissions; the best config (score, core, FP per target,
  gate; details folded); the answer set in one line (composition folded); the evaluator and its nodes; warnings
  grouped by message.
- **History**: one row per run, in sub-tabs Overview (true score vs the field's best and median, rank %, local - true),
  Answer set (accuracy vs the released truth, zoomed), Evaluator (configs per run, per minute, per GATK run) and
  Search (gain over baseline, pick gap, Spearman local ~ true).

- **Address**: `https://<main's IP>:8443` (the bare `http://<IP>` redirects there). Real Let's Encrypt certificate for the bare IP
  (no domain needed). A self-signed certificate is not enough: browsers then fail to load Streamlit's lazy JS chunks
  ("Failed to fetch dynamically imported module").
- **Access**: public, no password (read-only; node names anonymized).
- **Data**: main writes `<DATA_DIR>/status/live.json` and `rounds.jsonl` (`main/status.py`, `main/review.py`); the
  dashboard only reads them, so it can never slow down or break mining. Node names are anonymized (public page).
- **When the TRUE columns fill in**: once `main.review` has reviewed the run - about one round after it for the
  review, a few rounds for the rank; a replay is reviewed right after mining.
- **Service**: `systemctl status|restart minos-dashboard` (unit in `/etc/systemd/system/minos-dashboard.service`,
  runs `.venv/bin/streamlit run app.py` from this folder; config in `.streamlit/config.toml`).
- **Install** (main only, not the fleet nodes): `.venv/bin/pip install -r dashboard/requirements.txt`;
  certificate: see below.

## Port 80 and the certificate

- `redirect.py` (service `minos-dashboard-redirect`, plain python3) answers on port 80: it serves Let's Encrypt
  challenge files from `.certs/webroot/` and redirects everything else to `https://`.
- `renew_cert.sh` gets/renews a short-lived (~6 day) Let's Encrypt certificate for the IP with
  [lego](https://github.com/go-acme/lego) (`/opt/lego/lego`, account/state in `/opt/lego/data`) and installs it as
  `.certs/cert.pem` / `.certs/key.pem`, restarting the dashboard when it changed. The timer `minos-dashboard-cert.timer`
  runs it twice a day (`minos-dashboard-cert.service`); lego only renews at half-life.
- New machine: download the lego release binary to `/opt/lego/lego`, start the redirect service, run
  `DASHBOARD_IP=<IP> dashboard/renew_cert.sh` once, then enable the timer.
