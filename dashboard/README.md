# dashboard

The live dashboard on main: what the miner is doing right now, and an editor for the search's settings. It is
behind a login. The History tab lists finished runs from `status/rounds.jsonl`.

- **Address**: `https://<main's IP>:8443` (the bare `http://<IP>` redirects there).
- **First use**: set the password on main, then log in as `admin`:

  ```bash
  .venv/bin/python -m dashboard passwd            # asked twice, never echoed
  ```

## What it shows

- **Live tab**: round, phase, time left, best local score, submissions, configs scored, busy workers, and the answer set. The schedule is the current engine: opening, six climbing basins, no new basin in the last minutes, then one merge at the final. A domain card shows the score it opened at and its first scan. The scout sampler is shown as off. Scores, fleet nodes, and the event feed are unchanged.
- **Themes**: "clean white" for light and "midnight" for dark. The page follows the system setting until the
  Light/Dark button in the top bar is used; the choice is remembered in that browser.
- **History tab**: one block per finished run, newest first. Select a run and three groups open under it. Score: TRUE, local at submit, rank, field best, gap. Search: GATK workers (saved from the next round on; older runs have none), configs evaluated, how long the search ran. Answer set: donor, its score against the released truth (the same core the live Result panel shows), germline exact match, N, build time.
- **Settings tab**: every setting in `common/tuning.py` by group, with the value in force on main, the default and
  your value. Live settings apply on main within about a second; round settings wait for the next round. Values are
  checked before the file is written, and main checks them again; a rejected value keeps the one in force.

## How it works

- `server.py` (FastAPI + uvicorn) reads `status/live.json`, which main rewrites every few seconds (`main/status.py`,
  schema 2), and reads/writes `tuning.json` in the repo folder (git-ignored), which main re-reads when it changes. It
  never talks to main directly, so the dashboard can never slow down or break mining.
- `static/`: one page in plain JavaScript (no build step) that polls `/api/status` every 2 s; uPlot (MIT, in
  `static/vendor/`) draws the score chart. Nothing loads from the internet.
- **Login** (`auth.py`): one account in `~/.minos/dashboard/auth.json` (mode 600, outside the repo): a scrypt hash
  and a signing key. A session cookie lasts 12 hours (HttpOnly, Secure, SameSite=Strict). Five failed logins from one
  address lock it out for 15 minutes. A new password ends every session.
- **Hardening**: every page but the login needs the session; writes also need a same-origin request from the page
  itself (`X-Minos` header and Origin check); a strict Content-Security-Policy (no inline scripts or styles, no
  framing); no API docs routes.
- **Service**: `minos-dashboard` (`/etc/systemd/system/minos-dashboard.service`) runs
  `.venv/bin/python -m dashboard serve` from the repo folder. `systemctl status|restart minos-dashboard`.
- **Development**: `.venv/bin/python -m dashboard serve --local --port 8080` serves plain HTTP on 127.0.0.1 only;
  reach it through an SSH tunnel. `tests/test_dashboard.py` covers login, lockout, sessions and the settings checks.

## Port 80 and the certificate

- `redirect.py` (service `minos-dashboard-redirect`, plain python3) answers on port 80: it serves Let's Encrypt
  challenge files from `.certs/webroot/` and redirects everything else to `https://<IP>:8443`.
- `renew_cert.sh` gets/renews a short-lived (~6 day) Let's Encrypt certificate for the IP with
  [lego](https://github.com/go-acme/lego) (`/opt/lego/lego`, account/state in `/opt/lego/data`) and installs it as
  `.certs/cert.pem` / `.certs/key.pem`, restarting the dashboard if it is installed. The timer
  `minos-dashboard-cert.timer` runs it twice a day (`minos-dashboard-cert.service`); lego only renews at half-life.
- A real certificate for the bare IP is needed: a login page must not run on a self-signed one.
- New machine: download the lego release binary to `/opt/lego/lego`, start the redirect service, run
  `DASHBOARD_IP=<IP> dashboard/renew_cert.sh` once, then enable the timer, install the service below and set the
  password.

```ini
[Unit]
Description=Minos live dashboard (HTTPS on 8443, login)
After=network-online.target

[Service]
WorkingDirectory=/minos/minos-miner
Environment=DATA_DIR=/minos/data
ExecStart=/minos/minos-miner/.venv/bin/python -m dashboard serve
Restart=always
RestartSec=5
TimeoutStopSec=5

[Install]
WantedBy=multi-user.target
```
