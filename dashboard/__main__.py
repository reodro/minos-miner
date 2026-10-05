"""python -m dashboard serve [--port 8443]        the live dashboard over HTTPS (the certificate in dashboard/.certs)
python -m dashboard serve --local --port 8080    plain HTTP on 127.0.0.1 only (development, or an SSH tunnel)
python -m dashboard passwd [--user admin]        set the one account's password (asked twice, never echoed)
"""
from __future__ import annotations

import argparse
import getpass
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CERT, KEY = HERE / ".certs" / "cert.pem", HERE / ".certs" / "key.pem"


def main() -> None:
    ap = argparse.ArgumentParser(prog="python -m dashboard")
    sp = ap.add_subparsers(dest="cmd", required=True)
    s = sp.add_parser("serve")
    s.add_argument("--port", type=int, default=8443)
    s.add_argument("--local", action="store_true", help="plain HTTP on 127.0.0.1 (no certificate needed)")
    p = sp.add_parser("passwd")
    p.add_argument("--user", default="admin")
    args = ap.parse_args()

    if args.cmd == "passwd":
        from dashboard.auth import AUTH_FILE, set_password
        pw = getpass.getpass(f"new password for {args.user}: ")
        if pw != getpass.getpass("again: "):
            sys.exit("the two passwords differ; nothing changed")
        try:
            set_password(args.user, pw)
        except ValueError as e:
            sys.exit(f"{e}; nothing changed")
        print(f"saved to {AUTH_FILE}; every open session has ended")
        return

    import uvicorn
    from dashboard.server import create_app
    if args.local:
        uvicorn.run(create_app(secure_cookie=False), host="127.0.0.1", port=args.port, log_level="warning")
        return
    if not (CERT.exists() and KEY.exists()):
        sys.exit(f"no certificate in {CERT.parent}: run dashboard/renew_cert.sh, or use --local")
    uvicorn.run(create_app(), host="0.0.0.0", port=args.port, ssl_certfile=str(CERT), ssl_keyfile=str(KEY),
                log_level="warning", proxy_headers=False)


if __name__ == "__main__":
    main()
