"""Port 80: serves Let's Encrypt http-01 challenges (dashboard/renew_cert.sh) and redirects everything else to the
dashboard (https on DASHBOARD_PORT), so typing the bare address (5.196.64.61) still opens it."""
import os
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

WEBROOT = Path(__file__).resolve().parent / ".certs" / "webroot"
CHALLENGE = "/.well-known/acme-challenge/"
PORT = int(os.environ.get("DASHBOARD_PORT", "8443"))


class Redirect(BaseHTTPRequestHandler):
    def do_GET(self):
        if self.path.startswith(CHALLENGE):
            f = WEBROOT / ".well-known" / "acme-challenge" / Path(self.path[len(CHALLENGE):]).name
            if f.is_file():
                body = f.read_bytes()
                self.send_response(200)
                self.send_header("Content-Type", "text/plain")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
                return
            self.send_error(404)
            return
        host = (self.headers.get("Host") or "").split(":")[0]
        self.send_response(301)
        self.send_header("Location", f"https://{host}:{PORT}{self.path}")
        self.end_headers()

    do_HEAD = do_GET

    def log_message(self, *args):
        pass


if __name__ == "__main__":
    ThreadingHTTPServer(("0.0.0.0", 80), Redirect).serve_forever()
