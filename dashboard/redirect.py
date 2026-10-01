"""Port 80 -> https:// redirect, so typing the bare address (http://5.196.64.61) still opens the dashboard."""
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer


class Redirect(BaseHTTPRequestHandler):
    def do_GET(self):
        host = (self.headers.get("Host") or "").split(":")[0]
        self.send_response(301)
        self.send_header("Location", f"https://{host}{self.path}")
        self.end_headers()

    do_HEAD = do_GET

    def log_message(self, *args):
        pass


if __name__ == "__main__":
    ThreadingHTTPServer(("0.0.0.0", 80), Redirect).serve_forever()
