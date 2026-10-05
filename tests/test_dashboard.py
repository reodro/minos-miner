"""The dashboard server: login and lockout, the session cookie, status, and the settings editor's checks."""
import json
import tempfile
import unittest
from pathlib import Path

from fastapi.testclient import TestClient

from common.tuning import read_overrides
from dashboard.auth import MAX_FAILS, check_password, hash_password, set_password
from dashboard.server import COOKIE, create_app

W = {"X-Minos": "1"}
PASSWORD = "correct horse battery"


class DashboardTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        d = Path(self.tmp.name)
        self.status, self.tuning, self.auth = d / "status", d / "tuning.json", d / "auth" / "auth.json"
        self.status.mkdir()
        (self.status / "live.json").write_text(json.dumps({"schema": 2, "round": {"tag": "r1"}}))
        self.c = TestClient(create_app(self.status, self.tuning, self.auth, secure_cookie=False, fail_delay_s=0.0))

    def tearDown(self):
        self.tmp.cleanup()

    def login(self, password=PASSWORD, headers=W):
        return self.c.post("/api/login", json={"user": "admin", "password": password}, headers=headers)

    def test_everything_but_the_login_page_needs_a_session(self):
        r = self.c.get("/", follow_redirects=False)
        self.assertEqual((r.status_code, r.headers["location"]), (303, "/login"))
        self.assertEqual(self.c.get("/login").status_code, 200)
        self.assertEqual(self.c.get("/api/status").status_code, 401)
        self.assertEqual(self.c.get("/api/settings").status_code, 401)
        self.assertEqual(self.c.put("/api/settings", json={"set": {"scout.share": 0.2}}, headers=W).status_code, 401)
        self.assertIn("frame-ancestors 'none'", self.c.get("/login").headers["content-security-policy"])

    def test_no_account_no_login(self):
        self.assertEqual(self.login().status_code, 503)

    def test_login_status_and_logout(self):
        set_password("admin", PASSWORD, self.auth)
        self.assertEqual(oct(self.auth.stat().st_mode & 0o777), "0o600")
        self.assertEqual(self.login(headers={}).status_code, 400)          # no X-Minos header: not from our page
        r = self.login()
        self.assertEqual(r.status_code, 200)
        self.assertIn(COOKIE, r.cookies)
        self.assertEqual(self.c.get("/api/status").json()["live"]["round"]["tag"], "r1")
        self.assertEqual(self.c.get("/", follow_redirects=False).status_code, 200)
        self.c.post("/api/logout", headers=W)
        self.assertEqual(self.c.get("/api/status").status_code, 401)

    def test_a_new_password_ends_every_session(self):
        set_password("admin", PASSWORD, self.auth)
        self.login()
        set_password("admin", "another long password", self.auth)
        self.assertEqual(self.c.get("/api/status").status_code, 401)

    def test_failed_logins_lock_the_address_out(self):
        set_password("admin", PASSWORD, self.auth)
        for _ in range(MAX_FAILS):
            self.assertEqual(self.login("wrong password!!").status_code, 401)
        self.assertEqual(self.login().status_code, 429)                     # even the right password, for a while

    def test_settings_are_checked_before_the_file_is_written(self):
        set_password("admin", PASSWORD, self.auth)
        self.login()
        r = self.c.put("/api/settings", json={"set": {"scout.share": 0.2, "domains.max": 4}}, headers=W)
        self.assertEqual(r.status_code, 200)
        self.assertEqual(read_overrides(self.tuning)[0], {"domains.max": 4, "scout.share": 0.2})
        r = self.c.put("/api/settings", json={"set": {"scout.share": 5, "pick.band": 0.1}}, headers=W)
        self.assertEqual(r.status_code, 422)
        self.assertEqual([x["name"] for x in r.json()["rejected"]], ["scout.share"])
        self.assertEqual(read_overrides(self.tuning)[0], {"domains.max": 4, "scout.share": 0.2})   # unchanged
        r = self.c.put("/api/settings", json={"reset": ["scout.share"], "set": {"domains.max": 6}}, headers=W)
        self.assertEqual(read_overrides(self.tuning)[0], {})                # a default needs no override
        self.assertEqual(r.json()["overrides"], {})
        self.assertEqual(self.c.put("/api/settings", json={"set": {"pick.band": 0.1}}).status_code, 400)
        self.assertEqual(self.c.put("/api/settings", json={"set": {"pick.band": 0.1}},
                                    headers={**W, "Origin": "https://evil.example"}).status_code, 400)
        meta = self.c.get("/api/settings").json()
        self.assertEqual(len(meta["knobs"]), len(meta["defaults"]))

    def test_password_hashes(self):
        h = hash_password("x" * 12)
        self.assertTrue(check_password("x" * 12, h))
        self.assertFalse(check_password("y" * 12, h))
        self.assertFalse(check_password("x" * 12, "garbage"))
        with self.assertRaises(ValueError):
            set_password("admin", "", self.auth)


if __name__ == "__main__":
    unittest.main()
