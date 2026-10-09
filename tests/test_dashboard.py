"""The dashboard server: login and lockout, the session cookie, status, and the settings editor's checks."""
import json
import tempfile
import unittest
from pathlib import Path

from fastapi.testclient import TestClient

from common.tuning import read_overrides
from dashboard.field import FieldCache, place
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
        self.assertEqual(self.c.get("/api/history").status_code, 401)
        self.assertEqual(self.c.get("/api/settings").status_code, 401)
        self.assertEqual(self.c.put("/api/settings", json={"set": {"candidates.mutate_frac": 0.2}}, headers=W).status_code, 401)
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
        r = self.c.put("/api/settings", json={"set": {"candidates.mutate_frac": 0.2, "candidates.workers_per_candidate": 100}}, headers=W)
        self.assertEqual(r.status_code, 200)
        self.assertEqual(read_overrides(self.tuning)[0], {"candidates.workers_per_candidate": 100, "candidates.mutate_frac": 0.2})
        r = self.c.put("/api/settings", json={"set": {"candidates.mutate_frac": 5, "candidates.core_band": 0.1}}, headers=W)
        self.assertEqual(r.status_code, 422)
        self.assertEqual([x["name"] for x in r.json()["rejected"]], ["candidates.mutate_frac"])
        self.assertEqual(read_overrides(self.tuning)[0], {"candidates.workers_per_candidate": 100, "candidates.mutate_frac": 0.2})   # unchanged
        r = self.c.put("/api/settings", json={"reset": ["candidates.mutate_frac"], "set": {"candidates.workers_per_candidate": 300}}, headers=W)
        self.assertEqual(read_overrides(self.tuning)[0], {})                # a default needs no override
        self.assertEqual(r.json()["overrides"], {})
        self.assertEqual(self.c.put("/api/settings", json={"set": {"candidates.core_band": 0.1}}).status_code, 400)
        self.assertEqual(self.c.put("/api/settings", json={"set": {"candidates.core_band": 0.1}},
                                    headers={**W, "Origin": "https://evil.example"}).status_code, 400)
        meta = self.c.get("/api/settings").json()
        self.assertEqual(len(meta["knobs"]), len(meta["defaults"]))

    def test_history_is_one_row_per_run_newest_first(self):
        set_password("admin", PASSWORD, self.auth)
        self.login()
        rows = [
            {"type": "round", "run_id": "old", "tag": "20261002T092900Z", "region": "chr18:1-2", "kind": "replay",
             "started": 100.0, "ended": 200.0,
             "answer_set": {"donor": "HG001", "regime": "giab"},
             "search": {"best_local": 88.1}, "evaluator": {"evaluated": 10, "failed": 1, "workers": 162}},
            {"type": "review", "run_id": "old", "tag": "20261002T092900Z",
             "answer_set": {"germline_exact": {"match": 9, "truth": 10}},
             "search": {"final_true": 89.1, "final_local": 88.1, "baseline_true": 86.0, "spearman": 0.99,
                        "pick_gap": 0.01, "gain_over_baseline": 3.1, "configs_noise": [1, 2, 3]}},
            {"type": "rank", "run_id": "old", "tag": "20261002T092900Z", "true_score": 89.1, "rank": 23, "field": 66,
             "field_best": 89.9, "field_median": 89.4, "gap_to_best": 0.8},
            {"type": "round", "run_id": "new", "tag": "20261007T070500Z", "region": "chr20:3-4", "kind": "replay",
             "started": 300.0, "ended": 400.0, "search": {"best_local": 89.7}, "evaluator": {"evaluated": 5, "failed": 0}},
            {"type": "review", "run_id": "old",
             "answer_set": {"germline_exact": {"match": 9, "truth": 10},
                            "as_scored": {"core": 0.9912}},
             "search": {"final_true": 89.2, "spearman": 0.5, "pick_gap": 0.0}},
            "not json",
            {"type": "note", "run_id": "old"},
        ]
        (self.status / "rounds.jsonl").write_text("".join(json.dumps(r) + "\n" if isinstance(r, dict) else r + "\n" for r in rows))
        body = self.c.get("/api/history").json()
        self.assertEqual([r["run_id"] for r in body["runs"]], ["new", "old"])
        old = body["runs"][1]
        self.assertEqual(old["final_true"], 89.2)          # the later review line replaces the earlier one
        self.assertEqual(old["rank"], 23)
        self.assertEqual(old["gap_to_best"], 0.8)
        self.assertEqual(old["germline_match"], 9)
        self.assertEqual(old["workers"], 162)
        self.assertEqual(old["answer_score"], 0.9912)
        self.assertNotIn("configs_noise", old)
        new = body["runs"][0]
        self.assertFalse(new["reviewed"])
        self.assertFalse(new["ranked"])
        self.assertEqual(new["best_local"], 89.7)
        self.assertIsNone(new["true_score"])

    def test_a_replay_shows_its_live_rank_on_the_published_field(self):
        set_password("admin", PASSWORD, self.auth)
        self.login()
        (self.status / "live.json").write_text(json.dumps({
            "schema": 2,
            "round": {"tag": "r1", "id": "2026-10-09T01:06:00+00:00", "kind": "replay"},
            "search": {"best": 85.05},
        }))
        cache = FieldCache(lambda rid: [85.30, 85.07, 84.00], ttl_s=60)
        self.c = TestClient(create_app(self.status, self.tuning, self.auth, secure_cookie=False,
                                       fail_delay_s=0.0, field_cache=cache))
        self.login()
        rank = self.c.get("/api/status").json()["rank"]
        self.assertEqual(rank["rank"], 3)
        self.assertEqual(rank["field"], 4)
        self.assertEqual(rank["field_best"], 85.3)
        self.assertEqual(rank["gap_to_best"], 0.25)
        self.assertTrue(rank["filling"])

    def test_a_live_round_has_no_field_rank(self):
        set_password("admin", PASSWORD, self.auth)
        (self.status / "live.json").write_text(json.dumps({
            "schema": 2, "round": {"tag": "r1", "id": "round", "kind": "live"}, "search": {"best": 85.05},
        }))
        called = []
        cache = FieldCache(lambda rid: called.append(rid) or [90.0], ttl_s=60)
        self.c = TestClient(create_app(self.status, self.tuning, self.auth, secure_cookie=False,
                                       fail_delay_s=0.0, field_cache=cache))
        self.login()
        self.assertIsNone(self.c.get("/api/status").json()["rank"])
        self.assertEqual(called, [])

    def test_place_counts_a_replay_as_one_more_in_the_field(self):
        row = place(89.0, [89.3, 89.1, 88.0], replay=True)
        self.assertEqual((row["rank"], row["field"], row["gap_to_best"]), (3, 4, 0.3))

    def test_password_hashes(self):
        h = hash_password("x" * 12)
        self.assertTrue(check_password("x" * 12, h))
        self.assertFalse(check_password("y" * 12, h))
        self.assertFalse(check_password("x" * 12, "garbage"))
        with self.assertRaises(ValueError):
            set_password("admin", "", self.auth)


if __name__ == "__main__":
    unittest.main()
