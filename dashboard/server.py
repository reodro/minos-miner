"""The live dashboard's web server (FastAPI). It only reads main's status file and reads/writes the settings file,
so it can never slow down or break mining; main applies settings itself (common/tuning.py).

  GET  /                 the dashboard (login first)
  GET  /login            the login page
  POST /api/login        {"user", "password"} -> a session cookie; 5 failures from one address lock it for 15 min
  POST /api/logout
  GET  /api/status       main's live.json as it is, the server's clock, and the run's review/rank rows when they exist
  GET  /api/history      one row per finished run from rounds.jsonl (round, review, rank), newest first
  GET  /api/settings     every setting's metadata and default, the file's overrides and what main would reject
  PUT  /api/settings     {"set": {name: value}, "reset": [name]} -> checked, then the file is replaced atomically

Every request but the login page and static files needs the session cookie. Writes also need the header
`X-Minos: 1` and a same-origin Origin, so another site cannot make a logged-in browser change settings.
"""
from __future__ import annotations

import asyncio
import json
import re
import time
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, Request
from fastapi.responses import FileResponse, JSONResponse, RedirectResponse, Response
from fastapi.staticfiles import StaticFiles

from common.settings import SETTINGS
from common.tuning import DEFAULTS, GROUPS, TUNING_FILE, describe, read_overrides, validate, write_overrides
from dashboard.auth import AUTH_FILE, SESSION_S, Auth

HERE = Path(__file__).resolve().parent
STATIC = HERE / "static"
COOKIE = "minos_session"
CSP = ("default-src 'self'; script-src 'self'; style-src 'self'; img-src 'self' data:; connect-src 'self'; "
       "font-src 'self'; frame-ancestors 'none'; base-uri 'none'; form-action 'self'")


def _num(row: Optional[dict], *keys):
    cur = row or {}
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            return None
        cur = cur[k]
    return cur if isinstance(cur, (int, float)) and not isinstance(cur, bool) else None


def _text(row: Optional[dict], *keys):
    cur = row or {}
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            return None
        cur = cur[k]
    return cur if isinstance(cur, str) else None


def history_of(path: Path) -> list:
    """One compact row per finished run. Later lines of the same type replace earlier ones. Newest run first.
    The review's per-config list is left out: the panel only needs the summary."""
    found: dict = {}
    order: list = []
    try:
        lines = path.read_text().splitlines()
    except OSError:
        return []
    for line in lines:
        try:
            row = json.loads(line)
        except ValueError:
            continue
        rid, kind = row.get("run_id"), row.get("type")
        if not rid or kind not in ("round", "review", "rank"):
            continue
        if rid not in found:
            found[rid] = {}
            order.append(rid)
        found[rid][kind] = row
    out = []
    for rid in order:
        parts = found[rid]
        rnd, rev, rank = parts.get("round"), parts.get("review"), parts.get("rank")
        search = (rev or {}).get("search") or {}
        germ = ((rev or {}).get("answer_set") or {}).get("germline_exact") or {}
        out.append({
            "run_id": rid,
            "tag": _text(rnd, "tag") or _text(rev, "tag") or _text(rank, "tag"),
            "region": _text(rnd, "region"),
            "kind": _text(rnd, "kind"),
            "mode": _text(rnd, "mode"),
            "started": _num(rnd, "started"),
            "ended": _num(rnd, "ended"),
            "donor": _text(rnd, "answer_set", "donor"),
            "regime": _text(rnd, "answer_set", "regime"),
            "best_local": _num(rnd, "search", "best_local"),
            "evaluated": _num(rnd, "evaluator", "evaluated"),
            "failed": _num(rnd, "evaluator", "failed"),
            "nodes": _num(rnd, "evaluator", "nodes"),
            "workers": _num(rnd, "evaluator", "workers"),
            "build_s": _num(rnd, "answer_set", "build_s"),
            "n_ours": _num(rnd, "answer_set", "n_targets"),
            "final_true": search.get("final_true") if isinstance(search.get("final_true"), (int, float)) else None,
            "final_local": search.get("final_local") if isinstance(search.get("final_local"), (int, float)) else None,
            "baseline_true": search.get("baseline_true") if isinstance(search.get("baseline_true"), (int, float)) else None,
            "spearman": search.get("spearman") if isinstance(search.get("spearman"), (int, float)) else None,
            "pick_gap": search.get("pick_gap") if isinstance(search.get("pick_gap"), (int, float)) else None,
            "gain_over_baseline": search.get("gain_over_baseline") if isinstance(search.get("gain_over_baseline"), (int, float)) else None,
            "germline_match": germ.get("match") if isinstance(germ.get("match"), int) else None,
            "germline_truth": germ.get("truth") if isinstance(germ.get("truth"), int) else None,
            "n_truth": _num(rev, "answer_set", "n_targets", "truth"),
            "answer_score": _num(rev, "answer_set", "as_scored", "core"),
            "true_score": _num(rank, "true_score"),
            "rank": _num(rank, "rank"),
            "field": _num(rank, "field"),
            "field_best": _num(rank, "field_best"),
            "field_median": _num(rank, "field_median"),
            "gap_to_best": _num(rank, "gap_to_best"),
            "reviewed": rev is not None,
            "ranked": rank is not None,
        })
    out.sort(key=lambda r: (r.get("ended") or r.get("started") or 0), reverse=True)
    return out


def create_app(status_dir: Optional[Path] = None, tuning_file: Optional[Path] = None, auth_file: Optional[Path] = None,
               secure_cookie: bool = True, fail_delay_s: float = 1.0) -> FastAPI:
    status_file = Path(status_dir or SETTINGS.status_dir) / "live.json"
    history_file = status_file.with_name("rounds.jsonl")
    tuning_path = Path(tuning_file or TUNING_FILE)
    auth = Auth(Path(auth_file or AUTH_FILE))
    app = FastAPI(docs_url=None, redoc_url=None, openapi_url=None)
    app.mount("/static", StaticFiles(directory=STATIC), name="static")
    write_lock = asyncio.Lock()

    def addr(request: Request) -> str:
        return request.client.host if request.client else "?"

    def user_of(request: Request) -> Optional[str]:
        return auth.user(request.cookies.get(COOKIE))

    def denied() -> JSONResponse:
        return JSONResponse({"error": "log in first"}, status_code=401)

    def same_origin_write(request: Request) -> bool:
        if request.headers.get("x-minos") != "1":
            return False
        origin = request.headers.get("origin")
        return origin is None or origin.split("://", 1)[-1] == request.headers.get("host")

    @app.middleware("http")
    async def headers(request: Request, call_next):
        resp: Response = await call_next(request)
        resp.headers["Content-Security-Policy"] = CSP
        resp.headers["X-Content-Type-Options"] = "nosniff"
        resp.headers["Referrer-Policy"] = "no-referrer"
        resp.headers["X-Frame-Options"] = "DENY"
        if secure_cookie:
            resp.headers["Strict-Transport-Security"] = "max-age=31536000"
        if request.url.path.startswith("/api/") or request.url.path in ("/", "/login"):
            resp.headers["Cache-Control"] = "no-store"
        return resp

    # ---- pages ------------------------------------------------------------------------------------------------------
    @app.get("/")
    async def index(request: Request):
        if user_of(request) is None:
            return RedirectResponse("/login", status_code=303)
        return FileResponse(STATIC / "index.html")

    @app.get("/login")
    async def login_page():
        return FileResponse(STATIC / "login.html")

    # ---- session ------------------------------------------------------------------------------------------------------
    @app.get("/api/session")
    async def session(request: Request):
        return {"user": user_of(request), "account": auth.account() is not None}

    @app.post("/api/login")
    async def login(request: Request):
        if not same_origin_write(request):
            return JSONResponse({"error": "bad request"}, status_code=400)
        a = addr(request)
        wait = auth.locked_for(a)
        if wait > 0:
            return JSONResponse({"error": f"too many failed logins: try again in {int(wait // 60) + 1} min"},
                                status_code=429)
        if auth.account() is None:
            return JSONResponse({"error": "no password is set: on main run `.venv/bin/python -m dashboard passwd`"},
                                status_code=503)
        try:
            body = await request.json()
            user, password = str(body.get("user", "")), str(body.get("password", ""))
        except (ValueError, AttributeError):
            return JSONResponse({"error": "bad request"}, status_code=400)
        token = await asyncio.to_thread(auth.login, a, user, password)     # scrypt takes ~0.1 s: off the loop
        if token is None:
            await asyncio.sleep(fail_delay_s)              # slows guessing further
            return JSONResponse({"error": "wrong user name or password"}, status_code=401)
        resp = JSONResponse({"user": user})
        resp.set_cookie(COOKIE, token, max_age=SESSION_S, httponly=True, secure=secure_cookie, samesite="strict")
        return resp

    @app.post("/api/logout")
    async def logout(request: Request):
        resp = JSONResponse({"ok": True})
        resp.delete_cookie(COOKIE, httponly=True, secure=secure_cookie, samesite="strict")
        return resp

    # ---- data -------------------------------------------------------------------------------------------------------
    @app.get("/api/status")
    async def status(request: Request):
        if user_of(request) is None:
            return denied()
        try:
            live = json.loads(status_file.read_text())
        except FileNotFoundError:
            live = None
        except (OSError, ValueError):
            live = None                                   # main is replacing it right now: the next poll gets it
        fill_wave_starts(live)
        return {"now": time.time(), "live": live, "result": result_of(live)}

    @app.get("/api/history")
    async def history(request: Request):
        if user_of(request) is None:
            return denied()
        return {"runs": history_of(history_file)}

    _wave_at = re.compile(r"engine: wave (\d+) \([^)]*\) at t\+([0-9.]+)m")
    _wave_cache: dict = {}

    def fill_wave_starts(live: Optional[dict]) -> None:
        """The running search publishes each wave's start. An older process does not, so take it from that round's log."""
        try:
            waves = live["search"]["engine"].get("waves") or []
            tag = (live.get("round") or {}).get("tag") or ""
        except (AttributeError, KeyError, TypeError):
            return
        if not waves or not tag or all(w.get("started_s") is not None for w in waves):
            return
        log_dir = Path(__file__).resolve().parents[1] / "test" / "runs"
        logs = sorted(log_dir.glob(f"replay_{tag}_*.log"))
        if not logs:
            return
        path = logs[-1]
        try:
            stamp = path.stat().st_mtime
        except OSError:
            return
        key = (str(path), stamp)
        found = _wave_cache.get(key)
        if found is None:
            try:
                text = path.read_text(errors="replace")
            except OSError:
                return
            found = {int(m.group(1)): round(float(m.group(2)) * 60, 1) for m in _wave_at.finditer(text)}
            _wave_cache.clear()
            _wave_cache[key] = found
        for w in waves:
            if w.get("started_s") is None and w.get("wave") in found:
                w["started_s"] = found[w["wave"]]

    def result_of(live: Optional[dict]) -> Optional[dict]:
        """The round, review and rank rows main appended for the run on screen (rounds.jsonl), newest of each."""
        run_id = ((live or {}).get("round") or {}).get("run_id")
        if not run_id:
            return None
        out: dict = {}
        try:
            with open(history_file) as f:
                for line in f:
                    if run_id not in line:
                        continue
                    try:
                        row = json.loads(line)
                    except ValueError:
                        continue
                    if row.get("run_id") == run_id and row.get("type") in ("round", "review", "rank"):
                        out[row["type"]] = row
        except OSError:
            return None
        return out or None

    def settings_view() -> dict:
        raw, error = read_overrides(tuning_path)
        ok, rejected = validate(raw) if not error else ({}, [])
        return {"knobs": describe(), "groups": GROUPS, "defaults": DEFAULTS, "file": str(tuning_path),
                "overrides": raw, "accepted": ok, "rejected": rejected, "error": error}

    @app.get("/api/settings")
    async def get_settings(request: Request):
        if user_of(request) is None:
            return denied()
        return settings_view()

    @app.put("/api/settings")
    async def put_settings(request: Request):
        who = user_of(request)
        if who is None:
            return denied()
        if not same_origin_write(request):
            return JSONResponse({"error": "bad request"}, status_code=400)
        try:
            body = await request.json()
            changes = dict(body.get("set") or {})
            resets = list(body.get("reset") or [])
        except (ValueError, AttributeError, TypeError):
            return JSONResponse({"error": "bad request"}, status_code=400)
        async with write_lock:
            raw, error = read_overrides(tuning_path)
            if error:
                return JSONResponse({"error": f"{error}: fix or delete the file on main first"}, status_code=409)
            new = {k: v for k, v in raw.items() if k not in resets}
            new.update(changes)
            new = {k: v for k, v in new.items() if not (k in DEFAULTS and v == DEFAULTS[k])}   # a default needs no override
            ok, rejected = validate(new)
            if rejected:
                return JSONResponse({"error": "some values were rejected; nothing was saved", "rejected": rejected},
                                    status_code=422)
            write_overrides(tuning_path, ok)
        return settings_view()

    return app
