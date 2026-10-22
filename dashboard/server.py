"""The live dashboard's web server (FastAPI). It only reads main's status file and reads/writes the settings file,
so it can never slow down or break mining; main applies settings itself (common/tuning.py).

  GET  /                 the dashboard (login first)
  GET  /login            the login page
  POST /api/login        {"user", "password"} -> a session cookie; 5 failures from one address lock it for 15 min
  POST /api/logout
  GET  /api/status       main's live.json as it is, the server's clock, and the run's review/rank rows when they exist
  GET  /api/history      one row per finished run from rounds.jsonl (round, review, rank), newest first
  GET  /api/points       every config scored this round, projected to 3D (dashboard/points.py)
  GET  /api/curves       per-setting response curves over the stored rounds (test/bench/param_curves.py)
  GET  /api/settings     every setting's metadata and default, the file's overrides and what main would reject
  PUT  /api/settings     {"set": {name: value}, "reset": [name]} -> checked, then the file is replaced atomically

Every request but the login page and static files needs the session cookie. Writes also need the header
`X-Minos: 1` and a same-origin Origin, so another site cannot make a logged-in browser change settings.
"""
from __future__ import annotations

import asyncio
import json
import re
import sqlite3
import time
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, Request
from fastapi.middleware.gzip import GZipMiddleware
from fastapi.responses import FileResponse, JSONResponse, RedirectResponse, Response
from fastapi.staticfiles import StaticFiles

from common.settings import SETTINGS
from common.tuning import DEFAULTS, GROUPS, TUNING_FILE, describe, read_overrides, validate, write_overrides
from dashboard.auth import AUTH_FILE, SESSION_S, Auth
from dashboard.field import FieldCache, live_rank, published_scores, refresh_rank

HERE = Path(__file__).resolve().parent
STATIC = HERE / "static"
RERANK_S = 24 * 3600          # runs that ended within this long are placed again on the published field
CLAMPED = STATIC / "clamped_true.json"   # test/bench/clamped_true.py: final configs scored as the validators run them


def clamped_true() -> dict:
    """run id -> the final config's TRUE score after clamping to the platform's ranges (2026-10-18: the platform
    clamps a submission and the validators score the clamped config). Empty when the file is missing or unreadable."""
    try:
        return json.loads(CLAMPED.read_text())
    except (OSError, ValueError):
        return {}


def apply_clamped(row: dict, table: dict) -> dict:
    hit = table.get(row.get("run_id") or "")
    if hit and isinstance(hit.get("true_score"), (int, float)):
        row = {**row, "unclamped_true": row.get("true_score"), "true_score": hit["true_score"],
               "clamped": hit.get("clamped") or []}
    return row
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


def _core_timeline(evaluator: dict) -> list:
    """Busy-core samples for the history chart: [minutes into the round, busy, capacity]."""
    raw = evaluator.get("core_timeline")
    if not isinstance(raw, list):
        return []
    out = []
    for row in raw:
        if not isinstance(row, (list, tuple)) or len(row) < 2:
            continue
        try:
            t, busy = float(row[0]), int(row[1])
            cap = int(row[2]) if len(row) > 2 and row[2] is not None else None
        except (TypeError, ValueError):
            continue
        out.append([t, busy, cap])
    return out


SIMILAR_SETTINGS = 4          # same hill as search/climb.py: fewer than this many setting differences


def _cfg_params(raw) -> Optional[dict]:
    try:
        d = json.loads(raw) if isinstance(raw, str) else raw
    except (TypeError, ValueError):
        return None
    if not isinstance(d, dict):
        return None
    p = d.get("params")
    return p if isinstance(p, dict) else d


def _same_hill(a: Optional[dict], b: Optional[dict]) -> bool:
    if not a or not b:
        return False
    names = set(a) | set(b)
    return sum(a.get(n) != b.get(n) for n in names) < SIMILAR_SETTINGS


def attach_climber_found(live: Optional[dict], cache: dict) -> None:
    """Each climber's best returned score (its own B{id} runs), kept even before the centre moves there.

    Also the highest Step B GA (tag G) scores that are a different hill from every seated climber
    and from each other (fewer than four setting differences is the same hill), so the dashboard
    can show them while the live search process only publishes how many GA children are in flight.

    Only this run's rows count, from the engine's started_ts on: the results store keeps every replay of a round,
    and the same B{id} tags come back each run (044200Z, 2026-10-16: the live panel showed the last run's scores).
    """
    try:
        eng = live["search"]["engine"]
        basins = eng["basins"]
        tag = live["round"]["tag"]
    except (AttributeError, KeyError, TypeError):
        return
    if not isinstance(basins, list) or not tag:
        return
    path = Path(SETTINGS.cache_dir) / "results" / f"{tag}.sqlite"
    started = eng.get("started_ts")
    since = float(started) - 5.0 if isinstance(started, (int, float)) and not isinstance(started, bool) else None
    where, args = (" and created >= ?", (since,)) if since is not None else ("", ())
    now = time.time()
    if cache.get("path") != str(path) or cache.get("since") != since or now - cache.get("at", 0) > 2:
        by_tag: dict = {}
        ga_ranked: list = []
        params_by_key: dict = {}
        climber_keys = [b.get("best_key") for b in basins if isinstance(b, dict) and isinstance(b.get("best_key"), str)]
        try:
            con = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
            try:
                for row_tag, score in con.execute(
                        "select tag, max(score) from results where score is not null and tag like 'B%'" + where
                        + " group by tag", args):
                    if isinstance(score, (int, float)):
                        by_tag[row_tag] = float(score)
                for key, score in con.execute(
                        "select config_key, max(score) from results where tag='G' and score is not null" + where
                        + " group by config_key order by max(score) desc limit 80", args):
                    if key and isinstance(score, (int, float)):
                        ga_ranked.append({"key": key, "score": round(float(score), 4)})
                need = [row["key"] for row in ga_ranked] + climber_keys
                if need:
                    q = "select config_key, config from results where config_key in ({})".format(
                        ",".join("?" * len(need)))
                    for key, raw in con.execute(q, need):
                        if key in params_by_key:
                            continue
                        p = _cfg_params(raw)
                        if p is not None:
                            params_by_key[key] = p
            finally:
                con.close()
        except (OSError, sqlite3.Error):
            by_tag = cache.get("by_tag") or {}
            ga_ranked = cache.get("ga_ranked") or []
            params_by_key = cache.get("params_by_key") or {}
        cache["path"], cache["at"], cache["by_tag"], cache["ga_ranked"], cache["params_by_key"] = (
            str(path), now, by_tag, ga_ranked, params_by_key)
        cache["since"] = since
    by_tag = cache.get("by_tag") or {}
    params_by_key = cache.get("params_by_key") or {}
    held_keys, held_scores = set(), set()
    hills = []
    for basin in basins:
        if not isinstance(basin, dict):
            continue
        centre = basin.get("best")
        found = by_tag.get(f"B{basin.get('id')}")
        vals = [v for v in (found, centre) if isinstance(v, (int, float)) and not isinstance(v, bool)]
        if vals:
            basin["found"] = round(max(vals), 4)
        key = basin.get("best_key")
        if isinstance(key, str) and key:
            held_keys.add(key)
            p = params_by_key.get(key)
            if p is not None:
                hills.append(p)
        for v in (basin.get("found"), centre):
            if isinstance(v, (int, float)) and not isinstance(v, bool):
                held_scores.add(round(float(v), 4))
    ga = eng.get("ga")
    if not isinstance(ga, dict):
        ga = {}
        eng["ga"] = ga
    top = []
    for row in cache.get("ga_ranked") or []:
        if row["key"] in held_keys or row["score"] in held_scores:
            continue
        p = params_by_key.get(row["key"])
        if any(_same_hill(p, h) for h in hills):
            continue
        top.append({"key": row["key"], "score": row["score"]})
        if p is not None:
            hills.append(p)
        if len(top) == 3:
            break
    ga["top"] = top


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
            "round_id": _text(rank, "round_id") or _text(rnd, "round_id"),
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
            "core_timeline": _core_timeline((rnd or {}).get("evaluator") or {}),
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
               secure_cookie: bool = True, fail_delay_s: float = 1.0,
               field_cache: Optional[FieldCache] = None) -> FastAPI:
    status_file = Path(status_dir or SETTINGS.status_dir) / "live.json"
    history_file = status_file.with_name("rounds.jsonl")
    tuning_path = Path(tuning_file or TUNING_FILE)
    auth = Auth(Path(auth_file or AUTH_FILE))
    ranks = field_cache if field_cache is not None else FieldCache(published_scores)
    app = FastAPI(docs_url=None, redoc_url=None, openapi_url=None)
    app.mount("/static", StaticFiles(directory=STATIC), name="static")
    app.add_middleware(GZipMiddleware, minimum_size=50_000)      # /api/points carries every scored config
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
        # never cached: a cached page kept an old layout after an update (the island-mode Step A panel, 2026-10-21)
        return FileResponse(STATIC / "index.html", headers={"Cache-Control": "no-store"})

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
        attach_climber_found(live, found_cache)
        rank = await asyncio.to_thread(live_rank, live, ranks)
        result = result_of(live)
        if result and result.get("rank"):
            kind = ((live or {}).get("round") or {}).get("kind")
            row = apply_clamped(result["rank"], clamped_true())
            result["rank"] = await asyncio.to_thread(refresh_rank, row, ranks, replay=kind == "replay")
        return {"now": time.time(), "live": live, "result": result, "rank": rank}

    @app.get("/api/points")
    async def all_points(request: Request):
        if user_of(request) is None:
            return denied()
        try:
            live = json.loads(status_file.read_text())
        except (OSError, ValueError):
            live = None
        tag = ((live or {}).get("round") or {}).get("tag")
        if not tag:
            return {"round": None, "n": 0}
        from dashboard.points import points
        return await asyncio.to_thread(points, tag)

    @app.get("/api/curves")
    async def curves(request: Request):
        if user_of(request) is None:
            return denied()
        runs = HERE.parent / "test" / "runs"
        try:
            data = json.loads((runs / "param_curves.json").read_text())
            summary = json.loads((runs / "param_curves_summary.json").read_text())
        except (OSError, ValueError):
            return {"bins": {}, "rounds": {}, "summary": []}
        return {"bins": data.get("bins", {}), "rounds": data.get("rounds", {}), "summary": summary}

    @app.get("/api/history")
    async def history(request: Request):
        if user_of(request) is None:
            return denied()
        runs = history_of(history_file)

        def fresh() -> list:                      # recent runs' ranks on the field as published now (refresh_rank)
            cutoff, table = time.time() - RERANK_S, clamped_true()
            for r in runs:
                fixed = r.get("run_id") in table
                if fixed:
                    r.update(apply_clamped(r, table))
                if r.get("true_score") is not None and (fixed or (r.get("ended") or 0) > cutoff):
                    r.update(refresh_rank(r, ranks, replay=r.get("kind") == "replay"))
            return runs
        return {"runs": await asyncio.to_thread(fresh)}

    _wave_at = re.compile(r"engine: wave (\d+) \([^)]*\) at t\+([0-9.]+)m")
    _dom_at = re.compile(r"domain (\d+) (?:-> [\d.]+|born at [\d.]+).*t\+([0-9.]+)m")
    _wave_cache: dict = {}
    found_cache: dict = {}

    def fill_wave_starts(live: Optional[dict]) -> None:
        """An older search process does not publish wave starts or when a domain's best last moved. Take both from its log."""
        try:
            eng = live["search"]["engine"]
            waves = eng.get("waves") or []
            domains = eng.get("domains") or []
            tag = (live.get("round") or {}).get("tag") or ""
        except (AttributeError, KeyError, TypeError):
            return
        need_waves = any(w.get("started_s") is None for w in waves)
        need_domains = any(d.get("best_s") is None for d in domains)
        if not tag or not (need_waves or need_domains):
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
            wave_at = {int(m.group(1)): round(float(m.group(2)) * 60, 1) for m in _wave_at.finditer(text)}
            dom_at: dict = {}
            for m in _dom_at.finditer(text):
                dom_at[int(m.group(1))] = round(float(m.group(2)) * 60, 1)
            found = {"waves": wave_at, "domains": dom_at}
            _wave_cache.clear()
            _wave_cache[key] = found
        for w in waves:
            if w.get("started_s") is None and w.get("wave") in found["waves"]:
                w["started_s"] = found["waves"][w["wave"]]
        for d in domains:
            if d.get("best_s") is None and d.get("id") in found["domains"]:
                d["best_s"] = found["domains"][d["id"]]

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
