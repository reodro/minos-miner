"""Minos miner dashboard (read-only). Reads what main writes (main/status.py); never talks to the miner.

    run:  see dashboard/README.md (systemd service minos-dashboard, HTTPS on port 443)
    password: dashboard/.streamlit/secrets.toml  (password = "...", git-ignored)
"""
from __future__ import annotations

import hmac
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import streamlit as st

DATA_DIR = Path(os.environ.get("DATA_DIR", "/minos/data"))
STATUS = DATA_DIR / "status"

st.set_page_config(page_title="Minos miner", page_icon="🧬", layout="wide")


# ---- login (splash page) -----------------------------------------------------------------------------------------
def logged_in() -> bool:
    if st.session_state.get("ok"):
        return True
    st.title("🧬 Minos miner")
    with st.form("login"):
        pw = st.text_input("Password", type="password")
        if st.form_submit_button("Enter"):
            if hmac.compare_digest(pw.encode(), str(st.secrets.get("password", "")).encode()):
                st.session_state["ok"] = True
                st.rerun()
            else:
                time.sleep(1.0)                     # slows down guessing
                st.error("Wrong password")
    return False


def ago(ts) -> str:
    if not ts:
        return "-"
    s = int(time.time() - ts)
    return f"{s}s ago" if s < 120 else f"{s // 60} min ago"


def mmss(sec) -> str:
    sec = int(max(0, sec))
    return f"{sec // 60}:{sec % 60:02d}"


# ---- live round ----------------------------------------------------------------------------------------------------
@st.fragment(run_every="10s")
def live():
    path = STATUS / "live.json"
    if not path.exists():
        st.info("No round has run yet.")
        return
    s = json.loads(path.read_text())
    r, a, srch, ex = s.get("round", {}), s.get("answer_set") or {}, s.get("search", {}), s.get("executor") or {}
    started = r.get("started_ts") or time.time()

    c1, c2, c3 = st.columns([2, 1, 1])
    c1.subheader(f"Round {r.get('tag', '-')}  ·  {r.get('region', '')}")
    c2.metric("Status", r.get("mode", "-"))
    left = (r["deadline_ts"] - time.time()) if r.get("deadline_ts") else None
    c3.metric("Time left", mmss(left) if left is not None else "-")
    st.caption(f"updated {ago(s.get('updated'))}")

    ev = srch.get("eval_timeline") or []
    recent = [e for e in ev if e[0] > time.time() - 300]
    rate = (recent[-1][1] - recent[0][1]) / max((recent[-1][0] - recent[0][0]) / 60, 1e-9) if len(recent) > 1 else 0
    m = st.columns(5)
    m[0].metric("Configs evaluated", srch.get("evaluated", 0))
    m[1].metric("Per minute (last 5 min)", f"{rate:.0f}")
    m[2].metric("Best local score", f"{srch['best']:.3f}" if srch.get("best") is not None else "-")
    m[3].metric("Submissions", len(s.get("submissions") or []))
    m[4].metric("Failed / in flight", f"{srch.get('failed', 0)} / {srch.get('in_flight', 0)}")

    left_col, right_col = st.columns(2)
    with left_col:
        st.markdown("**Answer set**")
        if a:
            secs = a.get("seconds") or {}
            st.table(pd.DataFrame({"": ["regime / donor", "germline records", "synthetic calls", "N (fp denominator)",
                                        "build time"],
                                   "value": [f"{a.get('regime')} / {a.get('donor')}", str(a.get("germline")),
                                             str(a.get("synthetic")), f"{a.get('n_targets')} ({a.get('n_targets_source')})",
                                             f"{sum(secs.values()):.0f} s" if secs else "-"]}).set_index(""))
            if a.get("degraded"):
                st.warning(a["degraded"])
        else:
            st.write("building...")
    with right_col:
        st.markdown(f"**Evaluator** ({ex.get('kind', '-')})")
        if ex:
            st.write(f"workers busy {ex.get('running', 0)} / {ex.get('workers', 0)} · queued {ex.get('queued', 0)}")
            if ex.get("nodes"):
                st.dataframe(pd.DataFrame(ex["nodes"]), hide_index=True, use_container_width=True)

    g1, g2 = st.columns(2)
    best = srch.get("best_timeline") or []
    if best:
        g1.markdown("**Best local score over time** (minutes into the round)")
        g1.line_chart(pd.DataFrame({"minute": [(t - started) / 60 for t, _ in best], "best": [v for _, v in best]})
                      .set_index("minute"))
    if ev:
        g2.markdown("**Configs evaluated over time**")
        g2.line_chart(pd.DataFrame({"minute": [(t - started) / 60 for t, _ in ev], "evaluated": [v for _, v in ev]})
                      .set_index("minute"))

    subs = s.get("submissions") or []
    if subs:
        st.markdown("**Submissions** (the last one before the deadline counts)")
        st.dataframe(pd.DataFrame([{"minute": round((t - started) / 60, 1), "config": k, "local score": round(v, 4)}
                                   for t, k, v in subs]), hide_index=True, use_container_width=True)
    logs = s.get("log") or []
    if logs:
        with st.expander(f"Warnings and errors ({len(logs)})"):
            st.dataframe(pd.DataFrame([{"time": datetime.fromtimestamp(t, timezone.utc).strftime("%H:%M:%S"),
                                        "level": lv, "where": w, "message": msg} for t, lv, w, msg in reversed(logs)]),
                         hide_index=True, use_container_width=True)


# ---- history -----------------------------------------------------------------------------------------------------------
def history():
    path = STATUS / "rounds.jsonl"
    if not path.exists():
        st.info("No finished rounds yet.")
        return
    rounds, truth = {}, {}
    for line in path.read_text().splitlines():
        try:
            row = json.loads(line)
        except json.JSONDecodeError:
            continue
        (truth if row.get("type") == "true" else rounds)[row["round_id"]] = row
    rows = []
    for rid, r in sorted(rounds.items(), key=lambda kv: kv[1].get("ended") or 0, reverse=True):
        t = truth.get(rid, {})
        rows.append({"finished (UTC)": datetime.fromtimestamp(r["ended"], timezone.utc).strftime("%m-%d %H:%M"),
                     "round": r.get("tag"), "region": r.get("region"), "donor": f"{r.get('regime')}/{r.get('donor')}",
                     "mode": r.get("mode"), "configs": r.get("evaluated"), "best local": r.get("best_local"),
                     "true score": t.get("true_score"),
                     "rank": f"{t['rank']}/{t['field']}" if t.get("rank") else None, "note": r.get("note") or ""})
    df = pd.DataFrame(rows)
    st.dataframe(df, hide_index=True, use_container_width=True)
    ranked = [r for r in rows if r["true score"] is not None]
    if ranked:
        st.markdown("**True score of the final submission** (rounds with a known true score)")
        st.line_chart(pd.DataFrame(ranked[::-1]).set_index("round")[["true score"]])


if logged_in():
    top = st.columns([6, 1])
    top[0].title("🧬 Minos miner")
    if top[1].button("Log out"):
        st.session_state.clear()
        st.rerun()
    tab_live, tab_hist = st.tabs(["Live round", "History"])
    with tab_live:
        live()
    with tab_hist:
        history()
