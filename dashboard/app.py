"""Minos miner dashboard (read-only). Reads what main writes (main/status.py); never talks to the miner.

    run:  see dashboard/README.md (systemd service minos-dashboard, HTTPS on port 8443)
"""
from __future__ import annotations

import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.settings import SETTINGS  # noqa: E402 - the repo root is on the path only after the line above

STATUS = Path(SETTINGS.status_dir)                  # the same folder main/status.py writes (STATUS_DIR overrides)

st.set_page_config(page_title="Minos miner", page_icon="🧬", layout="wide")


def ago(ts) -> str:
    if not ts:
        return "-"
    s = int(time.time() - ts)
    return f"{s}s ago" if s < 120 else f"{s // 60} min ago"


def mmss(sec) -> str:
    sec = int(max(0, sec))
    return f"{sec // 60}:{sec % 60:02d}"


def _s(v) -> str:
    return "" if v is None else str(v)


def _n(v) -> str:
    return "-" if v is None else f"{v:,}"


def scores_chart(points, best_tl, subs, started):
    """Every scored config as a dot, the best-so-far as a step line, submissions as red diamonds.
    Configs that scored 0 are left out (counted as wasted above) so the y-axis can zoom to the real scores."""
    import altair as alt
    pts = pd.DataFrame([{"minute": (t - started) / 60, "score": v, "config": k} for t, v, k in points if v > 0])
    best = pd.DataFrame([{"minute": (t - started) / 60, "best": v} for t, v in best_tl])
    sub = pd.DataFrame([{"minute": (t - started) / 60, "score": v, "config": k} for t, k, v in subs])
    layers = []
    if not pts.empty:
        layers.append(alt.Chart(pts).mark_circle(size=40, opacity=0.5).encode(
            x=alt.X("minute:Q", title="minutes into the round"),
            y=alt.Y("score:Q", scale=alt.Scale(zero=False), title="local score"),
            tooltip=["config", alt.Tooltip("score:Q", format=".3f"), alt.Tooltip("minute:Q", format=".1f")]))
    if not best.empty:
        layers.append(alt.Chart(best).mark_line(interpolate="step-after", strokeWidth=2).encode(
            x="minute:Q", y=alt.Y("best:Q", scale=alt.Scale(zero=False))))
    if not sub.empty:
        layers.append(alt.Chart(sub).mark_point(shape="diamond", size=140, filled=True, color="#d62728").encode(
            x="minute:Q", y="score:Q", tooltip=["config", alt.Tooltip("score:Q", format=".3f")]))
    return alt.layer(*layers).properties(height=320) if layers else None


def grouped_log(logs):
    """Warnings/errors grouped by message (digits masked), most frequent first."""
    import re
    groups = {}
    for t, lv, where, msg in logs:
        key = (lv, where, re.sub(r"\d+(\.\d+)?", "#", msg))
        g = groups.setdefault(key, {"count": 0, "last": t, "level": lv, "where": where, "example": msg})
        g["count"] += 1
        if t >= g["last"]:
            g["last"], g["example"] = t, msg
    rows = sorted(groups.values(), key=lambda g: -g["count"])
    return pd.DataFrame([{"count": g["count"], "level": g["level"], "where": g["where"],
                          "last (UTC)": datetime.fromtimestamp(g["last"], timezone.utc).strftime("%H:%M:%S"),
                          "message": g["example"]} for g in rows])


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
    running = not str(r.get("mode", "")).startswith("finished")

    # header
    c1, c2, c3 = st.columns([2, 1, 1])
    c1.subheader(f"Round {r.get('tag', '-')}  ·  {r.get('region', '')}")
    c2.metric("Status", r.get("mode", "-"))
    left = (r["deadline_ts"] - time.time()) if r.get("deadline_ts") else None
    c3.metric("Time left", mmss(left) if left is not None and running else "-")
    st.caption(f"{r.get('kind') or ''} run {r.get('run_id') or ''} · updated {ago(s.get('updated'))}")

    # evaluations: measured counts only
    scored, zero = srch.get("scored", srch.get("evaluated", 0)), srch.get("zero", 0)
    runs_done, delivered = srch.get("gatk_runs_done") or 0, srch.get("delivered") or 0
    st.markdown("#### Evaluations")
    # evaluated = configs that actually ran through GATK (core time + budget); rejected ones never ran
    evaluated = delivered if delivered else scored
    gatk_failed = max(evaluated - scored, 0)
    rejected = max(srch.get("failed", 0) - gatk_failed, 0)
    minutes = max(((s.get("updated") or time.time()) - started) / 60, 1e-9)
    m = st.columns(7)
    m[0].metric("Evaluated", _n(evaluated), help="Configs that ran through GATK this round: each one used core time "
                "and eval budget (useful + wasted + failed in GATK/scoring).")
    m[1].metric("Useful", _n(scored - zero), help="Got a real score.")
    m[2].metric("Wasted (scored 0)", f"{zero} ({100 * zero / evaluated:.0f}%)" if evaluated else "0",
                help="Ran, but no usable calls / plausibility gate failed: core time that taught us nothing.")
    m[3].metric("Failed", gatk_failed, help=f"Ran but GATK or scoring crashed. Not counted: {rejected} proposals "
                "rejected before running (outside the platform's parameter ranges; no core time).")
    m[4].metric("Running", srch.get("in_flight", 0) if running else 0)
    m[5].metric("Evaluated / min", f"{evaluated / minutes:.1f}", help="Whole round so far.")
    m[6].metric("Configs per GATK run", f"{delivered / runs_done:.1f}" if runs_done and delivered else "-",
                help=f"Configs out of GATK / GATK runs finished this round ({delivered} / {runs_done}). Batching runs "
                     "several configs in one GATK run; above 1 is the gain.")

    # one chart: every score over time
    st.markdown("#### Scores over time")
    chart = scores_chart(srch.get("points") or [], srch.get("best_timeline") or [], s.get("submissions") or [],
                         started)
    if chart is not None:
        st.altair_chart(chart, width="stretch")
        st.caption("dots: every scored config (0-scores left out) · line: best so far · red diamonds: submissions")
    else:
        st.write("no scores yet")

    # best config
    st.markdown("#### Best config")
    d = srch.get("best_detail") or {}
    b = st.columns(5)
    b[0].metric("Best local score", f"{srch['best']:.3f}" if srch.get("best") is not None else "-")
    b[1].metric("Core (weighted F1)", f"{d['core']:.4f}" if d.get("core") is not None else "-")
    b[2].metric("FP per target", d.get("fp_per_target", "-"), help="region FP / N; score term exp(-fp/8)")
    b[3].metric("Plausibility gate", {True: "passed", False: "FAILED"}.get(d.get("gate_passed"), "-"))
    subs = s.get("submissions") or []
    b[4].metric("Submissions", len(subs))
    if d:
        with st.expander(f"Details of {d.get('config')}: per class, SNP ratios, parameters, submissions"):
            c1, c2 = st.columns([3, 2])
            c1.markdown("**Per class**, scored against our answer set")
            c1.dataframe(pd.DataFrame([{"class": k, **v} for k, v in (d.get("classes") or {}).items()]),
                         hide_index=True, width="stretch")
            ra = d.get("ratios") or {}
            c2.markdown("**SNP ratios** (plausibility gate)")
            c2.dataframe(pd.DataFrame([
                {"ratio": "het / hom", "our answer set": ra.get("hethom_truth_snp"), "this config": ra.get("hethom_query_snp")},
                {"ratio": "Ti / Tv", "our answer set": ra.get("titv_truth_snp"), "this config": ra.get("titv_query_snp")}]),
                hide_index=True, width="stretch")
            c2.caption(f"region FP {d.get('region_fp')} · N {d.get('n_targets')}")
            if subs:
                st.markdown("**Submissions** (the last one before the deadline counts)")
                st.dataframe(pd.DataFrame([{"minute": round((t - started) / 60, 1), "config": k,
                                            "local score": round(v, 4)} for t, k, v in subs]),
                             hide_index=True, width="stretch")
            if d.get("params"):
                st.markdown("**Parameters**")
                st.json(d["params"])

    # answer set, one line
    st.markdown("#### Answer set")
    if a:
        g, y = _g(a, "composition", "germline") or {}, _g(a, "composition", "synthetic") or {}
        secs = a.get("seconds") or {}
        parts = [f"**{a.get('regime')} / {a.get('donor')}**",
                 f"germline {_n(g.get('total', a.get('germline')))} (het:hom {g.get('het_hom', '-')})",
                 f"synthetic {_n(y.get('total', a.get('synthetic')))} ({_n(y.get('confident'))} confident)",
                 f"N {a.get('n_targets')} ({a.get('n_targets_source') or '-'})"]
        if secs:
            parts.append(f"built in {sum(secs.values()):.0f} s")
        st.markdown(" · ".join(parts))
        if a.get("degraded"):
            st.warning(a["degraded"])
        comp = a.get("composition") or {}
        if comp:
            with st.expander("Answer-set composition"):
                rows = [("records", "total"), ("SNP", "snp"), ("indel", "indel"), ("het (GT 0/1, 1/2)", "het"),
                        ("hom (GT 1/1)", "hom"), ("het : hom", "het_hom"), ("Ti / Tv (SNP)", "titv"),
                        ("confident (p ≥ 0.99)", "confident"), ("mean p", "mean_p")]
                st.table(pd.DataFrame({"": [x[0] for x in rows], "germline": [_s(g.get(k)) for _, k in rows],
                                       "synthetic": [_s(y.get(k)) for _, k in rows]}).set_index(""))
    else:
        st.write("building...")

    # evaluator
    st.markdown(f"#### Evaluator ({ex.get('kind', '-')})")
    if ex:
        rd = ex.get("round") or {}
        now = (f"workers busy {ex.get('running', 0)} / {ex.get('workers', 0)} · queued {ex.get('queued', 0)}"
               if running else f"{ex.get('workers', 0)} workers · round over")
        done = f", {rd['gatk_runs_done']} finished" if "gatk_runs_done" in rd else ""
        st.write(f"{now} · this round: {rd.get('gatk_runs', 0)} GATK runs started{done} "
                 f"({rd.get('batches_run', 0)} batched), {rd.get('repairs_run', 0)} divergence repairs")
        if ex.get("nodes"):
            st.dataframe(pd.DataFrame(ex["nodes"]), hide_index=True, width="stretch")

    logs = s.get("log") or []
    if logs:
        with st.expander(f"Warnings and errors ({len(logs)} lines, grouped)"):
            st.dataframe(grouped_log(logs), hide_index=True, width="stretch")


# ---- history -----------------------------------------------------------------------------------------------------------
def load_runs() -> list:
    """One dict per run: its round line + the latest review / rank lines (main.review). Old lines still load."""
    path = STATUS / "rounds.jsonl"
    runs, by_id = [], {}
    if not path.exists():
        return runs
    for line in path.read_text().splitlines():
        try:
            row = json.loads(line)
        except json.JSONDecodeError:
            continue
        t = row.get("type")
        if t == "round":
            run = {"round": row, "review": None, "rank": None}
            runs.append(run)
            by_id[row.get("run_id") or f"{row['round_id']}@{row.get('ended')}"] = run
        elif t in ("review", "rank") and row.get("run_id") in by_id:
            by_id[row["run_id"]][t] = row
        elif t == "true":                                   # old format: attach to the latest run of that round
            for run in reversed(runs):
                if run["round"]["round_id"] == row["round_id"]:
                    run["rank"] = {"true_score": row.get("true_score"), "rank": row.get("rank"),
                                   "field": row.get("field")}
                    break
    return sorted(runs, key=lambda r: r["round"].get("ended") or 0, reverse=True)


def _g(d, *keys):
    for k in keys:
        d = d.get(k) if isinstance(d, dict) else None
    return d


def _when(r) -> str:
    return datetime.fromtimestamp(r["round"]["ended"], timezone.utc).strftime("%m-%d %H:%M") \
        if r["round"].get("ended") else ""


def _state(r) -> str:
    if r["rank"]:
        return "done"
    if r["review"]:
        return "waiting for rank"
    return "waiting for truth" if r["round"].get("run_id") else "-"


def _trend(rows: list, cols: list, title: str, zero_line: bool = False, kind: str = "line"):
    """One chart over runs (oldest first): y-axis zoomed to the data; `zero_line` draws y = 0 for deltas."""
    import altair as alt
    pts = [{"run": i, "round": r.get("round"), "measure": c, "value": r[c]}
           for i, r in enumerate(rows[::-1], 1) for c in cols if r.get(c) is not None]
    if not pts:
        return
    st.markdown(title)
    df = pd.DataFrame(pts)
    lo, hi = float(df["value"].min()), float(df["value"].max())
    if zero_line:
        lo, hi = min(lo, 0.0), max(hi, 0.0)
    pad = (hi - lo) * 0.15 or max(abs(hi) * 0.01, 0.01)        # flat series still get a readable axis
    enc = dict(x=alt.X("run:O", title="run (oldest first)"),
               y=alt.Y("value:Q", scale=alt.Scale(domain=[lo - pad, hi + pad], nice=False),
                       axis=alt.Axis(format="~g"), title=None),
               tooltip=["round", "measure", alt.Tooltip("value:Q", format=".4f")])
    if len(cols) > 1:
        enc["color"] = alt.Color("measure:N", title=None)
    if kind == "bar":
        chart = alt.Chart(df).mark_bar().encode(**{**enc, "y": alt.Y("value:Q", title=None)})
    else:
        chart = alt.Chart(df).mark_line(point=True).encode(**enc)
    if zero_line:
        chart = chart + alt.Chart(pd.DataFrame({"y": [0]})).mark_rule(strokeDash=[4, 4], color="gray").encode(y="y:Q")
    st.altair_chart(chart.properties(height=220), width="stretch")


def history():
    runs = load_runs()
    if not runs:
        st.info("No finished rounds yet.")
        return
    st.caption("One row per run (a round mined once live, or replayed). main.review adds the TRUE side once the "
               "round's truth is released (about one round later) and the rank once it is finalized.")
    t_over, t_as, t_ev, t_se = st.tabs(["Overview", "Answer set", "Evaluator", "Search"])

    with t_over:
        rows = []
        for r in runs:
            rd, rk, rv = r["round"], r["rank"] or {}, r["review"] or {}
            rows.append({"finished (UTC)": _when(r), "round": rd.get("tag"), "kind": rd.get("kind"),
                         "donor": f"{_g(rd, 'answer_set', 'regime') or rd.get('regime')}/"
                                  f"{_g(rd, 'answer_set', 'donor') or rd.get('donor')}",
                         "mode": rd.get("mode"),
                         "configs": _g(rd, "evaluator", "evaluated") or rd.get("evaluated"),
                         "best local": _g(rd, "search", "best_local") or rd.get("best_local"),
                         "final true": _g(rv, "search", "final_true") or rk.get("true_score"),
                         "rank": f"{rk['rank']}/{rk['field']}" if rk.get("rank") else None,
                         "vs field best": -rk["gap_to_best"] if rk.get("gap_to_best") is not None else None,
                         "vs field median": rk.get("vs_median"),
                         "top %": round(100 * rk["rank"] / rk["field"], 1) if rk.get("rank") and rk.get("field")
                         else None,
                         "local - true": _g(rv, "search", "local_minus_true"),
                         "review": _state(r), "note": rd.get("note") or ""})
        st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
        st.caption("Raw scores change with every round's difficulty, so the charts show deltas: our final config's "
                   "TRUE score minus the field's best / median (0 = level with it), and our local minus true score.")
        _trend(rows, ["vs field best", "vs field median"], "**Our true score vs the field**", zero_line=True)
        _trend(rows, ["top %"], "**Rank as % of the field** (lower is better)")
        _trend(rows, ["local - true"], "**Local - true score of the final config** (0 = our answer set judged it "
                                       "exactly)", zero_line=True)

    with t_as:
        rows = []
        for r in runs:
            a, acc = r["round"].get("answer_set") or {}, (r["review"] or {}).get("answer_set") or {}
            g, y = _g(a, "composition", "germline") or {}, _g(a, "composition", "synthetic") or {}
            rows.append({"finished (UTC)": _when(r), "round": r["round"].get("tag"),
                         "donor": f"{a.get('regime') or r['round'].get('regime')}/{a.get('donor') or r['round'].get('donor')}",
                         "germline": g.get("total"),
                         "germ het:hom": g.get("het_hom"), "germ Ti/Tv": g.get("titv"),
                         "synthetic": y.get("total"), "syn het:hom": y.get("het_hom"),
                         "syn confident": y.get("confident"), "N ours": a.get("n_targets"),
                         "N truth": _g(acc, "n_targets", "truth"), "build s": a.get("build_s"),
                         "germline exact P": _g(acc, "germline_exact", "precision"),
                         "germline exact R": _g(acc, "germline_exact", "recall"),
                         "synthetic P": _g(acc, "synthetic", "precision"),
                         "synthetic R": _g(acc, "synthetic", "recall"),
                         "synthetic GT P": _g(acc, "synthetic_gt", "precision"),
                         "degraded": a.get("degraded") or ""})
        st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
        st.caption("P = share of ours that is in the released truth; R = share of the truth we have. Germline exact = "
                   "same position, alleles and genotype (what the score treats as exact).")
        _trend(rows, ["germline exact P", "germline exact R"], "**Germline accuracy** (zoomed)")
        _trend(rows, ["synthetic P", "synthetic R", "synthetic GT P"], "**Synthetic accuracy** (zoomed)")
        for x in rows:
            x["N ours / truth"] = round(x["N ours"] / x["N truth"], 4) if x.get("N ours") and x.get("N truth") else None
        _trend(rows, ["N ours / truth"], "**N estimate / real N** (1 = exact)")

    with t_ev:
        rows = []
        for r in runs:
            e = r["round"].get("evaluator") or {}
            rows.append({"finished (UTC)": _when(r), "round": r["round"].get("tag"), "executor": e.get("kind"),
                         "nodes": e.get("nodes"), "configs": e.get("evaluated"), "per min": e.get("per_min"),
                         "GATK runs": e.get("gatk_runs"), "configs per GATK run": e.get("per_gatk_run"),
                         "batches": e.get("batches"), "repairs": e.get("repairs"), "failed": e.get("failed"),
                         "worker utilization": e.get("utilization"),
                         "review re-runs": (r["review"] or {}).get("rerun")})
        st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
        _trend(rows, ["configs"], "**Configs evaluated per run** (total)")
        _trend(rows, ["per min"], "**Configs per minute**")
        _trend(rows, ["configs per GATK run"], "**Configs per GATK run** (batching gain)")

    with t_se:
        rows = []
        for r in runs:
            sr, rv = r["round"].get("search") or {}, _g(r["review"] or {}, "search") or {}
            rows.append({"finished (UTC)": _when(r), "round": r["round"].get("tag"),
                         "best local": sr.get("best_local"), "min to best": sr.get("time_to_best_min"),
                         "baseline local": sr.get("baseline_local"),
                         "scored 0": f"{sr['zero_scored']}/{sr['scored']}" if sr.get("scored") else None,
                         "submissions": sr.get("submissions"), "final true": rv.get("final_true"),
                         "baseline true": rv.get("baseline_true"), "gain over baseline": rv.get("gain_over_baseline"),
                         "best true (top-k)": rv.get("best_true"), "pick gap": rv.get("pick_gap"),
                         "spearman local~true": rv.get("spearman")})
        st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
        st.caption("Top-k = the final, the baseline and the best configs by local score, each re-scored against the "
                   "released truth. Pick gap = best true among them - true of the final (0 = we picked the best). "
                   "Spearman: does our local score order configs like the truth does?")
        _trend(rows, ["gain over baseline"], "**Gain over the baseline config** (TRUE score; above 0 = search helped)",
               zero_line=True)
        _trend(rows, ["pick gap"], "**Pick gap** (best TRUE among top-k - TRUE of our final; 0 = we picked the best)")
        _trend(rows, ["spearman local~true"], "**Spearman local ~ true** (1 = our local score orders configs like "
                                              "the truth)")
        reviewed = [r for r in runs if r["review"]]
        if reviewed:
            pick = st.selectbox("Configs of one reviewed run", [f"{_when(r)} {r['round'].get('tag')}" for r in reviewed])
            r = reviewed[[f"{_when(x)} {x['round'].get('tag')}" for x in reviewed].index(pick)]
            st.dataframe(pd.DataFrame(r["review"].get("configs") or []), hide_index=True, width="stretch")


st.title("🧬 Minos miner")
tab_live, tab_hist = st.tabs(["Live round", "History"])
with tab_live:
    live()
with tab_hist:
    history()
