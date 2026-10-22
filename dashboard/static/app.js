"use strict";
/* Minos live dashboard: polls /api/status every 2 s and draws it; the Settings tab edits tuning.json through
   /api/settings. No inline styles or scripts (the server's CSP forbids them): elements get styles through the DOM. */

const POLL_MS = 2000;
const STALE_S = 30;
const $ = (id) => document.getElementById(id);
const SVG = "http://www.w3.org/2000/svg";

let serverOffset = 0;          // server clock - browser clock, seconds
let last = null;               // the last /api/status payload
let lastResult = null;         // the run's round/review/rank rows (rounds.jsonl), once main has written them
let liveTrue = null;          // true score of the current best, written beside the dashboard while the round runs
let lastRank = null;           // where the live best sits on the published field, refreshed with each poll
let chart = null, chartKey = "";
let coresChart = null, coresChartKey = "";
let liveCoreTrace = [];        // [unix t, busy, capacity] while this page is open
let liveCoreRun = null;
let settingsMeta = null;       // /api/settings
const draft = {};              // setting name -> text typed but not saved
const resets = new Set();      // settings to put back to their defaults

// ---- helpers --------------------------------------------------------------------------------------------------------
function h(tag, props, ...kids) {
  const el = document.createElement(tag);
  for (const [k, v] of Object.entries(props || {})) {
    if (v == null || v === false) continue;
    if (k === "class") el.className = v;
    else if (k === "style") for (const [sk, sv] of Object.entries(v)) { if (sk.startsWith("--")) el.style.setProperty(sk, sv); else el.style[sk] = sv; }
    else if (k === "text") el.textContent = v;
    else if (k.startsWith("on")) el.addEventListener(k.slice(2), v);
    else el.setAttribute(k, v === true ? "" : v);
  }
  for (const c of kids.flat(Infinity)) if (c != null && c !== false) el.append(c instanceof Node ? c : String(c));
  return el;
}
function s(tag, attrs, ...kids) {
  const el = document.createElementNS(SVG, tag);
  for (const [k, v] of Object.entries(attrs || {})) if (v != null) el.setAttribute(k, v);
  for (const c of kids.flat()) if (c != null) el.append(c instanceof Node ? c : document.createTextNode(String(c)));
  return el;
}
const now = () => Date.now() / 1000 + serverOffset;
const fix = (x, d = 4) => (x == null || Number.isNaN(x) ? "–" : Number(x).toFixed(d));
const int = (x) => (x == null ? "–" : Math.round(x).toLocaleString("en-US"));
const pct = (a, b) => (b ? Math.max(0, Math.min(100, (100 * a) / b)) : 0);
function dur(sec) {
  if (sec == null || !isFinite(sec)) return "–";
  const neg = sec < 0; sec = Math.abs(Math.round(sec));
  const m = Math.floor(sec / 60), r = sec % 60;
  const out = m >= 60 ? `${Math.floor(m / 60)}h ${String(m % 60).padStart(2, "0")}m` : `${m}:${String(r).padStart(2, "0")}`;
  return neg ? `-${out}` : out;
}
const mins = (sec) => (sec == null ? "–" : `${(sec / 60).toFixed(1)} min`);
function clock(ts) {
  const d = new Date(ts * 1000);
  return d.toLocaleTimeString("en-GB", {hour12: false});
}
function ago(ts) {
  const a = now() - ts;
  if (a < 5) return "just now";
  if (a < 90) return `${Math.round(a)} s ago`;
  if (a < 5400) return `${Math.round(a / 60)} min ago`;
  return `${(a / 3600).toFixed(1)} h ago`;
}
function bar(frac, cls, mark) {
  const b = h("div", {class: `bar ${cls || ""}`}, h("span", {style: {width: `${Math.max(0, Math.min(100, frac * 100))}%`}}));
  if (mark != null) b.append(h("span", {class: "mark", style: {width: `${Math.max(0, Math.min(100, mark * 100))}%`}}));
  return b;
}
const chip = (text, cls) => h("span", {class: `chip ${cls || ""}`, text});
function kv(pairs) {
  return h("dl", {class: "kv"}, pairs.filter(Boolean).map(([k, v]) => [h("dt", {text: k}), h("dd", {}, v)]));
}
function css(name) { return getComputedStyle(document.documentElement).getPropertyValue(name).trim(); }
function toast(text) {
  const t = h("div", {class: "toast", role: "status", text});
  document.body.append(t);
  setTimeout(() => t.remove(), 4000);
}
async function api(path, opts = {}) {
  const r = await fetch(path, {...opts, headers: {"Content-Type": "application/json", "X-Minos": "1", ...(opts.headers || {})}});
  if (r.status === 401) { location.replace("/login"); throw new Error("logged out"); }
  return r;
}

// ---- tabs -----------------------------------------------------------------------------------------------------------
function showTab() {
  const hash = ["#settings", "#history", "#curves"].includes(location.hash) ? location.hash.slice(1) : "live";
  $("view-live").hidden = hash !== "live";
  $("view-history").hidden = hash !== "history";
  $("view-curves").hidden = hash !== "curves";
  $("view-settings").hidden = hash !== "settings";
  $("tab-live").classList.toggle("on", hash === "live");
  $("tab-history").classList.toggle("on", hash === "history");
  $("tab-curves").classList.toggle("on", hash === "curves");
  $("tab-settings").classList.toggle("on", hash === "settings");
  if (hash === "settings") loadSettings();
  else if (hash === "history") loadHistory();
  else if (hash === "curves") loadCurves();
  else if (chart) setTimeout(sizeChart, 0);
}
window.addEventListener("hashchange", showTab);

// the Settings bar sticks just under the top bar, whose height depends on how it wraps
const topBar = document.querySelector(".top");
function setTopHeight() { document.documentElement.style.setProperty("--top-h", `${topBar.offsetHeight}px`); }
setTopHeight();
if (window.ResizeObserver) new ResizeObserver(setTopHeight).observe(topBar);

// ---- theme: dark only (user, 2026-10-21) ----
function effectiveTheme() { return "dark"; }
$("drawer-close").addEventListener("click", closeDetail);
$("drawer").addEventListener("click", (e) => { if (e.target === $("drawer")) closeDetail(); });
document.addEventListener("keydown", (e) => { if (e.key === "Escape" && detail) closeDetail(); });
$("logout").addEventListener("click", async () => {
  await fetch("/api/logout", {method: "POST", headers: {"X-Minos": "1"}}).catch(() => {});
  location.replace("/login");
});

// ---- polling --------------------------------------------------------------------------------------------------------
async function poll() {
  try {
    const r = await api("/api/status");
    const data = await r.json();
    serverOffset = data.now - Date.now() / 1000;
    last = data.live;
    lastResult = data.result || null;
    lastRank = data.rank || null;
    noteLiveCores(last);
    try {
      const tr = await fetch("/static/true_live.json?t=" + Date.now());
      liveTrue = tr.ok ? await tr.json() : null;
    } catch (e) {
      liveTrue = null;
    }
    render();
  } catch (e) {
    if (e.message !== "logged out") $("updated").textContent = "server not answering";
  } finally {
    setTimeout(poll, POLL_MS);
  }
}

// ---- live view ------------------------------------------------------------------------------------------------------
function render() {
  const L = last;
  const banners = $("banners");
  banners.replaceChildren();
  if (!L) {
    $("round-name").textContent = "no round yet";
    banners.append(h("div", {class: "banner info", text: "Main has not written a status file yet. It appears once main starts a round."}));
    renderSettingsState();
    return;
  }
  const R = L.round || {}, S = L.search || {}, E = S.engine || null, X = L.executor || null;
  $("round-name").textContent = `${R.tag || R.id || "round"} · ${R.region || ""} · ${R.kind || ""}`;
  const finished = (R.mode || "").startsWith("finished");
  const phase = finished ? "finished" : E ? E.phase : (R.mode || "starting");
  const pc = $("phase-chip");
  pc.textContent = phase;
  pc.className = `chip ${finished ? "finished" : phaseClass(phase)}`;
  const age = now() - (L.updated || 0);
  $("updated").textContent = `updated ${ago(L.updated || 0)}`;
  if (!finished && age > STALE_S) {
    banners.append(h("div", {class: "banner bad", text: `Main has not updated the status for ${dur(age)}. It may have stopped, or the round ended without a final write.`}));
  }
  if (finished) banners.append(h("div", {class: "banner info", text: `Round ${R.mode}. The panels below show the round as it ended; the next round appears when main starts it.`}));
  if (L.schema !== 2) banners.append(h("div", {class: "banner", text: "This status file comes from an older main. Some panels stay empty until main runs the new code."}));
  const set = L.settings || {};
  if (set.error) banners.append(h("div", {class: "banner bad", text: `Settings file: ${set.error}`}));

  renderKpis(L, R, S, E, X, finished);
  renderResult(lastResult, finished);
  renderSchedule(E);
  renderLineage(E);
  renderGlobalA(S.global_a || null);
  renderIslands(E);
  renderClimbers(E, finished, X);
  renderDetail(E, finished);
  renderStepB(E, X);
  renderStepC(E);
  renderStepA(E, R);
  renderChart(S, L.submissions || []);
  renderCoresChart(S, X, R);
  renderNodes(X, finished);
  renderFeed(L, E);
  renderSettingsState();
}

function kpi(label, value, note, cls) {
  return h("div", {class: `kpi ${cls || ""}`}, h("div", {class: "label", text: label}), h("div", {class: "value", text: value, title: value}),
    note ? h("div", {class: "note", text: note}) : null);
}

function shownTrue(R, S) {
  // The reviewed rank is the pick scored on the released truth. Until that row exists, the current best's true score
  // comes from true_live.json. The local score is a different answer set, so it is not a place on the published field.
  const reviewed = lastResult && lastResult.rank;
  if (reviewed && reviewed.true_score != null) return reviewed;
  if (liveTrue && liveTrue.tag === (R && R.tag) && liveTrue.state === "ready"
      && liveTrue.config === (S && S.best_config) && liveTrue.true_score != null) return liveTrue;
  return null;
}

function rankKpi(row) {
  if (!row || row.rank == null) {
    const scoring = liveTrue && liveTrue.state === "scoring";
    return kpi("True rank", "–", scoring ? "scoring the current best" : "waiting for the true score", "k-rank");
  }
  const gap = row.gap_to_best == null ? "" : `gap ${fix(row.gap_to_best)}`;
  const best = row.field_best == null ? "" : `best ${fix(row.field_best)}`;
  const note = [gap, best].filter(Boolean).join(" · ");
  return kpi("True rank", `${int(row.rank)} / ${int(row.field)}`, note, "k-rank");
}

function phaseClass(phase) {
  // the chip's colour: the step the engine is in
  if (!phase) return "";
  if (phase.startsWith("A")) return "opening";
  if (phase.startsWith("B")) return "climb";
  if (phase.startsWith("C")) return "merge";
  return phase.split(" ")[0];
}

function climbers(E) { return (E && E.basins) || []; }

function renderKpis(L, R, S, E, X, finished) {
  const subs = L.submissions || [];
  const lastSub = subs.length ? subs[subs.length - 1] : null;
  const left = R.deadline_ts ? R.deadline_ts - now() : null;
  const perRun = S.gatk_runs_done ? (S.delivered || 0) / S.gatk_runs_done : null;
  const A = L.answer_set;
  const truth = shownTrue(R, S);
  $("kpis").replaceChildren(
    kpi("Time left", finished ? "done" : left != null ? dur(left) : "–",
      E && E.schedule ? `search ${dur(E.elapsed_s)} of ${dur(E.schedule.length_s)}` : (R.mode || ""), "k-time"),
    kpi("Best local", truth ? fix(truth.true_score) : fix(S.best),
      truth ? `true score · local ${fix(S.best)}` : (E && E.pick ? `pick ${fix(E.pick.score)}` : S.best_config ? `config ${S.best_config.slice(0, 8)}` : "local score"), "k-best"),
    rankKpi(truth),
    kpi("Submitted", lastSub ? fix(lastSub[2]) : "–", lastSub ? `${subs.length} so far · ${ago(lastSub[0])}` : "nothing yet", "k-sub"),
    kpi("Configs scored", int(S.scored), `${int(S.evaluated)} run · ${int(S.failed)} failed · ${int(S.zero)} scored 0`, "k-configs"),
    finished ? kpi("Workers busy", "–", "round over", "k-workers")
      : kpi("Workers busy", X ? int(X.running) : "–", X ? `of ${int(X.workers)} · ${int(X.queued)} queued · ${int(S.in_flight)} in flight` : "", "k-workers"),
    kpi("GATK runs", int(S.gatk_runs), perRun ? `${perRun.toFixed(2)} configs per run` : "", "k-gatk"),
    kpi("Answer set", A ? int(A.germline) : "building",
      A ? `germline + ${int(A.synthetic)} synthetic · ${A.regime || ""} · donor ${A.donor || "?"} · ${int(A.n_targets)} targets` : "", "k-answer"),
  );
}

function segText(cls, label, share) {
  // Step A sits in the first few minutes, so the full words do not always fit in the band.
  if (cls === "climb") return share >= 0.04 ? "Step B" : "B";
  if (share >= 0.1) return label;
  if (cls === "opening") return "W0";
  if (cls === "generations") return "Gen 1";
  if (cls === "stepc") return "C";
  if (cls === "finetune") return "FT";
  return share >= 0.05 ? label : "";
}

function segments(E) {
  // Wave 0 is capped by candidates.max_min. Generation 1 then gets candidates.generation_min of its own.
  // v9.29 skips generation 1 and the GA opening: Wave 0, the dummy peaks, and the valley test are Step A,
  // then Step B, then the fine-tune. A step that has finished uses its real end; until then the bar shows the cap.
  const sc = E.schedule || {};
  const len = sc.length_s || 0;
  const span = (a, b) => {
    const lo = Math.max(0, a);
    return [lo, Math.max(lo, Math.min(len, b))];
  };
  const pick = sc.pick_s != null ? sc.pick_s : len;
  const aDone = sc.a_done_s, bDone = sc.b_done_s;
  if ((E.version || "").startsWith("v9.34") || (E.version || "").startsWith("v9.33") || (E.version || "").startsWith("v9.32") || (E.version || "").startsWith("v9.31") || (E.version || "").startsWith("v9.30") || (E.version || "").startsWith("v9.29")) {
    const aEnd = aDone != null ? aDone : Math.max(E.elapsed_s || 0, sc.a_max_s || 0);
    const bEnd = bDone != null ? bDone : Math.max(aEnd, sc.b_end_s || aEnd);
    return [
      ["opening", "Wave 0", ...span(0, aEnd)],
      ["climb", "Step B", ...span(aEnd, bEnd)],
      ["finetune", "Fine-tune", ...span(bEnd, pick)],
      ["final", "pick", ...span(pick, len)],
    ].filter((g) => g[3] - g[2] > 1);
  }
  const gens = (E.waves || []).filter((w) => w.started_s != null && /^generation /.test(w.label || ""));
  const genStart = gens.length ? Math.min(...gens.map((w) => w.started_s)) : null;
  if (sc.gen_s == null) {
    const aEnd = aDone != null ? aDone : Math.max(E.elapsed_s || 0, sc.a_max_s || 0);
    const bEnd = bDone != null ? bDone : Math.max(aEnd, sc.b_end_s || aEnd);
    const parts = [];
    if (genStart != null && genStart > 1 && aEnd - genStart > 1) {
      parts.push(["opening", "Wave 0", ...span(0, genStart)]);
      parts.push(["generations", "Generation 1", ...span(genStart, aEnd)]);
    } else parts.push(["opening", "Step A", ...span(0, aEnd)]);
    parts.push(["climb", "Step B", ...span(aEnd, bEnd)]);
    parts.push(["stepc", "Step C", ...span(bEnd, pick)]);
    parts.push(["final", "pick", ...span(pick, len)]);
    return parts.filter((g) => g[3] - g[2] > 1);
  }
  const w0End = genStart != null ? genStart : (aDone != null ? aDone : (sc.a_max_s || 0));
  const genEnd = aDone != null ? aDone : (genStart != null ? genStart + sc.gen_s : w0End + sc.gen_s);
  const bStart = aDone != null ? aDone : genEnd;
  const bStop = bDone != null ? bDone : Math.max(bStart, sc.b_end_s || bStart);
  const parts = [
    ["opening", "Wave 0", ...span(0, w0End)],
    ["generations", "Generation 1", ...span(w0End, genEnd)],
  ];
  const gaOpenCap = sc.ga_open_s || 0;
  let climbStart = bStart;
  if (gaOpenCap > 1) {
    const gaEnd = sc.ga_open_done_s != null
      ? sc.ga_open_done_s
      : Math.min(bStop, bStart + gaOpenCap, Math.max(bStart, E.elapsed_s || bStart));
    parts.push(["generations", "GA open", ...span(bStart, gaEnd)]);
    climbStart = gaEnd;
  }
  parts.push(["climb", "Step B", ...span(climbStart, bStop)]);
  const ft = sc.c_reserve_s != null ? pick - sc.c_reserve_s : null;
  if (ft != null && ft < pick) {
    if (ft > bStop + 1) {
      parts.push(["stepc", "Step C GA", ...span(bStop, ft)]);
      parts.push(["finetune", "Fine-tune", ...span(ft, pick)]);
    } else parts.push(["finetune", "Fine-tune", ...span(bStop, pick)]);
  } else parts.push(["stepc", "Step C", ...span(bStop, pick)]);
  parts.push(["final", "pick", ...span(pick, len)]);
  return parts.filter((g) => g[3] - g[2] > 1);
}

function roundAxis(E) {
  // the round's time axis, shared by the schedule bar and the climber tree: from the round's start to its deadline.
  // The engine's times are seconds from its own start (`off` into the round); the search stops `margin` before
  // the deadline, which leaves main the time to submit.
  const R = (last && last.round) || {};
  const vals = (last && last.settings && last.settings.values) || {};
  const margin = vals["main.stop_before_deadline_s"] != null ? vals["main.stop_before_deadline_s"] : 60;
  const sc = (E && E.schedule) || null;
  const off = E && R.started_ts && E.started_ts ? Math.max(0, E.started_ts - R.started_ts) : 0;
  const total = R.started_ts && R.deadline_ts ? R.deadline_ts - R.started_ts : (sc ? off + sc.length_s + margin : 72 * 60);
  const searchEnd = sc ? Math.min(total, off + sc.length_s) : total - margin;
  return {total, off, searchEnd, start: R.started_ts};
}

function renderAnswerLane(ax, finished) {
  const lane = $("answer-lane");
  const A = last && last.answer_set;
  if (!ax.start) { lane.replaceChildren(); return; }
  if (A && !A.ready_ts) { lane.replaceChildren(); return; }               // a status file from before main recorded it
  const ready = A && A.ready_ts ? A.ready_ts - ax.start : null;
  const end = ready != null ? ready : Math.max(5, now() - ax.start);
  const frac = ax.total > 0 ? Math.max(0, Math.min(1, end / ax.total)) : 0;
  const pct = `${frac * 100}%`;
  const label = ready != null ? `answer set ${Math.round(ready)} s` : "answer set…";
  // The bar is the build's true length. A build of about a minute is a sliver of a 72-minute round,
  // so the words sit beside it instead of inside it.
  const wide = frac > 0.7;
  lane.replaceChildren(
    h("div", {class: "seg", style: {left: "0", width: pct}, title: label}, wide ? label : null),
    wide ? null : h("span", {class: "lane-label", style: {left: pct}, text: label}),
  );
}

function renderSchedule(E) {
  const el = $("timeline-sub");
  if (!el) return;
  if (!E || !E.schedule) { el.textContent = E ? E.version || "" : ""; return; }
  const sc = E.schedule;
  const plan = (E.version || "").startsWith("v9.34") || (E.version || "").startsWith("v9.33") || (E.version || "").startsWith("v9.32") || (E.version || "").startsWith("v9.31") || (E.version || "").startsWith("v9.30")
    ? `wave 0 cap ${mins(sc.a_max_s)}, a dummy on each apparent peak, two far scout waves, Step B until ${mins(sc.b_end_s)}. One compass on the highest score, then the joint cube. Free cores look 0.10 away. Fine-tune ${mins(sc.c_reserve_s)} stays on this mountain`
    : (E.version || "").startsWith("v9.29")
    ? `wave 0 cap ${mins(sc.a_max_s)}, a dummy on each apparent peak, no generation 1, a valley test, Step B until ${mins(sc.b_end_s)}. The GA breeds the view and the grid scores one setting at a time. Fine-tune ${mins(sc.c_reserve_s)}. A late Step B score replaces the maximum`
    : (E.version || "").startsWith("v9.28") || (E.version || "").startsWith("v9.27") || (E.version || "").startsWith("v9.26")
    ? `wave 0 cap ${mins(sc.a_max_s)}, generation 1 ${mins(sc.gen_s)} after it, Step B GA open ${mins(sc.ga_open_s || 0)}, Step B until ${mins(sc.b_end_s)}. Elite probes of the best, once. Then the grid keeps the fleet. Fine-tune ${mins(sc.c_reserve_s)} with no GA`
    : (E.version || "").startsWith("v9.25")
    ? `wave 0 cap ${mins(sc.a_max_s)}, generation 1 ${mins(sc.gen_s)} after it, Step B GA open ${mins(sc.ga_open_s || 0)}, Step B until ${mins(sc.b_end_s)}. While GA is moving a candidate the climbers send nothing. Then one rung. Fine-tune ${mins(sc.c_reserve_s)} with no GA`
    : (E.version || "").startsWith("v9.24") || (E.version || "").startsWith("v9.23")
    ? `wave 0 cap ${mins(sc.a_max_s)}, generation 1 ${mins(sc.gen_s)} after it, Step B GA open ${mins(sc.ga_open_s || 0)}, Step B until ${mins(sc.b_end_s)}, then a fine-tune of ${mins(sc.c_reserve_s)} with no GA`
    : sc.gen_s != null
    ? `wave 0 cap ${mins(sc.a_max_s)}, generation 1 ${mins(sc.gen_s)} after it, Step B GA open ${mins(sc.ga_open_s || 0)}, Step B until ${mins(sc.b_end_s)}, Step C keeps ${mins(sc.c_reserve_s)} to fine-tune`
    : `Step B until ${mins(sc.b_end_s)}, then Step C, then the pick`;
  el.textContent = `${E.version || "search"} · ${plan} · scaled ×${sc.scale.toFixed(2)}`;
}

function climberClass(c) {
  if (c.status === "climbing") return "live";
  if (c.status === "peak") return "final";
  if (c.status === "landed") return "landed";
  if (c.status === "joined") return "retired";
  return "stopped";
}

const STATUS_TEXT = {climbing: "climbing", peak: "at its peak", stopped: "stopped", joined: "joined another",
  landed: "short climb done"};

function renderLineage(E) {
  const box = $("lineage");
  const cl = climbers(E);
  if (!E || !E.schedule) {
    box.replaceChildren(h("div", {class: "empty", text: E ? "This round's search publishes no schedule." : "The engine has not started."}));
    return;
  }
  const rows = [...cl].sort((a, b) => a.born_s - b.born_s || a.id - b.id);
  const row = new Map(rows.map((c, i) => [c.id, i]));
  const W = Math.max(520, box.clientWidth || 640), top = 22, rh = 22;
  const left = parseFloat(css("--axis-l")) || 46, right = parseFloat(css("--axis-r")) || 66;   // same axis as the schedule
  const H = top + rows.length * rh + 8;
  const ax = roundAxis(E), L = ax.total, o = ax.off;
  const x = (t) => left + ((W - left - right) * Math.max(0, Math.min(L, t + o))) / L;
  const y = (i) => top + i * rh + rh / 2;
  const svg = s("svg", {width: W, height: H, viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Climbers over the round"});
  for (const [cls, label, a, b] of segments(E)) {
    const w = Math.max(0, x(b) - x(a));
    svg.append(s("rect", {class: `band ${cls}`, x: x(a), y: 0, width: w, height: H}));
    svg.append(s("text", {x: x(a) + 4, y: 13}, segText(cls, label, w / W)));
  }
  const tNow = E.elapsed_s;
  for (const c of rows) {                                   // a Step C landing hangs off the peak it came from
    if (c.parent == null || !row.has(c.parent)) continue;
    const xp = x(c.born_s);
    svg.append(s("path", {class: "link", d: `M ${xp} ${y(row.get(c.parent))} L ${xp} ${y(row.get(c.id))}`}));
  }
  for (const c of rows) {
    const i = row.get(c.id);
    const end = c.ended_s != null ? c.ended_s : tNow;
    const r = s("rect", {class: `row ${climberClass(c)}`, x: x(c.born_s), y: y(i) - 6, width: Math.max(3, x(end) - x(c.born_s)), height: 12, rx: 3});
    r.append(s("title", {}, `Climber ${c.id}: ${fix(c.opened)} → ${fix(c.best)}, ${STATUS_TEXT[c.status] || c.status}; ` +
      `from ${c.origin} at ${mins(c.born_s)}; ${c.iterations} grids`));
    svg.append(r);
    for (const [t, sc, why, part] of c.history || []) {     // each update in the colour of the part that found it
      const kind = why.startsWith("C") ? "c" : part && part !== "other" ? part
        : why === "ga" ? "ga" : why === "jump" ? "jump" : why === "returned" ? "returned"
        : why === "view" || why === "met" || why.endsWith("grid") || why.endsWith("line") ? "grid" : "";
      const m = s("circle", {class: kind ? `move ${kind}` : "move", cx: x(t), cy: y(i), r: 3.5});
      m.append(s("title", {}, `${mins(t)}: ${fix(sc)} · ${kind === "c" ? "Step C" : kind === "returned" ? "own run back after its wave (grid or jump, not recorded)" : PART_NAME[kind] || why} (${why})`));
      svg.append(m);
    }
    svg.append(s("text", {x: 4, y: y(i) + 4}, `C${c.id}`));
    svg.append(s("text", {class: "score", x: x(end) + 5, y: y(i) + 4}, fix(c.best, 3)));
  }
  if (tNow < L) svg.append(s("line", {class: "nowline", x1: x(tNow), x2: x(tNow), y1: 16, y2: H}));
  box.replaceChildren(svg);
}

function renderResult(res, finished) {
  const panel = $("result-panel");
  const rev = res && res.review, rank = res && res.rank;
  if (!finished && !rev && !rank) { panel.hidden = true; return; }
  panel.hidden = false;
  const box = $("result");
  if (!rev && !rank) {
    $("result-sub").textContent = "";
    box.replaceChildren(h("div", {class: "empty", text: "The truth review runs once the round's truth is released (a replay: right after the search; a live round: about one round later)."}));
    return;
  }
  const s = (rev && rev.search) || {}, a = (rev && rev.answer_set) || {};
  const syn = a.synthetic || {}, ger = a.germline_exact || {};
  $("result-sub").textContent = rank ? `reviewed ${ago(rank.at)}` : "review done; the rank follows";
  box.replaceChildren(
    kpi("True score", rank ? fix(rank.true_score) : s.final_true != null ? fix(s.final_true) : "–",
      s.final_local != null ? `local ${fix(s.final_local)} · pick gap ${fix(s.pick_gap)}` : "", "k-best"),
    kpi("Rank", rank ? `${rank.rank} / ${rank.field}` : "pending",
      rank ? `field best ${fix(rank.field_best)} · gap ${fix(rank.gap_to_best)}` : "", "k-sub"),
    kpi("Vs the field median", rank ? `${rank.vs_median >= 0 ? "+" : ""}${fix(rank.vs_median)}` : "–",
      rank ? `median ${fix(rank.field_median)}` : "", "k-configs"),
    kpi("Over the baseline", s.gain_over_baseline != null ? `+${fix(s.gain_over_baseline)}` : "–",
      s.spearman != null ? `local vs true order: Spearman ${fix(s.spearman, 2)}` : "", "k-gatk"),
    a.as_scored && a.as_scored.core != null
      ? kpi("Answer set", fix(a.as_scored.core, 4),
          `as scored · ${a.as_scored.missed} missed, ${a.as_scored.extra} extra (≤ ${fix(a.as_scored.core_points, 3)} pts) · strict synthetic ${syn.precision != null ? fix(syn.precision, 3) : "–"} / ${syn.recall != null ? fix(syn.recall, 3) : "–"}`, "k-answer")
      : kpi("Answer set", syn.precision != null ? `${fix(syn.precision, 3)} / ${fix(syn.recall, 3)}` : "–",
          `synthetic precision / recall · germline ${ger.precision != null ? fix(ger.precision, 3) : "–"} / ${ger.recall != null ? fix(ger.recall, 3) : "–"}`, "k-answer"),
  );
}

function bestUpdated(c, E) {
  if (c.best_s == null) return "updated –";
  const ts = E.started_ts != null ? E.started_ts + c.best_s : null;
  return ts != null ? `updated ${mins(c.best_s)} · ${ago(ts)}` : `updated ${mins(c.best_s)}`;
}

function renderClimbers(E, finished, X) {
  const box = $("climbers");
  const cl = climbers(E);
  if (!E) { box.replaceChildren(h("div", {class: "empty", text: "Waiting for the engine."})); $("climbers-sub").textContent = ""; return; }
  const main = cl.filter((c) => c.parent == null).sort((a, b) => a.id - b.id);
  const short = cl.filter((c) => c.parent != null);
  $("climbers-sub").textContent = short.length
    ? `${short.length} short climb${short.length === 1 ? "" : "s"} in the tree above`
    : (E.version || "").startsWith("v9.34") || (E.version || "").startsWith("v9.33") || (E.version || "").startsWith("v9.32") || (E.version || "").startsWith("v9.31") || (E.version || "").startsWith("v9.30")
      ? "one climber compasses the highest score. free cores look for a fingerprint 0.10 away, then wait"
    : (E.version || "").startsWith("v9.29")
      ? "one climber, on the highest fingerprint. the centre moves when its grid finishes"
    : (E.version || "").startsWith("v9.28")
      ? "the grid writes the setting it found. a new high in the pile does not move the centre"
    : (E.version || "").startsWith("v9.27")
      ? "elite probes go out together. a free worker walks further up the line"
    : (E.version || "").startsWith("v9.26")
      ? "elite probes of the best, then the grid keeps the fleet. a peaked centre is not taken"
    : (E.version || "").startsWith("v9.25")
      ? "while GA is still moving a candidate, the climbers send nothing and the fleet breeds. then one rung"
    : (E.version || "").startsWith("v9.24")
      ? "a config 0.10 from every pile, and above the lowest climber, takes that seat now. the fine-tune restarts from a higher score"
    : (E.version || "").startsWith("v9.23")
      ? "climbers keep their seats until the fine-tune. a peaked seat goes to the best candidate. the fine-tune climbs the best config and does not breed"
    : (E.version || "").startsWith("v9.22")
      ? "six piles. four climb. two candidates wait, and the best one takes a seat when a climber reaches the top"
    : (E.version || "").startsWith("v9.21")
      ? "one climber per callset pile. the highest score within 0.10 becomes that pile's reference"
    : (E.version || "").startsWith("v9.20")
      ? "one climber per callset pile. a peaked fingerprint is not climbed again. a higher returned run takes the centre"
    : (E.version || "").startsWith("v9.19")
      ? "one climber per callset pile, up to four. piles stay 0.10 apart. a run that comes back higher takes the centre"
    : (E.version || "").startsWith("v9.18")
      ? "four climbers. compass then a cube; leftover cores keep running GA; the far seat stays a far hill"
    : (E.version || "").startsWith("v9.17")
      ? "four climbers. compass then a cube on the live axes; one far seat; leftover GA is half elite-anchored"
    : (E.version || "").startsWith("v9.16")
      ? "four climbers. one seat is the best far hill even when it is behind 0.10"
    : (E.version || "").startsWith("v9.15")
      ? "four climbers. two GA children cannot share a hill; Step B opens with GA, then the lowest takes a free GA score"
      : (E.version || "").startsWith("v9.14")
      ? "four climbers. Step B opens with GA on the fleet; then the lowest takes the highest free GA score"
      : (E.version || "").startsWith("v9.13")
      ? "four climbers. the lowest takes the highest free GA score as soon as that score is higher"
      : (E.version || "").startsWith("v9.12")
      ? "four climbers. the lowest takes the highest GA score as soon as that score is higher"
      : (E.version || "").startsWith("v9.11")
      ? "four climbers. when a read or assembly grid ends, the lowest climber takes the highest GA score if it is behind"
      : (E.version || "").startsWith("v9.9")
      ? "five climbers. a GA best replaces the climber with the lowest best score"
      : (E.version || "").startsWith("v9.8")
      ? "a config that differs in fewer than four settings is the same hill"
      : (E.version || "").startsWith("v9.7")
      ? "an empty seat takes the best scored config that is far from every climber"
      : (E.version || "").startsWith("v9.6")
      ? "after the plateau probes, a child is sent only when it is far from every climber"
      : (E.version || "").startsWith("v9.5")
      ? "a value that has already tied a config's score is not sent again"
      : (E.version || "").startsWith("v9.4")
      ? "on a shared best score the GA probes one setting by one rung; a config within two settings of a climber is not climbed again"
      : (E.version || "").startsWith("v9.3")
      ? "seats stay within 0.10 of the best; a flat seat gives its workers away; the GA's best goes to the lowest climber"
      : (E.version || "").startsWith("v9.2")
      ? "main points first, genotyping rides along; PairHMM and the cheap list run on the centre only after those points are out; the GA is last"
      : "a free worker takes the next unsent value; the GA runs only on workers those updates leave free";
  if (!main.length) { box.replaceChildren(h("div", {class: "empty", text: (E.phase || "").startsWith("A") ? "Step A is still running." : "No climber."})); return; }
  const groups = new Map();
  for (const c of main) {
    const i = c.island == null ? 0 : c.island;
    if (!groups.has(i)) groups.set(i, []);
    groups.get(i).push(c);
  }
  const isl = islandRows(E);
  box.replaceChildren(...[...groups.keys()].sort((a, b) => a - b).map((i) => {
    const row = isl.find((r) => r.island === i) || {};
    return h("div", {class: "island-group", style: {"--isl": islandColor(i)}},
      h("div", {class: "island-head"}, h("i", {class: "dot"}), h("b", {text: `Island ${i}`}),
        h("span", {class: "muted num", text: `${int(row.capacity)} workers · best ${fix(row.best)} · ${groups.get(i).length} climber${groups.get(i).length === 1 ? "" : "s"}`})),
      h("div", {class: "cards"}, groups.get(i).map((c) => climberCard(c, finished, X, E))));
  }));
}

function climberCard(c, finished, X, E) {
  {
    const w = c.wave;
    const growth = Object.entries(c.growth || {});
    const best = c.found != null ? c.found : c.best;
    const i = c.island == null ? 0 : c.island;
    return h("div", {class: `card clickable ${climberClass(c)}`, role: "button", tabindex: "0",
        style: {borderTopColor: islandColor(i)}, title: "Open this climber's parts",
        onclick: () => openDetail(i, c.id), onkeydown: (e) => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); openDetail(i, c.id); } }},
      h("div", {class: "head"}, h("span", {class: "isl-tag", style: {background: islandColor(i)}, text: `I${i}`}), h("b", {text: `Climber ${c.id}`}),
        chip(finished && c.status === "climbing" ? "stopped" : STATUS_TEXT[c.status] || c.status, c.status === "climbing" && !finished ? "climbing" : c.status === "peak" ? "final" : "")),
      h("div", {class: "score-pair"},
        h("div", {},
          h("div", {class: "best-label", text: "best score"}),
          h("div", {class: "big", text: fix(best)})),
        h("div", {},
          h("div", {class: "best-label", text: "centre"}),
          h("div", {class: "big", text: fix(c.best)}),
          h("div", {class: "updated", text: bestUpdated(c, E)}),
          h("div", {class: "mono", text: (c.best_key || "").slice(0, 12)}))),
      w ? h("div", {class: "wave-line"}, h("span", {text: w.kind === "jumps" ? "random jumps" : `${w.kind} grid`}),
        h("span", {class: "num muted", text: `${int(w.back)} / ${int(w.sent)} back`})) : null,
      w ? bar(w.sent ? w.back / w.sent : 0, "", null) : null,
      kv([
        ["cores", climberCores(c, X)],
        ["grids", int(c.iterations)],
        ["configs sent", int(c.evaluated)],
        ["started at", mins(c.born_s)],
      ]),
      growth.length ? h("div", {class: "names"}, h("span", {class: "label", text: "growing steps "}),
        growth.map(([n, k]) => `${n} ${k > 0 ? "+" : ""}${k}`).join(" · ")) : null,
      partsLine(c),
      coreShare(c, E),
      h("div", {class: "origin", text: `From ${c.origin}. Click for grid, centre jumps, GA and RF-UCB.`}));
  }
}

// The share of the island's workers each part holds now: this climber's grid and jumps, the island's GA and RF-UCB.
function coreShare(c, E) {
  const i = c.island == null ? 0 : c.island;
  const row = islandRows(E).find((r) => r.island === i) || {};
  const cap = row.capacity || E.capacity || 0, P = c.parts || {}, sp = (row.sources || {}).split || {};
  if (!cap || !c.parts) return null;
  const vals = [["grid", (P.grid || {}).inflight || 0], ["jump", (P.jump || {}).inflight || 0], ["ga", sp.ga || 0], ["rf", sp.rf || 0]];
  return h("div", {class: "split"},
    h("div", {class: "best-label", text: "cores now (share of the island)"}),
    h("div", {class: "stack"}, vals.map(([k, v]) => h("span", {class: `seg ${k}`, style: {width: `${pct(v, cap)}%`}, title: `${PART_NAME[k]} ${int(v)}`}))),
    h("div", {class: "core-share"}, vals.map(([k, v]) => h("span", {class: `part-chip ${k}`},
      h("b", {text: PART_SHORT[k]}), ` ${pct(v, cap).toFixed(0)}%`, k === "ga" || k === "rf" ? h("span", {class: "muted", text: " isl"}) : null))));
}

function partsLine(c) {
  const P = c.parts;
  if (!P) return null;
  return h("div", {class: "parts-line"},
    ["grid", "jump", "ga", "rf"].map((k) => h("span", {class: `part-chip ${k}`, title: PART_NAME[k]},
      h("b", {text: PART_SHORT[k]}), ` ${P[k] && P[k].inflight != null ? int(P[k].inflight) + " live · " : ""}${int(P[k] ? P[k].moves : 0)} moves`)));
}

// ---- islands: 3D view, distances, and the climber detail ------------------------------------------------------------
const PART_NAME = {grid: "Grid", jump: "Centre jump", ga: "GA", rf: "RF-UCB"};
const PART_SHORT = {grid: "GRID", jump: "JUMP", ga: "GA", rf: "RF"};
const ISLAND_TOKENS = ["--accent", "--merge", "--final", "--bad", "--open", "--top-muted"];   // apart from the part colours
function islandColor(i) { return css(ISLAND_TOKENS[i % ISLAND_TOKENS.length]) || "#888"; }

// GA and RF-UCB runs are tagged G and R on every island, and their configs share batches: the executor's cores per
// tag, split between the islands by their configs in flight, is what each island's GA and RF-UCB really hold.
function realCores(rows, tag, key) {
  const X = (last && last.executor) || {}, cores = (X.cores || {})[tag];
  if (cores == null) return null;
  const n = rows.map((r) => (((r.sources || {}).split || {})[key]) || 0), tot = n.reduce((a, b) => a + b, 0);
  return rows.map((r, i) => (tot ? Math.round((cores * n[i]) / tot) : 0));
}
function withRealCores(rows) {
  const ga = realCores(rows, "G", "ga"), rf = realCores(rows, "R", "rf");
  return rows.map((r, i) => {
    const sp = (r.sources || {}).split;
    if (!sp || (ga == null && rf == null)) return r;
    return {...r, sources: {...r.sources, split: {...sp, ga: ga == null ? sp.ga : ga[i], rf: rf == null ? sp.rf : rf[i]}}};
  });
}

function islandRows(E) {
  if (!E) return [];
  if (Array.isArray(E.islands) && E.islands.length) return withRealCores(E.islands);
  return [{island: 0, capacity: E.capacity, best: E.best, phase: E.phase, best_vec: E.best_vec,
           sources: E.sources, cloud: E.cloud || [], climbers: climbers(E).length}];
}

// Principal axes of the encoded configs, recomputed at most every 15 s so the view does not spin with each poll.
let pcaBasis = null, pcaAt = 0, pcaDim = 0;
function pcaProject(vectors) {
  const V = vectors.filter((v) => Array.isArray(v) && v.length);
  if (!V.length) return null;
  const d = V[0].length;
  const mean = new Array(d).fill(0);
  for (const v of V) for (let j = 0; j < d; j++) mean[j] += v[j] / V.length;
  if (!pcaBasis || pcaDim !== d || Date.now() - pcaAt > 15000) {
    const C = Array.from({length: d}, () => new Array(d).fill(0));
    for (const v of V) for (let a = 0; a < d; a++) { const x = v[a] - mean[a]; for (let b = 0; b < d; b++) C[a][b] += x * (v[b] - mean[b]); }
    const basis = [];
    for (let k = 0; k < 3; k++) {
      let u = Array.from({length: d}, (_, j) => Math.sin(j * 1.7 + k * 3.1) + 1.1);
      for (let it = 0; it < 80; it++) {
        const w = C.map((row) => row.reduce((acc, x, j) => acc + x * u[j], 0));
        for (const b of basis) { const dot = w.reduce((acc, x, j) => acc + x * b[j], 0); for (let j = 0; j < d; j++) w[j] -= dot * b[j]; }
        const n = Math.hypot(...w) || 1;
        u = w.map((x) => x / n);
      }
      const big = u.reduce((m, x, j) => Math.abs(x) > Math.abs(u[m]) ? j : m, 0);
      if (u[big] < 0) u = u.map((x) => -x);              // a fixed sign, so the axes do not flip between refreshes
      basis.push(u);
    }
    pcaBasis = basis; pcaAt = Date.now(); pcaDim = d;
  }
  return (v) => (Array.isArray(v) && v.length === d ? pcaBasis.map((b) => b.reduce((acc, x, j) => acc + x * (v[j] - mean[j]), 0)) : null);
}

// ---- the 3D island view: canvas, no library (the page's CSP allows no eval) ----------------------------------------
// World points come from the configs' 25-setting encodings projected on 3 main axes (/api/points gives the axes and
// every scored config; until it answers, the axes come from the live clouds). Drag turns, the wheel zooms.
const v3 = {yaw: 0.7, pitch: 0.35, zoom: 1, auto: false, mode: "top", colour: "source", fit: "all",
            all: null, fetchedAt: 0, fetching: false, scene: null, screen: [], hover: null, raf: 0, wired: false, t0: performance.now()};
const SRC_TOKEN = {ga: "--warn", rf: "--rf", step_a: "--muted", released: "--faint", other: "--faint"};
const HEAT = [[44, 62, 145], [42, 157, 143], [233, 196, 106], [244, 162, 97], [230, 57, 70]];
function heat(t) {
  t = Math.max(0, Math.min(1, t)) * (HEAT.length - 1);
  const i = Math.min(HEAT.length - 2, Math.floor(t)), f = t - i, a = HEAT[i], b = HEAT[i + 1];
  return `rgb(${a.map((x, k) => Math.round(x + (b[k] - x) * f)).join(",")})`;
}
function rgba(col, a) {
  const c = col.trim();
  if (c.startsWith("#")) {
    const n = c.length === 4 ? c.slice(1).split("").map((x) => parseInt(x + x, 16)) : [1, 3, 5].map((k) => parseInt(c.slice(k, k + 2), 16));
    return `rgba(${n.join(",")},${a})`;
  }
  return c.startsWith("rgb(") ? c.replace("rgb(", "rgba(").replace(")", `,${a})`) : c;
}

async function fetchAllPoints() {
  if (v3.fetching || Date.now() - v3.fetchedAt < 60000) return;
  v3.fetching = true;
  try {
    const r = await api("/api/points");
    const d = await r.json();
    if (d && d.n) { v3.all = d; v3.scene = null; }
  } catch (e) { /* the next poll tries again */ } finally { v3.fetching = false; v3.fetchedAt = Date.now(); }
}

function projector(liveVectors) {
  const A = v3.all;
  // the server's axes, when the live vectors have its dimensions (island mode's banded spaces had fewer, and every
  // island point was dropped: an empty view); else the live points' own axes
  if (A && A.basis && A.mean && (!liveVectors.length || liveVectors.some((v) => Array.isArray(v) && v.length === A.mean.length))) {
    const d = A.mean.length;
    return (v) => (Array.isArray(v) && v.length === d ? A.basis.map((b) => b.reduce((acc, x, j) => acc + x * (v[j] - A.mean[j]), 0)) : null);
  }
  return pcaProject(liveVectors);
}

function settingDiff(a, b) {
  const names = (v3.all && v3.all.names) || [];
  if (!a || !b) return "";
  return a.map((x, j) => [j, b[j] - x]).filter(([, dx]) => Math.abs(dx) > 1e-9)
    .map(([j, dx]) => `${names[j] || `#${j}`} ${dx > 0 ? "+" : "−"}${Math.abs(dx * 100).toFixed(0)}%`).join(", ");
}

function buildScene(E) {
  const rows = islandRows(E), cl = climbers(E).filter((c) => c.parent == null);
  const live = [];
  for (const r of rows) { for (const [v] of r.cloud || []) live.push(v); if (r.best_vec) live.push(r.best_vec); }
  for (const c of cl) { if (c.vec) live.push(c.vec); for (const t of c.trail || []) live.push(t[0]); }
  const proj = projector(live);
  if (!proj) return null;
  const bg = [];
  if (v3.mode === "all" && v3.all && v3.all.n) {
    const A = v3.all;
    for (let i = 0; i < A.n; i++) bg.push({p: [A.xyz[3 * i], A.xyz[3 * i + 1], A.xyz[3 * i + 2]], s: A.score[i], src: A.sources[A.src[i]],
      island: A.island[i], climber: A.climber[i]});
  } else {
    for (const r of rows) for (const [v, sc] of r.cloud || []) { const p = proj(v); if (p) bg.push({p, s: sc, src: "climber", island: r.island}); }
  }
  const bests = rows.map((r) => ({p: proj(r.best_vec), island: r.island, best: r.best})).filter((b) => b.p);
  const paths = cl.map((c) => {
    const i = c.island == null ? 0 : c.island;
    const steps = (c.trail || []).map((t) => ({p: proj(t[0]), v: t[0], s: t[1], part: t[2]})).filter((x) => x.p);
    const here = proj(c.vec);
    return {c, island: i, steps, here};
  });
  const links = (E.island_distances || []).map((dd) => {
    const a = bests.find((b) => b.island === dd.a), b = bests.find((x) => x.island === dd.b);
    return a && b ? {a: a.p, b: b.p, dd} : null;
  }).filter(Boolean);
  const scores = bg.map((q) => q.s).sort((a, b) => a - b);
  const lo = scores.length ? scores[Math.floor(scores.length * 0.05)] : 0, hi = scores.length ? scores[scores.length - 1] : 1;
  // colour by rank (percentile among the shown points), in 16 steps: scores bunch near the top (round 08:34: 75% within
  // 0.44 of the best), so a linear scale painted almost every point orange-red
  const order = bg.map((q, i) => i).sort((a, b) => bg[a].s - bg[b].s);
  order.forEach((i, r) => { bg[i].t = order.length > 1 ? Math.round((r / (order.length - 1)) * 15) / 15 : 1; });
  return {bg, bests, paths, links, lo, hi};
}

function bgTip() {               // shared hover text of a background point (this = the hovered entry)
  const q = this.q;
  return `${q.src === "climber" ? `climber ${q.climber >= 0 ? q.climber : ""} · island ${q.island}` : q.src.replace("_", " ")} · ${fix(q.s)}`;
}

function draw3d(now) {
  const box = $("space3d"), sc = v3.scene;
  if (!sc) return;
  let cv = box.querySelector("canvas");
  if (!cv) {
    box.replaceChildren();
    cv = h("canvas", {class: "space-canvas"});
    box.append(cv, h("div", {class: "space-tip", id: "space-tip", hidden: true}), spaceControls());
  }
  const dpr = Math.min(1.5, window.devicePixelRatio || 1), W = Math.max(320, box.clientWidth || 600), H = 480;
  if (cv.width !== Math.round(W * dpr) || cv.height !== Math.round(H * dpr)) { cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr); cv.style.width = `${W}px`; cv.style.height = `${H}px`; }
  const g = cv.getContext("2d");
  g.setTransform(dpr, 0, 0, dpr, 0, 0);
  // the frame: what "fit" shows
  let fitPts = sc.bg.map((q) => q.p).concat(sc.bests.map((b) => b.p));
  for (const pa of sc.paths) for (const st of pa.steps) fitPts.push(st.p);
  if (v3.fit !== "all") {
    const pa = sc.paths.find((x) => `${x.island}:${x.c.id}` === v3.fit);
    if (pa && pa.steps.length) fitPts = pa.steps.map((st) => st.p).concat(pa.here ? [pa.here] : []);
  }
  if (!fitPts.length) return;
  const c = [0, 1, 2].map((k) => fitPts.reduce((a, p) => a + p[k], 0) / fitPts.length);
  const rad = fitPts.reduce((m, p) => Math.max(m, Math.hypot(p[0] - c[0], p[1] - c[1], p[2] - c[2])), 1e-6);   // no spread: 45k points
  const cy = Math.cos(v3.yaw), sy = Math.sin(v3.yaw), cp = Math.cos(v3.pitch), sp = Math.sin(v3.pitch);
  const scale = 0.40 * Math.min(W, H) * v3.zoom;
  const P = (p) => {
    const x0 = (p[0] - c[0]) / rad, y0 = (p[1] - c[1]) / rad, z0 = (p[2] - c[2]) / rad;
    const x1 = cy * x0 + sy * z0, z1 = -sy * x0 + cy * z0, y2 = cp * y0 - sp * z1, z2 = sp * y0 + cp * z1;
    const f = 3.2 / (3.2 + z2);
    return [W / 2 + scale * x1 * f, H / 2 - scale * y2 * f, z2, f];
  };
  const W0 = (q) => [c[0] + q[0] * rad, c[1] + q[1] * rad, c[2] + q[2] * rad];      // unit frame -> world
  const key = [v3.yaw, v3.pitch, v3.zoom, v3.fit, v3.mode, v3.colour, W, H, v3.sceneId, effectiveTheme()].join("|");
  if (v3.bgKey !== key) {                       // the background layer is redrawn only when the view or data moves
    const off = v3.bgCanvas || (v3.bgCanvas = document.createElement("canvas"));
    off.width = cv.width; off.height = cv.height;
    const bgc = off.getContext("2d");
    bgc.setTransform(dpr, 0, 0, dpr, 0, 0);
    const bgScreen = [];
    // background and floor
    const panel = css("--panel"), panel2 = css("--panel-2"), ink = css("--ink"), mut = css("--muted"), line = css("--line");
    const bgGrad = bgc.createRadialGradient(W / 2, H * 0.42, 10, W / 2, H / 2, Math.max(W, H) * 0.75);
    bgGrad.addColorStop(0, panel); bgGrad.addColorStop(1, panel2);
    bgc.fillStyle = bgGrad; bgc.fillRect(0, 0, W, H);
    bgc.lineWidth = 1;
    for (let k = -5; k <= 5; k++) {
      const t = k / 5, al = 0.35 - Math.abs(t) * 0.2;
      bgc.strokeStyle = rgba(line, al);
      let a = P(W0([t * 1.2, -1.05, -1.2])), b = P(W0([t * 1.2, -1.05, 1.2]));
      bgc.beginPath(); bgc.moveTo(a[0], a[1]); bgc.lineTo(b[0], b[1]); bgc.stroke();
      a = P(W0([-1.2, -1.05, t * 1.2])); b = P(W0([1.2, -1.05, t * 1.2]));
      bgc.beginPath(); bgc.moveTo(a[0], a[1]); bgc.lineTo(b[0], b[1]); bgc.stroke();
    }
    // axis gizmo
    const gz = [[1, 0, 0, "1"], [0, 1, 0, "2"], [0, 0, 1, "3"]], o = [44, H - 40];
    for (const [ax, ay, az, label] of gz) {
      const x1 = cy * ax + sy * az, z1 = -sy * ax + cy * az, y2 = cp * ay - sp * z1;
      bgc.strokeStyle = mut; bgc.beginPath(); bgc.moveTo(o[0], o[1]); bgc.lineTo(o[0] + 22 * x1, o[1] - 22 * y2); bgc.stroke();
      bgc.fillStyle = mut; bgc.font = "10px system-ui"; bgc.fillText(`axis ${label}`, o[0] + 26 * x1 - 8, o[1] - 26 * y2 + 3);
    }
    // background points, far first, batched: one fill per colour and depth step (2,500 single arcs per frame made the
    // auto-rotation stutter)
    const all = v3.mode === "all";
    const batches = new Map();
    for (const q of sc.bg) {
      const sp2 = P(q.p);
      const depth = Math.max(0.15, Math.min(1, (1.6 - sp2[2]) / 2.2));
      const dk = Math.round(depth * 4);
      if (q._col == null || q._colKey !== v3.colour) {
        q._colKey = v3.colour;
        q._col = v3.colour === "score" ? heat(q.t == null ? 1 : q.t)
          : q.src === "climber" && q.island >= 0 ? islandColor(q.island) : css(SRC_TOKEN[q.src] || "--faint");
      }
      const k = `${dk}|${q._col}`;
      let b = batches.get(k);
      if (!b) batches.set(k, (b = {dk, col: q._col, pts: []}));
      const r = (all ? 1.6 : 2.6) * sp2[3];
      b.pts.push(sp2[0], sp2[1], r);
      if (!v3.auto) bgScreen.push({x: sp2[0], y: sp2[1], r: Math.max(3, r + 1), q, tip: bgTip, w: 0});
    }
    for (const b of [...batches.values()].sort((x, y) => x.dk - y.dk)) {
      bgc.fillStyle = rgba(b.col, (all ? 0.55 : 0.6) * Math.max(0.15, b.dk / 4));
      bgc.beginPath();
      for (let i = 0; i < b.pts.length; i += 3) { const r = b.pts[i + 2]; bgc.rect(b.pts[i] - r, b.pts[i + 1] - r, 2 * r, 2 * r); }
      bgc.fill();
    }
    v3.bgKey = key; v3.bgScreen = bgScreen;
  }
  g.setTransform(1, 0, 0, 1, 0, 0);
  g.drawImage(v3.bgCanvas, 0, 0);
  g.setTransform(dpr, 0, 0, dpr, 0, 0);
  const screen = v3.bgScreen.slice();
  const panel = css("--panel"), ink = css("--ink"), mut = css("--muted");
  // island links
  for (const L of sc.links) {
    const a = P(L.a), b = P(L.b);
    g.setLineDash([6, 5]); g.strokeStyle = rgba(mut, 0.9); g.lineWidth = 1.5;
    g.beginPath(); g.moveTo(a[0], a[1]); g.lineTo(b[0], b[1]); g.stroke(); g.setLineDash([]);
    g.fillStyle = mut; g.font = "11px system-ui";
    g.fillText(`${fix(L.dd.distance, 2)} apart · ${L.dd.settings} settings`, (a[0] + b[0]) / 2 + 6, (a[1] + b[1]) / 2 - 6);
  }
  // climber paths with numbered steps
  const pulse = 0.5 + 0.5 * Math.sin(((now || performance.now()) - v3.t0) / 380);
  for (const pa of sc.paths) {
    const col = islandColor(pa.island), seq = pa.steps.concat(pa.here ? [{p: pa.here, v: pa.c.vec, s: pa.c.best, part: "now"}] : []);
    const sps = seq.map((st) => P(st.p));
    for (let k = 1; k < sps.length; k++) {
      g.strokeStyle = rgba(col, 0.25 + 0.65 * (k / sps.length)); g.lineWidth = 1.5 + 1.5 * (k / sps.length);
      g.beginPath(); g.moveTo(sps[k - 1][0], sps[k - 1][1]); g.lineTo(sps[k][0], sps[k][1]); g.stroke();
    }
    seq.forEach((st, k) => {
      const q = sps[k];
      if (st.part === "now") return;
      const pc = st.part === "grid" ? css("--good") : st.part === "jump" ? css("--scout") : st.part === "ga" ? css("--warn") : st.part === "rf" ? css("--rf") : col;
      g.fillStyle = pc; g.beginPath(); g.arc(q[0], q[1], 3 * q[3], 0, 6.2832); g.fill();
      if (v3.fit !== "all" || sps.length <= 12) { g.fillStyle = mut; g.font = "9px system-ui"; g.fillText(String(k + 1), q[0] + 4, q[1] - 4); }
      const prev = k ? seq[k - 1].v : null;
      screen.push({x: q[0], y: q[1], r: 6, w: 2, tip: () => `climber ${pa.c.id} · step ${k + 1} · ${PART_NAME[st.part] || st.part} · ${fix(st.s)}${prev ? `\n${settingDiff(prev, st.v)}` : ""}`});
    });
    if (pa.here) {
      const q = P(pa.here), climbing = pa.c.status === "climbing", r = (climbing ? 7 : 5.5) * q[3];
      if (climbing) { g.strokeStyle = rgba(col, 0.55 * (1 - pulse)); g.lineWidth = 2; g.beginPath(); g.arc(q[0], q[1], r + 4 + 8 * pulse, 0, 6.2832); g.stroke(); }
      g.shadowColor = col; g.shadowBlur = climbing ? 14 : 6;
      const grad = g.createRadialGradient(q[0] - r / 3, q[1] - r / 3, 1, q[0], q[1], r);
      grad.addColorStop(0, "#ffffff"); grad.addColorStop(0.35, col); grad.addColorStop(1, rgba(col, 0.85));
      g.fillStyle = grad; g.beginPath(); g.arc(q[0], q[1], r, 0, 6.2832); g.fill(); g.shadowBlur = 0;
      g.fillStyle = col; g.font = "bold 11px system-ui"; g.fillText(`c${pa.c.id}`, q[0] + r + 4, q[1] + 4);
      screen.push({x: q[0], y: q[1], r: r + 4, w: 3, tip: () => `climber ${pa.c.id} · island ${pa.island} · ${pa.c.status} · centre ${fix(pa.c.best)} · ${seq.length - 1} moves`});
    }
  }
  // island bests
  for (const b of sc.bests) {
    const q = P(b.p), col = islandColor(b.island), r = 8 * q[3];
    g.shadowColor = col; g.shadowBlur = 18;
    g.fillStyle = col; g.beginPath(); g.moveTo(q[0], q[1] - r); g.lineTo(q[0] + r, q[1]); g.lineTo(q[0], q[1] + r); g.lineTo(q[0] - r, q[1]); g.closePath(); g.fill();
    g.shadowBlur = 0; g.strokeStyle = ink; g.lineWidth = 1; g.stroke();
    g.fillStyle = col; g.font = "bold 12px system-ui"; g.fillText(`I${b.island} ${fix(b.best, 3)}`, q[0] + r + 5, q[1] - r);
    screen.push({x: q[0], y: q[1], r: r + 4, w: 4, tip: () => `island ${b.island} best ${fix(b.best)}`});
  }
  // legend
  g.font = "11px system-ui"; g.fillStyle = mut;
  const legend = v3.colour === "score" ? `colour: score, ${fix(sc.lo, 2)} (blue) → ${fix(sc.hi, 2)} (red)`
    : "colour: island (climbers) · amber GA · magenta RF-UCB · grey Step A";
  g.fillText(`${legend} · ${v3.mode === "all" ? `${int(sc.bg.length)} scored configs` : "top 40 per island"}`, 12, 18);
  v3.screen = screen;
}

function spaceControls() {
  const btn = (id, text, on) => h("button", {class: "btn small", id, type: "button", onclick: on, text});
  const sel = h("select", {class: "small-select", id: "space-fit", onchange: (e) => { v3.fit = e.target.value; draw3d(); }});
  return h("div", {class: "space-ctl"},
    btn("space-mode", "All scored", () => { v3.mode = v3.mode === "all" ? "top" : "all"; v3.scene = null; if (v3.mode === "all") { v3.fetchedAt = 0; fetchAllPoints().then(() => renderIslands(last && last.search && last.search.engine)); } renderIslands(last && last.search && last.search.engine); }),
    btn("space-colour", "Colour: score", () => { v3.colour = v3.colour === "score" ? "source" : "score"; draw3d(); syncControls(); }),
    btn("space-auto", "Auto-rotate", () => { v3.auto = !v3.auto; v3.bgKey = null; draw3d(); syncControls(); }),
    btn("space-reset", "Reset", () => { Object.assign(v3, {yaw: 0.7, pitch: 0.35, zoom: 1, fit: "all"}); draw3d(); syncControls(); }),
    sel);
}

function syncControls() {
  const n = v3.all && v3.all.n;
  const m = $("space-mode"); if (m) m.textContent = v3.mode === "all" ? "Top 40 per island" : `All scored${n ? ` (${int(n)})` : ""}`;
  const c = $("space-colour"); if (c) c.textContent = v3.colour === "score" ? "Colour: source" : "Colour: score";
  const a = $("space-auto"); if (a) a.classList.toggle("primary", v3.auto);
  const sel = $("space-fit");
  if (sel && v3.scene) {
    const opts = [["all", "Fit: everything"]].concat(v3.scene.paths.map((pa) => [`${pa.island}:${pa.c.id}`, `Fit: climber ${pa.c.id}`]));
    const key = opts.map((o) => o[0]).join("|");
    if (sel.dataset.key !== key) { sel.replaceChildren(...opts.map(([v, t]) => h("option", {value: v, text: t}))); sel.dataset.key = key; }
    sel.value = opts.some((o) => o[0] === v3.fit) ? v3.fit : "all";
  }
}

function wire3d() {
  if (v3.wired) return;
  v3.wired = true;
  const box = $("space3d");
  let drag = null;
  box.addEventListener("pointerdown", (e) => { if (e.target.closest(".space-ctl")) return; drag = [e.clientX, e.clientY]; box.setPointerCapture(e.pointerId); });
  box.addEventListener("pointermove", (e) => {
    if (drag) {
      v3.yaw += (e.clientX - drag[0]) * 0.01;
      v3.pitch = Math.max(-1.5, Math.min(1.5, v3.pitch + (e.clientY - drag[1]) * 0.01));
      drag = [e.clientX, e.clientY];
      draw3d();
      return;
    }
    const rect = box.getBoundingClientRect(), x = e.clientX - rect.left, y = e.clientY - rect.top;
    let best = null;
    for (const s2 of v3.screen) {
      const d = Math.hypot(s2.x - x, s2.y - y);
      if (d <= s2.r && (!best || s2.w > best.w || (s2.w === best.w && d < best.d))) best = {...s2, d};
    }
    const tip = $("space-tip");
    if (!tip) return;
    if (!best) { tip.hidden = true; return; }
    tip.textContent = best.tip();
    tip.hidden = false;
    tip.style.left = `${Math.min(rect.width - 260, x + 14)}px`;
    tip.style.top = `${y + 14}px`;
  });
  box.addEventListener("pointerup", () => { drag = null; });
  box.addEventListener("pointerleave", () => { const tip = $("space-tip"); if (tip) tip.hidden = true; });
  box.addEventListener("wheel", (e) => { e.preventDefault(); v3.zoom = Math.max(0.3, Math.min(8, v3.zoom * (e.deltaY < 0 ? 1.12 : 0.89))); draw3d(); }, {passive: false});
  // auto-rotation: every frame, turning at a fixed speed whatever the frame rate (the batched points keep a frame
  // cheap); climbers moving without auto-rotation: at most 4 redraws a second. A drag or the wheel redraws at once.
  const tick = (now) => {
    if (v3.scene && !document.hidden && $("view-live") && !$("view-live").hidden) {
      const dt = Math.min(0.1, (now - (v3.lastFrame || now)) / 1000);
      if (v3.auto) { v3.yaw += 0.35 * dt; draw3d(now); }
      else if (now - (v3.lastTick || 0) >= 250 && v3.scene.paths.some((pa) => pa.c.status === "climbing")) {
        v3.lastTick = now; draw3d(now);
      }
    }
    v3.lastFrame = now;
    v3.raf = requestAnimationFrame(tick);
  };
  v3.raf = requestAnimationFrame(tick);
}

function renderIslands(E) {
  const box = $("space3d"), dist = $("island-dist");
  const rows = islandRows(E);
  if (!E || !rows.length) { dist.replaceChildren(h("div", {class: "empty", text: "Waiting for the engine."})); return; }
  const A = v3.all, sampled = A && A.n_total > A.n ? ` · ${int(A.n)} of ${int(A.n_total)} points (the best 1,000 and an even sample)` : "";
  $("islands-sub").textContent = `${rows.length} island${rows.length === 1 ? "" : "s"} · the configs' 25 settings on their 3 main axes${sampled} · drag to turn, scroll to zoom, hover for details`;
  fetchAllPoints();
  v3.scene = buildScene(E);
  v3.sceneId = (v3.sceneId || 0) + 1;
  if (!v3.scene) {
    box.replaceChildren(h("div", {class: "empty", text: "No positions yet: they appear once the engine publishes scored configs."}));
  } else {
    wire3d();
    draw3d();
    syncControls();
  }
  const pairs = (E.island_distances || []).map((dd) => h("div", {class: "rowline"},
    h("div", {class: "top-line"}, h("span", {}, h("i", {class: "dot", style: {background: islandColor(dd.a)}}), ` Island ${dd.a} ↔ `,
      h("i", {class: "dot", style: {background: islandColor(dd.b)}}), ` Island ${dd.b}`),
      h("span", {class: "num", text: `${fix(dd.distance, 3)} apart · ${dd.settings} of 25 settings differ`})),
    bar(Math.min(1, dd.distance / 2), "", null)));
  const summary = rows.map((r) => {
    const S = r.sources || {}, ga = S.ga || {}, rf = S.rf || {};
    return h("div", {class: "rowline"},
      h("div", {class: "top-line"}, h("span", {}, h("i", {class: "dot", style: {background: islandColor(r.island)}}), h("b", {text: ` Island ${r.island}`}),
        h("span", {class: "muted", text: ` · ${int(r.capacity)} workers · ${r.climbers != null ? r.climbers : "–"} climbers`})),
        h("span", {class: "num", text: `best ${fix(r.best)}`})),
      h("div", {class: "muted num small", text: `GA ${int(ga.inflight)} live / ${int(ga.sent)} sent, best ${fix(ga.best)} · RF-UCB ${int(rf.inflight)} live / ${int(rf.sent)} sent, best ${fix(rf.best)}`}),
      splitBar(S.split));
  });
  dist.replaceChildren(...summary, ...(pairs.length ? pairs : [h("div", {class: "empty", text: rows.length > 1 ? "Distances appear once every island has a best config." : "One island: no distance to show."})]));
}

// Runs in flight on one island by part, against the grid and centre-jump quotas (climb.grid_cores, climb.jump_cores).
function splitBar(sp) {
  if (!sp || !sp.cap) return null;
  const parts = [["grid", sp.grid], ["jump", sp.jump], ["rf", sp.rf], ["ga", sp.ga]];
  const used = parts.reduce((a, [, v]) => a + (v || 0), 0);
  return h("div", {class: "split"},
    h("div", {class: "stack"}, parts.map(([k, v]) => h("span", {class: `seg ${k}`, style: {width: `${pct(v || 0, Math.max(sp.cap, used))}%`},
      title: `${PART_NAME[k]} ${int(v)}`}))),
    h("div", {class: "muted num small", text: `in flight: grid ${int(sp.grid)}${sp.grid_quota ? ` / ${int(sp.grid_quota)}` : ""} · jumps ${int(sp.jump)}${sp.jump_quota ? ` / ${int(sp.jump_quota)}` : ""} · RF ${int(sp.rf)} · GA ${int(sp.ga)} · of ${int(sp.cap)} workers`}));
}

let detail = null;                // {island, id} of the open climber card
function openDetail(island, id) { detail = {island, id}; $("drawer").hidden = false; render(); $("drawer-close").focus(); }
function closeDetail() { detail = null; $("drawer").hidden = true; }

function partPanel(key, title, sub, pairs, extra) {
  return h("section", {class: `part ${key}`}, h("h3", {}, h("span", {class: `part-chip ${key}`}, h("b", {text: PART_SHORT[key]})), ` ${title}`),
    sub ? h("div", {class: "muted small", text: sub}) : null, kv(pairs), extra || null);
}

function renderDetail(E, finished) {
  if (!detail) return;
  const c = climbers(E).find((x) => x.id === detail.id && (x.island == null ? 0 : x.island) === detail.island);
  if (!c) { $("drawer-title").textContent = `Climber ${detail.id} has left the round`; $("drawer-parts").replaceChildren(); return; }
  const i = c.island == null ? 0 : c.island;
  const row = islandRows(E).find((r) => r.island === i) || {};
  const S = row.sources || {}, ga = S.ga || {}, rf = S.rf || {}, P = c.parts || {};
  const p = (k) => P[k] || {moves: 0, gain: 0};
  const lastMove = (k) => (p(k).last_s == null ? "–" : `${mins(p(k).last_s)} into the round`);
  $("drawer-title").replaceChildren(h("span", {class: "isl-tag", style: {background: islandColor(i)}, text: `I${i}`}),
    ` Climber ${c.id} · ${finished && c.status === "climbing" ? "stopped" : STATUS_TEXT[c.status] || c.status} · centre ${fix(c.best)} · ${int(c.evaluated)} configs sent`);
  const w = c.wave;
  const steps = Object.entries(c.steps || {});
  const waves = (c.waves || []).slice(-6).reverse();
  const arms = Object.entries((E.leaps || {}).by_setting || {}).sort((a, b) => (b[1].wins || 0) - (a[1].wins || 0)).slice(0, 6);
  $("drawer-parts").replaceChildren(
    partPanel("grid", "Grid", w && w.kind !== "jumps" ? `${w.kind} wave: ${int(w.back)} / ${int(w.sent)} back` : "no grid wave out now", [
      ["in flight now", int(p("grid").inflight)], ["runs sent", int(p("grid").sent)], ["centre moves", int(p("grid").moves)],
      ["gain", `+${fix(p("grid").gain)}`], ["last move", lastMove("grid")], ["grid waves", int(c.iterations)], ["flat settings", `${int(c.flat)} of ${int(c.expensive)}`]],
      h("div", {}, steps.length ? h("div", {class: "names"}, h("span", {class: "label", text: "steps now "}),
          steps.map(([n, k]) => `${n} ±${k.join("/")}`).join(" · ")) : null,
        waves.length ? h("div", {class: "names"}, h("span", {class: "label", text: "last waves "}),
          waves.map((x) => `${x[0]} ${x[1]}/${x[2]} in ${x[3]}m`).join(" · ")) : null)),
    partPanel("jump", "Centre jump", w && w.kind === "jumps" ? `jump wave: ${int(w.back)} / ${int(w.sent)} back` : "jumps ride with every grid wave", [
      ["in flight now", int(p("jump").inflight)], ["runs sent", int(p("jump").sent)], ["centre moves", int(p("jump").moves)],
      ["gain", `+${fix(p("jump").gain)}`], ["last move", lastMove("jump")]],
      arms.length ? h("div", {class: "names"}, h("span", {class: "label", text: "island's best arms "}),
        arms.map(([n, a]) => `${n} ${int(a.wins)}/${int(a.tries)}`).join(" · ")) : null),
    partPanel("ga", "GA", `island ${i}: ${int(ga.inflight)} in flight, ${int(ga.back)} / ${int(ga.sent)} back, best ${fix(ga.best)}`, [
      ["moves into this climber", int(p("ga").moves)], ["gain", `+${fix(p("ga").gain)}`], ["last move", lastMove("ga")]]),
    partPanel("rf", "RF-UCB", `island ${i}: ${int(rf.inflight)} in flight, ${int(rf.back)} / ${int(rf.sent)} back, best ${fix(rf.best)}`, [
      ["moves into this climber", int(p("rf").moves)], ["gain", `+${fix(p("rf").gain)}`], ["last move", lastMove("rf")],
      ["forests fitted", int(rf.fits)], ["trained on", int(rf.trained_on)], ["share · kappa", `${fix(rf.share, 2)} · ${fix(rf.kappa, 2)}`]]));
}

function climberCores(c, X) {
  const table = X && X.cores;
  if (table) {
    const n = table[`B${c.id}`];
    return n == null ? "0" : int(n);
  }
  const w = c.wave;
  if (w && w.sent != null) return int(Math.max(0, (w.sent || 0) - (w.back || 0)));
  return "0";
}

function gaNote(E) {
  return (E.version || "").startsWith("v9.34")
    ? "The compass keeps running. The other cores run the pile GA: half crosses inside the piles and may land on the climber, half looks for a fingerprint 0.10 away. A linked group stays with one parent. The highest score becomes the fingerprint, and a pile it overlaps is merged into it. Sixteen cores stay free so the compass starts on its own heads."
    : (E.version || "").startsWith("v9.34") || (E.version || "").startsWith("v9.33")
    ? "One climber compasses until the fine-tune. Sixteen cores stay free so the compass starts on its own heads. The population mates on the other cores, and a higher score within 0.10 becomes the next centre when the compass returns. A linked group stays with one parent."
    : (E.version || "").startsWith("v9.32")
    ? "One climber compasses the highest score until the fine-tune. The hill lane crosses two piles. A verified linkage group stays with one parent. An unlinked knob comes from either. A child is sent only when it lands 0.10 away. For the last 12 minutes the compass also tries one rung further."
    : (E.version || "").startsWith("v9.31")
    ? "One climber compasses the highest score until the fine-tune. A face within 0.05 is not live. Two live faces send their joint cube. The hill lane crosses two piles and keeps a child only when it lands 0.10 away. For the last 12 minutes the compass also tries one rung further."
    : (E.version || "").startsWith("v9.30")
    ? "One climber compasses the highest score until the fine-tune. A face within 0.05 is not live. Two live faces send their joint cube, and the centre moves when that grid is back. For the last 12 minutes the compass also tries one rung further; that rung stays out of the cube. Free cores look for a fingerprint 0.10 away, then probe settings the compass has not tried. The fine-tune is the same compass and does not look for a new fingerprint."
    : (E.version || "").startsWith("v9.29")
    ? "One climber, on the highest fingerprint. The other piles wait. Step A places a dummy on each apparent peak and skips generation 1. The grid is a compass: one rung either side. The centre moves when that compass finishes. The fine-tune climbs the best."
    : (E.version || "").startsWith("v9.28")
    ? "Six piles. Four climb until the fine-tune. While a grid is walking, a new high in the pile does not replace the fingerprint or the centre. The grid writes the setting it found."
    : (E.version || "").startsWith("v9.27")
    ? "Six piles. Four climb until the fine-tune. Elite probes of the best go out together. A worker the planned grid did not fill walks further up the same line."
    : (E.version || "").startsWith("v9.26")
    ? "Six piles. Four climb until the fine-tune. While the GA owns the fleet, every child copies the best config and changes one to three settings. That window runs once. Then the grid keeps the fleet and the GA does not breed. A peaked centre is not taken."
    : (E.version || "").startsWith("v9.25")
    ? "Six piles. Four climb until the fine-tune. While a GA child or a candidate is still moving a climber, the climbers send nothing and the fleet breeds. Once that goes quiet, one rung past the compass walks to the peak. The line does not widen. The fine-tune climbs the best config, restarts when a higher score arrives, and does not breed."
    : (E.version || "").startsWith("v9.24")
    ? "Six piles. Four climb until the fine-tune. A config at least 0.10 from every fingerprint, and above the lowest climbing pile, takes that climber now. It does not move onto a centre another climber already holds. The fine-tune climbs the best config, restarts when a higher score arrives, and does not breed."
    : (E.version || "").startsWith("v9.23")
    ? "Six piles. Four climb until the fine-tune. A peaked seat goes to the best candidate. A candidate that scores higher than the lowest climbing pile takes that climber. The first half of Step B keeps a short grid. The second half takes five rungs along each side. The fine-tune climbs the best config and does not breed."
    : (E.version || "").startsWith("v9.22")
    ? "Six piles. Four climb. Two candidates wait, and the best one takes a seat when a climber reaches the top. The highest score within 0.10 of a pile becomes that climber's centre, even when it cannot become the fingerprint."
    : (E.version || "").startsWith("v9.21")
    ? "Up to four piles, one climber each. The highest score within 0.10 of a pile becomes its reference, and the climber develops from it. A fingerprint within 0.10 of the history is not explored again."
    : (E.version || "").startsWith("v9.20")
    ? "Up to four piles, one climber each. A pile breeds from up to 32 parents. Every fingerprint a pile stands on is kept, and a fingerprint within 0.10 of that history is not explored. A run that comes back above the centre takes the centre."
    : (E.version || "").startsWith("v9.19")
    ? "Up to four piles, one climber each. A pile is a callset fingerprint. Piles stay 0.10 apart. Half of leftover GA breeds inside the piles. Half looks for a new fingerprint. A run that comes back higher takes the centre."
    : (E.version || "").startsWith("v9.18")
    ? "Four climbers. Compass then a cube. Leftover cores run GA on the same pass. One far seat stays on a far hill: a near-hill GA goes to a near climber. Half of leftover children use the current best as a parent."
    : (E.version || "").startsWith("v9.17")
    ? "Four climbers. A grid is a compass of each setting at plus or minus one step, then a cube only on the axes that beat the centre. One seat is the best far hill. Half of leftover GA children use the current best as a parent."
    : (E.version || "").startsWith("v9.16")
    ? "Four climbers. One seat is the best far hill even when it is more than 0.10 behind. Scouts for that seat start from the far parents."
    : (E.version || "").startsWith("v9.15")
    ? "Four climbers. A GA child is not sent on the same hill as another GA child or a climber. Step B opens with the whole fleet on the GA. The lowest takes a free GA score."
    : (E.version || "").startsWith("v9.14")
    ? "Four climbers. Step B opens with the whole fleet on the GA for up to 8 minutes. Grids start once a free GA beats the lowest climber. The lowest takes the highest GA score that is not already a climber."
    : (E.version || "").startsWith("v9.13")
    ? "Four climbers. The 60 parents are far apart. The lowest climber takes the highest GA score that is not already a climber. Step B mutation starts at 0.3 and rises to 0.7."
    : (E.version || "").startsWith("v9.12")
    ? "Four climbers. The 60 parents are far apart. When the highest GA score is above the lowest climber's best, that climber takes it. Step B mutation starts at 0.3 and rises to 0.7."
    : (E.version || "").startsWith("v9.11")
    ? "Four climbers. The 60 parents are far apart. When a read or assembly grid ends, the lowest climber takes the highest GA score if it is behind. Step B mutation starts at 0.3 and rises to 0.7."
    : (E.version || "").startsWith("v9.10")
    ? "The 60 parents are far apart. When a read or assembly grid ends, the lowest climber takes the highest GA score if it is behind. Step B mutation starts at 0.3 and rises to 0.7."
    : (E.version || "").startsWith("v9.9")
    ? "The 60 parents are far apart. A GA best replaces the climber whose own runs have the lowest best. Step B mutation starts at 0.3 and rises to 0.7."
    : (E.version || "").startsWith("v9.8")
    ? "The 60 parents are far apart. Step B mutation starts at 0.3 and rises to 0.7. A child that differs in fewer than four settings from a climber is not sent."
    : (E.version || "").startsWith("v9.7")
    ? "The 60 parents are far apart. An empty seat is filled from scored configs that are far from the climbers."
    : (E.version || "").startsWith("v9.6")
    ? "After the one-rung probes, breeding continues. A child closer than the seat distance to a climber is not sent."
    : (E.version || "").startsWith("v9.5")
    ? "A child is not sent when every value it changes has already tied that score."
    : (E.version || "").startsWith("v9.4")
    ? "On a plateau, one-rung probes replace children. A child within two settings of a climber is not sent."
    : (E.version || "").startsWith("v9.3")
    ? "Lowest priority. When it holds the best score, the lowest climber takes that config."
    : "Lowest priority. Sent only while every ready candidate update is already out.";
}

function pilesOn(E) {
  // Callset piles start at v9.19. v9.16 seats four climbers and does not draw a pile board.
  const n = /^v9\.(\d+)/.exec((E && E.version) || "");
  return !!(n && Number(n[1]) >= 19);
}

function renderStepB(E, X) {
  const box = $("stepb");
  if (!box) return;
  if (!E) { box.replaceChildren(h("div", {class: "empty", text: "Waiting for the engine."})); return; }
  const ga = E.ga || {};
  const piles = pilesOn(E) ? (E.piles || []) : [];
  if (!ga.sent && !piles.length) {
    box.replaceChildren(h("div", {class: "empty", text: (E.phase || "").startsWith("A") ? "Step A is still running." : "The GA has not sent a child yet."}));
    return;
  }
  const heldKeys = new Set(climbers(E).map((c) => c.best_key).filter(Boolean));
  const heldScores = new Set(climbers(E).flatMap((c) => [c.found, c.best]).filter((s) => s != null).map((s) => +Number(s).toFixed(4)));
  const top = (ga.top || []).filter((row) => !heldKeys.has(row.key) && !heldScores.has(+Number(row.score).toFixed(4))).slice(0, 3);
  const gaCores = X && X.cores && X.cores.G != null ? X.cores.G : 0;
  box.replaceChildren(h("div", {class: "rowline"},
    h("div", {class: "top-line"}, h("b", {text: "GA on idle workers"}),
      h("span", {class: "num muted", text: `round best ${fix(E.best)} · ${int(ga.back)} / ${int(ga.sent)} back · ${int(gaCores)} cores`})),
    ga.sent ? bar(ga.sent ? ga.back / ga.sent : 0, "scout", null) : null,
    pilesOn(E) ? pileBoard(piles) : null,
    top.length ? h("div", {class: "ga-top"},
      top.map((row, i) => h("div", {class: "ga-top-item"},
        h("div", {class: "best-label", text: i === 0 ? "top" : `#${i + 1}`}),
        h("div", {class: "big", text: fix(row.score)}),
        h("div", {class: "mono", text: (row.key || "").slice(0, 12)}))))
      : (ga.sent ? h("div", {class: "empty", text: "No scored GA child besides the climbers."}) : null),
    h("div", {class: "muted", text: gaNote(E)})));
}

function pileBoard(piles) {
  if (!piles.length) return h("div", {class: "muted", text: "No pile yet."});
  return h("div", {class: "piles"}, piles.map((p) => h("div", {class: "pile"},
    h("div", {class: "best-label", text: p.candidate ? `Pile ${p.pile} · candidate` : `Pile ${p.pile}`}),
    h("div", {class: "ga-top"}, (p.top || []).slice(0, 3).map((row, i) => h("div", {class: "ga-top-item"},
      h("div", {class: "best-label", text: i === 0 ? "highest" : `#${i + 1}`}),
      h("div", {class: "big", text: fix(row.score)}),
      h("div", {class: "mono", text: (row.key || "").slice(0, 12)})))))));
}

function renderStepC(E) {
  const head = $("h-stepc");
  if (head) head.textContent = (E && ((E.version || "").startsWith("v9.34") || (E.version || "").startsWith("v9.33") || (E.version || "").startsWith("v9.32") || (E.version || "").startsWith("v9.31") || (E.version || "").startsWith("v9.30") || (E.version || "").startsWith("v9.29") || (E.version || "").startsWith("v9.28") || (E.version || "").startsWith("v9.27") || (E.version || "").startsWith("v9.26") || (E.version || "").startsWith("v9.25") || (E.version || "").startsWith("v9.24") || (E.version || "").startsWith("v9.23")))
    ? "Step C: the fine-tune"
    : "Step C: two generations, then the fine-tune";
  const box = $("stepc");
  const log = (E && E.step_c) || [];
  const waves = ((E && E.waves) || []).filter((w) => /^step C/i.test(w.label || ""));
  const fine = climbers(E).filter((c) => c.origin === "step C");
  const phase = (E && E.phase) || "";
  if (!E || (!log.length && !waves.length && !fine.length && !phase.startsWith("C"))) {
    const at = E && E.schedule ? mins(E.schedule.b_end_s) : "the fine-tune";
    const keep = E && E.schedule && E.schedule.c_reserve_s != null ? mins(E.schedule.c_reserve_s) : "8 min";
    const v929 = E && ((E.version || "").startsWith("v9.34") || (E.version || "").startsWith("v9.33") || (E.version || "").startsWith("v9.32") || (E.version || "").startsWith("v9.31") || (E.version || "").startsWith("v9.30") || (E.version || "").startsWith("v9.29"));
    const v923 = E && ((E.version || "").startsWith("v9.28") || (E.version || "").startsWith("v9.27") || (E.version || "").startsWith("v9.26") || (E.version || "").startsWith("v9.25") || (E.version || "").startsWith("v9.24") || (E.version || "").startsWith("v9.23"));
    box.replaceChildren(h("div", {class: "empty", text: v929
      ? `Starts at ${at}. Every climber except the highest leaves. The fine-tune keeps ${keep}. A Step B score that returns after this opens replaces the maximum when it is higher.`
      : v923
      ? `Starts at ${at}. Climbers keep their seats until then. The fine-tune keeps ${keep} and does not breed.`
      : `Starts at ${at}, after every climber has left Step B. Two generations, then the fine-tune keeps ${keep}.`}));
    return;
  }
  const byId = new Map(climbers(E).map((c) => [c.id, c]));
  const OUT = {running: "running", moved: "found a better point: Step B climbs on", "second highs": "second high points: short climbs",
    "landing won": "a short climb won: Step B climbs on", final: "the peak is final"};
  const rows = [];
  if (phase.startsWith("C")) rows.push(h("div", {class: "muted", text: phase}));
  for (const w of waves) {
    rows.push(h("div", {class: "rowline"},
      h("div", {class: "top-line"}, h("b", {text: w.label}),
        h("span", {class: "num muted", text: `${int(w.back)} / ${int(w.sent)} back · ${mins(w.started_s)}`})),
      bar(w.sent ? w.back / w.sent : 0, "merge", null)));
  }
  for (const c of fine) {
    rows.push(h("div", {class: "rowline"},
      h("div", {class: "top-line"}, h("b", {text: `Fine-tune C${c.id}`}),
        h("span", {class: "num muted", text: `${fix(c.opened)} → ${fix(c.best)}`})),
      h("div", {class: "muted", text: `${STATUS_TEXT[c.status] || c.status} · ${int(c.iterations)} grids`})));
  }
  for (const e of log) {
    const shorts = (e.landings || []).map((id) => byId.get(id)).filter(Boolean);
    rows.push(h("div", {class: "rowline"},
      h("div", {class: "top-line"}, h("b", {text: `C${e.peak} · pass ${e.pass}`}),
        h("span", {class: "num muted", text: `${int(e.points)} points · ${int(e.gaps)} gaps · ${mins(e.started_s)}`})),
      h("div", {}, chip(OUT[e.outcome] || e.outcome, e.outcome === "moved" || e.outcome === "landing won" ? "good" : e.outcome === "running" ? "climbing" : "")),
      h("div", {class: "muted", text: `from ${fix(e.from)}${e.detail && !(shorts.length && e.detail.startsWith("short climbs")) ? ` · ${e.detail}` : ""}`}),
      shorts.length ? h("div", {class: "muted", text: `short climbs: ${shorts.map((d) => `C${d.id} ${fix(d.opened)} → ${fix(d.best)}`).join(", ")}`}) : null));
  }
  box.replaceChildren(...rows);
}

const GA_SHORT = {base_quality_score_threshold: "BQST", active_probability_threshold: "APT",
  min_base_quality_score: "MBQ", phred_scaled_global_read_mismapping_rate: "mismapping", pcr_indel_model: "PCR",
  min_pruning: "min_pruning"};

function gaStageTitle(name) {
  if (name === "A1") return "A1 · every value of the 6 influential settings";
  if (name === "A3") return "A3 · island heads (and fill)";
  const m = /^A4_(\d+)_pm(.+)$/.exec(name);
  return m ? `Re-tune pass ${m[1]} (±${m[2].replace("+", ", ±")} rungs)` : name;
}

function renderGlobalA(G) {
  const panel = $("globala-panel"), box = $("globala");
  if (!G) { panel.hidden = true; return; }
  panel.hidden = false;
  $("globala-sub").textContent = globalASub(G);
  box.replaceChildren(...globalARows(G));
}

function globalASub(G) {
  const done = G.done_minute != null;
  return `${done ? `done at t+${fix(G.done_minute, 1)}m` : "running"} · stops by t+${fix(G.stop_minute, 1)}m` +
    ` · fill ${G.fill ? "on" : "off"}`;
}

function globalARows(G) {
  const rows = Object.entries(G.stages || {}).map(([name, st]) => h("div", {class: "rowline"},
    h("div", {class: "top-line"}, h("b", {text: gaStageTitle(name)}),
      h("span", {class: "num muted", text: `${int(st.scored)} / ${int(st.sent)} scored · t+${fix(st.minute, 1)}m`})),
    bar(st.sent ? st.scored / st.sent : 0, "open", null),
    h("div", {class: "muted", text: `best ${st.best != null ? fix(st.best) : "–"}${st.cut ? " · cut at the stop time" : ""}`})));
  const hills = Object.entries(G.hills || {});
  if (hills.length) rows.push(h("div", {class: "rowline"}, h("b", {text: "Hills (one island per hill and PCR model)"}),
    kv(hills.map(([n, hs]) => [GA_SHORT[n] || n, hs.map((x) => `${x[0]} (${x[1]}–${x[2]})`).join(" · ")]))));
  const climb = G.climbing || [];
  if (climb.length) rows.push(h("div", {class: "rowline"}, h("b", {text: `Islands climbing (${climb.length})`}),
    kv(climb.map((c) => [`Island ${c.island} · ${c.influential.pcr_indel_model} / mp ${c.influential.min_pruning}`,
      `${fix(c.score)} · band ${Object.entries(c.band || {}).map(([n, b]) => `${GA_SHORT[n] || n} ${b[0]}–${b[1]}`).join(", ") || "–"}`]))));
  const segs = G.segments || [];
  if (segs.length) rows.push(h("div", {class: "rowline"}, h("b", {text: `Segments found (${segs.length})`}),
    kv(segs.slice(0, 8).map((sg, i) => [`#${i + 1} · ${sg.influential.pcr_indel_model} / MBQ ${sg.influential.min_base_quality_score}`,
      `${fix(sg.score)} · BQST ${sg.influential.base_quality_score_threshold} · APT ${sg.influential.active_probability_threshold}`]))));
  return rows.length ? rows : [h("div", {class: "empty", text: "Step A is starting."})];
}

function renderStepA(E, R) {
  const box = $("waves");
  const G = last && last.search && last.search.global_a;
  if (G) {                         // island mode: its global Step A (the engine's own waves stay empty)
    $("h-stepa").textContent = `Step A: global (island mode) · ${globalASub(G)}`;
    box.replaceChildren(...globalARows(G));
    return;
  }
  $("h-stepa").textContent = "Step A: candidates";
  if (!E || !(E.waves || []).length) { box.replaceChildren(h("div", {class: "empty", text: R && R.mode ? R.mode : "Not started."})); return; }
  const stepA = E.waves.filter((w) => !/^step C/i.test(w.label || ""));
  const waves = stepA.map((w, i) => {
    const next = stepA[i + 1];
    const gain = next && w.best_at_launch != null && next.best_at_launch != null ? next.best_at_launch - w.best_at_launch : null;
    const title = w.wave === 0 ? "Wave 0: seeds, portfolio, hypercube" : (w.label || `Generation ${w.wave}`);
    return h("div", {class: "rowline"},
      h("div", {class: "top-line"}, h("b", {text: title}),
        h("span", {class: "num muted", text: `${int(w.back)} / ${int(w.sent)} back`})),
      bar(w.sent ? w.back / w.sent : 0, "open", null),
      gain != null ? h("div", {class: "muted", text: `the next wave started ${gain >= 0 ? "+" : ""}${fix(gain)} higher`}) : null);
  });
  const cands = E.candidates || [];
  const per = (last && last.settings && last.settings.values && last.settings.values["candidates.workers_per_candidate"]) || null;
  box.replaceChildren(...waves,
    cands.length ? h("div", {class: "rowline"}, h("b", {text: `Candidates (${cands.length})`}),
      kv(cands.map((c, i) => [`C${i} · ${c.group}`, `${fix(c.score)} · core ${fix(c.core)} · FP ${int(c.fp)}`])))
      : h("div", {class: "muted", text: (E.version || "").startsWith("v9.34") || (E.version || "").startsWith("v9.33") || (E.version || "").startsWith("v9.32") || (E.version || "").startsWith("v9.31") || (E.version || "").startsWith("v9.30")
        ? "Wave 0, one dummy on each apparent peak, then two scout waves from the far parents."
        : (E.version || "").startsWith("v9.29")
        ? "Wave 0, then one dummy on each apparent peak. One climber sits on the highest fingerprint."
        : (E.version || "").startsWith("v9.6") || (E.version || "").startsWith("v9.5") || (E.version || "").startsWith("v9.4") || (E.version || "").startsWith("v9.3")
        ? "Seats stay within 0.10 of the best score, picked when generation 1 is back."
        : per
          ? `One candidate per ${int(per)} workers, picked when generation 1 is back.`
          : "The candidates are picked when generation 1 is back."}));
}

function renderPick(E, finished) {
  const box = $("pick");
  if (!E) { box.replaceChildren(h("div", {class: "empty", text: "Waiting for the engine."})); return; }
  const cl = climbers(E);
  const best = cl.length ? cl.reduce((a, b) => (b.best > a.best ? b : a)) : null;
  box.replaceChildren(
    h("div", {class: "rowline"}, h("b", {text: "Named pick"}),
      E.pick ? kv([["config", h("span", {class: "mono", text: E.pick.key.slice(0, 12)})], ["its score", fix(E.pick.score)],
        ["named at", mins(E.pick.at_s)]])
        : h("div", {class: "muted", text: `the highest local score, named at ${E.schedule ? mins(E.schedule.pick_s) : "the end"}`})),
    h("div", {class: "rowline"}, h("b", {text: "Best climber"}),
      best ? kv([["climber", `C${best.id}`], ["score", fix(best.best)], ["from", best.origin],
        ["gain over its start", `${best.best - best.opened >= 0 ? "+" : ""}${fix(best.best - best.opened)}`]])
        : h("div", {class: "muted", text: "none yet"})),
    kv([["best of everything scored", fix(E.best)], ["workers", int(E.capacity)]]),
  );
}

// ---- chart ----------------------------------------------------------------------------------------------------------
// The live best (and scores near it) fill most of the height, still linear inside that band, so a gain of a few
// tenths is a visible step. Scores below the band stay on the chart, in a short strip at the bottom.
let FOCUS_LO = 70, FOCUS_HI = 100;
const FOCUS_SHARE = 0.84;
const CHART_H = 440;
function scoreWindow(S, subs) {
  const vals = [];
  for (const p of S.points || []) if (p[1] > 0) vals.push(p[1]);
  for (const p of S.best_timeline || []) if (p[1] > 0) vals.push(p[1]);
  for (const p of subs || []) if (p[2] > 0) vals.push(p[2]);
  if (S.best > 0) vals.push(S.best);
  if (!vals.length) return {lo: 70, hi: 100};
  const best = Math.max(...vals);
  const near = vals.filter((v) => v >= best - 8);
  const lo = Math.min(...near);
  let floor = Math.floor(lo * 2) / 2;
  let top = Math.ceil((best + 0.4) * 2) / 2;
  if (top - floor < 2) floor = Math.round((top - 2) * 2) / 2;
  if (floor < 0) floor = 0;
  return {lo: floor, hi: top};
}
function toFocus(v) {
  if (v == null || !isFinite(v)) return null;
  const low = 1 - FOCUS_SHARE;
  if (v <= 0) return 0;
  if (v < FOCUS_LO) return low * (v / Math.max(FOCUS_LO, 1e-9));
  return low + FOCUS_SHARE * ((v - FOCUS_LO) / Math.max(FOCUS_HI - FOCUS_LO, 1e-9));
}
function fromFocus(a) {
  const low = 1 - FOCUS_SHARE;
  if (a <= low) return FOCUS_LO * (a / Math.max(low, 1e-9));
  return FOCUS_LO + (FOCUS_HI - FOCUS_LO) * ((a - low) / FOCUS_SHARE);
}
function focusTicks() {
  const span = FOCUS_HI - FOCUS_LO;
  const step = span > 10 ? 5 : span > 4 ? 1 : 0.5;
  const ticks = [0];
  if (FOCUS_LO > 1) ticks.push(FOCUS_LO / 2);
  for (let x = FOCUS_LO; x <= FOCUS_HI + 1e-9; x += step) ticks.push(Math.round(x * 100) / 100);
  if (ticks[ticks.length - 1] < FOCUS_HI - 1e-9) ticks.push(FOCUS_HI);
  return ticks;
}
function chartData(S, subs) {
  const pts = (S.points || []).filter((p) => p[1] > 0);
  const best = S.best_timeline || [];
  const map = new Map();
  const put = (t, i, v) => { const k = Math.round(t * 10) / 10; if (!map.has(k)) map.set(k, [null, null, null]); map.get(k)[i] = v; };
  for (const [t, v] of pts) put(t, 0, v);
  for (const [t, v] of best) put(t, 1, v);
  for (const [t, , v] of subs) put(t, 2, v);
  const xs = [...map.keys()].sort((a, b) => a - b);
  const col = (i) => xs.map((t) => toFocus(map.get(t)[i]));
  return [xs, col(0), col(1), col(2)];
}
function sizeChart() {
  const el = $("chart");
  if (chart && el.clientWidth) chart.setSize({width: el.clientWidth, height: CHART_H});
  const coresEl = $("cores-chart");
  if (coresChart && coresEl && coresEl.clientWidth) coresChart.setSize({width: coresEl.clientWidth, height: 220});
}

function noteLiveCores(L) {
  if (!L) return;
  const rid = (L.round || {}).run_id;
  if (rid !== liveCoreRun) { liveCoreTrace = []; liveCoreRun = rid; }
  const X = L.executor || {};
  if (X.running == null || X.workers == null) return;
  const t = L.updated || now();
  const lastPt = liveCoreTrace[liveCoreTrace.length - 1];
  if (lastPt && Math.abs(lastPt[0] - t) < 1 && lastPt[1] === X.running && lastPt[2] === X.workers) return;
  liveCoreTrace.push([t, X.running, X.workers]);
}

function coreTraceFromSearch(S) {
  return (S && S.core_timeline) || [];
}

function mergeCoreTraces(saved, live) {
  const map = new Map();
  for (const row of saved || []) {
    if (!row || row[0] == null) continue;
    map.set(Math.round(row[0] * 10) / 10, [row[1], row[2]]);
  }
  for (const row of live || []) {
    if (!row || row[0] == null) continue;
    map.set(Math.round(row[0] * 10) / 10, [row[1], row[2]]);
  }
  return [...map.keys()].sort((a, b) => a - b).map((t) => [t, map.get(t)[0], map.get(t)[1]]);
}

function renderCoresChart(S, X, R) {
  const el = $("cores-chart");
  if (!el || typeof uPlot === "undefined") return;
  const saved = coreTraceFromSearch(S);
  const rows = mergeCoreTraces(saved, liveCoreTrace);
  if (!rows.length) {
    if (coresChart) { coresChart.destroy(); coresChart = null; }
    el.replaceChildren(h("div", {class: "empty", text: "Busy cores appear as the fleet reports in."}));
    return;
  }
  const xs = rows.map((r) => r[0]);
  const busy = rows.map((r) => r[1]);
  const cap = rows.map((r) => r[2]);
  const theme = effectiveTheme();
  if (!coresChart || coresChartKey !== theme) {
    if (coresChart) coresChart.destroy();
    el.replaceChildren();
    coresChartKey = theme;
    const axis = {stroke: css("--muted"), grid: {stroke: css("--line"), width: 1}, ticks: {stroke: css("--line")}};
    coresChart = new uPlot({
      width: el.clientWidth || 800, height: 220,
      scales: {x: {time: true}, y: {range: (u, lo, hi) => [0, Math.max(hi, 1)]}},
      axes: [{...axis}, {...axis, size: 48, values: (u, v) => v.map((x) => String(Math.round(x)))}],
      series: [
        {},
        {label: "busy", stroke: css("--good"), width: 2, spanGaps: true, points: {show: false},
         paths: uPlot.paths.stepped({align: 1})},
        {label: "fleet", stroke: css("--faint"), width: 1, spanGaps: true, points: {show: false},
         paths: uPlot.paths.stepped({align: 1})},
      ],
    }, [xs, busy, cap], el);
  } else {
    coresChart.setData([xs, busy, cap]);
  }
}
function renderChart(S, subs) {
  const el = $("chart");
  if (typeof uPlot === "undefined") return;
  const win = scoreWindow(S, subs);
  FOCUS_LO = win.lo;
  FOCUS_HI = win.hi;
  const sub = $("chart-sub");
  if (sub) {
    sub.textContent = S.best != null
      ? `best ${fix(S.best)} · ${FOCUS_LO.toFixed(1)}–${FOCUS_HI.toFixed(1)} fills the chart; scores below ${FOCUS_LO.toFixed(1)} stay in the bottom band`
      : `the live best fills this range; scores below it stay in the bottom band`;
  }
  const data = chartData(S, subs);
  if (!data[0].length) { if (chart) { chart.destroy(); chart = null; } el.replaceChildren(h("div", {class: "empty", text: "No scores yet."})); return; }
  const theme = effectiveTheme();
  const key = `${theme}:${FOCUS_LO}:${FOCUS_HI}`;
  if (!chart || chartKey !== key) {
    if (chart) chart.destroy();
    el.replaceChildren();
    chartKey = key;
    const axis = {stroke: css("--muted"), grid: {stroke: css("--line"), width: 1}, ticks: {stroke: css("--line")}};
    const scoreOf = (u, v) => (v == null ? "" : fromFocus(v).toFixed(4));
    const ticks = focusTicks();
    const label = (x) => {
      const s = fromFocus(x);
      return FOCUS_HI - FOCUS_LO <= 4 ? s.toFixed(1) : String(Math.round(s));
    };
    chart = new uPlot({
      width: el.clientWidth || 800, height: CHART_H,
      scales: {x: {time: true, range: (u, lo, hi) => (hi - lo < 120 ? [lo - 60, hi + 60] : [lo, hi])},
               y: {range: (u, lo, hi) => [0, Math.max(1, hi)]}},
      axes: [
        {...axis},
        {...axis, size: 52, splits: () => ticks.map(toFocus),
         values: (u, v) => v.map(label)},
      ],
      series: [
        {},
        {label: "scored config", stroke: css("--faint"), paths: () => null, points: {show: true, size: 3, fill: css("--faint"), stroke: css("--faint")}, value: scoreOf},
        {label: "best so far", stroke: css("--accent"), width: 2, spanGaps: true, points: {show: false}, paths: uPlot.paths.stepped({align: 1}), value: scoreOf},
        {label: "submitted", stroke: css("--final"), paths: () => null, points: {show: true, size: 9, fill: css("--final"), stroke: css("--final")}, value: scoreOf},
      ],
      hooks: {draw: [(u) => {
        const y70 = u.valToPos(toFocus(FOCUS_LO), "y", true);
        const yTop = u.valToPos(u.scales.y.max, "y", true);
        const ctx = u.ctx;
        ctx.save();
        ctx.fillStyle = css("--accent");
        ctx.globalAlpha = 0.07;
        ctx.fillRect(u.bbox.left, Math.min(y70, yTop), u.bbox.width, Math.abs(y70 - yTop));
        ctx.globalAlpha = 1;
        ctx.strokeStyle = css("--warn");
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(u.bbox.left, y70);
        ctx.lineTo(u.bbox.left + u.bbox.width, y70);
        ctx.stroke();
        ctx.restore();
      }]},
    }, data, el);
  } else {
    chart.setData(data);
  }
}
function sizeHistCharts() {
  for (const plot of Object.values(histCharts)) {
    const el = plot && plot.root && plot.root.parentElement;
    if (plot && el && el.clientWidth) plot.setSize({width: el.clientWidth, height: 200});
  }
}
window.addEventListener("resize", () => { sizeChart(); sizeHistCharts(); if (last) renderLineage((last.search || {}).engine); });

function renderNodes(X, finished) {
  const box = $("nodes");
  if (!X) { box.replaceChildren(h("div", {class: "empty", text: "No evaluator yet."})); $("nodes-sub").textContent = ""; return; }
  const extra = X.kind === "fleet" ? ` · ${int(X.speculated || 0)} slow runs duplicated, ${int(X.speculation_won || 0)} won` : "";
  $("nodes-sub").textContent = finished ? `${X.kind} · round over: no work on any node${extra}`
    : `${X.kind} · ${int(X.running)} of ${int(X.workers)} workers busy · ${int(X.queued)} queued${extra}`;
  const nodes = X.nodes && X.nodes.length ? X.nodes : [{node: "this machine", ready: true, workers: X.workers, busy: X.running, small: null, solo_s: X.solo_s}];
  box.replaceChildren(...nodes.map((n) => {
    const busy = finished ? 0 : n.busy;                    // a finished round has no work anywhere
    const state = n.rejected ? chip("rejected", "bad") : finished ? chip("idle") : !n.ready ? chip("not ready", "warn")
      : n.sick ? chip("very slow", "bad") : n.slow ? chip("slow", "warn") : chip("ready", "good");
    return h("div", {class: `node ${n.ready && !finished ? "" : "down"}`},
      h("div", {class: "head"}, h("span", {text: n.node}), state),
      bar(n.workers ? busy / n.workers : 0, "", null),
      h("div", {class: "num muted", text: `${int(busy)} / ${int(n.workers)} busy${!finished && n.small != null ? ` · ${int(n.small)} small` : ""}${n.solo_s ? ` · ${int(n.solo_s)} s per run` : ""}`}));
  }));
}

function eventTag(row) {
  if (row.kind === "setting") return "setting";
  if (row.kind === "error") return "error";
  if (row.kind === "warning") return "warning";
  const text = row.text || "";
  if (/^step B\b/i.test(text) || /^climber\b/i.test(text)) return "step B";
  if (/^GA\b/.test(text)) return "GA";
  if (/^step C\b/i.test(text) || /^GA cross\b/.test(text)) return "step C";
  if (/^step A\b/i.test(text) || /^wave\b/i.test(text) || /^generation\b/i.test(text)
      || /^valley\b/i.test(text) || /^no valley\b/i.test(text) || /^met\b/i.test(text)) return "step A";
  return "engine";
}

let feedOpen = null;

function renderFeed(L, E) {
  const rows = [];
  for (const [t, text] of (E && E.events) || []) rows.push({t, kind: "engine", text});
  for (const c of ((L.settings || {}).changes) || []) rows.push({t: c.t, kind: "setting", text: `${c.name}: ${JSON.stringify(c.old)} → ${JSON.stringify(c.new)}`});
  for (const [t, level, name, msg] of L.log || []) rows.push({t, kind: level === "ERROR" ? "error" : "warning", text: `${name}: ${msg}`});
  rows.sort((a, b) => b.t - a.t);
  const order = ["step B", "GA", "step A", "step C", "engine", "setting", "warning", "error"];
  const tone = { "step B": "climb", "GA": "scout", "step A": "opening", "step C": "merge", setting: "scout", warning: "warn", error: "bad", engine: "" };
  const groups = new Map(order.map((tag) => [tag, []]));
  for (const r of rows) {
    const tag = eventTag(r);
    if (!groups.has(tag)) groups.set(tag, []);
    groups.get(tag).push(r);
  }
  const box = $("feed");
  const tags = [...groups.keys()].filter((tag) => groups.get(tag).length);
  if (!tags.length) { box.replaceChildren(h("div", {class: "empty", text: "Nothing yet."})); return; }
  box.replaceChildren(...tags.flatMap((tag) => {
    const items = groups.get(tag);
    const open = feedOpen === tag;
    const head = h("article", {class: open ? "hist-run sel" : "hist-run",
      onclick: () => { feedOpen = open ? null : tag; if (last) renderFeed(last, (last.search || {}).engine); }},
      h("div", {class: "who"}, tag, h("span", {class: "sub", text: `${int(items.length)} events`})),
      h("span", {class: "spacer"}),
      chip(tag, tone[tag] || ""));
    if (!open) return head;
    const shown = items.slice(0, 80);
    return [head, h("div", {class: "hist-detail"},
      ...shown.map((r) => h("div", {class: "feed-line"},
        h("span", {class: "t", text: clock(r.t), title: E ? `${mins(r.t - E.started_ts)} into the search` : ""}),
        h("span", {class: "msg", text: r.text}))),
      items.length > shown.length ? h("div", {class: "muted", text: `${int(shown.length)} of ${int(items.length)}, newest first`}) : null)];
  }));
}

// ---- settings -------------------------------------------------------------------------------------------------------
async function loadSettings() {
  try {
    const r = await api("/api/settings");
    settingsMeta = await r.json();
    renderSettings();
  } catch (e) { /* the poll shows the server state */ }
}

function inForce() { return ((last && last.settings) || {}).values || null; }
function shown(v) { return Array.isArray(v) ? v.join(", ") : v == null ? "" : String(v); }
function parse(knob, text) {
  const t = text.trim();
  if (knob.kind === "list") return t.split(/[\s,]+/).filter(Boolean).map(Number);
  return Number(t);
}

function renderSettings() {
  const M = settingsMeta;
  if (!M) return;
  const groups = $("settings-groups");
  const byGroup = {};
  for (const k of M.knobs) (byGroup[k.group] = byGroup[k.group] || []).push(k);
  const rejected = Object.fromEntries((M.rejected || []).map((r) => [r.name, r.reason]));
  groups.replaceChildren(...Object.keys(M.groups).filter((g) => byGroup[g]).map((g) => h("section", {class: "panel group"},
    h("h3", {text: M.groups[g]}),
    h("div", {class: "stable-wrap"}, h("table", {class: "stable"},
      h("colgroup", {}, ["c-name", "c-force", "c-default", "c-value", "c-act"].map((c) => h("col", {class: c}))),
      h("thead", {}, h("tr", {}, ["Setting", "In force on main", "Default", "Your value", ""].map((c) => h("th", {text: c})))),
      h("tbody", {}, byGroup[g].map((k) => settingRow(k, M, rejected[k.name]))))))));
  const unknown = Object.keys(M.overrides || {}).filter((n) => !M.defaults.hasOwnProperty(n));
  const banner = $("settings-banner");
  banner.replaceChildren();
  if (M.error) banner.append(h("div", {class: "banner bad", text: `${M.error}. Fix or delete ${M.file} on main; until then nothing can be saved here.`}));
  if (unknown.length) banner.append(h("div", {class: "banner", text: `The file has settings main does not know: ${unknown.join(", ")}. Saving removes them only if you reset them.`},
    " ", h("button", {class: "btn small", type: "button", text: "Reset them", onclick: () => { unknown.forEach((n) => resets.add(n)); updateButtons(); }})));
  renderSettingsState();
  updateButtons();
}

function settingRow(k, M, reason) {
  const over = M.overrides.hasOwnProperty(k.name);
  const value = draft.hasOwnProperty(k.name) ? draft[k.name] : over ? shown(M.overrides[k.name]) : "";
  const input = h("input", {id: `in-${k.name}`, value, placeholder: shown(k.default), inputmode: k.kind === "list" ? "text" : "decimal",
    "aria-label": k.name, class: draft.hasOwnProperty(k.name) ? "changed" : ""});
  input.addEventListener("input", () => {
    const base = over ? shown(M.overrides[k.name]) : "";
    if (input.value === base) delete draft[k.name]; else draft[k.name] = input.value;
    resets.delete(k.name);
    input.classList.toggle("changed", draft.hasOwnProperty(k.name));
    input.classList.remove("invalid");
    updateButtons();
  });
  const unit = k.unit && k.unit !== "count" ? ` ${k.unit}` : "";
  return h("tr", {class: over ? "overridden" : ""},
    h("td", {class: "name"}, h("b", {class: "mono", text: k.name}), " ", k.when === "round" ? chip("next round", "warn") : chip("live", "live"),
      h("div", {class: "help", text: `${k.help} Allowed ${k.lo} to ${k.hi}${unit}.`})),
    h("td", {class: "val", id: `force-${k.name}`}),
    h("td", {class: "val muted", text: `${shown(k.default)}${unit}`}),
    h("td", {}, input, reason ? h("div", {class: "reason", text: `main rejects this: ${reason}`}) : null, h("div", {class: "reason", id: `err-${k.name}`})),
    h("td", {}, over ? h("button", {class: "btn small", type: "button", text: "Use default", onclick: () => {
      resets.add(k.name); delete draft[k.name]; input.value = ""; input.classList.add("changed"); updateButtons();
    }}) : null));
}

function renderSettingsState() {
  const st = (last && last.settings) || null;
  const el = $("settings-state");
  if (!st) { el.textContent = "Main's settings are unknown until it starts a round."; }
  else el.textContent = `${st.loaded_at ? `main read the file ${ago(st.loaded_at)}` : "main has not read the file yet"}${Object.keys(st.pending || {}).length ? ` · ${Object.keys(st.pending).length} waiting for the next round` : ""}`;
  if (!settingsMeta) return;
  const vals = inForce() || {}, pending = (st && st.pending) || {};
  for (const k of settingsMeta.knobs) {
    const cell = document.getElementById(`force-${k.name}`);
    if (!cell) continue;
    const v = vals.hasOwnProperty(k.name) ? shown(vals[k.name]) : "?";
    cell.replaceChildren(...[h("span", {text: v}), pending.hasOwnProperty(k.name) ? h("div", {class: "muted", text: `→ ${shown(pending[k.name])} next round`}) : null].filter(Boolean));
  }
}

function updateButtons() {
  const n = Object.keys(draft).length + resets.size;
  $("save").disabled = n === 0 || !!(settingsMeta && settingsMeta.error);
  $("discard").disabled = n === 0;
  $("save").textContent = n ? `Save ${n} change${n > 1 ? "s" : ""}` : "Save changes";
}

$("discard").addEventListener("click", () => {
  for (const k of Object.keys(draft)) delete draft[k];
  resets.clear();
  renderSettings();
});

$("save").addEventListener("click", async () => {
  const M = settingsMeta;
  const set = {};
  let bad = false;
  for (const [name, text] of Object.entries(draft)) {
    const k = M.knobs.find((x) => x.name === name);
    const err = document.getElementById(`err-${name}`);
    if (text.trim() === "") { resets.add(name); continue; }
    const v = parse(k, text);
    if ((Array.isArray(v) && (!v.length || v.some((x) => !isFinite(x)))) || (!Array.isArray(v) && !isFinite(v))) {
      bad = true;
      if (err) err.textContent = k.kind === "list" ? "write whole numbers separated by commas" : "write a number";
      document.getElementById(`in-${name}`).classList.add("invalid");
      continue;
    }
    set[name] = v;
  }
  if (bad) return;
  $("save").disabled = true;
  try {
    const r = await api("/api/settings", {method: "PUT", body: JSON.stringify({set, reset: [...resets]})});
    const body = await r.json();
    if (!r.ok) {
      for (const x of body.rejected || []) {
        const err = document.getElementById(`err-${x.name}`);
        if (err) err.textContent = x.reason;
        const inp = document.getElementById(`in-${x.name}`);
        if (inp) inp.classList.add("invalid");
      }
      toast(body.error || `Not saved (${r.status})`);
      updateButtons();
      return;
    }
    for (const k of Object.keys(draft)) delete draft[k];
    resets.clear();
    settingsMeta = body;
    renderSettings();
    toast("Saved. Main applies live settings within about a second; round settings at the next round.");
  } catch (e) {
    toast("Not saved: the server did not answer.");
    updateButtons();
  }
});

// ---- history: finished runs from rounds.jsonl. Loaded when the tab opens, not on the live poll. ----
let historyRuns = [];
let historySelected = null;
const histCharts = {};

function when(ts) {
  if (ts == null) return "";
  const d = new Date(ts * 1000);
  const p = (n) => String(n).padStart(2, "0");
  return `${d.getUTCFullYear()}-${p(d.getUTCMonth() + 1)}-${p(d.getUTCDate())} ${p(d.getUTCHours())}:${p(d.getUTCMinutes())} UTC`;
}
function rankOf(r) {
  if (r.rank == null || r.field == null) return "–";
  return `${int(r.rank)} / ${int(r.field)}`;
}
function searchedFor(r) {
  if (r.started == null || r.ended == null) return "–";
  const min = Math.max(0, Math.round((r.ended - r.started) / 60));
  const h = Math.floor(min / 60), m = min % 60;
  return h ? `${h} h ${m} min` : `${m} min`;
}
function reviewState(r) {
  if (r.ranked) return ["ranked", "good"];
  if (r.reviewed) return ["reviewed", "warn"];
  return ["waiting", ""];
}
function stat(tone, label, value) {
  return h("div", {class: `hist-stat ${tone}`}, h("div", {class: "label", text: label}), h("div", {class: "value", text: value}));
}
function group(title, tiles) {
  return h("div", {class: "hist-group"}, h("div", {class: "cap", text: title}), h("div", {class: "hist-stats"}, ...tiles));
}

async function loadHistory() {
  $("history-state").textContent = "Loading…";
  try {
    const r = await api("/api/history");
    if (!r.ok) { $("history-state").textContent = "Could not load history."; return; }
    const body = await r.json();
    historyRuns = body.runs || [];
    if (historySelected && !historyRuns.some((x) => x.run_id === historySelected)) historySelected = null;
    $("history-state").textContent = historyRuns.length ? `${historyRuns.length} finished runs` : "No finished runs yet.";
    renderHistoryCharts();
    renderHistory();
  } catch (e) { /* a 401 already sent the browser to the login page */ }
}

function historyStamp(tag) {
  return tag && tag.length >= 15 ? `${tag.slice(6, 8)} ${tag.slice(9, 11)}:${tag.slice(11, 13)}` : "";
}

function mountHist(id, runs, series, kind, empty) {
  const el = $(id);
  if (typeof uPlot === "undefined") return null;
  const cols = series.map((s) => runs.map(s.get));
  if (!cols.some((c) => c.some((v) => v != null))) {
    el.replaceChildren(h("div", {class: "empty", text: empty || "Nothing reviewed yet."}));
    return null;
  }
  el.replaceChildren();
  const axis = {stroke: css("--muted"), grid: {stroke: css("--line"), width: 1}, ticks: {stroke: css("--line")}};
  const yText = kind === "rho" ? (x) => x.toFixed(3) : kind === "n" ? (x) => String(Math.round(x)) : (x) => x.toFixed(2);
  return new uPlot({
    width: el.clientWidth || 320, height: 200,
    scales: {
      x: {time: false, range: (u, lo, hi) => hi - lo < 1 ? [lo - 0.6, hi + 0.6] : [lo - 0.4, hi + 0.4]},
      y: {range: (u, lo, hi) => { const span = hi - lo || Math.abs(hi) || 1; return [lo - span * 0.08, hi + span * 0.08]; }},
    },
    axes: [
      {...axis, values: (u, splits) => splits.map((v) => {
        const i = Math.round(v);
        if (Math.abs(v - i) > 0.05 || !runs[i]) return "";
        return historyStamp(runs[i].tag);
      })},
      {...axis, size: 48, values: (u, v) => v.map(yText)},
    ],
    series: [{label: "round"}, ...series.map((s) => ({label: s.label, stroke: css(s.color), width: 2, spanGaps: false,
      points: {show: true, size: 5, fill: css(s.color), stroke: css(s.color)}}))],
  }, [runs.map((_, i) => i), ...cols], el);
}

const CORE_COLORS = ["--accent", "--good", "--final", "--scout", "--open", "--merge", "--warn"];

function renderHistoryCoreChart() {
  if (histCharts.cores) { histCharts.cores.destroy(); histCharts.cores = null; }
  const el = $("hist-cores");
  if (!el || typeof uPlot === "undefined") return null;
  let shown = historyRuns.filter((r) => (r.core_timeline || []).length);
  if (shown.length > 10) {
    const keep = shown.slice(0, 9);
    const extra = historySelected && shown.find((r) => r.run_id === historySelected);
    if (extra && !keep.some((r) => r.run_id === extra.run_id)) keep.push(extra);
    shown = keep;
  }
  if (!shown.length) {
    el.replaceChildren(h("div", {class: "empty", text: "Busy-core traces start with the next finished round."}));
    return null;
  }
  const xsSet = new Set();
  for (const r of shown) for (const row of r.core_timeline) xsSet.add(Math.round(row[0] * 10) / 10);
  const xs = [...xsSet].sort((a, b) => a - b);
  const cols = [];
  const series = [];
  shown.forEach((r, i) => {
    const map = new Map(r.core_timeline.map((row) => [Math.round(row[0] * 10) / 10, row[1]]));
    const tmin = Math.min(...r.core_timeline.map((row) => row[0]));
    const tmax = Math.max(...r.core_timeline.map((row) => row[0]));
    let last = null;
    cols.push(xs.map((t) => {
      if (t < tmin - 1e-9 || t > tmax + 1e-9) { last = null; return null; }
      if (map.has(t)) last = map.get(t);
      return last;
    }));
    const on = r.run_id === historySelected;
    const color = on ? "--accent" : CORE_COLORS[i % CORE_COLORS.length];
    series.push({
      label: historyStamp(r.tag) || r.tag,
      stroke: css(color),
      width: on ? 3 : 1.5,
      spanGaps: false,
      points: {show: false},
      paths: uPlot.paths.stepped({align: 1}),
    });
  });
  el.replaceChildren();
  const axis = {stroke: css("--muted"), grid: {stroke: css("--line"), width: 1}, ticks: {stroke: css("--line")}};
  return new uPlot({
    width: el.clientWidth || 320, height: 220,
    scales: {
      x: {time: false, range: (u, lo, hi) => hi - lo < 1 ? [lo - 0.6, hi + 0.6] : [lo, hi]},
      y: {range: (u, lo, hi) => [0, Math.max(hi, 1)]},
    },
    axes: [
      {...axis, values: (u, v) => v.map((x) => `${Math.round(x)}m`)},
      {...axis, size: 48, values: (u, v) => v.map((x) => String(Math.round(x)))},
    ],
    series: [{label: "min"}, ...series],
  }, [xs, ...cols], el);
}

function renderHistoryCharts() {
  for (const k of Object.keys(histCharts)) {
    if (histCharts[k]) histCharts[k].destroy();
    histCharts[k] = null;
  }
  const runs = [...historyRuns].reverse();                          // oldest on the left
  const gapOf = (r) => {
    if (r.gap_to_best != null) return r.gap_to_best;
    const t = r.true_score != null ? r.true_score : r.final_true;
    return r.field_best != null && t != null ? r.field_best - t : null;
  };
  histCharts.score = mountHist("hist-score", runs, [
    {label: "gap", color: "--accent", get: gapOf},
  ], "score");
  histCharts.cores = renderHistoryCoreChart();
  histCharts.n = mountHist("hist-n", runs, [
    {label: "ours", color: "--accent", get: (r) => r.n_ours},
    {label: "truth", color: "--scout", get: (r) => r.n_truth},
  ], "n");
  histCharts.rho = mountHist("hist-rho", runs, [
    {label: "Spearman", color: "--final", get: (r) => r.spearman},
  ], "rho");
}

function renderHistory() {
  const list = $("history-list");
  if (!historyRuns.length) {
    list.replaceChildren(h("div", {class: "muted", text: "No finished runs yet."}));
  } else {
    list.replaceChildren(...historyRuns.flatMap((r) => {
      const sel = r.run_id === historySelected;
      const item = h("article", {class: sel ? "hist-run sel" : "hist-run",
        onclick: () => { historySelected = sel ? null : r.run_id; renderHistory(); histCharts.cores = renderHistoryCoreChart(); }},
        h("div", {class: "who"}, r.tag || r.run_id,
          h("span", {class: "sub", text: [r.region, r.kind, when(r.ended || r.started)].filter(Boolean).join(" · ")})),
        h("span", {class: "spacer"}),
        h("span", {class: `chip ${reviewState(r)[1]}`, text: reviewState(r)[0]}));
      if (!sel) return item;
      const trueScore = r.true_score != null ? r.true_score : r.final_true;
      const local = r.final_local != null ? r.final_local : r.best_local;
      const germ = r.germline_truth ? `${int(r.germline_match)} / ${int(r.germline_truth)}` : "–";
      const n = r.n_truth != null ? `${int(r.n_ours)} / ${int(r.n_truth)}` : int(r.n_ours);
      const who = [r.regime, r.donor].filter(Boolean).join(" ") || "–";
      return [item, h("div", {class: "hist-detail"},
        group("Score", [
          stat("score", "TRUE", fix(trueScore, 3)),
          stat("score", "Local at submit", fix(local, 3)),
          stat("score", "Rank", rankOf(r)),
          stat("score", "Field best", fix(r.field_best, 3)),
          stat("score", "Gap", fix(r.gap_to_best, 3)),
        ]),
        group("Search", [
          stat("work", "Cores", int(r.workers)),
          stat("work", "Evaluated", int(r.evaluated)),
          stat("work", "Searched", searchedFor(r)),
        ]),
        group("Answer set", [
          stat("set", "Donor", who),
          stat("set", "Score", r.answer_score != null ? fix(r.answer_score, 4) : "–"),
          stat("set", "Germline exact", germ),
          stat("set", "N", n),
          stat("set", "Build", r.build_s != null ? `${fix(r.build_s, 0)} s` : "–"),
        ]))];
    }));
  }
}

$("history-refresh").addEventListener("click", loadHistory);

// ---- start ----------------------------------------------------------------------------------------------------------
showTab();
poll();


// ---- curves: per-setting response curves over the stored rounds (GET /api/curves) ---------------------------------
const CURVE_GROUPS = [
  ["Stable peak", "The best value sits in nearly the same place in most rounds: start at the usual peak and search near it."],
  ["Depends on BAM type", "The peak follows the BAM's base qualities: low on binned BAMs (few quality levels), high on full-quality BAMs."],
  ["Moves per round", "The peak changes from round to round, so these need searching every round."],
  ["Flat", "The score barely changes across the whole range."],
];
const CURVE_YMIN = -3;
let curvesData = null, curvesLoading = false;

async function loadCurves() {
  if (curvesData) { renderCurves(); return; }
  if (curvesLoading) return;
  curvesLoading = true;
  try {
    const r = await api("/api/curves");
    curvesData = await r.json();
  } catch (e) {
    curvesData = null;
    $("curves-sections").replaceChildren(h("div", {class: "empty", text: "The curves could not be loaded."}));
  } finally {
    curvesLoading = false;
  }
  if (curvesData) renderCurves();
}

function curveMedians(rounds, name, k) {
  const out = [];
  for (let b = 0; b < k; b++) {
    const g = Object.values(rounds).map((r) => (r[name] ? r[name].gap[b] : null)).filter((x) => x != null).sort((a, c) => a - c);
    if (g.length < 10) continue;
    const m = g.length, med = m % 2 ? g[(m - 1) / 2] : (g[m / 2 - 1] + g[m / 2]) / 2;
    out.push({bin: b, med, n: m, best: g[m - 1]});
  }
  return out;
}

function renderCurves() {
  const box = $("curves-sections"), D = curvesData;
  const rounds = D.rounds || {}, bins = D.bins || {}, summary = D.summary || [];
  const nRounds = Object.keys(rounds).length;
  $("curves-sub").textContent = nRounds ? `${nRounds} stored rounds, ${summary.length} settings · local scores` : "";
  if (!summary.length) { box.replaceChildren(h("div", {class: "empty", text: "No curves yet: run test/bench/param_curves.py."})); return; }
  box.replaceChildren(...CURVE_GROUPS.map(([g, note]) => {
    const items = summary.filter((x) => x.group === g);
    if (!items.length) return null;
    return h("section", {class: "panel"},
      h("h2", {text: g}), h("div", {class: "muted", text: note}),
      h("div", {class: "curve-grid"}, items.map((x) => curveCard(x, rounds, bins[x.setting] || []))));
  }).filter(Boolean));
}

function curveCard(x, rounds, labels) {
  const k = labels.length, W = 320, H = 150, m = {t: 8, r: 10, b: 30, l: 28};
  const X = (b) => m.l + ((W - m.l - m.r) * b) / Math.max(1, k - 1);
  const y = (g) => m.t + (H - m.t - m.b) * (Math.max(CURVE_YMIN, Math.min(0, g)) / CURVE_YMIN);
  const svg = s("svg", {viewBox: `0 0 ${W} ${H}`, width: "100%", role: "img", "aria-label": `Response curve of ${x.setting}`});
  for (const t of [0, -1, -2, -3]) {
    svg.append(s("line", {x1: m.l, x2: W - m.r, y1: y(t), y2: y(t), class: "curve-grid-line"}));
    svg.append(s("text", {x: m.l - 5, y: y(t) + 3.5, "text-anchor": "end", class: "curve-tick"}, t === 0 ? "0" : "−" + Math.abs(t)));
  }
  const path = (pts) => pts.map((p, i) => `${i ? "L" : "M"}${X(p[0]).toFixed(1)},${y(p[1]).toFixed(1)}`).join("");
  for (const r of Object.values(rounds)) {
    const g = r[x.setting] ? r[x.setting].gap : [];
    const pts = g.map((v, b) => [b, v]).filter((p) => p[1] != null);
    if (pts.length > 1) svg.append(s("path", {d: path(pts), class: "curve-round"}));
  }
  const med = curveMedians(rounds, x.setting, k);
  if (med.length) {
    svg.append(s("path", {d: path(med.map((p) => [p.bin, p.med])), class: "curve-median"}));
    const top = med.reduce((a, c) => (c.med > a.med ? c : a));
    svg.append(s("circle", {cx: X(top.bin), cy: y(top.med), r: 4, class: "curve-peak"}));
  }
  const ticks = k <= 6 ? [...Array(k).keys()] : [0, Math.round((k - 1) / 2), k - 1];
  ticks.forEach((b, i) => svg.append(s("text", {x: X(b), y: H - m.b + 14, class: "curve-tick",
    "text-anchor": k <= 6 ? "middle" : i === 0 ? "start" : i === ticks.length - 1 ? "end" : "middle"}, String(labels[b]))));
  // hover: one invisible column per value with its numbers
  const step = (W - m.l - m.r) / Math.max(1, k - 1);
  for (const p of med) {
    const col = s("rect", {x: X(p.bin) - step / 2, y: m.t, width: Math.max(2, step), height: H - m.t - m.b, class: "curve-hit"});
    col.append(s("title", {}, `value ${labels[p.bin]}\nmedian ${p.med >= 0 ? "+" : ""}${p.med.toFixed(3)} (${p.n} rounds)\nbest round ${p.best >= 0 ? "+" : ""}${p.best.toFixed(3)}`));
    svg.append(col);
  }
  return h("div", {class: "curve-card"},
    h("div", {class: "curve-head"}, h("b", {text: x.setting}), chip(x.shape)),
    h("div", {class: "muted num curve-note"}, `Peak most often at ${x.top_peak} (${x.top_peak_rounds} of ${x.rounds} rounds) · middle half ${x.peak_low} to ${x.peak_high}`),
    svg);
}
