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
    else if (k === "style") Object.assign(el.style, v);
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
  const hash = location.hash === "#settings" ? "settings" : location.hash === "#history" ? "history" : "live";
  $("view-live").hidden = hash !== "live";
  $("view-history").hidden = hash !== "history";
  $("view-settings").hidden = hash !== "settings";
  $("tab-live").classList.toggle("on", hash === "live");
  $("tab-history").classList.toggle("on", hash === "history");
  $("tab-settings").classList.toggle("on", hash === "settings");
  if (hash === "settings") loadSettings();
  else if (hash === "history") loadHistory();
  else if (chart) setTimeout(sizeChart, 0);
}
window.addEventListener("hashchange", showTab);

// the Settings bar sticks just under the top bar, whose height depends on how it wraps
const topBar = document.querySelector(".top");
function setTopHeight() { document.documentElement.style.setProperty("--top-h", `${topBar.offsetHeight}px`); }
setTopHeight();
if (window.ResizeObserver) new ResizeObserver(setTopHeight).observe(topBar);

// ---- theme: the system setting until the toggle is used; the choice is remembered in this browser ----
const darkQuery = matchMedia("(prefers-color-scheme: dark)");
function effectiveTheme() {
  const t = document.documentElement.dataset.theme;
  return t === "dark" || t === "light" ? t : darkQuery.matches ? "dark" : "light";
}
function showThemeButton() {
  const t = effectiveTheme(), other = t === "dark" ? "light" : "dark";
  $("theme").textContent = other === "dark" ? "Dark" : "Light";
  $("theme").setAttribute("aria-label", `Switch to the ${other} theme`);
}
function themeChanged() {
  showThemeButton();
  if (last) render();                           // the chart takes its colors from the theme when it is drawn
  if (!$("view-history").hidden) renderHistoryCharts();
}
$("theme").addEventListener("click", () => {
  const next = effectiveTheme() === "dark" ? "light" : "dark";
  document.documentElement.dataset.theme = next;
  try { localStorage.setItem("minos-theme", next); } catch (e) { /* not remembered, still switched */ }
  themeChanged();
});
darkQuery.addEventListener("change", themeChanged);
showThemeButton();
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
  renderClimbers(E, finished, X);
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
    for (const [t, sc, why] of c.history || []) {           // GA children and the climber's own updates are different colors
      const kind = why === "ga" ? "ga" : why === "jump" ? "jump" : why.startsWith("C") ? "c"
        : why === "view" || why === "met" || why.endsWith("grid") || why.endsWith("line") ? "grid" : "";
      const m = s("circle", {class: kind ? `move ${kind}` : "move", cx: x(t), cy: y(i), r: kind === "grid" ? 2.5 : 3.5});
      m.append(s("title", {}, `${mins(t)}: ${fix(sc)} (${why})`));
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
  box.replaceChildren(...main.map((c) => {
    const w = c.wave;
    const growth = Object.entries(c.growth || {});
    const best = c.found != null ? c.found : c.best;
    return h("div", {class: `card ${climberClass(c)}`},
      h("div", {class: "head"}, h("b", {text: `Climber ${c.id}`}),
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
      h("div", {class: "origin", text: `From ${c.origin}.`}));
  }));
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

function renderStepA(E, R) {
  const box = $("waves");
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
