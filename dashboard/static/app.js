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
let chart = null, chartKey = "";
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
  const settings = location.hash === "#settings";
  $("view-live").hidden = settings;
  $("view-settings").hidden = !settings;
  $("tab-live").classList.toggle("on", !settings);
  $("tab-settings").classList.toggle("on", settings);
  if (settings) loadSettings();
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
  pc.className = `chip ${phase.split(" ")[0]}`;
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
  renderSchedule(E, finished);
  renderLineage(E);
  renderDomains(E, finished);
  renderMerges(E);
  renderWaves(E, R);
  renderScout(E, S, finished);
  renderChart(S, L.submissions || []);
  renderNodes(X, finished);
  renderFeed(L, E);
  renderSettingsState();
}

function kpi(label, value, note, cls) {
  return h("div", {class: `kpi ${cls || ""}`}, h("div", {class: "label", text: label}), h("div", {class: "value", text: value, title: value}),
    note ? h("div", {class: "note", text: note}) : null);
}

function renderKpis(L, R, S, E, X, finished) {
  const subs = L.submissions || [];
  const lastSub = subs.length ? subs[subs.length - 1] : null;
  const left = R.deadline_ts ? R.deadline_ts - now() : null;
  const live = E ? E.domains.filter((d) => d.status === "live").length : 0;
  const merging = E ? E.domains.filter((d) => d.status === "merging").length : 0;
  const perRun = S.gatk_runs_done ? (S.delivered || 0) / S.gatk_runs_done : null;
  const A = L.answer_set;
  $("kpis").replaceChildren(
    kpi("Time left", finished ? "done" : left != null ? dur(left) : "–",
      E ? `search ${dur(E.elapsed_s)} of ${dur(E.schedule.length_s)}` : (R.mode || ""), "k-time"),
    kpi("Best local", fix(S.best), E && E.pick ? `pick ${fix(E.pick.score)}` : S.best_config ? `config ${S.best_config.slice(0, 8)}` : "", "k-best"),
    kpi("Submitted", lastSub ? fix(lastSub[2]) : "–", lastSub ? `${subs.length} so far · ${ago(lastSub[0])}` : "nothing yet", "k-sub"),
    kpi("Configs scored", int(S.scored), `${int(S.evaluated)} run · ${int(S.failed)} failed · ${int(S.zero)} scored 0`, "k-configs"),
    finished ? kpi("Workers busy", "–", "round over", "k-workers")
      : kpi("Workers busy", X ? int(X.running) : "–", X ? `of ${int(X.workers)} · ${int(X.queued)} queued · ${int(S.in_flight)} in flight` : "", "k-workers"),
    finished ? kpi("Domains", E ? `${E.domains.filter((d) => d.status !== "retired").length}` : "–", E ? `at the end · ${E.domains.length} in the round` : "", "k-domains")
      : kpi("Domains", E ? `${live} / ${E.max_domains}` : "–",
        E ? `climbing · ${merging} merging in ${E.merges.length} merge${E.merges.length === 1 ? "" : "s"} · ${E.domains.length} so far` : "", "k-domains"),
    kpi("GATK runs", int(S.gatk_runs), perRun ? `${perRun.toFixed(2)} configs per run` : "", "k-gatk"),
    kpi("Answer set", A ? int(A.germline) : "building",
      A ? `germline + ${int(A.synthetic)} synthetic · ${A.regime || ""} · donor ${A.donor || "?"} · ${int(A.n_targets)} targets` : "", "k-answer"),
  );
}

function segments(E) {
  const sc = E.schedule;
  const founded = E.domains.length ? Math.min(...E.domains.map((d) => d.born_s)) : null;
  const openEnd = founded != null ? founded : Math.max(E.elapsed_s, Math.min(sc.opening_max_s, sc.reduce_to_3_start_s));
  const cut = (a, b) => [Math.max(0, a), Math.max(Math.max(0, a), Math.min(sc.length_s, b))];
  return [
    ["opening", "Opening", ...cut(0, openEnd)],
    ["climb", "Climb", ...cut(openEnd, sc.reduce_to_3_start_s)],
    ["merge", "6 → 3", ...cut(sc.reduce_to_3_start_s, sc.reduce_to_3_by_s)],
    ["climb", "Climb", ...cut(sc.reduce_to_3_by_s, sc.reduce_to_1_start_s)],
    ["merge", "3 → 1", ...cut(sc.reduce_to_1_start_s, sc.final_by_s)],
    ["final", "Final", ...cut(sc.final_by_s, sc.length_s)],
  ].filter((g) => g[3] - g[2] > 1);
}

function roundAxis(E) {
  // the round's time axis, shared by the schedule bar and the domain tree: from the round's start to its deadline.
  // The engine's times are seconds from its own start (`off` into the round); the search stops `margin` before
  // the deadline, which leaves main the time to submit.
  const R = (last && last.round) || {};
  const vals = (last && last.settings && last.settings.values) || {};
  const margin = vals["main.stop_before_deadline_s"] != null ? vals["main.stop_before_deadline_s"] : 60;
  const off = E && R.started_ts && E.started_ts ? Math.max(0, E.started_ts - R.started_ts) : 0;
  const total = R.started_ts && R.deadline_ts ? R.deadline_ts - R.started_ts : (E ? off + E.schedule.length_s + margin : 72 * 60);
  const searchEnd = E ? Math.min(total, off + E.schedule.length_s) : total - margin;
  return {total, off, searchEnd, start: R.started_ts};
}

function renderAnswerLane(ax, finished) {
  const lane = $("answer-lane");
  const A = last && last.answer_set;
  if (!ax.start) { lane.replaceChildren(); return; }
  if (A && !A.ready_ts) { lane.replaceChildren(); return; }               // a status file from before main recorded it
  const ready = A && A.ready_ts ? A.ready_ts - ax.start : null;
  const end = ready != null ? ready : Math.max(5, now() - ax.start);
  const pct = (t) => `${(100 * Math.max(0, Math.min(ax.total, t))) / ax.total}%`;
  const label = ready != null ? `answer set ${Math.round(ready)} s` : "answer set…";
  lane.replaceChildren(h("div", {class: "seg", style: {left: "0", width: pct(end)}, title: label, text: label}));
}

function renderSchedule(E, finished) {
  const box = $("schedule"), ticks = $("schedule-ticks");
  const ax = roundAxis(E);
  renderAnswerLane(ax, finished);
  if (!E) {
    const R = (last && last.round) || {};
    const what = R.mode === "building answer set" ? "building the answer set, then the opening" : "waiting for the fleet to settle, then the opening";
    box.replaceChildren(h("div", {class: "seg opening", style: {left: "0", width: "100%"}, text: what}));
    ticks.replaceChildren();
    $("schedule-sub").textContent = "";
    return;
  }
  const L = ax.total, o = ax.off;
  const x = (t) => `${(100 * Math.max(0, Math.min(L, t))) / L}%`;
  const segs = segments(E).map(([cls, label, a, b]) => [cls, label, a + o, b + o]);
  segs.push(["submit", "submit", ax.searchEnd, L]);
  const els = segs.filter((g) => g[3] > g[2]).map(([cls, label, a, b]) => h("div", {class: `seg ${cls}`,
    style: {left: x(a), width: `${(100 * (b - a)) / L}%`}, title: `${label}: ${mins(a)} to ${mins(b)}`}, label));
  const t = finished ? L : E.elapsed_s + o;
  els.push(h("div", {class: "done", style: {width: x(t)}}));
  els.push(h("div", {class: "pick", style: {left: x(E.schedule.pick_s + o)}, title: "plateau-centre pick"}));
  if (!finished) els.push(h("div", {class: "now", style: {left: x(t)}, title: `now: ${mins(t)}`}));
  box.replaceChildren(...els);
  const all = [...new Set([0, ...segs.map((g) => Math.round(g[3] / 60))])].sort((a, b) => a - b);
  const marks = [];
  for (const m of all) {                                    // labels closer than 2 minutes: keep the later one
    if (marks.length && m - marks[marks.length - 1] < 2) marks[marks.length - 1] = m; else marks.push(m);
  }
  ticks.replaceChildren(...marks.map((m) => h("span", {style: {left: x(m * 60)}, text: `${m}`})));
  $("schedule-sub").textContent = `minutes from the round's start; the search's stages scaled ×${E.schedule.scale.toFixed(2)} to this round`;
}

function rowClass(d) {
  if (d.status === "retired") return "retired";
  if (d.status === "merging") return "merging";
  return d.state === "climb" ? "live" : d.state;
}

function renderLineage(E) {
  const box = $("lineage");
  if (!E || !E.domains.length) {
    box.replaceChildren(h("div", {class: "empty", text: E ? "No domains yet: they are founded when the opening ends." : "The engine has not started."}));
    return;
  }
  const doms = [...E.domains].sort((a, b) => a.born_s - b.born_s || a.id - b.id);
  const row = new Map(doms.map((d, i) => [d.id, i]));
  const W = Math.max(520, box.clientWidth || 640), top = 22, rh = 22;
  const left = parseFloat(css("--axis-l")) || 46, right = parseFloat(css("--axis-r")) || 66;   // same axis as the schedule
  const H = top + doms.length * rh + 8;
  const ax = roundAxis(E), L = ax.total, o = ax.off;                   // the same axis as the schedule bar
  const x = (t) => left + ((W - left - right) * Math.max(0, Math.min(L, t + o))) / L;
  const y = (i) => top + i * rh + rh / 2;
  const svg = s("svg", {width: W, height: H, viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Domain tree over the round"});
  for (const [cls, label, a, b] of segments(E)) {
    svg.append(s("rect", {class: `band ${cls}`, x: x(a), y: 0, width: Math.max(0, x(b) - x(a)), height: H}));
    svg.append(s("text", {x: x(a) + 4, y: 13}, label));
  }
  const tNow = E.elapsed_s;
  for (const d of doms) {
    const i = row.get(d.id);
    for (const p of d.parents) {
      if (!row.has(p)) continue;
      const xp = x(d.born_s);
      svg.append(s("path", {class: "link", d: `M ${xp - 6} ${y(row.get(p))} L ${xp} ${y(row.get(p))} L ${xp} ${y(i)}`}));
    }
  }
  for (const d of doms) {
    const i = row.get(d.id);
    const end = d.ended_s != null ? d.ended_s : tNow;
    const r = s("rect", {class: `row ${rowClass(d)}`, x: x(d.born_s), y: y(i) - 6, width: Math.max(3, x(end) - x(d.born_s)), height: 12, rx: 3});
    r.append(s("title", {}, `Domain ${d.id}: best ${fix(d.best)}, ${d.status}, ${d.state}/${d.mode}; born at ${mins(d.born_s)} (${d.origin})` +
      (d.ended_s != null ? `; retired at ${mins(d.ended_s)}` : "")));
    svg.append(r);
    if (d.origin.startsWith("a scout find")) svg.append(s("circle", {class: "find", cx: x(d.born_s), cy: y(i), r: 4}));
    svg.append(s("text", {x: 4, y: y(i) + 4}, `D${d.id}`));
    svg.append(s("text", {class: "score", x: x(end) + 5, y: y(i) + 4}, fix(d.best, 3)));
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

function renderDomains(E, finished) {
  const box = $("domains");
  if (!E) { box.replaceChildren(h("div", {class: "empty", text: "Waiting for the engine."})); $("domains-sub").textContent = ""; return; }
  const live = E.domains.filter((d) => d.status !== "retired").sort((a, b) => a.id - b.id);   // stable order, as in the tree
  $("h-domains").firstChild.textContent = finished ? "Final domains " : "Live domains ";
  const gone = E.domains.filter((d) => d.status === "retired");
  $("domains-sub").textContent = gone.length ? `${gone.length} retired` : "";
  if (!live.length) { box.replaceChildren(h("div", {class: "empty", text: E.phase === "opening" ? "The opening is still running." : "No live domain."})); return; }
  box.replaceChildren(...live.map((d) => {
    const cls = d.status === "merging" ? "merging" : d.state === "final" ? "final" : "";
    return h("div", {class: `card ${cls}`},
      h("div", {class: "head"}, h("b", {text: `Domain ${d.id}`}),
        finished ? chip("finished", "final") : chip(d.status === "merging" ? "merging" : d.state, d.status === "merging" ? "merging" : d.state),
        finished ? null : chip(d.mode)),
      h("div", {class: "big", text: fix(d.best)}),
      bar(d.allotted ? d.inflight / Math.max(d.allotted, d.inflight) : 0, "", null),
      kv([
        ["climbing", `${d.progress_per_min >= 0 ? "+" : ""}${fix(d.progress_per_min)} / min`],
        finished ? null : ["in flight", d.status === "merging" ? `${d.inflight}, merging` : `${d.inflight} of ${fix(d.allotted, 0)} allotted`],
        ["configs tried", int(d.tried)],
        finished || !d.chains ? null : ["chains", `${d.chains}${d.perturb_rounds ? ` · perturbation round ${d.perturb_rounds}` : ""}`],
        ["born at", mins(d.born_s)],
        d.parents.length ? ["merged from", d.parents.map((p) => `D${p}`).join(", ")] : null,
        ["centre", h("span", {class: "mono", text: d.centre.slice(0, 12)})],
      ]),
      h("div", {class: "origin", text: `From ${d.origin}.`}));
  }));
}

function renderMerges(E) {
  const box = $("merges");
  if (!E || !E.merges.length) { box.replaceChildren(h("div", {class: "empty", text: "No merge running."})); return; }
  box.replaceChildren(...E.merges.map((m) => h("div", {class: "rowline"},
    h("div", {class: "top-line"}, h("span", {}, h("b", {text: `Merge ${m.id}`}), " · ", m.parents.map((p) => `D${p}`).join(" + ")),
      h("span", {class: "num muted", text: `${m.back} / ${m.sent} back`})),
    bar(m.sent ? m.back / m.sent : 0, "merge", 0.8),
    h("div", {class: "muted", text: `${m.why}; best so far ${fix(m.best)}; started at ${mins(m.started_s)}`}))));
}

function renderWaves(E, R) {
  const box = $("waves");
  if (!E || !E.waves.length) { box.replaceChildren(h("div", {class: "empty", text: R && R.mode ? R.mode : "Not started."})); return; }
  box.replaceChildren(...E.waves.map((w, i) => {
    const next = E.waves[i + 1];
    const gain = next && w.best_at_launch > 0 ? next.best_at_launch - w.best_at_launch : null;
    return h("div", {class: "rowline"},
      h("div", {class: "top-line"}, h("b", {text: w.wave === 0 ? "Wave 0: seeds + hypercube" : `Generation ${w.wave}`}),
        h("span", {class: "num muted", text: `${int(w.back)} / ${int(w.sent)} back`})),
      bar(w.sent ? w.back / w.sent : 0, "open", null),
      gain != null ? h("div", {class: "muted", text: `the next wave started ${gain >= 0 ? "+" : ""}${fix(gain)} higher`}) : null);
  }), h("div", {class: "muted", text: E.phase === "opening" ? "Domains are founded when the opening ends." : `Founded ${E.domains.filter((d) => d.origin === "the opening").length} domains.`}));
}

function renderScout(E, S, finished) {
  const box = $("scout");
  if (!E) { box.replaceChildren(h("div", {class: "empty", text: "Waiting for the engine."})); return; }
  const target = Math.max(1, Math.floor(E.capacity * E.scout.share));
  const finds = E.domains.filter((d) => d.origin.startsWith("a scout find")).length;
  box.replaceChildren(
    h("div", {class: "rowline"},
      h("div", {class: "top-line"}, h("b", {text: "Scout"}),
        h("span", {class: "num muted", text: finished ? "stopped" : `${int(E.scout.in_flight)} / ${int(target)} in flight`})),
      bar(finished ? 0 : E.scout.in_flight / target, "scout", null)),
    kv([
      ["finds tested", int(E.scout.finds_tested)],
      ["finds being tested", int(E.scout.finds_testing)],
      ["domains founded from finds", int(finds)],
      ["scored candidates kept", int(E.candidates)],
    ]),
    h("div", {class: "rowline"}, h("b", {text: "Pick"}),
      E.pick ? kv([["config", h("span", {class: "mono", text: E.pick.key.slice(0, 12)})], ["its score", fix(E.pick.score)],
        ["best score", fix(E.pick.best)], ["neighbour mean", fix(E.pick.neighbour_mean)]])
        : h("div", {class: "muted", text: `made ${mins(E.schedule.length_s - E.schedule.pick_s)} before the end`})),
  );
}

// ---- chart ----------------------------------------------------------------------------------------------------------
function chartData(S, subs) {
  const pts = (S.points || []).filter((p) => p[1] > 0);
  const best = S.best_timeline || [];
  const map = new Map();
  const put = (t, i, v) => { const k = Math.round(t * 10) / 10; if (!map.has(k)) map.set(k, [null, null, null]); map.get(k)[i] = v; };
  for (const [t, v] of pts) put(t, 0, v);
  for (const [t, v] of best) put(t, 1, v);
  for (const [t, , v] of subs) put(t, 2, v);
  const xs = [...map.keys()].sort((a, b) => a - b);
  return [xs, xs.map((t) => map.get(t)[0]), xs.map((t) => map.get(t)[1]), xs.map((t) => map.get(t)[2])];
}
function sizeChart() {
  const el = $("chart");
  if (chart && el.clientWidth) chart.setSize({width: el.clientWidth, height: 300});
}
function renderChart(S, subs) {
  const el = $("chart");
  if (typeof uPlot === "undefined") return;
  const data = chartData(S, subs);
  if (!data[0].length) { if (chart) { chart.destroy(); chart = null; } el.replaceChildren(h("div", {class: "empty", text: "No scores yet."})); return; }
  const best = Math.max(...data[2].filter((v) => v != null), ...data[1].filter((v) => v != null));
  const theme = effectiveTheme();
  if (!chart || chartKey !== theme) {
    if (chart) chart.destroy();
    el.replaceChildren();
    chartKey = theme;
    const axis = {stroke: css("--muted"), grid: {stroke: css("--line"), width: 1}, ticks: {stroke: css("--line")}};
    chart = new uPlot({
      width: el.clientWidth || 800, height: 300,
      scales: {x: {time: true, range: (u, lo, hi) => (hi - lo < 120 ? [lo - 60, hi + 60] : [lo, hi])},
               y: {range: (u, lo, hi) => {
        const top = Math.max(hi, best);
        return [Math.max(lo, top - 3), top + 0.15];
      }}},
      axes: [{...axis}, {...axis, size: 60, values: (u, v) => v.map((x) => x.toFixed(2))}],
      series: [
        {},
        {label: "scored config", stroke: css("--faint"), paths: () => null, points: {show: true, size: 3, fill: css("--faint"), stroke: css("--faint")}},
        {label: "best so far", stroke: css("--accent"), width: 2, spanGaps: true, points: {show: false}, paths: uPlot.paths.stepped({align: 1})},
        {label: "submitted", stroke: css("--final"), paths: () => null, points: {show: true, size: 9, fill: css("--final"), stroke: css("--final")}},
      ],
    }, data, el);
  } else {
    chart.setData(data);
  }
}
window.addEventListener("resize", () => { sizeChart(); if (last) renderLineage((last.search || {}).engine); });

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

function renderFeed(L, E) {
  const rows = [];
  for (const [t, text] of (E && E.events) || []) rows.push({t, kind: "engine", cls: "", text});
  for (const c of ((L.settings || {}).changes) || []) rows.push({t: c.t, kind: "setting", cls: "scout", text: `${c.name}: ${JSON.stringify(c.old)} → ${JSON.stringify(c.new)}`});
  for (const [t, level, name, msg] of L.log || []) rows.push({t, kind: level.toLowerCase(), cls: level === "ERROR" ? "bad" : "warn", text: `${name}: ${msg}`});
  rows.sort((a, b) => b.t - a.t);
  const box = $("feed");
  if (!rows.length) { box.replaceChildren(h("div", {class: "empty", text: "Nothing yet."})); return; }
  box.replaceChildren(...rows.slice(0, 120).map((r) => h("div", {class: "ev"},
    h("span", {class: "t", text: clock(r.t), title: E ? `${mins(r.t - E.started_ts)} into the search` : ""}),
    h("span", {}, chip(r.kind, r.cls)), h("span", {class: "msg", text: r.text}))));
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

// ---- start ----------------------------------------------------------------------------------------------------------
showTab();
poll();
