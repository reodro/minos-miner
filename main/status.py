"""What the dashboard shows: main writes a live snapshot and a per-round history; the dashboard only reads them.

    <status_dir>/live.json                  the current round, rewritten every few seconds (atomic replace);
                                            "schema" 2: + "settings" (common.tuning: in force, pending, rejected,
                                            changes) and the engine's notes under search.engine (phase, schedule,
                                            waves, every domain with origin and parents, merges, scout, events)
    <status_dir>/rounds.jsonl               one "round" line per finished run (sections answer_set / evaluator /
                                            search), later "review" and "rank" lines from main.review
    <status_dir>/runs/<run_id>/             what main.review needs once the truth is released:
        manifest.json                       final + baseline config, every config's params and local score,
                                            the settings at the start and every change made during the run
        answer_set.json.gz                  our answer set as used at the end of the run

A run is one mining of one round ("<tag>-<start epoch>"): a stored round can be replayed many times.

Nothing here can slow down or break mining: every method swallows its own errors. Node names are anonymized
("node 1", ...) because the dashboard is public.
"""
from __future__ import annotations

import collections
import gzip
import json
import logging
import os
import threading
import time
from pathlib import Path
from typing import Optional

from common.settings import SETTINGS
from common.tuning import TUNING
from common.types import CLASSES

SCHEMA = 2                                # live.json layout version (the dashboard checks it)
TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

log = logging.getLogger(__name__)
STATUS_DIR = Path(SETTINGS.status_dir)


class _Recent(logging.Handler):
    """Keeps the last warnings/errors for the dashboard's log panel."""

    def __init__(self, n: int = 40):
        super().__init__(logging.WARNING)
        self.lines = collections.deque(maxlen=n)

    def emit(self, record):
        try:
            self.lines.append([record.created, record.levelname, record.name, record.getMessage()[:300]])
        except Exception:  # noqa: BLE001
            pass


class StatusBoard:
    def __init__(self, directory: Path = STATUS_DIR, every_s: float = 5.0):
        self.dir, self.every_s = directory, every_s
        self.dir.mkdir(parents=True, exist_ok=True)
        self._lock = threading.Lock()
        self._last = 0.0
        self.state: dict = {}
        self.kind = "live"                                  # "live" | "replay" (set by main.__main__)
        self._client = None
        self._answer = None
        self._baseline: Optional[str] = None
        self._baseline_cfg: Optional[dict] = None
        self._util: list = []                               # running / workers samples for average utilization
        self._ex0: Optional[dict] = None                    # executor counters at the round's first tick
        self._seen: set = set()                             # config keys already plotted
        self._tuning0: Optional[dict] = None                # the settings in force when the round started
        self._recent = _Recent()
        logging.getLogger().addHandler(self._recent)

    # ---- round lifecycle (called by mine_round) -------------------------------------------------------------------
    def round_started(self, round, deadline_ts: Optional[float], baselines=()) -> None:
        self._safe(self._started, round, deadline_ts, baselines)

    def _started(self, round, deadline_ts, baselines):
        self._client, self._answer, self._util, self._ex0, self._seen = None, None, [], None, set()
        self._tuning0 = dict(TUNING.values)
        self._baseline = baselines[0].key if baselines else None
        self._baseline_cfg = baselines[0].to_dict() if baselines else None
        now = time.time()
        self.state = {"round": {"id": round.round_id, "tag": round.tag, "region": round.region,
                                "run_id": f"{round.tag}-{int(now)}", "kind": self.kind,
                                "deadline_ts": deadline_ts, "started_ts": now, "mode": "building answer set"},
                      "answer_set": None, "search": {"evaluated": 0, "failed": 0, "in_flight": 0, "best": None,
                                                     "scored": 0, "zero": 0, "gatk_runs": 0, "points": [],
                                                     "best_timeline": [], "eval_timeline": []},
                      "submissions": [], "executor": None, "schema": SCHEMA, "settings": TUNING.snapshot()}
        self._write(force=True)

    def answer_set(self, a) -> None:
        self._safe(self._answer_set, a)

    def _answer_set(self, a):
        self._answer = a
        n = a.notes or {}
        self.state["answer_set"] = {
            "regime": a.regime.value, "donor": a.donor, "germline": len(a.germline), "synthetic": len(a.synthetic),
            "n_targets": a.n_targets, "n_targets_source": n.get("n_targets_source"), "seconds": n.get("seconds"),
            "degraded": n.get("degraded"), "version": a.version,
            "composition": {"germline": _composition(a.germline),
                            "synthetic": _composition([c.variant for c in a.synthetic], [c.p for c in a.synthetic])}}
        self.state["round"]["mode"] = "searching" if a.germline else "static fallback"
        self._write(force=True)

    def tick(self, client, out, executor) -> None:
        """Called from the round loop; writes at most every `every_s` seconds."""
        now = time.time()
        if now - self._last < self.every_s:
            return
        self._safe(self._tick, client, out, executor, now)

    def _tick(self, client, out, executor, now):
        self._client = client
        s = self.state.setdefault("search", {})
        best = client.best()
        full = {k: r for (k, fid), r in client.latest.items() if fid == "full"}
        for k, r in full.items():                            # a dot per newly scored config (5 s resolution)
            if k not in self._seen and r.score is not None:
                self._seen.add(k)
                s.setdefault("points", []).append([now, round(r.score, 4), k])
        s["points"] = s["points"][-5000:]
        s["failure_reasons"] = sorted(out.failures.items(), key=lambda kv: -kv[1])[:8]
        s["engine"] = (getattr(client, "notes", None) or {}).get("engine")
        s.update(evaluated=out.evaluated, failed=out.failed, in_flight=client.in_flight,
                 scored=sum(1 for r in full.values() if r.score is not None),
                 zero=sum(1 for r in full.values() if r.score is not None and r.score <= 0),
                 best=best.score if best else None, best_config=best.config_key if best else None,
                 spent=round(client.spent, 2), sent=len(client.configs), scores=_scores(client),
                 best_detail=_detail(best, client))
        if best:
            s.setdefault("best_timeline", []).append([now, best.score])
        s.setdefault("eval_timeline", []).append([now, s["scored"]])
        for k in ("best_timeline", "eval_timeline"):
            s[k] = s[k][-2000:]
        self.state["submissions"] = [[t, key, score] for t, key, score in out.submissions]
        self.state["settings"] = TUNING.snapshot()
        stats = getattr(executor, "stats", None)
        ex = stats() if stats else None
        if ex is not None:
            if self._ex0 is None:                            # executor counters live as long as the process
                self._ex0 = {k: ex.get(k) or 0 for k in ("gatk_runs", "gatk_runs_done", "batches_run",
                                                          "repairs_run", "done")}
            ex["round"] = {k: (ex.get(k) or 0) - v for k, v in self._ex0.items()}
            s["gatk_runs"] = ex["round"]["gatk_runs"]
            s["gatk_runs_done"] = ex["round"]["gatk_runs_done"]
            s["delivered"] = ex["round"]["done"]
        self.state["executor"] = ex
        if ex and ex.get("workers"):
            self._util.append(ex.get("running", 0) / ex["workers"])
        self._write()

    def round_ended(self, mined) -> None:
        self._safe(self._ended, mined)

    def _ended(self, mined):
        r = self.state.get("round", {})
        out, client = mined.outcome, self._client
        r["mode"] = f"finished ({mined.mode})"
        srch = self.state.setdefault("search", {})
        if out:                                  # final numbers, not the last throttled tick
            srch.update(evaluated=out.evaluated, failed=out.failed, in_flight=0,
                        best=out.best.score if out.best else None)
            self.state["submissions"] = [[t, key, score] for t, key, score in out.submissions]
        if client:
            srch.update(spent=round(client.spent, 2), sent=len(client.configs), scores=_scores(client))
            full = {k: res for (k, fid), res in client.latest.items() if fid == "full" and res.score is not None}
            for k, res in full.items():                      # results that landed after the last tick
                if k not in self._seen:
                    self._seen.add(k)
                    srch.setdefault("points", []).append([time.time(), round(res.score, 4), k])
            srch.update(scored=len(full), zero=sum(1 for res in full.values() if res.score <= 0))
        a = self.state.get("answer_set") or {}
        final = out.submissions[-1][1] if out and out.submissions else (self._baseline if mined.mode == "static"
                                                                        else None)
        ex = self.state.get("executor") or {}
        started, ended = r.get("started_ts") or time.time(), time.time()
        evaluated = out.evaluated if out else 0
        scores = srch.get("scores") or []
        local = {k: res.score for (k, fid), res in (client.latest.items() if client else []) if fid == "full"}
        best_tl = srch.get("best_timeline") or []
        best = out.best.score if out and out.best else None
        row = {"type": "round", "run_id": r.get("run_id"), "kind": r.get("kind"), "round_id": mined.round_id,
               "tag": r.get("tag"), "region": r.get("region"), "mode": mined.mode, "note": mined.note,
               "started": started, "ended": ended,
               "answer_set": {"regime": a.get("regime"), "donor": a.get("donor"), "n_targets": a.get("n_targets"),
                              "n_targets_source": a.get("n_targets_source"),
                              "build_s": round(sum((a.get("seconds") or {}).values()), 1) if a.get("seconds") else None,
                              "degraded": a.get("degraded"), "composition": a.get("composition")},
               "evaluator": {"kind": ex.get("kind"), "evaluated": evaluated, "failed": out.failed if out else 0,
                             "per_min": round(evaluated / max((ended - started) / 60, 1e-9), 2),
                             "gatk_runs": srch.get("gatk_runs"),
                             "per_gatk_run": round(srch["delivered"] / srch["gatk_runs_done"], 2)
                             if srch.get("gatk_runs_done") and srch.get("delivered") else None,
                             "batches": (ex.get("round") or {}).get("batches_run"),
                             "repairs": (ex.get("round") or {}).get("repairs_run"),
                             "nodes": len(ex.get("nodes") or []) or None,
                             "utilization": round(sum(self._util) / len(self._util), 3) if self._util else None},
               "search": {"best_local": best, "final_config": final, "baseline_config": self._baseline,
                          "baseline_local": local.get(self._baseline),
                          "time_to_best_min": next((round((t - started) / 60, 1) for t, v in best_tl
                                                    if best is not None and v >= best - 1e-9), None),
                          "zero_scored": sum(1 for x in scores if x <= 0), "scored": len(scores),
                          "submissions": len(out.submissions) if out else 0}}
        self._append(row)
        self._save_run(r.get("run_id"), row, client, final, out)
        if ex:                                                # the round is over: nothing runs on any node for it
            ex.update(running=0, queued=0)
            for n in ex.get("nodes") or []:
                n.update(busy=0, jobs=0, small=0)
        self._write(force=True)

    def _save_run(self, run_id, row, client, final, out=None):
        """What main.review needs later: the configs (params + local score), the configs that failed (what to
        exclude from the space) and our answer set."""
        if not run_id:
            return
        d = self.dir / "runs" / run_id
        d.mkdir(parents=True, exist_ok=True)
        configs = {}
        if client:
            for (k, fid), res in client.latest.items():
                cfg = client.configs.get(k)
                if fid == "full" and cfg is not None:
                    configs[k] = {"config": cfg.to_dict(), "local": res.score, "gate_passed": res.gate_passed}
        if self._baseline and self._baseline not in configs and self._baseline_cfg:   # e.g. a static fallback
            configs[self._baseline] = {"config": self._baseline_cfg, "local": None, "gate_passed": None}
        failed = dict(getattr(out, "failed_configs", None) or {})
        manifest = {k: row[k] for k in ("run_id", "kind", "round_id", "tag", "region", "mode", "started", "ended")}
        manifest.update(final_config=final, baseline_config=self._baseline, configs=configs, failed=failed)
        manifest["tuning"] = {"at_start": self._tuning0,
                              "changes": [c for c in TUNING.snapshot()["changes"] if c["t"] >= (row.get("started") or 0)]}
        (d / "manifest.json").write_text(json.dumps(manifest))
        a = client.answer_set if client else self._answer
        if a is not None:
            with gzip.open(d / "answer_set.json.gz", "wt") as f:
                json.dump({"version": a.version, "regime": a.regime.value, "donor": a.donor, "n_targets": a.n_targets,
                           "germline": [[v.chrom, v.pos, v.ref, v.alt, v.gt] for v in a.germline],
                           "synthetic": [[c.variant.chrom, c.variant.pos, c.variant.ref, c.variant.alt, c.variant.gt,
                                          c.p] for c in a.synthetic]}, f)

    def append(self, row: dict) -> None:
        """Extra history lines (main.review: "review" and "rank")."""
        self._safe(self._append, row)

    # ---- io -----------------------------------------------------------------------------------------------------------
    def _write(self, force: bool = False):
        with self._lock:
            self.state["updated"] = time.time()
            self.state["log"] = list(self._recent.lines)
            tmp = self.dir / "live.json.tmp"
            tmp.write_text(json.dumps(self.state))
            os.replace(tmp, self.dir / "live.json")
            self._last = time.time()

    def _append(self, row: dict):
        with self._lock, open(self.dir / "rounds.jsonl", "a") as f:
            f.write(json.dumps(row) + "\n")

    @staticmethod
    def _safe(fn, *args):
        try:
            fn(*args)
        except Exception as e:  # noqa: BLE001 - the dashboard must never break mining
            log.debug("status board: %s", e)


# ---- summaries (plain python; the dashboard only formats them) -----------------------------------------------------
def _composition(variants, ps=None) -> dict:
    """SNP/indel split, genotype counts (het/hom), het:hom ratio, Ti/Tv of SNPs; for synthetic also the confidence."""
    snp = indel = het = hom = other = ti = tv = 0
    for v in variants:
        alts = v.alt.split(",")
        if len(v.ref) == 1 and all(len(x) == 1 for x in alts):
            snp += 1
            if len(alts) == 1:
                ti, tv = (ti + 1, tv) if (v.ref, v.alt) in TRANSITIONS else (ti, tv + 1)
        else:
            indel += 1
        a = v.gt.replace("|", "/").split("/")
        if len(a) == 2 and a[0] == a[1] and a[0] not in (".", "0"):
            hom += 1
        elif len(a) == 2 and a[0] != a[1]:
            het += 1
        else:
            other += 1
    c = {"total": len(variants), "snp": snp, "indel": indel, "het": het, "hom": hom, "other_gt": other,
         "het_hom": round(het / hom, 3) if hom else None, "titv": round(ti / tv, 3) if tv else None}
    if ps is not None:
        c["confident"] = sum(1 for p in ps if p >= 0.99)
        c["mean_p"] = round(sum(ps) / len(ps), 3) if ps else None
    return c


def _scores(client) -> list:
    """Every scored config's current score (full fidelity), for the distribution chart."""
    return sorted(round(r.score, 3) for (k, fid), r in client.latest.items()
                  if fid == "full" and r.score is not None)[-5000:]


def _detail(best, client) -> Optional[dict]:
    if best is None:
        return None
    per = {}
    for c in CLASSES:
        k = best.counts.get(c)
        if k is None:
            continue
        p = k.tp / (k.tp + k.fp) if k.tp + k.fp else None
        r = k.tp / (k.tp + k.fn) if k.tp + k.fn else None
        f1 = 2 * p * r / (p + r) if p and r else 0.0
        per[c] = {"tp": k.tp, "fp": k.fp, "fn": k.fn, "f1": round(f1, 4)}
    cfg = client.configs.get(best.config_key)
    return {"config": best.config_key, "score": best.score, "core": best.core, "gate_passed": best.gate_passed,
            "region_fp": best.region_fp, "n_targets": best.n_targets,
            "fp_per_target": round(best.region_fp / best.n_targets, 3) if best.n_targets else None,
            "score_sd": best.score_sd, "ratios": dict(best.ratios), "classes": per,
            "params": cfg.params if cfg else None}
