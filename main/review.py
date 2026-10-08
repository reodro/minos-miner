"""Post-round review: once a round's truth is released, judge every part of what main did in it. Owner: main.

Runs next to the live miner (a background thread on main, which is otherwise idle: the fleet does the evaluating)
and every `poll_s` (5 min) looks at the recent runs in <status_dir>/rounds.jsonl:

    truth released?  -> "review" line: answer-set accuracy (ours vs the released truth: strict record match, and
                        "as_scored" - our answer set scored like a callset, so only what the score can see counts),
                        and the TRUE score of
                        the final config, the baseline and the top-k configs by local score (their VCFs are reused
                        when on this machine, else GATK re-runs them here - exact, so the result is what the
                        validator saw), with Spearman(local, true), pick gap and local - true of the final.
    round finalized? -> "rank" line: the final config's rank among the real submissions, gap to the field's best.

A replay has both at once (the round is long finished), so `python -m main replay` reviews right after mining.

    python -m main review <run_id|tag>      review (and rank) one stored run now, e.g. to backfill
"""
from __future__ import annotations

import gzip
import hashlib
import json
import logging
import os
import statistics
import threading
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from common.settings import SETTINGS, default_workers
from common.types import FULL, CallSet, Config, EvalJob, Round, Variant

log = logging.getLogger(__name__)

STATUS_DIR = Path(SETTINGS.status_dir)
GIVE_UP_REVIEW_S = 24 * 3600            # truth should be out one round (72 min) after the deadline
GIVE_UP_RANK_S = 72 * 3600              # finalization takes a few rounds


# ---- pure parts (unit-tested) ----------------------------------------------------------------------------------------
def _gt(gt: str) -> str:
    return "/".join(sorted(gt.replace("|", "/").split("/"), key=lambda a: (a == ".", a)))


def _pr(ours: set, truth: set) -> dict:
    hit = len(ours & truth)
    return {"ours": len(ours), "truth": len(truth), "match": hit,
            "precision": round(hit / len(ours), 4) if ours else None,
            "recall": round(hit / len(truth), 4) if truth else None}


def answer_set_accuracy(germline: Iterable[Variant], synthetic: Iterable[Tuple[Variant, float]],
                        truth_germline: Iterable[Variant], truth_synthetic: Iterable[Variant],
                        n_targets: Optional[int], truth_n: Optional[int]) -> dict:
    """Ours vs the released truth.
    germline: exact records (chrom, pos, ref, alt, unphased GT) - what the score treats as exact - and sites alone;
    synthetic: by allele (chrom, pos, ref, alt), all candidates and the confident ones (p >= 0.99); N: estimate vs real."""
    site = lambda v: (v.chrom, v.pos, v.ref, v.alt)                       # noqa: E731
    exact = lambda v: (v.chrom, v.pos, v.ref, v.alt, _gt(v.gt))          # noqa: E731
    germline, synthetic = list(germline), list(synthetic)
    tg, ts = list(truth_germline), list(truth_synthetic)
    return {"germline_exact": _pr({exact(v) for v in germline}, {exact(v) for v in tg}),
            "germline_sites": _pr({site(v) for v in germline}, {site(v) for v in tg}),
            "synthetic": _pr({site(v) for v, _ in synthetic}, {site(v) for v in ts}),
            "synthetic_confident": _pr({site(v) for v, p in synthetic if p >= 0.99}, {site(v) for v in ts}),
            "synthetic_gt": _pr({exact(v) for v, _ in synthetic}, {exact(v) for v in ts}),
            "n_targets": {"ours": n_targets, "truth": truth_n,
                          "ratio": round(n_targets / truth_n, 4) if n_targets and truth_n else None}}


def write_answer_set_vcf(germline: Iterable[Variant], synthetic: Iterable[Tuple[Variant, float]],
                         contigs: Dict[str, int], path: Path) -> Path:
    """Our answer set as a bgzipped, indexed VCF: every record with its genotype, germline and synthetic alike - all
    the score sees of a truth file (the SYNTHETIC flag is written for reading; the score never uses it)."""
    import pysam
    recs = [(v, False) for v in germline] + [(v, True) for v, _ in synthetic]
    order = {c: i for i, c in enumerate(contigs)}
    recs.sort(key=lambda t: (order.get(t[0].chrom, len(order)), t[0].pos, t[0].ref, t[0].alt))
    path = Path(path)
    text = path.with_name(path.name[:-len(".gz")] if path.name.endswith(".gz") else path.name + ".txt")
    with open(text, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.writelines(f"##contig=<ID={c},length={n}>\n" for c, n in contigs.items())
        f.write('##INFO=<ID=SYNTHETIC,Number=0,Type=Flag,Description="Predicted injected variant">\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tANSWER\n")
        for v, syn in recs:
            f.write(f"{v.chrom}\t{v.pos}\t.\t{v.ref}\t{v.alt}\t.\tPASS\t{'SYNTHETIC' if syn else '.'}\tGT\t{v.gt}\n")
    pysam.tabix_compress(str(text), str(path), force=True)
    pysam.tabix_index(str(path), preset="vcf", force=True)
    text.unlink()
    return path


def answer_set_as_scored(result) -> dict:
    """Our answer set scored like a callset against the released truth, through the validator's own path: what its
    errors are worth to the score. Spellings of one variant match, labels are ignored, genotypes count. `core` is the
    difficulty-weighted F1 (1.0 = no error the score can see); `core_points` = what its gap is worth in score points."""
    from common.scoring import V2_CORE_WEIGHT
    classes = {}
    for k, c in result.counts.items():
        f1 = 2 * c.tp / (2 * c.tp + c.fp + c.fn) if c.tp + c.fp + c.fn else None
        classes[k] = {"tp": c.tp, "fn": c.fn, "fp": c.fp, "f1": None if f1 is None else round(f1, 5)}
    core = result.core
    return {"core": None if core is None else round(core, 5),
            "missed": sum(c.fn for c in result.counts.values()), "extra": sum(c.fp for c in result.counts.values()),
            "core_points": None if core is None else round(100 * V2_CORE_WEIGHT * (1 - core), 4),
            "classes": classes}


def pick_configs(configs: Dict[str, dict], final: Optional[str], baseline: Optional[str], top_k: int,
                 extra: Optional[List[str]] = None) -> List[str]:
    """final, baseline, the basin ledger, then the best `top_k` by local score (no duplicates, only configs we have;
    configs that scored 0 - no usable calls - say nothing about the answer set and would only tie in the Spearman)."""
    ranked = sorted((k for k, c in configs.items() if c.get("local")),
                    key=lambda k: configs[k]["local"], reverse=True)
    out = []
    for k in [final, baseline] + list(extra or []) + ranked[:top_k]:
        if k and k in configs and k not in out:
            out.append(k)
    return out


def summarize(local: Dict[str, float], true: Dict[str, float], final: Optional[str],
              baseline: Optional[str]) -> dict:
    """Search/answer-set judgement over the configs that got a true score."""
    keys = [k for k in true if true[k] is not None and local.get(k) is not None]
    out = {"final_true": true.get(final), "final_local": local.get(final), "baseline_true": true.get(baseline),
           "n": len(keys), "spearman": None, "best_true": None, "best_true_config": None, "pick_gap": None}
    if out["final_true"] is not None and out["final_local"] is not None:
        out["local_minus_true"] = round(out["final_local"] - out["final_true"], 4)
    if out["final_true"] is not None and out["baseline_true"] is not None:
        out["gain_over_baseline"] = round(out["final_true"] - out["baseline_true"], 4)
    if keys:
        best = max(keys, key=lambda k: true[k])
        out.update(best_true=true[best], best_true_config=best)
        if out["final_true"] is not None:
            out["pick_gap"] = round(true[best] - out["final_true"], 4)
    if len(keys) >= 3:
        from scipy.stats import spearmanr
        rho = spearmanr([local[k] for k in keys], [true[k] for k in keys]).correlation
        out["spearman"] = None if rho != rho else round(float(rho), 4)            # NaN when a side is constant
    return out


def history(path: Path) -> Tuple[List[dict], set, set]:
    """Round lines (oldest first) and the run ids already reviewed / ranked."""
    rounds, reviewed, ranked = [], set(), set()
    if not path.exists():
        return rounds, reviewed, ranked
    for line in path.read_text().splitlines():
        try:
            row = json.loads(line)
        except json.JSONDecodeError:
            continue
        t = row.get("type")
        if t == "round" and row.get("run_id"):
            rounds.append(row)
        elif t == "review":
            reviewed.add(row.get("run_id"))
        elif t == "rank":
            ranked.add(row.get("run_id"))
    return rounds, reviewed, ranked


# ---- the reviewer -----------------------------------------------------------------------------------------------------
class Reviewer:
    def __init__(self, board=None, status_dir: Path = STATUS_DIR, platform=None, workers: Optional[int] = None,
                 top_k: int = 8,
                 poll_s: float = 300.0):
        self.board, self.dir, self.workers, self.top_k, self.poll_s = board, status_dir, workers, top_k, poll_s
        self._platform = platform
        self._scorer = None
        self._stop = threading.Event()

    # ---- background loop (live) ----------------------------------------------------------------------------------
    def start(self) -> "Reviewer":
        threading.Thread(target=self._loop, name="review", daemon=True).start()
        return self

    def stop(self) -> None:
        self._stop.set()

    def _loop(self) -> None:
        while not self._stop.is_set():
            try:
                self.check()
            except Exception:  # noqa: BLE001 - reviewing must never take the miner down
                log.exception("review pass failed")
            self._stop.wait(self.poll_s)

    def check(self) -> None:
        """One pass: review every run whose truth is out, rank every reviewed run whose round is finalized."""
        rounds, reviewed, ranked = history(self.dir / "rounds.jsonl")
        now = time.time()
        for row in rounds:
            rid, age = row["run_id"], now - (row.get("ended") or now)
            if not (self.dir / "runs" / rid / "manifest.json").exists():
                continue
            if rid not in reviewed and age < GIVE_UP_REVIEW_S:
                if self.review(rid):
                    reviewed.add(rid)
            if rid in reviewed and rid not in ranked and age < GIVE_UP_RANK_S:
                self.rank(rid)

    # ---- one run ---------------------------------------------------------------------------------------------------
    def platform(self):
        if self._platform is None:
            from common.platform import Platform
            self._platform = Platform()
        return self._platform

    def _truth(self, m: dict) -> Optional[str]:
        """The round's released files on disk (downloaded once revealed), or None while not released."""
        from backtest.data import fetch, round_dir
        d = round_dir(m["tag"])
        if os.path.exists(os.path.join(d, "truth.vcf.gz")) and os.path.exists(os.path.join(d, "mutations.vcf.gz")):
            return d
        if "reveal" not in self.platform().verification(m["round_id"]):
            return None
        return fetch(m["round_id"], self.platform(), with_bam=not os.path.exists(os.path.join(d, "input.bam")))

    def review(self, run_id: str) -> Optional[dict]:
        d = self.dir / "runs" / run_id
        m = json.loads((d / "manifest.json").read_text())
        if self._truth(m) is None:
            return None
        log.info("review %s: truth released, reviewing", run_id)
        from answer_set import OracleProvider
        from backtest.data import load
        rnd = load(m["tag"], reference=_reference(m))
        oracle = OracleProvider().build(rnd)
        row = {"type": "review", "run_id": run_id, "round_id": m["round_id"], "tag": m["tag"], "at": time.time()}

        ours = _load_answer_set(d / "answer_set.json.gz")
        if ours is not None:
            row["answer_set"] = answer_set_accuracy(ours["germline"], ours["synthetic"], oracle.germline,
                                                    [c.variant for c in oracle.synthetic], ours["n_targets"],
                                                    oracle.n_targets)
            try:
                row["answer_set"]["as_scored"] = self._answer_set_as_scored(rnd, oracle, ours, d)
            except Exception as e:  # noqa: BLE001 - the strict numbers above still stand
                log.warning("review %s: scoring our answer set failed: %s", run_id, e)
        configs = m.get("configs") or {}
        basin_keys = [b.get(field) for b in (m.get("basins") or []) for field in ("first_key", "best_key")]
        keys = pick_configs(configs, m.get("final_config"), m.get("baseline_config"), self.top_k, basin_keys)
        true, gate, rerun = self._true_scores(rnd, oracle, keys, configs, m)
        local = {k: configs[k].get("local") for k in keys}
        basin_set = {k for k in basin_keys if k}
        row["configs"] = [{"config": k, "local": local.get(k), "true": true.get(k), "gate_passed": gate.get(k),
                           "role": "final" if k == m.get("final_config") else
                           "baseline" if k == m.get("baseline_config") else
                           "basin" if k in basin_set else "top"} for k in keys]
        row["search"] = summarize(local, true, m.get("final_config"), m.get("baseline_config"))
        row["rerun"] = rerun
        self._append(row)
        self._grow_portfolio(rnd, m)
        log.info("review %s: final true %s, pick gap %s, spearman %s", run_id, row["search"].get("final_true"),
                 row["search"].get("pick_gap"), row["search"].get("spearman"))
        return row

    def rank(self, run_id: str) -> Optional[dict]:
        """Once the round is finalized: where the final config's true score lands among the real submissions."""
        from backtest.data import field_scores, round_dir
        from backtest.kpi import rank_in_field
        m = json.loads((self.dir / "runs" / run_id / "manifest.json").read_text())
        rev = self._review_row(run_id)
        true = ((rev or {}).get("search") or {}).get("final_true")
        if true is None:
            return None
        scores = os.path.join(round_dir(m["tag"]), "scores.json")
        # fetch() writes scores.json when the round is downloaded, often before the field has filled.
        # A replay then ranked that early file (2 of 4) and never looked again. Always take the live
        # field when it is at least as large as the file on disk.
        try:
            live = self.platform().scores(m["round_id"])
        except Exception:  # noqa: BLE001 - a platform blip keeps the file already on disk
            live = None
        live_n = 0
        if isinstance(live, dict):
            live_n = len({r.get("submission_id") for r in (live.get("scores") or [])
                          if r.get("combined_final") is not None and r.get("submission_id")})
        disk_n = len(field_scores(m["tag"])) if os.path.exists(scores) else 0
        if live_n and live_n >= disk_n:
            os.makedirs(os.path.dirname(scores), exist_ok=True)
            with open(scores + ".part", "w") as f:
                json.dump(live, f)
            os.replace(scores + ".part", scores)
        elif not os.path.exists(scores):
            return None
        field = field_scores(m["tag"])
        # a replay was never submitted, so it joins the field; a live run is already one of its rows
        size = len(field) + (1 if m.get("kind") == "replay" else 0)
        row = {"type": "rank", "run_id": run_id, "round_id": m["round_id"], "tag": m["tag"], "at": time.time(),
               "true_score": true, "rank": rank_in_field(true, field), "field": size,
               "field_best": max(field) if field else None,
               "field_median": round(statistics.median(field), 4) if field else None,
               "gap_to_best": round(max(field) - true, 4) if field else None,
               "vs_median": round(true - statistics.median(field), 4) if field else None}
        self._append(row)
        log.info("rank %s: %s of %s (field best %s)", run_id, row["rank"], size, row["field_best"])
        return row

    # ---- helpers ---------------------------------------------------------------------------------------------------
    def _review_row(self, run_id: str) -> Optional[dict]:
        path = self.dir / "rounds.jsonl"
        found = None
        for line in path.read_text().splitlines():
            try:
                row = json.loads(line)
            except json.JSONDecodeError:
                continue
            if row.get("type") == "review" and row.get("run_id") == run_id:
                found = row
        return found

    def _grow_portfolio(self, rnd: Round, m: dict) -> None:
        """The run's diverse near-best configs join our part of the opening's portfolio (search.portfolio)."""
        try:
            from common.params import ParamRanges
            from search import portfolio
            from search.space import Space
            try:
                r = ParamRanges.from_platform()
            except Exception:  # noqa: BLE001 - offline: the checked-in snapshot
                r = ParamRanges.from_file(str(Path(__file__).resolve().parents[1] / "tests" / "fixtures" /
                                              "parameter_ranges.json"))
            profile = portfolio.bam_profile(rnd.bam, rnd.chrom, rnd.start, rnd.end) if rnd.bam else None
            n = portfolio.add_run(Space(r, "gatk"), m.get("configs") or {}, m["tag"], profile)
            log.info("review %s: %d configs added to the portfolio (profile %s)", m["tag"], n, profile)
        except Exception as e:  # noqa: BLE001 - the portfolio must never lose a review
            log.warning("review: portfolio update failed: %s", e)

    def _answer_set_as_scored(self, rnd: Round, oracle, ours: dict, run_dir: Path) -> dict:
        import evaluator
        with open(rnd.reference + ".fai") as f:
            contigs = {c: int(n) for c, n, *_ in (line.split("\t") for line in f) if c == rnd.chrom}
        vcf = write_answer_set_vcf(ours["germline"], ours["synthetic"], contigs, run_dir / "answer_set.vcf.gz")
        self._scorer = self._scorer or evaluator.scorer()
        res = self._scorer.score(CallSet(rnd.round_id, "answer-set", str(vcf), _sha(vcf), 0.0), oracle)
        return answer_set_as_scored(res)

    def _append(self, row: dict) -> None:
        if self.board is not None:
            self.board.append(row)
            return
        with open(self.dir / "rounds.jsonl", "a") as f:
            f.write(json.dumps(row) + "\n")

    def _true_scores(self, rnd: Round, oracle, keys: List[str], configs: Dict[str, dict], m: dict):
        """True score per config: its VCF from this machine if there, else GATK runs it here (exact)."""
        import evaluator
        callsets, missing = {}, []
        work = evaluator.callset_dir(rnd.tag)
        for k in keys:
            vcf = work / f"{k}_{FULL}.vcf.gz"
            if vcf.exists():
                callsets[k] = CallSet(rnd.round_id, k, str(vcf), _sha(vcf), 0.0)
            elif configs[k].get("config"):
                missing.append(k)
        if missing:
            callsets.update(self._run(rnd, {k: Config.from_dict(configs[k]["config"]) for k in missing}))
        self._scorer = self._scorer or evaluator.scorer()
        true, gate = {}, {}
        for k, cs in callsets.items():
            try:
                res = self._scorer.score(cs, oracle)
                true[k], gate[k] = res.score, res.gate_passed
            except Exception as e:  # noqa: BLE001 - one bad callset must not lose the review
                log.warning("review: scoring %s failed: %s", k, e)
        return true, gate, len(missing)

    def _run(self, rnd: Round, configs: Dict[str, Config]) -> Dict[str, CallSet]:
        import evaluator
        # solo and in parallel: one config per worker, no batching. A batch runs its variants one after another in
        # one JVM thread (~2.5-3x a solo run for assembly-level variants), which made a 9-config review take ~10 min.
        workers = min(self.workers or len(configs), len(configs), default_workers())
        executor = evaluator.local_executor(workers, batching=False)
        log.info("review %s: running %d configs here, solo on %d workers (their VCFs are not on this machine)",
                 rnd.tag, len(configs), workers)
        jobs = {f"review-{rnd.tag}-{i}": k for i, k in enumerate(configs)}
        executor.submit([EvalJob(j, rnd, configs[k], FULL) for j, k in jobs.items()])
        out, deadline = {}, time.time() + 3600
        while len(out) < len(jobs) and time.time() < deadline:
            for cs in executor.poll(5.0):
                if cs.job_id in jobs and cs.ok and cs.vcf:
                    out[jobs[cs.job_id]] = cs
                elif cs.job_id in jobs:
                    log.warning("review: %s failed: %s", jobs[cs.job_id], cs.error)
                    out.setdefault(jobs[cs.job_id], None)
        executor._pool.shutdown(wait=False)
        return {k: cs for k, cs in out.items() if cs is not None}


def _reference(m: dict) -> str:
    chrom = Round.parse_region(m["region"])[0]
    return SETTINGS.reference_fasta(chrom)


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _load_answer_set(path: Path) -> Optional[dict]:
    if not path.exists():
        return None
    with gzip.open(path, "rt") as f:
        a = json.load(f)
    return {"n_targets": a.get("n_targets"),
            "germline": [Variant(c, int(p), r, al, g) for c, p, r, al, g in a["germline"]],
            "synthetic": [(Variant(c, int(p), r, al, g), float(pp)) for c, p, r, al, g, pp in a["synthetic"]]}


def find_run(ref: str, status_dir: Path = STATUS_DIR) -> str:
    """A run id, or a round tag (its newest run)."""
    if (status_dir / "runs" / ref).exists():
        return ref
    runs = sorted(p.name for p in (status_dir / "runs").glob(f"{ref}-*"))
    if not runs:
        raise SystemExit(f"no stored run for {ref!r} in {status_dir / 'runs'}")
    return runs[-1]
