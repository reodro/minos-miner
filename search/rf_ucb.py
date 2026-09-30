"""RF_UCB: random-forest surrogate + upper-confidence-bound picks, streaming, family-aware. Simple first version.

Loop: keep every worker busy (the executor's idle_slots; `in_flight` configs if it cannot tell). Whenever one frees, refit the forest on every result the client
has (all strategies share them), score `n_candidates` proposals - mostly local moves around the best configs so
far, some uniform random - by mean + kappa * spread across trees, and submit the best one together with a small
family of cheap variants of it (same config except parameters the evaluator can batch almost for free), so the
executor runs them as one GATK job. Failed configs count as the worst score seen, so the model learns to avoid them.
The surrogate is scikit-learn's ExtraTreesRegressor; uncertainty is the spread of its trees' predictions.

Tuning knobs are the arguments; nothing here is fitted to real rounds yet - that is the search part's work.
"""
from __future__ import annotations

import asyncio
import random
from typing import List, Optional, Set

import numpy as np
from sklearn.ensemble import ExtraTreesRegressor

from common.params import ParamRanges
from common.types import FULL, Config
from search.base import Strategy
from search.client import EvalClient, Outcome
from search.space import Space


def rf_ucb(ranges: ParamRanges, *, tool: str = "gatk", n_init: int = 16, in_flight: int = 8, n_candidates: int = 300,
           kappa: float = 1.0, top_k: int = 8, explore_frac: float = 0.25, family: int = 4,
           family_max_cost: float = 0.15, n_trees: int = 30, max_pending: int = 256, seed: int = 0,
           tag: str = "rf_ucb") -> Strategy:
    space = Space(ranges, tool)
    rng = random.Random(seed)

    async def run(client: EvalClient) -> None:
        failed: Set[str] = set()
        proposed: Set[str] = set()
        model: List[Optional[ExtraTreesRegressor]] = [None]
        fitted_on = [0]

        def observations():
            ok = [(client.configs[k], r.score) for (k, fid), r in client.latest.items()
                  if fid == FULL and r.score is not None and k in client.configs and client.configs[k].tool == tool]
            if not ok:
                return [], []
            worst = min(s for _, s in ok)
            rows = ok + [(client.configs[k], worst) for k in failed if k in client.configs]
            return [c for c, _ in rows], [s for _, s in rows]

        def refit_if_due():
            configs, y = observations()
            if len(y) < 4:
                return
            if model[0] is None or len(y) - fitted_on[0] >= max(1, len(y) // 20):
                model[0] = ExtraTreesRegressor(n_estimators=n_trees, min_samples_leaf=3, max_features="sqrt",
                                               random_state=rng.randrange(1 << 30)).fit(
                    np.array([space.encode(c) for c in configs]), np.array(y))
                fitted_on[0] = len(y)

        def ucb_many(cs: List[Config]) -> List[float]:
            """mean + kappa * spread of the trees' predictions, for many configs at once."""
            X = np.array([space.encode(c) for c in cs])
            per_tree = np.stack([t.predict(X) for t in model[0].estimators_])
            return list(per_tree.mean(axis=0) + kappa * per_tree.std(axis=0))

        def fresh(c: Optional[Config]) -> bool:
            return c is not None and c.key not in proposed and (c.key, FULL) not in client.latest and c.key not in failed

        def propose() -> List[Config]:
            if model[0] is None:
                c = space.random(rng)
                return [c] if fresh(c) else []
            configs, y = observations()
            ranked = [c for _, c in sorted(zip(y, configs), key=lambda t: -t[0])][:top_k]
            cands = []
            for _ in range(n_candidates):
                c = space.random(rng) if (rng.random() < explore_frac or not ranked) else space.neighbour(rng.choice(ranked), rng)
                if fresh(c):
                    cands.append(c)
            if not cands:
                return []
            scores = ucb_many(cands)
            best = cands[max(range(len(cands)), key=scores.__getitem__)]
            fam = []
            if family > 0:
                near = [space.neighbour(best, rng, n_params=2, scale=0.2) for _ in range(40)]
                cheap = [c for c in near if fresh(c) and c.key != best.key
                         and client.variant_cost(best, c) <= family_max_cost]
                seen = {best.key}
                ranked_cheap = [c for _, c in sorted(zip(ucb_many(cheap), cheap), key=lambda t: -t[0])] if cheap else []
                for c in ranked_cheap:
                    if c.key not in seen:
                        seen.add(c.key); fam.append(c)
                    if len(fam) == family:
                        break
            return [best] + fam

        def note(o: Outcome):
            if not o.ok and o.config is not None and not (o.error or "").startswith("budget"):
                failed.add(o.config.key)

        def room() -> bool:
            """Propose more while workers are free (a family runs as one job), else keep `in_flight` configs
            running when the executor cannot report free workers. Capped so the queue cannot grow unbounded."""
            if len(pending) >= max_pending:
                return False
            idle = client.idle_slots()
            return (len(pending) < in_flight) if idle is None else (idle > 0 or not pending)

        pending = set()
        init = []
        while len(init) < n_init:
            c = space.random(rng)
            if fresh(c) and c.key not in {x.key for x in init}:
                init.append(c)
        proposed.update(c.key for c in init)
        pending |= set(client.submit_many(init, tag=tag))
        while True:
            while not client.exhausted and room():
                refit_if_due()
                batch = propose()
                if not batch:
                    break
                proposed.update(c.key for c in batch)
                pending |= set(client.submit_many(batch, tag=tag))
            if not pending:
                return
            done, pending = await asyncio.wait(pending, return_when=asyncio.FIRST_COMPLETED)
            for f in done:
                note(f.result())
            if client.exhausted and not pending:
                return

    return run
