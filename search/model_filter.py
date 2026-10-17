"""A per-class score model of the round, to choose which GA children are worth a run (2026-10-17).

test/bench/model_feasibility.py on 84,975 configs of 20261016T130600Z: gradient boosting per class (tp, fp, fn of the six
classes and the region FP, the score rebuilt with the validator's own formula) ranks near-top configs with Spearman 0.74
(trained on one run, tested on later runs) to 0.92 (random split); a direct score model 0.28-0.81. It cannot name the last
call (error 0.02-0.04, one call ~0.005), so it filters, it does not pick: with climb.model_filter = k the GA breeds k
times as many children and sends the 1/k the model ranks highest. 1 (the default) changes nothing.

    fit     once 2,000 configs are scored, then every 4,000 more; features: each searched setting's ladder position
    rank    one batch prediction per breeding call
"""
from __future__ import annotations

import logging
from typing import List, Optional, Sequence

from common.types import CLASSES, ClassCounts

log = logging.getLogger(__name__)
FIRST_FIT, REFIT_EVERY = 2000, 4000
TARGETS = [(c, x) for c in CLASSES for x in ("tp", "fp", "fn")]


class ScoreModel:
    def __init__(self, space, ladders):
        self.names = list(space.names)
        self.ladders = ladders
        self.models = None
        self.fitted_at = 0
        self.n_targets = None

    def features(self, params) -> List[float]:
        out = []
        for n in self.names:
            try:
                out.append(float(self.ladders.index(n, params.get(n))))
            except Exception:  # noqa: BLE001
                out.append(-1.0)
        return out

    def maybe_fit(self, results: Sequence) -> bool:
        """results: (config, ScoreResult) of every scored config. Refits when enough new ones came in."""
        n = len(results)
        if n < FIRST_FIT or (self.models is not None and n - self.fitted_at < REFIT_EVERY):
            return False
        try:
            import numpy as np
            from sklearn.ensemble import HistGradientBoostingRegressor
        except ImportError:  # pragma: no cover - numpy and scikit-learn are in requirements.txt
            return False
        X = np.array([self.features(c.params) for c, _ in results])
        Y = np.array([[getattr(r.counts[c], x) for c, x in TARGETS] + [r.region_fp] for _, r in results], dtype=float)
        try:
            self.models = [HistGradientBoostingRegressor(max_iter=150, learning_rate=0.1, random_state=0).fit(X, Y[:, j])
                           for j in range(Y.shape[1])]
        except Exception:  # noqa: BLE001 - a failed fit keeps the previous model
            log.exception("model filter: fit failed")
            return False
        self.n_targets = results[-1][1].n_targets
        self.fitted_at = n
        return True

    def predict(self, configs: Sequence) -> Optional[List[float]]:
        if self.models is None or not configs:
            return None
        import numpy as np
        from common.scoring import v2_score
        X = np.array([self.features(c.params) for c in configs])
        P = np.column_stack([m.predict(X) for m in self.models])
        out = []
        for row in P:
            counts = {c: ClassCounts(tp=max(0, int(round(row[i * 3]))), fp=max(0, int(round(row[i * 3 + 1]))),
                                     fn=max(0, int(round(row[i * 3 + 2])))) for i, c in enumerate(CLASSES)}
            out.append(v2_score(counts, max(0, int(round(row[-1]))), self.n_targets or 1, None)[1] or 0.0)
        return out


def model_of(ctx) -> ScoreModel:
    m = getattr(ctx, "_score_model", None)
    if m is None:
        m = ScoreModel(ctx.space, ctx.ladders)
        ctx._score_model = m
    return m


def scored_results(ctx) -> list:
    out = []
    client = ctx.client
    for (key, fid), r in client.latest.items():
        if r.score is None or r.gate_passed is False or not getattr(r, "counts", None):
            continue
        c = client.configs.get(key)
        if c is not None:
            out.append((c, r))
    return out


def keep_best(ctx, children: List, n: int) -> List:
    """The n children the model ranks highest (all of them while there is no model yet)."""
    if len(children) <= n:
        return children
    m = model_of(ctx)
    if m.maybe_fit(scored_results(ctx)):
        ctx.say("model filter: per-class model fitted on %d scored configs", m.fitted_at)
    pred = m.predict(children)
    if pred is None:
        return children[:n]
    order = sorted(range(len(children)), key=lambda i: -pred[i])
    return [children[i] for i in order[:n]]
