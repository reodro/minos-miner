"""RF_UCB's model: a random forest of the score over the whole space, fitted on every scored config of the round.

Encoding: each searched setting is its position on its ladder (search/v9/steps.py), 0..1. A forest splits on
thresholds, so a setting whose score jumps from one value to the next (min_mapping_quality_score, the BAM groups) is
fitted as easily as a smooth one.

Each tree is grown on its own random rf.tree_sample share of ALL the round's configs (no region preferred): at 40,000
configs a fit takes 3.4 s on 14 threads instead of 14 s on 8 with every config in every tree, and the trees differ
more, which is the spread UCB reads. Every config is encoded once and cached.

Weights: a config's weight is 1/sqrt(the configs within rf.density_radius of it, itself included), so a cluster (the
portfolio's members and their neighbours, a hill sampled densely by local moves) does not pull the whole forest toward
itself.

UCB = mean + kappa_eff x spread, where spread is the standard deviation across the trees. kappa_eff is set at every
ranking so that, among the 5% best predicted, the bonus is rf.kappa times the spread of their predicted means: a fixed
kappa against a tree spread of 0.058 and moves worth 0.01 ranked proposals by noise (legacy FINDINGS 2026-09-10).
"""
from __future__ import annotations

import logging
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from common.types import Config

log = logging.getLogger(__name__)


NOV_CAP = 3.0          # the novelty bonus never exceeds this many typical (median-candidate) bonuses


class Forest:
    def __init__(self, space, ladders, trees: int = 200, min_leaf: int = 2, seed: int = 0, jobs: int = 8,
                 sample: float = 1.0):
        self.space, self.ladders = space, ladders
        self.names = [n for n in space.names if ladders.movable(n)]
        self.trees, self.min_leaf, self.seed, self.jobs = trees, min_leaf, seed, jobs
        self.sample = sample                        # each tree's bootstrap: this share of ALL configs (1.0: all)
        self.model = None
        self.n_fit = 0
        self._rows: Dict[str, List[float]] = {}     # config key -> its encoding (configs are encoded once)
        self._tree = None                           # KD-tree of the configs of the last fit: distance to the data
        self.bad_model = None                       # P(a config crashes or scores zero), when enough are known

    # ---- encoding -------------------------------------------------------------------------------------------------
    def row(self, c: Config, keep: bool = False) -> List[float]:
        """A config's encoding. Only the scored configs are kept (keep=True, at fitting): candidates are encoded once
        per ranking and dropped, or the cache would grow by every candidate of every refit (millions a round)."""
        r = self._rows.get(c.key)
        if r is None:
            r = []
            for n in self.names:
                k = len(self.ladders.rungs[n]) - 1
                r.append(self.ladders.index(n, c.params.get(n)) / k if k > 0 else 0.0)
            if keep:
                self._rows[c.key] = r
        return r

    def matrix(self, configs: Sequence[Config], keep: bool = False) -> np.ndarray:
        return np.ascontiguousarray([self.row(c, keep) for c in configs], dtype=np.float32)

    # ---- fitting --------------------------------------------------------------------------------------------------
    @staticmethod
    def density_weights(X: np.ndarray, radius: float) -> np.ndarray:
        """1/sqrt(configs within `radius` of each, by mean absolute difference per setting, itself included)."""
        if radius <= 0 or len(X) < 2:
            return np.ones(len(X))
        from sklearn.neighbors import KDTree
        counts = KDTree(X, metric="manhattan").query_radius(X, r=radius * X.shape[1], count_only=True)
        return 1.0 / np.sqrt(np.maximum(counts, 1.0))

    def fit(self, configs: Sequence[Config], scores: Sequence[float], radius: float) -> None:
        from sklearn.ensemble import RandomForestRegressor
        if len(configs) < 8:
            self.model = None
            return
        X = self.matrix(configs, keep=True)
        y = np.asarray(scores, dtype=np.float64)
        w = self.density_weights(X, radius)
        sample = self.sample if 0 < self.sample < 1 else None
        model = RandomForestRegressor(n_estimators=self.trees, min_samples_leaf=self.min_leaf, n_jobs=self.jobs,
                                      max_samples=sample, random_state=self.seed + self.n_fit)
        model.fit(X, y, sample_weight=w)
        from sklearn.neighbors import KDTree
        tree = KDTree(X, metric="manhattan")
        self.model, self._tree, self.n_fit = model, tree, self.n_fit + 1

    def fit_bad(self, good: Sequence[Config], bad: Sequence[Config], min_bad: int = 20) -> None:
        """A classifier of crashes and zero scores, apart from the score forest: taught as a floor inside the score
        forest (second fleet replay) they pulled every prediction down by 1.4-2.1 and lowered its Spearman (sweep 0.91
        -> 0.76)."""
        if len(bad) < min_bad or len(good) < min_bad:
            self.bad_model = None
            return
        from sklearn.ensemble import RandomForestClassifier
        X = self.matrix(list(good) + list(bad), keep=True)
        y = np.array([0] * len(good) + [1] * len(bad))
        sample = self.sample if 0 < self.sample < 1 else None
        m = RandomForestClassifier(n_estimators=max(20, self.trees // 2), min_samples_leaf=self.min_leaf, n_jobs=self.jobs,
                                   max_samples=sample, class_weight="balanced", random_state=self.seed + 3)
        m.fit(X, y)
        self.bad_model = m

    def p_bad(self, configs: Sequence[Config]) -> np.ndarray:
        if self.bad_model is None or not configs:
            return np.zeros(len(configs))
        return self.bad_model.predict_proba(self.matrix(configs))[:, 1]

    @property
    def ready(self) -> bool:
        return self.model is not None

    # ---- prediction -----------------------------------------------------------------------------------------------
    def per_tree(self, configs: Sequence[Config]) -> np.ndarray:
        """trees x configs. One float32 array through the raw trees (sklearn's predict re-validates per tree)."""
        X = self.matrix(configs)
        return np.stack([t.tree_.predict(X)[:, 0] for t in self.model.estimators_])

    def mean_sd(self, configs: Sequence[Config]) -> Tuple[np.ndarray, np.ndarray]:
        if not self.ready or not configs:
            return np.zeros(len(configs)), np.full(len(configs), np.inf)
        pt = self.per_tree(configs)
        return pt.mean(axis=0), pt.std(axis=0, ddof=1)

    def novelty(self, configs: Sequence[Config]) -> np.ndarray:
        """Each config's distance to the nearest config of the last fit: mean ladder distance per setting."""
        if self._tree is None or not configs:
            return np.zeros(len(configs))
        X = self.matrix(configs)
        d, _ = self._tree.query(X, k=1)
        return d[:, 0] / max(1, X.shape[1])

    def ucb(self, configs: Sequence[Config], kappa: float, kappa_dist: float = 0.0
            ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """(ucb, mean, spread, novelty). UCB = mean + k x spread + kd x novelty, both bonuses weighed against the spread
        of the predicted means among the 5% best predicted: the spread bonus's median there is kappa times it, a typical
        candidate's novelty bonus kappa_dist times it.

        The novelty term is there because a forest does not know where it has not looked: far from the data each
        tree carries over its nearest box's value, the trees agree, and the spread stays low exactly where an unseen
        hill could be (at the start, the portfolio's boxes would fill the empty space)."""
        mu, sd = self.mean_sd(configs)
        nov = self.novelty(configs)
        if not self.ready or not len(mu):
            return mu, mu, sd, nov
        top = np.argsort(-mu)[:max(2, len(mu) // 20)]
        mu_spread = float(np.std(mu[top]))
        sd_top = float(np.median(sd[top]))
        nov_all = float(np.median(nov))
        k = kappa * mu_spread / sd_top if sd_top > 0 else 0.0
        # The novelty bonus of a typical candidate (the median over all of them, not over the top, which is mostly
        # local moves with almost no novelty: scaling by that blew the term up and ranked "most novel first") is
        # kappa_dist x the predicted spread at the top; a config twice as far from the data gets twice that.
        kd = kappa_dist * mu_spread / nov_all if nov_all > 0 else 0.0
        # Capped at NOV_CAP typical bonuses: on the first fleet replay a restart reset QGA to uniform, its draws were
        # far from everything, and the uncapped term let them take 52% of the picks with nothing to climb from.
        bonus = np.minimum(kd * nov, NOV_CAP * kappa_dist * mu_spread)
        u = mu + k * sd + bonus
        pb = self.p_bad(configs)
        if pb.any():                                # a likely crash or zero sinks: by its probability, across the pool
            u = u - pb * float(u.max() - u.min())
        return u, mu, sd, nov

    def flat(self, centre: Config, name: str, points: int, resolution: float) -> bool:
        """QLock's test: sweeping `name` alone across its ladder from `centre` moves the predicted score less than
        `resolution`, counting the trees' spread (an axis the forest barely knows is never flat)."""
        if not self.ready:
            return False
        rungs = self.ladders.rungs[name]
        idx = sorted({int(round(i * (len(rungs) - 1) / max(1, points - 1))) for i in range(points)})
        cs = [c for c in (self.space._assemble({**centre.params, name: rungs[i]}) for i in idx) if c is not None]
        if len(cs) < 2:
            return False
        mu, sd = self.mean_sd(cs)
        return float((mu + sd).max() - (mu - sd).min()) < resolution
