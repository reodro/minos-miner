"""The (core, FP) forest: ranks candidate configs by predicted score gain per cost when there are more candidates
than free workers. Owner: search.

Two ExtraTrees, one on core and one on the FP rate, combined through the real v2 formula; the spread between trees
is the uncertainty bonus. It orders scouts and pair moves; it never accepts a move - that is always an exact
evaluation (on a stepped landscape the model cannot resolve the 0.005-0.05 differences that decide a step).
"""
from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

import numpy as np

from common.scoring import score_from_parts
from common.types import Config
from search.space import Space


class Ranker:
    def __init__(self, space: Space, n_trees: int = 40, min_points: int = 30, seed: int = 0):
        self.space, self.n_trees, self.min_points, self.seed = space, n_trees, min_points, seed
        self._x: List[List[float]] = []
        self._core: List[float] = []
        self._fp: List[float] = []
        self._models = None
        self._fitted_on = 0
        self._worst = (1.0, 0.0)

    def add(self, config: Config, core: Optional[float], fp_per_target: Optional[float], failed: bool = False) -> None:
        if failed or core is None or fp_per_target is None:
            core, fp_per_target = self._worst                 # a crash / no calls: as bad as the worst seen
        else:
            self._worst = (min(self._worst[0], core), max(self._worst[1], fp_per_target))
        self._x.append(self.space.encode(config))
        self._core.append(float(core))
        self._fp.append(float(fp_per_target))

    @property
    def ready(self) -> bool:
        return self._models is not None

    def stale(self, every: int = 100) -> bool:
        return len(self._x) >= self.min_points and len(self._x) - self._fitted_on >= every

    def fit(self) -> None:
        """Refit on everything seen (a few seconds at 10k points; call it from a thread)."""
        from sklearn.ensemble import ExtraTreesRegressor
        n = min(len(self._x), len(self._core), len(self._fp))      # add() may append while this thread runs
        if n < self.min_points:
            return
        x = np.array(self._x[:n])
        core = ExtraTreesRegressor(n_estimators=self.n_trees, min_samples_leaf=2, max_features=0.5,
                                   random_state=self.seed, n_jobs=1).fit(x, np.array(self._core[:n]))
        fp = ExtraTreesRegressor(n_estimators=self.n_trees, min_samples_leaf=2, max_features=0.5,
                                 random_state=self.seed + 1, n_jobs=1).fit(x, np.array(self._fp[:n]))
        self._models, self._fitted_on = (core, fp), n

    def predict(self, configs: Sequence[Config]) -> Tuple[np.ndarray, np.ndarray]:
        """(mean predicted score, spread between trees) per config."""
        core_m, fp_m = self._models
        x = np.array([self.space.encode(c) for c in configs])
        core_t = np.stack([t.predict(x) for t in core_m.estimators_])
        fp_t = np.stack([t.predict(x) for t in fp_m.estimators_])
        per_tree = np.vectorize(score_from_parts)(core_t, fp_t)
        return per_tree.mean(axis=0), per_tree.std(axis=0)

    def rank(self, configs: Sequence[Config], costs: Optional[Sequence[float]] = None, best: float = 0.0,
             kappa: float = 0.5) -> List[int]:
        """Indices of `configs`, best first, by (predicted gain over `best` + kappa x spread) per cost. Without a
        fitted model the order is unchanged."""
        if not configs or not self.ready:
            return list(range(len(configs)))
        mean, spread = self.predict(configs)
        gain = mean - best + kappa * spread
        cost = np.array(costs if costs is not None else [1.0] * len(configs), float)
        util = np.where(gain > 0, gain / np.maximum(cost, 0.02), gain)
        return [int(i) for i in np.argsort(-util)]
