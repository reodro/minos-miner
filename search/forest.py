"""A small random forest regressor in pure Python (no numpy in this project): extremely randomized trees.

Each tree is grown on a bootstrap sample; at each node it tries `max_features` random features with one random
threshold each and keeps the split that most reduces variance. predict() returns the mean and the spread across
trees - the spread is the uncertainty that UCB rewards. Good enough for a few thousand points and ~25 features;
swap for scikit-learn if the search part adds dependencies.
"""
from __future__ import annotations

import math
import random
from typing import List, Optional, Sequence, Tuple


class _Node:
    __slots__ = ("feature", "threshold", "left", "right", "value")

    def __init__(self, value: float, feature: int = -1, threshold: float = 0.0,
                 left: Optional["_Node"] = None, right: Optional["_Node"] = None):
        self.value, self.feature, self.threshold, self.left, self.right = value, feature, threshold, left, right


class ExtraTreesRegressor:
    def __init__(self, n_trees: int = 30, min_leaf: int = 3, max_depth: int = 12, max_features: Optional[int] = None,
                 seed: int = 0):
        self.n_trees, self.min_leaf, self.max_depth, self.max_features = n_trees, min_leaf, max_depth, max_features
        self.rng = random.Random(seed)
        self.trees: List[_Node] = []

    def fit(self, X: Sequence[Sequence[float]], y: Sequence[float]) -> "ExtraTreesRegressor":
        n, d = len(X), len(X[0])
        m = self.max_features or max(1, int(math.sqrt(d)) + 1)
        self.trees = []
        for _ in range(self.n_trees):
            idx = [self.rng.randrange(n) for _ in range(n)]
            self.trees.append(self._grow(X, y, idx, 0, d, m))
        return self

    def _grow(self, X, y, idx: List[int], depth: int, d: int, m: int) -> _Node:
        ys = [y[i] for i in idx]
        mean = sum(ys) / len(ys)
        if depth >= self.max_depth or len(idx) < 2 * self.min_leaf or max(ys) - min(ys) < 1e-12:
            return _Node(mean)
        best: Optional[Tuple[float, int, float, List[int], List[int]]] = None
        for f in self.rng.sample(range(d), min(m, d)):
            vals = [X[i][f] for i in idx]
            lo, hi = min(vals), max(vals)
            if hi - lo < 1e-12:
                continue
            t = self.rng.uniform(lo, hi)
            left = [i for i in idx if X[i][f] <= t]
            right = [i for i in idx if X[i][f] > t]
            if len(left) < self.min_leaf or len(right) < self.min_leaf:
                continue
            sse = _sse([y[i] for i in left]) + _sse([y[i] for i in right])
            if best is None or sse < best[0]:
                best = (sse, f, t, left, right)
        if best is None:
            return _Node(mean)
        _, f, t, left, right = best
        return _Node(mean, f, t, self._grow(X, y, left, depth + 1, d, m), self._grow(X, y, right, depth + 1, d, m))

    def predict(self, x: Sequence[float]) -> Tuple[float, float]:
        """(mean, standard deviation) of the trees' predictions."""
        preds = []
        for node in self.trees:
            while node.left is not None:
                node = node.left if x[node.feature] <= node.threshold else node.right
            preds.append(node.value)
        mu = sum(preds) / len(preds)
        return mu, math.sqrt(sum((p - mu) ** 2 for p in preds) / len(preds))


def _sse(v: List[float]) -> float:
    m = sum(v) / len(v)
    return sum((a - m) ** 2 for a in v)
