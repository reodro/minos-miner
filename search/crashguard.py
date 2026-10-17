"""Crash guard: learn which configs make GATK crash in this round, and stop sending them.

044200Z (2026-10-16, v4): 4.9% of runs failed, every one with GATK's "Graph must have ref source and sink vertices".
By source: GA children 56%, riders lost with a crashed head 18%, centre jumps 15%. The crash depends on combinations
(a short reference window: low assembly_region_padding with small or large region sizes, the active threshold), so a
per-setting edge (Climber.beyond_crash) only catches the jumps. This guard fits a shallow decision tree (scikit-learn)
on the round's crashed configs against its completed ones, over each searched setting's ladder position, and refits as
outcomes arrive. Its probabilities are the raw crash rates of the tree's leaves (no class weighting: with a 5% base
rate, balanced weights would put a leaf that crashes 7% of the time at 0.6). Context.fresh asks it before anything is sent: a config the tree gives SKIP_P or more is not sent; a
climber sends a head the tree gives SOLO_P or more without riders, so a crash does not take them down too.

Audits (2026-10-17): the guard withheld 5,224 configs on 20261016T130600Z while 349 runs crashed, and a withheld config
was never checked. AUDIT of them are now sent anyway; a leaf whose audited configs crash less than AUDIT_RATE of the time
(after AUDIT_MIN audits) stops blocking until the next fit, and every audited outcome is training data for that fit.
A withheld config is never marked proposed, so the search can send it again once its leaf is off.
"""
from __future__ import annotations

import collections
import logging
import random
from typing import Dict, List, Optional

log = logging.getLogger(__name__)

MIN_CRASHES = 40            # no verdicts before this many crashes
REFIT_EVERY = 150           # outcomes between refits
SKIP_P = 0.6                # crash probability at which a config is not sent
SOLO_P = 0.3                # at which a head goes without riders
KEEP_OK = 6000              # completed configs kept for fitting (most recent)
AUDIT = 0.1                 # share of the configs it would withhold that are sent anyway
AUDIT_MIN, AUDIT_RATE = 5, 0.3   # a leaf whose audits crash less often than this (after AUDIT_MIN) stops blocking


class CrashGuard:
    def __init__(self, space, ladders):
        self.names = list(space.names)
        self.ladders = ladders
        self.crashed: List[List[float]] = []
        self.ok: collections.deque = collections.deque(maxlen=KEEP_OK)
        self.model = None
        self.pending = 0
        self.skipped: set = set()
        self._cache: Dict[str, float] = {}
        self.rng = random.Random(5)
        self.audits: Dict[str, int] = {}                # audited config key -> its leaf
        self.leaf_audit: Dict[int, List[int]] = {}      # leaf -> [crashed, audited]
        self.off: set = set()                           # leaves that stopped blocking (this fit)
        self.n_audits = self.n_audit_crashed = self.n_off = 0

    def features(self, params) -> List[float]:
        out = []
        for n in self.names:
            try:
                out.append(float(self.ladders.index(n, params.get(n))))
            except Exception:  # noqa: BLE001 - a value the ladder cannot place
                out.append(-1.0)
        return out

    def note(self, params, crashed: bool, key: Optional[str] = None) -> None:
        """One outcome: a crash, or a run that came back (scored or gate-failed)."""
        (self.crashed if crashed else self.ok).append(self.features(params))
        leaf = self.audits.pop(key, None) if key is not None else None
        if leaf is not None:
            st = self.leaf_audit.setdefault(leaf, [0, 0])
            st[0] += int(crashed)
            st[1] += 1
            self.n_audit_crashed += int(crashed)
            if leaf not in self.off and st[1] >= AUDIT_MIN and st[0] / st[1] < AUDIT_RATE:
                self.off.add(leaf)
                self.n_off += 1
                log.info("crash guard: a rule stops blocking: %d of %d audited configs crashed", st[0], st[1])
        self.pending += 1
        if self.pending >= REFIT_EVERY and len(self.crashed) >= MIN_CRASHES and len(self.ok) >= MIN_CRASHES:
            self.fit()

    def fit(self) -> None:
        try:
            from sklearn.tree import DecisionTreeClassifier
        except ImportError:  # pragma: no cover - scikit-learn is in requirements.txt
            return
        x = list(self.crashed) + list(self.ok)
        y = [1] * len(self.crashed) + [0] * len(self.ok)
        model = DecisionTreeClassifier(max_depth=5, min_samples_leaf=10, random_state=0)
        try:
            model.fit(x, y)
        except Exception:  # noqa: BLE001 - a fit that fails leaves the previous verdicts in place
            log.exception("crash guard: fit failed")
            return
        self.model = model
        self.pending = 0
        self._cache.clear()
        self.off.clear()                                # a new tree: its leaves start blocking again
        self.leaf_audit.clear()
        self.audits.clear()

    def p_crash(self, config) -> float:
        if self.model is None:
            return 0.0
        key = getattr(config, "key", None)
        if key is not None and key in self._cache:
            return self._cache[key]
        p = self._leaf_rate(self.features(config.params))
        if key is not None:
            self._cache[key] = p
        return p

    def _leaf(self, x: List[float]) -> int:
        t = self.model.tree_
        node = 0
        while t.children_left[node] != -1:
            node = t.children_left[node] if x[t.feature[node]] <= t.threshold[node] else t.children_right[node]
        return node

    def _leaf_rate(self, x: List[float]) -> float:
        """The crash rate of x's leaf: predict_proba's answer, walked in Python (one sample at a time through
        predict_proba costs ~0.2 ms of input checks; fresh() asks for every GA child)."""
        t = self.model.tree_
        node = 0
        while t.children_left[node] != -1:
            node = t.children_left[node] if x[t.feature[node]] <= t.threshold[node] else t.children_right[node]
        v = t.value[node][0]
        total = float(v[0] + v[1])
        return float(v[1]) / total if total > 0 else 0.0

    def blocks(self, config) -> bool:
        if self.model is None or self.p_crash(config) < SKIP_P:
            return False
        key = getattr(config, "key", None)
        leaf = self._leaf(self.features(config.params))
        if leaf in self.off:
            return False
        if key not in self.skipped and self.rng.random() < AUDIT:
            self.audits[key] = leaf                     # sent anyway: its outcome checks the rule
            self.n_audits += 1
            return False
        self.skipped.add(key)
        return True

    def risky(self, config) -> bool:
        return self.model is not None and self.p_crash(config) >= SOLO_P

    def summary(self) -> Optional[str]:
        if not self.crashed and not self.skipped:
            return None
        audited = self.n_audits - len(self.audits)
        return ("crash guard: %d crashes and %d completed runs learnt, %d configs not sent; %d audited (%d back, %d "
                "crashed), %d rules stopped blocking" % (len(self.crashed), len(self.ok), len(self.skipped),
                                                         self.n_audits, audited, self.n_audit_crashed, self.n_off))
