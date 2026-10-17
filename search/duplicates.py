"""Duplicate callsets: park a config predicted to repeat a callset, audit the prediction, and send it the moment the
prediction is contradicted or could matter. Nothing is dropped on a prediction alone.

044200Z v4 (2026-10-16): 9,688 of 42,025 scored configs (141 of 752 core-hours) returned a callset the same climber had
already returned, mostly a one-setting change of max_num_haplotypes_in_population, max_alternate_alleles,
indel_heterozygosity or contamination_fraction_to_filter. A change that tied on one centre gave a different score on
another for 882 changes, so a remembered tie is a prediction, not a fact.

    rule      (setting, a, b): a change between two values of one setting. It is established once it has tied
              dups.min_ties times between scored configs that differ only there, and never given a different score.
    park      a config one established change away from a config already sent (its twin, scored or running) is
              predicted to return the twin's callset. It is kept under its rule, not dropped.
    audit     dups.audit of the predicted configs are sent anyway. Any scored pair that contradicts a rule (an audit,
              or any other pair) breaks the rule for the rest of the round, and every config parked under it is sent.
    near      a config whose twin scores within dups.near of the island's best (when it is parked, or when the twin
              comes back) is sent at once: only those could have moved the round's best.
    exact     a config already scored or running is never sent again (Context.fresh); that is not a prediction.
    never     nothing is parked while the round's pick is certified (local_opt.certify_round).

Replay of v4's own send order through this policy: search/DESIGN.md (savings, and configs that stayed parked).
"""
from __future__ import annotations

import collections
import random
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from common.tuning import tune

PARKED = "parked: predicted duplicate"       # the error text of a parked config's outcome
TIE = 1e-9


class DuplicateGuard:
    def __init__(self, names: Sequence[str], rng: Optional[random.Random] = None):
        self.names = list(names)
        self.rng = rng or random.Random(0)
        n = len(self.names)
        self.vec: Dict[str, tuple] = {}                                              # key -> setting values
        self.sent: List[Dict[int, Dict[object, str]]] = [collections.defaultdict(dict) for _ in range(n)]
        self.scored: List[Dict[int, Dict[object, str]]] = [collections.defaultdict(dict) for _ in range(n)]
        self.score_of: Dict[str, float] = {}
        self.rules: Dict[tuple, List[int]] = collections.defaultdict(lambda: [0, 0])  # rule -> [tied, differed]
        self.broken: set = set()
        self.parked: Dict[tuple, Dict[str, tuple]] = collections.defaultdict(dict)  # rule -> key -> (config, tag)
        self.parked_rule: Dict[str, tuple] = {}
        self.waiting_on: Dict[str, List[str]] = collections.defaultdict(list)       # twin key -> parked keys
        self.n_parked = self.n_audits = self.n_broken = self.n_released = self.n_near = 0

    # ---- indexing ----------------------------------------------------------------------------------------------------
    @staticmethod
    def rule(name: str, a, b) -> tuple:
        lo, hi = sorted((a, b), key=repr)
        return (name, lo, hi)

    def _vec(self, params) -> tuple:
        return tuple(params.get(n) for n in self.names)

    @staticmethod
    def _proj(v: tuple, d: int) -> int:
        return hash(v[:d] + v[d + 1:])

    def _same_but(self, v: tuple, key: str, d: int) -> bool:
        """`key` differs from `v` in setting d only (a hash match is checked against the values themselves)."""
        u = self.vec.get(key)
        return u is not None and all(u[i] == v[i] for i in range(len(v)) if i != d) and u[d] != v[d]

    def _established(self, r: tuple) -> bool:
        if r in self.broken:
            return False
        tied, differed = self.rules.get(r, (0, 0))
        return differed == 0 and tied >= int(tune("dups.min_ties"))

    # ---- deciding ----------------------------------------------------------------------------------------------------
    def predicted_twin(self, params) -> Optional[Tuple[tuple, str]]:
        """(rule, twin key) when one established change separates `params` from a config already sent."""
        v = self._vec(params)
        for d, name in enumerate(self.names):
            b = self.sent[d].get(self._proj(v, d))
            if not b:
                continue
            for a, key in b.items():
                if a == v[d] or key in self.parked_rule:
                    continue
                r = self.rule(name, a, v[d])
                if self._established(r) and self._same_but(v, key, d):
                    return r, key
        return None

    def decide(self, config, tag: str, best: Optional[float]) -> bool:
        """True: park `config` (it is not sent now). False: send it (not predicted, near the best, or an audit)."""
        if not int(tune("dups.enabled")):
            return False
        hit = self.predicted_twin(config.params)
        if hit is None:
            return False
        r, twin = hit
        s = self.score_of.get(twin)
        if s is not None and best is not None and s >= best - float(tune("dups.near")):
            self.n_near += 1
            return False
        if self.rng.random() < float(tune("dups.audit")):
            self.n_audits += 1
            return False
        self.parked[r][config.key] = (config, tag)
        self.parked_rule[config.key] = r
        if s is None:
            self.waiting_on[twin].append(config.key)
        self.n_parked += 1
        return True

    # ---- learning ----------------------------------------------------------------------------------------------------
    def note_sent(self, config) -> None:
        if config.key in self.parked_rule:
            return
        v = self._vec(config.params)
        self.vec.setdefault(config.key, v)
        for d in range(len(self.names)):
            self.sent[d][self._proj(v, d)].setdefault(v[d], config.key)

    def note_score(self, config, score: float, best: Optional[float]) -> List[tuple]:
        """A scored config. Returns the (config, tag) pairs to send now: parked under a rule this score breaks, or
        waiting on this config as their twin when it scores within dups.near of the best."""
        if config.key in self.score_of:
            return []
        self.score_of[config.key] = score
        v = self._vec(config.params)
        self.vec.setdefault(config.key, v)
        out: List[tuple] = []
        for d, name in enumerate(self.names):
            b = self.scored[d][self._proj(v, d)]
            for a, key in b.items():
                if a == v[d] or not self._same_but(v, key, d):
                    continue
                r = self.rule(name, a, v[d])
                if abs(self.score_of[key] - score) < TIE:
                    self.rules[r][0] += 1
                else:
                    self.rules[r][1] += 1
                    if r not in self.broken:
                        self.broken.add(r)
                        self.n_broken += 1
                        out.extend(self._release_rule(r))
            b.setdefault(v[d], config.key)
        near = best is not None and score >= best - float(tune("dups.near"))
        for k in self.waiting_on.pop(config.key, []):
            if near and k in self.parked_rule:
                out.append(self._release_one(k))
                self.n_near += 1
        return out

    def _release_rule(self, r: tuple) -> List[tuple]:
        got = list(self.parked.pop(r, {}).values())
        for c, _ in got:
            self.parked_rule.pop(c.key, None)
        self.n_released += len(got)
        return got

    def _release_one(self, key: str) -> tuple:
        r = self.parked_rule.pop(key)
        item = self.parked[r].pop(key)
        self.n_released += 1
        return item

    def still_parked(self) -> int:
        return len(self.parked_rule)

    def summary(self) -> Optional[str]:
        if not (self.n_parked or self.n_audits):
            return None
        return ("duplicates: %d configs parked as predicted repeats of a callset; %d audited, %d sent because their twin "
                "is near the best, %d rules broken, %d sent after a broken rule or a near twin, %d still parked"
                % (self.n_parked, self.n_audits, self.n_near, self.n_broken, self.n_released, self.still_parked()))
