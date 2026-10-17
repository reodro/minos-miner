"""The second score: the base score's order, with no flat (2026-10-16). Owner: search.

Three rules, from the search designer:
  1. where the base score S rises or falls, the second score rises or falls the same way;
  2. S(a) > S(b) means second(a) > second(b);
  3. two different configs never share a second score.
Rules 2 and 3 leave one choice: the order inside a set of configs that share S. So the second score is S followed by
tie-breakers, compared left to right: a tuple, not one float. Two different scores can be 1.6e-12 apart (2.9M stored
rows), and a float near 90 moves in steps of 1.4e-14, so one float holds about 100 tie-break levels while one flat
holds up to a thousand configs.

    second(c) = (S, G, H)
    S  the base score exactly as the validator computes it (common.scoring); 0 for a config that fails the gate.
    G  (0, 0, 0) when the gate passes. When it fails: minus the ratios the callset cannot measure, minus how far ti/tv
       and het/hom sit outside their limits (each in units of its limit), then the score without the gate. The
       validator's flat 0 becomes a slope toward the gate: a violated constraint ordered by its violation (Deb 2000).
    H  the config's key, then its settings: two different configs never tie.

What the stored rounds say about the order inside a flat (4 rounds, 880K configs, 2026-10-16):
  - a tie is always the same TP/FP/FN in every class; 35% of tied pairs have byte-identical VCFs, 63% differ in QUAL;
  - near the top a strict peak is rare (under 1% of configs): the climb stops on flats;
  - the same move from two tied configs gives the same score 96-99.7% of the time: a flat's members act as one point;
  - no tie-breaker ranked a flat's members better than chance: QUAL soft score 0.51-0.53, Gaussian process mean
    0.44-0.53, nearest higher config 0.46-0.50, random 0.46-0.52.
So H only makes the order strict; it says nothing about where a flat's exit is. A climber should not count an H-only
gain as progress (`Second.informative` drops H): a hash order on a flat makes about 1/(d+1) of its points false peaks
(Neidhart, Szendro and Krug 2014), and a tie-breaker that does not point at the exit can slow a search down
exponentially (Brockhoff et al. 2009). On a flat the exit is an untried move, not a sideways one.
"""
from __future__ import annotations

import json
from typing import NamedTuple, Optional, Tuple

from common.scoring import score_from_parts
from common.types import FULL, Config, ScoreResult
from utils.scoring import GATE_HETHOM_MAX_DELTA, GATE_TITV_MAX_DELTA, _usable_ratio  # upstream, on sys.path via common.scoring

GATE_CHECKS = (("titv_truth_snp", "titv_query_snp", GATE_TITV_MAX_DELTA),
               ("hethom_truth_snp", "hethom_query_snp", GATE_HETHOM_MAX_DELTA))
PASSED = (0, 0.0, 0.0)


class Second(NamedTuple):
    s: float                              # the base score; 0 when the gate fails
    gate: Tuple[int, float, float]        # PASSED, or (-unmeasured ratios, -distance outside the limits, ungated score)
    code: Tuple[str, str]                 # config key, then its settings

    @property
    def informative(self) -> tuple:
        """S and the gate slope: the part that carries evidence. H only makes the order strict."""
        return self.s, self.gate


def gate_slope(result: ScoreResult) -> Tuple[int, float, float]:
    """PASSED for a config that passes the gate; otherwise how far it is from passing, then its ungated score.

    Same reading of the ratios as upstream `plausibility_gate`: a truth ratio that cannot be measured skips its check,
    and a query ratio that cannot be measured while the truth's can is a failure (counted apart, it is worse than any
    distance)."""
    if result.gate_passed is not False:
        return PASSED
    ratios = result.ratios or {}
    unmeasured, outside = 0, 0.0
    for truth_key, query_key, limit in GATE_CHECKS:
        truth = _usable_ratio(ratios.get(truth_key))
        if truth is None:
            continue
        query = _usable_ratio(ratios.get(query_key))
        if query is None:
            unmeasured += 1
        else:
            outside += max(0.0, abs(query - truth) - limit) / limit
    ungated = score_from_parts(result.core or 0.0, result.region_fp / max(float(result.n_targets), 1.0))
    return -unmeasured, -outside, ungated


def code(config: Config) -> Tuple[str, str]:
    """The key, then the settings themselves (two configs with different settings never share this)."""
    return config.key, json.dumps(config.params, sort_keys=True, default=str)


def second(result: Optional[ScoreResult], config: Config) -> Optional[Second]:
    """The second score of a scored config; None when it has no score (crashed, or still running)."""
    if result is None or result.score is None:
        return None
    s = float(result.score) if result.gate_passed is not False else 0.0
    return Second(s, gate_slope(result), code(config))


def of(ctx, config: Config) -> Optional[Second]:
    """The second score from the round's latest full-region result."""
    return second(ctx.client.latest.get((config.key, FULL)), config)
