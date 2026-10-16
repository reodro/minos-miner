"""Facts a finished run gives the search, applied to later proposals. Owner: search.

max_num_haplotypes_in_population: the assembly reports the smallest cap that reproduces it. A later proposal
above that cap is rewritten to a cap that already ran at or above it. Every cap at or below it stays, because a
tighter cap is a different assembly and can score higher. A later report never lowers that floor: lowering it
would start rewriting caps that were still distinct. A later report that raises the floor replaces a stored cap
that now sits below the floor, so a higher cap is never rewritten to a smaller, different assembly.

max_alternate_alleles is not rewritten. On the stored VCFs a higher cap above one that had left alleles
unused was a different callset 1,137 times out of 16,700. A count of ALT alleles can be too small for the
cap that produced them.

The calling threshold is not rewritten from equal scores. The score can fall and rise again as calls drop, so
two equal ends do not mean the values between them are the same run.
"""
from __future__ import annotations

from typing import Dict, Optional, Tuple

from common.types import Config

MNH = "max_num_haplotypes_in_population"


def _sig(params: dict, drop: str) -> str:
    return Config("gatk", {k: v for k, v in params.items() if k != drop}).key


class Reductions:
    def __init__(self):
        self.hap: Dict[str, Tuple[int, int]] = {}          # signature without the cap -> (needed, cap that ran)

    def note(self, config: Config, score: Optional[float] = None, haplotype_cap: Optional[int] = None) -> Optional[str]:
        """Record what this run taught. A short phrase when the haplotype floor rises, else None.

        `score` is accepted from the engine and deliberately unused. Equal scores are not proof that every
        value between two settings is the same run."""
        del score
        params = config.params
        if haplotype_cap is None or isinstance(params.get(MNH), bool) or not isinstance(params.get(MNH), int):
            return None
        ran = int(params[MNH])
        if ran < haplotype_cap:
            return None
        sig = _sig(params, MNH)
        prev = self.hap.get(sig)
        if prev is not None and haplotype_cap <= prev[0]:
            return None
        # The stored cap is one that already ran at or above this floor. A report of "needed 0" can name a
        # small cap; a later real floor must not keep rewriting higher caps down to that small one.
        rep = ran if prev is None or prev[1] < haplotype_cap else prev[1]
        self.hap[sig] = (haplotype_cap, rep)
        return f"haplotype cap {haplotype_cap} (ran {rep})"

    def adjust(self, params: dict) -> dict:
        """A copy of `params`. A haplotype cap proven identical to a finished run is rewritten to that run.
        Caps at or below the reported floor are returned unchanged."""
        p = dict(params)
        hap = p.get(MNH)
        known = self.hap.get(_sig(p, MNH)) if isinstance(hap, int) and not isinstance(hap, bool) else None
        if known is not None and hap > known[0]:
            p[MNH] = known[1]
        return p
