"""Multivariate TPE beside RF-UCB. Optuna samples the mixed GATK space; each suggestion is expanded into the
same designed cheap family, and the trial is told the best score in that family.

TPE is the exploration lane: it is built for mixed int / float / enum spaces and for pairs of parameters that
have to move together. It reads nothing the forest does not also see, because both propose through one client.
"""
from __future__ import annotations

import warnings
from typing import List, Optional, Set, Tuple

import optuna
from optuna.exceptions import ExperimentalWarning
from optuna.trial import TrialState

warnings.filterwarnings("ignore", category=ExperimentalWarning)

from common.params import ParamRanges
from common.types import FULL, Config
from search.base import Strategy
from search.client import EvalClient
from search.family import designed_family
from search.schedule import run_streaming
from search.space import Space

optuna.logging.set_verbosity(optuna.logging.WARNING)


def tpe(ranges: ParamRanges, *, tool: str = "gatk", in_flight: int = 4, family: int = 6,
        family_max_cost: float = 0.15, n_startup: int = 10, late_s: float = 180.0, max_pending: int = 256,
        seed: int = 0, tag: str = "tpe") -> Strategy:
    space = Space(ranges, tool)

    async def run(client: EvalClient) -> None:
        proposed: Set[str] = set()
        study = optuna.create_study(
            direction="maximize",
            sampler=optuna.samplers.TPESampler(multivariate=True, constant_liar=True, n_startup_trials=n_startup,
                                               seed=seed))
        open_batches: List[Tuple[optuna.trial.Trial, List]] = []

        def fresh(c: Optional[Config]) -> bool:
            return c is not None and c.key not in proposed and (c.key, FULL) not in client.latest

        def late() -> bool:
            left = client.time_left()
            return left is not None and left < late_s

        def suggest() -> Optional[Config]:
            for _ in range(20):
                trial = study.ask()
                params = {}
                for name in space.names:
                    spec = space.specs[name]
                    if name == "base_quality_score_threshold":
                        params[name] = trial.suggest_int(name, 6, int(spec.max))
                    elif spec.allowed:
                        params[name] = trial.suggest_categorical(name, list(spec.allowed))
                    elif spec.type == "bool":
                        params[name] = trial.suggest_categorical(name, [False, True])
                    elif spec.type == "int":
                        params[name] = trial.suggest_int(name, int(spec.min), int(spec.max))
                    elif space._log(spec):
                        params[name] = trial.suggest_float(name, float(spec.min), float(spec.max), log=True)
                    else:
                        params[name] = trial.suggest_float(name, float(spec.min), float(spec.max))
                config = space._assemble(params)
                # Repair can change a value (assembly bounds). Tell the failure instead of scoring a different
                # config than the one TPE thinks it asked for.
                unchanged = config is not None and all(config.params.get(n) == v for n, v in params.items())
                if unchanged and fresh(config):
                    open_batches.append((trial, []))
                    return config
                study.tell(trial, state=TrialState.FAIL)
            return None

        def fill():
            if late():
                return []
            config = suggest()
            if config is None:
                return []
            fam = designed_family(
                space, config, family,
                lambda c: fresh(c) and c.key != config.key and client.variant_cost(config, c) <= family_max_cost)
            batch = [config] + fam
            proposed.update(c.key for c in batch)
            futs = list(client.propose_many(batch, tag=tag))
            trial, _ = open_batches[-1]
            open_batches[-1] = (trial, futs)
            return futs

        def on_done(_fut):
            still = []
            for trial, futs in open_batches:
                if not futs:
                    study.tell(trial, state=TrialState.FAIL)
                    continue
                if not all(f.done() for f in futs):
                    still.append((trial, futs))
                    continue
                scores = []
                for f in futs:
                    outcome = f.result()
                    if outcome.ok and outcome.result.gate_passed is not False and outcome.score is not None:
                        scores.append(outcome.score)
                if scores:
                    study.tell(trial, max(scores))
                else:
                    study.tell(trial, state=TrialState.FAIL)
            open_batches[:] = still

        def room(pending) -> bool:
            cap = max(max_pending, (family + 1) * (client.capacity() or 0))
            if len(pending) >= cap:
                return False
            idle = client.idle_slots()
            return (len(pending) < in_flight) if idle is None else (idle > 0 or not pending)

        await run_streaming(client, room, fill, on_done)

    return run
