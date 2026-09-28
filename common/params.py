"""Validate configs against the platform's published parameter ranges.

The platform rejects (or defaults) any out-of-range value, which silently turns a
searched config into a different one. Every Config main submits or evaluator runs
is validated here first.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Dict, List

from common.types import Config


@dataclass(frozen=True)
class ParamSpec:
    name: str
    type: str                   # "int" | "float" | "bool" | "enum"
    default: object
    min: object = None
    max: object = None
    allowed: tuple = ()

    def check(self, v) -> str:
        """Return "" if valid, else a reason."""
        if self.type == "bool":
            if not isinstance(v, bool):
                return f"{self.name}: expected bool, got {v!r}"
        elif self.type in ("int", "float"):
            if isinstance(v, bool) or not isinstance(v, (int, float)):
                return f"{self.name}: expected {self.type}, got {v!r}"
            if self.type == "int" and not float(v).is_integer():
                return f"{self.name}: expected int, got {v!r}"
            if self.min is not None and v < self.min or self.max is not None and v > self.max:
                return f"{self.name}: {v!r} outside [{self.min}, {self.max}]"
        if self.allowed and v not in self.allowed:
            return f"{self.name}: {v!r} not in {list(self.allowed)}"
        return ""


# Combinations the platform's per-parameter ranges accept but GATK 4.5.0.0 rejects (the validator then scores
# nothing). Each verified by a real GATK run on 2026-09-28.
CROSS_CONSTRAINTS = {
    "gatk": [
        ("min_assembly_region_size <= max_assembly_region_size",
         lambda p: p.get("min_assembly_region_size", 50) <= p.get("max_assembly_region_size", 300)),
        ("base_quality_score_threshold >= 6",
         lambda p: p.get("base_quality_score_threshold", 18) >= 6),
    ],
}


class ParamRanges:
    def __init__(self, raw: dict):
        self.raw = raw
        self.tools: Dict[str, Dict[str, ParamSpec]] = {}
        for tool, spec in raw["tools"].items():
            self.tools[tool] = {
                p["name"]: ParamSpec(p["name"], p["type"], p.get("default"), p.get("min"), p.get("max"),
                                     tuple(p.get("allowed_values") or ()))
                for p in spec["parameters"]
            }

    @classmethod
    def from_file(cls, path: str) -> "ParamRanges":
        with open(path) as f:
            return cls(json.load(f))

    @classmethod
    def from_platform(cls) -> "ParamRanges":
        from common.platform import Platform
        return cls(Platform().parameter_ranges())

    def errors(self, config: Config) -> List[str]:
        specs = self.tools.get(config.tool)
        if specs is None:
            return [f"unknown tool {config.tool!r}"]
        errs = [f"unknown parameter {k!r}" for k in config.params if k not in specs]
        errs += [e for k, v in config.params.items() if k in specs for e in [specs[k].check(v)] if e]
        if not errs:
            params = {**{n: s.default for n, s in specs.items()}, **config.params}
            errs += [f"GATK rejects: needs {desc}" for desc, ok in CROSS_CONSTRAINTS.get(config.tool, []) if not ok(params)]
        return errs

    def defaults(self, tool: str) -> Config:
        return Config(tool, {n: s.default for n, s in self.tools[tool].items()})
