#!/usr/bin/env python3
"""Write search/engine_tag.py so the next fleet_replay syncs this engine tag."""
from __future__ import annotations

import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
PATH = ROOT / "search" / "engine_tag.py"

TAGS = {
    "v9.16": (False, "far seat only (B27)"),
    "v9.17": (True, "compass-then-cube, far seat, elite leftover GA (B28+B27+B26)"),
    "v9.18": (True, "v9.16+v9.17, leftover fill, far seat stays far (B29+B30)"),
    "v9.19": (False, "four callset piles, gap 0.10, half the GA looks for a new fingerprint"),
    "v9.20": (False, "32 parents per pile; a higher returned run takes the centre; fingerprint history blocks a similar pile"),
    "v9.21": (False, "the highest score within 0.10 of a pile becomes its reference"),
    "v9.22": (False, "six piles, four climb, two candidates take a free seat; the highest score within 0.10 becomes the climber's centre"),
    "v9.23": (False, "climbers keep their seats until the fine-tune; the fine-tune does not breed; a higher candidate takes the lowest climber; the second half of Step B widens the grid to five rungs"),
    "v9.24": (False, "a far config above the lowest climbing pile takes that climber immediately; the fine-tune restarts from a higher score"),
    "v9.25": (False, "while GA is still moving a candidate the climbers send nothing and the fleet breeds; then one rung"),
    "v9.26": (False, "a peaked centre is not taken; elite probes of the best run once, then the grid keeps the fleet"),
    "v9.27": (False, "elite probes go out even when close; a free worker then walks further up the line"),
    "v9.28": (False, "while a grid is walking, the pile does not move the centre; the grid writes the setting"),
    "v9.29": (False, "one climber on the highest fingerprint, a quadratic dummy, and a grid that writes the centre when it finishes"),
    "v9.30": (False, "one compass on the highest score, a joint cube, and a hill lane that only looks 0.10 away"),
    "v9.31": (False, "the hill lane crosses two piles; the compass and the cube stay as in v9.30"),
    "v9.32": (False, "the hill-lane cross keeps each verified linkage group with one parent"),
    "v9.33": (False, "the compass keeps 16 cores free and the population mates while it runs"),
    "v9.34": (False, "the compass stays, and the pile GA breeds beside it"),
}


def main() -> None:
    if len(sys.argv) != 2 or sys.argv[1] not in TAGS:
        sys.exit(f"usage: {sys.argv[0]} v9.16|v9.17|v9.18|v9.19|v9.20|v9.21|v9.22|v9.23|v9.24|v9.25|v9.26|v9.27|v9.28|v9.29|v9.30|v9.31|v9.32|v9.33|v9.34")
    ver = sys.argv[1]
    elite, why = TAGS[ver]
    PATH.write_text(
        '"""Which engine tag the next fleet_replay syncs."""\n'
        f"VERSION = {ver!r}\n"
        f"ELITE_LEFTOVER = {elite}\n"
    )
    print(f"{ver}: ELITE_LEFTOVER={elite} ({why})")


if __name__ == "__main__":
    main()
