"""True scores (the released truth) for a run's top results, through main.review's own path: the config's VCF on
main if it is there, else GATK reruns it here, then the validator's scorer. Answers "did we submit the right
config?" beyond the review's top configs (which are often ties: configs with identical calls).

    .venv/bin/python scripts/true_scores.py <run_id> [--top N]      the top N distinct local scores (default 7)
    .venv/bin/python scripts/true_scores.py <run_id> <config_key>...

run_id: <tag>-<epoch>, as in status/runs/ (main.review.find_run accepts a round tag for its newest run)."""
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from answer_set import OracleProvider  # noqa: E402
from backtest.data import load  # noqa: E402
from main.review import Reviewer, _reference, find_run  # noqa: E402

args = sys.argv[1:]
run_id = find_run(args[0])
rv = Reviewer()
m = json.loads((rv.dir / "runs" / run_id / "manifest.json").read_text())
configs = {k: v for k, v in (m.get("configs") or {}).items() if v.get("local") is not None}
if len(args) > 1 and args[1] != "--top":
    keys = args[1:]
else:
    n = int(args[2]) if len(args) > 2 else 7
    keys, seen = [], set()
    for k in sorted(configs, key=lambda x: -configs[x]["local"]):
        s = round(configs[k]["local"], 9)                  # identical calls give identical scores: one of each
        if s not in seen:
            seen.add(s)
            keys.append(k)
        if len(keys) == n:
            break
rnd = load(m["tag"], reference=_reference(m))
oracle = OracleProvider().build(rnd)
t0 = time.time()
true, gate, rerun = rv._true_scores(rnd, oracle, keys, configs, m)
print(f"{run_id}: {len(true)} of {len(keys)} scored ({rerun} rerun here) in {time.time() - t0:.0f} s; "
      f"submitted {m.get('final_config')}")
print(f"{'#':>2}  {'config':16}  {'local':>9}  {'true':>9}  note")
for i, k in enumerate(sorted(keys, key=lambda x: -configs[x]["local"]), 1):
    t = true.get(k)
    note = "submitted" if k == m.get("final_config") else ""
    print(f"{i:>2}  {k:16}  {configs[k]['local']:9.4f}  {t if t is None else f'{t:9.4f}'}  {note}")
