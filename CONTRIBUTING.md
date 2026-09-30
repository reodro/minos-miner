# Working on minos-miner

Four parts, three people. These rules keep us from breaking each other - and, above all, from breaking
**exactness**: every number we act on must equal what the validator would compute.

## 1. Exactness first

- Anything that touches GATK execution or scoring must pass a **bit-identical test against separate stock runs**
  (or the validator-exact scoring path) before it merges. "Looks the same" is not a test; compare every call
  (position, alleles, QUAL, genotype) and, where possible, the complete VCF lines.
- **Never reimplement the score formula.** Use `common/scoring.py`, which wraps upstream
  `minos_subnet/utils/scoring.py`.
- A fast path that can be wrong must **detect** when it is wrong and fall back (example: batch variants that
  diverge are flagged and repaired, never silently returned).

## 2. Measure it ourselves

- Numbers from the previous project (`/minos/legacy`) are **ideas, not evidence**. Several of its "verified"
  results were wrong when we measured. Reproduce before relying on anything.
- Noise floor: identical full search runs differed by ~0.08 points. Do not claim an effect below ~0.1 from a
  single pair; replicate.
- Backtest KPIs must not use information the live miner will not have (example: `N` is estimated, not read from
  the released `mutations.vcf`).

## 3. Where things go

| What | Where | In git? |
|---|---|---|
| Part code | `common/`, `answer_set/`, `evaluator/`, `search/`, `main/`, `backtest/` | yes |
| Unit tests (contracts, fast, no network/GATK unless `RUN_SLOW=1`) | `tests/` | yes |
| Experiments, benchmarks, exactness tests, scratch output | `test/` | **no** (ignored) |
| Data: BAMs, references, callsets, caches | `$DATA_DIR` (default `/minos/data`) | **no** |
| Secrets (tokens, wallets, keys) | never in the repo; the fleet's SSH key passphrase lives only in `~/.minos/fleet.json` on main (mode 600) | **never** |

## 4. How the parts may depend on each other

- `common/types.py` is the **contract**. Changing it needs agreement from all owners (it breaks the others).
- `answer_set`, `evaluator` and `search` do not import each other's internals. `main` wires them together through
  their public APIs (`__init__.py`). `search` talks to evaluation only through `EvalClient`.
- If you need something from another part, add it to that part's public API (ask its owner) instead of reaching in.

## 5. Git

- One branch per task, then a pull request. Another member reviews before merging to `main`.
- `.venv/bin/python -m unittest discover -s tests -t .` must pass.
- Commit messages: first line says what changed; the body says why and what was measured (numbers, test names).
- `CODEOWNERS` routes reviews to the part's owner.

## 6. Libraries

- If a standard library covers it (numpy, scikit-learn, scipy, optuna, pymoo/DEAP, ...), **propose adding it**
  instead of hand-writing the algorithm. Add it to `requirements.txt` (pinned) so every machine installs the same.

## 7. Machines

- Every machine (main VPS and fleet nodes) runs **floor(cores x 0.95 - 1)** GATK workers, 1 thread each
  (`common.settings.default_workers()`, override with `MINOS_WORKERS`).
- Same GATK image (`broadinstitute/gatk:4.5.0.0`) and the same patch jar everywhere; a node that fails its
  self-test must not be used.
- Long jobs run in the background; never block a shared session waiting on them.

## 8. Before you change something, ask

- Does it keep validator parity (calls and score identical)?
- Does it change the **order of random draws** inside GATK? (Batching depends on it - see evaluator/README.md.)
- Does it cost round time? Anything slow at round start delays the first evaluation.
- Does it work on a fleet node, not just on the main VPS?
- Does the search choose configs by **score**? Never restrict the config space for the evaluator's convenience -
  only exclude settings that can never score (see search/README.md).
