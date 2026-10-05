"""Pull one round's run history to main before the fleet is turned off: every node's callsets (VCFs, crash
stderr, JVM logs - they exist only on the nodes), plus main's run records, result store, fleet logs and the run's
log, into $DATA_DIR/history/<label>/ (<node>/ per node, main/ for main's side).

    .venv/bin/python scripts/pull_history.py <round_tag> <label> [run_log]
    e.g. scripts/pull_history.py 20261002T092900Z 20261002T092900Z_v4 test/runs/replay_20261002T092900Z_1924.log

A 72-minute run on ~800 workers is ~30-40 GB and takes ~5-7 minutes (nodes in parallel). Re-running it resumes."""
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.settings import SETTINGS  # noqa: E402
from evaluator.fleet.inventory import REPO, ensure_agent, load  # noqa: E402

tag, label = sys.argv[1], sys.argv[2]
run_log = sys.argv[3] if len(sys.argv) > 3 else None
root = Path(SETTINGS.data_dir) / "history" / label
root.mkdir(parents=True, exist_ok=True)
hosts = load()
ensure_agent(hosts)


def pull(h) -> str:
    dst = root / h.name.replace(":", "_")
    dst.mkdir(parents=True, exist_ok=True)
    ssh = " ".join(["ssh", *h.ssh_opts()])
    t0 = time.time()
    p = subprocess.run(["rsync", "-a", "--partial", "-e", ssh, f"{h.user}@{h.host}:{h.repo}/cache/callsets/{tag}/",
                        f"{dst}/"], capture_output=True, text=True, timeout=3600)
    n = sum(len(f) for _, _, f in os.walk(dst))
    return f"{h.name}: rc {p.returncode}, {n} files, {time.time() - t0:.0f} s {p.stderr.strip()[-120:]}"


with ThreadPoolExecutor(max(1, len(hosts))) as ex:
    for line in ex.map(pull, hosts):
        print(line, flush=True)
main = root / "main"
for d in (REPO / "status" / "runs").glob(f"{tag}-*"):
    shutil.copytree(d, main / "runs" / d.name, dirs_exist_ok=True)
(main / "results").mkdir(parents=True, exist_ok=True)
for f in (REPO / "cache" / "results").glob(f"{tag}*"):
    shutil.copy(f, main / "results")
if run_log:
    shutil.copy(run_log, main)
logs = Path.home() / ".minos" / "fleet" / "logs"
if logs.is_dir():
    shutil.copytree(logs, main / "fleet_logs", dirs_exist_ok=True)
subprocess.run(["du", "-sh", str(root)])
print("done - the nodes can be turned off", flush=True)
