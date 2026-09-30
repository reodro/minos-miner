"""Control every sub-VPS from main, over SSH, using the inventory (~/.minos/fleet.json).

    python -m evaluator.fleet.control status              reachable? cores, workers, docker, patch jar, genomes, code, disk
    python -m evaluator.fleet.control provision [HOST..]  copy the repo, run scripts/setup_node.sh (fresh machines)
    python -m evaluator.fleet.control sync [HOST..]       copy the repo + self-test bundle, rebuild the patch jar
    python -m evaluator.fleet.control selftest [HOST..]   run the node self-test everywhere
    python -m evaluator.fleet.control reap [HOST..]       kill stray GATK containers left by a crashed node

Hosts run in parallel. Nothing here stores secrets: keys come from key_filename / ssh-agent.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from typing import Callable, List

from common.settings import SETTINGS
from evaluator.fleet import selftest
from evaluator.fleet.coordinator import push
from evaluator.fleet.inventory import REPO, Host, code_digest, load

EXCLUDES = ["--exclude", ".git/", "--exclude", ".venv/", "--exclude", "test/", "--exclude", "__pycache__/",
            "--exclude", "evaluator/gatk_patch/build/", "--exclude", "*.pyc"]


def sh(h: Host, cmd: str, timeout: int = 600) -> subprocess.CompletedProcess:
    return subprocess.run(h.ssh(cmd), capture_output=True, text=True, timeout=timeout)


def copy_repo(h: Host) -> None:
    sh(h, f"mkdir -p {h.repo}")
    p = subprocess.run(h.rsync([f"{REPO}/"], f"{h.repo}/", "--delete", *EXCLUDES), capture_output=True, text=True,
                       timeout=1800)
    if p.returncode != 0:
        raise RuntimeError(p.stderr[-300:])


def status(h: Host) -> dict:
    py = f"{h.repo}/.venv/bin/python"
    script = (f"nproc; free -g | awk '/Mem/{{print $2}}'; df -BG --output=avail {h.data_dir} 2>/dev/null | tail -1; "
              f"docker image inspect broadinstitute/gatk:4.5.0.0 >/dev/null 2>&1 && echo gatk || echo no-gatk; "
              f"test -f {h.repo}/evaluator/gatk_patch/build/minos-patch.jar && echo jar || echo no-jar; "
              f"ls -d {h.data_dir}/reference/chr*/ 2>/dev/null | wc -l; "
              f"cd {h.repo} 2>/dev/null && {py} -c 'from evaluator.fleet.inventory import code_digest; "
              f"from common.settings import default_workers; print(code_digest()); print(default_workers())' 2>/dev/null")
    p = sh(h, script, timeout=60)
    if p.returncode not in (0, 1) and not p.stdout:
        return {"host": h.name, "reachable": False, "error": p.stderr.strip()[-120:]}
    v = p.stdout.split("\n")
    get = lambda i: v[i].strip() if len(v) > i else "?"  # noqa: E731
    return {"host": h.name, "reachable": True, "cores": get(0), "ram_gb": get(1), "disk_free": get(2),
            "gatk_image": get(3) == "gatk", "patch_jar": get(4) == "jar", "genomes": get(5),
            "code_same_as_main": get(6) == code_digest(), "workers": h.workers or get(7)}


def provision(h: Host) -> dict:
    copy_repo(h)
    p = sh(h, f"cd {h.repo} && DATA_DIR={h.data_dir} bash scripts/setup_node.sh", timeout=7200)
    sync(h)
    return {"host": h.name, "ok": p.returncode == 0, "tail": (p.stdout + p.stderr).strip()[-300:]}


def sync(h: Host) -> dict:
    copy_repo(h)
    b = selftest.make_bundle(SETTINGS.data_dir, SETTINGS.rounds_dir)
    push(h, [b["bam"], b["bai"]], f"{h.data_dir}/selftest/")
    p = sh(h, f"cd {h.repo} && bash evaluator/gatk_patch/build.sh", timeout=1200)
    return {"host": h.name, "ok": p.returncode == 0, "tail": (p.stdout + p.stderr).strip()[-200:]}


def run_selftest(h: Host) -> dict:
    p = sh(h, f"cd {h.repo} && DATA_DIR={h.data_dir} {h.repo}/.venv/bin/python -m evaluator.fleet.node --selftest-only",
           timeout=900)
    try:
        return {"host": h.name, **json.loads(p.stdout.strip().splitlines()[-1])}
    except (IndexError, json.JSONDecodeError):
        return {"host": h.name, "ok": False, "error": (p.stderr or p.stdout).strip()[-200:]}


def reap(h: Host) -> dict:
    p = sh(h, "ids=$(docker ps -q --filter name=minos-node-); [ -n \"$ids\" ] && docker kill $ids | wc -l || echo 0")
    return {"host": h.name, "killed": p.stdout.strip()}


def on_hosts(fn: Callable[[Host], dict], names: List[str]) -> List[dict]:
    hosts = [h for h in load() if not names or h.name in names or h.host in names]
    if not hosts:
        sys.exit("no matching enabled hosts in the inventory")

    def safe(h):
        try:
            return fn(h)
        except Exception as e:  # noqa: BLE001 - report per host, keep going
            return {"host": h.name, "ok": False, "error": str(e)[-200:]}
    with ThreadPoolExecutor(max(1, len(hosts))) as pool:
        return list(pool.map(safe, hosts))


def main() -> None:
    ap = argparse.ArgumentParser(prog="python -m evaluator.fleet.control")
    ap.add_argument("cmd", choices=["status", "provision", "sync", "selftest", "reap"])
    ap.add_argument("hosts", nargs="*", help="limit to these hosts (default: every enabled host)")
    a = ap.parse_args()
    fn = {"status": status, "provision": provision, "sync": sync, "selftest": run_selftest, "reap": reap}[a.cmd]
    for row in on_hosts(fn, a.hosts):
        print(json.dumps(row))


if __name__ == "__main__":
    main()
