"""Control every sub-VPS from main, over SSH, using the inventory (~/.minos/fleet.json).

    python -m evaluator.fleet.control status              reachable? cores, workers, docker, patch jar, genomes, code, disk
    python -m evaluator.fleet.control provision [HOST..]  copy the repo + main's genomes, run scripts/setup_node.sh
    python -m evaluator.fleet.control sync [HOST..]       copy the repo + self-test bundle, rebuild the patch jar
    python -m evaluator.fleet.control selftest [HOST..]   run the node self-test everywhere
    python -m evaluator.fleet.control reap [HOST..]       kill GATK containers left by a crashed node (refuses while an
                                                          agent is running there)

Hosts run in parallel. Keys come from key_filename; a key_passphrase in the inventory is unlocked in a private
ssh-agent (see inventory.ensure_agent).
"""
from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from typing import Callable, List

from common.settings import SETTINGS
from evaluator.fleet import selftest
from evaluator.fleet.coordinator import push
from evaluator.fleet.inventory import REPO, Host, code_digest, ensure_agent, load

EXCLUDES = ["--exclude", ".git/", "--exclude", ".venv/", "--exclude", "test/", "--exclude", "__pycache__/",
            "--exclude", "evaluator/gatk_patch/build/", "--exclude", "*.pyc"]


def sh(h: Host, cmd: str, timeout: int = 600) -> subprocess.CompletedProcess:
    return subprocess.run(h.ssh(cmd), capture_output=True, text=True, timeout=timeout)


def copy_repo(h: Host) -> None:
    """This repo and main's copy of the upstream minos_subnet (score formula, templates), so a node can never run a
    different upstream version than main."""
    up = SETTINGS.minos_subnet_dir
    sh(h, f"mkdir -p {shlex.quote(h.repo)} {shlex.quote(up)}")
    for src, dst, extra in ((REPO, h.repo, EXCLUDES),
                            (up, up, ["--exclude", ".git/", "--exclude", ".venv/", "--exclude", "datasets/",
                                      "--exclude", "__pycache__/"])):
        p = subprocess.run(h.rsync([f"{src}/"], f"{dst}/", "--delete", *extra), capture_output=True, text=True,
                           timeout=1800)
        if p.returncode != 0:
            raise RuntimeError(p.stderr[-300:])


def status(h: Host) -> dict:
    py = f"{h.repo}/.venv/bin/python"
    script = (f"echo cores=$(nproc); echo ram_gb=$(free -g | awk '/Mem/{{print $2}}'); "
              f"echo disk_free=$(df -BG --output=avail {h.data_dir} 2>/dev/null | tail -1 | tr -d ' '); "
              f"echo gatk_image=$(docker image inspect broadinstitute/gatk:4.5.0.0 >/dev/null 2>&1 && echo yes || echo no); "
              f"echo patch_jar=$(test -f {h.repo}/evaluator/gatk_patch/build/minos-patch.jar && echo yes || echo no); "
              f"echo genomes=$(ls -d {h.data_dir}/reference/chr*/ 2>/dev/null | wc -l); "
              f"cd {h.repo} 2>/dev/null && {py} -c 'from evaluator.fleet.inventory import code_digest; "
              f"from common.settings import default_workers; print(\"code=\" + code_digest()); "
              f"print(\"workers=\" + str(default_workers()))' 2>/dev/null; true")
    p = sh(h, script, timeout=60)
    if "cores=" not in p.stdout:
        return {"host": h.name, "reachable": False, "error": p.stderr.strip()[-160:]}
    v = dict(line.split("=", 1) for line in p.stdout.splitlines() if "=" in line)
    return {"host": h.name, "reachable": True, "cores": v.get("cores"), "ram_gb": v.get("ram_gb"),
            "disk_free": v.get("disk_free") or "?", "gatk_image": v.get("gatk_image") == "yes",
            "patch_jar": v.get("patch_jar") == "yes", "genomes": v.get("genomes"),
            "code_same_as_main": v.get("code") == code_digest(), "workers": h.workers or v.get("workers", "?")}


def push_genomes(h: Host) -> None:
    """Main's reference genomes (FASTA, index, dict, SDF) to the host, every chromosome main has; the node's setup
    then downloads only what is still missing. No dependence on the platform's endpoint for what main already has."""
    sh(h, f"mkdir -p {shlex.quote(h.data_dir + '/reference')}")
    p = subprocess.run(h.rsync([f"{SETTINGS.data_dir}/reference/"], f"{h.data_dir}/reference/"), capture_output=True,
                       text=True, timeout=7200)
    if p.returncode != 0:
        raise RuntimeError(f"pushing genomes: {p.stderr[-300:]}")


def provision(h: Host) -> dict:
    copy_repo(h)
    push_genomes(h)
    p = sh(h, f"cd {shlex.quote(h.repo)} && DATA_DIR={shlex.quote(h.data_dir)} bash scripts/setup_node.sh", timeout=7200)
    s = sync(h)
    ok = p.returncode == 0 and s["ok"]
    return {"host": h.name, "ok": ok, "tail": (p.stdout + p.stderr).strip()[-300:] if p.returncode else s["tail"]}


def sync(h: Host) -> dict:
    copy_repo(h)
    b = selftest.make_bundle(SETTINGS.data_dir, SETTINGS.rounds_dir)
    push(h, [b["bam"], b["bai"]], f"{h.data_dir}/selftest/")
    p = sh(h, f"cd {shlex.quote(h.repo)} && bash evaluator/gatk_patch/build.sh", timeout=1200)
    return {"host": h.name, "ok": p.returncode == 0, "tail": (p.stdout + p.stderr).strip()[-200:]}


def run_selftest(h: Host) -> dict:
    repo = shlex.quote(h.repo)
    p = sh(h, f"cd {repo} && DATA_DIR={shlex.quote(h.data_dir)} {repo}/.venv/bin/python -m evaluator.fleet.node "
              f"--selftest-only", timeout=900)
    try:
        return {"host": h.name, **json.loads(p.stdout.strip().splitlines()[-1])}
    except (IndexError, json.JSONDecodeError):
        return {"host": h.name, "ok": False, "error": (p.stderr or p.stdout).strip()[-200:]}


def reap(h: Host) -> dict:
    """Kill every node GATK container on the host - only when no node agent is running there (a running agent's
    containers are live jobs; stop main first)."""
    busy = sh(h, "pgrep -f 'evaluator.fleet.node --stdi[o]' >/dev/null && echo busy || echo idle").stdout.strip()
    if busy == "busy":
        return {"host": h.name, "killed": 0, "skipped": "a node agent is running there (stop main first)"}
    p = sh(h, "ids=$(docker ps -q --filter name=minos-node-); [ -n \"$ids\" ] && docker kill $ids | wc -l || echo 0")
    return {"host": h.name, "killed": p.stdout.strip()}


def on_hosts(fn: Callable[[Host], dict], names: List[str]) -> List[dict]:
    hosts = [h for h in load() if not names or h.name in names or h.host in names]
    if not hosts:
        sys.exit("no matching enabled hosts in the inventory")
    ensure_agent(hosts)

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
