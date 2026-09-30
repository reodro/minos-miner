"""The fleet inventory (~/.minos/fleet.json on main) and how main reaches each sub-VPS over SSH.

    [{"host": "1.2.3.4", "user": "root", "port": 22, "key_filename": "~/.ssh/id_ed25519",
      "workers": null, "enabled": true}]

- `workers`: null = the node's own rule, floor(cores x 0.95 - 1); a number overrides it.
- `port`: several machines may share one IP behind different ports (e.g. vast.ai).
- optional `repo` / `data_dir`: where the code and data live on that host (default /minos/minos-miner, /minos/data).
- `host: "local"`: run the node as a local subprocess instead of over SSH (tests, or this machine as a node).
- an entry without `user` / `key_filename` inherits them from the first entry that has them.
There is deliberately no passphrase field: unlock keys with ssh-agent. The file is re-read at every round start,
so hosts can be added or disabled between rounds.
"""
from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

INVENTORY = Path(os.environ.get("MINOS_FLEET_INVENTORY", str(Path.home() / ".minos" / "fleet.json")))
CONTROL_DIR = Path.home() / ".minos" / "ssh"                 # ssh ControlMaster sockets (one connection per host)
REPO = Path(__file__).resolve().parents[2]
CODE_DIRS = ("common", "answer_set", "evaluator", "search", "main", "backtest", "scripts")
CODE_SUFFIXES = (".py", ".java", ".json", ".sh")


@dataclass(frozen=True)
class Host:
    host: str
    user: str = "root"
    port: int = 22
    key_filename: Optional[str] = None
    workers: Optional[int] = None
    enabled: bool = True
    repo: str = "/minos/minos-miner"
    data_dir: str = "/minos/data"
    extra: dict = field(default_factory=dict, compare=False, hash=False)

    @property
    def name(self) -> str:
        return self.host if self.port == 22 else f"{self.host}:{self.port}"

    @property
    def local(self) -> bool:
        return self.host == "local"

    def ssh_opts(self) -> List[str]:
        CONTROL_DIR.mkdir(parents=True, exist_ok=True, mode=0o700)
        opts = ["-p", str(self.port), "-o", "BatchMode=yes", "-o", "StrictHostKeyChecking=accept-new",
                "-o", "ServerAliveInterval=15", "-o", "ServerAliveCountMax=4",
                "-o", "ControlMaster=auto", "-o", f"ControlPath={CONTROL_DIR}/%r@%h:%p", "-o", "ControlPersist=600"]
        if self.key_filename:
            opts += ["-i", os.path.expanduser(self.key_filename)]
        return opts

    def ssh(self, remote_cmd: str) -> List[str]:
        """argv running `remote_cmd` on the host (a login-less, non-interactive shell)."""
        if self.local:
            return ["bash", "-c", remote_cmd]
        return ["ssh", *self.ssh_opts(), f"{self.user}@{self.host}", remote_cmd]

    def rsync(self, srcs: List[str], dst: str, *extra: str) -> List[str]:
        """argv copying local `srcs` to `dst` on the host (rsync over the same ssh connection)."""
        if self.local:
            return ["rsync", "-a", *extra, *srcs, dst]
        ssh = " ".join(["ssh", *self.ssh_opts()])
        return ["rsync", "-a", "--partial", "-e", ssh, *extra, *srcs, f"{self.user}@{self.host}:{dst}"]


def load(path: Path = INVENTORY, enabled_only: bool = True) -> List[Host]:
    """Hosts from the inventory. An entry without `user` / `key_filename` takes them from the first entry that has
    them, so a new host can be just {"host": ..., "port": ...}."""
    if not path.exists():
        return []
    entries = json.loads(path.read_text())
    donor = next((e for e in entries if e.get("key_filename")), {})
    hosts = []
    for e in entries:
        e = {**{k: donor[k] for k in ("user", "key_filename") if k in donor}, **e}
        known = {k: e[k] for k in ("host", "user", "port", "key_filename", "workers", "enabled", "repo", "data_dir")
                 if k in e and e[k] is not None}
        h = Host(**known, extra={k: v for k, v in e.items() if k not in known})
        if h.enabled or not enabled_only:
            hosts.append(h)
    return hosts


def code_digest(root: Path = REPO) -> str:
    """Hash of the tracked-kind source files: main refuses a node whose code differs from its own."""
    h = hashlib.sha256()
    for d in CODE_DIRS:
        for p in sorted((root / d).rglob("*")):
            if p.is_file() and p.suffix in CODE_SUFFIXES and "build" not in p.parts and "__pycache__" not in p.parts:
                h.update(str(p.relative_to(root)).encode())
                h.update(p.read_bytes())
    return h.hexdigest()[:16]
