"""The fleet inventory (~/.minos/fleet.json on main) and how main reaches each sub-VPS over SSH.

    [{"host": "1.2.3.4", "user": "root", "port": 22, "key_filename": "~/.ssh/id_ed25519",
      "workers": null, "enabled": true}]

- `workers`: null = the node's own rule, floor(cores x 0.95 - 1); a number overrides it.
- `port`: several machines may share one IP behind different ports (e.g. vast.ai).
- optional `repo` / `data_dir`: where the code and data live on that host (default /minos/minos-miner, /minos/data).
- `host: "local"`: run the node as a local subprocess instead of over SSH (tests, or this machine as a node).
- an entry without `user` / `key_filename` / `key_passphrase` inherits them from the first entry that has them.
- `key_passphrase`: main starts its own ssh-agent and unlocks the key with it; the passphrase is handed to ssh-add
  through an environment variable and never logged or written anywhere else. Main keeps this file at mode 600 itself.
The file is re-read at every round start, so hosts can be added or disabled between rounds.
"""
from __future__ import annotations

import atexit
import hashlib
import json
import os
import stat
import subprocess
import tempfile
import threading
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
    key_passphrase: Optional[str] = field(default=None, repr=False, compare=False)   # never printed
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
    if path.stat().st_mode & (stat.S_IRWXG | stat.S_IRWXO):
        try:
            os.chmod(path, 0o600)                               # keep it private; never a reason to stop mining
        except OSError:
            pass
    donor = next((e for e in entries if e.get("key_filename")), {})
    hosts = []
    for e in entries:
        e = {**{k: donor[k] for k in ("user", "key_filename", "key_passphrase") if k in donor}, **e}
        known = {k: e[k] for k in ("host", "user", "port", "key_filename", "key_passphrase", "workers", "enabled",
                                   "repo", "data_dir") if k in e and e[k] is not None}
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


# ---- ssh-agent for passphrase-protected keys --------------------------------------------------------------------
_AGENT_LOCK = threading.Lock()
_UNLOCKED: set = set()


def ensure_agent(hosts: List[Host]) -> None:
    """Start one ssh-agent for this process (its children inherit SSH_AUTH_SOCK) and add every key that has a
    passphrase. The passphrase goes to ssh-add through an environment variable read by a tiny askpass script; neither
    the script nor anything else on disk contains it. The agent is killed when this process exits."""
    need = {(os.path.expanduser(h.key_filename), h.key_passphrase) for h in hosts if h.key_filename and h.key_passphrase}
    need = {k for k in need if k[0] not in _UNLOCKED}
    if not need:
        return
    with _AGENT_LOCK:
        if "MINOS_AGENT_PID" not in os.environ:
            out = subprocess.run(["ssh-agent", "-s"], capture_output=True, text=True, check=True).stdout
            for line in out.splitlines():
                if line.startswith(("SSH_AUTH_SOCK=", "SSH_AGENT_PID=")):
                    k, _, v = line.partition(";")[0].partition("=")
                    os.environ[k] = v
            os.environ["MINOS_AGENT_PID"] = os.environ["SSH_AGENT_PID"]
            atexit.register(lambda pid=os.environ["SSH_AGENT_PID"]: subprocess.run(["kill", pid], capture_output=True))
        d = tempfile.mkdtemp(prefix="minos-askpass-")
        ask = os.path.join(d, "askpass.sh")
        with open(ask, "w") as f:
            f.write('#!/bin/sh\nprintf "%s" "$MINOS_KEY_PASSPHRASE"\n')
        os.chmod(ask, 0o700)
        try:
            for key, passphrase in need:
                env = {**os.environ, "SSH_ASKPASS": ask, "SSH_ASKPASS_REQUIRE": "force", "DISPLAY": ":0",
                       "MINOS_KEY_PASSPHRASE": passphrase}
                p = subprocess.run(["ssh-add", key], stdin=subprocess.DEVNULL, capture_output=True, text=True, env=env)
                if p.returncode != 0:
                    raise RuntimeError(f"ssh-add {key} failed (wrong passphrase?)")   # no passphrase in the message
                _UNLOCKED.add(key)
        finally:
            os.remove(ask)
            os.rmdir(d)

