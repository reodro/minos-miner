"""Run a `docker run` command inside the image's own files, with no container.

The unpacked rootfs is shared and read-only (overlay lower dir). Each job gets a private upper
dir, its own mount and pid and net namespace, and the cgroup limits Docker was given (--cpus,
--memory). A job can stay up: the warm JVM returns while Java is still serving requests.
`MINOS_RUNTIME=docker` leaves this unused.
"""
from __future__ import annotations

import json
import logging
import os
import shlex
import shutil
import signal
import subprocess
import threading
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

log = logging.getLogger(__name__)

_VALUE_FLAGS = {
    "--name", "--cpus", "--memory", "--memory-swap", "--network", "--entrypoint",
    "-v", "--volume", "-e", "--env", "-w", "--workdir", "-u", "--user",
}
_BOOL_FLAGS = {"--rm", "-d", "--detach", "--init", "--network"}


def use_sandbox() -> bool:
    return os.environ.get("MINOS_RUNTIME", "sandbox") != "docker"


@dataclass
class DockerSpec:
    image: str
    argv: List[str] = field(default_factory=list)   # arguments after the image name (they replace the image CMD)
    binds: List[Tuple[str, str, bool]] = field(default_factory=list)   # host, guest, read-only
    env: Dict[str, str] = field(default_factory=dict)
    name: str = ""
    cpus: Optional[float] = None
    memory: Optional[str] = None
    network_none: bool = False
    entrypoint: Optional[List[str]] = None


def parse_docker_run(argv: List[str]) -> DockerSpec:
    """The subset of `docker run` this evaluator builds. Anything else is refused."""
    if not argv or argv[0] != "docker" or "run" not in argv[:3]:
        raise ValueError("not a docker run command")
    rest = argv[argv.index("run") + 1:]
    spec = DockerSpec(image="")
    i = 0
    while i < len(rest):
        token = rest[i]
        if token == "--":
            i += 1
            break
        if token in ("--network",) or token.startswith("--network="):
            value = token.split("=", 1)[1] if "=" in token else rest[i + 1]
            if "=" not in token:
                i += 1
            spec.network_none = value == "none"
            i += 1
            continue
        parsed = _flag(rest, i)
        if parsed is None:
            break
        name, value = parsed
        i = value
        flag, arg = name
        if flag in ("--name",):
            spec.name = arg
        elif flag in ("--cpus",):
            spec.cpus = float(arg)
        elif flag in ("--memory", "--memory-swap"):
            if flag == "--memory":
                spec.memory = arg
        elif flag in ("-v", "--volume"):
            spec.binds.append(_bind(arg))
        elif flag in ("-e", "--env"):
            key, _, val = arg.partition("=")
            spec.env[key] = val
        elif flag == "--entrypoint":
            spec.entrypoint = arg.split()
        i += 1
    if i >= len(rest):
        raise ValueError("docker run has no image")
    spec.image = rest[i]
    spec.argv = rest[i + 1:]
    return spec


def _flag(argv: List[str], i: int):
    token = argv[i]
    if not token.startswith("-"):
        return None
    for flag in _VALUE_FLAGS:
        if token == flag:
            if i + 1 >= len(argv):
                raise ValueError(f"{flag} needs a value")
            return (flag, argv[i + 1]), i + 1
        if token.startswith(flag + "="):
            return (flag, token.split("=", 1)[1]), i
    if token in _BOOL_FLAGS or token.startswith("--"):
        return (token, ""), i
    raise ValueError(f"unsupported docker flag {token}")


def _memory_bytes(text: str) -> int:
    """Docker's 4g / 512m into the integer systemd's MemoryMax accepts."""
    raw = text.strip().lower()
    mult = {"k": 1024, "m": 1024 ** 2, "g": 1024 ** 3, "b": 1}
    if raw[-1] in mult and not raw[:-1].isdigit():
        return int(float(raw[:-1]) * mult[raw[-1]])
    if raw[-1] in mult:
        return int(raw[:-1]) * mult[raw[-1]]
    return int(raw)


def _bind(arg: str) -> Tuple[str, str, bool]:
    ro = arg.endswith(":ro")
    body = arg[:-3] if ro else arg
    host, guest = body.rsplit(":", 1)
    return host, guest, ro


def command_for(spec: DockerSpec, image_config: dict) -> List[str]:
    """What Docker would execute: entrypoint, then the args (or the image CMD when there are none)."""
    entry = spec.entrypoint if spec.entrypoint is not None else list(image_config.get("Entrypoint") or [])
    tail = spec.argv if spec.argv else list(image_config.get("Cmd") or [])
    return entry + tail


def image_env(image_config: dict, extra: Dict[str, str]) -> Dict[str, str]:
    env = {}
    for item in image_config.get("Env") or []:
        key, _, val = item.partition("=")
        env[key] = val
    env.update(extra)
    return env


# ---- rootfs ---------------------------------------------------------------------------------------------------------

def _cache_dir() -> Path:
    from common.settings import SETTINGS
    return Path(SETTINGS.cache_dir) / "rootfs"


def _digest_of(image: str) -> Optional[str]:
    """The directory name of an image pinned as name@sha256:hex. A tag has no digest in the string."""
    marker = "@sha256:"
    if marker not in image:
        return None
    digest = image.split(marker, 1)[1]
    if len(digest) >= 12 and all(c in "0123456789abcdef" for c in digest):
        return digest
    return None


def _names_path() -> Path:
    return _cache_dir() / "names.json"


def _load_names() -> Dict[str, str]:
    try:
        data = json.loads(_names_path().read_text())
    except (OSError, json.JSONDecodeError):
        return {}
    return data if isinstance(data, dict) else {}


def _remember(image: str, ident: str) -> None:
    """Record which unpacked directory `image` uses, so the next job does not ask Docker."""
    path = _names_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    names = _load_names()
    if names.get(image) == ident:
        return
    names[image] = ident
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(names))
    os.replace(tmp, path)


def _cached(image: str) -> Optional[Tuple[Path, dict]]:
    """The unpacked image, when this machine already has it. No Docker."""
    idents = []
    digest = _digest_of(image)
    if digest:
        idents.append(digest)
    mapped = _load_names().get(image)
    if mapped and mapped not in idents:
        idents.append(mapped)
    for ident in idents:
        dest = _cache_dir() / ident
        config_path = dest / "config.json"
        if (dest / ".ready").exists() and config_path.exists():
            return dest, json.loads(config_path.read_text())
    return None


def _image_id(image: str) -> str:
    out = subprocess.run(["docker", "image", "inspect", "--format", "{{.Id}}", image],
                         capture_output=True, text=True, timeout=60)
    if out.returncode != 0:
        raise RuntimeError(f"image {image} is not present ({out.stderr.strip()})")
    return out.stdout.strip().split(":")[-1]


def _image_config(image: str) -> dict:
    out = subprocess.run(["docker", "image", "inspect", "--format", "{{json .Config}}", image],
                         capture_output=True, text=True, timeout=60)
    if out.returncode != 0:
        raise RuntimeError(out.stderr.strip())
    cfg = json.loads(out.stdout)
    return {k: cfg.get(k) for k in ("Env", "Entrypoint", "Cmd", "WorkingDir")}


def ensure_rootfs(image: str) -> Tuple[Path, dict]:
    """Unpack the image once into the cache. Later jobs mount it read-only and do not call Docker."""
    hit = _cached(image)
    if hit:
        return hit
    cache = _cache_dir()
    cache.mkdir(parents=True, exist_ok=True)
    lock_path = cache / ".lock"
    lock_path.touch()
    import fcntl
    with open(lock_path) as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        hit = _cached(image)
        if hit:
            return hit
        ident = _image_id(image)
        dest = cache / ident
        ready = dest / ".ready"
        if ready.exists():
            _remember(image, ident)
            return dest, json.loads((dest / "config.json").read_text())
        staging = dest.parent / f".partial-{ident[:12]}-{uuid.uuid4().hex[:6]}"
        if staging.exists():
            subprocess.run(["rm", "-rf", str(staging)], check=False)
        staging.mkdir(parents=True)
        name = f"minos-unpack-{uuid.uuid4().hex[:10]}"
        try:
            created = subprocess.run(["docker", "create", "--name", name, image],
                                     capture_output=True, text=True, timeout=120)
            if created.returncode != 0:
                raise RuntimeError(created.stderr.strip())
            extract = subprocess.run(f"docker export {name} | tar -C {shlex.quote(str(staging))} -xf -",
                                     shell=True, timeout=1800)
            if extract.returncode != 0:
                raise RuntimeError(f"exporting {image} failed")
            (staging / "config.json").write_text(json.dumps(_image_config(image)))
            (staging / ".ready").write_text(ident + "\n")
            if dest.exists():
                subprocess.run(["rm", "-rf", str(dest)], check=False)
            staging.rename(dest)
        except Exception:
            subprocess.run(["rm", "-rf", str(staging)], check=False)
            raise
        finally:
            subprocess.run(["docker", "rm", name], capture_output=True)
        _remember(image, ident)
        return dest, json.loads((dest / "config.json").read_text())


# ---- running --------------------------------------------------------------------------------------------------------

class _Job:
    def __init__(self, proc: subprocess.Popen, unit: str, name: str, script: Path, dirs: List[str]):
        self.proc, self.unit, self.name, self.script, self.dirs = proc, unit, name, script, dirs

    def kill(self) -> None:
        if self.unit:
            subprocess.run(["systemctl", "kill", "-s", "KILL", "--kill-who=all", self.unit], capture_output=True)
        try:
            os.killpg(self.proc.pid, signal.SIGKILL)
        except (ProcessLookupError, PermissionError):
            pass
        try:
            self.proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            self.proc.kill()
        self.script.unlink(missing_ok=True)
        for path in self.dirs:
            subprocess.run(["rm", "-rf", path], capture_output=True)


_jobs: Dict[str, _Job] = {}
_jobs_lock = threading.Lock()


def _chroot_bin() -> str:
    """Absolute path of chroot. The script runs it under `env -i`, which clears PATH, so a later lookup cannot work.
    Ubuntu 22.04 installs it at /usr/sbin/chroot; newer Ubuntu installs the same program at /usr/bin/chroot."""
    for candidate in ("/usr/bin/chroot", "/usr/sbin/chroot"):
        if os.access(candidate, os.X_OK):
            return candidate
    found = shutil.which("chroot")
    if found:
        return found
    raise RuntimeError("chroot not found")


def _script(rootfs: Path, spec: DockerSpec, config: dict, command: List[str],
            mnt: str, upper: str, work: str) -> str:
    lines = [
        "set -e",
        f"mkdir -p {shlex.quote(mnt)} {shlex.quote(upper)} {shlex.quote(work)}",
        "mount --make-rprivate /",
        (f"mount -t overlay overlay -o lowerdir={shlex.quote(str(rootfs))},"
         f"upperdir={shlex.quote(upper)},workdir={shlex.quote(work)} {shlex.quote(mnt)}"),
    ]
    for host, guest, ro in spec.binds:
        target = mnt + guest
        lines.append(f"mkdir -p {shlex.quote(target)}")
        lines.append(f"mount --bind {shlex.quote(host)} {shlex.quote(target)}")
        if ro:
            lines.append(f"mount -o remount,bind,ro {shlex.quote(target)}")
    lines.append(f"mkdir -p {shlex.quote(mnt + '/proc')} {shlex.quote(mnt + '/dev')} {shlex.quote(mnt + '/tmp')}")
    lines.append(f"mount -t proc proc {shlex.quote(mnt + '/proc')}")
    lines.append(f"mount --rbind /dev {shlex.quote(mnt + '/dev')}")
    lines.append(f"mount -t tmpfs tmpfs {shlex.quote(mnt + '/tmp')}")
    lines.append("ulimit -n 1048576 || ulimit -n 65536 || true")
    lines.append("ip link set lo up >/dev/null 2>&1 || true")
    wd = config.get("WorkingDir") or "/"
    lines.append(f"if [ -d {shlex.quote(mnt + wd)} ]; then cd {shlex.quote(mnt + wd)}; else cd {shlex.quote(mnt)}; fi")
    env = image_env(config, spec.env)
    assigns = " ".join(f"{shlex.quote(k)}={shlex.quote(v)}" for k, v in env.items())
    quoted = " ".join(shlex.quote(a) for a in command)
    lines.append(f"exec /usr/bin/env -i {assigns} {shlex.quote(_chroot_bin())} {shlex.quote(mnt)} {quoted}")
    return "\n".join(lines) + "\n"


def _launcher(spec: DockerSpec, script_path: Path, unit: str) -> List[str]:
    unshare = ["unshare", "--mount", "--pid", "--fork", "--kill-child", "--propagation", "private"]
    if spec.network_none:
        unshare.append("--net")
    unshare += ["sh", str(script_path)]
    limits = []
    if spec.cpus:
        limits += ["-p", f"CPUQuota={max(1, int(round(spec.cpus * 100)))}%"]
    if spec.memory:
        limits += ["-p", f"MemoryMax={_memory_bytes(spec.memory)}"]
    if not limits:
        return unshare
    return ["systemd-run", "--scope", f"--unit={unit}", *limits, "--", *unshare]


def start(argv: List[str], *, stdout=None) -> subprocess.Popen:
    """Start the command. The Popen stays alive until the image process exits (a warm JVM stays up)."""
    spec = parse_docker_run(list(argv))
    rootfs, config = ensure_rootfs(spec.image)
    return start_in(rootfs, config, spec, stdout=stdout)


def start_in(rootfs: Path, config: dict, spec: DockerSpec, *, stdout=None) -> subprocess.Popen:
    """Start `spec` inside an already-unpacked rootfs. Tests use this with a tiny root."""
    command = command_for(spec, config)
    if not command:
        raise RuntimeError(f"{spec.image} has no command")
    token = uuid.uuid4().hex[:12]
    script_path = Path(f"/tmp/minos-sandbox-{token}.sh")
    dirs = [f"/tmp/minos-mnt-{token}", f"/tmp/minos-upper-{token}", f"/tmp/minos-work-{token}"]
    script_path.write_text(_script(rootfs, spec, config, command, *dirs))
    unit = f"minos-{token}"
    proc = subprocess.Popen(_launcher(spec, script_path, unit), stdout=stdout, stderr=subprocess.STDOUT,
                            start_new_session=True)
    scope = f"{unit}.scope" if spec.cpus or spec.memory else ""
    job = _Job(proc, scope, spec.name, script_path, dirs)
    proc.minos_kill = job.kill          # type: ignore[attr-defined]
    with _jobs_lock:
        if spec.name:
            _jobs[spec.name] = job
        _jobs[f"pid:{proc.pid}"] = job
    return proc


def run(argv: List[str], *, timeout: Optional[float] = None, capture_output: bool = True,
        text: bool = True, **_ignored) -> subprocess.CompletedProcess:
    """Run to completion. A timeout kills the scope, not just the client."""
    proc = start(argv, stdout=subprocess.PIPE if capture_output else None)
    try:
        out, err = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        proc.minos_kill()               # type: ignore[attr-defined]
        raise
    finally:
        release(proc)
    if text and isinstance(out, bytes):
        out = out.decode("utf-8", "replace")
    # GATK and hap.py write their errors to stderr. The scope merges both streams; callers read either one.
    body = out or ""
    return subprocess.CompletedProcess(argv, proc.returncode, body, body)


def kill_name(name: str) -> None:
    with _jobs_lock:
        job = _jobs.get(name)
    if job is not None:
        job.kill()
        _drop(job.proc)


def kill_prefix(prefix: str) -> int:
    """Kill every sandbox job whose docker --name starts with prefix. Returns how many were still running."""
    with _jobs_lock:
        jobs = [j for n, j in _jobs.items() if not n.startswith("pid:") and n.startswith(prefix)]
    killed = 0
    for job in jobs:
        if job.proc.poll() is None:
            killed += 1
        job.kill()
        _drop(job.proc)
    return killed


def _drop(proc: subprocess.Popen) -> None:
    with _jobs_lock:
        for key in [k for k, j in _jobs.items() if j.proc is proc]:
            _jobs.pop(key, None)


def release(proc: subprocess.Popen) -> None:
    """Drop a finished job and remove its private overlay directories."""
    with _jobs_lock:
        jobs = [j for j in _jobs.values() if j.proc is proc]
    seen = set()
    for job in jobs:
        if id(job) in seen:
            continue
        seen.add(id(job))
        job.script.unlink(missing_ok=True)
        for path in job.dirs:
            subprocess.run(["rm", "-rf", path], capture_output=True)
    _drop(proc)


def score_cpus() -> Optional[float]:
    """CPU limit for one hap.py scoring (MINOS_SCORE_CPUS, default 1; 0 = none). Upstream runs hap.py with
    --threads 4 and no limit: on a big node with one scorer per worker, bursts of scorings then oversubscribed the CPU
    ~5x (load ~1000 on 180 cores). At 1 core a full scoring takes ~10 s instead of ~5 s idle; results are identical."""
    raw = (os.environ.get("MINOS_SCORE_CPUS") or "").strip()
    value = float(raw) if raw else 1.0
    return value if value > 0 else None


def score_threads() -> Optional[int]:
    """hap.py's own thread count for one scoring (MINOS_SCORE_THREADS, default 1; 0 = upstream's --threads 4).
    Under the 1-CPU limit four threads only contend: measured 2026-10-05, a cold full scoring 13.6 s -> 10.5 s with
    one thread, identical results."""
    raw = (os.environ.get("MINOS_SCORE_THREADS") or "").strip()
    value = int(raw) if raw else 1
    return value if value > 0 else None


def with_score_limit(args: List[str]) -> List[str]:
    """`docker run ... <hap.py image>` with the scoring CPU limit added (unless it already has one), and hap.py's
    --threads replaced by MINOS_SCORE_THREADS when set."""
    cpus = score_cpus()
    if "run" not in args[:3] or not any("hap-py" in a or "hap.py" in a for a in args):
        return args
    args = list(args)
    threads = score_threads()
    if threads is not None:
        for j in range(len(args) - 1):
            if args[j] == "--threads":
                args[j + 1] = str(threads)
    if cpus is None or any(a == "--cpus" or a.startswith("--cpus=") for a in args):
        return args
    i = args.index("run") + 1
    return args[:i] + [f"--cpus={cpus:g}"] + args[i:]


class _ScoringSubprocess:
    """Stands in for utils.scoring.subprocess. docker run goes through the sandbox (hap.py with the scoring CPU
    limit); everything else is unchanged."""

    _minos_sandbox = True

    def run(self, args, **kwargs):
        if isinstance(args, (list, tuple)) and args and args[0] == "docker":
            return run(with_score_limit(list(args)), timeout=kwargs.get("timeout"),
                       capture_output=kwargs.get("capture_output", False), text=kwargs.get("text", False))
        return subprocess.run(args, **kwargs)

    def __getattr__(self, name):
        return getattr(subprocess, name)


def install_scoring_shim() -> None:
    """Point upstream scoring's subprocess at the sandbox. Call it in every process that scores."""
    if not use_sandbox():
        return
    import utils.scoring as scoring
    if getattr(scoring.subprocess, "_minos_sandbox", False):
        return
    scoring.subprocess = _ScoringSubprocess()
