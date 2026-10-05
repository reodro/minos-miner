"""The dashboard's one account: a scrypt password hash, signed session cookies, and a lockout after failed logins.

The account lives in AUTH_FILE (~/.minos/dashboard/auth.json, mode 600, outside the repo; MINOS_DASHBOARD_AUTH
overrides the path), written by `python -m dashboard passwd`. Nothing secret is in the repo. Setting a new password
also makes a new signing key, so every open session ends.
"""
from __future__ import annotations

import base64
import hashlib
import hmac
import json
import os
import secrets
import threading
import time
from pathlib import Path
from typing import Dict, List, Optional

from itsdangerous import BadSignature, URLSafeTimedSerializer

AUTH_FILE = Path(os.environ.get("MINOS_DASHBOARD_AUTH") or Path.home() / ".minos" / "dashboard" / "auth.json")
SESSION_S = 12 * 3600
MAX_FAILS = 5                       # failed logins from one address within LOCK_S lock it out for LOCK_S
LOCK_S = 15 * 60
SCRYPT = {"n": 2 ** 15, "r": 8, "p": 1}


def _b64(b: bytes) -> str:
    return base64.b64encode(b).decode()


def hash_password(password: str, salt: Optional[bytes] = None, n: int = SCRYPT["n"], r: int = SCRYPT["r"],
                  p: int = SCRYPT["p"]) -> str:
    salt = salt or secrets.token_bytes(16)
    h = hashlib.scrypt(password.encode(), salt=salt, n=n, r=r, p=p, maxmem=128 * 1024 * 1024, dklen=32)
    return f"scrypt${n}${r}${p}${_b64(salt)}${_b64(h)}"


def check_password(password: str, stored: str) -> bool:
    try:
        kind, n, r, p, salt, h = stored.split("$")
        if kind != "scrypt":
            return False
        got = hash_password(password, base64.b64decode(salt), int(n), int(r), int(p)).split("$")[-1]
        return hmac.compare_digest(got, h)
    except (ValueError, TypeError):
        return False


def set_password(user: str, password: str, path: Path = AUTH_FILE) -> None:
    if not password:
        raise ValueError("the password cannot be empty")
    path.parent.mkdir(parents=True, exist_ok=True)
    os.chmod(path.parent, 0o700)
    data = {"user": user, "hash": hash_password(password), "secret": secrets.token_hex(32), "set_at": time.time()}
    tmp = path.with_name(f".{path.name}.tmp")
    fd = os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    with os.fdopen(fd, "w") as f:
        json.dump(data, f)
    os.replace(tmp, path)


class Auth:
    """Reads the account file on each check (a new password applies without a restart)."""

    def __init__(self, path: Path = AUTH_FILE):
        self.path = Path(path)
        self._fails: Dict[str, List[float]] = {}
        self._lock = threading.Lock()

    def account(self) -> Optional[dict]:
        try:
            data = json.loads(self.path.read_text())
            return data if {"user", "hash", "secret"} <= set(data) else None
        except (OSError, ValueError):
            return None

    # ---- lockout --------------------------------------------------------------------------------------------------
    def locked_for(self, addr: str) -> float:
        """Seconds this address stays locked out (0: not locked)."""
        now = time.time()
        with self._lock:
            fails = [t for t in self._fails.get(addr, []) if now - t < LOCK_S]
            self._fails[addr] = fails
            if len(fails) >= MAX_FAILS:
                return LOCK_S - (now - fails[-MAX_FAILS])
            return 0.0

    def _failed(self, addr: str) -> None:
        with self._lock:
            self._fails.setdefault(addr, []).append(time.time())

    # ---- login and sessions -------------------------------------------------------------------------------------------
    def login(self, addr: str, user: str, password: str) -> Optional[str]:
        """A session token, or None. The caller checks `locked_for` first."""
        acc = self.account()
        ok = acc is not None and hmac.compare_digest(user.encode(), acc["user"].encode()) and \
            check_password(password, acc["hash"])
        if not ok:
            self._failed(addr)
            return None
        with self._lock:
            self._fails.pop(addr, None)
        return self._serializer(acc).dumps({"u": acc["user"]})

    def user(self, token: Optional[str]) -> Optional[str]:
        acc = self.account()
        if not token or acc is None:
            return None
        try:
            data = self._serializer(acc).loads(token, max_age=SESSION_S)
        except BadSignature:
            return None
        return data.get("u") if data.get("u") == acc["user"] else None

    @staticmethod
    def _serializer(acc: dict) -> URLSafeTimedSerializer:
        return URLSafeTimedSerializer(acc["secret"], salt="minos-dashboard-session")
