"""Accounts for a shared server install.

The desktop app has no logins. When PRIMERFORGE_AUTH is on (the Linux installer
sets it), every page and API call needs a signed-in user:

  admin   everything, including Settings, species installs and user accounts
  user    design primers, run BLAST/BLAT, browse transcripts; can delete only
          their own runs and tickets

Everyone sees everyone's runs and tickets, as a lab group sharing one server
would expect. Passwords are stored as salted hashes (Werkzeug's scrypt/pbkdf2),
a session ends when its password changes, and repeated failed sign-ins to one
account from one address are paused.
"""
from __future__ import annotations

import os
import re
import secrets
import threading
import time

from werkzeug.security import check_password_hash, generate_password_hash

from . import config, store

ROLES = ("admin", "user")
MIN_PASSWORD = 10
USERNAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._@-]{1,63}$")

SCHEMA = """
CREATE TABLE IF NOT EXISTS users (
    id            INTEGER PRIMARY KEY AUTOINCREMENT,
    username      TEXT NOT NULL UNIQUE COLLATE NOCASE,
    full_name     TEXT,
    password_hash TEXT NOT NULL,
    role          TEXT NOT NULL DEFAULT 'user',
    active        INTEGER NOT NULL DEFAULT 1,
    must_change   INTEGER NOT NULL DEFAULT 0,
    created_at    REAL NOT NULL,
    last_login    REAL
);
"""


class AuthError(ValueError):
    pass


def init() -> None:
    with store.db() as conn:
        conn.executescript(SCHEMA)


# ---------------------------------------------------------------------------
# Secret key for signing session cookies
# ---------------------------------------------------------------------------
def secret_key() -> str:
    """From the environment, else a random key kept in the data directory."""
    env = os.environ.get("PRIMERFORGE_SECRET_KEY", "").strip()
    if env:
        return env
    path = config.DATA_DIR / "secret_key"
    try:
        value = path.read_text(encoding="utf-8").strip()
        if len(value) >= 32:
            return value
    except OSError:
        pass
    config.DATA_DIR.mkdir(parents=True, exist_ok=True)
    value = secrets.token_hex(32)
    path.write_text(value, encoding="utf-8")
    try:
        os.chmod(path, 0o600)
    except OSError:
        pass
    return value


# ---------------------------------------------------------------------------
# Users
# ---------------------------------------------------------------------------
def generate_password() -> str:
    """A readable temporary password: four groups of letters and digits."""
    alphabet = "abcdefghjkmnpqrstuvwxyzABCDEFGHJKLMNPQRSTUVWXYZ23456789"
    return "-".join("".join(secrets.choice(alphabet) for _ in range(4)) for _ in range(4))


def _check_password(password: str) -> None:
    if len(password) < MIN_PASSWORD:
        raise AuthError(f"Passwords need at least {MIN_PASSWORD} characters.")


def _public(row) -> dict:
    user = dict(row)
    user.pop("password_hash", None)
    user["is_admin"] = user["role"] == "admin"
    return user


def list_users() -> list[dict]:
    init()
    with store.db() as conn:
        return [_public(r) for r in conn.execute(
            "SELECT * FROM users ORDER BY role != 'admin', username COLLATE NOCASE")]


def get_user(user_id: int) -> dict | None:
    with store.db() as conn:
        row = conn.execute("SELECT * FROM users WHERE id=?", (user_id,)).fetchone()
    return dict(row) if row else None


def find_user(username: str) -> dict | None:
    init()
    with store.db() as conn:
        row = conn.execute("SELECT * FROM users WHERE username=?", (username.strip(),)).fetchone()
    return dict(row) if row else None


def count_admins(exclude_id: int | None = None) -> int:
    with store.db() as conn:
        return conn.execute("SELECT COUNT(*) AS n FROM users WHERE role='admin' AND active=1 "
                            "AND id != ?", (exclude_id or -1,)).fetchone()["n"]


def create_user(username: str, password: str, role: str = "user", full_name: str = "",
                must_change: bool = False) -> dict:
    init()
    username = username.strip()
    if not USERNAME_RE.match(username):
        raise AuthError("Usernames are 2-64 characters: letters, digits, '.', '_', '-' or '@'.")
    if role not in ROLES:
        raise AuthError(f"Role must be one of: {', '.join(ROLES)}.")
    _check_password(password)
    if find_user(username):
        raise AuthError(f"There is already a user called '{username}'.")
    with store.db() as conn:
        conn.execute("INSERT INTO users(username, full_name, password_hash, role, must_change,"
                     " created_at) VALUES (?,?,?,?,?,?)",
                     (username, full_name.strip() or None, generate_password_hash(password),
                      role, int(must_change), time.time()))
    return _public(find_user(username))


def set_password(user_id: int, password: str, must_change: bool = False) -> None:
    _check_password(password)
    with store.db() as conn:
        conn.execute("UPDATE users SET password_hash=?, must_change=? WHERE id=?",
                     (generate_password_hash(password), int(must_change), user_id))


def set_role(user_id: int, role: str) -> None:
    if role not in ROLES:
        raise AuthError(f"Role must be one of: {', '.join(ROLES)}.")
    user = get_user(user_id)
    if not user:
        raise AuthError("No such user.")
    if user["role"] == "admin" and role != "admin" and count_admins(exclude_id=user_id) == 0:
        raise AuthError("Keep at least one administrator.")
    with store.db() as conn:
        conn.execute("UPDATE users SET role=? WHERE id=?", (role, user_id))


def set_active(user_id: int, active: bool) -> None:
    user = get_user(user_id)
    if not user:
        raise AuthError("No such user.")
    if not active and user["role"] == "admin" and count_admins(exclude_id=user_id) == 0:
        raise AuthError("Keep at least one active administrator.")
    with store.db() as conn:
        conn.execute("UPDATE users SET active=? WHERE id=?", (int(active), user_id))


def delete_user(user_id: int) -> None:
    user = get_user(user_id)
    if not user:
        raise AuthError("No such user.")
    if user["role"] == "admin" and count_admins(exclude_id=user_id) == 0:
        raise AuthError("Keep at least one administrator.")
    with store.db() as conn:
        conn.execute("DELETE FROM users WHERE id=?", (user_id,))


def session_token(user: dict) -> str:
    """Changes whenever the password does, which signs out other sessions."""
    return user["password_hash"][-16:]


# ---------------------------------------------------------------------------
# Signing in, with a brake on password guessing
# ---------------------------------------------------------------------------
FAIL_WINDOW = 15 * 60
FAIL_LIMIT = 8
_failures: dict[str, list[float]] = {}
_fail_lock = threading.Lock()


def locked_out(key: str) -> int:
    """Seconds until sign-in may be tried again for this address and username (0 = now)."""
    now = time.time()
    with _fail_lock:
        recent = [t for t in _failures.get(key, []) if now - t < FAIL_WINDOW]
        _failures[key] = recent
        if len(recent) < FAIL_LIMIT:
            return 0
        return int(FAIL_WINDOW - (now - recent[0])) + 1


def authenticate(username: str, password: str, address: str) -> dict:
    # Counted per address and username: colleagues behind one institutional NAT address
    # are not locked out by each other's typos, and guessing stays slow for each account.
    key = f"{address}|{(username or '').strip().lower()}"
    wait = locked_out(key)
    if wait:
        raise AuthError(f"Too many failed attempts. Try again in {wait // 60 + 1} minutes.")
    user = find_user(username or "")
    # Hash something even for unknown users so response time does not reveal them.
    ok = check_password_hash(user["password_hash"] if user else _DUMMY_HASH, password or "")
    if not user or not ok or not user["active"]:
        with _fail_lock:
            _failures.setdefault(key, []).append(time.time())
        if user and ok and not user["active"]:
            raise AuthError("This account has been disabled. Ask an administrator.")
        raise AuthError("Wrong username or password.")
    with _fail_lock:
        _failures.pop(key, None)
    with store.db() as conn:
        conn.execute("UPDATE users SET last_login=? WHERE id=?", (time.time(), user["id"]))
    return user


_DUMMY_HASH = generate_password_hash(secrets.token_hex(16))
