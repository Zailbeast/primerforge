"""Where BLAT runs: inside WSL on Windows, or directly on a Linux server.

UCSC only ships BLAT for Linux and macOS, so on Windows the BLAT binaries, their
2bit genomes and the gfServer indexes all live inside a WSL distribution and
everything here shells out through ``wsl.exe``; nothing assumes a particular
distro name, home directory or drive mount root. On Linux the same calls run the
programs directly ("native" mode), so callers never need to know which it is.
"""
from __future__ import annotations

import os
import platform
import re
import shlex
import shutil
import subprocess
import threading
from pathlib import Path, PureWindowsPath

from . import config

NO_WINDOW = getattr(subprocess, "CREATE_NO_WINDOW", 0)
# Distros that exist only to back other tools and have no userland worth using.
IGNORED_DISTROS = {"docker-desktop", "docker-desktop-data", "rancher-desktop",
                   "rancher-desktop-data"}


# On Linux there is no bridge: commands run as ordinary child processes.
NATIVE = os.name != "nt"


class WslError(RuntimeError):
    pass


_cache: dict[str, object] = {}
_lock = threading.Lock()


def _env() -> dict:
    env = dict(os.environ)
    env["WSL_UTF8"] = "1"          # otherwise wsl.exe's own messages are UTF-16
    return env


def exe() -> str | None:
    if os.name != "nt":
        return None
    found = shutil.which("wsl.exe") or shutil.which("wsl")
    if found:
        return found
    system = Path(os.environ.get("SystemRoot", r"C:\Windows")) / "System32" / "wsl.exe"
    return str(system) if system.exists() else None


def native_name() -> str:
    """A label for the Linux system BLAT runs on directly, e.g. 'Ubuntu 24.04.1 LTS'."""
    try:
        for line in Path("/etc/os-release").read_text(encoding="utf-8").splitlines():
            if line.startswith("PRETTY_NAME="):
                return line.split("=", 1)[1].strip().strip('"')
    except OSError:
        pass
    return f"{platform.system()} {platform.release()}".strip()


def distros() -> list[dict]:
    """Installed distributions as [{name, default, state, version}]."""
    if NATIVE:
        return []
    wsl = exe()
    if not wsl:
        return []
    try:
        proc = subprocess.run([wsl, "-l", "-v"], capture_output=True, timeout=30,
                              env=_env(), creationflags=NO_WINDOW)
    except (OSError, subprocess.TimeoutExpired):
        return []
    raw = proc.stdout
    text = raw.decode("utf-16-le", errors="ignore") if b"\x00" in raw else \
        raw.decode("utf-8", errors="ignore")
    out = []
    for line in text.splitlines()[1:]:
        line = line.replace("\x00", "").rstrip()
        if not line.strip():
            continue
        default = line.lstrip().startswith("*")
        parts = line.replace("*", " ").split()
        if len(parts) >= 3:
            out.append({"name": " ".join(parts[:-2]), "state": parts[-2],
                        "version": parts[-1], "default": default})
    return out


def distro() -> str | None:
    """The distribution BLAT runs in: the configured one, else the WSL default."""
    if NATIVE:
        return native_name()
    wanted = (config.settings().get("wsl_distro") or "").strip()
    names = distros()
    usable = [d for d in names if d["name"].lower() not in IGNORED_DISTROS]
    if wanted:
        return wanted if any(d["name"] == wanted for d in names) else None
    for d in usable:
        if d["default"]:
            return d["name"]
    return usable[0]["name"] if usable else None


def available() -> bool:
    return distro() is not None


def _base(distro_name: str | None = None) -> list[str]:
    if NATIVE:
        return []
    wsl = exe()
    name = distro_name or distro()
    if not wsl or not name:
        raise WslError("WSL is not installed or has no Linux distribution. "
                       "Install one with 'wsl --install -d Ubuntu'.")
    return [wsl, "-d", name]


def run(args: list[str], cwd: str | None = None, timeout: float | None = None,
        check: bool = False) -> subprocess.CompletedProcess:
    """Run a Linux program directly (no shell). Paths in args must be Linux paths."""
    if NATIVE:
        cmd, cwd_arg = list(map(str, args)), cwd
    else:
        cmd, cwd_arg = _base() + (["--cd", cwd] if cwd else []) + ["--exec", *map(str, args)], None
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, cwd=cwd_arg,
                          env=_env(), encoding="utf-8", errors="replace",
                          creationflags=NO_WINDOW)
    if check and proc.returncode != 0:
        raise WslError(f"{Path(str(args[0])).name} failed{'' if NATIVE else ' in WSL'}: "
                       f"{(proc.stderr or proc.stdout).strip()[:500]}")
    return proc


def bash(script: str, timeout: float | None = None, check: bool = False
         ) -> subprocess.CompletedProcess:
    return run(["bash", "-c", script], timeout=timeout, check=check)


def popen(args: list[str], cwd: str | None = None, **kw) -> subprocess.Popen:
    if NATIVE:
        return subprocess.Popen(list(map(str, args)), cwd=cwd, env=_env(), **kw)
    cmd = _base() + (["--cd", cwd] if cwd else []) + ["--exec", *map(str, args)]
    return subprocess.Popen(cmd, env=_env(), creationflags=NO_WINDOW, **kw)


def home() -> str:
    if NATIVE:
        # A server keeps BLAT with the rest of its data; the service account may have no $HOME.
        return str(config.DATA_DIR)
    with _lock:
        if "home" in _cache:
            return str(_cache["home"])
    proc = bash('printf %s "$HOME"', timeout=60, check=True)
    value = proc.stdout.strip()
    if not value.startswith("/"):
        raise WslError("Could not determine the WSL home directory.")
    with _lock:
        _cache["home"] = value
    return value


def linux_path(path: str) -> str:
    """Resolve a setting such as 'primerforge/blat' against the WSL home."""
    return path if path.startswith("/") else f"{home().rstrip('/')}/{path}"


def to_wsl(win_path: str | Path) -> str:
    """Translate a Windows path to the path WSL sees it at (/mnt/d/... by default)."""
    if NATIVE:
        return str(Path(win_path).resolve())
    p = PureWindowsPath(str(Path(win_path).resolve()))
    drive = p.drive.rstrip(":").lower()
    if not drive or len(drive) != 1:
        raise WslError(f"Cannot map '{win_path}' into WSL (not on a lettered drive).")
    with _lock:
        root = _cache.get("mount_root")
    if root is None:
        proc = run(["wslpath", "-u", "C:\\"], timeout=60)
        m = re.match(r"^(.*)/c/?$", proc.stdout.strip())
        root = m.group(1) if m else "/mnt"
        with _lock:
            _cache["mount_root"] = root
    rest = "/".join(p.parts[1:])
    return f"{root}/{drive}/{rest}"


def quote(value: str) -> str:
    return shlex.quote(str(value))


def reset_cache() -> None:
    with _lock:
        _cache.clear()
