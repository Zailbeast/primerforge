"""gfServer lifecycle inside WSL.

Ensembl answers BLAT searches from long-running gfServer processes that hold each
genome index in memory, with gfClient doing the per-query alignment. This module
does the same: one gfServer per species, started on first use from its
precomputed index (seconds rather than minutes), kept alive by the wsl.exe
process that launched it, and stopped cleanly when the app exits.
"""
from __future__ import annotations

import atexit
import subprocess
import threading
import time

from . import species, wsl

HOST = "127.0.0.1"


class BlatServerError(RuntimeError):
    pass


class _Server:
    def __init__(self, name: str):
        self.name = name
        self.proc: subprocess.Popen | None = None
        self.state = "stopped"
        self.message = ""
        self.started_at: float | None = None
        self.lock = threading.Lock()


_servers: dict[str, _Server] = {}
_registry_lock = threading.Lock()


def _server(name: str) -> _Server:
    with _registry_lock:
        if name not in _servers:
            _servers[name] = _Server(name)
        return _servers[name]


def _config(name: str) -> dict:
    m = species.load(name) or {}
    b = m.get("blat") or {}
    if not species.blat_installed(m):
        raise BlatServerError(f"BLAT is not installed for {name}. Install it from Settings.")
    return b


def ping(name: str, timeout: float = 25) -> dict | None:
    """gfServer status output as a dict, or None if nothing answers on the port."""
    b = _config(name)
    port = str(b["port"])
    # Without a listening server, `gfServer status` waits ten seconds for its connect
    # to time out, so first make sure a gfServer process for this port exists at all.
    # "gfServe[r]" matches the server's command line but not this bash -c line itself.
    script = (f"pgrep -f {wsl.quote('gfServe[r] start ' + HOST + ' ' + port)} >/dev/null && "
              f"{wsl.quote(b['bindir'] + '/gfServer')} status {HOST} {port}")
    try:
        proc = wsl.bash(script, timeout=timeout)
    except (subprocess.TimeoutExpired, wsl.WslError, OSError):
        return None
    if proc.returncode != 0:
        return None
    info = {}
    for line in proc.stdout.splitlines():          # "key value" lines
        key, _, value = line.strip().partition(" ")
        if key:
            info[key] = value.strip()
    return info if info.get("port") == port else None


def status(name: str) -> dict:
    srv = _server(name)
    m = species.load(name) or {}
    b = m.get("blat") or {}
    out = {"species": name, "installed": species.blat_installed(m), "port": b.get("port"),
           "distro": b.get("distro"), "dir": b.get("dir"), "step_size": b.get("step_size"),
           "state": srv.state, "message": srv.message, "started_at": srv.started_at}
    if not out["installed"]:
        out["state"] = "not installed"
        return out
    if srv.state in ("running", "stopped", "error"):
        info = ping(name, timeout=15)
        if info:
            out.update(state="running", info=info)
            if srv.state != "running":
                srv.state, srv.message = "running", "Answering on its port"
        elif srv.state == "running":
            srv.state, srv.message = "stopped", "No longer answering"
            out["state"] = "stopped"
    return out


def ensure_running(name: str, wait: float = 900) -> None:
    """Start the species' gfServer if needed and block until it answers."""
    srv = _server(name)
    with srv.lock:
        if srv.state == "running" and srv.proc and srv.proc.poll() is None:
            return
        if ping(name):
            srv.state, srv.message = "running", "Already running"
            return
        _start_locked(srv, wait)


def start(name: str) -> None:
    threading.Thread(target=_safe_start, args=(name,), daemon=True).start()


def _safe_start(name: str) -> None:
    try:
        ensure_running(name)
    except Exception as exc:                                     # noqa: BLE001
        srv = _server(name)
        srv.state, srv.message = "error", str(exc)


def _start_locked(srv: _Server, wait: float) -> None:
    b = _config(srv.name)
    srv.state, srv.message, srv.started_at = "starting", "Loading genome index", time.time()
    args = [f"{b['bindir']}/gfServer", "start", HOST, str(b["port"]),
            f"-stepSize={b.get('step_size') or 5}", "-canStop",
            f"-indexFile={b['index']}", f"-log={b['dir']}/gfServer.log", b["twobit"]]
    try:
        srv.proc = wsl.popen(args, cwd=b["dir"], stdin=subprocess.DEVNULL,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    except (OSError, wsl.WslError) as exc:
        srv.state, srv.message = "error", f"Could not launch gfServer in WSL: {exc}"
        raise BlatServerError(srv.message) from exc

    deadline = time.time() + wait
    while time.time() < deadline:
        if srv.proc.poll() is not None:
            output = (srv.proc.stdout.read() if srv.proc.stdout else "") or b""
            text = output.decode("utf-8", "replace") if isinstance(output, bytes) else output
            srv.state = "error"
            srv.message = f"gfServer exited: {text.strip()[-400:] or 'no output'}"
            raise BlatServerError(srv.message)
        if ping(srv.name, timeout=15):
            srv.state, srv.message = "running", f"Ready in {time.time() - srv.started_at:.0f}s"
            threading.Thread(target=_drain, args=(srv.proc,), daemon=True).start()
            return
        time.sleep(2)
    srv.state, srv.message = "error", "gfServer did not become ready in time"
    raise BlatServerError(srv.message)


def _drain(proc: subprocess.Popen) -> None:
    """Keep the pipe empty so a chatty gfServer can never block on a full buffer."""
    try:
        if proc.stdout:
            for _ in proc.stdout:
                pass
    except (OSError, ValueError):
        pass


def stop(name: str) -> None:
    srv = _server(name)
    with srv.lock:
        try:
            b = _config(name)
            wsl.run([f"{b['bindir']}/gfServer", "stop", HOST, str(b["port"])], timeout=30)
        except Exception:                                        # noqa: BLE001
            pass
        if srv.proc and srv.proc.poll() is None:
            try:
                srv.proc.wait(timeout=15)
            except subprocess.TimeoutExpired:
                srv.proc.kill()
        srv.proc = None
        srv.state, srv.message = "stopped", "Stopped"


def stop_all() -> None:
    for name, srv in list(_servers.items()):
        if srv.proc is not None:
            stop(name)


atexit.register(stop_all)
