"""Fill this folder with everything the Linux server installation needs.

Run on any computer with internet access (Windows, macOS or Linux, Python 3.9+):

    python installation/prepare.py              # app + dependencies
    python installation/prepare.py --archive    # also pack the folder into one .tar.gz

Afterwards the installation folder holds:

    install.sh, uninstall.sh, primerforge.sh, services/   the installer
    app/                    the PrimerForge application (copied from this project)
    dependencies/python/    Python packages for Python 3.9-3.14 on x86_64 Linux
    dependencies/blast/     NCBI BLAST+ for x86_64 Linux
    dependencies/blat/      UCSC BLAT programs (gfServer, gfClient, blat, faToTwoBit, twoBitInfo)
    dependencies/caddy/     Caddy web server, used for HTTPS
    dependencies/SHA256SUMS checked by install.sh, so a damaged copy is caught

Copy the whole folder to the server and run:  sudo bash install.sh

Downloads already present are kept, so running this again is quick; --refresh
fetches everything anew. Re-run it after changing the application so app/ is current.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
import tarfile
import time
import urllib.error
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
APP = HERE / "app"
DEPS = HERE / "dependencies"
USER_AGENT = "PrimerForge-installer/1.0"

APP_ITEMS = ["app.py", "serve.py", "requirements.txt", "README.md", "run.bat",
             "primerforge", "templates", "static", "tools"]
SKIP_DIRS = {"__pycache__", ".git"}
SKIP_SUFFIXES = {".pyc", ".pyo", ".log"}
TEXT_SUFFIXES = {".sh", ".service", ".py", ".html", ".js", ".css", ".txt", ".md", ".bat"}
INSTALLER_FILES = ["install.sh", "uninstall.sh", "primerforge.sh", "services"]

PYTHON_VERSIONS = ["3.9", "3.10", "3.11", "3.12", "3.13", "3.14"]
LINUX_PLATFORMS = ["manylinux2014_x86_64", "manylinux_2_17_x86_64", "manylinux_2_28_x86_64",
                   "manylinux_2_34_x86_64"]
# pip evaluates environment markers against the Python running this script, not the
# target, so dependencies that only older Pythons need are listed here explicitly.
EXTRA_FOR_OLD_PYTHON = {"3.9": ["importlib-metadata>=3.6", "zipp>=3.20"]}

BLAST_LATEST = "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
UCSC_MIRRORS = ("https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64",
                "https://hgdownload-euro.soe.ucsc.edu/admin/exe/linux.x86_64")
BLAT_BINARIES = ("blat/gfServer", "blat/gfClient", "blat/blat", "faToTwoBit", "twoBitInfo")
CADDY_FALLBACK = "2.11.4"


def say(message: str) -> None:
    print(message, flush=True)


def human(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} GB"


# ---------------------------------------------------------------------------
# Network
# ---------------------------------------------------------------------------
def fetch_text(url: str, timeout: int = 60) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", errors="replace")


def download(url: str, dest: Path, refresh: bool = False) -> Path:
    """Resumable download with retries; an existing complete file is kept."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and not refresh:
        say(f"  kept     {dest.relative_to(HERE)} ({human(dest.stat().st_size)})")
        return dest
    part = dest.with_name(dest.name + ".part")
    if refresh:
        part.unlink(missing_ok=True)
    last_error = None
    for attempt in range(5):
        have = part.stat().st_size if part.exists() else 0
        req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        if have:
            req.add_header("Range", f"bytes={have}-")
        try:
            with urllib.request.urlopen(req, timeout=120) as resp:
                if have and resp.status != 206:
                    have = 0
                total = int(resp.headers.get("Content-Length", 0)) + have
                got, shown = have, 0.0
                with open(part, "ab" if have else "wb") as out:
                    while True:
                        chunk = resp.read(1 << 20)
                        if not chunk:
                            break
                        out.write(chunk)
                        got += len(chunk)
                        if time.time() - shown > 3 and total:
                            shown = time.time()
                            print(f"\r  download {dest.name}: {human(got)} of {human(total)}   ",
                                  end="", flush=True)
            if total and part.stat().st_size < total:
                raise IOError("download ended early")
            part.replace(dest)
            print("\r", end="")
            say(f"  fetched  {dest.relative_to(HERE)} ({human(dest.stat().st_size)})")
            return dest
        except urllib.error.HTTPError as exc:
            if exc.code == 416 and part.exists():
                part.replace(dest)
                return dest
            if exc.code in (403, 404):
                raise
            last_error = exc
        except Exception as exc:                               # noqa: BLE001
            last_error = exc
        time.sleep(3 * (attempt + 1))
    raise RuntimeError(f"could not download {url}: {last_error}")


# ---------------------------------------------------------------------------
# Steps
# ---------------------------------------------------------------------------
def copy_app() -> None:
    say("\nApplication")
    if APP.exists():
        shutil.rmtree(APP)
    count = 0
    for item in APP_ITEMS:
        src = PROJECT / item
        paths = [src] if src.is_file() else sorted(p for p in src.rglob("*") if p.is_file())
        for path in paths:
            rel = path.relative_to(PROJECT)
            if SKIP_DIRS & set(rel.parts) or path.suffix in SKIP_SUFFIXES:
                continue
            target = APP / rel
            target.parent.mkdir(parents=True, exist_ok=True)
            data = path.read_bytes()
            if path.suffix in TEXT_SUFFIXES and path.suffix != ".bat":
                data = data.replace(b"\r\n", b"\n")        # Linux line endings
            target.write_bytes(data)
            count += 1
    version = re.search(r'__version__ = "([^"]+)"',
                        (APP / "primerforge" / "__init__.py").read_text(encoding="utf-8")).group(1)
    say(f"  copied   app/ ({count} files, PrimerForge {version})")


def normalise_installer() -> None:
    """Shell scripts and unit files must have Unix line endings to run on the server."""
    for name in INSTALLER_FILES:
        path = HERE / name
        for p in ([path] if path.is_file() else path.rglob("*")):
            if p.is_file():
                data = p.read_bytes()
                if b"\r\n" in data:
                    p.write_bytes(data.replace(b"\r\n", b"\n"))


def python_packages(refresh: bool) -> None:
    say("\nPython packages (for Python " + ", ".join(PYTHON_VERSIONS) + " on x86_64 Linux)")
    wheels = DEPS / "python"
    if refresh and wheels.exists():
        shutil.rmtree(wheels)
    wheels.mkdir(parents=True, exist_ok=True)
    requirements = PROJECT / "requirements.txt"
    for version in PYTHON_VERSIONS:
        cmd = [sys.executable, "-m", "pip", "download", "--quiet", "--disable-pip-version-check",
               "--dest", str(wheels), "--only-binary=:all:", "--implementation", "cp",
               "--python-version", version]
        for platform in LINUX_PLATFORMS:
            cmd += ["--platform", platform]
        cmd += ["-r", str(requirements), *EXTRA_FOR_OLD_PYTHON.get(version, [])]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"pip could not download packages for Python {version}:\n"
                               f"{(proc.stderr or proc.stdout).strip()[-1500:]}")
        say(f"  ready    Python {version}")
    total = sum(p.stat().st_size for p in wheels.glob("*.whl"))
    say(f"  {len(list(wheels.glob('*.whl')))} packages in dependencies/python ({human(total)})")


def blast(refresh: bool) -> str:
    say("\nNCBI BLAST+")
    listing = fetch_text(BLAST_LATEST)
    names = sorted(set(re.findall(r'href="(ncbi-blast-[\d.]+\+-x64-linux\.tar\.gz)"', listing)))
    if not names:
        raise RuntimeError(f"no x64 Linux BLAST+ package listed at {BLAST_LATEST}")
    name = names[-1]
    folder = DEPS / "blast"
    for old in folder.glob("ncbi-blast-*.tar.gz"):
        if old.name != name:
            old.unlink()                        # an older release than NCBI now offers
    dest = download(BLAST_LATEST + name, folder / name, refresh)
    expected = fetch_text(BLAST_LATEST + name + ".md5").split()[0].lower()
    actual = hashlib.md5(dest.read_bytes()).hexdigest()
    if actual != expected:
        dest.unlink()
        raise RuntimeError(f"{name} failed NCBI's MD5 check; run prepare.py again")
    say("  verified NCBI MD5 checksum")
    return re.search(r"ncbi-blast-([\d.]+)\+", name).group(1)


def reachable(url: str) -> bool:
    """A quick check, so an unreachable mirror costs seconds rather than a round of retries."""
    try:
        req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(req, timeout=20):
            return True
    except Exception:                                          # noqa: BLE001
        return False


def blat(refresh: bool) -> None:
    say("\nUCSC BLAT")
    folder = DEPS / "blat"
    mirrors = list(UCSC_MIRRORS)
    for rel in BLAT_BINARIES:
        name = rel.split("/")[-1]
        errors = []
        for mirror in list(mirrors):
            if not (folder / name).exists() or refresh:
                if not reachable(f"{mirror}/{rel}"):
                    errors.append(f"{mirror} unreachable")
                    continue
            try:
                download(f"{mirror}/{rel}", folder / name, refresh)
                mirrors.remove(mirror)
                mirrors.insert(0, mirror)               # try the working mirror first next time
                break
            except Exception as exc:                           # noqa: BLE001
                errors.append(f"{mirror}: {exc}")
        else:
            raise RuntimeError(f"could not download {name}: " + "; ".join(errors))
        if (folder / name).stat().st_size < 100_000:
            (folder / name).unlink()
            raise RuntimeError(f"{name} from UCSC is too small to be the program; try again")


def caddy(refresh: bool) -> str:
    say("\nCaddy (HTTPS)")
    try:
        release = json.loads(fetch_text("https://api.github.com/repos/caddyserver/caddy/releases/latest",
                                        timeout=30))
        version = release["tag_name"].lstrip("v")
    except Exception:                                          # noqa: BLE001
        version = CADDY_FALLBACK
    base = f"https://github.com/caddyserver/caddy/releases/download/v{version}/"
    name = f"caddy_{version}_linux_amd64.tar.gz"
    folder = DEPS / "caddy"
    for old in folder.glob("caddy_*_linux_amd64.tar.gz"):
        if old.name != name:
            old.unlink()
    dest = download(base + name, folder / name, refresh)
    try:
        sums = fetch_text(base + f"caddy_{version}_checksums.txt")
        expected = next(line.split()[0] for line in sums.splitlines() if line.strip().endswith(name))
        if hashlib.sha512(dest.read_bytes()).hexdigest() != expected.lower():
            dest.unlink()
            raise RuntimeError(f"{name} failed Caddy's checksum; run prepare.py again")
        say("  verified Caddy SHA-512 checksum")
    except StopIteration:
        say("  (Caddy's checksum list did not include this file; not verified)")
    return version


def write_checksums(versions: dict) -> None:
    lines = []
    for path in sorted(p for p in DEPS.rglob("*") if p.is_file()
                       and p.name not in ("SHA256SUMS", "VERSIONS.txt") and not p.name.endswith(".part")):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        lines.append(f"{digest}  {path.relative_to(DEPS).as_posix()}")
    (DEPS / "SHA256SUMS").write_bytes(("\n".join(lines) + "\n").encode())
    info = [f"Prepared {time.strftime('%Y-%m-%d %H:%M')} on {sys.platform} with Python "
            f"{sys.version.split()[0]}"] + [f"{k}: {v}" for k, v in versions.items()]
    (DEPS / "VERSIONS.txt").write_bytes(("\n".join(info) + "\n").encode())


def archive() -> Path:
    version = re.search(r'__version__ = "([^"]+)"',
                        (APP / "primerforge" / "__init__.py").read_text(encoding="utf-8")).group(1)
    dist = PROJECT / "dist"
    dist.mkdir(exist_ok=True)
    target = dist / f"primerforge-installation-{version}.tar.gz"
    executable = {"install.sh", "uninstall.sh", "primerforge.sh"}
    with tarfile.open(target, "w:gz", format=tarfile.PAX_FORMAT) as tar:
        for path in sorted(HERE.rglob("*")):
            rel = path.relative_to(HERE)
            if SKIP_DIRS & set(rel.parts) or path.name.endswith(".part") or path.suffix in SKIP_SUFFIXES:
                continue
            info = tar.gettarinfo(str(path), arcname=f"installation/{rel.as_posix()}")
            info.uid = info.gid = 0
            info.uname = info.gname = "root"
            if path.is_dir():
                info.mode = 0o755
                tar.addfile(info)
            else:
                is_program = rel.as_posix() in executable or rel.parts[:2] == ("dependencies", "blat")
                info.mode = 0o755 if is_program else 0o644
                with open(path, "rb") as fh:
                    tar.addfile(info, fh)
    return target


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refresh", action="store_true", help="download everything again")
    ap.add_argument("--archive", action="store_true",
                    help="also write dist/primerforge-installation-<version>.tar.gz")
    ap.add_argument("--no-dependencies", action="store_true",
                    help="only refresh app/ (the server then downloads its dependencies)")
    args = ap.parse_args()

    copy_app()
    normalise_installer()
    versions = {}
    if not args.no_dependencies:
        try:
            python_packages(args.refresh)
            versions["NCBI BLAST+"] = blast(args.refresh)
            blat(args.refresh)
            versions["Caddy"] = caddy(args.refresh)
        except Exception as exc:                               # noqa: BLE001
            say(f"\nStopped: {exc}")
            return 1
        write_checksums(versions)
    size = sum(p.stat().st_size for p in HERE.rglob("*") if p.is_file())
    say(f"\nThe installation folder is ready: {HERE} ({human(size)})")
    if args.archive:
        target = archive()
        say(f"Archive: {target} ({human(target.stat().st_size)})")
    say("Copy the whole folder to the server, then run:  sudo bash install.sh")
    return 0


if __name__ == "__main__":
    sys.exit(main())
