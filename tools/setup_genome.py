"""One-shot installer for the local specificity stack.

Downloads NCBI BLAST+ and the soft-masked GRCh38 primary assembly, builds a
FASTA index and a BLAST nucleotide database. Resumable: every stage is skipped
if its output already exists. Progress is mirrored to data/setup_status.json so
the web UI can show a live status.

    python tools/setup_genome.py            # run everything
    python tools/setup_genome.py --status   # print current state and exit
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import subprocess
import sys
import tarfile
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from primerforge import config                                    # noqa: E402
from primerforge.refseq import FastaIndex                         # noqa: E402

STATUS_PATH = config.DATA_DIR / "setup_status.json"
GENOME_GZ = config.GENOME_DIR / "genome.fa.gz"
BLAST_TGZ = config.BLAST_DIR / "blast.tar.gz"


def set_status(stage: str, message: str, pct: float | None = None,
               done: bool = False, error: str | None = None) -> None:
    config.ensure_dirs()
    payload = {"stage": stage, "message": message, "pct": pct, "done": done,
               "error": error, "updated_at": time.time(), "pid": os.getpid()}
    STATUS_PATH.write_text(json.dumps(payload), encoding="utf-8")
    tail = f" {pct:.1f}%" if pct is not None else ""
    print(f"[{stage}]{tail} {message}", flush=True)


def read_status() -> dict:
    if STATUS_PATH.exists():
        try:
            return json.loads(STATUS_PATH.read_text(encoding="utf-8"))
        except ValueError:
            pass
    return {"stage": "idle", "message": "Not started", "done": False}


def human(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} GB"


def download(url: str, dest: Path, stage: str, label: str) -> None:
    """Resumable HTTP download with byte-range continuation."""
    if dest.exists() and dest.stat().st_size > 0:
        head = urllib.request.Request(url, method="HEAD")
        try:
            with urllib.request.urlopen(head, timeout=60) as r:
                total = int(r.headers.get("Content-Length", 0))
            if total and dest.stat().st_size >= total:
                set_status(stage, f"{label} already downloaded", 100.0)
                return
        except Exception:
            pass

    have = dest.stat().st_size if dest.exists() else 0
    req = urllib.request.Request(url, headers={"User-Agent": config.USER_AGENT})
    if have:
        req.add_header("Range", f"bytes={have}-")
    mode = "ab" if have else "wb"
    started = time.time()
    with urllib.request.urlopen(req, timeout=120) as resp:
        total = int(resp.headers.get("Content-Length", 0)) + have
        with open(dest, mode) as out:
            got = have
            last = 0.0
            while True:
                chunk = resp.read(1 << 20)
                if not chunk:
                    break
                out.write(chunk)
                got += len(chunk)
                now = time.time()
                if now - last > 2:
                    last = now
                    rate = (got - have) / max(now - started, 1e-6)
                    pct = 100.0 * got / total if total else None
                    eta = (total - got) / rate if total and rate else 0
                    set_status(stage,
                               f"{label}: {human(got)} of {human(total)} "
                               f"at {human(rate)}/s, ~{eta / 60:.0f} min left", pct)
    set_status(stage, f"{label} downloaded ({human(dest.stat().st_size)})", 100.0)


def blast_url() -> str:
    """The current BLAST+ build for this platform; NCBI's LATEST directory only keeps
    the newest release, so a pinned file name eventually disappears."""
    import re
    suffix = config.blast_package_suffix()
    try:
        req = urllib.request.Request(config.BLAST_LATEST, headers={"User-Agent": config.USER_AGENT})
        with urllib.request.urlopen(req, timeout=60) as r:
            listing = r.read().decode("utf-8", errors="ignore")
        names = sorted(set(re.findall(r'href="(ncbi-blast-[\d.]+\+-' + re.escape(suffix) + ')"',
                                      listing)))
        if names:
            return config.BLAST_LATEST + names[-1]
    except Exception:                                  # noqa: BLE001
        pass
    return config.BLAST_URL


def install_blast() -> None:
    if config.blast_bin("blastn"):
        set_status("blast", "BLAST+ already installed", 100.0)
        return
    config.ensure_dirs()
    download(blast_url(), BLAST_TGZ, "blast", "BLAST+ toolkit")
    set_status("blast", "Extracting BLAST+ ...")
    with tarfile.open(BLAST_TGZ, "r:gz") as tar:
        try:
            tar.extractall(config.BLAST_DIR, filter="data")
        except TypeError:                     # Python < 3.12
            tar.extractall(config.BLAST_DIR)
    BLAST_TGZ.unlink(missing_ok=True)
    if not config.blast_bin("blastn"):
        raise RuntimeError("BLAST+ extracted but blastn was not found.")
    set_status("blast", "BLAST+ installed", 100.0)


def install_genome() -> None:
    config.ensure_dirs()
    if not config.GENOME_FASTA.exists():
        download(config.GENOME_URL, GENOME_GZ, "genome", "GRCh38 genome")
        set_status("genome", "Decompressing genome (this takes a few minutes) ...")
        tmp = config.GENOME_FASTA.with_suffix(".part")
        total = GENOME_GZ.stat().st_size
        with gzip.open(GENOME_GZ, "rb") as src, open(tmp, "wb") as dst:
            written = 0
            last = 0.0
            while True:
                chunk = src.read(1 << 22)
                if not chunk:
                    break
                dst.write(chunk)
                written += len(chunk)
                now = time.time()
                if now - last > 3:
                    last = now
                    set_status("genome", f"Decompressing: {human(written)} written")
        tmp.replace(config.GENOME_FASTA)
        GENOME_GZ.unlink(missing_ok=True)
        set_status("genome", f"Genome ready ({human(config.GENOME_FASTA.stat().st_size)})",
                   100.0)
    else:
        set_status("genome", "Genome FASTA already present", 100.0)

    if not config.GENOME_FAI.exists():
        set_status("index", "Building FASTA index ...")
        FastaIndex.build(config.GENOME_FASTA, config.GENOME_FAI,
                         progress=lambda n: set_status("index", f"Indexing contig {n}"))
        set_status("index", "FASTA index built", 100.0)
    else:
        set_status("index", "FASTA index already present", 100.0)


def build_blastdb() -> None:
    if config.blast_ready():
        set_status("blastdb", "BLAST database already built", 100.0)
        return
    makeblastdb = config.blast_bin("makeblastdb")
    if not makeblastdb:
        raise RuntimeError("makeblastdb not found; install BLAST+ first.")
    config.BLAST_DB.parent.mkdir(parents=True, exist_ok=True)
    set_status("blastdb", "Building BLAST database (10-25 minutes) ...")
    proc = subprocess.Popen(
        [makeblastdb, "-in", str(config.GENOME_FASTA), "-dbtype", "nucl",
         "-out", str(config.BLAST_DB), "-title", "GRCh38", "-max_file_sz", "3GB"],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
        creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
    assert proc.stdout
    for line in proc.stdout:
        line = line.strip()
        if line:
            set_status("blastdb", line[:180])
    if proc.wait() != 0:
        raise RuntimeError("makeblastdb failed; see the log above.")
    set_status("blastdb", "BLAST database built", 100.0)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--status", action="store_true", help="print status and exit")
    ap.add_argument("--skip-blast", action="store_true")
    args = ap.parse_args()

    if args.status:
        print(json.dumps(read_status(), indent=2))
        print(json.dumps({"genome_ready": config.genome_ready(),
                          "blast_ready": config.blast_ready()}, indent=2))
        return 0

    try:
        if not args.skip_blast:
            install_blast()
        install_genome()
        build_blastdb()
    except Exception as exc:                      # surfaced in the UI
        set_status("error", f"Setup failed: {exc}", error=str(exc))
        print(f"FAILED: {exc}", file=sys.stderr)
        return 1

    free = shutil.disk_usage(config.DATA_DIR).free
    set_status("ready", f"Local specificity checking is enabled. "
                        f"{human(free)} disk free.", 100.0, done=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
