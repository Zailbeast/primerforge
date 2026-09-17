"""Central configuration and on-disk layout for PrimerForge."""
from __future__ import annotations

import os
import platform
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = Path(os.environ.get("PRIMERFORGE_DATA", ROOT / "data"))
GENOME_DIR = DATA_DIR / "genome"
BLAST_DIR = DATA_DIR / "blast"
DB_PATH = DATA_DIR / "primerforge.sqlite"

# Reference assembly. Ensembl REST serves GRCh38 at rest.ensembl.org and
# GRCh37 at grch37.rest.ensembl.org; the local genome is always GRCh38.
ASSEMBLY = "GRCh38"
ENSEMBL_REST = "https://rest.ensembl.org"
ENSEMBL_REST_GRCH37 = "https://grch37.rest.ensembl.org"

# Local reference files, created by tools/setup_genome.py.
GENOME_FASTA = GENOME_DIR / "Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
GENOME_FAI = GENOME_DIR / "Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.fai"
BLAST_DB = GENOME_DIR / "blastdb" / "GRCh38"

GENOME_URL = (
    "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/"
    "Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
)
BLAST_LATEST = "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
BLAST_VERSION = "2.17.0"             # used when the LATEST listing cannot be read


def blast_package_suffix() -> str:
    """The BLAST+ build for this machine, as named in NCBI's download directory."""
    if os.name == "nt":
        return "x64-win64.tar.gz"
    machine = platform.machine().lower()
    arch = "aarch64" if machine in ("aarch64", "arm64") else "x64"
    return f"{arch}-{'macosx' if platform.system() == 'Darwin' else 'linux'}.tar.gz"


BLAST_URL = f"{BLAST_LATEST}ncbi-blast-{BLAST_VERSION}+-{blast_package_suffix()}"

# Server mode (a shared Linux install) is configured through the environment, which
# the installer writes to /etc/primerforge/primerforge.env. The desktop app leaves
# these unset: no logins, listening on this computer only.
def _flag(name: str, default: bool = False) -> bool:
    value = os.environ.get(name)
    return default if value is None else value.strip().lower() in ("1", "true", "yes", "on")


AUTH_ENABLED = _flag("PRIMERFORGE_AUTH")
HOST = os.environ.get("PRIMERFORGE_HOST", "127.0.0.1")
PORT = int(os.environ.get("PRIMERFORGE_PORT", "5000"))
TRUST_PROXY = _flag("PRIMERFORGE_TRUST_PROXY")          # behind Caddy/nginx on this host
SECURE_COOKIES = _flag("PRIMERFORGE_SECURE_COOKIES")    # site is served over HTTPS
THREADS = int(os.environ.get("PRIMERFORGE_THREADS", "24"))

HTTP_TIMEOUT = 60
USER_AGENT = "PrimerForge/1.0 (local primer design tool)"

# BLAST/BLAT search tool. Each installed species gets its own directory holding
# its FASTA files, BLAST databases, annotation and a manifest.json.
SPECIES_DIR = DATA_DIR / "species"
BLAST_JOBS_DIR = DATA_DIR / "blast_jobs"
BLAT_BIN_CACHE = DATA_DIR / "blat_bin"
SETTINGS_PATH = DATA_DIR / "settings.json"
ENSEMBL_FTP = "https://ftp.ensembl.org/pub"
ENSEMBL_WEB = "https://www.ensembl.org"
ENSEMBL_WEB_GRCH37 = "https://grch37.ensembl.org"
NCBI_WINDOWMASKER = "https://ftp.ncbi.nlm.nih.gov/blast/windowmasker_files"
# UCSC's main download host is intermittently unreachable from some networks;
# the European mirror carries identical binaries.
UCSC_EXE_MIRRORS = (
    "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64",
    "https://hgdownload-euro.soe.ucsc.edu/admin/exe/linux.x86_64",
)
BLAT_BASE_PORT = 17779

# On Windows BLAT runs inside WSL. Its binaries, 2bit genomes and gfServer indexes live
# on the Linux filesystem (relative paths are resolved against the WSL user's $HOME)
# because gfServer memory-maps its index and that is slow across the /mnt/<drive>
# bridge. On Linux BLAT runs directly and relative paths are resolved against DATA_DIR.
DEFAULT_SETTINGS = {
    "wsl_distro": os.environ.get("PRIMERFORGE_WSL_DISTRO", ""),   # "" = WSL default
    "blat_dir": os.environ.get("PRIMERFORGE_BLAT_DIR",
                               "primerforge/blat" if os.name == "nt" else "blat"),
    "blat_step_size": 5,             # what UCSC uses for its public BLAT servers
    "blast_workers": 2,
    "blast_threads": 0,              # 0 = share the CPUs between workers
    "blat_autostart": True,
}


def settings() -> dict:
    import json

    out = dict(DEFAULT_SETTINGS)
    if SETTINGS_PATH.exists():
        try:
            out.update(json.loads(SETTINGS_PATH.read_text(encoding="utf-8")))
        except ValueError:
            pass
    return out


def save_settings(values: dict) -> dict:
    import json

    current = settings()
    current.update({k: v for k, v in values.items() if k in DEFAULT_SETTINGS})
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    SETTINGS_PATH.write_text(json.dumps(current, indent=2), encoding="utf-8")
    return current


def blast_bin(name: str) -> str | None:
    """Locate a BLAST+ executable, preferring the bundled copy over PATH."""
    exe = name + (".exe" if os.name == "nt" else "")
    for candidate in BLAST_DIR.rglob(exe):
        if candidate.is_file():
            return str(candidate)
    from shutil import which

    return which(name)


def genome_ready() -> bool:
    return GENOME_FASTA.exists() and GENOME_FAI.exists()


def blast_ready() -> bool:
    return blast_bin("blastn") is not None and bool(
        list(BLAST_DB.parent.glob("GRCh38*.nsq")) or list(BLAST_DB.parent.glob("GRCh38*.nal"))
    )


def ensure_dirs() -> None:
    for d in (DATA_DIR, GENOME_DIR, BLAST_DIR, BLAST_DB.parent, SPECIES_DIR, BLAST_JOBS_DIR):
        d.mkdir(parents=True, exist_ok=True)
