"""Installed species and the BLAST/BLAT databases available for each.

Every species lives in data/species/<production_name>/ with a manifest.json that
records where each component was installed. The human genome and FASTA index that
PrimerForge's primer design already uses are adopted in place rather than
downloaded a second time.
"""
from __future__ import annotations

import json
import re
import threading
import time
from pathlib import Path

from . import config
from .blastconf import GENOMIC_SOURCES, SOURCES

COMPONENTS = {
    "karyotype": "Chromosomes and cytogenetic bands",
    "genome":    "Genome sequence (soft-masked FASTA)",
    "blastdb":   "BLAST genome database with repeat masks",
    "cdna":      "cDNA BLAST database",
    "ncrna":     "Non-coding RNA BLAST database",
    "pep":       "Protein BLAST database",
    "gtf":       "Gene annotation (genes, transcripts, exons)",
    "xrefs":     "Transcript identifiers (MANE Select, RefSeq NM/NP)",
    "refseq":    "RefSeq transcript BLAST database (NM/NR)",
    "repeats":   "WindowMasker repeat filter for queries",
    "blat":      "BLAT index (2bit genome + gfServer index)",
}
# BLAT only needs the genome, so it is built before the slower BLAST databases.
INSTALL_ORDER = ["karyotype", "genome", "blat", "blastdb", "gtf", "xrefs", "cdna", "ncrna",
                 "pep", "refseq", "repeats"]

_lock = threading.Lock()


def species_dir(name: str) -> Path:
    return config.SPECIES_DIR / name


def manifest_path(name: str) -> Path:
    return species_dir(name) / "manifest.json"


def _read(path: Path) -> dict | None:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def load(name: str) -> dict | None:
    return _read(manifest_path(name))


def save(manifest: dict) -> None:
    path = manifest_path(manifest["name"])
    path.parent.mkdir(parents=True, exist_ok=True)
    with _lock:
        tmp = path.with_suffix(".tmp")
        tmp.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        tmp.replace(path)


def update(name: str, **fields) -> dict:
    with _lock:
        manifest = _read(manifest_path(name)) or {"name": name}
    manifest.update(fields)
    manifest["updated_at"] = time.time()
    save(manifest)
    return manifest


def web_name(name: str) -> str:
    """Ensembl URL form: homo_sapiens -> Homo_sapiens."""
    return name[:1].upper() + name[1:]


def ensure_human() -> None:
    """Register the genome PrimerForge already installed as the human species."""
    existing = load("homo_sapiens")
    if existing is not None:
        # Registered before the primer-design genome finished downloading: adopt it now
        # rather than letting the species installer download it a second time.
        if not _exists((existing.get("genome") or {}).get("fai")) and config.genome_ready():
            update("homo_sapiens", genome={"fasta": str(config.GENOME_FASTA),
                                           "fai": str(config.GENOME_FAI),
                                           "file": config.GENOME_FASTA.name})
        return
    manifest = {
        "name": "homo_sapiens", "display_name": "Human", "scientific_name": "Homo sapiens",
        "common_name": "human", "assembly": config.ASSEMBLY, "taxon_id": 9606,
        "division": "EnsemblVertebrates", "created_at": time.time(),
    }
    if config.genome_ready():
        manifest["genome"] = {"fasta": str(config.GENOME_FASTA),
                              "fai": str(config.GENOME_FAI),
                              "file": config.GENOME_FASTA.name}
    save(manifest)


def installed_names() -> list[str]:
    if not config.SPECIES_DIR.exists():
        return []
    return sorted(p.parent.name for p in config.SPECIES_DIR.glob("*/manifest.json"))


def all_species() -> list[dict]:
    out = []
    for name in installed_names():
        m = load(name)
        if m:
            out.append(summary(m))
    out.sort(key=lambda s: (s["name"] != "homo_sapiens", s["display_name"]))
    return out


def _exists(path: str | None) -> bool:
    return bool(path) and Path(path).exists()


def _blastdb_exists(prefix: str | None, kind: str = "n") -> bool:
    if not prefix:
        return False
    p = Path(prefix)
    return any(p.parent.glob(f"{p.name}*.{kind}sq")) or any(p.parent.glob(f"{p.name}*.{kind}al"))


def available_sources(m: dict) -> dict[str, bool]:
    db = m.get("blastdb") or {}
    genome_db = _blastdb_exists(db.get("prefix"))
    masked = genome_db and db.get("mask_id") is not None
    return {
        "LATESTGP": genome_db,
        "LATESTGP_MASKED": masked,
        "LATESTGP_SOFT": masked,
        "CDNA_ALL": _blastdb_exists((m.get("cdna") or {}).get("prefix")),
        "NCRNA": _blastdb_exists((m.get("ncrna") or {}).get("prefix")),
        "PEP_ALL": _blastdb_exists((m.get("pep") or {}).get("prefix"), "p"),
        "REFSEQ_RNA": _blastdb_exists((m.get("refseq") or {}).get("prefix")),
        "MANE_SELECT": _exists((m.get("refseq") or {}).get("mane_alias")),
    }


def blat_installed(m: dict) -> bool:
    b = m.get("blat") or {}
    return bool(b.get("ready") and b.get("dir") and b.get("port"))


def component_state(m: dict) -> dict[str, bool]:
    sources = available_sources(m)
    return {
        "karyotype": bool(m.get("chromosomes")),
        "genome": _exists((m.get("genome") or {}).get("fai")),
        "blastdb": sources["LATESTGP_MASKED"],
        "cdna": sources["CDNA_ALL"],
        "ncrna": sources["NCRNA"],
        "pep": sources["PEP_ALL"],
        "gtf": _exists((m.get("gtf") or {}).get("sqlite")),
        "xrefs": _exists((m.get("transcripts") or {}).get("sqlite")),
        "refseq": sources["REFSEQ_RNA"],
        "repeats": _exists((m.get("repeats") or {}).get("windowmasker")),
        "blat": blat_installed(m),
    }


def summary(m: dict) -> dict:
    return {
        "name": m["name"],
        "display_name": m.get("display_name") or m["name"].replace("_", " ").capitalize(),
        "scientific_name": m.get("scientific_name") or m["name"].replace("_", " ").capitalize(),
        "assembly": m.get("assembly", ""),
        "release": m.get("release"),
        "taxon_id": m.get("taxon_id"),
        "sources": available_sources(m),
        "components": component_state(m),
        "blat": blat_installed(m),
        "has_karyotype": bool(m.get("karyotype")),
    }


def blast_db(m: dict, source: str) -> tuple[str, str | None]:
    """(database prefix, mask algorithm id) for a data source."""
    if source in GENOMIC_SOURCES:
        db = m.get("blastdb") or {}
        mask = str(db["mask_id"]) if source != "LATESTGP" and db.get("mask_id") is not None \
            else None
        return db.get("prefix", ""), mask
    if source == "MANE_SELECT":
        return (m.get("refseq") or {}).get("mane_prefix", ""), None
    key = {"CDNA_ALL": "cdna", "NCRNA": "ncrna", "PEP_ALL": "pep", "REFSEQ_RNA": "refseq"}[source]
    return (m.get(key) or {}).get("prefix", ""), None


def source_label(source: str) -> str:
    return SOURCES.get(source, {}).get("label", source)


# ---------------------------------------------------------------------------
# Karyotype
# ---------------------------------------------------------------------------
_CHROM_RE = re.compile(r"^(chr)?([0-9]+|[XYZWM]|MT|[IVX]+)$", re.I)


def karyotype(m: dict) -> list[dict]:
    """Chromosomes to draw, in karyotype order, with lengths and bands."""
    chroms = m.get("chromosomes") or {}
    order = m.get("karyotype") or []
    return [{"name": c, "length": chroms[c]["length"], "bands": chroms[c].get("bands", [])}
            for c in order if c in chroms]


def karyotype_from_fai(fai: Path) -> tuple[list[str], dict]:
    """Fallback karyotype when Ensembl's assembly endpoint is unavailable."""
    names, chroms = [], {}
    for line in fai.read_text(encoding="utf-8").splitlines():
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        name, length = parts[0], int(parts[1])
        if _CHROM_RE.match(name):
            names.append(name)
            chroms[name] = {"length": length, "bands": []}

    def key(n: str):
        core = n.lower().removeprefix("chr")
        if core.isdigit():
            return (0, int(core), "")
        return (1, 0, {"x": "0", "y": "1", "w": "2", "z": "3"}.get(core, "9" + core))

    names.sort(key=key)
    return names, chroms


def contig_lengths(m: dict) -> dict[str, int]:
    fai = (m.get("genome") or {}).get("fai")
    out: dict[str, int] = {}
    if fai and Path(fai).exists():
        for line in Path(fai).read_text(encoding="utf-8").splitlines():
            parts = line.split("\t")
            if len(parts) >= 2:
                out[parts[0]] = int(parts[1])
    return out


def ensembl_site(m: dict) -> str:
    return config.ENSEMBL_WEB_GRCH37 if m.get("assembly") == "GRCh37" else config.ENSEMBL_WEB
