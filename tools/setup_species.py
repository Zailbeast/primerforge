"""Install BLAST/BLAT databases for an Ensembl species.

    python tools/setup_species.py homo_sapiens                   # everything
    python tools/setup_species.py mus_musculus --components genome,blastdb,blat
    python tools/setup_species.py homo_sapiens --status

Components (each skipped when already installed, so re-running is resumable):

  karyotype  chromosome lengths and cytogenetic bands (Ensembl REST)
  genome     soft-masked genome FASTA + index (Ensembl FTP)
  blastdb    BLAST nucleotide database carrying the soft-mask as a mask track, so
             one database serves "genomic", "hard masked" and "soft masked" searches
  cdna/ncrna/pep   transcript, non-coding RNA and protein BLAST databases
  gtf        gene annotation index (overlapping genes, exons, cDNA->genome mapping)
  repeats    NCBI WindowMasker statistics used to repeat-mask queries
  blat       UCSC BLAT (in WSL on Windows, directly on Linux): binaries, a 2bit
             genome and a gfServer index

Progress is mirrored to data/species/<name>/install_status.json for the web UI.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from primerforge import annotation, config, species, transcripts, wsl   # noqa: E402
from primerforge.refseq import FastaIndex                        # noqa: E402

NO_WINDOW = getattr(subprocess, "CREATE_NO_WINDOW", 0)
BLAT_BINARIES = ("blat/gfServer", "blat/gfClient", "blat/blat", "faToTwoBit", "twoBitInfo")


class SetupError(RuntimeError):
    pass


# ---------------------------------------------------------------------------
# Status
# ---------------------------------------------------------------------------
class Status:
    def __init__(self, name: str, components: list[str]):
        self.path = species.species_dir(name) / "install_status.json"
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.data = {"species": name, "started_at": time.time(), "done": False,
                     "pid": os.getpid(),
                     "error": None, "component": None, "message": "Starting",
                     "pct": None,
                     "components": {c: {"state": "pending", "message": ""} for c in components}}
        self.write()

    def write(self) -> None:
        self.data["updated_at"] = time.time()
        tmp = self.path.with_suffix(".tmp")
        tmp.write_text(json.dumps(self.data), encoding="utf-8")
        tmp.replace(self.path)

    def set(self, component: str, message: str, pct: float | None = None,
            state: str = "running") -> None:
        self.data.update(component=component, message=message, pct=pct)
        self.data["components"].setdefault(component, {})
        self.data["components"][component].update(state=state, message=message)
        self.write()
        tail = f" {pct:.0f}%" if pct is not None else ""
        print(f"[{component}]{tail} {message}", flush=True)


def read_status(name: str) -> dict:
    path = species.species_dir(name) / "install_status.json"
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def human(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} GB"


# ---------------------------------------------------------------------------
# Network helpers
# ---------------------------------------------------------------------------
def list_dir(url: str) -> list[str]:
    req = urllib.request.Request(url, headers={"User-Agent": config.USER_AGENT})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                html = r.read().decode("utf-8", errors="ignore")
            return re.findall(r'href="([^"?/][^"]*)"', html)
        except urllib.error.HTTPError as exc:
            if exc.code == 404 or attempt == 3:
                raise SetupError(f"Could not list {url}: {exc}") from exc
            time.sleep(3 * (attempt + 1))
        except Exception as exc:                               # noqa: BLE001
            if attempt == 3:
                raise SetupError(f"Could not list {url}: {exc}") from exc
            time.sleep(3 * (attempt + 1))
    return []


def download(url: str, dest: Path, status: Status, component: str, label: str) -> None:
    """Resumable download with byte-range continuation and retries."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    total = 0
    for attempt in range(6):
        have = part.stat().st_size if part.exists() else 0
        req = urllib.request.Request(url, headers={"User-Agent": config.USER_AGENT})
        if have:
            req.add_header("Range", f"bytes={have}-")
        try:
            with urllib.request.urlopen(req, timeout=120) as resp:
                if have and resp.status != 206:          # server ignored the range
                    have = 0
                total = int(resp.headers.get("Content-Length", 0)) + have
                started, last, got = time.time(), 0.0, have
                with open(part, "ab" if have else "wb") as out:
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
                            status.set(component, f"{label}: {human(got)} of {human(total)} "
                                                  f"at {human(rate)}/s", pct)
            if total and part.stat().st_size < total:
                raise SetupError("download ended early")
            part.replace(dest)
            status.set(component, f"{label} downloaded ({human(dest.stat().st_size)})", 100)
            return
        except urllib.error.HTTPError as exc:
            if exc.code == 416 and part.exists():            # already complete
                part.replace(dest)
                return
            if exc.code in (403, 404):
                raise SetupError(f"{label} not found at {url} (HTTP {exc.code})") from exc
            err = exc
        except Exception as exc:                               # noqa: BLE001
            err = exc
        status.set(component, f"{label}: retrying after error ({err})")
        time.sleep(min(30, 4 * (attempt + 1)))
    raise SetupError(f"Failed to download {label} from {url}")


def gunzip(src: Path, dest: Path, status: Status, component: str, label: str) -> None:
    tmp = dest.with_name(dest.name + ".part")
    written, last = 0, 0.0
    with gzip.open(src, "rb") as fin, open(tmp, "wb") as fout:
        while True:
            chunk = fin.read(1 << 22)
            if not chunk:
                break
            fout.write(chunk)
            written += len(chunk)
            if time.time() - last > 3:
                last = time.time()
                status.set(component, f"Decompressing {label}: {human(written)}")
    tmp.replace(dest)


def rest_json(path: str):
    req = urllib.request.Request(config.ENSEMBL_REST + path,
                                 headers={"User-Agent": config.USER_AGENT,
                                          "Accept": "application/json"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.loads(r.read().decode())
        except Exception:                                      # noqa: BLE001
            if attempt == 3:
                return None
            time.sleep(4 * (attempt + 1))
    return None


def run(cmd: list, status: Status, component: str, label: str, cwd=None) -> str:
    status.set(component, label)
    proc = subprocess.Popen([str(c) for c in cmd], stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, cwd=cwd,
                            encoding="utf-8", errors="replace", creationflags=NO_WINDOW)
    lines = []
    assert proc.stdout
    for line in proc.stdout:
        line = line.strip()
        if line:
            lines.append(line)
            status.set(component, f"{label} {line[:160]}")
    if proc.wait() != 0:
        raise SetupError(f"{Path(str(cmd[0])).name} failed: {' | '.join(lines[-4:])}")
    return "\n".join(lines)


def ensembl_species_info(name: str) -> dict:
    cache = config.SPECIES_DIR / "ensembl_species.json"
    data = None
    if cache.exists() and time.time() - cache.stat().st_mtime < 7 * 86400:
        try:
            data = json.loads(cache.read_text(encoding="utf-8"))
        except ValueError:
            data = None
    if data is None:
        data = rest_json("/info/species?content-type=application/json")
        if data:
            cache.parent.mkdir(parents=True, exist_ok=True)
            cache.write_text(json.dumps(data), encoding="utf-8")
    for s in (data or {}).get("species", []):
        if s.get("name") == name:
            return s
    return {}


# ---------------------------------------------------------------------------
# Components
# ---------------------------------------------------------------------------
def c_karyotype(m: dict, st: Status) -> None:
    name = m["name"]
    st.set("karyotype", "Fetching assembly information from Ensembl")
    data = rest_json(f"/info/assembly/{name}?bands=1;content-type=application/json")
    if data and data.get("top_level_region"):
        regions = {r["name"]: r for r in data["top_level_region"]}
        karyo = [c for c in data.get("karyotype") or [] if c in regions]
        chroms = {}
        for c in karyo:
            bands = sorted(({"start": b["start"], "end": b["end"], "stain": b.get("stain"),
                             "id": b.get("id")} for b in regions[c].get("bands") or []),
                           key=lambda b: b["start"])
            chroms[c] = {"length": regions[c]["length"], "bands": bands}
        species.update(name, karyotype=karyo, chromosomes=chroms,
                       assembly=data.get("default_coord_system_version") or m.get("assembly"),
                       assembly_full=data.get("assembly_name"))
        st.set("karyotype", f"{len(karyo)} chromosomes", 100, state="done")
        return
    fai = (m.get("genome") or {}).get("fai")
    if not fai or not Path(fai).exists():
        raise SetupError("Ensembl assembly information unavailable and no genome index yet.")
    karyo, chroms = species.karyotype_from_fai(Path(fai))
    species.update(name, karyotype=karyo, chromosomes=chroms)
    st.set("karyotype", f"{len(karyo)} chromosomes (from the genome index)", 100, state="done")


def _ftp_species_dir(m: dict, kind: str) -> str:
    return f"{config.ENSEMBL_FTP}/current_fasta/{m['name']}/{kind}/"


def c_genome(m: dict, st: Status) -> None:
    g = m.get("genome") or {}
    if g.get("fai") and Path(g["fai"]).exists():
        st.set("genome", "Genome already installed", 100, state="done")
        return
    listing = list_dir(_ftp_species_dir(m, "dna"))
    pick = next((f for f in listing if f.endswith(".dna_sm.primary_assembly.fa.gz")), None) or \
        next((f for f in listing if f.endswith(".dna_sm.toplevel.fa.gz")), None)
    if not pick:
        raise SetupError("No soft-masked genome FASTA found on the Ensembl FTP site.")
    gdir = species.species_dir(m["name"]) / "genome"
    gz = gdir / pick
    fa = gdir / pick[:-3]
    if not fa.exists():
        download(_ftp_species_dir(m, "dna") + pick, gz, st, "genome", "Genome")
        gunzip(gz, fa, st, "genome", "genome")
        gz.unlink(missing_ok=True)
    fai = fa.with_name(fa.name + ".fai")
    st.set("genome", "Indexing genome FASTA")
    FastaIndex.build(fa, fai, progress=lambda n: st.set("genome", f"Indexing contig {n}"))
    species.update(m["name"], genome={"fasta": str(fa), "fai": str(fai), "file": fa.name})
    st.set("genome", f"Genome ready ({human(fa.stat().st_size)})", 100, state="done")


def c_blastdb(m: dict, st: Status) -> None:
    if species.available_sources(m)["LATESTGP_MASKED"]:
        st.set("blastdb", "Genome BLAST database already built", 100, state="done")
        return
    fasta = (m.get("genome") or {}).get("fasta")
    if not fasta or not Path(fasta).exists():
        raise SetupError("Install the genome component first.")
    bindir = _blast_bindir()
    out = species.species_dir(m["name"]) / "blastdb"
    out.mkdir(parents=True, exist_ok=True)
    masks = out / "repeats.asnt"
    prefix = out / "dna"
    run([bindir / "convert2blastmask", "-in", fasta, "-masking_algorithm", "repeat",
         "-masking_options", "repeatmasker,default", "-parse_seqids",
         "-outfmt", "maskinfo_asn1_text", "-out", masks],
        st, "blastdb", "Extracting soft-masked repeats (a few minutes)...")
    title = f"{m.get('scientific_name', m['name'])} {m.get('assembly', '')} genome".strip()
    run([bindir / "makeblastdb", "-in", fasta, "-dbtype", "nucl", "-parse_seqids",
         "-mask_data", masks, "-out", prefix, "-title", title, "-max_file_sz", "4GB"],
        st, "blastdb", "Building BLAST database...")
    info = run([bindir / "blastdbcmd", "-db", prefix, "-info"], st, "blastdb",
               "Reading database masks")
    match = re.search(r"^\s*(\d+)\s+repeat\b", info, re.M)
    if not match:
        raise SetupError("The database was built without its repeat mask track.")
    masks.unlink(missing_ok=True)
    species.update(m["name"], blastdb={"prefix": str(prefix), "mask_id": int(match.group(1)),
                                       "built_at": time.time()})
    st.set("blastdb", "Genome BLAST database ready (unmasked, hard- and soft-masked)",
           100, state="done")


def _blast_bindir() -> Path:
    exe = config.blast_bin("makeblastdb")
    if not exe:
        # A species can be installed before (or without) the human primer-design genome,
        # which is what normally brings BLAST+ along.
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from setup_genome import install_blast
        try:
            install_blast()
        except Exception as exc:                               # noqa: BLE001
            raise SetupError(f"NCBI BLAST+ is not installed and could not be downloaded: "
                             f"{exc}") from exc
        exe = config.blast_bin("makeblastdb")
    if not exe:
        raise SetupError("NCBI BLAST+ is not installed. Run Settings > Install local genome, "
                         "or python tools/setup_genome.py first.")
    return Path(exe).parent


def c_transcripts(m: dict, st: Status, kind: str) -> None:
    key = {"cdna": "CDNA_ALL", "ncrna": "NCRNA", "pep": "PEP_ALL"}[kind]
    if species.available_sources(m)[key]:
        st.set(kind, "Already installed", 100, state="done")
        return
    suffix = {"cdna": ".cdna.all.fa.gz", "ncrna": ".ncrna.fa.gz", "pep": ".pep.all.fa.gz"}[kind]
    listing = list_dir(_ftp_species_dir(m, kind))
    pick = next((f for f in listing if f.endswith(suffix)), None)
    if not pick:
        st.set(kind, f"Ensembl has no {kind} file for this species", state="skipped")
        return
    d = species.species_dir(m["name"]) / kind
    gz, fa = d / pick, d / pick[:-3]
    download(_ftp_species_dir(m, kind) + pick, gz, st, kind, kind)
    gunzip(gz, fa, st, kind, kind)
    gz.unlink(missing_ok=True)
    prefix = d / kind
    title = f"{m.get('scientific_name', m['name'])} {m.get('assembly', '')} {kind}".strip()
    run([_blast_bindir() / "makeblastdb", "-in", fa, "-dbtype",
         "prot" if kind == "pep" else "nucl", "-parse_seqids", "-out", prefix,
         "-title", title], st, kind, "Building BLAST database...")
    fa.unlink(missing_ok=True)
    species.update(m["name"], **{kind: {"prefix": str(prefix), "file": pick,
                                        "built_at": time.time()}})
    st.set(kind, "BLAST database ready", 100, state="done")


def c_gtf(m: dict, st: Status) -> None:
    g = m.get("gtf") or {}
    if g.get("sqlite") and Path(g["sqlite"]).exists():
        st.set("gtf", "Annotation already installed", 100, state="done")
        return
    # Ensembl moved the GTF tree from current_gtf/ to current/gtf/; accept either.
    listing, base = [], ""
    for base in (f"{config.ENSEMBL_FTP}/current/gtf/{m['name']}/",
                 f"{config.ENSEMBL_FTP}/current_gtf/{m['name']}/"):
        try:
            listing = list_dir(base)
            break
        except SetupError:
            continue
    pick = next((f for f in listing if re.search(r"\.\d+\.gtf\.gz$", f)
                 and ".chr." not in f and "abinitio" not in f and "patch" not in f), None)
    if not pick:
        raise SetupError("No GTF annotation found on the Ensembl FTP site.")
    d = species.species_dir(m["name"])
    gz = d / pick
    download(base + pick, gz, st, "gtf", "Annotation")
    release = re.search(r"\.(\d+)\.gtf\.gz$", pick)
    sqlite_path = d / "annotation.sqlite"
    st.set("gtf", "Indexing genes, transcripts and exons (1-3 minutes)...")
    annotation.reset()
    counts = annotation.build(gz, sqlite_path, progress=lambda c: st.set(
        "gtf", f"Indexed {c['genes']:,} genes, {c['transcripts']:,} transcripts, "
               f"{c['exons']:,} exons"))
    gz.unlink(missing_ok=True)
    species.update(m["name"], gtf={"sqlite": str(sqlite_path), "file": pick, **counts},
                   release=int(release.group(1)) if release else m.get("release"))
    st.set("gtf", f"{counts['genes']:,} genes, {counts['transcripts']:,} transcripts",
           100, state="done")


MANE_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/"
REFSEQ_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/"


def _refseq_dir(m: dict) -> str | None:
    """RefSeq's per-organism directory, e.g. Homo sapiens -> H_sapiens."""
    parts = (m.get("scientific_name") or m["name"].replace("_", " ")).split()
    if len(parts) < 2:
        return None
    return f"{parts[0][0].upper()}_{parts[1].lower()}"


def c_xrefs(m: dict, st: Status) -> None:
    """MANE Select records and Ensembl's RefSeq cross-references."""
    d = species.species_dir(m["name"])
    sqlite_path = transcripts.db_path(d)
    counts = {}

    tsv = None
    for base in (f"{config.ENSEMBL_FTP}/current/tsv/{m['name']}/",
                 f"{config.ENSEMBL_FTP}/current_tsv/{m['name']}/"):
        try:
            listing = list_dir(base)
        except SetupError:
            continue
        pick = next((f for f in listing if f.endswith(".refseq.tsv.gz")), None)
        if pick:
            tsv = d / pick
            if not tsv.exists():
                download(base + pick, tsv, st, "xrefs", "RefSeq cross-references")
            break
    if tsv is None:
        st.set("xrefs", "Ensembl publishes no RefSeq cross-reference file for this species")

    mane = mane_name = None
    if str(m.get("taxon_id")) == "9606":
        mane_name = next((f for f in list_dir(MANE_URL) if f.endswith(".summary.txt.gz")), None)
        if mane_name:
            mane = d / mane_name
            if not mane.exists():
                download(MANE_URL + mane_name, mane, st, "xrefs", "MANE summary")
    else:
        st.set("xrefs", "MANE covers human only; installing cross-references alone")

    if tsv is None and mane is None:
        raise SetupError("No transcript identifier tables are available for this species.")
    st.set("xrefs", "Indexing transcript identifiers...")
    counts = transcripts.build_xrefs(sqlite_path, tsv, mane,
                                     progress=lambda c: st.set("xrefs", f"{c['xrefs']:,} cross-references"))
    for f in (tsv, mane):
        if f:
            f.unlink(missing_ok=True)
    species.update(m["name"], transcripts={"sqlite": str(sqlite_path), **counts,
                                           "mane_release": mane_name})
    st.set("xrefs", f"{counts['xrefs']:,} RefSeq cross-references, {counts['mane']:,} MANE "
                    "transcripts", 100, state="done")


def c_refseq(m: dict, st: Status) -> None:
    """All RefSeq transcripts for the species as a BLAST database (plus a MANE-only alias)."""
    sources = species.available_sources(m)
    if sources["REFSEQ_RNA"] and sources["MANE_SELECT"]:
        st.set("refseq", "RefSeq transcripts already installed", 100, state="done")
        return
    if sources["REFSEQ_RNA"]:
        # Database is there but its MANE Select subset is not (or MANE arrived later).
        _mane_alias(m, st)
        return
    organism = _refseq_dir(m)
    base = f"{REFSEQ_URL}{organism}/mRNA_Prot/"
    try:
        listing = list_dir(base)
    except SetupError:
        st.set("refseq", f"RefSeq has no per-organism transcript set for {organism}",
               state="skipped")
        return
    files = sorted({f for f in listing if re.match(r"^[a-z_]+\.\d+\.rna\.fna\.gz$", f)},
                   key=lambda f: int(f.split(".")[1]))
    if not files:
        st.set("refseq", "No RefSeq transcript files found", state="skipped")
        return
    d = species.species_dir(m["name"]) / "refseq"
    d.mkdir(parents=True, exist_ok=True)
    combined = d / "refseq_rna.fa"
    if not combined.exists():
        parts = []
        for i, name in enumerate(files, 1):
            gz = d / name
            if not gz.exists():
                download(base + name, gz, st, "refseq",
                         f"RefSeq transcripts {i}/{len(files)}")
            parts.append(gz)
        st.set("refseq", "Combining RefSeq transcript files...")
        tmp = combined.with_suffix(".part")
        with open(tmp, "wb") as out:
            for gz in parts:
                with gzip.open(gz, "rb") as src:
                    shutil.copyfileobj(src, out, 1 << 22)
        tmp.replace(combined)
        for gz in parts:
            gz.unlink(missing_ok=True)

    prefix = d / "refseq"
    run([_blast_bindir() / "makeblastdb", "-in", combined, "-dbtype", "nucl", "-parse_seqids",
         "-out", prefix, "-title", f"RefSeq transcripts {m.get('scientific_name', '')}".strip()],
        st, "refseq", "Building RefSeq BLAST database...")

    st.set("refseq", "Indexing RefSeq accessions, gene symbols and descriptions...")
    sqlite_path = transcripts.db_path(species.species_dir(m["name"]))
    n = transcripts.build_refseq_info(sqlite_path, [combined],
                                      progress=lambda t: st.set("refseq", f"{t:,} accessions"))
    combined.unlink(missing_ok=True)

    species.update(m["name"], refseq={"prefix": str(prefix), "n_transcripts": n,
                                      "files": len(files), "built_at": time.time()})
    st.set("refseq", f"{n:,} RefSeq transcripts", 100)
    _mane_alias(species.load(m["name"]) or m, st)


def _mane_alias(m: dict, st: Status) -> None:
    """A MANE Select subset of the RefSeq database, searchable on its own."""
    refseq_db = (m.get("refseq") or {}).get("prefix")
    sqlite_path = transcripts.db_path(species.species_dir(m["name"]))
    accessions = transcripts.mane_accessions(sqlite_path) if sqlite_path.exists() else []
    n = (m.get("refseq") or {}).get("n_transcripts", 0)
    if not accessions:
        st.set("refseq", f"{n:,} RefSeq transcripts (no MANE set for this species)",
               100, state="done")
        return
    d = Path(refseq_db).parent
    st.set("refseq", f"Building the MANE Select subset ({len(accessions):,} transcripts)...")
    ids = d / "mane_ids.txt"
    ids.write_text("\n".join(accessions) + "\n", encoding="utf-8")
    binary = d / "mane_ids.bsl"
    run([_blast_bindir() / "blastdb_aliastool", "-seqid_file_in", ids,
         "-seqid_file_out", binary], st, "refseq", "Preparing the MANE accession list")
    mane_prefix = d / "mane"
    run([_blast_bindir() / "blastdb_aliastool", "-db", refseq_db, "-dbtype", "nucl",
         "-seqidlist", binary, "-out", mane_prefix, "-title", "MANE Select transcripts"],
        st, "refseq", "Building the MANE Select database")
    ids.unlink(missing_ok=True)
    alias = str(mane_prefix) + ".nal"
    species.update(m["name"], refseq={**(m.get("refseq") or {}),
                                      "mane_prefix": str(mane_prefix) if Path(alias).exists() else None,
                                      "mane_alias": alias if Path(alias).exists() else None})
    st.set("refseq", f"{n:,} RefSeq transcripts, {len(accessions):,} MANE Select",
           100, state="done")


def c_repeats(m: dict, st: Status) -> None:
    r = m.get("repeats") or {}
    if r.get("windowmasker") and Path(r["windowmasker"]).exists():
        st.set("repeats", "Repeat filter already installed", 100, state="done")
        return
    d = species.species_dir(m["name"]) / "repeats"
    dest = d / "wmasker.obinary"
    taxid = m.get("taxon_id")
    source = "NCBI"
    try:
        if not taxid:
            raise SetupError("no taxonomy id")
        download(f"{config.NCBI_WINDOWMASKER}/{taxid}/wmasker.obinary", dest, st, "repeats",
                 "WindowMasker statistics")
    except SetupError as exc:
        fasta = (m.get("genome") or {}).get("fasta")
        if not fasta:
            raise SetupError(f"NCBI has no WindowMasker file ({exc}) and no genome is "
                             "installed to compute one.") from exc
        source = "computed"
        wm = _blast_bindir() / "windowmasker"
        d.mkdir(parents=True, exist_ok=True)
        run([wm, "-mk_counts", "-in", fasta, "-infmt", "fasta", "-sformat", "obinary",
             "-out", dest, "-mem", "4096"], st, "repeats",
            "Computing WindowMasker statistics from the genome (can take 30+ minutes)...")
    species.update(m["name"], repeats={"windowmasker": str(dest), "source": source})
    st.set("repeats", f"Repeat filter ready ({source})", 100, state="done")


def _blat_binaries_to_cache(st: Status) -> dict[str, Path]:
    config.BLAT_BIN_CACHE.mkdir(parents=True, exist_ok=True)
    out = {}
    for rel in BLAT_BINARIES:
        name = rel.split("/")[-1]
        dest = config.BLAT_BIN_CACHE / name
        if not dest.exists() or dest.stat().st_size < 100_000:
            last_error = None
            for mirror in config.UCSC_EXE_MIRRORS:
                try:
                    _probe(f"{mirror}/{rel}")
                    download(f"{mirror}/{rel}", dest, st, "blat", f"BLAT {name}")
                    last_error = None
                    break
                except Exception as exc:                      # noqa: BLE001
                    last_error = exc
                    st.set("blat", f"{mirror} unavailable for {name}; trying next mirror")
            if last_error:
                raise SetupError(f"Could not download {name} from UCSC: {last_error}")
        out[name] = dest
    return out


def _probe(url: str) -> None:
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": config.USER_AGENT})
    with urllib.request.urlopen(req, timeout=20):
        pass


def c_blat(m: dict, st: Status) -> None:
    if species.blat_installed(m):
        st.set("blat", "BLAT index already built", 100, state="done")
        return
    fasta = (m.get("genome") or {}).get("fasta")
    if not fasta or not Path(fasta).exists():
        raise SetupError("Install the genome component first.")
    st.set("blat", "Preparing BLAT" if wsl.NATIVE else "Connecting to WSL")
    distro = wsl.distro()
    if not distro:
        raise SetupError("WSL has no Linux distribution. Install one with "
                         "'wsl --install -d Ubuntu' and try again.")
    s = config.settings()
    root = wsl.linux_path(s["blat_dir"])
    bindir, sdir = f"{root}/bin", f"{root}/{m['name']}"
    step = int(s.get("blat_step_size") or 5)

    cached = _blat_binaries_to_cache(st)
    st.set("blat", f"Installing BLAT binaries into {bindir}" if wsl.NATIVE
           else f"Installing BLAT binaries into {distro}:{bindir}")
    copies = " && ".join(f"cp -f {wsl.quote(wsl.to_wsl(p))} {wsl.quote(bindir)}/{name}"
                         for name, p in cached.items())
    wsl.bash(f"mkdir -p {wsl.quote(bindir)} {wsl.quote(sdir)} && {copies} && "
             f"chmod +x {wsl.quote(bindir)}/*", timeout=600, check=True)
    usage = wsl.run([f"{bindir}/gfServer"], timeout=60)
    version = (usage.stdout or usage.stderr).splitlines()[:1]

    twobit = f"{sdir}/genome.2bit"
    have = wsl.bash(f"test -s {wsl.quote(twobit)} && echo yes", timeout=60).stdout.strip()
    if have != "yes":
        st.set("blat", "Converting genome to 2bit (a few minutes)...")
        proc = wsl.run([f"{bindir}/faToTwoBit", wsl.to_wsl(fasta), twobit + ".part"])
        if proc.returncode != 0:
            raise SetupError(f"faToTwoBit failed: {proc.stderr.strip()[:400]}")
        wsl.run(["mv", "-f", twobit + ".part", twobit], check=True)

    index = f"{sdir}/genome.untrans.gfidx"
    marker = f"{index}.complete"
    have = wsl.bash(f"test -s {wsl.quote(index)} && test -e {wsl.quote(marker)} && echo yes",
                    timeout=60).stdout.strip()
    if have != "yes":
        st.set("blat", f"Building gfServer index (stepSize={step}); large genomes take "
                       "10-20 minutes...")
        # gfServer requires the index to be named <genome>.untrans.gfidx and records
        # sequence file names as given (gfClient resolves them against its seqDir), so
        # build in place from the species directory and mark completion separately.
        wsl.run(["rm", "-f", index, marker], check=True)
        proc = wsl.run([f"{bindir}/gfServer", "index", f"-stepSize={step}",
                        "genome.untrans.gfidx", "genome.2bit"], cwd=sdir)
        if proc.returncode != 0:
            wsl.run(["rm", "-f", index])
            raise SetupError(f"gfServer index failed: {(proc.stderr or proc.stdout)[:400]}")
        wsl.run(["touch", marker], check=True)

    sizes = wsl.bash(f"du -h {wsl.quote(twobit)} {wsl.quote(index)} | cut -f1 | paste -sd/",
                     timeout=60).stdout.strip()
    used = {(species.load(n) or {}).get("blat", {}).get("port")
            for n in species.installed_names() if n != m["name"]}
    port = (m.get("blat") or {}).get("port") or next(
        p for p in range(config.BLAT_BASE_PORT, config.BLAT_BASE_PORT + 500) if p not in used)
    species.update(m["name"], blat={
        "distro": distro, "dir": sdir, "bindir": bindir, "twobit": "genome.2bit",
        "index": "genome.untrans.gfidx", "port": port, "step_size": step,
        "version": version[0] if version else "", "ready": True, "built_at": time.time()})
    st.set("blat", f"BLAT ready{'' if wsl.NATIVE else ' in WSL'} ({sizes}), port {port}",
           100, state="done")


RUNNERS = {
    "karyotype": c_karyotype, "genome": c_genome, "blastdb": c_blastdb,
    "cdna": lambda m, st: c_transcripts(m, st, "cdna"),
    "ncrna": lambda m, st: c_transcripts(m, st, "ncrna"),
    "pep": lambda m, st: c_transcripts(m, st, "pep"),
    "gtf": c_gtf, "xrefs": c_xrefs, "refseq": c_refseq, "repeats": c_repeats, "blat": c_blat,
}


def install(name: str, components: list[str]) -> int:
    config.ensure_dirs()
    if name == "homo_sapiens":
        species.ensure_human()
    if not species.load(name):
        info = ensembl_species_info(name)
        if not info:
            print(f"Unknown Ensembl species '{name}'.", file=sys.stderr)
            return 2
        species.save({
            "name": name, "display_name": info.get("display_name") or name,
            "scientific_name": name.replace("_", " ").capitalize(),
            "common_name": info.get("common_name"), "assembly": info.get("assembly"),
            "taxon_id": int(info["taxon_id"]) if info.get("taxon_id") else None,
            "release": info.get("release"), "division": info.get("division"),
            "created_at": time.time()})
    else:
        info = ensembl_species_info(name)
        m = species.load(name) or {}
        if not m.get("release") and info.get("release"):
            species.update(name, release=info["release"])

    order = [c for c in species.INSTALL_ORDER if c in components]
    st = Status(name, order)
    failures = []
    for comp in order:
        manifest = species.load(name) or {"name": name}
        try:
            RUNNERS[comp](manifest, st)
        except Exception as exc:                               # noqa: BLE001
            failures.append(comp)
            st.set(comp, f"Failed: {exc}", state="error")
    st.data.update(done=True, component=None, pct=100,
                   error=(f"{len(failures)} component(s) failed: {', '.join(failures)}"
                          if failures else None),
                   message="Finished" if not failures else "Finished with errors")
    st.write()
    free = shutil.disk_usage(config.DATA_DIR).free
    print(f"Done. {human(free)} free on the data drive.", flush=True)
    return 1 if failures else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("species", help="Ensembl production name, e.g. homo_sapiens")
    ap.add_argument("--components", default=",".join(species.INSTALL_ORDER))
    ap.add_argument("--status", action="store_true")
    args = ap.parse_args()
    if args.status:
        print(json.dumps(read_status(args.species), indent=2))
        m = species.load(args.species)
        print(json.dumps(species.summary(m) if m else {}, indent=2))
        return 0
    comps = [c.strip() for c in args.components.split(",") if c.strip()]
    bad = [c for c in comps if c not in RUNNERS]
    if bad:
        print(f"Unknown component(s): {', '.join(bad)}", file=sys.stderr)
        return 2
    return install(args.species, comps)


if __name__ == "__main__":
    raise SystemExit(main())
