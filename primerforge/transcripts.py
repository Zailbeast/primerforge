"""Transcript identifiers: MANE Select, RefSeq accessions and Ensembl IDs.

Primer design is usually specified against one agreed transcript, so the tool needs
to answer "which transcript is the MANE Select one, and what is its NM accession?".
Two small tables give that offline:

  MANE summary (NCBI/EMBL-EBI)  one matched RefSeq/Ensembl transcript per gene,
                                flagged MANE Select or MANE Plus Clinical
  Ensembl cross-references      every ENST/ENSP with its RefSeq NM/NP accession

RefSeq transcript FASTA headers add each accession's gene symbol and description.
Everything lives in data/species/<name>/transcripts.sqlite next to the GTF index.
"""
from __future__ import annotations

import gzip
import re
import sqlite3
import threading
from pathlib import Path

from . import annotation

SCHEMA = """
CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT);
CREATE TABLE IF NOT EXISTS xref (
    transcript_id TEXT, protein_id TEXT, gene_id TEXT,
    refseq_mrna TEXT, refseq_peptide TEXT);
CREATE TABLE IF NOT EXISTS mane (
    symbol TEXT, name TEXT, hgnc TEXT, ncbi_gene TEXT, ensembl_gene TEXT,
    ensembl_nuc TEXT, ensembl_nuc_v TEXT, ensembl_prot TEXT,
    refseq_nuc TEXT, refseq_nuc_v TEXT, refseq_prot TEXT,
    status TEXT, chrom TEXT, start INTEGER, end INTEGER, strand INTEGER);
CREATE TABLE IF NOT EXISTS refseq_info (
    accession TEXT PRIMARY KEY, version TEXT, symbol TEXT, description TEXT,
    moltype TEXT, length INTEGER);
CREATE INDEX IF NOT EXISTS idx_xref_tx ON xref(transcript_id);
CREATE INDEX IF NOT EXISTS idx_xref_nm ON xref(refseq_mrna);
CREATE INDEX IF NOT EXISTS idx_mane_enst ON mane(ensembl_nuc);
CREATE INDEX IF NOT EXISTS idx_mane_nm ON mane(refseq_nuc);
CREATE INDEX IF NOT EXISTS idx_mane_symbol ON mane(symbol);
CREATE INDEX IF NOT EXISTS idx_refseq_symbol ON refseq_info(symbol);
"""

MANE_SELECT = "MANE Select"
MANE_PLUS = "MANE Plus Clinical"
_ACC_RE = re.compile(r"^([A-Z]{2}_\d+|ENS[A-Z]*[GTP]\d+)(?:\.(\d+))?$", re.I)
# ">NM_001008216.2 Homo sapiens UDP-galactose-4-epimerase (GALE), transcript variant 2, mRNA"
_HEADER_RE = re.compile(r"^(\S+)\s+(.*?)\s*\(([^()]+)\),\s*([^,]+(?:,\s*[^,]+)*)$")


def db_path(species_dir: Path) -> Path:
    return Path(species_dir) / "transcripts.sqlite"


def split_version(identifier: str) -> tuple[str, str]:
    m = _ACC_RE.match((identifier or "").strip())
    if not m:
        base, _, version = (identifier or "").strip().partition(".")
        return base, version
    return m.group(1), m.group(2) or ""


# ---------------------------------------------------------------------------
# Building
# ---------------------------------------------------------------------------
def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if str(path).endswith(".gz") \
        else open(path, "r", encoding="utf-8", errors="replace")


def build_xrefs(sqlite_path: Path, ensembl_tsv: Path | None, mane_summary: Path | None,
                progress=None) -> dict:
    """Load Ensembl's RefSeq cross-reference table and the MANE summary."""
    counts = {"xrefs": 0, "mane": 0}
    with sqlite3.connect(sqlite_path) as conn:
        conn.executescript(SCHEMA)
        if ensembl_tsv:
            conn.execute("DELETE FROM xref")
            by_tx: dict[str, dict] = {}
            with _open_text(ensembl_tsv) as fh:
                next(fh, None)
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 5 or not f[1]:
                        continue
                    entry = by_tx.setdefault(f[1], {"gene": f[0], "protein": f[2],
                                                    "mrna": None, "peptide": None})
                    if f[4] == "RefSeq_mRNA" and not entry["mrna"]:
                        entry["mrna"] = f[3]
                    elif f[4] == "RefSeq_peptide" and not entry["peptide"]:
                        entry["peptide"] = f[3]
            conn.executemany("INSERT INTO xref VALUES (?,?,?,?,?)",
                             [(tx, e["protein"], e["gene"], e["mrna"], e["peptide"])
                              for tx, e in by_tx.items()])
            counts["xrefs"] = len(by_tx)
            if progress:
                progress(counts)
        if mane_summary:
            conn.execute("DELETE FROM mane")
            rows = []
            with _open_text(mane_summary) as fh:
                header = next(fh, "").lstrip("#").rstrip("\n").split("\t")
                col = {name: i for i, name in enumerate(header)}
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) < len(header):
                        continue
                    enst, enst_v = split_version(f[col["Ensembl_nuc"]])
                    nm, nm_v = split_version(f[col["RefSeq_nuc"]])
                    rows.append((
                        f[col["symbol"]], f[col["name"]], f[col["HGNC_ID"]],
                        f[col["NCBI_GeneID"]].replace("GeneID:", ""),
                        split_version(f[col["Ensembl_Gene"]])[0],
                        enst, f[col["Ensembl_nuc"]], f[col["Ensembl_prot"]],
                        nm, f[col["RefSeq_nuc"]], f[col["RefSeq_prot"]],
                        f[col["MANE_status"]], f[col["GRCh38_chr"]],
                        int(f[col["chr_start"]] or 0), int(f[col["chr_end"]] or 0),
                        -1 if f[col["chr_strand"]] == "-" else 1))
            conn.executemany("INSERT INTO mane VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", rows)
            counts["mane"] = len(rows)
        conn.executescript(SCHEMA)
        conn.commit()
    reset()
    return counts


def build_refseq_info(sqlite_path: Path, fastas: list[Path], progress=None) -> int:
    """Index accession -> gene symbol, description and length from RefSeq FASTA headers."""
    with sqlite3.connect(sqlite_path) as conn:
        conn.executescript(SCHEMA)
        conn.execute("DELETE FROM refseq_info")
        rows, total = [], 0
        for path in fastas:
            acc = desc = symbol = moltype = ""
            length = 0
            with _open_text(path) as fh:
                for line in fh:
                    if line.startswith(">"):
                        if acc:
                            rows.append(_refseq_row(acc, symbol, desc, moltype, length))
                        head = line[1:].rstrip("\n")
                        acc, symbol, desc, moltype, length = head.split(" ", 1)[0], "", head, "", 0
                        m = _HEADER_RE.match(head)
                        if m:
                            symbol, desc, moltype = m.group(3), m.group(2), m.group(4)
                        else:
                            parts = head.split(" ", 1)
                            desc = parts[1] if len(parts) > 1 else ""
                    else:
                        length += len(line.strip())
                    if len(rows) >= 20000:
                        conn.executemany(
                            "INSERT OR REPLACE INTO refseq_info VALUES (?,?,?,?,?,?)", rows)
                        total += len(rows)
                        rows = []
                        if progress:
                            progress(total)
            if acc:
                rows.append(_refseq_row(acc, symbol, desc, moltype, length))
        conn.executemany("INSERT OR REPLACE INTO refseq_info VALUES (?,?,?,?,?,?)", rows)
        total += len(rows)
        conn.executescript(SCHEMA)
        conn.commit()
    reset()
    return total


def _refseq_row(acc: str, symbol: str, desc: str, moltype: str, length: int):
    base, version = split_version(acc)
    return (base, version, symbol, desc, moltype, length)


def mane_accessions(sqlite_path: Path, status: str | None = MANE_SELECT) -> list[str]:
    with sqlite3.connect(sqlite_path) as conn:
        sql = "SELECT refseq_nuc_v FROM mane"
        args: tuple = ()
        if status:
            sql += " WHERE status=?"
            args = (status,)
        return [r[0] for r in conn.execute(sql, args) if r[0]]


# ---------------------------------------------------------------------------
# Queries
# ---------------------------------------------------------------------------
_conns: dict[str, sqlite3.Connection] = {}
_lock = threading.Lock()


def _db(manifest: dict) -> sqlite3.Connection | None:
    path = (manifest.get("transcripts") or {}).get("sqlite")
    if not path or not Path(path).exists():
        return None
    key = f"{path}:{threading.get_ident()}"
    with _lock:
        conn = _conns.get(key)
        if conn is None:
            conn = sqlite3.connect(f"file:{path}?mode=ro", uri=True, check_same_thread=False)
            conn.row_factory = sqlite3.Row
            _conns[key] = conn
    return conn


def reset() -> None:
    with _lock:
        for c in _conns.values():
            try:
                c.close()
            except sqlite3.Error:
                pass
        _conns.clear()


def available(manifest: dict) -> bool:
    return _db(manifest) is not None


def mane_for(manifest: dict, identifier: str) -> dict | None:
    """MANE record for an Ensembl or RefSeq transcript/protein accession."""
    conn = _db(manifest)
    if conn is None or not identifier:
        return None
    base = split_version(identifier)[0]
    row = conn.execute(
        "SELECT * FROM mane WHERE ensembl_nuc=? OR refseq_nuc=? OR ensembl_prot LIKE ?"
        " OR refseq_prot LIKE ? LIMIT 1", (base, base, base + ".%", base + ".%")).fetchone()
    return dict(row) if row else None


def refseq_for(manifest: dict, ensembl_id: str) -> dict | None:
    """RefSeq NM/NP accessions for an Ensembl transcript (or its protein)."""
    conn = _db(manifest)
    if conn is None or not ensembl_id:
        return None
    base = split_version(ensembl_id)[0]
    row = conn.execute("SELECT * FROM xref WHERE transcript_id=? OR protein_id=? LIMIT 1",
                       (base, base)).fetchone()
    return dict(row) if row else None


def ensembl_for(manifest: dict, refseq_id: str) -> str | None:
    """Ensembl transcript matching a RefSeq accession, preferring the MANE pairing."""
    conn = _db(manifest)
    if conn is None or not refseq_id:
        return None
    base = split_version(refseq_id)[0]
    row = conn.execute("SELECT ensembl_nuc FROM mane WHERE refseq_nuc=? OR refseq_prot LIKE ?",
                       (base, base + ".%")).fetchone()
    if row:
        return row["ensembl_nuc"]
    row = conn.execute("SELECT transcript_id FROM xref WHERE refseq_mrna=? OR refseq_peptide=?"
                       " LIMIT 1", (base, base)).fetchone()
    return row["transcript_id"] if row else None


def span_for(manifest: dict, accession: str) -> dict | None:
    """Gene span from the MANE table, as a fallback location for a RefSeq accession."""
    from .variants import NC_ACCESSIONS

    m = mane_for(manifest, accession)
    if not m or not m.get("chrom"):
        return None
    by_nc = {nc: chrom for chrom, nc in NC_ACCESSIONS.items()}
    chrom = by_nc.get(m["chrom"], m["chrom"])
    return {"chrom": chrom, "start": m["start"], "end": m["end"], "strand": m["strand"]}


def refseq_info(manifest: dict, accession: str) -> dict | None:
    conn = _db(manifest)
    if conn is None or not accession:
        return None
    row = conn.execute("SELECT * FROM refseq_info WHERE accession=?",
                       (split_version(accession)[0],)).fetchone()
    return dict(row) if row else None


def annotate(manifest: dict, identifier: str) -> dict:
    """Identifier labels for a hit: matching accession and MANE status."""
    out: dict = {}
    if not identifier:
        return out
    base = split_version(identifier)[0]
    if base.upper().startswith("ENS"):
        x = refseq_for(manifest, base)
        if x:
            out["refseq_mrna"] = x.get("refseq_mrna")
            out["refseq_peptide"] = x.get("refseq_peptide")
    else:
        info = refseq_info(manifest, base)
        if info:
            out["symbol"] = info["symbol"]
            out["description"] = info["description"]
            out["moltype"] = info["moltype"]
        enst = ensembl_for(manifest, base)
        if enst:
            out["ensembl_transcript"] = enst
    m = mane_for(manifest, base)
    if m:
        out["mane"] = m["status"]
        out["mane_refseq"] = m["refseq_nuc_v"]
        out["mane_ensembl"] = m["ensembl_nuc_v"]
        out.setdefault("symbol", m["symbol"])
    return out


# ---------------------------------------------------------------------------
# Gene / transcript search
# ---------------------------------------------------------------------------
def search(manifest: dict, query: str, limit: int = 25) -> list[dict]:
    """Find genes by symbol, name, or any transcript/protein accession."""
    q = (query or "").strip()
    if len(q) < 2:
        return []
    conn = _db(manifest)
    gene_ids: list[str] = []
    seen = set()

    def add(gene_id: str | None):
        if gene_id and gene_id not in seen:
            seen.add(gene_id)
            gene_ids.append(gene_id)

    base = split_version(q)[0]
    if conn is not None:
        for row in conn.execute(
                "SELECT ensembl_gene, symbol FROM mane WHERE symbol=? COLLATE NOCASE"
                " OR ensembl_nuc=? OR refseq_nuc=? OR ensembl_gene=? OR ncbi_gene=?"
                " OR hgnc=? COLLATE NOCASE LIMIT ?", (q, base, base, base, q, q, limit)):
            add(row["ensembl_gene"])
        if not gene_ids:
            row = conn.execute("SELECT transcript_id FROM xref WHERE refseq_mrna=?"
                               " OR refseq_peptide=? LIMIT 1", (base, base)).fetchone()
            if row:
                tx = annotation.transcript(manifest, row["transcript_id"])
                if tx:
                    add(tx["gene_id"])
    genes = _genes_by_query(manifest, q, base, limit)
    for g in genes:
        add(g["id"])
    return [detail for gid in gene_ids[:limit]
            if (detail := gene_detail(manifest, gid)) is not None]


def _genes_by_query(manifest: dict, q: str, base: str, limit: int) -> list[dict]:
    conn = annotation.connection(manifest)
    if conn is None:
        return []
    rows = conn.execute(
        "SELECT * FROM genes WHERE name=? COLLATE NOCASE OR id=? LIMIT ?",
        (q, base, limit)).fetchall()
    if not rows:
        tx = conn.execute("SELECT gene_id FROM transcripts WHERE id=? OR protein_id=?"
                          " OR name=? COLLATE NOCASE LIMIT 1", (base, base, q)).fetchone()
        if tx:
            rows = conn.execute("SELECT * FROM genes WHERE id=?", (tx["gene_id"],)).fetchall()
    if not rows:
        rows = conn.execute(
            "SELECT * FROM genes WHERE name LIKE ? COLLATE NOCASE ORDER BY length(name) LIMIT ?",
            (q + "%", limit)).fetchall()
    return [dict(r) for r in rows]


def gene_detail(manifest: dict, gene_id: str) -> dict | None:
    """A gene with all of its transcripts, their accessions and MANE status."""
    conn = annotation.connection(manifest)
    if conn is None:
        return None
    gene = conn.execute("SELECT * FROM genes WHERE id=?",
                        (split_version(gene_id)[0],)).fetchone()
    if not gene:
        return None
    rows = conn.execute(
        "SELECT t.*, (SELECT COUNT(*) FROM exons x WHERE x.transcript_rid=t.rid) AS n_exons,"
        " (SELECT SUM(x.end - x.start + 1) FROM exons x WHERE x.transcript_rid=t.rid) AS cdna_len,"
        " (SELECT SUM(c.end - c.start + 1) FROM cds c WHERE c.transcript_rid=t.rid) AS cds_len"
        " FROM transcripts t WHERE t.gene_id=? ORDER BY t.canonical DESC, cdna_len DESC",
        (gene["id"],)).fetchall()
    transcripts = []
    for r in rows:
        tx = dict(r)
        ids = annotate(manifest, tx["id"])
        cds_len = tx.get("cds_len") or 0
        transcripts.append({
            "id": tx["id"], "version": tx["version"], "name": tx["name"],
            "biotype": tx["biotype"], "chrom": tx["chrom"], "start": tx["start"],
            "end": tx["end"], "strand": tx["strand"], "protein_id": tx["protein_id"],
            "canonical": bool(tx["canonical"]), "n_exons": tx["n_exons"],
            "cdna_len": tx["cdna_len"], "cds_len": cds_len,
            # Ensembl's CDS features exclude the stop codon, so every codon is a residue.
            "protein_len": cds_len // 3 if cds_len else 0,
            "refseq_mrna": ids.get("refseq_mrna"), "refseq_peptide": ids.get("refseq_peptide"),
            "mane": ids.get("mane"),
        })
    order = {MANE_SELECT: 0, MANE_PLUS: 1}
    transcripts.sort(key=lambda t: (order.get(t["mane"], 2), not t["canonical"],
                                    -(t["cdna_len"] or 0)))
    m = None
    for t in transcripts:
        if t["mane"] == MANE_SELECT:
            m = mane_for(manifest, t["id"])
            break
    return {"id": gene["id"], "version": gene["version"], "name": gene["name"],
            "biotype": gene["biotype"], "chrom": gene["chrom"], "start": gene["start"],
            "end": gene["end"], "strand": gene["strand"], "transcripts": transcripts,
            "hgnc": (m or {}).get("hgnc"), "ncbi_gene": (m or {}).get("ncbi_gene"),
            "description": (m or {}).get("name")}
