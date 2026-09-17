"""Gene annotation for BLAST/BLAT hits.

Ensembl's results table lists the genes each hit overlaps, maps cDNA and protein
hits back onto the genome, and marks exons and variants in the sequence views.
With the species GTF installed all of that is answered locally from an SQLite
index (an R*Tree over genes and exons); otherwise the Ensembl REST API is used
and its answers cached.
"""
from __future__ import annotations

import gzip
import re
import sqlite3
import threading
import urllib.parse
from contextlib import closing
from pathlib import Path

from . import config, net, store

SCHEMA = """
CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT);
CREATE TABLE IF NOT EXISTS chroms (idx INTEGER PRIMARY KEY, name TEXT UNIQUE);
CREATE TABLE IF NOT EXISTS genes (
    rid INTEGER PRIMARY KEY, id TEXT, version TEXT, name TEXT, biotype TEXT,
    chrom TEXT, start INTEGER, end INTEGER, strand INTEGER);
CREATE TABLE IF NOT EXISTS transcripts (
    rid INTEGER PRIMARY KEY, id TEXT, version TEXT, gene_id TEXT, gene_name TEXT,
    name TEXT, biotype TEXT, chrom TEXT, start INTEGER, end INTEGER, strand INTEGER,
    protein_id TEXT, canonical INTEGER DEFAULT 0);
CREATE TABLE IF NOT EXISTS exons (
    rid INTEGER PRIMARY KEY, transcript_rid INTEGER, number INTEGER,
    start INTEGER, end INTEGER);
CREATE TABLE IF NOT EXISTS cds (
    transcript_rid INTEGER, start INTEGER, end INTEGER, phase INTEGER);
CREATE VIRTUAL TABLE IF NOT EXISTS gene_rt USING rtree_i32(rid, c0, c1, s, e);
CREATE VIRTUAL TABLE IF NOT EXISTS exon_rt USING rtree_i32(rid, c0, c1, s, e);
"""
INDEXES = """
CREATE INDEX IF NOT EXISTS idx_tx_id ON transcripts(id);
CREATE INDEX IF NOT EXISTS idx_tx_protein ON transcripts(protein_id);
CREATE INDEX IF NOT EXISTS idx_gene_id ON genes(id);
CREATE INDEX IF NOT EXISTS idx_exon_tx ON exons(transcript_rid);
CREATE INDEX IF NOT EXISTS idx_cds_tx ON cds(transcript_rid);
"""

_ATTR_RE = re.compile(r'(\w+) "([^"]*)"')


class AnnotationError(RuntimeError):
    pass


# ---------------------------------------------------------------------------
# Building the index
# ---------------------------------------------------------------------------
def build(gtf_path: Path, sqlite_path: Path, progress=None) -> dict:
    """Stream an Ensembl GTF (optionally gzipped) into an SQLite annotation index."""
    tmp = sqlite_path.with_suffix(".building")
    tmp.unlink(missing_ok=True)
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    counts = {"genes": 0, "transcripts": 0, "exons": 0, "cds": 0}
    chrom_idx: dict[str, int] = {}
    tx_rid: dict[str, int] = {}
    protein_by_rid: dict[int, str] = {}     # protein ids only appear on CDS lines

    with closing(sqlite3.connect(tmp)) as conn:
        conn.execute("PRAGMA journal_mode=OFF")
        conn.execute("PRAGMA synchronous=OFF")
        conn.executescript(SCHEMA)
        genes, gene_rt, txs, exons, exon_rt, cds = [], [], [], [], [], []
        next_gene = next_tx = next_exon = 1

        def flush(force=False):
            nonlocal genes, gene_rt, txs, exons, exon_rt, cds
            if not force and len(exons) < 50_000:
                return
            conn.executemany("INSERT INTO genes VALUES (?,?,?,?,?,?,?,?,?)", genes)
            conn.executemany("INSERT INTO gene_rt VALUES (?,?,?,?,?)", gene_rt)
            conn.executemany("INSERT INTO transcripts VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)", txs)
            conn.executemany("INSERT INTO exons VALUES (?,?,?,?,?)", exons)
            conn.executemany("INSERT INTO exon_rt VALUES (?,?,?,?,?)", exon_rt)
            conn.executemany("INSERT INTO cds VALUES (?,?,?,?)", cds)
            genes, gene_rt, txs, exons, exon_rt, cds = [], [], [], [], [], []
            if progress:
                progress(counts)

        with opener(gtf_path, "rt", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9:
                    continue
                chrom, _src, feature, start, end, _score, strand, phase, attrs = cols
                if feature not in ("gene", "transcript", "exon", "CDS"):
                    continue
                a = dict(_ATTR_RE.findall(attrs))
                s, e = int(start), int(end)
                st = 1 if strand == "+" else -1
                if chrom not in chrom_idx:
                    chrom_idx[chrom] = len(chrom_idx) + 1
                ci = chrom_idx[chrom]

                if feature == "gene":
                    genes.append((next_gene, a.get("gene_id"), a.get("gene_version"),
                                  a.get("gene_name") or a.get("gene_id"),
                                  a.get("gene_biotype"), chrom, s, e, st))
                    gene_rt.append((next_gene, ci, ci, s, e))
                    next_gene += 1
                    counts["genes"] += 1
                elif feature == "transcript":
                    tid = a.get("transcript_id")
                    tx_rid[tid] = next_tx
                    txs.append((next_tx, tid, a.get("transcript_version"), a.get("gene_id"),
                                a.get("gene_name") or a.get("gene_id"),
                                a.get("transcript_name"), a.get("transcript_biotype"),
                                chrom, s, e, st, None,
                                1 if 'tag "Ensembl_canonical"' in attrs else 0))
                    next_tx += 1
                    counts["transcripts"] += 1
                elif feature == "exon":
                    rid = tx_rid.get(a.get("transcript_id"))
                    if rid is None:
                        continue
                    exons.append((next_exon, rid, int(a.get("exon_number") or 0), s, e))
                    exon_rt.append((next_exon, ci, ci, s, e))
                    next_exon += 1
                    counts["exons"] += 1
                else:  # CDS
                    rid = tx_rid.get(a.get("transcript_id"))
                    if rid is None:
                        continue
                    cds.append((rid, s, e, int(phase) if phase.isdigit() else 0))
                    if a.get("protein_id"):
                        protein_by_rid[rid] = a["protein_id"]
                    counts["cds"] += 1
                flush()
        flush(force=True)

        conn.executemany("UPDATE transcripts SET protein_id=? WHERE rid=?",
                         [(pid, rid) for rid, pid in protein_by_rid.items()])
        conn.executemany("INSERT INTO chroms VALUES (?,?)",
                         [(i, n) for n, i in chrom_idx.items()])
        conn.executemany("INSERT INTO meta VALUES (?,?)",
                         [(k, str(v)) for k, v in counts.items()] +
                         [("source", Path(gtf_path).name)])
        conn.executescript(INDEXES)
        conn.commit()
    tmp.replace(sqlite_path)
    return counts


# ---------------------------------------------------------------------------
# Queries
# ---------------------------------------------------------------------------
_conns: dict[str, sqlite3.Connection] = {}
_conn_lock = threading.Lock()


def _db(manifest: dict) -> sqlite3.Connection | None:
    path = (manifest.get("gtf") or {}).get("sqlite")
    if not path or not Path(path).exists():
        return None
    key = f"{path}:{threading.get_ident()}"
    with _conn_lock:
        conn = _conns.get(key)
        if conn is None:
            conn = sqlite3.connect(f"file:{path}?mode=ro", uri=True, check_same_thread=False)
            conn.row_factory = sqlite3.Row
            _conns[key] = conn
    return conn


def reset() -> None:
    with _conn_lock:
        for c in _conns.values():
            try:
                c.close()
            except sqlite3.Error:
                pass
        _conns.clear()


def has_local(manifest: dict) -> bool:
    return _db(manifest) is not None


def connection(manifest: dict) -> sqlite3.Connection | None:
    """Read-only handle on the annotation index, for callers that query it directly."""
    return _db(manifest)


def _chrom_idx(conn: sqlite3.Connection, chrom: str) -> int | None:
    for cand in (chrom, chrom.removeprefix("chr"), "chr" + chrom):
        row = conn.execute("SELECT idx FROM chroms WHERE name=?", (cand,)).fetchone()
        if row:
            return row["idx"]
    return None


def genes_in_region(manifest: dict, chrom: str, start: int, end: int) -> list[dict]:
    conn = _db(manifest)
    if conn is None:
        return _rest_genes(manifest, chrom, start, end)
    ci = _chrom_idx(conn, chrom)
    if ci is None:
        return []
    rows = conn.execute(
        "SELECT g.* FROM gene_rt r JOIN genes g ON g.rid=r.rid "
        "WHERE r.c0<=? AND r.c1>=? AND r.s<=? AND r.e>=? ORDER BY g.start",
        (ci, ci, end, start)).fetchall()
    return [{"id": r["id"], "name": r["name"], "biotype": r["biotype"], "chrom": r["chrom"],
             "start": r["start"], "end": r["end"], "strand": r["strand"]} for r in rows]


def exons_in_region(manifest: dict, chrom: str, start: int, end: int,
                    cached_only: bool = False) -> list[dict]:
    conn = _db(manifest)
    if conn is None:
        return _rest_exons(manifest, chrom, start, end, cached_only)
    ci = _chrom_idx(conn, chrom)
    if ci is None:
        return []
    rows = conn.execute(
        "SELECT x.number, x.start, x.end, t.id AS transcript_id, t.name AS transcript_name, "
        "t.gene_name, t.strand, t.biotype, t.canonical FROM exon_rt r "
        "JOIN exons x ON x.rid=r.rid JOIN transcripts t ON t.rid=x.transcript_rid "
        "WHERE r.c0<=? AND r.c1>=? AND r.s<=? AND r.e>=? ORDER BY x.start",
        (ci, ci, end, start)).fetchall()
    return [dict(r) for r in rows]


def transcript(manifest: dict, stable_id: str) -> dict | None:
    """Structure of a transcript (or of the transcript encoding a protein)."""
    conn = _db(manifest)
    if conn is None:
        return None
    base = stable_id.split(".")[0]
    row = conn.execute("SELECT * FROM transcripts WHERE id=?", (base,)).fetchone() or \
        conn.execute("SELECT * FROM transcripts WHERE protein_id=?", (base,)).fetchone()
    if not row:
        return None
    exons = conn.execute("SELECT number, start, end FROM exons WHERE transcript_rid=? "
                         "ORDER BY number", (row["rid"],)).fetchall()
    cds = conn.execute("SELECT start, end, phase FROM cds WHERE transcript_rid=?",
                       (row["rid"],)).fetchall()
    strand = row["strand"]
    cds_sorted = sorted((dict(c) for c in cds), key=lambda c: c["start"] * strand)
    return {"id": row["id"], "gene_id": row["gene_id"], "gene_name": row["gene_name"],
            "name": row["name"], "biotype": row["biotype"], "chrom": row["chrom"],
            "start": row["start"], "end": row["end"], "strand": strand,
            "protein_id": row["protein_id"],
            "exons": [dict(e) for e in exons], "cds": cds_sorted}


def _walk(segments: list[tuple[int, int]], strand: int, start: int, end: int,
          offset: int = 0) -> list[tuple[int, int]]:
    """Map a 1-based range along spliced segments (in transcription order) to genome."""
    out = []
    pos = 1 - offset
    for s, e in segments:
        length = e - s + 1
        lo, hi = max(start, pos), min(end, pos + length - 1)
        if lo <= hi:
            if strand > 0:
                out.append((s + (lo - pos), s + (hi - pos)))
            else:
                out.append((e - (hi - pos), e - (lo - pos)))
        pos += length
        if pos > end:
            break
    return out


def map_to_genome(manifest: dict, target_id: str, kind: str, start: int, end: int,
                  allow_rest: bool = True) -> dict | None:
    """Genomic span of a cDNA/ncRNA ("cdna") or protein ("pep") hit range.

    With a local annotation index installed the mapping is answered from it; the
    REST fallback is for species without one, and callers mapping many hits should
    turn it off so one unmappable accession cannot stall a job.
    """
    tx = transcript(manifest, target_id)
    if tx:
        if kind == "pep":
            if not tx["cds"]:
                return None
            segs = [(c["start"], c["end"]) for c in tx["cds"]]
            phase = tx["cds"][0]["phase"] or 0
            pieces = _walk(segs, tx["strand"], (start - 1) * 3 + 1, end * 3, offset=phase)
        else:
            segs = [(e["start"], e["end"]) for e in tx["exons"]]
            pieces = _walk(segs, tx["strand"], start, end)
        if not pieces:
            return None
        return {"chrom": tx["chrom"], "start": min(p[0] for p in pieces),
                "end": max(p[1] for p in pieces), "strand": tx["strand"],
                "blocks": sorted(pieces), "gene_id": tx["gene_id"],
                "gene_name": tx["gene_name"], "transcript_id": tx["id"]}
    if not allow_rest:
        return None
    return _rest_map(manifest, target_id, kind, start, end)


# ---------------------------------------------------------------------------
# Ensembl REST fallbacks (cached)
# ---------------------------------------------------------------------------
def _rest_base(manifest: dict) -> str:
    return config.ENSEMBL_REST_GRCH37 if manifest.get("assembly") == "GRCh37" \
        else config.ENSEMBL_REST


class NotCached(LookupError):
    """Raised by cache-only lookups so a page can render now and fetch later."""


def _rest_get(manifest: dict, path: str, max_age: float = 30 * 86400,
              cached_only: bool = False, timeout: int = 30, attempts: int = 2):
    url = _rest_base(manifest) + path
    key = f"rest:{url}"
    hit = store.cache_get(key, max_age=max_age)
    if hit is not None:
        return hit
    if cached_only:
        raise NotCached(url)
    try:
        data = net.fetch_json(url, attempts=attempts, timeout=timeout)
    except net.HttpError:
        return None
    store.cache_put(key, data)
    return data


def _region(chrom: str, start: int, end: int) -> str:
    return f"{urllib.parse.quote(chrom.removeprefix('chr'))}:{max(1, start)}-{end}"


def _rest_genes(manifest, chrom, start, end) -> list[dict]:
    if end - start > 5_000_000:
        return []
    data = _rest_get(manifest, f"/overlap/region/{manifest['name']}/"
                               f"{_region(chrom, start, end)}?feature=gene;"
                               "content-type=application/json") or []
    return [{"id": g.get("id"), "name": g.get("external_name") or g.get("id"),
             "biotype": g.get("biotype"), "chrom": chrom, "start": g.get("start"),
             "end": g.get("end"), "strand": g.get("strand")}
            for g in data if isinstance(g, dict)]


def _rest_exons(manifest, chrom, start, end, cached_only=False) -> list[dict]:
    if end - start > 1_000_000:
        return []
    base = f"/overlap/region/{manifest['name']}/{_region(chrom, start, end)}"
    exons = _rest_get(manifest, base + "?feature=exon;content-type=application/json",
                      cached_only=cached_only) or []
    txs = _rest_get(manifest, base + "?feature=transcript;content-type=application/json",
                    cached_only=cached_only) or []
    tx_by_id = {t.get("id"): t for t in txs if isinstance(t, dict)}
    out = []
    for x in exons:
        if not isinstance(x, dict):
            continue
        t = tx_by_id.get(x.get("Parent"), {})
        out.append({"number": x.get("rank"), "start": x.get("start"), "end": x.get("end"),
                    "transcript_id": x.get("Parent"),
                    "transcript_name": t.get("external_name"),
                    "gene_name": (t.get("external_name") or "").rsplit("-", 1)[0],
                    "strand": x.get("strand"), "biotype": t.get("biotype"),
                    "canonical": 1 if t.get("is_canonical") else 0})
    return out


def _rest_map(manifest, target_id, kind, start, end) -> dict | None:
    endpoint = "translation" if kind == "pep" else "cdna"
    data = _rest_get(manifest, f"/map/{endpoint}/{target_id.split('.')[0]}/"
                               f"{start}..{end}?content-type=application/json")
    maps = [m for m in (data or {}).get("mappings", []) if m.get("seq_region_name")]
    if not maps:
        return None
    return {"chrom": maps[0]["seq_region_name"], "start": min(m["start"] for m in maps),
            "end": max(m["end"] for m in maps), "strand": maps[0].get("strand", 1),
            "blocks": sorted((m["start"], m["end"]) for m in maps)}


def variants_in_region(manifest: dict, chrom: str, start: int, end: int,
                       cached_only: bool = False, somatic: bool = False) -> list[dict]:
    """Short variants from Ensembl Variation (online only; views are small regions).

    Germline variants by default; somatic=True gives the somatic mutations (COSMIC),
    which Ensembl reports without alleles. Ensembl can take a minute to answer for a
    variant-dense region it has not served recently, so views call this with
    cached_only=True (raising NotCached) and let the browser trigger the slow fetch
    in the background.
    """
    if end - start > 50_000:
        return []
    feature = "somatic_variation" if somatic else "variation"
    data = _rest_get(manifest, f"/overlap/region/{manifest['name']}/"
                               f"{_region(chrom, start, end)}?feature={feature};"
                               "content-type=application/json", max_age=90 * 86400,
                     cached_only=cached_only, timeout=120, attempts=3) or []
    out = []
    for v in data:
        if not isinstance(v, dict) or v.get("start") is None:
            continue
        out.append({"id": v.get("id"), "start": v["start"], "end": v["end"],
                    "alleles": "/".join(v.get("alleles") or []),
                    "allele_list": list(v.get("alleles") or []),
                    "consequence": (v.get("consequence_type") or "").replace("_", " "),
                    "consequence_type": v.get("consequence_type") or "",
                    "source": v.get("source")})
    return out
