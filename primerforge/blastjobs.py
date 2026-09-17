"""BLAST/BLAT tickets and jobs.

A ticket is one form submission. As on Ensembl it expands into one job per query
sequence per species; jobs run on a small worker pool, survive an app restart
(unfinished jobs are queued again), can be cancelled, and keep their hits in
SQLite and their raw output files under data/blast_jobs/<ticket>/<job>/.
"""
from __future__ import annotations

import json
import queue
import re
import shutil
import threading
import time
from pathlib import Path

from . import annotation, blastconf, config, engines, species, store
from .blastconf import (MAX_NUM_SEQUENCES, MAX_SEQUENCE_LENGTH, SEARCH_TYPE_BY_VALUE,
                        SOURCES, SPECIES_SELECTION_LIMIT, ConfigError)

SCHEMA = """
CREATE TABLE IF NOT EXISTS blast_tickets (
    id            TEXT PRIMARY KEY,
    created_at    REAL NOT NULL,
    description   TEXT,
    search_type   TEXT NOT NULL,
    query_type    TEXT NOT NULL,
    db_type       TEXT NOT NULL,
    source        TEXT NOT NULL,
    config_set    TEXT,
    configs_json  TEXT NOT NULL,
    species_json  TEXT NOT NULL,
    sequences_json TEXT NOT NULL,
    n_jobs        INTEGER NOT NULL DEFAULT 0
);
CREATE TABLE IF NOT EXISTS blast_jobs (
    id            TEXT PRIMARY KEY,
    ticket_id     TEXT NOT NULL REFERENCES blast_tickets(id) ON DELETE CASCADE,
    job_number    INTEGER NOT NULL,
    species       TEXT NOT NULL,
    assembly      TEXT,
    description   TEXT,
    summary       TEXT,
    seq_desc      TEXT,
    sequence      TEXT NOT NULL,
    seq_type      TEXT,
    status        TEXT NOT NULL,
    message       TEXT,
    notes_json    TEXT,
    created_at    REAL NOT NULL,
    started_at    REAL,
    finished_at   REAL,
    n_hits        INTEGER NOT NULL DEFAULT 0,
    n_found       INTEGER NOT NULL DEFAULT 0,
    command       TEXT,
    files_json    TEXT,
    engine        TEXT,
    role          TEXT
);
CREATE TABLE IF NOT EXISTS blast_hits (
    job_id    TEXT NOT NULL REFERENCES blast_jobs(id) ON DELETE CASCADE,
    idx       INTEGER NOT NULL,
    score     REAL, evalue REAL, pident REAL,
    gid       TEXT, gstart INTEGER, gend INTEGER,
    data_json TEXT NOT NULL,
    PRIMARY KEY (job_id, idx)
);
CREATE INDEX IF NOT EXISTS idx_bjobs_ticket ON blast_jobs(ticket_id);
CREATE INDEX IF NOT EXISTS idx_bjobs_status ON blast_jobs(status);
CREATE INDEX IF NOT EXISTS idx_btickets_created ON blast_tickets(created_at DESC);
"""

ACTIVE = ("queued", "running")

# A query can be a PCR primer. Primers sent from a design run carry their role
# explicitly; pasted ones are recognised by names such as "GALE_F" or "reverse primer".
PRIMER_ROLES = ("forward", "reverse")
PRIMER_MAX_LENGTH = 60
THREE_PRIME_WINDOW = 5          # bases at the 3' end that must pair for extension
_FORWARD_RE = re.compile(r"\b(forward|fwd|left)\b|(^|[\s_.\-])(f|fw)$", re.I)
_REVERSE_RE = re.compile(r"\b(reverse|rev|right)\b|(^|[\s_.\-])(r|rv)$", re.I)


class SubmitError(ValueError):
    pass


def primer_role(description: str | None, sequence: str, seq_type: str = "dna",
                explicit: str | None = None) -> str | None:
    if explicit in PRIMER_ROLES:
        return explicit
    if seq_type != "dna" or len(sequence) > PRIMER_MAX_LENGTH:
        return None
    text = (description or "").strip()
    fwd, rev = bool(_FORWARD_RE.search(text)), bool(_REVERSE_RE.search(text))
    return "forward" if fwd and not rev else "reverse" if rev and not fwd else None


_queue: "queue.Queue[str]" = queue.Queue()
_procs: dict[str, object] = {}
_cancelled: set[str] = set()
_lock = threading.Lock()
_started = False


def init() -> None:
    global _started
    with store.db() as conn:
        conn.executescript(SCHEMA)
        columns = {r["name"] for r in conn.execute("PRAGMA table_info(blast_jobs)")}
        if "role" not in columns:                  # databases created before primer roles
            conn.execute("ALTER TABLE blast_jobs ADD COLUMN role TEXT")
        columns = {r["name"] for r in conn.execute("PRAGMA table_info(blast_tickets)")}
        if "owner" not in columns:                 # who submitted it, on a shared server
            conn.execute("ALTER TABLE blast_tickets ADD COLUMN owner TEXT")
    with _lock:
        if _started:
            return
        _started = True
    with store.db() as conn:
        pending = [r["id"] for r in conn.execute(
            "SELECT id FROM blast_jobs WHERE status IN ('queued','running') "
            "ORDER BY created_at, job_number")]
        conn.execute("UPDATE blast_jobs SET status='queued', started_at=NULL "
                     "WHERE status='running'")
    for jid in pending:
        _queue.put(jid)
    for i in range(max(1, int(config.settings().get("blast_workers") or 2))):
        threading.Thread(target=_worker, name=f"blast-worker-{i}", daemon=True).start()


# ---------------------------------------------------------------------------
# Submission
# ---------------------------------------------------------------------------
def _species_summary(m: dict, search: dict, source: str) -> str:
    return (f"{search['method']} against {m.get('display_name', m['name'])} "
            f"{m.get('assembly', '')} ({species.source_label(source)})")


def validate(payload: dict) -> dict:
    """Check a form submission the way Ensembl's ticket does; return a clean ticket."""
    raw_seqs = payload.get("sequences")
    if isinstance(raw_seqs, str) or raw_seqs is None:
        parsed = blastconf.parse_sequences(raw_seqs or payload.get("sequence_text", ""))
        sequences = parsed["sequences"]
    else:
        sequences = []
        for s in raw_seqs:
            seq = "".join(str(s.get("sequence", "")).split()).upper()
            if not seq:
                continue
            if not blastconf.is_valid_sequence(seq) or not seq.replace("*", "").isalpha():
                raise SubmitError(f"Sequence '{s.get('description') or seq[:20]}' contains "
                                  "invalid characters.")
            if len(seq) > MAX_SEQUENCE_LENGTH:
                raise SubmitError(f"Sequence '{s.get('description') or seq[:20]}' is longer "
                                  f"than {MAX_SEQUENCE_LENGTH:,} characters.")
            reparsed = blastconf.parse_sequences(seq)["sequences"]
            seq_type = reparsed[0]["type"] if reparsed else "dna"
            if seq_type == "dna":
                seq = seq.replace("U", "T")
            sequences.append({"description": str(s.get("description") or "").strip()[:200],
                              "sequence": seq, "type": seq_type, "length": len(seq),
                              "role": s.get("role")})
    for s in sequences:
        s["role"] = primer_role(s["description"], s["sequence"], s["type"], s.get("role"))
    if not sequences:
        raise SubmitError("Enter at least one valid query sequence.")
    if len(sequences) > MAX_NUM_SEQUENCES:
        raise SubmitError(f"A maximum of {MAX_NUM_SEQUENCES} sequences can be searched at once.")

    names = payload.get("species") or []
    if isinstance(names, str):
        names = [names]
    names = list(dict.fromkeys(names))
    if not names:
        raise SubmitError("Select at least one species to search against.")
    if len(names) > SPECIES_SELECTION_LIMIT:
        raise SubmitError(f"Select at most {SPECIES_SELECTION_LIMIT} species.")
    manifests = []
    for n in names:
        m = species.load(n)
        if not m:
            raise SubmitError(f"Species '{n}' is not installed.")
        manifests.append(m)

    query_type = payload.get("query_type") or blastconf.guess_query_type(sequences)
    db_type = payload.get("db_type")
    source = payload.get("source")
    search_type = payload.get("search_type")
    if query_type not in blastconf.QUERY_TYPES:
        raise SubmitError("Choose DNA or protein as the query type.")
    if db_type not in blastconf.DB_TYPES:
        raise SubmitError("Choose a DNA or protein database to search against.")
    if source not in SOURCES or SOURCES[source]["db_type"] != db_type:
        raise SubmitError("Choose a database source that matches the database type.")
    search = SEARCH_TYPE_BY_VALUE.get(search_type)
    if not search:
        raise SubmitError("Choose a search tool.")
    if search["query_type"] != query_type or search["db_type"] != db_type \
            or source not in search["sources"]:
        raise SubmitError(f"{search['method']} cannot search a "
                          f"{blastconf.QUERY_TYPES[query_type]} query against "
                          f"{species.source_label(source)}.")
    if source in blastconf.RESTRICTIONS.get(search_type, []):
        raise SubmitError(f"{search['method']} cannot be used with "
                          f"{species.source_label(source)}.")

    mismatched = [s for s in sequences if s["type"] != query_type]
    if mismatched:
        kind = blastconf.QUERY_TYPES[mismatched[0]["type"]]
        raise SubmitError(
            f"{len(mismatched)} sequence(s) look like {kind} but the query type is "
            f"{blastconf.QUERY_TYPES[query_type]}. Change the query type or remove "
            f"'{mismatched[0]['description'] or mismatched[0]['sequence'][:20]}'.")
    if search["min_length"]:
        short = [s for s in sequences if s["length"] <= search["min_length"]]
        if short:
            raise SubmitError(f"{search['method']} needs sequences longer than "
                              f"{search['min_length']} bases; "
                              f"'{short[0]['description'] or short[0]['sequence']}' is "
                              f"{short[0]['length']}. Use BLASTN for short sequences.")

    for m in manifests:
        label = m.get("display_name", m["name"])
        if search_type == blastconf.BLAT_VALUE:
            if not species.blat_installed(m):
                raise SubmitError(f"BLAT is not installed for {label}. Install it in Settings "
                                  "or choose a BLAST program.")
        elif not species.available_sources(m).get(source):
            raise SubmitError(f"{species.source_label(source)} is not installed for {label}.")

    config_set = payload.get("config_set") or ""
    if config_set and config_set not in blastconf.CONFIG_SETS.get(search_type, {}):
        config_set = ""
    try:
        configs = blastconf.normalise_configs(search_type, payload.get("configs") or {},
                                              config_set=config_set)
    except ConfigError as exc:
        raise SubmitError(str(exc)) from exc

    return {"description": str(payload.get("description") or "").strip()[:300],
            "search_type": search_type, "query_type": query_type, "db_type": db_type,
            "source": source, "config_set": config_set, "configs": configs,
            "species": names, "sequences": sequences, "manifests": manifests}


def submit(payload: dict) -> str:
    t = validate(payload)
    ticket_id = store.new_id()
    now = time.time()
    search = SEARCH_TYPE_BY_VALUE[t["search_type"]]
    jobs = []
    number = 0
    for m in t["manifests"]:
        summary = _species_summary(m, search, t["source"])
        for s in t["sequences"]:
            number += 1
            jobs.append((store.new_id(), ticket_id, number, m["name"], m.get("assembly"),
                         t["description"] or s["description"] or summary, summary,
                         s["description"], s["sequence"], s["type"], "queued", now, s["role"]))
    with store.db() as conn:
        conn.execute(
            "INSERT INTO blast_tickets(id, created_at, description, search_type, query_type,"
            " db_type, source, config_set, configs_json, species_json, sequences_json, n_jobs)"
            " VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
            (ticket_id, now, t["description"], t["search_type"], t["query_type"], t["db_type"],
             t["source"], t["config_set"], json.dumps(t["configs"]), json.dumps(t["species"]),
             json.dumps([{k: s[k] for k in ("description", "sequence", "type", "role")}
                         for s in t["sequences"]]), len(jobs)))
        conn.executemany(
            "INSERT INTO blast_jobs(id, ticket_id, job_number, species, assembly, description,"
            " summary, seq_desc, sequence, seq_type, status, created_at, role)"
            " VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)", jobs)
    for j in jobs:
        _queue.put(j[0])
    return ticket_id


# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------
def job_dir(ticket_id: str, job_id: str) -> Path:
    return config.BLAST_JOBS_DIR / ticket_id / job_id


def _worker() -> None:
    while True:
        jid = _queue.get()
        try:
            _run(jid)
        except Exception as exc:                                  # noqa: BLE001
            _finish(jid, "failed", message=f"Unexpected error: {exc}")
        finally:
            _queue.task_done()


def _set(jid: str, **fields) -> None:
    cols = ", ".join(f"{k}=?" for k in fields)
    with store.db() as conn:
        conn.execute(f"UPDATE blast_jobs SET {cols} WHERE id=?", (*fields.values(), jid))


def _finish(jid: str, status: str, message: str | None = None, **extra) -> None:
    _set(jid, status=status, message=message, finished_at=time.time(), **extra)


def _run(jid: str) -> None:
    with store.db() as conn:
        row = conn.execute("SELECT * FROM blast_jobs WHERE id=?", (jid,)).fetchone()
        if not row or row["status"] != "queued":
            return
        trow = conn.execute("SELECT * FROM blast_tickets WHERE id=?",
                            (row["ticket_id"],)).fetchone()
    job, ticket = dict(row), _ticket_dict(trow)
    if jid in _cancelled:
        _finish(jid, "cancelled", "Cancelled before it started.")
        return
    manifest = species.load(job["species"])
    if not manifest:
        _finish(jid, "failed", f"Species {job['species']} is no longer installed.")
        return

    workdir = job_dir(job["ticket_id"], jid)
    workdir.mkdir(parents=True, exist_ok=True)
    safe_id = f"query_{job['job_number']}"
    (workdir / "input.fa").write_text(blastconf.fasta(safe_id, job["sequence"]),
                                      encoding="utf-8", newline="\n")
    _set(jid, status="running", started_at=time.time(), message=None)

    def on_proc(proc):
        with _lock:
            if proc is None:
                _procs.pop(jid, None)
            else:
                _procs[jid] = proc

    runner = engines.run_blat if ticket["search_type"] == blastconf.BLAT_VALUE \
        else engines.run_ncbi
    try:
        result = runner(job, ticket, manifest, workdir, on_proc=on_proc,
                        is_cancelled=lambda: jid in _cancelled)
    except engines.Cancelled:
        _finish(jid, "cancelled", "Cancelled.")
        return
    except engines.EngineError as exc:
        _finish(jid, "cancelled" if jid in _cancelled else "failed", str(exc))
        return

    if not get_job(jid):            # deleted while the search was running
        return
    hits = result["hits"]
    found = len(hits)
    notes = list(result["notes"])
    if found > blastconf.MAX_STORED_HITS:
        notes.append(f"{found:,} alignments were found; the best "
                     f"{blastconf.MAX_STORED_HITS:,} are shown.")
        hits = hits[:blastconf.MAX_STORED_HITS]
    if ticket["source"] in blastconf.GENOMIC_SOURCES and annotation.has_local(manifest):
        for h in hits:
            h["genes"] = [{"id": g["id"], "name": g["name"], "biotype": g["biotype"],
                           "strand": g["strand"]}
                          for g in annotation.genes_in_region(manifest, h["gid"], h["gstart"],
                                                              h["gend"])]
    rows = [(jid, i, h["score"], h["evalue"], h["pident"], h.get("gid"), h.get("gstart"),
             h.get("gend"), json.dumps(h)) for i, h in enumerate(hits)]
    with store.db() as conn:
        conn.execute("DELETE FROM blast_hits WHERE job_id=?", (jid,))
        conn.executemany("INSERT INTO blast_hits(job_id, idx, score, evalue, pident, gid,"
                         " gstart, gend, data_json) VALUES (?,?,?,?,?,?,?,?,?)", rows)
    _finish(jid, "done", None, n_hits=len(hits), n_found=found, notes_json=json.dumps(notes),
            command=result["command"], files_json=json.dumps(result["files"]),
            engine=result["engine"])


def cancel(jid: str) -> bool:
    job = get_job(jid)
    if not job or job["status"] not in ACTIVE:
        return False
    _cancelled.add(jid)
    with _lock:
        proc = _procs.get(jid)
    if proc is not None:
        try:
            proc.kill()
        except OSError:
            pass
    if job["status"] == "queued":
        _finish(jid, "cancelled", "Cancelled before it started.")
    return True


def resubmit_job(jid: str) -> bool:
    """Run a failed or cancelled job again with its original settings."""
    job = get_job(jid)
    if not job or job["status"] in ACTIVE:
        return False
    _cancelled.discard(jid)
    _set(jid, status="queued", message=None, started_at=None, finished_at=None, n_hits=0,
         n_found=0, notes_json=None)
    with store.db() as conn:
        conn.execute("DELETE FROM blast_hits WHERE job_id=?", (jid,))
    _queue.put(jid)
    return True


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------
def _ticket_dict(row) -> dict:
    t = dict(row)
    t["configs"] = json.loads(t.pop("configs_json"))
    t["species"] = json.loads(t.pop("species_json"))
    t["sequences"] = json.loads(t.pop("sequences_json"))
    t["search_caption"] = blastconf.search_caption(t["search_type"])
    t["method"] = SEARCH_TYPE_BY_VALUE.get(t["search_type"], {}).get("method", t["search_type"])
    t["source_label"] = species.source_label(t["source"])
    return t


def _job_dict(row) -> dict:
    j = dict(row)
    j["notes"] = json.loads(j.pop("notes_json") or "[]")
    j["files"] = json.loads(j.pop("files_json") or "{}")
    if not j.get("role") and j.get("sequence"):
        j["role"] = primer_role(j.get("seq_desc"), j["sequence"], j.get("seq_type") or "dna")
    return j


def get_ticket(ticket_id: str, with_jobs: bool = True) -> dict | None:
    with store.db() as conn:
        row = conn.execute("SELECT * FROM blast_tickets WHERE id=?", (ticket_id,)).fetchone()
        if not row:
            return None
        t = _ticket_dict(row)
        if with_jobs:
            t["jobs"] = [_job_dict(r) for r in conn.execute(
                "SELECT * FROM blast_jobs WHERE ticket_id=? ORDER BY job_number", (ticket_id,))]
    return t


def set_ticket_owner(ticket_id: str, owner: str | None) -> None:
    if owner:
        with store.db() as conn:
            conn.execute("UPDATE blast_tickets SET owner=? WHERE id=?", (owner, ticket_id))


def get_job(jid: str) -> dict | None:
    with store.db() as conn:
        row = conn.execute("SELECT * FROM blast_jobs WHERE id=?", (jid,)).fetchone()
    return _job_dict(row) if row else None


def list_tickets(limit: int = 100) -> list[dict]:
    with store.db() as conn:
        tickets = [_ticket_dict(r) for r in conn.execute(
            "SELECT * FROM blast_tickets ORDER BY created_at DESC LIMIT ?", (limit,))]
        ids = [t["id"] for t in tickets]
        jobs: dict[str, list] = {i: [] for i in ids}
        if ids:
            marks = ",".join("?" * len(ids))
            for r in conn.execute(
                    "SELECT id, ticket_id, job_number, species, assembly, description, summary,"
                    " seq_desc, length(sequence) AS seq_len, seq_type, status, message,"
                    " created_at, started_at, finished_at, n_hits, n_found, role"
                    f" FROM blast_jobs WHERE ticket_id IN ({marks}) ORDER BY job_number", ids):
                jobs[r["ticket_id"]].append(dict(r))
    for t in tickets:
        t.pop("sequences", None)
        t["jobs"] = jobs.get(t["id"], [])
    return tickets


def get_hits(jid: str, full: bool = False) -> list[dict]:
    with store.db() as conn:
        rows = conn.execute("SELECT idx, data_json FROM blast_hits WHERE job_id=? ORDER BY idx",
                            (jid,)).fetchall()
    out = []
    for r in rows:
        h = json.loads(r["data_json"])
        h["idx"] = r["idx"]
        _primer_flags(h)
        if not full:
            for k in ("qseq", "sseq", "aln"):
                h.pop(k, None)
        out.append(h)
    return out


def _primer_flags(h: dict) -> None:
    """Whether a hit could prime: a polymerase only extends from a paired 3' end.

    three_prime_ok -- the alignment reaches the query's last base and the final
                      five bases pair with no gap (the rule the in-silico PCR uses)
    full_match     -- additionally the whole query matches without mismatches
    """
    q, s = (h.get("qseq") or "").upper(), (h.get("sseq") or "").upper()
    nucleotide = h.get("program") in ("blastn", "blat")
    tail_q, tail_s = q[-THREE_PRIME_WINDOW:], s[-THREE_PRIME_WINDOW:]
    ok = (nucleotide and (h.get("qori") or 1) > 0 and h.get("qend") == h.get("qlen")
          and len(tail_q) == THREE_PRIME_WINDOW and "-" not in tail_q + tail_s
          and tail_q == tail_s)
    h["three_prime_ok"] = bool(ok)
    h["full_match"] = bool(ok and h.get("qstart") == 1 and h.get("pident", 0) >= 100
                           and not h.get("gapopen"))


def get_hit(jid: str, idx: int) -> dict | None:
    with store.db() as conn:
        r = conn.execute("SELECT idx, data_json FROM blast_hits WHERE job_id=? AND idx=?",
                         (jid, idx)).fetchone()
    if not r:
        return None
    h = json.loads(r["data_json"])
    h["idx"] = r["idx"]
    return h


def update_hit_genes(jid: str, genes_by_idx: dict[int, list]) -> None:
    with store.db() as conn:
        for idx, genes in genes_by_idx.items():
            r = conn.execute("SELECT data_json FROM blast_hits WHERE job_id=? AND idx=?",
                             (jid, idx)).fetchone()
            if r:
                h = json.loads(r["data_json"])
                h["genes"] = genes
                conn.execute("UPDATE blast_hits SET data_json=? WHERE job_id=? AND idx=?",
                             (json.dumps(h), jid, idx))


def delete_ticket(ticket_id: str) -> None:
    t = get_ticket(ticket_id)
    if not t:
        return
    for j in t["jobs"]:
        if j["status"] in ACTIVE:
            cancel(j["id"])
    with store.db() as conn:
        conn.execute("DELETE FROM blast_hits WHERE job_id IN "
                     "(SELECT id FROM blast_jobs WHERE ticket_id=?)", (ticket_id,))
        conn.execute("DELETE FROM blast_jobs WHERE ticket_id=?", (ticket_id,))
        conn.execute("DELETE FROM blast_tickets WHERE id=?", (ticket_id,))
    shutil.rmtree(config.BLAST_JOBS_DIR / ticket_id, ignore_errors=True)


def delete_job(jid: str) -> str | None:
    """Delete one job; the ticket goes too once it has no jobs left."""
    job = get_job(jid)
    if not job:
        return None
    if job["status"] in ACTIVE:
        cancel(jid)
    with store.db() as conn:
        conn.execute("DELETE FROM blast_hits WHERE job_id=?", (jid,))
        conn.execute("DELETE FROM blast_jobs WHERE id=?", (jid,))
        left = conn.execute("SELECT COUNT(*) AS n FROM blast_jobs WHERE ticket_id=?",
                            (job["ticket_id"],)).fetchone()["n"]
        if left:
            conn.execute("UPDATE blast_tickets SET n_jobs=? WHERE id=?", (left, job["ticket_id"]))
    shutil.rmtree(job_dir(job["ticket_id"], jid), ignore_errors=True)
    if not left:
        delete_ticket(job["ticket_id"])
    return job["ticket_id"]


def queue_position(jid: str) -> int | None:
    with store.db() as conn:
        rows = [r["id"] for r in conn.execute(
            "SELECT id FROM blast_jobs WHERE status='queued' ORDER BY created_at, job_number")]
    return rows.index(jid) + 1 if jid in rows else None
