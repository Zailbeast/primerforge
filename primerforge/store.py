"""SQLite-backed persistence: design runs, their results, and an API cache.

Every design is recorded so the dashboard can replay and compare past work.
Results are stored denormalised enough that the run pages never need to hit
the network again.
"""
from __future__ import annotations

import json
import sqlite3
import time
import uuid
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator

from .config import DB_PATH, ensure_dirs

SCHEMA = """
CREATE TABLE IF NOT EXISTS runs (
    id            TEXT PRIMARY KEY,
    created_at    REAL NOT NULL,
    finished_at   REAL,
    label         TEXT,
    assay         TEXT NOT NULL,
    input_raw     TEXT NOT NULL,
    params_json   TEXT NOT NULL,
    status        TEXT NOT NULL,
    n_input       INTEGER NOT NULL DEFAULT 0,
    n_ok          INTEGER NOT NULL DEFAULT 0,
    n_failed      INTEGER NOT NULL DEFAULT 0,
    specificity   INTEGER NOT NULL DEFAULT 0,
    error         TEXT
);

CREATE TABLE IF NOT EXISTS variants (
    id            TEXT PRIMARY KEY,
    run_id        TEXT NOT NULL REFERENCES runs(id) ON DELETE CASCADE,
    idx           INTEGER NOT NULL,
    input_hgvs    TEXT NOT NULL,
    status        TEXT NOT NULL,
    error         TEXT,
    gene          TEXT,
    chrom         TEXT,
    pos           INTEGER,
    locus_json    TEXT,
    template_json TEXT,
    n_pairs       INTEGER NOT NULL DEFAULT 0,
    warnings_json TEXT
);

CREATE TABLE IF NOT EXISTS pairs (
    id            TEXT PRIMARY KEY,
    variant_id    TEXT NOT NULL REFERENCES variants(id) ON DELETE CASCADE,
    run_id        TEXT NOT NULL,
    rank          INTEGER NOT NULL,
    left_seq      TEXT NOT NULL,
    right_seq     TEXT NOT NULL,
    left_tm       REAL, right_tm REAL,
    left_gc       REAL, right_gc REAL,
    product_size  INTEGER,
    penalty       REAL,
    metrics_json  TEXT,
    spec_json     TEXT,
    spec_status   TEXT,
    n_amplicons   INTEGER
);

CREATE TABLE IF NOT EXISTS cache (
    key        TEXT PRIMARY KEY,
    value      TEXT NOT NULL,
    created_at REAL NOT NULL
);

CREATE INDEX IF NOT EXISTS idx_variants_run ON variants(run_id);
CREATE INDEX IF NOT EXISTS idx_pairs_variant ON pairs(variant_id);
CREATE INDEX IF NOT EXISTS idx_runs_created ON runs(created_at DESC);
"""


def _connect() -> sqlite3.Connection:
    ensure_dirs()
    conn = sqlite3.connect(DB_PATH, timeout=30)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA foreign_keys=ON")
    return conn


@contextmanager
def db() -> Iterator[sqlite3.Connection]:
    conn = _connect()
    try:
        yield conn
        conn.commit()
    finally:
        conn.close()


def init() -> None:
    with db() as conn:
        conn.executescript(SCHEMA)
        columns = {r["name"] for r in conn.execute("PRAGMA table_info(runs)")}
        if "owner" not in columns:                 # who designed it, on a shared server
            conn.execute("ALTER TABLE runs ADD COLUMN owner TEXT")


def set_run_owner(run_id: str, owner: str | None) -> None:
    if owner:
        with db() as conn:
            conn.execute("UPDATE runs SET owner=? WHERE id=?", (owner, run_id))


def new_id() -> str:
    return uuid.uuid4().hex[:12]


# --------------------------------------------------------------------------
# Cache (Ensembl responses are slow and highly repeatable)
# --------------------------------------------------------------------------
def cache_get(key: str, max_age: float = 30 * 86400) -> Any | None:
    with db() as conn:
        row = conn.execute("SELECT value, created_at FROM cache WHERE key=?", (key,)).fetchone()
    if row and (time.time() - row["created_at"]) < max_age:
        return json.loads(row["value"])
    return None


def cache_put(key: str, value: Any) -> None:
    with db() as conn:
        conn.execute(
            "INSERT OR REPLACE INTO cache(key, value, created_at) VALUES (?,?,?)",
            (key, json.dumps(value), time.time()),
        )


# Nothing reads a cached Ensembl response older than 90 days (the longest max_age
# any caller asks for), so those rows are dead weight; on a shared server they
# would otherwise grow without limit.
CACHE_MAX_AGE = 90 * 86400


def prune_cache(max_age: float = CACHE_MAX_AGE) -> int:
    """Delete cache rows too old for any caller to use. Returns the row count."""
    with db() as conn:
        cur = conn.execute("DELETE FROM cache WHERE created_at < ?", (time.time() - max_age,))
        return cur.rowcount


def backup(target: Path) -> Path:
    """Copy the database to `target` using SQLite's online backup.

    Safe while the server is running and while WAL writes are in flight, which
    plain file copying is not: the result is always a consistent database.
    """
    target.parent.mkdir(parents=True, exist_ok=True)
    with db() as conn:
        dest = sqlite3.connect(target)
        try:
            conn.backup(dest)
        finally:
            dest.close()
    return target


# --------------------------------------------------------------------------
# Runs
# --------------------------------------------------------------------------
def create_run(label: str, assay: str, input_raw: str, params: dict, n_input: int,
               specificity: bool) -> str:
    run_id = new_id()
    with db() as conn:
        conn.execute(
            "INSERT INTO runs(id, created_at, label, assay, input_raw, params_json, status,"
            " n_input, specificity) VALUES (?,?,?,?,?,?,?,?,?)",
            (run_id, time.time(), label, assay, input_raw, json.dumps(params), "running",
             n_input, int(specificity)),
        )
    return run_id


def abandon_running_runs(message: str) -> int:
    """Close runs whose process died, and report how many there were.

    A design run lives on a worker thread, so a restart leaves its row at
    'running' with nothing left to finish it. Unlike a BLAST job it cannot be
    requeued: the variants designed before the restart are already stored, and
    repeating them would duplicate rows. The run is closed with what it managed
    to produce instead, so the page says what happened rather than hanging.
    """
    with db() as conn:
        rows = conn.execute("SELECT id FROM runs WHERE status='running'").fetchall()
        for row in rows:
            counts = conn.execute(
                "SELECT COALESCE(SUM(status='ok'), 0) AS ok,"
                " COALESCE(SUM(status<>'ok'), 0) AS failed"
                " FROM variants WHERE run_id=?", (row["id"],)).fetchone()
            conn.execute(
                "UPDATE runs SET status=?, finished_at=?, n_ok=?, n_failed=?, error=?"
                " WHERE id=?",
                ("partial" if counts["ok"] else "failed", time.time(),
                 counts["ok"], counts["failed"], message, row["id"]))
    return len(rows)


def finish_run(run_id: str, status: str, n_ok: int, n_failed: int,
               error: str | None = None) -> None:
    with db() as conn:
        conn.execute(
            "UPDATE runs SET status=?, finished_at=?, n_ok=?, n_failed=?, error=? WHERE id=?",
            (status, time.time(), n_ok, n_failed, error, run_id),
        )


def add_variant(run_id: str, idx: int, input_hgvs: str, status: str, error: str | None = None,
                locus: dict | None = None, template: dict | None = None,
                warnings: list | None = None) -> str:
    vid = new_id()
    with db() as conn:
        conn.execute(
            "INSERT INTO variants(id, run_id, idx, input_hgvs, status, error, gene, chrom, pos,"
            " locus_json, template_json, warnings_json) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
            (vid, run_id, idx, input_hgvs, status, error,
             (locus or {}).get("gene"), (locus or {}).get("chrom"), (locus or {}).get("start"),
             json.dumps(locus) if locus else None,
             json.dumps(template) if template else None,
             json.dumps(warnings or [])),
        )
    return vid


def add_pairs(variant_id: str, run_id: str, pairs: list[dict]) -> None:
    rows = []
    for p in pairs:
        spec = p.get("specificity") or {}
        rows.append((
            new_id(), variant_id, run_id, p["rank"], p["left"]["seq"], p["right"]["seq"],
            p["left"]["tm"], p["right"]["tm"], p["left"]["gc"], p["right"]["gc"],
            p["product_size"], p["penalty"], json.dumps(p),
            json.dumps(spec) if spec else None, spec.get("status"), spec.get("n_amplicons"),
        ))
    with db() as conn:
        conn.executemany(
            "INSERT INTO pairs(id, variant_id, run_id, rank, left_seq, right_seq, left_tm,"
            " right_tm, left_gc, right_gc, product_size, penalty, metrics_json, spec_json,"
            " spec_status, n_amplicons) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", rows)
        conn.execute("UPDATE variants SET n_pairs=? WHERE id=?", (len(pairs), variant_id))


def get_run(run_id: str) -> dict | None:
    with db() as conn:
        row = conn.execute("SELECT * FROM runs WHERE id=?", (run_id,)).fetchone()
        if not row:
            return None
        run = dict(row)
        run["params"] = json.loads(run["params_json"])
        variants = []
        for vrow in conn.execute(
                "SELECT * FROM variants WHERE run_id=? ORDER BY idx", (run_id,)):
            v = dict(vrow)
            v["locus"] = json.loads(v["locus_json"]) if v["locus_json"] else None
            v["template"] = json.loads(v["template_json"]) if v["template_json"] else None
            v["warnings"] = json.loads(v["warnings_json"] or "[]")
            v["pairs"] = [
                json.loads(p["metrics_json"])
                for p in conn.execute(
                    "SELECT metrics_json FROM pairs WHERE variant_id=? ORDER BY rank",
                    (v["id"],))
            ]
            variants.append(v)
        run["variants"] = variants
    return run


def list_runs(limit: int = 500, search: str = "", assay: str = "") -> list[dict]:
    sql = ("SELECT r.*, (SELECT group_concat(DISTINCT gene) FROM variants v"
           " WHERE v.run_id=r.id) AS genes FROM runs r WHERE 1=1")
    args: list = []
    if search:
        sql += (" AND (r.label LIKE ? OR r.input_raw LIKE ? OR r.id LIKE ? OR EXISTS"
                " (SELECT 1 FROM variants v WHERE v.run_id=r.id AND (v.gene LIKE ?"
                " OR v.input_hgvs LIKE ?)))")
        args += [f"%{search}%"] * 5
    if assay:
        sql += " AND r.assay=?"
        args.append(assay)
    sql += " ORDER BY r.created_at DESC LIMIT ?"
    args.append(limit)
    with db() as conn:
        return [dict(r) for r in conn.execute(sql, args)]


def delete_run(run_id: str) -> None:
    with db() as conn:
        conn.execute("DELETE FROM pairs WHERE run_id=?", (run_id,))
        conn.execute("DELETE FROM variants WHERE run_id=?", (run_id,))
        conn.execute("DELETE FROM runs WHERE id=?", (run_id,))


def stats() -> dict:
    """Aggregates that drive the dashboard charts."""
    with db() as conn:
        out: dict[str, Any] = {}
        out["totals"] = dict(conn.execute(
            "SELECT COUNT(*) AS runs, COALESCE(SUM(n_input),0) AS variants,"
            " COALESCE(SUM(n_ok),0) AS ok, COALESCE(SUM(n_failed),0) AS failed FROM runs"
        ).fetchone())
        out["pairs_total"] = conn.execute("SELECT COUNT(*) AS n FROM pairs").fetchone()["n"]
        out["by_day"] = [dict(r) for r in conn.execute(
            "SELECT date(created_at, 'unixepoch', 'localtime') AS day, COUNT(*) AS n"
            " FROM runs GROUP BY day ORDER BY day DESC LIMIT 60")]
        out["by_assay"] = [dict(r) for r in conn.execute(
            "SELECT assay, COUNT(*) AS n FROM runs GROUP BY assay ORDER BY n DESC")]
        out["by_gene"] = [dict(r) for r in conn.execute(
            "SELECT gene, COUNT(*) AS n FROM variants WHERE gene IS NOT NULL"
            " GROUP BY gene ORDER BY n DESC LIMIT 15")]
        out["pairs"] = [dict(r) for r in conn.execute(
            "SELECT left_tm, right_tm, product_size, penalty, spec_status FROM pairs"
            " WHERE left_tm IS NOT NULL LIMIT 5000")]
        out["spec"] = [dict(r) for r in conn.execute(
            "SELECT COALESCE(spec_status, 'not checked') AS status, COUNT(*) AS n"
            " FROM pairs GROUP BY status")]
    return out
