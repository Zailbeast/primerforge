"""Specificity checking by local BLAST plus in-silico PCR.

A primer only primes where its 3' end anneals, so raw BLAST hits are filtered
on 3'-terminal complementarity before being paired up into predicted amplicons.
Any two binding sites pointing at each other within the product size range make
a product -- including forward/forward and reverse/reverse combinations, which
are a common cause of unexpected bands.
"""
from __future__ import annotations

import os
import subprocess
import tempfile
from dataclasses import dataclass

from .config import BLAST_DB, blast_bin, blast_ready

OUTFMT = ("6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send"
          " evalue bitscore sstrand qlen qseq sseq")

# A hit is only amplifiable if the 3'-terminal bases pair perfectly.
THREE_PRIME_WINDOW = 5
MAX_TOTAL_MISMATCH = 5


class SpecificityError(RuntimeError):
    pass


@dataclass
class Hit:
    query: str
    chrom: str
    direction: int          # +1 extends towards higher coordinates
    outer: int              # genomic position of the primer 5' end
    three_prime: int        # genomic position of the primer 3' end
    mismatches: int
    identity: float
    three_prime_ok: bool


def _parse_hit(cols: list[str]) -> Hit | None:
    (qseqid, sseqid, pident, _length, mismatch, gapopen, qstart, qend,
     sstart, send, _evalue, _bits, sstrand, qlen, qseq, sseq) = cols[:16]
    qstart, qend = int(qstart), int(qend)
    sstart, send, qlen = int(sstart), int(send), int(qlen)

    # Require the alignment to reach the primer's 3'-most base.
    if qend != qlen:
        return None
    if int(gapopen) > 0:
        return None

    if sstrand == "plus":
        direction, outer, three_prime = 1, sstart - (qstart - 1), send
    else:
        direction, outer, three_prime = -1, sstart + (qstart - 1), send

    # The final few aligned bases must match exactly for extension to occur.
    tail_q, tail_s = qseq[-THREE_PRIME_WINDOW:], sseq[-THREE_PRIME_WINDOW:]
    three_prime_ok = tail_q.upper() == tail_s.upper()

    return Hit(query=qseqid, chrom=sseqid, direction=direction, outer=outer,
               three_prime=three_prime, mismatches=int(mismatch),
               identity=float(pident), three_prime_ok=three_prime_ok)


def _run_blast(primers: dict[str, str], threads: int = 0) -> list[Hit]:
    exe = blast_bin("blastn")
    if not exe or not blast_ready():
        raise SpecificityError("Local BLAST database is not installed.")
    threads = threads or max(1, (os.cpu_count() or 4) - 1)

    fasta = "".join(f">{name}\n{seq}\n" for name, seq in primers.items())
    tmp = tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False, encoding="utf-8")
    try:
        tmp.write(fasta)
        tmp.close()
        cmd = [exe, "-task", "blastn-short", "-db", str(BLAST_DB), "-query", tmp.name,
               "-outfmt", OUTFMT, "-evalue", "1000", "-word_size", "7",
               "-dust", "no", "-soft_masking", "false", "-max_target_seqs", "20000",
               "-num_threads", str(threads), "-penalty", "-2", "-reward", "1"]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=900,
                              creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
        if proc.returncode != 0:
            raise SpecificityError(f"blastn failed: {proc.stderr.strip()[:300]}")
        hits: list[Hit] = []
        for line in proc.stdout.splitlines():
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 16:
                hit = _parse_hit(cols)
                if hit and hit.mismatches <= MAX_TOTAL_MISMATCH:
                    hits.append(hit)
        return hits
    finally:
        try:
            os.unlink(tmp.name)
        except OSError:
            pass


def _amplicons(hits: list[Hit], min_size: int, max_size: int,
               slack: float = 1.5) -> list[dict]:
    """Pair opposing binding sites into predicted products."""
    by_chrom: dict[str, list[Hit]] = {}
    for h in hits:
        if h.three_prime_ok:
            by_chrom.setdefault(h.chrom, []).append(h)

    limit = int(max_size * slack)
    out: list[dict] = []
    for chrom, group in by_chrom.items():
        fwd = sorted((h for h in group if h.direction == 1), key=lambda h: h.outer)
        rev = sorted((h for h in group if h.direction == -1), key=lambda h: h.outer)
        if not fwd or not rev:
            continue
        for f in fwd:
            for r in rev:
                if r.outer <= f.outer:
                    continue
                size = r.outer - f.outer + 1
                if size > limit:
                    break
                if size < min_size:
                    continue
                out.append({
                    "chrom": chrom, "start": f.outer, "end": r.outer, "size": size,
                    "fwd": f.query, "rev": r.query,
                    "mismatches": f.mismatches + r.mismatches,
                    "identity": round((f.identity + r.identity) / 2, 1),
                })
    out.sort(key=lambda a: (a["mismatches"], -a["size"]))
    return out


def check_pairs(pairs: list[dict], locus_chrom: str, min_size: int, max_size: int,
                threads: int = 0) -> None:
    """Annotate each pair in place with its in-silico PCR result."""
    if not pairs:
        return
    if not blast_ready():
        for p in pairs:
            p["specificity"] = {
                "status": "unavailable",
                "message": "Local genome/BLAST database not installed.",
            }
        return

    primers: dict[str, str] = {}
    for p in pairs:
        primers[f"p{p['rank']}_F"] = p["left"]["seq"]
        primers[f"p{p['rank']}_R"] = p["right"]["seq"]

    try:
        hits = _run_blast(primers, threads=threads)
    except SpecificityError as exc:
        for p in pairs:
            p["specificity"] = {"status": "error", "message": str(exc)}
        return

    for p in pairs:
        names = {f"p{p['rank']}_F", f"p{p['rank']}_R"}
        mine = [h for h in hits if h.query in names]
        amps = _amplicons(mine, min_size, max_size)

        intended_start = p["left"]["start"]
        on_target, off_target = [], []
        for a in amps:
            same = (a["chrom"].replace("chr", "") == locus_chrom.replace("chr", "")
                    and abs(a["start"] - intended_start) <= 10
                    and a["fwd"].endswith("_F") and a["rev"].endswith("_R"))
            (on_target if same else off_target).append(a)

        n_sites = {
            "left": sum(1 for h in mine if h.query.endswith("_F") and h.three_prime_ok),
            "right": sum(1 for h in mine if h.query.endswith("_R") and h.three_prime_ok),
        }

        if not on_target:
            status, message = "no product", (
                "The intended amplicon was not recovered from the genome. "
                "Check the assembly build.")
        elif not off_target:
            status, message = "unique", "Single predicted product genome-wide."
        elif len(off_target) <= 2:
            status, message = "minor off-target", (
                f"{len(off_target)} additional predicted product(s); usually "
                "tolerable but worth checking on a gel.")
        else:
            status, message = "non-specific", (
                f"{len(off_target)} additional predicted products - expect "
                "extra bands.")

        p["specificity"] = {
            "status": status, "message": message,
            "n_amplicons": len(on_target) + len(off_target),
            "on_target": on_target[:1],
            "off_target": off_target[:20],
            "n_off_target": len(off_target),
            "binding_sites": n_sites,
        }
        if status in ("non-specific", "no product"):
            p["warnings"].append(message)
