"""Reference sequence retrieval.

Prefers the local soft-masked GRCh38 FASTA (instant, offline, and lower-case
letters mark repeats) and falls back to the Ensembl REST API so the tool is
fully usable before the genome has been downloaded.
"""
from __future__ import annotations

import threading
from pathlib import Path

from . import net, store
from .config import (ENSEMBL_REST, ENSEMBL_REST_GRCH37, GENOME_FAI, GENOME_FASTA,
                     genome_ready)


class SequenceError(RuntimeError):
    pass


# --------------------------------------------------------------------------
# Minimal samtools-compatible FASTA index
# --------------------------------------------------------------------------
class FastaIndex:
    """Random access into a plain FASTA using a .fai sidecar index."""

    def __init__(self, fasta: Path, fai: Path):
        self.fasta = fasta
        self.index: dict[str, tuple[int, int, int, int]] = {}
        self._lock = threading.Lock()
        self._fh = None
        with open(fai, "r", encoding="utf-8") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 5:
                    name, length, offset, linebases, linewidth = parts[:5]
                    self.index[name] = (int(length), int(offset), int(linebases),
                                        int(linewidth))

    @staticmethod
    def build(fasta: Path, fai: Path, progress=None) -> None:
        """Write a .fai for a FASTA with uniform line lengths."""
        entries: list[str] = []
        name = None
        length = offset = linebases = linewidth = 0
        with open(fasta, "rb") as fh:
            pos = 0
            for raw in fh:
                if raw.startswith(b">"):
                    if name is not None:
                        entries.append(
                            f"{name}\t{length}\t{offset}\t{linebases}\t{linewidth}")
                    name = raw[1:].split()[0].decode()
                    length = linebases = linewidth = 0
                    offset = pos + len(raw)
                    if progress:
                        progress(name)
                else:
                    stripped = len(raw.rstrip(b"\r\n"))
                    if linebases == 0:
                        linebases, linewidth = stripped, len(raw)
                    length += stripped
                pos += len(raw)
        if name is not None:
            entries.append(f"{name}\t{length}\t{offset}\t{linebases}\t{linewidth}")
        fai.write_text("\n".join(entries) + "\n", encoding="utf-8")

    def contigs(self) -> dict[str, int]:
        return {k: v[0] for k, v in self.index.items()}

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """1-based inclusive slice; clamps to contig bounds."""
        key = self._resolve(chrom)
        length, offset, linebases, linewidth = self.index[key]
        start = max(1, start)
        end = min(length, end)
        if end < start:
            return ""
        begin = offset + (start - 1) // linebases * linewidth + (start - 1) % linebases
        stop = offset + (end - 1) // linebases * linewidth + (end - 1) % linebases
        with self._lock:
            if self._fh is None:
                self._fh = open(self.fasta, "rb")
            self._fh.seek(begin)
            raw = self._fh.read(stop - begin + 1)
        return raw.replace(b"\n", b"").replace(b"\r", b"").decode()

    def _resolve(self, chrom: str) -> str:
        for cand in (chrom, chrom.replace("chr", ""), "chr" + chrom,
                     chrom.upper(), chrom.replace("M", "MT")):
            if cand in self.index:
                return cand
        raise SequenceError(f"Contig '{chrom}' is not in the local genome index.")


_index: FastaIndex | None = None
_index_lock = threading.Lock()


def local_index() -> FastaIndex | None:
    global _index
    if not genome_ready():
        return None
    with _index_lock:
        if _index is None:
            _index = FastaIndex(GENOME_FASTA, GENOME_FAI)
    return _index


def reset_index() -> None:
    """Drop the cached handle after the genome is (re)installed."""
    global _index
    with _index_lock:
        _index = None


# --------------------------------------------------------------------------
# Public API
# --------------------------------------------------------------------------
def source(assembly: str = "GRCh38") -> str:
    if assembly == "GRCh38" and local_index() is not None:
        return "local"
    return "ensembl"


def fetch(chrom: str, start: int, end: int, assembly: str = "GRCh38") -> str:
    """Reference sequence for a 1-based inclusive region, soft-masked."""
    if start > end:
        raise SequenceError(f"Invalid region {chrom}:{start}-{end}.")
    if end - start + 1 > 5_000_000:
        raise SequenceError("Requested region is too large (limit 5 Mb).")

    idx = local_index() if assembly == "GRCh38" else None
    if idx is not None:
        seq = idx.fetch(chrom, start, end)
        if seq:
            return seq
        raise SequenceError(f"No sequence at {chrom}:{start}-{end} in the local genome.")

    base = ENSEMBL_REST if assembly == "GRCh38" else ENSEMBL_REST_GRCH37
    region = f"{chrom.replace('chr', '')}:{start}..{end}:1"
    url = f"{base}/sequence/region/human/{region}?mask=soft&content-type=application/json"
    cached = store.cache_get(f"seq:{assembly}:{region}", max_age=90 * 86400)
    if cached:
        return cached
    try:
        seq = (net.fetch_json(url) or {}).get("seq", "")
    except net.HttpError as exc:
        raise SequenceError(
            f"Could not fetch reference sequence from Ensembl: {exc}. "
            "Install the local genome from Settings to work offline.") from exc
    if not seq:
        raise SequenceError(f"Ensembl returned no sequence for {region}.")
    store.cache_put(f"seq:{assembly}:{region}", seq)
    return seq


def repeat_fraction(seq: str) -> float:
    """Share of soft-masked (repeat) bases, 0-1."""
    if not seq:
        return 0.0
    return sum(1 for c in seq if c.islower()) / len(seq)


def gc_fraction(seq: str) -> float:
    if not seq:
        return 0.0
    up = seq.upper()
    return (up.count("G") + up.count("C")) / len(up)
