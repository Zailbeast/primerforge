"""Search engines: NCBI BLAST+ (Windows) and UCSC BLAT (WSL).

Both produce hits in the shape Ensembl stores for its results table:

  qid qstart qend qori qframe   query name, 1-based span, orientation, frame
  tid tstart tend tori tframe   subject (chromosome, transcript or protein)
  gid gstart gend gori          genomic location of the hit
  score evalue pident len       bit score, E-value, % identity, alignment length
  aln                           BTOP alignment string
  qseq sseq                     aligned query and subject, gaps as '-'

Starts are always <= ends; orientation lives in the *ori fields. Aligned strings
are written query-forward, with the subject reverse-complemented for minus-strand
hits, which is the convention BLAST+ itself uses.
"""
from __future__ import annotations

import math
import os
import re
import subprocess
from pathlib import Path

from . import annotation, blastconf, blatserver, config, species, transcripts, wsl
from .blastconf import GENOMIC_SOURCES, SEARCH_TYPE_BY_VALUE

NO_WINDOW = getattr(subprocess, "CREATE_NO_WINDOW", 0)
TAB_FIELDS = ("qseqid sseqid pident length mismatch gapopen qstart qend sstart send "
              "evalue bitscore score qframe sframe qlen slen btop qseq sseq stitle")

# Download formats. Keys are what the UI offers; BLAST+ ones are rendered from the
# ASN.1 archive with blast_formatter on request.
NCBI_FORMATS = {
    "txt":  ("0", "BLAST pairwise text", "text/plain", "txt"),
    "tab":  ("7 std btop", "BLAST tabular with comment lines", "text/tab-separated-values", "tsv"),
    "xml":  ("5", "BLAST XML", "application/xml", "xml"),
    "json": ("15", "BLAST JSON", "application/json", "json"),
    "csv":  ("10", "BLAST CSV", "text/csv", "csv"),
    "sam":  ("17", "SAM (BLASTN only)", "text/plain", "sam"),
    "asn":  (None, "BLAST archive (ASN.1)", "text/plain", "asn"),
}
BLAT_FORMATS = {
    "blast": ("blat.blast.txt", "BLAT output, NCBI BLAST style", "text/plain", "txt"),
    "blast8": ("blat.blast8.tsv", "BLAT tabular (blast8)", "text/tab-separated-values", "tsv"),
    "psl": ("blat.psl", "PSL", "text/plain", "psl"),
    "axt": ("blat.axt", "AXT alignments", "text/plain", "axt"),
}


class EngineError(RuntimeError):
    pass


class Cancelled(RuntimeError):
    pass


_COMP = str.maketrans("ACGTUNacgtunRYKMSWBDHVrykmswbdhv-", "TGCAANtgcaanYRMKSWVHDByrmkswvhdb-")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def btop(qseq: str, sseq: str) -> str:
    """BLAST trace-back operations: match run lengths and query/subject pairs."""
    out, run = [], 0
    for q, s in zip(qseq.upper(), sseq.upper()):
        if q == s and q != "-":
            run += 1
        else:
            if run:
                out.append(str(run))
                run = 0
            out.append(q + s)
    if run:
        out.append(str(run))
    return "".join(out)


def threads() -> int:
    s = config.settings()
    if int(s.get("blast_threads") or 0) > 0:
        return int(s["blast_threads"])
    workers = max(1, int(s.get("blast_workers") or 1))
    return max(1, (os.cpu_count() or 4) // workers)


def _popen_wait(cmd: list[str], on_proc, cwd: str | None = None, timeout: float | None = None,
                via_wsl: bool = False, wsl_cwd: str | None = None) -> tuple[int, str, str]:
    if via_wsl:
        proc = wsl.popen(cmd, cwd=wsl_cwd, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, text=True, encoding="utf-8",
                         errors="replace")
    else:
        proc = subprocess.Popen([str(c) for c in cmd], cwd=cwd, stdin=subprocess.DEVNULL,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                                encoding="utf-8", errors="replace", creationflags=NO_WINDOW)
    if on_proc:
        on_proc(proc)
    try:
        out, err = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        proc.kill()
        out, err = proc.communicate()
        raise EngineError("The search took too long and was stopped.")
    finally:
        if on_proc:
            on_proc(None)
    return proc.returncode, out or "", err or ""


# ---------------------------------------------------------------------------
# NCBI BLAST+
# ---------------------------------------------------------------------------
def run_ncbi(job: dict, ticket: dict, manifest: dict, workdir: Path, on_proc=None,
             is_cancelled=lambda: False) -> dict:
    search = SEARCH_TYPE_BY_VALUE[ticket["search_type"]]
    program = search["program"]
    exe = config.blast_bin(program)
    if not exe:
        raise EngineError("NCBI BLAST+ is not installed (Settings > Install local genome).")
    db, mask = species.blast_db(manifest, ticket["source"])
    if not db:
        raise EngineError(f"{species.source_label(ticket['source'])} is not installed for "
                          f"{manifest.get('display_name', manifest['name'])}.")
    configs = ticket["configs"]
    notes: list[str] = []
    archive = workdir / "blast.asn"
    query = workdir / "input.fa"
    extra: list[str] = []
    if configs.get("repeat_mask") == "1" and program in ("blastn", "blastx", "tblastx"):
        query, note = _repeat_mask_query(query, manifest, on_proc)
        if query.name != "input.fa":
            extra.append("-lcase_masking")
        if note:
            notes.append(note)
    if is_cancelled():
        raise Cancelled()
    cmd = [exe, "-db", db, "-query", str(query), "-outfmt", "11",
           "-out", str(archive), "-num_threads", str(threads())]
    cmd += blastconf.blast_args(ticket["search_type"], configs, ticket.get("config_set", ""))
    cmd += extra
    if mask:
        cmd += ["-db_hard_mask" if ticket["source"] == "LATESTGP_MASKED" else "-db_soft_mask",
                mask]

    code, _out, err = _popen_wait(cmd, on_proc)
    if is_cancelled():
        raise Cancelled()
    if code != 0:
        raise EngineError(_blast_error(err))
    for line in err.splitlines():
        if line.startswith("Warning"):
            notes.append(line.strip())

    tsv = workdir / "hits.tsv"
    code, _out, err = _popen_wait([config.blast_bin("blast_formatter"), "-archive", str(archive),
                                   "-outfmt", f"6 {TAB_FIELDS}", "-out", str(tsv)], on_proc)
    if code != 0:
        raise EngineError(_blast_error(err))

    hits = []
    genomic = ticket["source"] in GENOMIC_SOURCES
    kind = "pep" if ticket["source"] == "PEP_ALL" else "cdna"
    with open(tsv, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 20:
                continue
            hits.append(_ncbi_hit(cols, job, program))
    hits.sort(key=lambda h: (-h["score"], h["evalue"]))
    # Ensembl REST is only consulted for species with no local annotation, and then
    # only for the best hits, so a batch of unmappable accessions cannot stall a job.
    rest_budget = 0 if annotation.has_local(manifest) else 25
    for h in hits:
        if genomic:
            h.update(gid=h["tid"], gstart=h["tstart"], gend=h["tend"], gori=h["tori"])
        else:
            if _map_transcript_hit(h, manifest, kind, allow_rest=rest_budget > 0):
                rest_budget -= 1
            _label_identifiers(h, manifest)
    tsv.unlink(missing_ok=True)
    return {"hits": hits, "notes": notes, "command": _cmdline(cmd), "engine": "ncbi",
            "files": {"asn": archive.name}}


def _repeat_mask_query(query: Path, manifest: dict, on_proc) -> tuple[Path, str | None]:
    """Soft-mask repeats in the query with WindowMasker ("Filter query sequences using
    RepeatMasker"). Masking the query up front, then searching with -lcase_masking,
    works for BLASTN, BLASTX and TBLASTX alike; BLAST+ itself only accepts a
    WindowMasker database for some programs."""
    wm = (manifest.get("repeats") or {}).get("windowmasker")
    tool = config.blast_bin("windowmasker")
    if not wm or not Path(wm).exists() or not tool:
        return query, ("Repeat masking of the query was requested but the species' repeat "
                       "filter is not installed, so the query was searched unmasked.")
    masked = query.with_name("input.masked.fa")
    code, _out, err = _popen_wait([tool, "-ustat", wm, "-in", str(query), "-outfmt", "fasta",
                                   "-out", str(masked)], on_proc, timeout=600)
    if code != 0 or not masked.exists():
        return query, f"Repeat masking failed ({_blast_error(err)}); the query was searched unmasked."
    body = "".join(ln.strip() for ln in masked.read_text(encoding="utf-8").splitlines()
                   if not ln.startswith(">"))
    lower = sum(1 for ch in body if ch.islower())
    pct = 100.0 * lower / len(body) if body else 0.0
    if lower == len(body) and body:
        return query, ("The whole query is repetitive, so repeat masking was skipped to "
                       "avoid an empty search.")
    return masked, (f"Repeats were masked in {pct:.0f}% of the query (WindowMasker)."
                    if lower else None)


def _blast_error(err: str) -> str:
    lines = [ln.strip() for ln in err.splitlines() if ln.strip()]
    useful = [ln for ln in lines if "error" in ln.lower()] or lines
    return " ".join(useful)[:600] or "BLAST failed without an error message."


def _cmdline(cmd: list) -> str:
    return " ".join(f'"{c}"' if " " in str(c) else str(c) for c in cmd)


def _ncbi_hit(cols: list[str], job: dict, program: str) -> dict:
    (qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send,
     evalue, bitscore, score, qframe, sframe, qlen, slen, bt, qseq, sseq) = cols[:20]
    stitle = cols[20] if len(cols) > 20 else ""
    qs, qe, ss, se = int(qstart), int(qend), int(sstart), int(send)
    qf, sf = int(qframe or 0), int(sframe or 0)
    qori = -1 if (qf < 0 or qs > qe) else 1
    tori = -1 if (sf < 0 or ss > se) else 1
    hit = {
        "qid": job.get("seq_desc") or qseqid, "qstart": min(qs, qe), "qend": max(qs, qe),
        "qori": qori, "qframe": qf, "qlen": int(qlen),
        "tid": _clean_seqid(sseqid), "tstart": min(ss, se), "tend": max(ss, se),
        "tori": tori, "tframe": sf, "slen": int(slen),
        "score": float(bitscore), "raw_score": int(float(score)),
        "evalue": float(evalue), "pident": round(float(pident), 2), "len": int(length),
        "mismatch": int(mismatch), "gapopen": int(gapopen),
        "aln": bt, "qseq": qseq, "sseq": sseq, "tdesc": stitle,
        "program": program,
    }
    hit.update(_title_meta(stitle))
    return hit


_DB_TAGS = {"ref", "gb", "emb", "dbj", "sp", "tr", "pir", "prf", "tpg", "tpe", "tpd", "pdb"}


def _clean_seqid(sid: str) -> str:
    """Plain accession from a BLAST sequence id ('ref|NM_001008216.2|' -> 'NM_001008216.2')."""
    if "|" not in sid:
        return sid
    parts = [p for p in sid.split("|") if p]
    if not parts:
        return sid
    if parts[0] in _DB_TAGS and len(parts) > 1:
        return parts[1]
    if parts[0] in ("lcl", "gnl"):
        return parts[-1]
    return parts[-1]


_TITLE_KV = re.compile(r"(\w+):(\S+)")


def _title_meta(stitle: str) -> dict:
    """Gene, transcript and span annotations embedded in Ensembl FASTA headers."""
    if not stitle:
        return {}
    meta = dict(_TITLE_KV.findall(stitle.split(" description:")[0]))
    out = {}
    if meta.get("gene"):
        out["gene_hit"] = {"id": meta["gene"].split(".")[0],
                           "name": meta.get("gene_symbol") or meta["gene"].split(".")[0],
                           "biotype": meta.get("gene_biotype")}
    if meta.get("transcript"):
        out["transcript_id"] = meta["transcript"]
    m = re.search(r"(?:chromosome|scaffold|primary_assembly|supercontig|contig):[^:\s]+:"
                  r"([^:\s]+):(\d+):(\d+):(-?1)", stitle)
    if m:
        out["tspan"] = {"chrom": m.group(1), "start": int(m.group(2)), "end": int(m.group(3)),
                        "strand": int(m.group(4))}
    desc = stitle.split(" description:", 1)
    if len(desc) == 2:
        out["tdescription"] = re.sub(r"\s*\[Source:.*$", "", desc[1]).strip()
    return out


def _label_identifiers(h: dict, manifest: dict) -> None:
    """Add the matching RefSeq/Ensembl accession and MANE status to a transcript hit."""
    try:
        ids = transcripts.annotate(manifest, h["tid"])
    except Exception:                                            # noqa: BLE001
        return
    if ids:
        h["ids"] = ids
        if not h.get("gene_hit") and ids.get("symbol"):
            h["gene_hit"] = {"id": ids.get("gene_id") or "", "name": ids["symbol"],
                             "biotype": ids.get("moltype")}
        if not h.get("tdescription") and ids.get("description"):
            h["tdescription"] = ids["description"]


def _map_transcript_hit(h: dict, manifest: dict, kind: str, allow_rest: bool = True) -> bool:
    """Place a transcript or protein hit on the genome. Returns whether REST was used."""
    mapped = None
    used_rest = False
    target = h["tid"]
    if not target.upper().startswith("ENS"):
        # RefSeq accession: map through the Ensembl transcript it is matched to.
        try:
            target = transcripts.ensembl_for(manifest, target) or target
        except Exception:                                        # noqa: BLE001
            pass
    try:
        mapped = annotation.map_to_genome(manifest, target, kind, h["tstart"], h["tend"],
                                          allow_rest=allow_rest)
        used_rest = allow_rest and bool(mapped) and not annotation.has_local(manifest)
    except Exception:                                            # noqa: BLE001
        mapped = None
    if mapped:
        h.update(gid=mapped["chrom"], gstart=mapped["start"], gend=mapped["end"],
                 gori=mapped["strand"] * h["tori"], gblocks=mapped.get("blocks"))
        if not h.get("gene_hit") and mapped.get("gene_id"):
            h["gene_hit"] = {"id": mapped["gene_id"], "name": mapped.get("gene_name") or "",
                             "biotype": None}
    else:
        # No exon mapping: fall back to the span in the FASTA header (Ensembl) or the
        # gene span from MANE (RefSeq), marked as approximate.
        span = h.get("tspan")
        if not span:
            try:
                span = transcripts.span_for(manifest, h["tid"])
            except Exception:                                    # noqa: BLE001
                span = None
        if span:
            h.update(gid=span["chrom"], gstart=span["start"], gend=span["end"],
                     gori=span["strand"] * h["tori"], gapprox=True)
        else:
            h.update(gid=None, gstart=None, gend=None, gori=None)
    return used_rest


# ---------------------------------------------------------------------------
# BLAT via WSL
# ---------------------------------------------------------------------------
def run_blat(job: dict, ticket: dict, manifest: dict, workdir: Path, on_proc=None,
             is_cancelled=lambda: False) -> dict:
    name = manifest["name"]
    if not species.blat_installed(manifest):
        raise EngineError(f"BLAT is not installed for {manifest.get('display_name', name)}.")
    try:
        blatserver.ensure_running(name)
    except blatserver.BlatServerError as exc:
        raise EngineError(f"BLAT server could not start: {exc}") from exc
    if is_cancelled():
        raise Cancelled()

    b = manifest["blat"]
    c = ticket["configs"]
    wdir = wsl.to_wsl(workdir)
    common = [f"-minScore={c.get('min_score', '0')}",
              f"-minIdentity={c.get('min_identity', '0')}",
              f"-maxIntron={c.get('max_intron', '750000')}"]
    client = f"{b['bindir']}/gfClient"
    commands = []
    for fmt, fname in (("axt", BLAT_FORMATS["axt"][0]), ("psl", BLAT_FORMATS["psl"][0]),
                       ("blast", BLAT_FORMATS["blast"][0])):
        cmd = [client, f"-out={fmt}", *common, blatserver.HOST, str(b["port"]), b["dir"],
               "input.fa", fname]
        commands.append(" ".join(cmd))
        code, out, err = _popen_wait(cmd, on_proc, via_wsl=True, wsl_cwd=wdir, timeout=3600)
        if is_cancelled():
            raise Cancelled()
        if code != 0:
            raise EngineError(f"gfClient failed: {(err or out).strip()[:500]}")

    query_len = len(job["sequence"])
    hits = parse_axt((workdir / BLAT_FORMATS["axt"][0]).read_text(encoding="utf-8",
                                                                 errors="replace"),
                     job.get("seq_desc") or "query", query_len)
    hits.sort(key=lambda h: (-h["raw_score"], h["evalue"]))
    try:
        max_e = float(c.get("evalue", "1e-1"))
    except ValueError:
        max_e = 0.1
    notes = []
    kept = [h for h in hits if h["evalue"] <= max_e]
    if len(kept) < len(hits):
        notes.append(f"{len(hits) - len(kept)} BLAT alignment(s) above the E-value cut-off "
                     f"of {c.get('evalue')} were not reported.")
    limit = int(c.get("max_target_seqs", "100") or 100)
    if len(kept) > limit:
        notes.append(f"Showing the best {limit} of {len(kept)} alignments "
                     "(maximum number of hits to report).")
        kept = kept[:limit]
    for h in kept:
        h.update(gid=h["tid"], gstart=h["tstart"], gend=h["tend"], gori=h["tori"])
    _write_blast8(workdir / BLAT_FORMATS["blast8"][0], kept)
    return {"hits": kept, "notes": notes, "command": "\n".join(commands), "engine": "blat",
            "files": {k: v[0] for k, v in BLAT_FORMATS.items()}}


def blat_score_to_bits(score: int) -> int:
    return round(score * 0.0205)


def blat_score_to_evalue(score: int) -> float:
    """BLAT's own conversion for its NCBI-style output (kent src/lib/blastOut.c)."""
    return 3.0e9 * math.exp(-(score * 0.0205) * math.log(2))


def parse_axt(text: str, query_name: str, query_len: int) -> list[dict]:
    hits = []
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        head = lines[i].split()
        if len(head) < 9 or not head[0].isdigit():
            i += 1
            continue
        # Line 2 is the target (genome), line 3 the query; BLAT lower-cases repeats.
        tseq = lines[i + 1].strip().upper() if i + 1 < len(lines) else ""
        qseq = lines[i + 2].strip().upper() if i + 2 < len(lines) else ""
        i += 3
        _n, tname, ts, te, _qname, qs, qe, strand, score = head[:9]
        ts, te, qs, qe, score = int(ts), int(te), int(qs), int(qe), int(score)
        if strand == "-":
            # axt gives minus-strand query coordinates on the reverse complement;
            # flip to query-forward with the genome reverse-complemented instead.
            qs, qe = query_len - qe + 1, query_len - qs + 1
            qseq, tseq = revcomp(qseq), revcomp(tseq)
        size = len(qseq)
        matches = sum(1 for a, b in zip(qseq.upper(), tseq.upper()) if a == b and a != "-")
        gaps = qseq.count("-") + tseq.count("-")
        gap_opens = len(re.findall(r"-+", "".join(
            "-" if a == "-" or b == "-" else "x" for a, b in zip(qseq, tseq))))
        hits.append({
            "qid": query_name, "qstart": qs, "qend": qe, "qori": 1, "qframe": 0,
            "qlen": query_len, "tid": tname, "tstart": ts, "tend": te,
            "tori": -1 if strand == "-" else 1, "tframe": 0,
            "score": float(blat_score_to_bits(score)), "raw_score": score,
            "evalue": blat_score_to_evalue(score),
            "pident": round(100.0 * matches / size, 2) if size else 0.0, "len": size,
            "mismatch": size - matches - gaps, "gapopen": gap_opens,
            "aln": btop(qseq, tseq), "qseq": qseq, "sseq": tseq, "program": "blat",
        })
    return hits


def _write_blast8(path: Path, hits: list[dict]) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for h in hits:
            ss, se = (h["tstart"], h["tend"]) if h["tori"] == 1 else (h["tend"], h["tstart"])
            fh.write("\t".join(map(str, [
                h["qid"].split()[0] if h["qid"] else "query", h["tid"], f"{h['pident']:.2f}",
                h["len"], h["mismatch"], h["gapopen"], h["qstart"], h["qend"], ss, se,
                f"{h['evalue']:.1e}", f"{int(h['score'])}.0"])) + "\n")


# ---------------------------------------------------------------------------
# Downloads
# ---------------------------------------------------------------------------
def render_ncbi_format(workdir: Path, fmt: str) -> Path:
    spec = NCBI_FORMATS.get(fmt)
    if not spec:
        raise EngineError(f"Unknown format '{fmt}'.")
    archive = workdir / "blast.asn"
    if not archive.exists():
        raise EngineError("The BLAST archive for this job is missing.")
    if spec[0] is None:
        return archive
    out = workdir / f"download_{fmt}.{spec[3]}"
    if out.exists() and out.stat().st_mtime >= archive.stat().st_mtime:
        return out
    code, _o, err = _popen_wait([config.blast_bin("blast_formatter"), "-archive", str(archive),
                                 "-outfmt", spec[0], "-out", str(out)], None, timeout=900)
    if code != 0:
        out.unlink(missing_ok=True)
        raise EngineError(_blast_error(err))
    return out
