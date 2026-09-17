"""Text sequence views for a single BLAST/BLAT hit.

Ensembl offers three views per hit, each with a "Configure this page" panel:

  Alignment         query / match line / subject, with exon and variant markup
  Query sequence    the whole query with this and other HSPs highlighted
  Genomic sequence  the hit's genomic region with flanks, HSPs, exons, variants

Option names and values follow Ensembl's ViewConfig::Blast::* modules. Views are
rendered to HTML here: long sequences (200 kb queries, 10 kb flanks) are cheaper
to build as run-length-merged spans in Python than as per-base DOM in the browser.
"""
from __future__ import annotations

import functools
import html
import threading
from pathlib import Path

from . import annotation, blastconf, species
from .blastconf import GENOMIC_SOURCES
from .refseq import FastaIndex

MISMATCH, GAP, HSP_SEL, HSP_OTHER, EXON, VARIANT, REPEAT = 1, 2, 4, 8, 16, 32, 64
CLASS_NAMES = [(MISMATCH, "sv-mm"), (GAP, "sv-gap"), (HSP_SEL, "sv-hsp-sel"),
               (HSP_OTHER, "sv-hsp"), (EXON, "sv-exon"), (VARIANT, "sv-var"),
               (REPEAT, "sv-rep")]

WIDTHS = [(w, f"{w} bps per line") for w in range(30, 160, 10)]
FLANKS = [(f, f"{f:,} bp") for f in (0, 50, 100, 200, 300, 400, 500, 1000, 2000, 5000, 10000)]
LINE_NUMBERING = [("slice", "Relative to coordinate systems"),
                  ("sequence", "Relative to this sequence"), ("off", "None")]

VIEW_OPTIONS = {
    "alignment": [
        ("display_width", "Number of base pairs per row", WIDTHS, "60"),
        ("align_display", "Alignments display", [("off", "Off"),
                                                 ("line", "Mark matching bp with lines"),
                                                 ("dot", "Mark matching bp with dots")], "line"),
        ("exon_display", "Show exons", [("core", "Ensembl exons"), ("off", "None")], "core"),
        ("exon_ori", "Orientation of exons", [("fwd", "Forward only"), ("rev", "Reverse only"),
                                              ("all", "Both orientations")], "all"),
        ("snp_display", "Show variants", [("on", "Yes"), ("off", "No")], "off"),
        ("line_numbering", "Line numbering", LINE_NUMBERING, "slice"),
        ("title_display", "Display pop-up information on mouseover",
         [("yes", "Yes"), ("off", "No")], "yes"),
    ],
    "query": [
        ("display_width", "Number of base pairs per row", WIDTHS, "60"),
        ("hsp_display", "Alignment Markup", [("all", "All alignments"),
                                             ("sel", "Selected alignments only"),
                                             ("off", "No alignment markup")], "all"),
        ("line_numbering", "Line numbering", [("slice", "Relative to this sequence"),
                                              ("off", "None")], "slice"),
    ],
    "genomic": [
        ("display_width", "Number of base pairs per row", WIDTHS, "60"),
        ("flank5_display", "5' Flanking sequence", FLANKS, "300"),
        ("flank3_display", "3' Flanking sequence", FLANKS, "300"),
        ("orientation", "Orientation", [("fc", "Forward relative to coordinate system"),
                                        ("rc", "Reverse relative to coordinate system"),
                                        ("fa", "Forward relative to selected alignment")], "fa"),
        ("hsp_display", "Alignment Markup", [("all", "All alignments"),
                                             ("sel", "Selected alignments only"),
                                             ("off", "No alignment markup")], "all"),
        ("exon_display", "Show exons", [("core", "Ensembl exons"), ("off", "None")], "core"),
        ("exon_ori", "Orientation of exons", [("fwd", "Forward only"), ("rev", "Reverse only"),
                                              ("all", "Both orientations")], "all"),
        ("snp_display", "Show variants", [("on", "Yes"), ("off", "No")], "off"),
        ("repeat_display", "Show repeats (soft-masked) in lower case",
         [("on", "Yes"), ("off", "No")], "off"),
        ("line_numbering", "Line numbering", LINE_NUMBERING, "slice"),
        ("title_display", "Display pop-up information on mouseover",
         [("yes", "Yes"), ("off", "No")], "yes"),
    ],
}


class ViewError(RuntimeError):
    pass


def lookup(kind: str, manifest: dict, chrom: str, start: int, end: int,
           pending: list[dict]) -> list[dict]:
    """Annotation from the local index or the REST cache; slow network lookups are
    recorded in `pending` for the browser to fetch in the background."""
    try:
        return lookup_function(kind)(manifest, chrom, start, end, cached_only=True)
    except annotation.NotCached:
        pending.append({"kind": kind, "chrom": chrom, "start": start, "end": end})
        return []


def lookup_function(kind: str):
    """The annotation call behind a lookup kind: exons, variants or somatic."""
    if kind == "exons":
        return annotation.exons_in_region
    if kind == "somatic":
        return functools.partial(annotation.variants_in_region, somatic=True)
    return annotation.variants_in_region


def options(view: str, params: dict) -> dict:
    out = {}
    for name, _label, values, default in VIEW_OPTIONS[view]:
        value = str(params.get(name, default))
        out[name] = value if value in {str(v) for v, _ in values} else default
    return out


def option_form(view: str, current: dict) -> list[dict]:
    return [{"name": n, "label": label, "value": current[n],
             "values": [{"value": str(v), "caption": c} for v, c in values]}
            for n, label, values, _d in VIEW_OPTIONS[view]]


def _ori(v) -> str:
    return "+" if (v or 1) > 0 else "-"


def hit_summary(ticket: dict, hit: dict) -> list[tuple[str, str]]:
    method = blastconf.SEARCH_TYPE_BY_VALUE[ticket["search_type"]]["method"]
    kind = "BLAT" if ticket["search_type"] == blastconf.BLAT_VALUE else "BLAST"
    rows = [(f"{kind} type", method),
            ("Query location", f"{hit['qid']} {hit['qstart']:,} to {hit['qend']:,} "
                               f"({_ori(hit['qori'])})"),
            ("Database location", f"{hit['tid']} {hit['tstart']:,} to {hit['tend']:,} "
                                  f"({_ori(hit['tori'])})")]
    if hit.get("gid"):
        approx = " (transcript span)" if hit.get("gapprox") else ""
        rows.append(("Genomic location", f"{hit['gid']} {hit['gstart']:,} to "
                                         f"{hit['gend']:,} ({_ori(hit['gori'])}){approx}"))
    rows += [("Alignment score", f"{hit['score']:g}"), ("E-value", f"{hit['evalue']:.3g}"),
             ("Alignment length", f"{hit['len']:,}"),
             ("Percentage identity", f"{hit['pident']:.2f}")]
    return rows


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------
def _classes(mask: int) -> str:
    return " ".join(name for bit, name in CLASS_NAMES if mask & bit)


def _segments(seq: str, marks: list[int], titles: list[str | None] | None,
              lo: int, hi: int) -> str:
    out = []
    i = lo
    while i < hi:
        m = marks[i] if marks else 0
        t = titles[i] if titles else None
        j = i + 1
        while j < hi and (marks[j] if marks else 0) == m and (titles[j] if titles else None) == t:
            j += 1
        text = html.escape(seq[i:j])
        if m or t:
            attr = f' class="{_classes(m)}"' if m else ""
            if t:
                attr += f' title="{html.escape(t, quote=True)}"'
            out.append(f"<span{attr}>{text}</span>")
        else:
            out.append(text)
        i = j
    return "".join(out)


def _num(value, width: int) -> str:
    return f'<span class="sv-ln">{"" if value is None else f"{value:,}":>{width}}</span>'


class _Titles:
    """Per-position pop-up text, built lazily (most positions have none)."""

    def __init__(self, n: int, enabled: bool):
        self.enabled = enabled
        self.parts: dict[int, list[str]] = {}
        self.n = n

    def add(self, pos: int, text: str) -> None:
        if self.enabled and 0 <= pos < self.n:
            lst = self.parts.setdefault(pos, [])
            if text not in lst and len(lst) < 4:
                lst.append(text)

    def array(self) -> list[str | None] | None:
        if not self.enabled or not self.parts:
            return None
        arr: list[str | None] = [None] * self.n
        for p, texts in self.parts.items():
            arr[p] = "\n".join(texts)
        return arr


def _exon_filter(exons: list[dict], ori: str, display_strand: int) -> list[dict]:
    if ori == "fwd":
        return [e for e in exons if (e["strand"] or 1) == display_strand]
    if ori == "rev":
        return [e for e in exons if (e["strand"] or 1) != display_strand]
    return exons


def _exon_labels(exons: list[dict]) -> dict[tuple[int, int], str]:
    grouped: dict[tuple[int, int], list[dict]] = {}
    for e in exons:
        grouped.setdefault((e["start"], e["end"]), []).append(e)
    out = {}
    for span, group in grouped.items():
        genes = sorted({g["gene_name"] or "" for g in group} - {""})
        first = group[0]
        label = f"Exon {first.get('number') or ''} of {first.get('transcript_name') or first['transcript_id']}"
        if len(group) > 1:
            label += f" (+{len(group) - 1} more transcript{'s' if len(group) > 2 else ''})"
        if genes:
            label = f"{', '.join(genes)}: {label}"
        out[span] = label
    return out


# ---------------------------------------------------------------------------
# Alignment view
# ---------------------------------------------------------------------------
def alignment_view(ticket: dict, job: dict, hit: dict, manifest: dict, opts: dict) -> dict:
    qseq, sseq = hit.get("qseq") or "", hit.get("sseq") or ""
    if not qseq:
        raise ViewError("This hit has no stored alignment.")
    program = hit.get("program") or blastconf.SEARCH_TYPE_BY_VALUE[ticket["search_type"]]["program"]
    is_protein = ticket["query_type"] == "peptide" or ticket["db_type"] == "peptide" \
        or program == "tblastx"
    qstep = 3 if program in ("blastx", "tblastx") else 1
    sstep = 3 if program in ("tblastn", "tblastx") else 1
    width = int(opts["display_width"])
    n = len(qseq)
    titles_on = opts["title_display"] == "yes"

    qmarks = [0] * n
    smarks = [0] * n
    mid = []
    for i, (a, b) in enumerate(zip(qseq.upper(), sseq.upper())):
        if a == "-" or b == "-":
            qmarks[i] |= GAP if a == "-" else 0
            smarks[i] |= GAP if b == "-" else 0
            mid.append(" ")
        elif a == b:
            mid.append("|" if opts["align_display"] != "dot" else ".")
        else:
            qmarks[i] |= MISMATCH
            smarks[i] |= MISMATCH
            mid.append(" ")
    midline = "".join(mid)

    # Coordinates of every alignment column on each sequence (None in gaps).
    def coords(seq: str, first: int, ori: int, step: int) -> list[int | None]:
        out, pos = [], first
        for ch in seq:
            if ch == "-":
                out.append(None)
            else:
                out.append(pos)
                pos += step * ori
        return out

    qori, tori = hit["qori"] or 1, hit["tori"] or 1
    qc = coords(qseq, hit["qstart"] if qori > 0 else hit["qend"], qori, qstep)
    sc = coords(sseq, hit["tstart"] if tori > 0 else hit["tend"], tori, sstep)

    notes = []
    pending: list[dict] = []
    stitles = _Titles(n, titles_on)
    genomic = ticket["source"] in GENOMIC_SOURCES
    legend = [("sv-mm", "Mismatch"), ("sv-gap", "Gap")]
    if genomic and opts["exon_display"] != "off":
        exons = _exon_filter(lookup("exons", manifest, hit["tid"], hit["tstart"], hit["tend"],
                                     pending), opts["exon_ori"], tori)
        labels = _exon_labels(exons)
        spans = sorted(labels.items())
        for i, g in enumerate(sc):
            if g is None:
                continue
            centre = g + (tori if sstep == 3 else 0)
            for (s, e), label in spans:
                if s <= centre <= e:
                    smarks[i] |= EXON
                    stitles.add(i, label)
                    break
        legend.append(("sv-exon", "Exon"))
    elif not genomic and opts["exon_display"] != "off":
        notes.append("Exon and variant markup is shown for hits against the genome.")
    if genomic and opts["snp_display"] == "on":
        variants = lookup("variants", manifest, hit["tid"], hit["tstart"], hit["tend"], pending)
        by_pos = {}
        for v in variants:
            for p in range(v["start"], v["end"] + 1):
                by_pos.setdefault(p, v)
        for i, g in enumerate(sc):
            if g is not None and g in by_pos:
                v = by_pos[g]
                smarks[i] |= VARIANT
                stitles.add(i, f"{v['id']} {v['alleles']} {v['consequence']}".strip())
        legend.append(("sv-var", "Variant"))

    st_arr = stitles.array()
    numbering = opts["line_numbering"]
    label_w = 8
    num_w = max(len(f"{max(hit['qend'], hit['tend']):,}"), len(f"{n:,}")) + 1
    lines = []
    qcount = scount = 0
    for lo in range(0, n, width):
        hi = min(n, lo + width)

        def span(clist, count_before, seg, step=1, ori=1):
            real = [c for c in clist[lo:hi] if c is not None]
            if numbering == "slice":
                if not real:
                    return None, None
                # A translated residue covers a codon: end on its last nucleotide.
                return real[0], real[-1] + (step - 1) * ori
            if numbering == "sequence":
                k = sum(1 for ch in seg if ch != "-")
                return (count_before + 1 if k else count_before), count_before + k
            return None, None

        qs_, qe_ = span(qc, qcount, qseq[lo:hi], qstep, qori)
        ss_, se_ = span(sc, scount, sseq[lo:hi], sstep, tori)
        qcount += sum(1 for ch in qseq[lo:hi] if ch != "-")
        scount += sum(1 for ch in sseq[lo:hi] if ch != "-")
        nums = numbering != "off"
        q_line = (f'<span class="sv-lab">{"Query":<{label_w}}</span>'
                  + (_num(qs_, num_w) + " " if nums else "")
                  + _segments(qseq, qmarks, None, lo, hi)
                  + (" " + _num(qe_, 0) if nums else ""))
        rows = [q_line]
        if opts["align_display"] != "off":
            rows.append(" " * (label_w + (num_w + 1 if nums else 0))
                        + f'<span class="sv-mid">{html.escape(midline[lo:hi])}</span>')
        s_line = (f'<span class="sv-lab">{"Subject":<{label_w}}</span>'
                  + (_num(ss_, num_w) + " " if nums else "")
                  + _segments(sseq, smarks, st_arr, lo, hi)
                  + (" " + _num(se_, 0) if nums else ""))
        rows.append(s_line)
        lines.append("\n".join(rows))
    title = f"Alignment of {hit['qid']} with {hit['tid']}:{hit['tstart']:,}-{hit['tend']:,}"
    return {"html": "\n\n".join(lines), "legend": legend, "notes": notes, "title": title,
            "is_protein": is_protein, "pending": pending}


# ---------------------------------------------------------------------------
# Query sequence view
# ---------------------------------------------------------------------------
def query_view(ticket: dict, job: dict, hit: dict, all_hits: list[dict], opts: dict) -> dict:
    seq = job["sequence"]
    n = len(seq)
    marks = [0] * n
    titles = _Titles(n, True)
    legend = []
    if opts["hsp_display"] != "off":
        shown = [hit] if opts["hsp_display"] == "sel" else \
            [h for h in all_hits[:1000] if h["idx"] != hit["idx"]] + [hit]
        for h in shown:
            bit = HSP_SEL if h["idx"] == hit["idx"] else HSP_OTHER
            where = f"{h['gid']}:{h['gstart']:,}-{h['gend']:,}" if h.get("gid") else \
                f"{h['tid']}:{h['tstart']:,}-{h['tend']:,}"
            label = (f"{'Selected hit' if bit == HSP_SEL else 'Hit'} {h['idx'] + 1}: {where} "
                     f"({_ori(h.get('gori') or h['tori'])}), {h['pident']:.1f}% ID")
            for p in range(max(0, h["qstart"] - 1), min(n, h["qend"])):
                if bit == HSP_SEL:
                    marks[p] = (marks[p] & ~HSP_OTHER) | HSP_SEL
                elif not marks[p] & HSP_SEL:
                    marks[p] |= HSP_OTHER
            titles.add(max(0, h["qstart"] - 1), label)
        legend = [("sv-hsp-sel", "Selected alignment")]
        if opts["hsp_display"] == "all":
            legend.append(("sv-hsp", "Other alignments"))
    # Titles only at HSP starts would split spans oddly; apply them across each HSP run.
    arr = _run_titles(marks, titles)
    width = int(opts["display_width"])
    num_w = len(f"{n:,}") + 1
    lines = []
    for lo in range(0, n, width):
        hi = min(n, lo + width)
        body = _segments(seq, marks, arr, lo, hi)
        if opts["line_numbering"] == "slice":
            lines.append(f"{_num(lo + 1, num_w)} {body} {_num(hi, 0)}")
        else:
            lines.append(body)
    unit = "residues" if ticket["query_type"] == "peptide" else "bases"
    return {"html": "\n".join(lines), "legend": legend, "notes": [],
            "title": f"{job.get('seq_desc') or 'Query'} ({n:,} {unit})"}


def _run_titles(marks: list[int], titles: _Titles) -> list[str | None] | None:
    """Spread each HSP's start title over its highlighted run so hovering anywhere works."""
    if not titles.parts:
        return None
    arr: list[str | None] = [None] * len(marks)
    current = None
    for i, m in enumerate(marks):
        if i in titles.parts:
            current = "\n".join(titles.parts[i])
        if not m & (HSP_SEL | HSP_OTHER):
            current = None if i not in titles.parts else current
        arr[i] = current if m & (HSP_SEL | HSP_OTHER) else None
    return arr


# ---------------------------------------------------------------------------
# Genomic sequence view
# ---------------------------------------------------------------------------
_fai_cache: dict[str, FastaIndex] = {}
_fai_lock = threading.Lock()


def fasta_index(manifest: dict) -> FastaIndex | None:
    g = manifest.get("genome") or {}
    if not g.get("fai") or not Path(g["fai"]).exists():
        return None
    with _fai_lock:
        idx = _fai_cache.get(g["fasta"])
        if idx is None:
            idx = FastaIndex(Path(g["fasta"]), Path(g["fai"]))
            _fai_cache[g["fasta"]] = idx
    return idx


def genomic_region(hit: dict, opts: dict, contig_len: int | None = None) -> tuple[int, int, int]:
    f5, f3 = int(opts["flank5_display"]), int(opts["flank3_display"])
    gori = hit.get("gori") or 1
    if gori > 0:
        start, end = hit["gstart"] - f5, hit["gend"] + f3
    else:
        start, end = hit["gstart"] - f3, hit["gend"] + f5
    start = max(1, start)
    if contig_len:
        end = min(contig_len, end)
    ori = opts["orientation"]
    strand = 1 if ori == "fc" else -1 if ori == "rc" else gori * (hit.get("qori") or 1)
    return start, end, strand


def genomic_view(ticket: dict, job: dict, hit: dict, all_hits: list[dict], manifest: dict,
                 opts: dict) -> dict:
    if not hit.get("gid"):
        raise ViewError("This hit could not be placed on the genome.")
    idx = fasta_index(manifest)
    if idx is None:
        raise ViewError("The genome sequence for this species is not installed.")
    try:
        contig_len = idx.index[idx._resolve(hit["gid"])][0]
    except Exception as exc:                                    # noqa: BLE001
        raise ViewError(f"{hit['gid']} is not in the installed genome.") from exc
    start, end, strand = genomic_region(hit, opts, contig_len)
    seq = idx.fetch(hit["gid"], start, end)
    if opts["repeat_display"] != "on":
        seq = seq.upper()
    n = len(seq)
    if strand < 0:
        from .engines import revcomp
        seq = revcomp(seq)

    def to_i(g: int) -> int:
        return g - start if strand > 0 else end - g

    marks = [0] * n
    titles = _Titles(n, opts["title_display"] == "yes")
    legend = []

    if opts["repeat_display"] == "on":
        for i, ch in enumerate(seq):
            if ch.islower():
                marks[i] |= REPEAT
        legend.append(("sv-rep", "Repeat (lower case)"))

    if opts["hsp_display"] != "off":
        pool = [hit] if opts["hsp_display"] == "sel" else all_hits
        for h in pool:
            if h.get("gid") != hit["gid"] or h["gend"] < start or h["gstart"] > end:
                continue
            bit = HSP_SEL if h["idx"] == hit["idx"] else HSP_OTHER
            blocks = h.get("gblocks") or [(h["gstart"], h["gend"])]
            label = (f"{'Selected hit' if bit == HSP_SEL else 'Hit'} {h['idx'] + 1}: query "
                     f"{h['qstart']:,}-{h['qend']:,}, {h['pident']:.1f}% ID, "
                     f"E-value {h['evalue']:.2g}")
            for bs, be in blocks:
                for g in range(max(start, bs), min(end, be) + 1):
                    i = to_i(g)
                    if bit == HSP_SEL:
                        marks[i] = (marks[i] & ~HSP_OTHER) | HSP_SEL
                    elif not marks[i] & HSP_SEL:
                        marks[i] |= HSP_OTHER
                    if g in (bs, be) or g in (start, end):
                        titles.add(i, label)
        legend.append(("sv-hsp-sel", "Selected alignment"))
        if opts["hsp_display"] == "all":
            legend.append(("sv-hsp", "Other alignments"))

    notes = []
    pending: list[dict] = []
    if opts["exon_display"] != "off":
        exons = _exon_filter(lookup("exons", manifest, hit["gid"], start, end, pending),
                             opts["exon_ori"], strand)
        for (s, e), label in _exon_labels(exons).items():
            for g in range(max(start, s), min(end, e) + 1):
                i = to_i(g)
                marks[i] |= EXON
                titles.add(i, label)
        legend.append(("sv-exon", "Exon"))
        if not annotation.has_local(manifest):
            notes.append("Exons come from Ensembl REST; install the species annotation in "
                         "Settings to show them offline.")
    if opts["snp_display"] == "on":
        for v in lookup("variants", manifest, hit["gid"], start, end, pending):
            for g in range(max(start, v["start"]), min(end, v["end"]) + 1):
                i = to_i(g)
                marks[i] |= VARIANT
                titles.add(i, f"{v['id']} {v['alleles']} {v['consequence']}".strip())
        legend.append(("sv-var", "Variant"))

    arr = titles.array()
    width = int(opts["display_width"])
    numbering = opts["line_numbering"]
    num_w = len(f"{end:,}") + 1
    lines = []
    for lo in range(0, n, width):
        hi = min(n, lo + width)
        body = _segments(seq, marks, arr, lo, hi)
        if numbering == "slice":
            a = start + lo if strand > 0 else end - lo
            b = start + hi - 1 if strand > 0 else end - hi + 1
            lines.append(f"{_num(a, num_w)} {body} {_num(b, 0)}")
        elif numbering == "sequence":
            lines.append(f"{_num(lo + 1, num_w)} {body} {_num(hi, 0)}")
        else:
            lines.append(body)
    name = f"{hit['gid']}:{start}-{end}:{strand}"
    return {"html": "\n".join(lines), "legend": legend, "notes": notes,
            "title": f"{manifest.get('assembly', '')} {hit['gid']}:{start:,}-{end:,} "
                     f"({'forward' if strand > 0 else 'reverse'} strand)",
            "region": {"chrom": hit["gid"], "start": start, "end": end, "strand": strand},
            "fasta": blastconf.fasta(name, seq.upper()), "sequence": seq.upper(),
            "pending": pending}

