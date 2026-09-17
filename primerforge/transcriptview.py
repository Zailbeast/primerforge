"""Aligned transcript view: the entire genomic sequence with cDNA and protein in register.

A transcript is laid out as a list of columns, one per genomic base across the whole
gene (5' flank, exons, introns, 3' flank), each knowing its chromosome coordinate,
which exon, intron or flank it lies in, its cDNA position and HGVS c. number, and
the amino acid it carries. Rows are rendered from those columns, so the tracks
always line up and share one colour code:

    genomic   121 tggtaagTGGCAGAGAAGGTGCTGGTAACAGGT   exons upper case, introns lower case
    cDNA      c.1        ATGGCAGAGAAGGTGCTGGTAACAGGT   same exon colours, CDS/UTR, start/stop
    protein     1         M--A--E--K--V--L--V--T--G   residue under the codon's middle base

Sequence comes from the installed genome and exon/CDS structure from the gene
annotation, so the view works offline. When a RefSeq accession is linked, its own
sequence is compared base by base and differences from the genome are marked.
Variants are added from Ensembl.
"""
from __future__ import annotations

import html
import subprocess

from . import annotation, blastconf, config, seqviews, transcripts
from .engines import revcomp

# The standard genetic code, in TCAG order.
CODONS = {a + b + c: aa for (a, b, c), aa in zip(
    [(a, b, c) for a in "TCAG" for b in "TCAG" for c in "TCAG"],
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
THREE_LETTER = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "E": "Glu",
                "Q": "Gln", "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
                "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
                "Y": "Tyr", "V": "Val", "*": "Ter", "X": "Xaa"}

WIDTHS = [(w, f"{w} bases per line") for w in (30, 45, 60, 75, 90, 120, 150)]
FLANKS = [(f, "None" if f == 0 else f"{f:,} bp") for f in (0, 100, 250, 500, 1000, 2000, 5000)]
INTRON_ENDS = 100                    # bases kept at each end of an intron when shortened
MAX_FULL_REGION = 1_000_000          # longer regions shorten their introns automatically
VARIANT_REGION_LIMIT = 50_000        # Ensembl overlap requests are made per region below this

TABS = [("genomic", "Genomic DNA"), ("cdna", "cDNA"), ("aligned", "Aligned")]
VIEW_OPTIONS = [
    ("tab", "View", TABS, "genomic"),
    ("layout", "Sequence shown", [("genomic", "Entire genomic sequence (flanks, exons, introns)"),
                                  ("spliced", "Exons only (cDNA)")], "genomic"),
    ("flank5", "5′ flanking sequence", FLANKS, "1000"),
    ("flank3", "3′ flanking sequence", FLANKS, "1000"),
    ("intron_display", "Introns", [("full", "Show in full"),
                                   ("ends", f"Shorten to {INTRON_ENDS} bp at each end")], "full"),
    ("display_width", "Number of bases per row", WIDTHS, "60"),
    ("numbering", "Genomic numbering", [("chromosome", "Chromosome coordinates"),
                                        ("region", "Relative to the displayed sequence")],
     "chromosome"),
    ("cdna_numbering", "cDNA numbering", [("hgvs", "HGVS c. (c.1 = A of the ATG)"),
                                          ("transcript", "Transcript position")], "hgvs"),
    ("protein_display", "Protein translation", [("on", "Show"), ("off", "Hide")], "on"),
    ("snp_display", "Show variants", [("on", "Yes"), ("off", "No")], "off"),
    # cDNA tab (Ensembl's transcript cDNA sequence page)
    ("cdna_snp_display", "Show variants", [("on", "Yes, coloured by consequence"),
                                           ("off", "No")], "on"),
    ("coding_display", "Coding sequence line", [("on", "Show"), ("off", "Hide")], "on"),
    ("exon_display", "Exons", [("on", "Alternate black / blue"),
                               ("off", "Not marked")], "on"),
    ("codon_display", "Codons", [("on", "Shade alternate codons"), ("off", "Not shaded")], "on"),
    ("title_display", "Display pop-up information on mouseover",
     [("yes", "Yes"), ("off", "No")], "yes"),
]


class ViewError(RuntimeError):
    pass


def options(params: dict) -> dict:
    out = {}
    for name, _label, values, default in VIEW_OPTIONS:
        value = str(params.get(name, default))
        out[name] = value if value in {str(v) for v, _ in values} else default
    return out


# Options that apply to each tab; the rest are hidden from its Configure panel.
CDNA_ONLY = {"cdna_snp_display", "coding_display", "exon_display", "codon_display"}
TAB_FIELDS = {
    "genomic": {"flank5", "flank3", "intron_display", "display_width", "numbering",
                "snp_display", "title_display"},
    "cdna": {"display_width", "cdna_snp_display", "coding_display", "protein_display",
             "exon_display", "codon_display", "title_display"},
    "aligned": {n for n, *_ in VIEW_OPTIONS} - {"tab"} - CDNA_ONLY,
}


def option_form(current: dict) -> list[dict]:
    shown = TAB_FIELDS.get(current.get("tab", "aligned"), TAB_FIELDS["aligned"])
    return [{"name": n, "label": label, "value": current[n],
             "values": [{"value": str(v), "caption": c} for v, c in values]}
            for n, label, values, _d in VIEW_OPTIONS if n in shown]


def effective_layout(opts: dict) -> str:
    """Genomic DNA always shows the whole region; cDNA is always spliced."""
    return {"genomic": "genomic", "cdna": "spliced"}.get(opts.get("tab"), opts["layout"])


def translate(seq: str) -> str:
    return "".join(CODONS.get(seq[i:i + 3].upper(), "X")
                   for i in range(0, len(seq) - len(seq) % 3, 3))


# ---------------------------------------------------------------------------
# Reference cDNA from RefSeq
# ---------------------------------------------------------------------------
def refseq_sequence(manifest: dict, accession: str) -> tuple[str, str] | None:
    """(versioned accession, sequence) of a RefSeq transcript from the local database."""
    prefix = (manifest.get("refseq") or {}).get("prefix")
    tool = config.blast_bin("blastdbcmd")
    if not prefix or not tool or not accession:
        return None
    base, version = transcripts.split_version(accession)
    if not version:
        info = transcripts.refseq_info(manifest, base) or {}
        version = info.get("version") or ""
    entry = f"{base}.{version}" if version else base
    try:
        proc = subprocess.run([tool, "-db", prefix, "-entry", entry, "-outfmt", "%s"],
                              capture_output=True, text=True, timeout=60,
                              creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
    except (OSError, subprocess.TimeoutExpired):
        return None
    seq = "".join(proc.stdout.split()).upper()
    return (entry, seq) if proc.returncode == 0 and seq else None


# ---------------------------------------------------------------------------
# Column model
# ---------------------------------------------------------------------------
def _layout(tx: dict, index, opts: dict) -> tuple[list[dict], dict]:
    strand = tx["strand"]
    exons = sorted(tx["exons"], key=lambda e: e["start"] * strand)
    if not exons:
        raise ViewError("This transcript has no exons in the annotation.")
    n_exons = len(exons)
    ex_lo = min(e["start"] for e in exons)
    ex_hi = max(e["end"] for e in exons)
    genomic = effective_layout(opts) == "genomic"
    f5 = int(opts["flank5"]) if genomic else 0
    f3 = int(opts["flank3"]) if genomic else 0
    lo, hi = (ex_lo - f5, ex_hi + f3) if strand > 0 else (ex_lo - f3, ex_hi + f5)
    contig_len = index.index[index._resolve(tx["chrom"])][0]           # noqa: SLF001
    lo, hi = max(1, lo), min(contig_len, hi)
    raw = index.fetch(tx["chrom"], lo, hi).upper()
    if len(raw) < hi - lo + 1:
        raise ViewError(f"The genome has no sequence at {tx['chrom']}:{lo}-{hi}.")

    notes = []
    shorten = opts["intron_display"] == "ends"
    if genomic and not shorten and hi - lo + 1 > MAX_FULL_REGION:
        shorten = True
        notes.append(f"This gene spans {hi - lo + 1:,} bp, so introns are shortened to "
                     f"{INTRON_ENDS} bp at each end. Use exons only, or download the entire "
                     "genomic sequence as FASTA.")

    # Transcript-order coordinate t = pos * strand increases 5' -> 3'.
    ranges = [(min(e["start"] * strand, e["end"] * strand),
               max(e["start"] * strand, e["end"] * strand)) for e in exons]
    cds_positions = set()
    for c in tx["cds"]:
        cds_positions.update(range(c["start"], c["end"] + 1))
    phase = tx["cds"][0]["phase"] if tx["cds"] else 0

    if genomic:
        t_positions = range(min(lo * strand, hi * strand), max(lo * strand, hi * strand) + 1)
    else:
        t_positions = [t for a, b in ranges for t in range(a, b + 1)]

    columns: list[dict] = []
    ei = 0
    tx_pos = 0
    cds_pos = 0
    last_intron_skip = None
    previous = None
    for t in t_positions:
        pos = t * strand
        while ei < n_exons and t > ranges[ei][1]:
            ei += 1
        if ei < n_exons and ranges[ei][0] <= t <= ranges[ei][1]:
            region, number = "exon", ei + 1
        elif ei == 0:
            region, number = "flank5", None
        elif ei == n_exons:
            region, number = "flank3", None
        else:
            region, number = "intron", ei
            start_t, end_t = ranges[ei - 1][1] + 1, ranges[ei][0] - 1
            length = end_t - start_t + 1
            if shorten and length > 2 * INTRON_ENDS + 20 \
                    and start_t + INTRON_ENDS <= t <= end_t - INTRON_ENDS:
                if last_intron_skip != number:
                    last_intron_skip = number
                    columns.append({"skip": length - 2 * INTRON_ENDS, "intron": number})
                continue
        base = raw[pos - lo]
        if strand < 0:
            base = revcomp(base)
        col = {"pos": pos, "base": base, "region": region, "number": number,
               "tx": None, "cds": None, "aa": None, "aa_pos": None}
        if (region, number) != previous:
            previous = (region, number)
            col["feature"] = {"exon": f"exon {number}", "intron": f"intron {number}",
                              "flank5": "5′ flank", "flank3": "3′ flank"}[region]
        if region == "exon":
            tx_pos += 1
            col["tx"] = tx_pos
            if pos in cds_positions:
                cds_pos += 1
                col["cds"] = cds_pos
        elif region == "intron":
            start_t, end_t = ranges[number - 1][1] + 1, ranges[number][0] - 1
            if t - start_t < 2:
                col["splice"] = "donor"
            elif end_t - t < 2:
                col["splice"] = "acceptor"
        columns.append(col)

    exonic = [c for c in columns if "skip" not in c and c["tx"]]
    coding = [c for c in exonic if c["cds"]]
    meta = {"n_exons": n_exons, "region": {"chrom": tx["chrom"], "start": lo, "end": hi,
                                           "strand": strand},
            "notes": notes, "cdna_length": len(exonic), "cds_length": len(coding),
            "protein": "", "coding": bool(coding), "coding_columns": [], "phase": 0,
            "exonic_columns": exonic}
    if coding:
        cds_seq = "".join(c["base"] for c in coding)[phase:]
        protein = translate(cds_seq)
        # Ensembl's CDS excludes the stop codon; take it from the next exonic bases.
        after = exonic[exonic.index(coding[-1]) + 1:][:3]
        stop = "".join(c["base"] for c in after)
        if len(after) == 3 and CODONS.get(stop) == "*":
            for i, c in enumerate(after):
                c["stop"] = True
                c["cds"] = len(coding) + 1 + i
            protein += "*"
            coding = coding + after
        for c in coding[:3]:
            if not phase and CODONS.get("".join(x["base"] for x in coding[:3])) == "M":
                c["start"] = True
        # Residues sit under the middle base of their codon.
        for k, c in enumerate(coding):
            offset = k - (phase or 0)
            if offset >= 0 and offset % 3 == 1 and offset // 3 < len(protein):
                c["aa"] = protein[offset // 3]
                c["aa_pos"] = offset // 3 + 1
        meta["protein"] = protein.rstrip("*")
        meta["cds_length"] = len(coding)
        meta["coding_columns"] = coding
        meta["phase"] = phase or 0
        first_tx, last_tx = coding[0]["tx"], coding[-1]["tx"]
        for c in exonic:
            if c["tx"] < first_tx:
                c["c"] = f"-{first_tx - c['tx']}"
            elif c["tx"] > last_tx:
                c["c"] = f"*{c['tx'] - last_tx}"
            else:
                c["c"] = str(c["tx"] - first_tx + 1)
    else:
        for c in exonic:
            c["c"] = str(c["tx"])                     # n. numbering for non-coding RNA
    meta["cdna"] = "".join(c["base"] for c in exonic)
    meta["genomic"] = "".join(c["base"] for c in columns if "skip" not in c)
    meta["complete_genomic"] = not any("skip" in c for c in columns)
    return columns, meta


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------
def _num(value, width: int) -> str:
    text = "" if value is None else (f"{value:,}" if isinstance(value, int) else str(value))
    return f'<span class="sv-ln">{text:>{width}}</span>'


def _region_label(col: dict, n_exons: int) -> str:
    region, number = col["region"], col["number"]
    if region == "exon":
        kind = ("start codon" if col.get("start") else "stop codon" if col.get("stop")
                else "coding" if col["cds"] else "UTR" if col.get("c", "").startswith(("-", "*"))
                else "non-coding")
        return f"Exon {number} of {n_exons} · {kind}"
    if region == "intron":
        if col.get("splice"):
            return f"Intron {number} · splice {col['splice']} site"
        return f"Intron {number} of {n_exons - 1}"
    return "5′ flanking sequence (upstream)" if region == "flank5" else \
        "3′ flanking sequence (downstream)"


def _classes(col: dict, row: str) -> str:
    """Shared colour code for the genomic and cDNA rows."""
    region = col["region"]
    out = []
    if region == "exon":
        out.append("tv-e1" if col["number"] % 2 else "tv-e2")
        if col.get("start"):
            out.append("tv-start")
        elif col.get("stop"):
            out.append("tv-stop")
        elif col["cds"]:
            out.append("tv-cds")
        else:
            out.append("tv-utr")
        if col.get("mismatch"):
            out.append("tv-mm")
    elif row == "genomic":
        out.append("tv-intron" if region == "intron" else "tv-flank")
        if col.get("splice"):
            out.append("tv-splice")
    if col.get("variant"):
        out.append("sv-var")
    return " ".join(out)


def _row(label: str, start, cells, end, num_w: int, label_w: int) -> str:
    body = []
    i = 0
    while i < len(cells):
        text, cls, title = cells[i]
        j = i + 1
        while j < len(cells) and cells[j][1] == cls and cells[j][2] == title:
            text += cells[j][0]
            j += 1
        escaped = html.escape(text).replace(" ", "&nbsp;")
        if cls or title:
            attrs = f' class="{cls}"' if cls else ""
            if title:
                attrs += f' title="{html.escape(title, quote=True)}"'
            body.append(f"<span{attrs}>{escaped}</span>")
        else:
            body.append(escaped)
        i = j
    return (f'<span class="sv-lab">{label:<{label_w}}</span>{_num(start, num_w)} '
            f'{"".join(body)} {_num(end, 0)}')


def _render_block(block: list[dict], opts: dict, meta: dict, num_w: int, label_w: int,
                  region_start: int, mode: str = "aligned") -> str:
    titles = opts["title_display"] == "yes"
    n_exons = meta["n_exons"]
    genomic_cells, cdna_cells, protein_cells, variant_cells = [], [], [], []
    for col in block:
        label = _region_label(col, n_exons) if titles else None
        v = col.get("variant")
        vtitle = f"{v['id']} {v['alleles']} {v['consequence']}".strip() if (v and titles) else None
        base = col["base"] if col["region"] == "exon" else col["base"].lower()
        genomic_cells.append((base, _classes(col, "genomic"), vtitle or label))
        if col["region"] == "exon":
            ctitle = None
            if titles:
                ctitle = f"{label} · {'c.' if meta['coding'] else 'n.'}{col['c']}"
                if col.get("mismatch"):
                    ctitle += f" · RefSeq {col['cdna_base']}, genome {col['base']}"
            cdna_cells.append((col.get("cdna_base", col["base"]), _classes(col, "cdna"),
                               vtitle or (label if not col.get("mismatch") else ctitle)))
        else:
            cdna_cells.append((" ", "", None))
        variant_cells.append((_iupac(v, col["base"]) if v else " ", "sv-var" if v else "",
                              vtitle))
        if col["aa"]:
            name = THREE_LETTER.get(col["aa"], col["aa"])
            protein_cells.append((col["aa"], "tv-aa" + (" tv-stop-aa" if col["aa"] == "*" else ""),
                                  f"p.{name}{col['aa_pos']}" if titles else None))
        elif col["cds"]:
            protein_cells.append(("-", "tv-aa-link", None))
        else:
            protein_cells.append((" ", "", None))

    exonic = [c for c in block if c["region"] == "exon"]
    residues = [c for c in block if c["aa"]]
    if opts["numbering"] == "region":
        g_from = abs(block[0]["pos"] - region_start) + 1
        g_to = abs(block[-1]["pos"] - region_start) + 1
    else:
        g_from, g_to = block[0]["pos"], block[-1]["pos"]

    def cnum(c):
        if opts["cdna_numbering"] == "transcript":
            return c["tx"]
        return c["c"]

    # A marker at the end of the line names each exon, intron or flank that begins in it.
    features = [c["feature"] for c in block if c.get("feature")]
    marker = (f'  <span class="tv-feature">▸ {html.escape(", ".join(features))}</span>'
              if features else "")

    rows = []
    if any(cell[0] != " " for cell in variant_cells):
        rows.append(_row("", None, variant_cells, None, num_w, label_w))
    rows.append(_row("genomic", g_from, genomic_cells, g_to, num_w, label_w) + marker)
    if mode == "aligned" and exonic:     # the cDNA tab has its own renderer
        rows.append(_row("cDNA", cnum(exonic[0]), cdna_cells, cnum(exonic[-1]), num_w, label_w))
        if opts["protein_display"] == "on" and (residues or any(c["cds"] for c in exonic)):
            rows.append(_row("protein", residues[0]["aa_pos"] if residues else None,
                             protein_cells, residues[-1]["aa_pos"] if residues else None,
                             num_w, label_w))
    return "\n".join(rows)


# ---------------------------------------------------------------------------
# Variant consequences on this transcript (Ensembl/VEP terms and colours)
# ---------------------------------------------------------------------------
# Most severe first, with the colours Ensembl uses for each consequence.
CONSEQUENCES = [
    ("splice_acceptor_variant", "#ff581a", "Splice acceptor"),
    ("splice_donor_variant", "#ff581a", "Splice donor"),
    ("stop_gained", "#ff0000", "Stop gained"),
    ("frameshift_variant", "#9400d3", "Frameshift"),
    ("stop_lost", "#ff0000", "Stop lost"),
    ("start_lost", "#ffd700", "Start lost"),
    ("inframe_insertion", "#ff69b4", "In-frame insertion"),
    ("inframe_deletion", "#ff69b4", "In-frame deletion"),
    ("missense_variant", "#ffd700", "Missense"),
    ("protein_altering_variant", "#ff0080", "Protein altering"),
    ("splice_region_variant", "#ff7f50", "Splice region"),
    ("start_retained_variant", "#76ee00", "Start retained"),
    ("stop_retained_variant", "#76ee00", "Stop retained"),
    ("synonymous_variant", "#76ee00", "Synonymous"),
    ("coding_sequence_variant", "#458b00", "Coding sequence"),
    ("5_prime_UTR_variant", "#7ac5cd", "5′ UTR"),
    ("3_prime_UTR_variant", "#7ac5cd", "3′ UTR"),
    ("non_coding_transcript_exon_variant", "#32cd32", "Non-coding exon"),
]
SEVERITY = {term: i for i, (term, _c, _l) in enumerate(CONSEQUENCES)}
_DARK_BACKGROUNDS = {"#ff581a", "#ff0000", "#9400d3", "#ff69b4", "#ff0080", "#458b00"}
# Variants with these consequences turn their residues red in the protein line.
PROTEIN_CHANGING = {"missense_variant", "stop_gained", "stop_lost", "frameshift_variant"}
_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def consequence_css() -> str:
    """Background colour per consequence, shared by the page and its legend."""
    rules = []
    for term, colour, _label in CONSEQUENCES:
        cls = _csq_class(term)
        if colour in _DARK_BACKGROUNDS:
            rules.append(f".{cls}{{background:{colour};color:#fff}}")
        else:
            # Fixed backgrounds get fixed inks, so they read the same in either theme.
            rules.append(f".{cls}{{background:{colour};color:#333}}"
                         f".{cls}.tvc-exon-alt{{color:#1044ee}}")
    return "\n".join(rules)


def _csq_class(term: str) -> str:
    return "csq-" + term.lower().replace("_", "-")


def _worst(terms) -> str | None:
    terms = [t for t in terms if t in SEVERITY]
    return min(terms, key=SEVERITY.get) if terms else None


def _trim(ref: str, alt: str) -> tuple[int, str, str]:
    """Drop the shared leading/trailing bases VCF-style alleles carry; returns (offset, ref, alt)."""
    offset = 0
    while ref and alt and ref[0] == alt[0]:
        ref, alt, offset = ref[1:], alt[1:], offset + 1
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    return offset, ref, alt


def annotate_variants(meta: dict, variants: list[dict], strand: int) -> None:
    """Mark exonic columns with the variants over them, the way Ensembl's cDNA page does.

    Every variant gets its consequence on this transcript and covers its whole span
    (an insertion covers the bases either side). Each base shows one of its variants: the
    IUPAC code of a substitution or * for anything else, coloured by consequence. Codons
    hit by a missense, stop-gained, stop-lost or frameshift variant are flagged so the
    protein line shows their residue in red.
    """
    exonic = meta["exonic_columns"]
    by_pos = {c["pos"]: c for c in exonic}
    order = {id(c): i for i, c in enumerate(exonic)}
    coding = meta["coding_columns"]
    phase = meta["phase"]
    coding_index = {id(c): i for i, c in enumerate(coding)}
    # Exonic bases within 3 bp of a splice junction count as splice region.
    near_junction = set()
    for i in range(len(exonic) - 1):
        if exonic[i]["number"] != exonic[i + 1]["number"]:
            for j in range(max(0, i - 2), min(len(exonic), i + 4)):
                near_junction.add(id(exonic[j]))

    def codon_of(col):
        k = coding_index.get(id(col))
        if k is None or k < phase:
            return None, []
        ci = (k - phase) // 3
        cols = coding[phase + 3 * ci: phase + 3 * ci + 3]
        return (ci, cols) if len(cols) == 3 else (None, [])

    def base(col):
        return col.get("cdna_base", col["base"])

    for v in variants:
        alleles = [a.upper().replace("-", "") for a in v.get("allele_list") or []]
        valid = len(alleles) >= 2 and all(set(a) <= set("ACGTN") for a in alleles)
        snv = valid and all(len(a) == 1 for a in alleles)
        lo, hi = sorted((v["start"], v["end"]))       # insertions cover both flanking bases
        span = [by_pos[p] for p in range(lo, hi + 1) if p in by_pos]
        if not span:
            continue
        terms = []
        if not valid:                                   # e.g. HGMD_MUTATION, COSMIC_MUTATION
            term = "coding_sequence_variant" if any(c["cds"] for c in span) \
                else _utr_term(span[0], meta)
            if any(id(c) in near_junction for c in span):
                term = "splice_region_variant"
            terms.append(term)
        for alt in (alleles[1:] if valid else ()):
            offset, r, a = _trim(alleles[0], alt)
            start = v["start"] + offset
            if r:
                cols = [by_pos[p] for p in range(start, start + len(r)) if p in by_pos]
            else:                                       # pure insertion: the bases either side
                cols = [by_pos[p] for p in (start - 1, start) if p in by_pos]
            if not cols:
                continue
            coding_cols = [c for c in cols if id(c) in coding_index]
            if not coding_cols:
                term = _utr_term(cols[0], meta)
            elif len(r) != len(a):
                diff = len(a) - len(r)
                if diff % 3:
                    term = "frameshift_variant"
                elif any(c.get("stop") for c in coding_cols):
                    term = "stop_lost"
                elif any(c.get("start") for c in coding_cols):
                    term = "start_lost"
                else:
                    term = "inframe_insertion" if diff > 0 else "inframe_deletion"
            else:
                # Substitution: translate every codon it touches.
                alt_bases = a if strand > 0 else a[::-1].translate(_COMPLEMENT)
                ordered = sorted(cols, key=lambda c: order[id(c)])
                changed = {id(c): b for c, b in zip(ordered, alt_bases)}
                codon_terms, seen = [], set()
                for col in coding_cols:
                    ci, codon_cols = codon_of(col)
                    if ci is None or ci in seen:
                        continue
                    seen.add(ci)
                    ra = CODONS.get("".join(base(c) for c in codon_cols), "X")
                    aa = CODONS.get("".join(changed.get(id(c), base(c)) for c in codon_cols), "X")
                    if ci == 0 and ra == "M" and aa != "M":
                        codon_terms.append("start_lost")
                    elif ra == "*" and aa != "*":
                        codon_terms.append("stop_lost")
                    elif ra != "*" and aa == "*":
                        codon_terms.append("stop_gained")
                    elif ra == aa:
                        codon_terms.append("stop_retained_variant" if ra == "*" else
                                           "start_retained_variant" if ci == 0 and ra == "M"
                                           else "synonymous_variant")
                    else:
                        codon_terms.append("missense_variant")
                term = _worst(codon_terms) or "coding_sequence_variant"
            if term in PROTEIN_CHANGING:
                # Only the codons the trimmed change sits in, not the whole allele span; an
                # insertion changes the protein from the codon after it.
                hit = coding_cols if r else [max(coding_cols, key=lambda c: c["tx"])]
                for col in hit:
                    for c in codon_of(col)[1]:
                        c["aa_changed"] = True
            if any(id(c) in near_junction for c in cols) and \
                    SEVERITY[term] > SEVERITY["splice_region_variant"]:
                term = "splice_region_variant"
            terms.append(term)
        term = _worst(terms)
        if not term:
            continue
        bases = {x if strand > 0 else x.translate(_COMPLEMENT) for x in alleles} if snv else set()
        record = {**v, "term": term, "snv": snv, "bases": bases,
                  "length": v["end"] - v["start"] + 1}  # 0 for an insertion
        for col in span:
            col.setdefault("variants", []).append(record)

    for col in exonic:
        found = col.get("variants")
        if not found:
            continue
        # Ensembl draws the longest variants first and the most severe last among equals,
        # so a base shows its shortest, then most severe, variant: an SNV inside a long
        # deletion keeps its own code and colour, and an insertion beside it wins.
        shown = min(found, key=lambda x: (x["length"], SEVERITY[x["term"]]))
        col["shown"] = shown
        col["csq"] = shown["term"]
        col["code"] = _iupac({"alleles": "/".join(shown["bases"])}, base(col)) \
            if shown["snv"] else "*"


def _utr_term(col: dict, meta: dict) -> str:
    if not meta["coding"]:
        return "non_coding_transcript_exon_variant"
    return "5_prime_UTR_variant" if str(col.get("c", "")).startswith("-") else "3_prime_UTR_variant"


def _cells_html(cells) -> str:
    """(text, class, title, href) cells, merged into runs where nothing differs."""
    out = []
    i = 0
    while i < len(cells):
        text, cls, title, href = cells[i]
        j = i + 1
        while j < len(cells) and cells[j][1:] == (cls, title, href):
            text += cells[j][0]
            j += 1
        escaped = html.escape(text)
        # Newlines in pop-ups are encoded so every sequence line stays one line of HTML.
        attrs = (f' class="{cls}"' if cls else "") + \
            (f' title="{html.escape(title, quote=True).replace(chr(10), "&#10;")}"'
             if title else "")
        if href:
            out.append(f'<a href="{html.escape(href, quote=True)}" target="_blank" '
                       f'rel="noopener"{attrs}>{escaped}</a>')
        elif attrs:
            out.append(f"<span{attrs}>{escaped}</span>")
        else:
            out.append(escaped)
        i = j
    return "".join(out)


CDNA_BASE_TITLES_LIMIT = 20_000      # per-base pop-ups above this length would bloat the page


def _render_cdna_block(block: list[dict], opts: dict, meta: dict, num_w: int,
                       variant_url: str) -> str:
    """One block of Ensembl's transcript cDNA sequence page:

            *R   K KW *K*****K**                          variant codes, linked to Ensembl
        121 TGGCAGAGAAGGTGCTGGTAACAGG ... 180             cDNA, coloured by consequence
          2 TGGCAGAGAAGGTGCTGGTAACAGG ...  61             coding sequence
          1 M--A--E--K--V--L--V--T--G ...  20             protein, changed residues in red
    """
    titles = opts["title_display"] == "yes"
    base_titles = titles and meta["cdna_length"] <= CDNA_BASE_TITLES_LIMIT
    show_variants = opts["cdna_snp_display"] == "on"
    exons_on = opts["exon_display"] == "on"
    codons_on = opts["codon_display"] == "on"
    coding_index = meta["coding_index"]
    phase = meta["phase"]
    variant_cells, cdna_cells, coding_cells, protein_cells = [], [], [], []
    for col in block:
        variants = col.get("variants") if show_variants else None
        classes = []
        k = coding_index.get(id(col))
        if variants:
            classes += [_csq_class(col["csq"]), "tvc-var"]
        elif codons_on and k is not None and k >= phase and ((k - phase) // 3) % 2:
            classes.append("tvc-codon")
        if exons_on and col["number"] % 2 == 0:
            classes.append("tvc-exon-alt")
        if col.get("mismatch"):
            classes.append("tvc-mm")
        title = None
        if base_titles or (titles and (variants or col.get("mismatch"))):
            where = f"cDNA {col['tx']}"
            if meta["coding"]:
                where += f" · c.{col['c']}"
            lines = [f"{where} · exon {col['number']} of {meta['n_exons']}"]
            if col.get("mismatch"):
                lines.append(f"RefSeq {col['cdna_base']}, genome {col['base']}")
            ranked = sorted(variants or [], key=lambda x: (x["length"], SEVERITY[x["term"]]))
            for x in ranked[:6]:
                lines.append(f"{x['id']} {x['alleles']} · {x['term'].replace('_', ' ')}")
            if variants and len(variants) > 6:
                lines.append(f"…and {len(variants) - 6} more")
            title = "\n".join(lines)
        base = col.get("cdna_base", col["base"])
        cdna_cells.append((base, " ".join(classes), title, None))
        if variants:
            ident = col["shown"].get("id") or ""
            href = variant_url + ident if ident else None
            vtitle = "\n".join(title.split("\n")[1:]) if title else None
            variant_cells.append((col["code"], "tvc-code", vtitle, href))
        else:
            variant_cells.append((" ", "", None, None))
        coding_cells.append((base if col["cds"] else " ", "", None, None))
        if col["cds"]:
            cls = "tvc-pep-hit" if show_variants and col.get("aa_changed") else "tvc-pep"
            ptitle = f"p.{THREE_LETTER.get(col['aa'], col['aa'])}{col['aa_pos']}" \
                if titles and col["aa"] else None
            protein_cells.append((col["aa"] or "-", cls, ptitle, None))
        else:
            protein_cells.append((" ", "", None, None))

    width = int(opts["display_width"])

    def line(start, cells, end):
        left = "" if start is None else str(start)
        right = "" if end is None else str(end)
        pad = " " * (width - len(cells))
        return (f'<span class="tvc-num">{left:>{num_w}}</span> {_cells_html(cells)}{pad}  '
                f'<span class="tvc-num">{right:>{num_w}}</span>').rstrip()

    rows = []
    if show_variants and any(cell[0] != " " for cell in variant_cells):
        rows.append(line(None, variant_cells, None))
    rows.append(line(block[0]["tx"], cdna_cells, block[-1]["tx"]))
    coding = [c for c in block if c["cds"]]
    if coding and opts["coding_display"] == "on":
        rows.append(line(coding[0]["cds"], coding_cells, coding[-1]["cds"]))
    residues = [c for c in block if c["aa"]]
    if coding and opts["protein_display"] == "on":
        rows.append(line(residues[0]["aa_pos"] if residues else None, protein_cells,
                         residues[-1]["aa_pos"] if residues else None))
    return "\n".join(rows)


def build(manifest: dict, transcript_id: str, opts: dict, compare_with: str | None = None
          ) -> dict:
    tx = annotation.transcript(manifest, transcript_id)
    if not tx:
        raise ViewError(f"{transcript_id} is not in the installed gene annotation.")
    index = seqviews.fasta_index(manifest)
    if index is None:
        raise ViewError("The genome sequence for this species is not installed.")

    columns, meta = _layout(tx, index, opts)
    exonic = [c for c in columns if "skip" not in c and c["region"] == "exon"]
    notes = list(meta["notes"])

    # Compare the RefSeq transcript's own sequence with the genome.
    ids = transcripts.annotate(manifest, tx["id"])
    accession = compare_with or ids.get("mane_refseq") or ids.get("refseq_mrna")
    reference = None
    mismatches = 0
    # The Genomic DNA tab has no cDNA row, so there is nothing to compare there.
    if accession and not accession.upper().startswith("ENS") and opts.get("tab") != "genomic":
        reference = refseq_sequence(manifest, accession)
    if reference:
        ref_acc, ref_seq = reference
        if len(ref_seq) == len(exonic):
            for col, ref_base in zip(exonic, ref_seq):
                if ref_base != col["base"]:
                    col["mismatch"] = True
                    col["cdna_base"] = ref_base
                    mismatches += 1
            notes.append(f"cDNA row is {ref_acc}; it matches the genome at every base."
                         if not mismatches else
                         f"cDNA row is {ref_acc}; {mismatches} base(s) differ from the "
                         "reference genome and are marked in red.")
        else:
            notes.append(f"{ref_acc} is {len(ref_seq):,} bp but the genome-based transcript is "
                         f"{len(exonic):,} bp (an indel or a different UTR), so the cDNA row "
                         f"shows the {tx['id']} sequence from the genome.")

    mode = opts.get("tab") or "aligned"
    pending: list[dict] = []
    show_variants = opts["cdna_snp_display" if mode == "cdna" else "snp_display"] == "on"
    variants: list[dict] = []
    if show_variants:
        region = meta["region"]
        if mode != "cdna" and region["end"] - region["start"] + 1 <= VARIANT_REGION_LIMIT:
            spans = [(region["start"], region["end"])]
        elif mode == "cdna" and max(e["end"] for e in tx["exons"]) - \
                min(e["start"] for e in tx["exons"]) + 1 <= VARIANT_REGION_LIMIT:
            spans = [(min(e["start"] for e in tx["exons"]), max(e["end"] for e in tx["exons"]))]
        else:
            spans = [(e["start"] - 20, e["end"] + 20) for e in tx["exons"]]
            if mode != "cdna":
                notes.append("This region is too long to fetch every variant, so variants are "
                             "shown for exons and their splice regions.")
        seen = set()
        # Ensembl's cDNA page shows somatic mutations (COSMIC) alongside germline variants.
        kinds = ("variants", "somatic") if mode == "cdna" else ("variants",)
        for s, e in spans:
            for kind in kinds:
                for v in seqviews.lookup(kind, manifest, tx["chrom"], s, e, pending):
                    key = (v["id"], v["start"], v["end"])
                    if key not in seen:
                        seen.add(key)
                        variants.append(v)
    if mode == "cdna":
        return _build_cdna(manifest, tx, opts, meta, variants, notes, pending, reference,
                           mismatches)
    by_pos = {}
    for v in variants:
        for p in range(v["start"], v["end"] + 1):
            by_pos.setdefault(p, v)
    for col in columns:
        if "skip" not in col and col["pos"] in by_pos:
            col["variant"] = by_pos[col["pos"]]

    width = int(opts["display_width"])
    region_start = meta["region"]["start"] if tx["strand"] > 0 else meta["region"]["end"]
    num_w = max(len(f"{meta['region']['end']:,}"), 7) + 1
    label_w = 9
    parts = []
    segment: list[dict] = []

    def flush():
        for lo in range(0, len(segment), width):
            parts.append(_render_block(segment[lo:lo + width], opts, meta, num_w, label_w,
                                       region_start, mode))
        segment.clear()

    for col in columns:
        if "skip" in col:
            flush()
            parts.append(f'<span class="tv-skip">{" " * (label_w + num_w + 1)}··· '
                         f'{col["skip"]:,} bp of intron {col["intron"]} not shown ···</span>')
        else:
            segment.append(col)
    flush()

    legend = [("tv-e1 tv-cds", "Coding exon"), ("tv-e2 tv-cds", "Next exon"),
              ("tv-e1 tv-utr", "UTR"), ("tv-e1 tv-start", "Start codon"),
              ("tv-e1 tv-stop", "Stop codon")]
    if effective_layout(opts) == "genomic":
        legend += [("tv-intron", "Intron (lower case)"), ("tv-intron tv-splice", "Splice site"),
                   ("tv-flank", "Flanking sequence")]
    if mismatches and mode != "genomic":           # differences live in the cDNA row
        legend.append(("tv-e1 tv-mm", "RefSeq differs from genome"))
    has_variants = any(c.get("variant") for c in columns if "skip" not in c)
    if has_variants:
        legend.append(("sv-var", "Variant"))

    return {
        # One track per line reads as continuous sequence; stacked tracks need a gap.
        "html": ("\n" if mode == "genomic" and not has_variants else "\n\n").join(parts),
        "legend": legend, "pending": pending, "notes": notes,
        "protein": meta["protein"], "cdna_length": meta["cdna_length"],
        "cds_length": meta["cds_length"], "exons": meta["n_exons"], "transcript": tx,
        "region": meta["region"], "cdna": meta["cdna"], "genomic": meta["genomic"],
        "complete_genomic": meta["complete_genomic"], "reference": reference[0] if reference else None,
        "mismatches": mismatches, "mode": mode, "layout": effective_layout(opts),
    }


def _build_cdna(manifest: dict, tx: dict, opts: dict, meta: dict, variants: list[dict],
                notes: list[str], pending: list[dict], reference, mismatches: int) -> dict:
    """The cDNA tab, laid out like Ensembl's transcript cDNA sequence page."""
    from . import species
    exonic = meta["exonic_columns"]
    annotate_variants(meta, variants, tx["strand"])
    variant_url = (f"{species.ensembl_site(manifest)}/{species.web_name(manifest['name'])}"
                   "/Variation/Explore?v=")
    meta["coding_index"] = {id(c): i for i, c in enumerate(meta["coding_columns"])}
    width = int(opts["display_width"])
    num_w = len(str(meta["cdna_length"]))
    parts = [_render_cdna_block(exonic[i:i + width], opts, meta, num_w, variant_url)
             for i in range(0, len(exonic), width)]

    legend = []
    if opts["exon_display"] == "on" and meta["n_exons"] > 1:
        legend += [("tvc-exon-key", "Odd exons"), ("tvc-exon-key tvc-exon-alt", "Even exons")]
    if opts["codon_display"] == "on" and meta["coding"]:
        legend.append(("tvc-codon", "Alternate codons"))
    if opts["cdna_snp_display"] == "on":
        present = {c["csq"] for c in exonic if c.get("csq")}
        legend += [(_csq_class(term), label) for term, _colour, label in CONSEQUENCES
                   if term in present]
        if any(c.get("aa_changed") for c in exonic) and opts["protein_display"] == "on":
            legend.append(("tvc-pep-hit", "Residue changed by a variant"))
    if mismatches:
        legend.append(("tvc-mm", "RefSeq differs from genome"))

    n_variants = len({(x["id"], x["start"]) for c in exonic for x in c.get("variants") or ()})
    species_label = (manifest.get("display_name")
                     or manifest["name"].replace("_", " ").capitalize())
    return {
        "html": "\n\n".join(parts),
        "legend": legend, "pending": pending, "notes": notes,
        "protein": meta["protein"], "cdna_length": meta["cdna_length"],
        "cds_length": meta["cds_length"], "exons": meta["n_exons"], "transcript": tx,
        "region": meta["region"], "cdna": "".join(c.get("cdna_base", c["base"]) for c in exonic),
        "genomic": meta["genomic"], "complete_genomic": meta["complete_genomic"],
        "reference": reference[0] if reference else None, "mismatches": mismatches,
        "mode": "cdna", "layout": "spliced",
        "title": f"{tx.get('gene_name') or tx['id']} {species_label} cDNA",
        "variant_count": n_variants,
        "variants_shown": opts["cdna_snp_display"] == "on",
    }


def _iupac(variant: dict, reference: str) -> str:
    """Ambiguity code for the alleles at a variant position, as sequence viewers show."""
    codes = {frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("GC"): "S",
             frozenset("AT"): "W", frozenset("GT"): "K", frozenset("AC"): "M",
             frozenset("CGT"): "B", frozenset("AGT"): "D", frozenset("ACT"): "H",
             frozenset("ACG"): "V", frozenset("ACGT"): "N"}
    alleles = {a for a in (variant.get("alleles") or "").upper().split("/") if len(a) == 1}
    alleles.add(reference.upper())
    alleles = {a for a in alleles if a in "ACGT"}
    if len(alleles) < 2:
        return "*"                      # an indel or other change with no single-base code
    return codes.get(frozenset(alleles), "N")


def summary(manifest: dict, tx: dict, view: dict) -> list[tuple[str, str]]:
    ids = transcripts.annotate(manifest, tx["id"])
    rows = [("Transcript", tx["id"] + (f".{tx.get('version')}" if tx.get("version") else "")),
            ("Gene", f"{tx.get('gene_name') or ''} ({tx.get('gene_id') or ''})"),
            ("Biotype", (tx.get("biotype") or "").replace("_", " "))]
    if ids.get("mane"):
        rows.append(("MANE", f"{ids['mane']}: {ids.get('mane_refseq', '')} = "
                             f"{ids.get('mane_ensembl', '')}"))
    if ids.get("refseq_mrna"):
        rows.append(("RefSeq", f"{ids['refseq_mrna']}"
                               + (f" / {ids['refseq_peptide']}" if ids.get("refseq_peptide") else "")))
    exon_lo = min(e["start"] for e in tx["exons"])
    exon_hi = max(e["end"] for e in tx["exons"])
    rows += [("Location", f"{tx['chrom']}:{exon_lo:,}-{exon_hi:,} ({'+' if tx['strand'] > 0 else '-'})"),
             ("Gene length", f"{exon_hi - exon_lo + 1:,} bp"),
             ("Exons", str(view["exons"])),
             ("cDNA length", f"{view['cdna_length']:,} bp")]
    if view["protein"]:
        rows.append(("Protein", f"{len(view['protein']):,} aa"))
    return rows


def fasta(name: str, seq: str) -> str:
    return blastconf.fasta(name, seq)


def genomic_sequence(manifest: dict, transcript_id: str, flank5: int, flank3: int) -> dict:
    """The complete gene region with flanks, 5'->3' on the transcript strand (no shortening)."""
    tx = annotation.transcript(manifest, transcript_id)
    if not tx:
        raise ViewError(f"{transcript_id} is not in the installed gene annotation.")
    index = seqviews.fasta_index(manifest)
    if index is None:
        raise ViewError("The genome sequence for this species is not installed.")
    ex_lo = min(e["start"] for e in tx["exons"])
    ex_hi = max(e["end"] for e in tx["exons"])
    strand = tx["strand"]
    lo, hi = (ex_lo - flank5, ex_hi + flank3) if strand > 0 else (ex_lo - flank3, ex_hi + flank5)
    contig_len = index.index[index._resolve(tx["chrom"])][0]           # noqa: SLF001
    lo, hi = max(1, lo), min(contig_len, hi)
    if hi - lo + 1 > 5_000_000:
        raise ViewError("The region is longer than 5 Mb; reduce the flanks.")
    seq = index.fetch(tx["chrom"], lo, hi).upper()
    if strand < 0:
        seq = revcomp(seq)
    return {"transcript": tx, "chrom": tx["chrom"], "start": lo, "end": hi, "strand": strand,
            "sequence": seq}
