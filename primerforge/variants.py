"""Turn what a user types into a precise genomic locus.

Accepts HGVS (c./g./n./m. on RefSeq, Ensembl or LRG references), plain VCF-style
coordinates and rsIDs. Resolution goes through the Ensembl VEP REST endpoint,
which returns coordinates *and* gene/exon/protein annotation in a single call;
the Variant Recoder is used as a fallback for inputs VEP rejects.
"""
from __future__ import annotations

import json
import re
import urllib.parse
from dataclasses import dataclass, field, asdict

from . import net, store
from .config import ENSEMBL_REST, ENSEMBL_REST_GRCH37


class VariantError(ValueError):
    """Raised when an input cannot be parsed or resolved to a locus."""


# --------------------------------------------------------------------------
# Input shapes
# --------------------------------------------------------------------------
HGVS_RE = re.compile(
    r"^(?P<ref>[A-Za-z0-9_.\-]+)\s*:\s*(?P<type>[cgnmr])\s*\.\s*(?P<desc>.+)$", re.I)
VCF_RE = re.compile(
    r"^(?:chr)?(?P<chrom>[0-9]{1,2}|[XYMT]{1,2})[\s:\-_|]+(?P<pos>[0-9,]+)"
    r"[\s:\-_|>/]+(?P<ref>[ACGTacgt]+)[\s:\-_|>/]+(?P<alt>[ACGTacgt]+|-)$")
RSID_RE = re.compile(r"^rs[0-9]+$", re.I)
# Substitution / del / dup / ins / delins / inv, with optional intronic offsets.
DESC_RE = re.compile(
    r"^(?:\*?-?\d+[+\-]?\d*)(?:_\*?-?\d+[+\-]?\d*)?"
    r"(?:[ACGT]+>[ACGT]+|del[ACGT]*|dup[ACGT]*|ins[ACGT]+|delins[ACGT]+|inv|[ACGT]+\[\d+\])$", re.I)

CHROM_ORDER = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass
class Locus:
    """A resolved variant, in 1-based inclusive genomic coordinates."""
    input: str
    kind: str                       # hgvs_c | hgvs_g | vcf | rsid
    assembly: str
    chrom: str
    start: int
    end: int
    ref: str
    alt: str
    vtype: str                      # SNV | deletion | insertion | duplication | delins | inv
    hgvs_g: str | None = None
    gene: str | None = None
    gene_id: str | None = None
    transcript: str | None = None
    mane: str | None = None
    hgvs_c: str | None = None
    hgvs_p: str | None = None
    exon: str | None = None
    intron: str | None = None
    consequence: str | None = None
    strand: int | None = None
    rsid: str | None = None
    transcript_region: dict | None = None   # {start, end, exons: [...]}
    notes: list[str] = field(default_factory=list)

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    def label(self) -> str:
        bits = [self.gene] if self.gene else []
        bits.append(self.hgvs_c or self.hgvs_g or self.input)
        if self.hgvs_p:
            bits.append(self.hgvs_p.split(":")[-1])
        return "  ".join(b for b in bits if b)

    def to_dict(self) -> dict:
        d = asdict(self)
        d["length"] = self.length
        d["label"] = self.label()
        return d


# --------------------------------------------------------------------------
# HTTP helper
# --------------------------------------------------------------------------
def _rest(path: str, assembly: str = "GRCh38", cache: bool = True) -> object:
    base = ENSEMBL_REST if assembly == "GRCh38" else ENSEMBL_REST_GRCH37
    url = f"{base}{path}"
    key = f"rest:{url}"
    if cache:
        hit = store.cache_get(key)
        if hit is not None:
            return hit
    try:
        data = net.fetch_json(url)
    except net.HttpError as exc:
        raise VariantError(str(exc)) from exc
    if cache:
        store.cache_put(key, data)
    return data


# --------------------------------------------------------------------------
# Syntax checking -- fast, offline, and with actionable messages
# --------------------------------------------------------------------------
def classify(raw: str) -> tuple[str, str]:
    """Return (kind, cleaned) or raise VariantError with a helpful message."""
    text = " ".join(raw.split()).strip().strip(",;")
    if not text:
        raise VariantError("Empty input.")
    if RSID_RE.match(text):
        return "rsid", text.lower()
    m = VCF_RE.match(text)
    if m:
        return "vcf", text
    m = HGVS_RE.match(text)
    if m:
        ref, vtype, desc = m.group("ref"), m.group("type").lower(), m.group("desc")
        desc = desc.replace(" ", "")
        if not DESC_RE.match(desc):
            raise VariantError(
                f"'{desc}' is not a recognised HGVS change. Expected something like "
                "215C>G, 1521_1523del, 35dup, 35_36insACGT or 68_70delinsAC.")
        if vtype == "c" and not re.match(r"^(NM_|ENST|LRG_|XM_|NR_)", ref, re.I):
            raise VariantError(
                f"'{ref}' does not look like a transcript. A c. description needs a "
                "transcript reference such as NM_000546.6 or ENST00000269305.9.")
        if vtype == "g" and not re.match(r"^(NC_|ENSG|LRG_|chr)?[0-9XYMT]", ref, re.I):
            raise VariantError(
                f"'{ref}' does not look like a genomic reference for a g. description. "
                "Use NC_000017.11 or a chromosome name such as 17.")
        if vtype == "c" and "." not in ref:
            raise VariantError(
                f"'{ref}' has no version number. Use a versioned transcript such as "
                f"{ref}.1 so the coordinates are unambiguous.")
        return ("hgvs_" + vtype), f"{ref}:{vtype}.{desc}"
    if ":" in text and re.search(r"[CGNMR]\s*\.", text):
        raise VariantError(
            "HGVS variant types are lower case. Use 'c.' for coding, 'g.' for genomic.")
    if re.match(r"^[A-Za-z0-9_.]+$", text):
        raise VariantError(
            f"'{text}' is missing a variant description. Expected REFERENCE:TYPE.CHANGE, "
            "for example NM_000546.6:c.215C>G.")
    raise VariantError(
        f"Could not interpret '{text}'. Supported: HGVS (NM_000546.6:c.215C>G), "
        "genomic HGVS (NC_000017.11:g.7676154G>C), coordinates (17-7676154-G-C) or rsID.")


def _vtype_from_alleles(ref: str, alt: str, raw: str) -> str:
    low = raw.lower()
    if "delins" in low:
        return "delins"
    if "dup" in low:
        return "duplication"
    if "inv" in low:
        return "inv"
    if "ins" in low and "del" not in low:
        return "insertion"
    if "del" in low:
        return "deletion"
    if ref == "-" or not ref:
        return "insertion"
    if alt == "-" or not alt:
        return "deletion"
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    return "delins"


# --------------------------------------------------------------------------
# Resolution
# --------------------------------------------------------------------------
def _vep_path(kind: str, cleaned: str) -> str:
    q = urllib.parse.quote(cleaned, safe="")
    opts = "content-type=application/json&numbers=1&hgvs=1&canonical=1&mane=1&xref_refseq=1"
    if kind == "rsid":
        return f"/vep/human/id/{q}?{opts}"
    if kind == "vcf":
        m = VCF_RE.match(cleaned)
        assert m
        chrom = m.group("chrom").upper().replace("CHR", "")
        pos = int(m.group("pos").replace(",", ""))
        ref, alt = m.group("ref").upper(), m.group("alt").upper()
        end = pos + len(ref) - 1
        region = f"{chrom}:{pos}-{end}"
        return f"/vep/human/region/{region}/{alt}?{opts}"
    return f"/vep/human/hgvs/{q}?{opts}"


def _pick_transcript(consequences: list[dict], prefer: str | None) -> dict | None:
    """Prefer the transcript the user actually named, then MANE, then canonical."""
    if not consequences:
        return None
    if prefer:
        base = prefer.split(".")[0].upper()
        for t in consequences:
            ids = {str(t.get("transcript_id", "")).split(".")[0].upper(),
                   str(t.get("refseq_transcript_ids", "")).upper()}
            refseqs = t.get("refseq_transcript_ids") or []
            if isinstance(refseqs, str):
                refseqs = [refseqs]
            ids |= {r.split(".")[0].upper() for r in refseqs}
            mane = t.get("mane_select") or ""
            ids.add(mane.split(".")[0].upper())
            if base in ids:
                return t
    for key in ("mane_select", "canonical"):
        for t in consequences:
            if t.get(key):
                return t
    coding = [t for t in consequences if t.get("biotype") == "protein_coding"]
    return (coding or consequences)[0]


def resolve(raw: str, assembly: str = "GRCh38") -> Locus:
    """Resolve free-text variant input into a fully annotated Locus."""
    kind, cleaned = classify(raw)
    prefer = cleaned.split(":")[0] if kind.startswith("hgvs") else None

    data = _rest(_vep_path(kind, cleaned), assembly)
    if not isinstance(data, list) or not data:
        raise VariantError(f"Ensembl could not resolve '{cleaned}'.")
    rec = data[0]

    chrom = str(rec.get("seq_region_name") or "")
    start, end = rec.get("start"), rec.get("end")
    if not chrom or start is None:
        raise VariantError(f"Ensembl returned no genomic position for '{cleaned}'.")

    notes: list[str] = []
    allele = str(rec.get("allele_string") or "")
    parts = [p for p in allele.split("/") if p]
    ref = parts[0] if parts else ""
    alts = parts[1:]
    if len(alts) > 1:
        notes.append(
            f"{cleaned} is multi-allelic ({'/'.join(alts)}); primers flank the position "
            f"so they work for every allele. Showing {alts[0]}.")
    alt = alts[0] if alts else ""

    # VEP reports alleles in the orientation of the *input*: a c. description on a
    # minus-strand gene comes back transcript-oriented (top-level strand == -1).
    # Everything downstream works in genomic orientation, so flip when needed.
    if rec.get("strand") == -1:
        ref = revcomp(ref) if ref not in ("", "-") else ref
        alt = revcomp(alt) if alt not in ("", "-") else alt

    vtype = _vtype_from_alleles(ref, alt, cleaned)

    # Ensembl reports pure insertions as start = end + 1; normalise to a span.
    if start > end:
        start, end = end, start

    tc = _pick_transcript(rec.get("transcript_consequences") or [], prefer)
    if prefer and tc:
        named = prefer.split(".")[0].upper()
        got = str(tc.get("transcript_id", "")).split(".")[0].upper()
        mane = str(tc.get("mane_select", "")).split(".")[0].upper()
        refseqs = tc.get("refseq_transcript_ids") or []
        if isinstance(refseqs, str):
            refseqs = [refseqs]
        known = {got, mane} | {r.split(".")[0].upper() for r in refseqs}
        if named not in known:
            notes.append(
                f"Annotation shown is for {tc.get('transcript_id')}; Ensembl did not "
                f"expose {prefer} directly. Genomic coordinates are unaffected.")

    rsid = None
    for cid in ([rec.get("id")] if rec.get("id") else []) + (rec.get("colocated_variants") or []):
        cid = cid.get("id") if isinstance(cid, dict) else cid
        if isinstance(cid, str) and RSID_RE.match(cid):
            rsid = cid
            break

    loc = Locus(
        input=raw.strip(), kind=kind, assembly=assembly, chrom=chrom,
        start=int(start), end=int(end), ref=ref or "-", alt=alt or "-", vtype=vtype,
        hgvs_g=f"{_nc_for(chrom, assembly)}:g."
               f"{_g_desc(start, end, ref, alt, vtype)}" if chrom else None,
        gene=(tc or {}).get("gene_symbol"), gene_id=(tc or {}).get("gene_id"),
        transcript=(tc or {}).get("transcript_id"), mane=(tc or {}).get("mane_select"),
        hgvs_c=(tc or {}).get("hgvsc"), hgvs_p=(tc or {}).get("hgvsp"),
        exon=(tc or {}).get("exon"), intron=(tc or {}).get("intron"),
        consequence=rec.get("most_severe_consequence"),
        strand=(tc or {}).get("strand"), rsid=rsid, notes=notes)

    # In repeats an indel has several equally valid coordinates: VEP/SPDI left-align,
    # while HGVS mandates 3'-most shifting. SNVs are unambiguous, so only pay for the
    # (slower) Variant Recoder call when it can actually change the answer.
    if vtype != "SNV" and assembly == "GRCh38":
        canonical = _canonical_hgvs_g(cleaned, chrom)
        if canonical and canonical != loc.hgvs_g:
            loc.notes.append(
                f"Lies in a repeat: HGVS 3'-shifts this to {canonical}, while the "
                f"left-aligned position is {loc.hgvs_g}. Both describe the same change.")
            loc.hgvs_g = canonical

    if loc.transcript:
        try:
            loc.transcript_region = transcript_structure(loc.transcript, assembly)
        except VariantError:
            pass                     # the gene track is optional decoration
    return loc


def _canonical_hgvs_g(cleaned: str, chrom: str) -> str | None:
    """Authoritative genomic HGVS from the Ensembl Variant Recoder, if available."""
    want = _nc_for(chrom, "GRCh38")
    try:
        data = _rest(
            f"/variant_recoder/human/{urllib.parse.quote(cleaned, safe='')}"
            "?content-type=application/json")
    except VariantError:
        return None
    if not isinstance(data, list) or not data:
        return None
    for entry in data:
        if not isinstance(entry, dict):
            continue
        for allele in entry.values():
            if not isinstance(allele, dict):
                continue
            for g in allele.get("hgvsg") or []:
                if g.startswith(want):
                    return g
    return None


NC_ACCESSIONS = {
    "1": "NC_000001.11", "2": "NC_000002.12", "3": "NC_000003.12", "4": "NC_000004.12",
    "5": "NC_000005.10", "6": "NC_000006.12", "7": "NC_000007.14", "8": "NC_000008.11",
    "9": "NC_000009.12", "10": "NC_000010.11", "11": "NC_000011.10", "12": "NC_000012.12",
    "13": "NC_000013.11", "14": "NC_000014.9", "15": "NC_000015.10", "16": "NC_000016.10",
    "17": "NC_000017.11", "18": "NC_000018.10", "19": "NC_000019.10", "20": "NC_000020.11",
    "21": "NC_000021.9", "22": "NC_000022.11", "X": "NC_000023.11", "Y": "NC_000024.10",
    "MT": "NC_012920.1",
}


def _nc_for(chrom: str, assembly: str) -> str:
    if assembly != "GRCh38":
        return f"chr{chrom}"
    return NC_ACCESSIONS.get(chrom.upper().replace("CHR", ""), f"chr{chrom}")


def _g_desc(start: int, end: int, ref: str, alt: str, vtype: str) -> str:
    span = f"{start}" if start == end else f"{start}_{end}"
    if vtype == "SNV":
        return f"{start}{ref}>{alt}"
    if vtype == "deletion":
        return f"{span}del"
    if vtype == "duplication":
        return f"{span}dup"
    if vtype == "insertion":
        return f"{start}_{end}ins{alt}"
    if vtype == "inv":
        return f"{span}inv"
    return f"{span}delins{alt}"


def transcript_structure(transcript_id: str, assembly: str = "GRCh38") -> dict:
    """Exon coordinates for the gene track drawn under each amplicon."""
    tid = transcript_id.split(".")[0]
    data = _rest(f"/lookup/id/{tid}?expand=1&content-type=application/json", assembly)
    if not isinstance(data, dict):
        raise VariantError("Unexpected transcript payload.")
    exons = sorted(
        ({"start": e["start"], "end": e["end"]} for e in data.get("Exon", [])),
        key=lambda e: e["start"])
    tl = data.get("Translation") or {}
    return {
        "id": data.get("id"), "display_name": data.get("display_name"),
        "start": data.get("start"), "end": data.get("end"),
        "strand": data.get("strand"), "biotype": data.get("biotype"),
        "cds_start": tl.get("start"), "cds_end": tl.get("end"),
        "exons": exons,
    }


def parse_input_block(text: str) -> list[str]:
    """Split a textarea into individual variant strings (newline or comma separated)."""
    items: list[str] = []
    for line in text.replace(";", "\n").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        # Commas separate variants only when they are not inside an HGVS description.
        chunks = [line] if ":" in line and line.count(":") == 1 and "," not in line \
            else [c for c in line.split(",") if c.strip()]
        items.extend(c.strip() for c in chunks if c.strip())
    return items
