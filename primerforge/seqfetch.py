"""Fetch a query sequence by identifier, like Ensembl's BLAST form.

Ensembl stable IDs come from the Ensembl REST API (ENSG -> genomic, ENST -> cDNA,
ENSP -> protein); RefSeq accessions from NCBI E-utilities; UniProt accessions from
UniProt; anything else is tried as an ENA/EMBL accession.
"""
from __future__ import annotations

import re
import urllib.parse

from . import blastconf, config, net, store

ENSEMBL_RE = re.compile(r"^ENS[A-Z]*([GTPE])\d{6,}(\.\d+)?$", re.I)
REFSEQ_NUC_RE = re.compile(r"^(NM|NR|XM|XR|NC|NG|NT|NW)_\d+(\.\d+)?$", re.I)
REFSEQ_PROT_RE = re.compile(r"^(NP|XP|YP|WP)_\d+(\.\d+)?$", re.I)
UNIPROT_RE = re.compile(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})"
                        r"(-\d+)?$", re.I)
ID_RE = re.compile(r"^[A-Za-z][A-Za-z0-9._\-]{3,40}$")


class FetchError(ValueError):
    pass


def looks_like_id(text: str) -> bool:
    t = text.strip()
    return bool(ID_RE.match(t) and re.search(r"\d", t))


def fetch(identifier: str) -> list[dict]:
    ident = identifier.strip()
    if not looks_like_id(ident):
        raise FetchError(f"'{ident}' does not look like a sequence identifier.")
    cache_key = f"seqfetch:{ident.upper()}"
    cached = store.cache_get(cache_key, max_age=30 * 86400)
    if cached:
        return cached

    m = ENSEMBL_RE.match(ident)
    if m:
        kind = {"G": "genomic", "T": "cdna", "P": "protein", "E": "genomic"}[m.group(1).upper()]
        url = (f"{config.ENSEMBL_REST}/sequence/id/{urllib.parse.quote(ident.split('.')[0])}"
               f"?type={kind};content-type=text/x-fasta")
        source = "Ensembl"
    elif REFSEQ_NUC_RE.match(ident) or REFSEQ_PROT_RE.match(ident):
        db = "protein" if REFSEQ_PROT_RE.match(ident) else "nuccore"
        url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
               f"db={db}&id={urllib.parse.quote(ident)}&rettype=fasta&retmode=text")
        source = "NCBI RefSeq"
    elif UNIPROT_RE.match(ident):
        url = f"https://rest.uniprot.org/uniprotkb/{urllib.parse.quote(ident.upper())}.fasta"
        source = "UniProt"
    else:
        url = f"https://www.ebi.ac.uk/ena/browser/api/fasta/{urllib.parse.quote(ident)}"
        source = "ENA"

    try:
        raw = net.fetch(url, accept="text/plain, text/x-fasta, */*", attempts=2, timeout=60)
    except net.HttpError as exc:
        raise FetchError(f"Could not find any sequence for {ident} in {source} ({exc}).") from exc
    text = raw.decode("utf-8", errors="replace")
    if not text.lstrip().startswith(">"):
        raise FetchError(f"Could not find any sequence for {ident} in {source}.")
    if len(text) > blastconf.MAX_SEQUENCE_LENGTH * 1.2:
        raise FetchError(f"{ident} is longer than the {blastconf.MAX_SEQUENCE_LENGTH:,}-character "
                         "limit for a query sequence.")
    parsed = blastconf.parse_sequences(text)
    if not parsed["sequences"]:
        raise FetchError(f"{source} returned no usable sequence for {ident}.")
    out = [{"description": s["description"], "sequence": s["sequence"], "type": s["type"],
            "source": source} for s in parsed["sequences"]]
    store.cache_put(cache_key, out)
    return out
