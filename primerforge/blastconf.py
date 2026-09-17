"""BLAST/BLAT search configuration, mirroring the Ensembl BLAST/BLAT tool.

The option lists, per-program defaults, sensitivity presets, data sources and
valid query/database/program combinations follow Ensembl's own configuration
(public-plugins: tools/modules/EnsEMBL/Web/BlastConstants.pm and MULTI.ini), so a
search set up here behaves like the same search on ensembl.org. Where Ensembl
offers a value NCBI BLAST+ rejects on the command line (for example a word size
of 2 for BLASTN) it is filtered out here rather than left to fail at run time.
"""
from __future__ import annotations

import re

MAX_SEQUENCE_LENGTH = 200_000
MAX_NUM_SEQUENCES = 30
DNA_THRESHOLD_PERCENT = 85
SEQUENCE_VALID_CHARS = r"A-Za-z\-\.\*\?=~"
SPECIES_SELECTION_LIMIT = 25
MAX_STORED_HITS = 20_000

QUERY_TYPES = {"dna": "DNA", "peptide": "Protein"}
DB_TYPES = {"dna": "DNA database", "peptide": "Protein database"}

# Ensembl's ENSEMBL_BLAST_DATASOURCES, in display order.
SOURCES = {
    "LATESTGP":        {"db_type": "dna", "label": "Genomic sequence"},
    "LATESTGP_MASKED": {"db_type": "dna", "label": "Genomic sequence (hard masked)"},
    "LATESTGP_SOFT":   {"db_type": "dna", "label": "Genomic sequence (soft masked)"},
    "CDNA_ALL":        {"db_type": "dna", "label": "cDNAs (transcripts/splice variants)"},
    "NCRNA":           {"db_type": "dna", "label": "Ensembl Non-coding RNA genes"},
    "REFSEQ_RNA":      {"db_type": "dna", "label": "RefSeq transcripts (NM/NR)"},
    "MANE_SELECT":     {"db_type": "dna", "label": "MANE Select transcripts (NM)"},
    "PEP_ALL":         {"db_type": "peptide", "label": "Proteins (Ensembl)"},
}
GENOMIC_SOURCES = ("LATESTGP", "LATESTGP_MASKED", "LATESTGP_SOFT")
TRANSCRIPT_SOURCES = ("CDNA_ALL", "NCRNA", "REFSEQ_RNA", "MANE_SELECT")
DNA_SOURCES = ("LATESTGP", "LATESTGP_MASKED", "LATESTGP_SOFT", "CDNA_ALL", "NCRNA",
               "REFSEQ_RNA", "MANE_SELECT")

BLAST_TYPES = {"NCBIBLAST": "NCBI Blast", "BLAT": "BLAT"}

SEARCH_TYPES = [
    {"value": "NCBIBLAST_BLASTN", "blast_type": "NCBIBLAST", "method": "BLASTN",
     "program": "blastn", "query_type": "dna", "db_type": "dna",
     "sources": list(DNA_SOURCES), "min_length": 0},
    {"value": "NCBIBLAST_BLASTX", "blast_type": "NCBIBLAST", "method": "BLASTX",
     "program": "blastx", "query_type": "dna", "db_type": "peptide",
     "sources": ["PEP_ALL"], "min_length": 0},
    {"value": "NCBIBLAST_BLASTP", "blast_type": "NCBIBLAST", "method": "BLASTP",
     "program": "blastp", "query_type": "peptide", "db_type": "peptide",
     "sources": ["PEP_ALL"], "min_length": 0},
    {"value": "NCBIBLAST_TBLASTN", "blast_type": "NCBIBLAST", "method": "TBLASTN",
     "program": "tblastn", "query_type": "peptide", "db_type": "dna",
     "sources": list(DNA_SOURCES), "min_length": 0},
    {"value": "NCBIBLAST_TBLASTX", "blast_type": "NCBIBLAST", "method": "TBLASTX",
     "program": "tblastx", "query_type": "dna", "db_type": "dna",
     "sources": list(DNA_SOURCES), "min_length": 0},
    {"value": "BLAT_BLAT", "blast_type": "BLAT", "method": "BLAT",
     "program": "blat", "query_type": "dna", "db_type": "dna",
     "sources": ["LATESTGP"], "min_length": 26},
]
SEARCH_TYPE_BY_VALUE = {s["value"]: s for s in SEARCH_TYPES}
BLAT_VALUE = "BLAT_BLAT"

# Sources a program cannot search even though its db type matches.
RESTRICTIONS = {"NCBIBLAST_TBLASTN": ["NCRNA"], "NCBIBLAST_BLASTP": ["NCRNA"]}

# ---------------------------------------------------------------------------
# Scoring option tables (value "AnB" = gap opening A, extension B)
# ---------------------------------------------------------------------------
MATRIX_GAPS = {
    "BLOSUM45": ["14n2", "13n3", "12n3", "11n3", "10n3", "16n2", "15n2", "13n2", "12n2",
                 "19n1", "18n1", "17n1", "16n1"],
    "BLOSUM50": ["13n2", "13n3", "12n3", "11n3", "10n3", "9n3", "16n2", "15n2", "14n2",
                 "12n2", "19n1", "18n1", "17n1", "16n1", "15n1"],
    "BLOSUM62": ["11n1", "10n2", "9n2", "8n2", "7n2", "6n2", "13n1", "12n1", "10n1", "9n1",
                 "11n2"],
    "BLOSUM80": ["10n1", "25n2", "13n2", "9n2", "8n2", "7n2", "6n2", "11n1", "9n1"],
    "BLOSUM90": ["10n1", "9n2", "8n2", "7n2", "6n2", "11n1", "9n1"],
    "PAM30":    ["9n1", "7n2", "6n2", "5n2", "10n1", "8n1"],
    "PAM70":    ["10n1", "8n2", "7n2", "6n2", "11n1", "9n1"],
    "PAM250":   ["14n2", "15n3", "14n3", "13n3", "12n3", "11n3", "17n2", "16n2", "15n2",
                 "13n2", "21n1", "20n1", "19n1", "18n1", "17n1"],
}

REWARD_LABELS = {
    "1_5": "1,-5", "1_4": "1,-4", "2_7": "2,-7", "1_3": "1,-3", "2_5": "2,-5",
    "1_2": "1,-2", "2_3": "2,-3", "3_4": "3,-4", "4_5": "4,-5", "1_1": "1,-1",
    "3_2": "3,-2", "5_4": "5,-4",
}
# Gap costs BLAST+ accepts for each reward/penalty pair. Ensembl's defaults use
# 5n2 with 1,-2 and 1,-3, which BLAST+ supports, so it is listed first there.
REWARD_GAPS = {
    "1_5": ["3n3"],
    "1_4": ["1n2", "0n2", "2n1", "1n1"],
    "2_7": ["2n4", "0n4", "4n2", "2n2"],
    "1_3": ["5n2", "2n2", "1n2", "0n2", "2n1", "1n1"],
    "2_5": ["2n4", "0n4", "4n2", "2n2"],
    "1_2": ["5n2", "2n2", "1n2", "0n2", "3n1", "2n1", "1n1"],
    "2_3": ["4n4", "2n4", "0n4", "3n3", "6n2", "5n2", "4n2", "2n2"],
    "3_4": ["6n3", "5n3", "4n3", "6n2", "5n2", "4n2"],
    "4_5": ["6n5", "5n5", "4n5", "3n5"],
    "1_1": ["3n2", "2n2", "1n2", "0n2", "4n1", "3n1", "2n1"],
    "3_2": ["5n5"],
    "5_4": ["10n6", "8n6"],
}


def gap_caption(value: str) -> str:
    a, b = value.split("n")
    return f"Opening: {a}, Extension: {b}"


def _vals(items, caption=None):
    return [{"value": str(v), "caption": caption(v) if caption else str(v)} for v in items]


# Field groups in Ensembl's order. Each field: name, type (dropdown|checklist),
# label, values. Checklists are on/off.
CONFIG_FIELDS = [
    ("general", "General options", [
        {"name": "max_target_seqs", "type": "dropdown",
         "label": "Maximum number of hits to report",
         "values": _vals([10, 50, 100, 250, 500, 1000, 5000])},
        {"name": "evalue", "type": "dropdown",
         "label": "Maximum E-value for reported alignments",
         "values": _vals(["1e-200", "1e-100", "1e-50", "1e-10", "1e-5", "1e-4", "1e-3",
                          "1e-2", "1e-1", "1.0", "10", "100", "1000", "10000", "100000"])},
        {"name": "word_size", "type": "dropdown",
         "label": "Word size for seeding alignments", "values": _vals(range(2, 16))},
        {"name": "max_hsps", "type": "dropdown", "label": "Maximum HSPs per hit",
         "values": _vals([1, 2, 5, 10, 50, 100])},
    ]),
    ("scoring", "Scoring options", [
        {"name": "matrix", "type": "dropdown", "label": "Scoring matrix to use",
         "values": _vals(sorted(MATRIX_GAPS))},
        {"name": "score", "type": "dropdown", "label": "Match/Mismatch scores",
         "values": [{"value": k, "caption": REWARD_LABELS[k]} for k in sorted(REWARD_GAPS)]},
        {"name": "gap_dna", "type": "dropdown", "label": "Gap penalties",
         "depends_on": "score",
         "values_by": {k: _vals(v, gap_caption) for k, v in REWARD_GAPS.items()}},
        {"name": "gappenalty", "type": "dropdown", "label": "Gap penalties",
         "depends_on": "matrix",
         "values_by": {k: _vals(v, gap_caption) for k, v in MATRIX_GAPS.items()}},
        {"name": "ungapped", "type": "checklist", "label": "Disallow gaps in Alignment"},
        {"name": "comp_based_stats", "type": "dropdown", "label": "Compositional adjustments",
         "values": [
             {"value": "0", "caption": "No adjustment"},
             {"value": "1", "caption": "Composition-based statistics"},
             {"value": "2", "caption": "Conditional compositional score matrix adjustment"},
             {"value": "3", "caption": "Universal compositional score matrix adjustment"}]},
        {"name": "threshold", "type": "dropdown",
         "label": "Minimium score to add a word to the BLAST lookup table",
         "values": _vals([11, 12, 13, 14, 15, 16, 20, 999])},
    ]),
    ("filters_and_masking", "Filters and masking options", [
        {"name": "dust", "type": "checklist", "label": "Filter low complexity regions"},
        {"name": "seg", "type": "checklist", "label": "Filter low complexity regions"},
        {"name": "repeat_mask", "type": "checklist",
         "label": "Filter query sequences using RepeatMasker",
         "help": "Query repeats are masked with NCBI WindowMasker for the species."},
    ]),
    ("blat", "BLAT options", [
        {"name": "min_score", "type": "dropdown", "label": "Minimum alignment score",
         "values": _vals([0, 20, 30, 50, 100, 200])},
        {"name": "min_identity", "type": "dropdown", "label": "Minimum sequence identity (%)",
         "values": _vals([0, 25, 50, 75, 80, 85, 90, 95, 98, 100])},
        {"name": "max_intron", "type": "dropdown", "label": "Maximum intron size (bp)",
         "values": _vals([0, 1000, 10000, 50000, 100000, 750000, 2000000])},
    ]),
]
FIELD_BY_NAME = {f["name"]: f for _, _, fields in CONFIG_FIELDS for f in fields}

CONFIG_DEFAULTS = {
    "all": {"evalue": "1e-1", "max_target_seqs": "100"},
    "NCBIBLAST_BLASTN": {"word_size": "11", "score": "1_2", "ungapped": "0",
                         "gap_dna": "5n2", "dust": "1", "repeat_mask": "1",
                         "max_hsps": "100"},
    "NCBIBLAST_BLASTP": {"word_size": "3", "ungapped": "0", "matrix": "BLOSUM62",
                         "gappenalty": "11n1", "threshold": "11", "comp_based_stats": "2",
                         "seg": "1", "max_hsps": "100"},
    "NCBIBLAST_BLASTX": {"word_size": "3", "matrix": "BLOSUM62", "gappenalty": "11n1",
                         "threshold": "11", "seg": "1", "repeat_mask": "1",
                         "max_hsps": "100"},
    "NCBIBLAST_TBLASTN": {"word_size": "3", "ungapped": "0", "matrix": "BLOSUM62",
                          "gappenalty": "11n1", "threshold": "13", "comp_based_stats": "2",
                          "seg": "1", "max_hsps": "100"},
    "NCBIBLAST_TBLASTX": {"word_size": "3", "matrix": "BLOSUM62", "threshold": "13",
                          "seg": "1", "max_hsps": "100"},
    # Ensembl runs gfClient with -minScore=0 -minIdentity=0.
    "BLAT_BLAT": {"min_score": "0", "min_identity": "0", "max_intron": "750000"},
}

SENSITIVITY_OPTIONS = [
    {"value": "near", "caption": "Near match"},
    {"value": "near_oligo", "caption": "Short sequences"},
    {"value": "normal", "caption": "Normal", "selected": True},
    {"value": "distant", "caption": "Distant homologies"},
]
_DNA_SETS = {
    "near":       {"word_size": "15", "dust": "1", "evalue": "10", "score": "1_3", "gap_dna": "5n2"},
    "near_oligo": {"word_size": "7", "dust": "0", "evalue": "1000", "score": "1_3", "gap_dna": "5n2"},
    "normal":     {"word_size": "11", "dust": "1", "evalue": "10", "score": "1_3", "gap_dna": "5n2"},
    "distant":    {"word_size": "9", "dust": "1", "evalue": "10", "score": "1_1", "gap_dna": "2n1"},
}
_PROTEIN_SETS = {
    "near":    {"matrix": "BLOSUM90", "gappenalty": "10n1"},
    "normal":  {"matrix": "BLOSUM62", "gappenalty": "11n1"},
    "distant": {"matrix": "BLOSUM45", "gappenalty": "14n2"},
}
CONFIG_SETS = {
    "NCBIBLAST_BLASTN": _DNA_SETS,
    "NCBIBLAST_BLASTP": _PROTEIN_SETS,
    "NCBIBLAST_BLASTX": _PROTEIN_SETS,
    "NCBIBLAST_TBLASTN": _PROTEIN_SETS,
    "NCBIBLAST_TBLASTX": _PROTEIN_SETS,
}

# Command-line limits BLAST+ enforces on word size.
WORD_SIZE_RANGE = {"blastn": (4, 15), "blastp": (2, 7), "blastx": (2, 7),
                   "tblastn": (2, 7), "tblastx": (2, 3)}


class ConfigError(ValueError):
    pass


def search_caption(search_type: str) -> str:
    s = SEARCH_TYPE_BY_VALUE.get(search_type)
    if not s:
        return search_type
    label = BLAST_TYPES[s["blast_type"]]
    return s["method"] if s["method"] == label else f"{s['method']} ({label})"


def fields_for(search_type: str) -> list[str]:
    """Config field names shown for a search type, in display order."""
    defaults = {**CONFIG_DEFAULTS["all"], **CONFIG_DEFAULTS.get(search_type, {})}
    return [f["name"] for _, _, fields in CONFIG_FIELDS for f in fields
            if f["name"] in defaults]


def defaults_for(search_type: str) -> dict:
    return {k: v for k, v in {**CONFIG_DEFAULTS["all"],
                              **CONFIG_DEFAULTS.get(search_type, {})}.items()}


def allowed_values(search_type: str, name: str, configs: dict) -> list[str] | None:
    """Valid values for a field given the other selections (None = on/off field)."""
    field = FIELD_BY_NAME[name]
    if field["type"] == "checklist":
        return None
    if "values_by" in field:
        parent = configs.get(field["depends_on"])
        return [v["value"] for v in field["values_by"].get(parent, [])]
    values = [v["value"] for v in field["values"]]
    if name == "word_size":
        program = SEARCH_TYPE_BY_VALUE[search_type]["program"]
        lo, hi = WORD_SIZE_RANGE.get(program, (2, 15))
        values = [v for v in values if lo <= int(v) <= hi]
    return values


def normalise_configs(search_type: str, submitted: dict | None,
                      config_set: str = "") -> dict:
    """Validate submitted options against the search type; fill gaps with defaults."""
    if search_type not in SEARCH_TYPE_BY_VALUE:
        raise ConfigError(f"Unknown search tool '{search_type}'.")
    submitted = dict(submitted or {})
    configs = defaults_for(search_type)
    preset = CONFIG_SETS.get(search_type, {}).get(config_set or "", {})
    for key, value in preset.items():
        if key in configs:
            configs[key] = str(value)
    for key in fields_for(search_type):
        if key in submitted and submitted[key] is not None:
            value = submitted[key]
            if isinstance(value, bool):
                value = "1" if value else "0"
            configs[key] = str(value)

    for key in list(configs):
        allowed = allowed_values(search_type, key, configs)
        if allowed is None:
            configs[key] = "1" if configs[key] in ("1", "yes", "true", "on") else "0"
        elif configs[key] not in allowed:
            label = FIELD_BY_NAME[key]["label"]
            if not allowed:
                raise ConfigError(f"No valid values for '{label}'.")
            # A dependent field (gap costs) left over from a different score or
            # matrix silently snaps to the first valid choice, as the form does.
            if "values_by" in FIELD_BY_NAME[key]:
                configs[key] = allowed[0]
            else:
                raise ConfigError(f"'{configs[key]}' is not a valid value for '{label}'.")
    return configs


def describe_configs(search_type: str, configs: dict) -> list[dict]:
    """Human-readable config summary grouped like Ensembl's ticket details."""
    groups = []
    for _key, title, fields in CONFIG_FIELDS:
        rows = []
        for f in fields:
            if f["name"] not in configs:
                continue
            value = configs[f["name"]]
            if f["type"] == "checklist":
                text = "Yes" if value == "1" else "No"
            else:
                pool = f.get("values") or f.get("values_by", {}).get(
                    configs.get(f.get("depends_on", ""), ""), [])
                text = next((v["caption"] for v in pool if v["value"] == value), value)
            rows.append({"label": f["label"], "value": text})
        if rows:
            groups.append({"title": title, "rows": rows})
    return groups


def form_options() -> dict:
    """Static option data for the input form's JavaScript."""
    return {
        "query_types": QUERY_TYPES, "db_types": DB_TYPES,
        "sources": [{"value": k, **v} for k, v in SOURCES.items()],
        "search_types": [{**s, "caption": s["method"]} for s in SEARCH_TYPES],
        "restrictions": RESTRICTIONS,
        "config_fields": [{"key": k, "title": t, "fields": f} for k, t, f in CONFIG_FIELDS],
        "config_defaults": CONFIG_DEFAULTS,
        "config_sets": CONFIG_SETS,
        "sensitivity_options": SENSITIVITY_OPTIONS,
        "word_size_range": WORD_SIZE_RANGE,
        "limits": {"max_sequence_length": MAX_SEQUENCE_LENGTH,
                   "max_num_sequences": MAX_NUM_SEQUENCES,
                   "dna_threshold_percent": DNA_THRESHOLD_PERCENT,
                   "species": SPECIES_SELECTION_LIMIT},
        "blat_value": BLAT_VALUE,
    }


def _gap(value: str) -> tuple[str, str]:
    a, b = value.split("n")
    return a, b


def blast_args(search_type: str, configs: dict, config_set: str = "") -> list[str]:
    """NCBI BLAST+ arguments for the scoring, filtering and reporting options."""
    program = SEARCH_TYPE_BY_VALUE[search_type]["program"]
    c = configs
    args: list[str] = []
    if program == "blastn":
        args += ["-task", "blastn-short" if config_set == "near_oligo" else "blastn"]
    elif program in ("blastp", "blastx", "tblastn"):
        args += ["-task", program]

    args += ["-evalue", c["evalue"], "-max_target_seqs", c["max_target_seqs"]]
    if "max_hsps" in c:
        args += ["-max_hsps", c["max_hsps"]]
    if "word_size" in c:
        args += ["-word_size", c["word_size"]]

    if program == "blastn":
        reward, penalty = c["score"].split("_")
        opening, extension = _gap(c["gap_dna"])
        args += ["-reward", reward, "-penalty", f"-{penalty}",
                 "-gapopen", opening, "-gapextend", extension,
                 "-dust", "yes" if c.get("dust") == "1" else "no"]
    else:
        args += ["-matrix", c["matrix"]]
        if "gappenalty" in c and program != "tblastx":
            opening, extension = _gap(c["gappenalty"])
            args += ["-gapopen", opening, "-gapextend", extension]
        if "threshold" in c:
            args += ["-threshold", c["threshold"]]
        if "comp_based_stats" in c and program in ("blastp", "tblastn", "blastx"):
            args += ["-comp_based_stats", c["comp_based_stats"]]
        args += ["-seg", "yes" if c.get("seg") == "1" else "no"]

    if c.get("ungapped") == "1" and program != "tblastx":
        args.append("-ungapped")
    return args


# ---------------------------------------------------------------------------
# Query sequences
# ---------------------------------------------------------------------------
_VALID_RE = re.compile(rf"^[{SEQUENCE_VALID_CHARS}]*$")


def parse_sequences(raw: str, max_num: int = MAX_NUM_SEQUENCES) -> dict:
    """Split pasted text into sequences using the rules of Ensembl's BLAST form.

    A '>' or a blank line starts a new sequence; digits and whitespace are
    dropped; duplicates and sequences with invalid characters are ignored; RNA
    is converted to DNA; each sequence is typed as DNA when at least 85% of its
    letters are nucleotide codes.
    """
    text = (raw or "").replace("\r\n", "\n").replace("\r", "\n").strip()
    if not text:
        return {"sequences": [], "invalids": 0, "is_rna": False, "errors": []}
    chunks = [c for c in re.split(r"(?=>)|\n[ \t]*\n+", text) if c and c.strip()]
    sequences: list[dict] = []
    seen: set[str] = set()
    invalids = 0
    is_rna = False
    errors: list[str] = []

    for chunk in chunks:
        lines = [ln.strip() for ln in chunk.strip().split("\n")]
        lines = _space_to_newline(lines)
        description = ""
        body = []
        for j, line in enumerate(lines):
            if line.startswith((">", ";")):
                if j == 0:
                    description = line[1:].strip()
                continue
            body.append(re.sub(r"[\d\s]+", "", line.upper()))
        seq = "".join(body)
        if not seq:
            continue
        if not re.fullmatch(r"[A-Z*]+", seq):
            invalids += 1
            bad = "".join(sorted(set(re.sub(r"[A-Z*]", "", seq))))
            errors.append(f"{description or seq[:20]}: invalid characters ({bad})")
            continue
        if len(seq) > MAX_SEQUENCE_LENGTH:
            invalids += 1
            errors.append(f"{description or seq[:20]}: longer than "
                          f"{MAX_SEQUENCE_LENGTH:,} characters")
            continue
        dna_chars = sum(1 for ch in seq if ch in "ACTUGNX")
        seq_type = "dna" if 100 * dna_chars / len(seq) >= DNA_THRESHOLD_PERCENT else "peptide"
        if seq_type == "dna" and "U" in seq:
            is_rna = True
            seq = seq.replace("U", "T")
        if seq in seen:
            invalids += 1
            continue
        seen.add(seq)
        sequences.append({"description": description, "sequence": seq, "type": seq_type,
                          "length": len(seq)})
        if len(sequences) >= max_num:
            break
    return {"sequences": sequences, "invalids": invalids, "is_rna": is_rna,
            "errors": errors}


def _space_to_newline(lines: list[str]) -> list[str]:
    """A FASTA record copied from a web page can arrive with newlines turned to spaces."""
    if len(lines) != 1 or not lines[0].startswith(">"):
        return lines
    words = lines[0].split()
    seq_words: list[str] = []
    only60 = False
    while words:
        w = words.pop()
        if re.fullmatch(r"[A-Za-z*\-]+", w) and (not only60 or len(w) == 60):
            seq_words.insert(0, w)
            only60 = True
        else:
            words.append(w)
            break
    if not seq_words or not words:
        return lines
    return [" ".join(words)] + seq_words


def guess_query_type(sequences: list[dict]) -> str:
    dna = sum(1 for s in sequences if s["type"] == "dna")
    return "dna" if dna >= len(sequences) - dna else "peptide"


def fasta(description: str, sequence: str, width: int = 60) -> str:
    head = ">" + (description or "query")
    return head + "\n" + "\n".join(sequence[i:i + width]
                                   for i in range(0, len(sequence), width)) + "\n"


def is_valid_sequence(seq: str) -> bool:
    return bool(_VALID_RE.match(seq))
