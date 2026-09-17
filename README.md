# PrimerForge

A local web app that turns HGVS variant nomenclature into ordered-ready PCR
primers, and keeps every design so you can search, compare and re-export it later.

Paste `NM_000546.6:c.215C>G`, get validated primer pairs with an amplicon map,
thermodynamics, common-SNP warnings and a genome-wide specificity check.

---

## Quick start

Double-click **`run.bat`** (or run `python app.py`). The app opens at
<http://127.0.0.1:5000>.

Nothing else is required — it works immediately using Ensembl web services.
Installing the local genome (Settings → Install local genome) makes it faster,
lets it work offline, and turns on specificity checking.

```
python -m pip install -r requirements.txt
python app.py
```

## Install on a Linux server (shared, with logins)

Everything for the server lives in the **`installation`** folder. One command
installs it on an x86_64 Linux server with systemd (tested on Ubuntu 24.04, Debian
12 and Rocky Linux 9; also written for other Debian/Ubuntu, RHEL/Rocky/Alma/Fedora
and openSUSE releases).

**1. Fill the folder** on a computer with internet access (this one is fine):

```
python installation/prepare.py
```

This copies the application into `installation/app/` and downloads the
dependencies into `installation/dependencies/`, about 360 MB:

| Folder | Contents |
|---|---|
| `dependencies/python/` | Flask, Waitress, primer3-py and their dependencies for Python 3.9–3.14 on x86_64 Linux |
| `dependencies/blast/` | NCBI BLAST+ for x86_64 Linux (checked against NCBI's MD5) |
| `dependencies/blat/` | UCSC `gfServer`, `gfClient`, `blat`, `faToTwoBit`, `twoBitInfo` |
| `dependencies/caddy/` | Caddy, for HTTPS (checked against Caddy's SHA-512) |
| `dependencies/SHA256SUMS` | checksums the installer verifies on the server |

Downloads already present are kept (`--refresh` fetches them again). Run it again
after changing the application so `app/` is current. `--archive` also packs the
folder into `dist/primerforge-installation-<version>.tar.gz`, one file that is
easier to copy.

**2. Copy the whole `installation` folder to the server** (e.g. with WinSCP, or
`scp -r installation you@server:`), then run inside it:

```
cd installation
sudo bash install.sh                                # http://SERVER:8080
sudo bash install.sh --domain primers.mylab.org     # https, Let's Encrypt
sudo bash install.sh --self-signed                  # https, internal network
```

The installer:

- checks the bundled files against `SHA256SUMS`, so a damaged copy stops the
  install instead of failing later;
- installs Python 3.9+ and the libraries BLAST+ and BLAT need from the system's
  package manager, skipping this when they are already present;
- installs the Python packages, BLAST+, the BLAT programs and Caddy from the
  bundled dependencies, downloading only what is missing (so the folder also works
  without `prepare.py`'s downloads, given internet access);
- puts the app in `/opt/primerforge` (its own Python environment), data in
  `/var/lib/primerforge` and settings in `/etc/primerforge/primerforge.env`;
- runs it as the `primerforge` service (Waitress web server, starts at boot,
  restarts on failure) under an unprivileged `primerforge` account;
- turns on logins and creates an `admin` account, printing its password (also
  kept in `/etc/primerforge/initial-admin-password`, readable by root only);
- with `--domain` or `--self-signed`, adds a Caddy HTTPS proxy
  (`primerforge-proxy` service) and opens ports 80/443 in ufw or firewalld, or the
  app port without HTTPS;
- starts downloading the GRCh38 genome and building every human BLAST/BLAT
  database in the background (`primerforge-setup` service). Genome data is not
  bundled: the server needs internet access for it, and for the Ensembl web
  services the app uses anyway. In testing this took
  47 minutes on a fast connection (it depends mostly on download speed from
  Ensembl) and uses 11 GB, about 15 GB at its peak. Primer design works
  immediately through Ensembl's web services; searches become available as each
  database finishes, and the web service restarts itself at the end to pick
  everything up.

Options: `--port N`, `--data-dir DIR`, `--species homo_sapiens,mus_musculus`,
`--no-data` (download later), `--admin NAME`, `--email ADDRESS` (Let's Encrypt
notices), `--no-firewall`. Running `install.sh` from a newer installation folder
upgrades in place and keeps data, accounts and settings. Recommended server: 4+ CPU cores,
16 GB memory (BLAT maps its 2 GB human index into memory and BLAST reads the
genome databases), and 25 GB free disk for human (40 GB with room for more species).

On Linux BLAT runs directly (no WSL): binaries and indexes live in
`<data-dir>/blat`, and gfServer runs as a child of the web service.

### Accounts

With logins on, everyone signs in. **Users** design primers, run BLAST/BLAT and
browse transcripts; they see the whole group's runs and tickets but can delete only
their own. **Administrators** also manage Settings, species installs and accounts.
Administrators add people on the **Users** page, which shows a one-time temporary
password; the person picks their own at first sign-in. Passwords need 10+
characters and are stored as salted hashes. Resetting or changing a password signs
out that account's other sessions, and after 8 failed sign-ins to one account from
one address, further attempts there are paused for 15 minutes.

### Managing the server

```
sudo primerforge status                 services, installed data, accounts
sudo primerforge logs                   follow the web service log
sudo primerforge data-progress          follow the genome/database download
sudo primerforge setup-data             (re)start the download for the configured species
sudo primerforge setup-data --species mus_musculus    add a species now
sudo primerforge users add alice [--admin]            prints a temporary password
sudo primerforge users passwd alice                   e.g. a forgotten password
sudo primerforge users list | role | disable | enable | delete
sudo primerforge restart
sudo primerforge uninstall [--purge]    --purge also deletes data, accounts and settings
```

The desktop app (`run.bat`, `python app.py`) is unchanged: no logins, and it
listens on this computer only.

## What it accepts

| Form | Example |
|---|---|
| Coding HGVS | `NM_000546.6:c.215C>G` |
| Intronic / splice | `NM_000546.6:c.375+1G>A` |
| Deletions, duplications, insertions | `NM_000492.4:c.1521_1523del`, `NM_007294.4:c.5266dup` |
| Genomic HGVS | `NC_000017.11:g.7676154G>C` |
| Plain coordinates | `17-7676154-G-C`, `chr17:7676154 G>C` |
| dbSNP | `rs334` |

One per line, or paste a whole worklist (up to 200 per run). Syntax is checked
as you type, offline, with a message that says what to fix.

## Assays

**Sanger confirmation** — 400–900 bp products with the variant held at least
120 bp from either primer, so it lands in clean sequence well past the noisy
start of the read.

**Long-range / gap-filling** — 1–5 kb products with longer, hotter primers
(25 nt, ~66 °C) suited to long-range polymerases and a two-step protocol.

Every threshold in either preset can be overridden under *Advanced parameters*.

## What each run gives you

- **Resolved variant** — GRCh38 coordinates, gene, MANE transcript, exon or
  intron number, protein consequence and rsID, from Ensembl VEP.
- **Primer pairs** — Primer3 (SantaLucia 1998 thermodynamics) with Tm, GC,
  hairpin and dimer melting points, and each oligo mapped to genomic coordinates.
- **Amplicon map** — exon track, variant position, both primers and the distance
  from the variant to each primer.
- **Common-variant check** — gnomAD allele frequencies for anything under a
  primer, flagged hard when it falls within 5 bp of a 3' end, where it causes
  allele dropout. dbSNP presence alone is not used: it lists a variant at almost
  every position, so only frequency carries signal.
- **Specificity** — every primer is BLASTed genome-wide; hits whose 3' ends
  anneal are paired into predicted amplicons, including forward/forward and
  reverse/reverse combinations. Pairs with a single predicted product rank first.

Export any run as CSV (an order sheet), FASTA (the oligos) or JSON (everything).

## BLAST/BLAT

The **BLAST/BLAT** page is a local version of Ensembl's BLAST/BLAT tool, with the
same inputs, options and result views.

**Input.** Paste up to 30 sequences (plain text or FASTA, up to 200,000 characters
each), upload a FASTA file, or type an identifier: Ensembl stable IDs, RefSeq,
UniProt and ENA accessions are fetched for you. Each sequence is typed as DNA or
protein automatically (RNA is converted to DNA); duplicates and invalid characters
are reported. Choose one or more species (up to 25); one job runs per sequence per
species.

**Databases.** Genomic sequence (unmasked, hard-masked or soft-masked), cDNAs,
Ensembl non-coding RNAs and Ensembl proteins.

**Programs.** BLASTN, BLASTX, BLASTP, TBLASTN, TBLASTX (NCBI BLAST+) and BLAT
(UCSC BLAT in WSL). Only valid query/database/program combinations are offered;
BLAT is chosen by default when it is installed and needs sequences longer than 26 bases.

**Options.** Ensembl's sensitivity presets (Near match, Short sequences, Normal,
Distant homologies) and every option from its configuration panel: maximum hits,
E-value, word size, HSPs per hit, match/mismatch scores, scoring matrix, gap
penalties, ungapped alignment, compositional adjustments, lookup threshold,
low-complexity filtering (DUST/SEG) and query repeat masking (WindowMasker). BLAT
also has minimum score, minimum identity and maximum intron size.

**Tickets.** Recent tickets list every job with its status and hit count. Jobs
queue and run in the background, survive an app restart, and can be cancelled, run
again, edited and resubmitted, or deleted.

**Results.**
- Job details with every setting used.
- Hits drawn on a banded karyotype, coloured by %ID with the best hit boxed.
- HSP distribution along the query.
- A sortable, filterable and pageable results table with selectable columns:
  genomic location, overlapping genes (or gene hit), orientation, query
  coordinates, length, score, E-value and %ID. Clicking a pointer opens a hit
  pop-up with links.
- Per hit, three views, each with a *Configure this page* panel:
  - **Alignment** — match lines or dots, exon and variant markup, line numbering.
  - **Query sequence** — this and other HSPs highlighted.
  - **Genomic sequence** — adjustable flanks and orientation, HSP, exon, variant
    and repeat markup, FASTA download and *BLAST this sequence*.
- Downloads: BLAST pairwise text, tabular, XML, JSON, CSV and SAM, or BLAT's BLAST-
  style text, blast8, PSL and AXT, plus the results table as CSV or TSV.

cDNA and protein hits are mapped back onto the genome through the gene annotation.
Links open the region, gene or transcript on ensembl.org.

**Primer pairs.** Primer pairs on the design results page have *BLAST primers* and
*BLAST/BLAT amplicon* buttons. *BLAST primers* searches each primer as its own
BLASTN job with the "Short sequences" settings, since BLAT needs more than 26 bases.
The results page then shows a **Forward primer** division with a **Reverse primer**
division below it. In each division a *3′ end* column marks where the primer's last
5 bases pair (so it can prime) or the whole primer matches, and a filter shows only
those sites. The karyotype plots only those binding sites, and a summary warns when
a primer matches perfectly at more than one place. Pasted primers named like
`GALE_F` / `GALE_R` or "forward primer" / "reverse primer" get the same layout.

## Transcripts

The **Transcripts** page answers "which transcript should I design against, and what
is its NM accession?". Search a gene symbol, an Ensembl ID (ENSG/ENST/ENSP), a RefSeq
accession (NM/NP), an HGNC ID or an NCBI Gene ID, and every transcript of the gene is
listed with the **MANE Select** one first, alongside its RefSeq and Ensembl accessions,
exon count, cDNA length and protein length.

Opening a transcript gives three tabs, each with its own options and copy/download
buttons; the chosen tab is remembered:

- **Genomic DNA** — the entire genomic sequence (5′ flank, exons, introns, 3′ flank)
  as one continuous colour-coded block, with a marker at the end of each line where
  an exon, intron or flank begins.
- **cDNA** — laid out like Ensembl's transcript cDNA page (see below).
- **Aligned** — genomic, cDNA and protein in register, in the layout sequence viewers
  use for primer design:

```
genomic   23,799,054 aggcagcctcagccacctctgagactctgtatcctgggccagGTGCCATGGCAGAGAAGG 23,798,995
cDNA              -5                                           GTGCCATGGCAGAGAAGG 13
protein            1                                                -M--A--E--K-- 4
```

The genomic and cDNA rows share one colour code: exons alternate blue and violet,
coding bases are bold, UTRs italic, the start codon green and the stop codon red.
Introns are lower case with their splice donor and acceptor sites highlighted, and
flanking sequence is dimmed. The cDNA row carries the RefSeq transcript's own
sequence, and any base where it differs from the reference genome is marked in red.
cDNA is numbered in HGVS c. notation (c.1 is the A of the ATG, c.-N in the 5′ UTR,
c.*N in the 3′ UTR) and each residue sits under the middle base of its codon.
Line numbers give the coordinate and c. position at each end of a row; hovering
shows the exon or intron and splice site, and each residue's p. position.

The cDNA tab reproduces Ensembl's transcript cDNA sequence page, titled e.g.
"GALE Human cDNA". Each block has four lines: variant codes, the cDNA (the RefSeq
accession's own sequence, numbered by transcript position), the coding sequence
(numbered from the A of the ATG) and the protein, one residue under the middle base
of each codon:

```
    *R   K KW *K*****K*******KR**DM*RS***Y*YK****R*Y*Y*YN*******
121 TGGCAGAGAAGGTGCTGGTAACAGGTGGGGCTGGCTACATTGGCAGCCACACGGTGCTGG   180
  2 TGGCAGAGAAGGTGCTGGTAACAGGTGGGGCTGGCTACATTGGCAGCCACACGGTGCTGG    61
  1 M--A--E--K--V--L--V--T--G--G--A--G--Y--I--G--S--H--T--V--L--    20
```

It uses Ensembl's colours and rules. Germline variants and COSMIC somatic mutations
come from Ensembl. Each variant base is underlined and shaded by its consequence
(missense gold, synonymous green, stop gained red, frameshift purple, in-frame indel
pink, splice region coral, coding-sequence/somatic dark green, UTR teal). The code
above it is the IUPAC ambiguity code of a substitution or `*` for anything else,
and it links to the variant in Ensembl. Where variants overlap, the shortest is
shown (and the most severe of equal length), so an SNV inside a long deletion keeps
its own code. Alternate codons are shaded pale yellow, even-numbered exons are blue,
and residues changed by a missense, stop or frameshift variant are red. Hovering a
base lists its cDNA and c. position, exon and every variant over it. Options switch
variants, the coding and protein lines, exon and codon colouring, and the width. On
GALE (NM_001008216) the view matches Ensembl's page base for base.

Options set the 5′ and 3′ flank lengths (up to 5 kb), show introns in full or
shortened, switch to exons only, and change the row width. Numbering can be
chromosome coordinates or relative to the displayed sequence (genomic), and c.
notation or transcript position (the Aligned tab's cDNA row). The protein track and variants (fetched
from Ensembl in the background, marked with IUPAC ambiguity codes) can be turned
on or off. Genes over 1 Mb shorten their introns automatically. The entire genomic
sequence with flanks (never shortened), the cDNA and the protein can be copied or
downloaded as FASTA, the cDNA sent to BLAST/BLAT, or the accession carried into
primer design.

MANE Select and RefSeq accessions also appear beside every transcript and protein hit
in BLAST results, and **RefSeq transcripts (NM/NR)** and **MANE Select transcripts (NM)**
are searchable databases in their own right, so hits can be reported as NM accessions.

Identifiers come from the MANE summary (NCBI/EMBL-EBI), Ensembl's RefSeq
cross-reference table and RefSeq's own transcript FASTA, all held locally.

### Installing databases

Settings → *BLAST/BLAT species databases* installs any Ensembl species from its
current release. The same can be done from a terminal:

```
python tools/setup_species.py homo_sapiens
python tools/setup_species.py mus_musculus --components karyotype,genome,blat,blastdb,gtf
```

Components: `karyotype`, `genome`, `blat`, `blastdb` (one database carrying the
repeat mask, serving unmasked, hard- and soft-masked searches), `gtf` (gene
annotation index), `xrefs` (MANE Select and RefSeq NM/NP identifiers), `cdna`,
`ncrna`, `pep`, `refseq` (all RefSeq transcripts, plus a MANE Select subset) and
`repeats` (WindowMasker statistics). Each is skipped once installed, so re-running
resumes. Human needs about 11 GB for everything, plus 3 GB inside WSL for BLAT.
MANE covers human only; `refseq` needs a per-organism RefSeq set and is skipped for
species that have none.

### BLAT runs in WSL (Windows) or directly (Linux)

On a Linux server BLAT runs natively; see *Install on a Linux server*. On Windows,
UCSC only builds BLAT for Linux and macOS, so it runs inside the Windows Subsystem
for Linux (any distribution; `wsl --install -d Ubuntu` if you have none). The
installer copies `gfServer`, `gfClient`, `blat`, `faToTwoBit` and `twoBitInfo` into
`~/primerforge/blat/bin` in WSL, converts the genome to 2bit and builds a
`gfServer` index (`-stepSize=5`, as UCSC uses) on the Linux filesystem. That is
about a minute for human.

Searches work the way Ensembl runs BLAT: one `gfServer` per species holds the index
in memory and `gfClient` aligns each query. A server starts from its prebuilt index
in about a second, the first time it is needed or when the app starts. It is
listed with start/stop controls under Settings → *BLAT server (WSL)* and stops when
the app exits. The WSL distribution can be chosen there too. BLAT scores and
E-values are computed with BLAT's own formulas, so they match its BLAST-style output.

## Run history

The **Runs** page keeps every design: activity over time, specificity outcomes,
Tm and product-size distributions, most-designed genes, and a searchable,
sortable table. Click any row to reopen the full result.

## Local genome (optional)

Settings → *Install local genome* downloads the soft-masked GRCh38 primary
assembly from Ensembl and NCBI BLAST+, then builds a BLAST database. About 7 GB
on disk. The download is resumable and the app stays usable while it runs.

Soft-masking is what lets the tool warn you when a primer falls in a repeat.

## Layout

```
app.py                  Flask routes, API and exports
serve.py                production server for Linux installs (Waitress)
installation/           everything for a Linux server install:
  install.sh            one-command installer
  prepare.py            copies the app in and downloads the dependencies
  uninstall.sh          removal (also 'sudo primerforge uninstall')
  primerforge.sh        the 'primerforge' management command
  services/             systemd units: web service, data setup, HTTPS proxy
  app/, dependencies/   filled by prepare.py
primerforge/
  config.py             paths, URLs, readiness checks, server settings
  auth.py               accounts, password hashing, sign-in throttling
  cli.py                server administration: users, data setup, status
  net.py                HTTP with retries and backoff
  variants.py           HGVS parsing and resolution to a genomic locus
  refseq.py             local FASTA index, Ensembl fallback
  presets.py            assay presets and parameter validation
  design.py             Primer3 wrapper, coordinate mapping, gnomAD annotation
  specificity.py        BLAST and in-silico PCR
  pipeline.py           run orchestration and progress
  store.py              SQLite runs, results and API cache
  blastconf.py          BLAST/BLAT options, presets and combinations (Ensembl's)
  blastjobs.py          BLAST/BLAT tickets, job queue and persistence
  engines.py            BLAST+ and BLAT runners, parsers and download formats
  blatserver.py         gfServer lifecycle (in WSL or native)
  wsl.py                wsl.exe bridge on Windows; runs commands directly on Linux
  species.py            installed species and their databases
  annotation.py         GTF index: overlapping genes, exons, cDNA/protein -> genome
  transcripts.py        MANE Select, RefSeq/Ensembl accessions, gene and transcript search
  transcriptview.py     cDNA aligned to the genome with its protein translation
  seqviews.py           Alignment, Query sequence and Genomic sequence views
  seqfetch.py           fetch query sequences by identifier
tools/setup_genome.py   genome + BLAST installer
tools/setup_species.py  BLAST/BLAT database installer for any Ensembl species
data/                   SQLite database, genome, BLAST databases, species (generated)
```

## Notes and limits

- Coordinates are GRCh38. The local genome and specificity checking are
  GRCh38-only.
- In repeats, HGVS 3'-shifts an indel while VEP/SPDI left-align it. Where these
  differ the tool shows the HGVS-canonical string and says so; the few-base
  difference does not affect primer placement.
- In-silico PCR predicts products from sequence alone. It does not model
  polymerase, cycling conditions or secondary structure in the template — treat
  it as a strong filter, not a guarantee.
- Designs are not a substitute for wet-lab validation.

## Data sources

Ensembl REST (VEP, Variant Recoder, sequence, assembly, overlap, variation),
Ensembl FTP (genomes, cDNA, ncRNA, peptides, GTF), gnomAD v4 (allele frequencies),
NCBI BLAST+ and WindowMasker statistics, UCSC BLAT, NCBI E-utilities, UniProt and
ENA (sequence fetch), and Primer3 via `primer3-py`.
