"""Tests for PrimerForge.

    python -m unittest discover -s tests -t .

Everything here runs offline and in a throwaway data directory: no Ensembl
request, no genome, no BLAST+ install. PRIMERFORGE_DATA is redirected before
any part of the app is imported, because config.py reads it at import time and
the whole on-disk layout is derived from it. Importing a test module without
going through this package would write into the real data directory.
"""
from __future__ import annotations

import os
import tempfile
from pathlib import Path

_TMP = Path(tempfile.mkdtemp(prefix="primerforge-tests-"))
os.environ["PRIMERFORGE_DATA"] = str(_TMP)
os.environ.pop("PRIMERFORGE_AUTH", None)          # tests decide their own auth mode

DATA_DIR = _TMP
